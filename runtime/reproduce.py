"""Run the complete E071 reanalysis from a selected-source candidate.

Default: integrity check and execution plan only. --execute creates a NEW,
external workspace. Inputs are hash checked; download needs --fetch-inputs.
This version supports the verified Windows / Python 3.12.10 runtime only.
The original scientific scripts and expected results are not modified.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path, PurePosixPath
import platform
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
PACKAGE = "rescience_offline_candidate_20260927"
OUTPUTS = {"verification.json", "itpr3_claim_verification.csv", "itpr3_sample_values.csv",
           "source_manifest.json", "original_itpr3_replay.json", "original_nine_rule_replay.json"}
RESULT_HASHES = {"core_seed1": "a52f57efdc4c82a08264ce4b17ce9c3eed45ff079bec1d2ed21df5d93c1cb221",
                 "core_seed77": "a52f57efdc4c82a08264ce4b17ce9c3eed45ff079bec1d2ed21df5d93c1cb221",
                 "gsea_current": "3bba0847db786c4a591d73263c2d95355768146ebdd83600374a7241c2d2d790",
                 "loo_frozen": "29e349d0ef9198d15a1b085e07cb78790ce3e471e321e1b74cb708677c1a3aa6"}
SCOPE = {"loo_deletions": 25, "cross_model_state_pairs": 257,
         "label_assignments": {"CEP290_LCA": 252, "P23H": 126, "rd10": 20},
         "gsea_models": 3, "gsea_tested_pathways_per_model": 277,
         "gsea_permutations_per_model": 2000, "overlap_draws_per_test": 10000}
VERSIONS = {"python": "3.12.10", "numpy": "2.5.3", "pandas": "2.3.3",
            "scipy": "1.18.1", "requests": "2.34.2", "gseapy": "1.3.1", "matplotlib": "3.11.2"}


class ReplayError(RuntimeError):
    pass


def require(ok, reason):
    if not ok:
        raise ReplayError(reason)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def load(path):
    def invalid(value):
        raise ReplayError("Non-finite JSON: " + value)
    return json.loads(Path(path).read_text(encoding="utf-8"), parse_constant=invalid)


def relative(root, name):
    p = PurePosixPath(name)
    require(bool(name) and not p.is_absolute() and ".." not in p.parts and
            "\\" not in name and ":" not in name, "Unsafe relative path")
    root = Path(root).resolve()
    target = root / name
    for node in [target, *target.parents]:
        if node == root:
            break
        require(not node.is_symlink() and not node.is_junction(), "Linked path is not allowed")
    require(root in target.resolve().parents, "Path escapes root")
    return target


def checked(data, row):
    require(len(data) == row["bytes"] and sha(data) == row["sha256"],
            "File differs from contract: " + row.get("path", row.get("id", "input")))
    return data


def write_new(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write(data)
    require(path.read_bytes() == data, "Write readback failed")


def write_json(path, value):
    write_new(path, (json.dumps(value, indent=2, allow_nan=False) + "\n").encode("utf-8"))


def preflight(root):
    root = Path(root).resolve()
    contract = load(root / "release_manifest.json")
    require(contract["schema_version"] == 1 and contract["package_directory"] == PACKAGE,
            "Unsupported package contract")
    seen = set()
    for row in contract["files"]:
        name = row["path"]
        require(name not in seen, "Duplicate manifest file")
        seen.add(name)
        require(not any(part in {"deps", "runs", "evidence", "__pycache__"}
                        for part in PurePosixPath(name).parts), "Private file in release contract")
        require(not name.lower().endswith((".gmt", ".pyd", ".exe", "tested_gene_sets.json")),
                "Private input or platform binary in release contract")
        checked(relative(root, name).read_bytes(), row)
    actual_files = {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file()}
    require(actual_files == seen | {"release_manifest.json"}, "Unlisted or missing candidate files")
    inputs = load(root / "input_manifest.json")
    require(sha((root / "input_manifest.json").read_bytes()) == contract["input_manifest_sha256"],
            "Input contract changed")
    require(len(inputs["inputs"]) == 15 and len({x["id"] for x in inputs["inputs"]}) == 15 and
            sum(len(x["targets"]) for x in inputs["inputs"]) == 16 and
            len({p for x in inputs["inputs"] for p in x["targets"]}) == 16,
            "Incomplete scientific input scope")
    require({"reproduce.py", "acquire_inputs.py", "input_manifest.json",
             "build_rescience_figures.py", PACKAGE + "/run.py",
             PACKAGE + "/full_replay/replay.py", PACKAGE + "/full_replay/compare.py"} <= seen,
            "Missing runtime entry point")
    bounded_seen = set()
    for row in contract["bounded_manifest"]:
        relative(root / PACKAGE, row["path"])
        require(row["path"] not in bounded_seen, "Duplicate bounded manifest file")
        bounded_seen.add(row["path"])
    require({"run.py", "verify_e071.py", "replay_existing.py"} <= bounded_seen and
            {"reference_outputs/" + x for x in OUTPUTS} <= bounded_seen,
            "Bounded contract omits required source or result files")
    return contract, inputs


def new_workspace(root, work):
    work = Path(work)
    require(work.is_absolute() and not work.exists(), "Workspace must be a NEW absolute path")
    require(work.parent.is_dir() and len(work.name) >= 3, "Workspace parent must already exist")
    require(not work.is_symlink() and not work.is_junction(), "Linked workspace refused")
    resolved = work.resolve()
    require(root != resolved and root not in resolved.parents, "Keep inputs outside the source package")
    require(resolved != Path(resolved.anchor), "Drive root is not a workspace")
    return resolved


def runtime(python):
    code = ("import json,platform,importlib.metadata as m; "
            "print(json.dumps(dict(system=platform.system(),python=platform.python_version(),"
            "**{n:m.version(n) for n in ['numpy','pandas','scipy','requests','gseapy','matplotlib']})))")
    result = subprocess.run([python, "-B", "-c", code], check=True, capture_output=True, text=True)
    values = json.loads(result.stdout)
    require(values.get("system") == "Windows", "Only the verified Windows serialization is supported")
    require(all(values.get(k) == v for k, v in VERSIONS.items()), "Pinned runtime does not match")
    return values


def prepare(root, work, contract, inputs, cache=None, fetch=False, allow_mutable=False):
    work = new_workspace(root, work)
    work.mkdir()
    write_new(work / ".codex-task-temp", b"E071 isolated reproducibility run. Preserve until reviewed.\n")
    for row in contract["files"]:
        write_new(relative(work, row["path"]), checked(relative(root, row["path"]).read_bytes(), row))
    pack = work / PACKAGE
    provided_cache_ids = []
    for item in inputs["inputs"]:
        available = None
        for name in item["targets"]:
            target = relative(pack, name)
            if target.is_file():
                available = checked(target.read_bytes(), item)
                break
        if available is None and cache:
            for name in item["targets"]:
                source = relative(cache, name)
                if source.is_file():
                    available = checked(source.read_bytes(), item)
                    provided_cache_ids.append(item["id"])
                    break
        if available is not None:
            for name in item["targets"]:
                target = relative(pack, name)
                if target.exists():
                    checked(target.read_bytes(), item)
                else:
                    write_new(target, available)
    spec = importlib.util.spec_from_file_location("e071_selected_acquisition", root / "acquire_inputs.py")
    acquire = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(acquire)
    ready, origins = acquire.resolve_inputs(inputs, "download" if fetch else "cache-check",
                                            cache_root=pack, allow_mutable=allow_mutable)
    require(len(ready) == 15, "Not every required input was acquired")
    for item in inputs["inputs"]:
        data = checked(ready[item["id"]], item)
        for name in item["targets"]:
            target = relative(pack, name)
            if target.exists():
                checked(target.read_bytes(), item)
            else:
                write_new(target, data)
    write_json(pack / "MANIFEST.json", {"scope": "Selected-source bounded stage within the full pipeline",
                                        "files": contract["bounded_manifest"]})
    for row in contract["bounded_manifest"]:
        checked(relative(pack, row["path"]).read_bytes(), row)
    return pack, {"origins": origins, "provided_cache_ids": provided_cache_ids,
                  "input_count": 15, "input_targets": 16}


def graph(work, python):
    pack = Path(work) / PACKAGE
    full = pack / "full_replay"
    result = [("bounded", [python, "-B", str(pack / "run.py"), "--run-name", "bounded_current"], "1")]
    for action, name, seed in [("core", "core_seed1", "1"), ("core", "core_seed77", "77"),
                               ("gsea", "gsea_current", "1"), ("loo", "loo_frozen", "1")]:
        result.append((name, [python, "-B", str(full / "replay.py"), action, "--run-name", name], seed))
    result += [("compare", [python, "-B", str(full / "compare.py")], "1"),
               ("figures", [python, "-B", str(Path(work) / "build_rescience_figures.py"),
                            "--out", str(Path(work) / "figures"), "--preview-dir", str(work)], "1")]
    return result


def validate_compare(report):
    require(report.get("status") == "PASS_RETAINED_DETERMINISTIC_METRICS", "Comparison failed")
    for key, n, flag in [("historical_checks", 26, "match"),
                          ("loo_and_permutation_rounded_checks", 27, "match"),
                          ("highlighted_pathway_display_checks", 18, "match"),
                          ("overlap_historical_display_comparison", 5, "rounds_to_historical")]:
        rows = report.get(key, [])
        require(len(rows) == n and all(x.get(flag) is True for x in rows), "Incomplete or failed " + key)
    seeds = report.get("hashseed_determinism", {})
    require(seeds.get("pythonhashseed_values") == ["1", "77"] and
            seeds.get("full_core_json_equal") is True and seeds.get("byte_identical") is True,
            "Seed reproducibility failed")
    require(report.get("scope") == SCOPE, "Scientific scope differs")


def validate_stage(name, work, contract):
    pack = work / PACKAGE
    full = pack / "full_replay"
    if name == "bounded":
        out = pack / "runs/bounded_current"
        receipt = load(out / "REPLAY_RECEIPT.json")
        require(receipt.get("status") == "PASS_OFFLINE_SIX_OUTPUT_BYTE_IDENTITY" and
                receipt.get("analysis_network_attempts") == 0, "Bounded stage failed")
        rows = receipt.get("comparison", [])
        require(len(rows) == 6 and {x["file"] for x in rows} == OUTPUTS, "Missing bounded output")
        for row in rows:
            data = (out / row["file"]).read_bytes()
            require(row.get("byte_identical") is True and sha(data) == row["sha256"] == row["expected_sha256"]
                    and data == (pack / "reference_outputs" / row["file"]).read_bytes(), "Bounded byte mismatch")
    elif name in RESULT_HASHES:
        out = full / "runs" / name
        rec = load(out / "RUN_RECEIPT.json")
        actual = sha((out / "result.json").read_bytes())
        require(actual == RESULT_HASHES[name] == rec["result_sha256"] and
                rec["status"] == "EXECUTED_OFFLINE" and rec["analysis_network_attempts"] == 0 and
                rec["calculation_rules_modified"] is False and
                rec["pythonhashseed"] == ("77" if name == "core_seed77" else "1"), "Full stage differs")
    elif name == "compare":
        validate_compare(load(full / "FULL_REPLAY_COMPARISON.json"))
        write_new(work / "results/itpr3_claim_verification.csv",
                  (pack / "runs/bounded_current/itpr3_claim_verification.csv").read_bytes())
    elif name == "figures":
        rec = load(work / "figures/FIGURE_PROVENANCE.json")
        expected = {"core": full / "runs/core_seed1/result.json", "loo": full / "runs/loo_frozen/result.json",
                    "gsea": full / "runs/gsea_current/result.json",
                    "itpr3": pack / "runs/bounded_current/itpr3_claim_verification.csv"}
        require(rec["input_hashes"] == {k: sha(p.read_bytes()) for k, p in expected.items()}, "Figure input mismatch")
        require(len(rec["figures"]) == 5, "Missing figure")
        for row in rec["figures"]:
            require(sha(relative(work / "figures", row["svg"]).read_bytes()) == row["sha256"], "Figure mismatch")
    for row in contract["files"]:
        checked(relative(work, row["path"]).read_bytes(), row)
    for row in contract["bounded_manifest"]:
        checked(relative(pack, row["path"]).read_bytes(), row)


def validate_inputs(pack, inputs):
    count = 0
    for item in inputs["inputs"]:
        for target in item["targets"]:
            checked(relative(pack, target).read_bytes(), item)
            count += 1
    require(count == 16, "Incomplete post-stage input validation")
    return count


def execute(root, work, contract, inputs, python, cache=None, fetch=False, allow_mutable=False):
    versions = runtime(python)
    pack, acquisition = prepare(root, work, contract, inputs, cache, fetch, allow_mutable)
    work = pack.parent
    completed = []
    for name, command, seed in graph(work, python):
        env = os.environ.copy()
        env.update(PYTHONHASHSEED=seed, PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1",
                   MPLCONFIGDIR=str(work / "matplotlib_cache"))
        env.pop("PYTHONPATH", None)
        with (work / (name + ".log")).open("x", encoding="utf-8") as log:
            subprocess.run(command, cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT,
                           check=True, timeout=600)
        validate_stage(name, work, contract)
        validate_inputs(pack, inputs)
        completed.append(name)
    preflight(root)
    result = {"status": "PASS_SELECTED_SOURCE_FULL_PIPELINE", "completed": completed,
              "runtime": versions, "scope": SCOPE, "acquisition": acquisition,
              "source_hashes": {x["path"]: x["sha256"] for x in contract["files"]},
              "release_manifest_sha256": sha((root / "release_manifest.json").read_bytes()),
              "result_hashes": RESULT_HASHES, "scientific_sources_modified": False,
              "license_changed_by_this_run": False, "repository_published_by_this_run": False, "submission_performed_by_this_run": False,
              "privacy": "This receipt uses relative source paths. Execution logs, figures' private provenance and generated KEGG memberships stay in the local workspace and are not release files.",
              "historical_gmt_identity": "NOT_ESTABLISHED", "visual_review": "NOT_PERFORMED_BY_RUNNER"}
    write_json(work / "PORTABLE_EXECUTION_RECEIPT.json", result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--input-cache", type=Path)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--fetch-inputs", action="store_true")
    parser.add_argument("--allow-mutable", action="store_true")
    args = parser.parse_args(argv)
    contract, inputs = preflight(HERE)
    work = new_workspace(HERE, args.work_root)
    if args.execute:
        result = execute(HERE, work, contract, inputs, args.python,
                         args.input_cache, args.fetch_inputs, args.allow_mutable)
    else:
        result = {"status": "DRY_RUN_NO_WRITES_NO_NETWORK_NO_ANALYSIS",
                  "contract_files_verified": len(contract["files"]),
                  "input_assets": 15, "input_targets": 16, "commands": graph(work, args.python),
                  "required_runtime": {"system": "Windows", **VERSIONS}, "source_scope": SCOPE,
                  "repository_published_by_this_run": False, "license_changed_by_this_run": False}
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    try:
        main()
    except (ReplayError, OSError, ValueError, KeyError, subprocess.SubprocessError) as error:
        print(json.dumps({"status": "FAIL_CLOSED", "error": str(error)}), file=sys.stderr)
        sys.exit(2)
