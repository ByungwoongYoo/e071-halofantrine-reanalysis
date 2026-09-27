"""Integrity-checked offline runner for the bounded 2026-09-27 E071 replay.

Preserves the existing calculations, never overwrites a run directory, and
denies Python socket operations before executing the unchanged scripts.
"""
from pathlib import Path
import argparse
import contextlib
import datetime as dt
import hashlib
import io
import json
import runpy
import shutil
import socket
import sys

ROOT = Path(__file__).resolve().parent
BASE = Path("upstream/02_COORDINATOR/BIOASSET_DEEP_IP_20260919")
OUTPUTS = ["verification.json", "itpr3_claim_verification.csv",
           "itpr3_sample_values.csv", "source_manifest.json",
           "original_itpr3_replay.json", "original_nine_rule_replay.json"]
AUTHORED = ["run.py", "README.md", "requirements.txt", "LICENSE_STATUS.json"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def prepare(source):
    """Copy only the explicit frozen-input contract; do not alter originals."""
    source = source.resolve()
    manifest_path = ROOT / "MANIFEST.json"
    if manifest_path.exists():
        raise FileExistsError("Package already sealed; refusing to replace its manifest")
    source_manifest = json.loads((source / "results/source_manifest.json").read_text())
    names = {entry[role]["filename"] for entry in source_manifest["sources"].values()
             for role in ("deg", "cpm")}
    names.update([f"GSE{n}_suppl_index.html" for n in (327431, 327432, 327434)])
    names.update(["ensembl_all_homolog_query_20260927.xml",
                  "ensembl_all_homolog_query_result_20260927.tsv"])
    pairs = [(Path(n), Path(n)) for n in ("verify_e071.py", "replay_existing.py")]
    pairs += [(Path("source_data") / n, Path("source_data") / n) for n in sorted(names)]
    pairs += [(BASE / n, BASE / n) for n in
              ("e071_itpr3_source_diagnostic.py", "e071_venn_rule_diagnostic.py",
               "SREP_SUBMISSION_PACKAGE_V1/ortholog_mapping/ensembl_one2one_core.tsv")]
    pairs += [(Path("results") / n, Path("reference_outputs") / n) for n in OUTPUTS]
    pairs += [(Path("results/source_manifest.json"), Path("results/source_manifest.json"))]
    pairs += [(Path("results") / n, Path("evidence") / n) for n in
              ("verification_audit.json", "E071_RECEIPT.json", "decisive_claim_table.json")]
    assert len(pairs) == len({str(dst) for _, dst in pairs})
    for src, dst in pairs:
        assert (source / src).is_file(), src
        assert not (ROOT / dst).exists(), dst
    copied = []
    for src, dst in pairs:
        target = ROOT / dst
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source / src, target)
        assert digest(source / src) == digest(target), dst
        copied.append({"path": dst.as_posix(), "source_relative_path": src.as_posix(),
                       "bytes": target.stat().st_size, "sha256": digest(target)})
    files = copied + [{"path": n, "bytes": (ROOT / n).stat().st_size,
                       "sha256": digest(ROOT / n), "origin": "new packaging only"}
                      for n in AUTHORED]
    write_json(manifest_path, {"source_commit": "cdbda819817b7f1005c7b0dcfa2489468347b10e",
                              "scope": "bounded September 27 replay, not entire V29",
                              "files": files, "copied_file_count": len(copied)})
    print(f"SEALED {len(files)} files; {len(copied)} exact copies")


def verify_integrity():
    manifest = json.loads((ROOT / "MANIFEST.json").read_text())
    for item in manifest["files"]:
        path = ROOT / item["path"]
        assert path.is_file(), item["path"]
        assert path.stat().st_size == item["bytes"], item["path"]
        assert digest(path) == item["sha256"], item["path"]
    return len(manifest["files"])


def execute(out_name):
    count = verify_integrity()
    out = (ROOT / "runs" / out_name).resolve()
    assert out.parent == (ROOT / "runs").resolve(), "Use a single run name"
    out.mkdir(parents=True, exist_ok=False)
    sys.dont_write_bytecode = True
    blocked = []
    loopback_probes = []

    def deny_network(event, args):
        # urllib3 checks IPv6 support by binding an ephemeral loopback port.
        # That local feature probe sends no traffic and is not a remote request.
        if event == "socket.bind" and args[1] in (("::1", 0), ("127.0.0.1", 0)):
            loopback_probes.append(event)
            return
        if event in {"socket.connect", "socket.getaddrinfo", "socket.gethostbyname",
                     "socket.gethostbyaddr", "socket.sendto", "socket.bind"}:
            blocked.append(event)
            raise RuntimeError("Offline replay prohibits network: " + event)

    sys.addaudithook(deny_network)
    try:
        socket.getaddrinfo("offline-guard-self-test.invalid", 443)
    except RuntimeError as exc:
        assert "Offline replay prohibits network" in str(exc)
    else:
        raise AssertionError("Network-denial self-test failed")
    assert blocked == ["socket.getaddrinfo"]
    for name in ("verify_e071.py", "replay_existing.py"):
        capture = io.StringIO()
        sys.argv = [str(ROOT / name), "--offline", "--output-dir", str(out)]
        with contextlib.redirect_stdout(capture):
            runpy.run_path(str(ROOT / name), run_name="__main__")
        (out / (name + ".log")).write_text(capture.getvalue(), encoding="utf-8")
    assert blocked == ["socket.getaddrinfo"], "A calculation attempted network access"
    comparison = []
    for name in OUTPUTS:
        actual, expected = out / name, ROOT / "reference_outputs" / name
        equal = actual.read_bytes() == expected.read_bytes()
        comparison.append({"file": name, "byte_identical": equal,
                           "sha256": digest(actual), "expected_sha256": digest(expected)})
        assert equal, f"Output differs from frozen September 27 reference: {name}"
    assert verify_integrity() == count
    result = {"status": "PASS_OFFLINE_SIX_OUTPUT_BYTE_IDENTITY",
              "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "python": sys.version, "sealed_files_verified_before_and_after": count,
              "manifest_sha256": digest(ROOT / "MANIFEST.json"),
              "network_guard_self_test": "PASS", "analysis_network_attempts": 0,
              "allowed_local_ipv6_feature_probes": len(loopback_probes),
              "comparison": comparison,
              "not_claimed": ["fresh raw-count limma refit", "exact source S7A/S7B reproduction",
                              "fresh GSEA/overlap-null/full LOO replay", "license clearance",
                              "full V29 replay", "independent reviewer", "submission readiness"]}
    write_json(out / "REPLAY_RECEIPT.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-from", type=Path)
    parser.add_argument("--run-name", default="replay_1")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    if args.prepare_from:
        prepare(args.prepare_from)
    elif args.verify_only:
        print("INTEGRITY_PASS", verify_integrity())
    else:
        execute(args.run_name)
