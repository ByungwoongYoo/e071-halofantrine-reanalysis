"""Prepare private E071 inputs; do not redistribute database bytes or run analyses.

Default is dry-run. Download is explicit, HTTPS-only, size-bounded and hash-strict.
The manifest describes observed snapshots, not immutable upstream availability or
permission to redistribute. No historical or sealed file is modified.
"""
from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import sys
import urllib.parse
import urllib.request

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
# The public BioMart URL returned this exact official archive in an observed
# HTTPS 308 redirect on 2026-09-27; no wildcard subdomain permission is used.
ALLOWED_HOSTS = {"ftp.ncbi.nlm.nih.gov", "www.ensembl.org",
                 "jun2026.archive.ensembl.org", "maayanlab.cloud"}


class InputError(RuntimeError):
    pass


def sha(data):
    return hashlib.sha256(data).hexdigest()


def relative_path(root, value):
    p = PurePosixPath(value)
    if p.is_absolute() or ".." in p.parts or "\\" in value or ":" in value:
        raise InputError("Unsafe relative input path: " + value)
    root = Path(root).resolve()
    target = (root / value).resolve()
    if target == root or root not in target.parents:
        raise InputError("Input path escapes its root")
    return target


def verify_bytes(item, data):
    actual = sha(data)
    if len(data) != item["bytes"] or actual != item["sha256"]:
        raise InputError("INPUT_HASH_MISMATCH " + item["id"] +
                         " expected=" + item["sha256"] + " actual=" + actual)
    if item.get("format") == "gzip":
        # Reading to EOF validates the gzip CRC; expansion is independently bounded.
        with gzip.GzipFile(fileobj=io.BytesIO(data)) as stream:
            total = 0
            while chunk := stream.read(1024 * 1024):
                total += len(chunk)
                if total > 64 * 1024 * 1024:
                    raise InputError("Unexpected gzip expansion")
    return data


def checked_url(url):
    parsed = urllib.parse.urlsplit(url)
    if (parsed.scheme != "https" or parsed.hostname not in ALLOWED_HOSTS or
            parsed.username or parsed.password or parsed.port not in (None, 443)):
        raise InputError("Non-official or unsafe acquisition URL")
    return url


class OfficialRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        checked_url(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def download(url, maximum_bytes):
    checked_url(url)
    request = urllib.request.Request(url, headers={"User-Agent": "E071-input-replay/1.0"})
    opener = urllib.request.build_opener(OfficialRedirect())
    with opener.open(request, timeout=45) as response:
        checked_url(response.url)
        if response.status != 200:
            raise InputError("Unexpected response status")
        body = response.read(maximum_bytes + 1)
        if len(body) > maximum_bytes:
            raise InputError("Input exceeds declared download bound")
        # Exact expected digest rejects challenges; three legacy index assets are HTML.
        return body


def query_xml(kind):
    attrs = ["ensembl_gene_id"]
    if kind == "gsea":
        attrs.append("external_gene_name")
    attrs += ["mmusculus_homolog_ensembl_gene"]
    if kind != "all":
        attrs.append("mmusculus_homolog_orthology_type")
    attrs.append("rnorvegicus_homolog_ensembl_gene")
    if kind != "all":
        attrs.append("rnorvegicus_homolog_orthology_type")
    count = '' if kind == "core" else ' count=""'
    return ('<?xml version="1.0" encoding="UTF-8"?>\n'
            '<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1"' +
            count + ' datasetConfigVersion="0.6">\n'
            '<Dataset name="hsapiens_gene_ensembl" interface="default">\n' +
            ''.join('<Attribute name="' + a + '"/>\n' for a in attrs) +
            '</Dataset></Query>')


def map_response(data, kind):
    """Port only the existing BioMart filter/serialization, not analysis arithmetic."""
    rows = list(csv.reader(io.StringIO(data.decode("utf-8")), delimiter="\t"))
    names = (["human", "symbol", "mouse", "mouse_type", "rat", "rat_type"]
             if kind == "gsea" else ["human", "mouse", "mouse_type", "rat", "rat_type"])
    if not rows or len(rows[0]) != len(names):
        raise InputError("Unexpected BioMart schema")
    if any(len(r) != len(names) for r in rows[1:]):
        raise InputError("Ragged BioMart response")
    retained = [dict(zip(names, r)) for r in rows[1:]]
    retained = [r for r in retained if r["mouse_type"] == r["rat_type"] == "ortholog_one2one"]
    keys = ["human", "symbol", "mouse", "rat"] if kind == "gsea" else ["human", "mouse", "rat"]
    for r in retained:
        for k in keys:
            r[k] = re.sub(r"\.\d+$", "", r[k]).strip()
    retained = [r for r in retained if all(r[k] for k in keys)]
    # Sequential uniqueness filtering matches e071_freeze_ortholog_maps.py.
    for k in ("human", "mouse", "rat"):
        counts = collections.Counter(r[k] for r in retained)
        retained = [r for r in retained if counts[r[k]] == 1]
    values = sorted(set(tuple(r[k] for k in keys) for r in retained))
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(keys)
    writer.writerows(values)
    return output.getvalue().encode("utf-8")


def normalized_gmt(raw):
    rows = []
    seen = set()
    for line in raw.decode("utf-8").splitlines():
        fields = line.strip().split("\t")
        if len(fields) < 2:
            continue
        if fields[0] in seen:
            raise InputError("Duplicate GMT term")
        seen.add(fields[0])
        genes = [g.split(",")[0] for g in fields[2:] if g.split(",")[0]]
        rows.append("\t".join(fields[:2] + genes))
    # Existing Windows replay used write_text newline translation. Make it explicit.
    return ("\r\n".join(rows) + "\r\n").encode("utf-8")


def load_manifest(path):
    manifest = json.loads(Path(path).read_text(encoding="utf-8"))
    ids = [x["id"] for x in manifest["inputs"]]
    if len(ids) != len(set(ids)):
        raise InputError("Duplicate input IDs")
    for item in manifest["inputs"]:
        if not re.fullmatch(r"[0-9a-f]{64}", item["sha256"]):
            raise InputError("Invalid manifest digest")
        if item.get("url"):
            checked_url(item["url"])
        for p in item["targets"]:
            relative_path(HERE, p)
    return manifest


def resolve_inputs(manifest, mode, cache_root=None, selected=None, allow_mutable=False, fetch=download):
    items = {x["id"]: x for x in manifest["inputs"]}
    requested = list(items) if not selected else selected
    if any(x not in items for x in requested):
        raise InputError("Unknown selected input")
    ready = {}
    origins = {}

    def resolve(key):
        if key in ready:
            return ready[key]
        item = items[key]
        if cache_root:
            for value in item["targets"]:
                path = relative_path(cache_root, value)
                if path.is_file():
                    data = verify_bytes(item, path.read_bytes())
                    ready[key], origins[key] = data, "verified-local-cache"
                    return data
        if mode == "cache-check":
            raise InputError("CACHE_INPUT_MISSING " + key)
        if item.get("mutable") and not allow_mutable:
            raise InputError("MUTABLE_SOURCE_REQUIRES_EXPLICIT_FLAG " + key)
        recipe = item["recipe"]
        if recipe == "query-xml":
            data = query_xml("all").replace("\n", "\r\n").encode("utf-8")
        elif recipe == "normalize-gmt":
            data = normalized_gmt(resolve(item["depends_on"]))
        elif recipe in ("biomart-all", "biomart-core", "biomart-gsea"):
            kind = recipe.split("-")[1]
            url = item["url"] + "?" + urllib.parse.urlencode({"query": query_xml(kind)})
            raw = fetch(url, item["max_response_bytes"])
            if kind == "all":
                data = raw
            else:
                canonical = map_response(raw, kind)
                if sha(canonical) != item["canonical_lf_sha256"]:
                    raise InputError("MAPPING_CONTENT_DRIFT " + key)
                data = canonical.replace(b"\n", b"\r\n")
        elif recipe == "official-bytes":
            data = fetch(item["url"], item["bytes"] + 4096)
        else:
            raise InputError("Unsupported acquisition recipe")
        data = verify_bytes(item, data)
        ready[key], origins[key] = data, recipe
        return data

    for key in requested:
        resolve(key)
    return ready, origins


def materialize(manifest, ready, output_root):
    root = Path(output_root).resolve()
    # Inputs must remain outside this potentially distributable source directory.
    if root == HERE or HERE in root.parents:
        raise InputError("Private inputs must be outside public_candidate")
    entries = [(item, relative_path(root, target)) for item in manifest["inputs"]
               if item["id"] in ready for target in item["targets"]]
    for item, target in entries:
        if target.exists():
            raise InputError("Refusing to overwrite existing input: " + str(target))
        verify_bytes(item, ready[item["id"]])
    created = []
    try:
        for item, target in entries:
            target.parent.mkdir(parents=True, exist_ok=True)
            with target.open("xb") as out:
                created.append(target)
                out.write(ready[item["id"]])
            verify_bytes(item, target.read_bytes())
    except Exception:
        # Only files exclusively created by this call; never user originals.
        for target in created:
            if root in target.resolve().parents:
                target.unlink(missing_ok=True)
        raise
    return len(entries)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=HERE / "input_manifest.json")
    parser.add_argument("--mode", choices=["dry-run", "cache-check", "download"], default="dry-run")
    parser.add_argument("--cache-root", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--only", nargs="+")
    parser.add_argument("--allow-mutable", action="store_true")
    args = parser.parse_args(argv)
    manifest = load_manifest(args.manifest)
    if args.mode == "dry-run":
        result = {"status": "PLAN_ONLY_NO_NETWORK_NO_WRITES", "inputs": manifest["inputs"]}
    else:
        if args.mode == "cache-check" and args.cache_root is None:
            raise InputError("cache-check requires --cache-root")
        ready, origins = resolve_inputs(manifest, args.mode, args.cache_root, args.only, args.allow_mutable)
        written = materialize(manifest, ready, args.output_dir) if args.output_dir else 0
        result = {"status": "EXACT_INPUT_BYTES_VERIFIED", "input_count": len(ready),
                  "complete_input_set": len(ready) == len(manifest["inputs"]),
                  "origins": origins, "files_written": written,
                  "manifest_sha256": sha(args.manifest.read_bytes()),
                  "analysis_executed": False, "release_rights_cleared": False}
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    try:
        main()
    except (InputError, OSError, ValueError) as exc:
        print(json.dumps({"status": "FAIL_CLOSED", "error": str(exc)}), file=sys.stderr)
        sys.exit(2)
