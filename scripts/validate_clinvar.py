"""Validate Phase 2 ClinVar retrieval.

Checks:
 1. raw file exists and its size matches metadata
 2. sha256(raw) == metadata.sha256
 3. md5(raw) == metadata.md5 AND == the official ClinVar .md5
 4. parser reproduces the intermediate dataset (re-parse to temp, byte-compare)
 5. raw file was not modified (its current checksum still equals the recorded one)

Exit 0 = all checks pass.
"""
import hashlib
import json
import os
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.clinvar import parse_gene_subset


def _hash(path, algo):
    h = getattr(hashlib, algo)()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024 * 16), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    cfg = yaml.safe_load(open("config/config.yaml"))
    raw_dir = os.path.join(cfg["output_dirs"]["raw"], "clinvar")
    int_dir = cfg["output_dirs"]["intermediate"]
    gene = cfg["study"]["gene"]
    raw = os.path.join(raw_dir, cfg["clinvar"]["file"])
    meta_path = os.path.join(raw_dir, "metadata.json")
    md5_official_path = raw + ".md5"

    problems = []
    meta = json.load(open(meta_path))

    # 1. existence + size
    if not os.path.exists(raw):
        problems.append("raw file missing")
    elif os.path.getsize(raw) != meta["file_size_bytes"]:
        problems.append("raw size mismatch vs metadata")

    # 2. sha256
    if os.path.exists(raw):
        if _hash(raw, "sha256") != meta["sha256"]:
            problems.append("sha256 mismatch (raw modified or corrupt)")

    # 3. md5 (vs metadata + official)
    if os.path.exists(raw):
        md5 = _hash(raw, "md5")
        if md5 != meta["md5"]:
            problems.append("md5 mismatch vs metadata")
        official = open(md5_official_path).read().strip().split()[0]
        if md5 != official:
            problems.append("md5 mismatch vs official ClinVar .md5")

    # 4. reproducibility: re-parse to temp and byte-compare
    out_tsv = os.path.join(int_dir, f"clinvar_{gene.lower()}_raw.tsv")
    tmp_tsv = os.path.join(int_dir, ".reparse_tmp.tsv")
    tmp_summary = os.path.join(int_dir, ".reparse_tmp_summary.json")
    if os.path.exists(raw) and os.path.exists(out_tsv):
        parse_gene_subset(raw, gene, tmp_tsv, tmp_summary)
        if _hash(tmp_tsv, "sha256") != _hash(out_tsv, "sha256"):
            problems.append("re-parse output differs from stored intermediate (not reproducible)")
        os.remove(tmp_tsv)
        os.remove(tmp_summary)

    if problems:
        print("VALIDATION FAILED:")
        for p in problems:
            print("  -", p)
        return 1
    print("VALIDATION PASSED: raw intact, checksums match, parser reproducible.")
    print(f"  raw: {raw}")
    print(f"  sha256: {meta['sha256'][:16]}...  size: {meta['file_size_bytes']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
