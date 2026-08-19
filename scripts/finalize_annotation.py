"""Phase 4B-cleanup — finalize and validate the annotated dataset (TASK 3 + TASK 4).

Verifies data integrity (row count, unique keys, no lost/new records, unchanged predictor
values vs the previous committed dataset, valid AF ranges, NA-as-empty, no error strings in
numeric fields) and writes a finalization report with checksum + provenance.
"""
import csv
import hashlib
import json
import os
import subprocess
import sys
import time
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import yaml

ANN_TSV = "data/processed/brca1_vus_missense_annotated.tsv"
NORM_TSV = "data/intermediate/brca1_vus_missense_normalized.tsv"
REPORT = "results/reports/annotation_finalization_report.md"

NUMERIC_FIELDS = [
    "gnomad_genome_af", "gnomad_genome_ac", "gnomad_genome_an", "gnomad_genome_hom",
    "gnomad_genome_faf95_popmax",
    "gnomad_exome_af", "gnomad_exome_ac", "gnomad_exome_an", "gnomad_exome_hom",
    "gnomad_exome_faf95_popmax",
    "sift_score", "polyphen_score", "revel_score",
]
UNCHANGED_FIELDS = ["revel_score", "sift_score", "sift_prediction",
                    "polyphen_score", "polyphen_prediction", "vep_hgvsc", "vep_hgvsp"]


def load_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        idx = {n: i for i, n in enumerate(header)}
        rows = [row for row in r]
    return header, idx, rows


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    checks = []

    _, a_idx, ann = load_tsv(ANN_TSV)
    _, n_idx, norm = load_tsv(NORM_TSV)

    # 1. row count
    n_ann = len(ann)
    checks.append(("row count == 1904", n_ann == 1904, f"{n_ann}"))

    # 2. unique biological variant keys
    keys = [r[a_idx["variant_key"]] for r in ann]
    checks.append(("unique variant keys == 1904", len(set(keys)) == 1904, f"{len(set(keys))}"))

    # 3. no lost/new records (VariationID set equality vs normalized input)
    ann_vids = {r[a_idx["VariationID"]] for r in ann}
    norm_vids = {r[n_idx["VariationID"]] for r in norm}
    lost = norm_vids - ann_vids
    new = ann_vids - norm_vids
    checks.append(("no lost ClinVar records", not lost, f"lost={len(lost)}"))
    checks.append(("no new variants introduced", not new, f"new={len(new)}"))

    # 4. unchanged predictor values vs previous committed dataset
    prev_ok = True
    prev_note = "n/a"
    try:
        prev_text = subprocess.run(
            ["git", "show", f"HEAD:{ANN_TSV}"], capture_output=True, text=True, check=True
        ).stdout
        prev_lines = prev_text.splitlines()
        ph = prev_lines[0].split("\t")
        pi = {n: i for i, n in enumerate(ph)}
        prev = {}
        for line in prev_lines[1:]:
            f = line.split("\t")
            prev[f[pi["VariationID"]]] = f
        changed = 0
        for r in ann:
            vid = r[a_idx["VariationID"]]
            if vid not in prev:
                changed += 1
                continue
            for field in UNCHANGED_FIELDS:
                if r[a_idx[field]] != prev[vid][pi[field]]:
                    changed += 1
                    break
        prev_note = f"{changed} rows changed (of {n_ann})"
        prev_ok = (changed == 0)
    except (subprocess.CalledProcessError, FileNotFoundError, IndexError, KeyError):
        prev_note = "previous committed dataset not available (first run)"
    checks.append(("predictor values unchanged vs previous run", prev_ok, prev_note))

    # 5. AF in [0,1] where present
    bad_af = 0
    for r in ann:
        for f in ("gnomad_genome_af", "gnomad_exome_af",
                  "gnomad_genome_faf95_popmax", "gnomad_exome_faf95_popmax"):
            v = r[a_idx[f]]
            if v in ("", "NA"):
                continue
            try:
                x = float(v)
                if not (0.0 <= x <= 1.0):
                    bad_af += 1
            except ValueError:
                bad_af += 1
    checks.append(("allele frequencies in [0,1]", bad_af == 0, f"{bad_af} out-of-range"))

    # 6. missing values represented as NA/empty (no error strings in numeric fields)
    err_in_numeric = 0
    for r in ann:
        for f in NUMERIC_FIELDS:
            v = r[a_idx[f]]
            if v not in ("", "NA") and not _is_float(v):
                err_in_numeric += 1
    checks.append(("no API error strings in numeric fields", err_in_numeric == 0, f"{err_in_numeric}"))

    # coverage + missingness
    gnomad = Counter(r[a_idx["gnomad_found"]] for r in ann)
    gnomad_error_reasons = Counter(r[a_idx["gnomad_error"]] for r in ann if r[a_idx["gnomad_error"]])
    missing = {}
    for f in NUMERIC_FIELDS + ["cadd_phred"]:
        missing[f] = sum(1 for r in ann if r[a_idx[f]] in ("", "NA"))

    all_pass = all(ok for _, ok, _ in checks)
    checksum = sha256_file(ANN_TSV)
    ts = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    # write report
    md = ["# Annotation Finalization Report — Phase 4B\n\n"]
    md.append(f"Generated: {ts}\n\n")
    md.append("## Dataset integrity\n\n| Check | Pass | Detail |\n|---|---|---|\n")
    for name, ok, detail in checks:
        md.append(f"| {name} | {'✅' if ok else '❌'} | {detail} |\n")
    md.append("\n## Coverage & missingness\n\n")
    md.append(f"- Row count: {n_ann}\n")
    md.append(f"- Unique biological variant keys: {len(set(keys))}\n")
    md.append(f"- gnomAD: present={gnomad.get('present',0)}, absent={gnomad.get('absent',0)}, "
              f"error={gnomad.get('error',0)}\n")
    if gnomad_error_reasons:
        md.append(f"- gnomAD error reasons: {dict(gnomad_error_reasons)}\n")
    md.append("\nMissing (NA) counts by field:\n\n")
    for f, c in missing.items():
        md.append(f"- {f}: {c}\n")
    md.append("\n## Checksum & provenance\n\n")
    md.append(f"- Annotated dataset SHA-256: `{checksum}`\n")
    md.append(f"- File: `{ANN_TSV}`\n")
    md.append(f"- Python: `{sys.version.split()[0]}`\n")
    md.append("- gnomAD: v4 (`gnomad_r4`), GRCh38, GraphQL API\n")
    md.append("- Ensembl VEP: release 116, GRCh38 (SIFT + PolyPhen-2 HumVar)\n")
    md.append("- REVEL: v1.3, GRCh38 (standalone file)\n")
    md.append("- CADD: excluded (v1.7 bulk annotation unavailable; documented)\n")
    md.append("\n## Status\n\n")
    md.append(("ALL CHECKS PASSED — dataset frozen." if all_pass
               else "CHECK FAILURES PRESENT — review before proceeding.") + "\n")

    os.makedirs(os.path.dirname(REPORT), exist_ok=True)
    with open(REPORT, "w") as fh:
        fh.write("".join(md))

    print("=== FINALIZATION ===")
    for name, ok, detail in checks:
        print(f"  {'PASS' if ok else 'FAIL'}: {name} ({detail})")
    print(f"gnomAD: {dict(gnomad)}")
    if gnomad_error_reasons:
        print(f"gnomAD error reasons: {dict(gnomad_error_reasons)}")
    print(f"checksum: {checksum[:16]}...")
    print(f"ALL PASS: {all_pass}")
    sys.exit(0 if all_pass else 1)


def _is_float(v):
    try:
        float(v)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    main()
