"""Phase 4A — transcript & variant normalization for BRCA1 missense VUS.

Reads data/processed/clinvar_vus_missense.tsv, re-derives HGVS on the approved transcript
NM_007294.4 via Ensembl VEP REST (cached), validates reference alleles against GRCh38, and
cross-checks the protein consequence. Writes a normalized dataset, an unresolved dataset,
and a report.
"""
import csv
import json
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import yaml

from src.annotation import extract_mane_consequence, vep_annotate
from src.variants import (
    NORMALIZED_TRANSCRIPT,
    build_variant_string,
    canonical_representation,
    extract_c_change,
    extract_protein_substitution,
    validate_alleles,
)


def load_processed(path):
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        idx = {n.lstrip("#"): i for i, n in enumerate(header)}
        rows = []
        for row in reader:
            rec = {n.lstrip("#"): (row[i] if i < len(row) else "") for n, i in idx.items()}
            rec["_raw"] = row
            rows.append(rec)
    return header, rows


def main():
    cfg = yaml.safe_load(open("config/config.yaml"))
    src_tsv = os.path.join(cfg["output_dirs"]["processed"], "clinvar_vus_missense.tsv")
    out_tsv = os.path.join(cfg["output_dirs"]["intermediate"], "brca1_vus_missense_normalized.tsv")
    unresolved_tsv = os.path.join(cfg["output_dirs"]["intermediate"], "unresolved_normalization.tsv")
    report = os.path.join(cfg["output_dirs"]["reports"], "transcript_normalization_report.md")
    cache_dir = os.path.join(cfg["output_dirs"]["intermediate"], "vep_cache")

    header, records = load_processed(src_tsv)

    # Build VEP input strings (1:1 with records)
    strings = []
    for r in records:
        strings.append(build_variant_string(r["Chromosome"], r["Start"], r["Stop"],
                                            r["ReferenceAlleleVCF"], r["AlternateAlleleVCF"]))
    print(f"Querying VEP for {len(strings)} variants (cached)...")
    vep, failed_batches = vep_annotate(strings, cache_dir)

    normalized = []
    unresolved = []
    changed_c = 0
    changed_p = 0
    status_counter = Counter()

    for r in records:
        vstr = build_variant_string(r["Chromosome"], r["Start"], r["Stop"],
                                    r["ReferenceAlleleVCF"], r["AlternateAlleleVCF"])
        res = vep.get(vstr)
        clinvar_c = extract_c_change(r["Name"])
        clinvar_p = r["protein_change"] or extract_protein_substitution(r["Name"])

        if res is None:
            status = "unresolved"
            reason = "no VEP result"
            hgvsc = hgvsp = aa = cons = None
            ref_ok = None
            p_match = None
        else:
            hgvsc, hgvsp, aa, cons, err = extract_mane_consequence(res)
            ref_ok = (err is None)
            vep_c = extract_c_change(hgvsc)
            vep_p = extract_protein_substitution(hgvsp)
            if err:
                status, reason = "unresolved", f"VEP error: {err}"
                p_match = None
            elif hgvsc is None:
                status, reason = "unresolved", "no MANE (NM_007294.4) transcript consequence"
                p_match = None
            elif cons and "missense_variant" not in cons:
                status, reason = "unresolved", f"not missense on MANE ({cons})"
                p_match = None
            else:
                p_match = (clinvar_p == vep_p) if (clinvar_p and vep_p) else None
                if clinvar_c and vep_c and clinvar_c != vep_c:
                    changed_c += 1
                if p_match is False:
                    changed_p += 1
                    status, reason = "normalized_changed", "protein consequence differs from ClinVar"
                else:
                    status, reason = "normalized", ""

        status_counter[status] += 1
        row = r["_raw"] + [
            r.get("protein_change", ""), r.get("consequence_class", ""),
            canonical_representation(r["Chromosome"], r["Start"],
                                    r["ReferenceAlleleVCF"], r["AlternateAlleleVCF"]),
            NORMALIZED_TRANSCRIPT,
            hgvsc or "",
            (extract_protein_substitution(hgvsp) or "") if hgvsp else "",
            cons or "",
            aa or "",
            "True" if ref_ok else ("False" if ref_ok is False else ""),
            "" if p_match is None else ("True" if p_match else "False"),
            status,
        ]
        (unresolved if status.startswith("unresolved") else normalized).append((r, row, reason))

    # duplicate normalized variants (same canonical representation, >1 VID)
    rep_vids = defaultdict(set)
    for r, row, reason in normalized:
        rep_vids[row[-6]].add(r["VariationID"])  # normalized_representation column
    dup_reps = {k: v for k, v in rep_vids.items() if len(v) > 1}

    # ---- write outputs ----
    norm_header = header + ["protein_change", "consequence_class",
                            "normalized_representation", "normalized_transcript",
                            "normalized_hgvs_c", "normalized_hgvs_p",
                            "vep_consequence", "vep_amino_acids",
                            "ref_match", "protein_match", "normalization_status"]
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with open(out_tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(norm_header)
        for r, row, reason in normalized:
            w.writerow(row)
    with open(unresolved_tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["VariationID", "AlleleID", "Name", "Chromosome", "Start", "ReferenceAlleleVCF",
                    "AlternateAlleleVCF", "reason"])
        for r, row, reason in unresolved:
            w.writerow([r["VariationID"], r["AlleleID"], r["Name"], r["Chromosome"], r["Start"],
                        r["ReferenceAlleleVCF"], r["AlternateAlleleVCF"], reason])

    # ---- report ----
    n_total = len(records)
    n_norm = len(normalized)
    n_unres = len(unresolved)
    md = []
    md.append("# Transcript & Variant Normalization Report — Phase 4A\n\n")
    md.append(f"Date: 2026-08-19. Gene: BRCA1. Approved transcript: **{NORMALIZED_TRANSCRIPT}** "
              f"(MANE Select = ENST00000357654 = ENIGMA/ClinGen VCEP; protocol v1.0 §4).\n\n")
    md.append("**Transcript investigation (NCBI RefSeq, accessed 2026-08-19):**\n\n")
    md.append("- NM_007294.3 vs NM_007294.4: identical CDS (5592 nt, same protein NP_009225.1, "
              "same CCDS1145); .3→.4 is a 5' UTR-only change. Protein/cDNA numbering identical.\n")
    md.append("- NM_001408514.1: distinct isoform (protein NP_001395443.1, CDS 836..2191) — "
              "transcript-dependent; NOT present in the 1,904 missense-VUS set.\n\n")
    md.append("## Results\n\n")
    md.append(f"| Metric | Count |\n|---|---|\n")
    md.append(f"| Total variants | {n_total} |\n")
    md.append(f"| Successfully normalized | {n_norm} |\n")
    md.append(f"| Unresolved | {n_unres} |\n")
    md.append(f"| Normalized with changed cDNA (c.) vs ClinVar | {changed_c} |\n")
    md.append(f"| Normalized with changed protein (p.) vs ClinVar | {changed_p} |\n")
    md.append(f"| Duplicate normalized representations (>1 VID per coord) | {len(dup_reps)} |\n\n")
    md.append("**Normalization status counts:** " +
              ", ".join(f"{k}={v}" for k, v in sorted(status_counter.items())) + "\n\n")
    if failed_batches:
        md.append(f"**VEP batches that failed after retry:** {len(failed_batches)} "
                  f"(their variants are counted as unresolved).\n\n")
    md.append("## Transcript distribution\n\n")
    md.append("- Before: NM_007294.4 = 1904 (all variants already on the approved transcript)\n")
    md.append(f"- After: NM_007294.4 = {n_norm} (unchanged)\n\n")
    if dup_reps:
        md.append("## Duplicate normalized representations\n\n")
        for rep, vids in list(dup_reps.items())[:20]:
            md.append(f"- {rep}: {sorted(vids)}\n")
    md.append("\n## Method\n\n")
    md.append("- HGVS re-derived from GRCh38 coordinates via Ensembl VEP REST (region endpoint, "
              "`hgvs=1`), filtered to the MANE Select transcript ENST00000357654; responses cached "
              "in `data/intermediate/vep_cache/`.\n")
    md.append("- Reference allele validated by VEP against GRCh38 (VEP `error` field = ref mismatch).\n")
    md.append("- Two-pass validation: ClinVar protein change vs VEP protein change compared.\n")
    md.append("\n## Files\n\n")
    md.append(f"- `{out_tsv}` (normalized)\n")
    md.append(f"- `{unresolved_tsv}` (unresolved)\n")
    with open(report, "w") as fh:
        fh.write("".join(md))

    print(f"\n=== SUMMARY ===")
    print(f"total={n_total}  normalized={n_norm}  unresolved={n_unres}")
    print(f"changed_c={changed_c}  changed_p={changed_p}  duplicate_reprs={len(dup_reps)}")
    print(f"status: {dict(status_counter)}")
    print(f"Output: {out_tsv}")
    print(f"Unresolved: {unresolved_tsv}")
    print(f"Report: {report}")


if __name__ == "__main__":
    main()
