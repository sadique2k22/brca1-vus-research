"""Phase 3 — ClinVar data audit and VUS/missense filtering.

Reads data/intermediate/clinvar_brca1_raw.tsv, audits the dataset (duplicates,
submissions, conditions, conflicts, coordinates, transcripts, builds, types), applies the
FINAL inclusion criteria from protocol.md v1.0 via src/variants.filter_vus_missense, and
writes:
  - data/processed/clinvar_vus_missense.tsv
  - results/reports/clinvar_audit.md
"""
import csv
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import yaml

from src.variants import classify_protein_change, filter_vus_missense

TRANSCRIPT_RE = re.compile(r'(N[MR]_\d+\.\d+)')


def load_records(tsv_path):
    header = None
    records = []
    idx = {}
    with open(tsv_path, newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if header is None:
                header = row
                idx = {name.lstrip("#"): i for i, name in enumerate(row)}
                continue
            rec = {name.lstrip("#"): (row[i] if i < len(row) else "") for name, i in idx.items()}
            rec["_raw"] = row
            records.append(rec)
    return header, records


def main():
    cfg = yaml.safe_load(open("config/config.yaml"))
    gene = cfg["study"]["gene"]
    int_tsv = os.path.join(cfg["output_dirs"]["intermediate"], f"clinvar_{gene.lower()}_raw.tsv")
    out_tsv = os.path.join(cfg["output_dirs"]["processed"], "clinvar_vus_missense.tsv")
    out_report = os.path.join(cfg["output_dirs"]["reports"], "clinvar_audit.md")

    header, records = load_records(int_tsv)

    # ---- Audit investigations ----
    n_rows = len(records)
    vids = Counter(r["VariationID"] for r in records if r["VariationID"])
    empty_vid = sum(1 for r in records if not r["VariationID"])

    assembly_by_vid = defaultdict(set)
    for r in records:
        if r["VariationID"]:
            assembly_by_vid[r["VariationID"]].add(r["Assembly"])

    n_grch38_only = sum(1 for s in assembly_by_vid.values() if s == {"GRCh38"})
    n_grch37_only = sum(1 for s in assembly_by_vid.values() if s == {"GRCh37"})
    n_both = sum(1 for s in assembly_by_vid.values() if "GRCh38" in s and "GRCh37" in s)
    n_na_involved = sum(1 for s in assembly_by_vid.values() if "na" in s)

    # submissions
    nsub = Counter(r["NumberSubmitters"] for r in records if r["Assembly"] == "GRCh38")
    multi_submit = sum(c for k, c in nsub.items() if k.isdigit() and int(k) > 1)

    # conditions (GRCh38 rows)
    multi_cond = 0
    for r in records:
        if r["Assembly"] == "GRCh38" and "|" in (r.get("PhenotypeList") or ""):
            multi_cond += 1

    # transcripts
    transcripts = Counter()
    no_transcript = 0
    for r in records:
        if r["Assembly"] == "GRCh38":
            m = TRANSCRIPT_RE.search(r["Name"] or "")
            if m:
                transcripts[m.group(1)] += 1
            else:
                no_transcript += 1

    # multiple representations (GRCh38 coordinate groups with >1 VID)
    coord_vids = defaultdict(set)
    for r in records:
        if r["Assembly"] == "GRCh38" and r["VariationID"]:
            key = (r["Chromosome"], r["PositionVCF"], r["ReferenceAlleleVCF"], r["AlternateAlleleVCF"])
            coord_vids[key].add(r["VariationID"])
    multi_repr_groups = sum(1 for v in coord_vids.values() if len(v) > 1)

    # types
    types = Counter(r["Type"] for r in records if r["Assembly"] == "GRCh38")

    # missing coords / protein / hgvs (GRCh38 rows)
    grch38_rows = [r for r in records if r["Assembly"] == "GRCh38"]
    missing_coord = sum(1 for r in grch38_rows if not r["Chromosome"] or not r["Start"])
    missing_protein = sum(1 for r in grch38_rows if "p." not in (r["Name"] or ""))
    missing_hgvs = sum(1 for r in grch38_rows if ":c." not in (r["Name"] or ""))

    # conflicting significance (aggregate, unique VIDs)
    sig = Counter(r["ClinicalSignificance"] for r in grch38_rows)

    # ---- Apply FINAL inclusion criteria ----
    kept, steps = filter_vus_missense(records)

    # ---- Write processed TSV (original 43 cols + derived) ----
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with open(out_tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header + ["protein_change", "consequence_class"])
        for r in kept:
            name = r.get("Name", "")
            w.writerow(r["_raw"] + [extract_pc(name), classify_protein_change(name)])

    # ---- Write report ----
    md = []
    md.append("# ClinVar Data Audit — Phase 3\n")
    md.append(f"Date: 2026-08-19. Gene: {gene}. Source: {int_tsv} (derived from "
              f"`variant_summary.txt.gz`, build {cfg['clinvar']['build']}).\n")
    md.append("## Audit findings\n")
    md.append("| # | Investigation | Finding |\n|---|---|---|\n")
    md.append(f"| 1 | Duplicate variants (per-assembly rows) | {n_rows} rows → {len(vids)} unique Variation IDs |\n")
    md.append(f"| 2 | Duplicate Variation IDs | each VID has GRCh37 and/or GRCh38 rows; "
              f"{n_both} VIDs have both, {n_grch37_only} GRCh37-only, {n_grch38_only} GRCh38-only |\n")
    md.append(f"| 3 | Multiple submissions | {multi_submit} GRCh38 variants have >1 submitter "
              f"(NumberSubmitters distribution: {dict(sorted(nsub.items()))}) |\n")
    md.append(f"| 4 | Multiple conditions | {multi_cond} GRCh38 variants have >1 condition (pipe-delimited PhenotypeList) |\n")
    md.append(f"| 5 | Conflicting significance | {sig.get('Conflicting classifications of pathogenicity', 0)} GRCh38 variants "
              f"have aggregate 'Conflicting classifications of pathogenicity' |\n")
    md.append(f"| 6 | Review-status differences | see distribution below; "
              f"{Counter(r['ReviewStatus'] for r in grch38_rows).get('criteria provided, conflicting classifications', 0)} 'conflicting classifications' |\n")
    md.append(f"| 7 | Missing genomic coordinates | {missing_coord} GRCh38 rows missing chr/start; "
              f"{n_na_involved} VIDs have an Assembly='na' row |\n")
    md.append(f"| 8 | Missing transcript info | {no_transcript} GRCh38 rows have no NM_/NR_ accession in Name |\n")
    md.append(f"| 9 | Missing protein change | {missing_protein} GRCh38 rows have no p. notation |\n")
    md.append(f"| 10 | Multiple transcripts | transcript accessions present: {dict(sorted(transcripts.items()))} |\n")
    md.append(f"| 11 | Multiple representations (same coord, >1 VID) | {multi_repr_groups} GRCh38 coordinate groups map to >1 VID |\n")
    md.append(f"| 12 | Genome-build inconsistencies | per-variant: both={n_both}, GRCh37-only={n_grch37_only}, "
              f"GRCh38-only={n_grch38_only}, 'na' involved={n_na_involved} |\n")
    md.append(f"| 13 | Unexpected variant types | {dict(sorted(types.items()))} |\n")
    if empty_vid:
        md.append(f"\n**Note:** {empty_vid} rows have an empty VariationID (flagged).\n")

    md.append("\n## Clinical significance (unique GRCh38 rows)\n\n")
    for k, v in sig.most_common():
        md.append(f"- {k}: {v}\n")

    md.append("\n## Review status (unique GRCh38 rows)\n\n")
    for k, v in Counter(r["ReviewStatus"] for r in grch38_rows).most_common():
        md.append(f"- {k}: {v}\n")

    md.append("\n## Filtering steps (FINAL inclusion criteria)\n\n")
    md.append("| Step | Before | Removed | Remaining | Reason |\n|---|---|---|---|---|\n")
    for name, before, removed, remaining, reason in steps:
        md.append(f"| {name} | {before} | {removed} | {remaining} | {reason} |\n")

    md.append(f"\n## Result\n\n")
    md.append(f"**Final candidate VUS (missense, GRCh38): {len(kept)} variants.**\n\n")
    md.append("Written to `data/processed/clinvar_vus_missense.tsv` (original 43 columns + "
              "`protein_change` + `consequence_class`; traceable via `VariationID`/`#AlleleID`).\n\n")

    md.append("## Limitations of this classification\n\n")
    md.append("- Missense is inferred from ClinVar p. notation + `Type='single nucleotide variant'`; "
              "variant_summary.txt has **no molecular-consequence (MCNS) column**, so exonic variants that are "
              "missense *and* splice-affecting are not distinguished here (deferred to Phase 4 VEP).\n")
    md.append("- Splice-donor/acceptor SNVs (e.g. c.4096+1G>A) have no p. notation and are correctly excluded.\n")
    md.append("- Transcript normalization to NM_007294.4 (VCEP) is deferred to Phase 4; transcripts are reported above.\n")

    with open(out_report, "w") as fh:
        fh.write("".join(md))

    # concise stdout summary
    print("=== FILTERING STEPS ===")
    for name, before, removed, remaining, reason in steps:
        print(f"  {name}: {before} -> {remaining} (removed {removed}: {reason})")
    print(f"\nFINAL candidate VUS (missense, GRCh38): {len(kept)}")
    print(f"Output: {out_tsv}")
    print(f"Report: {out_report}")


def extract_pc(name):
    m = re.search(r'p\.([^)\s]+)', name or "")
    return m.group(1) if m else ""


if __name__ == "__main__":
    main()
