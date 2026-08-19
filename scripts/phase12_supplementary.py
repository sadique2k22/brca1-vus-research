"""Phase 12 — generate supplementary tables (S1-S4) and the figure package.

Derived outputs only; no values are modified. Reads frozen datasets read-only.
"""
import csv
import os
import shutil

SUP = "manuscript/supplementary"
FIGDIR = "manuscript/figures"


def load(path):
    with open(path, newline="", encoding="utf-8") as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        idx = {n.lstrip("#"): i for i, n in enumerate(header)}
        rows = [row for row in r]
    return header, idx, rows


def write_tsv(path, cols, idx, rows, extra_header=None):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for r in rows:
            w.writerow([r[idx[c]] if c in idx and idx[c] < len(r) else "" for c in cols])


def main():
    os.makedirs(SUP, exist_ok=True)

    # ---- S1: all 1,904 variants ----
    header, idx, rows = load("data/processed/brca1_vus_missense_annotated.tsv")
    s1_cols = ["VariationID", "AlleleID", "variant_key", "Chromosome", "Start",
               "ReferenceAlleleVCF", "AlternateAlleleVCF", "normalized_hgvs_c",
               "normalized_hgvs_p", "protein_change", "ClinicalSignificance",
               "normalized_transcript", "gnomad_found", "gnomad_exome_af",
               "gnomad_genome_af", "gnomad_exome_faf95_popmax", "revel_score",
               "sift_score", "polyphen_score"]
    write_tsv(os.path.join(SUP, "S1_all_1904_variants.tsv"), s1_cols, idx, rows)
    # xlsx for S1
    try:
        from openpyxl import Workbook
        wb = Workbook(); ws = wb.active; ws.title = "S1"
        ws.append(s1_cols)
        for r in rows:
            ws.append([r[idx[c]] if c in idx and idx[c] < len(r) else "" for c in s1_cols])
        wb.save(os.path.join(SUP, "S1_all_1904_variants.xlsx"))
    except ImportError:
        pass

    # ---- S2: 373 variants with Findlay measurements ----
    _, idx2, rows2 = load("data/processed/brca1_vus_missense_with_functional.tsv")
    s2_cols = ["VariationID", "variant_key", "protein_change", "normalized_hgvs_c",
               "findlay_score", "findlay_score_rep1", "findlay_score_rep2",
               "findlay_score_rna", "revel_score", "sift_score", "polyphen_score"]
    s2_rows = [r for r in rows2 if r[idx2["findlay_available"]] == "present"]
    with open(os.path.join(SUP, "S2_findlay_373.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(s2_cols)
        for r in s2_rows:
            w.writerow([r[idx2[c]] if c in idx2 and idx2[c] < len(r) else "" for c in s2_cols])

    # ---- S3: 41-variant evidence matrix ----
    with open("results/tables/phase9_evidence_matrix.tsv", newline="", encoding="utf-8") as fh:
        data = fh.read()
    with open(os.path.join(SUP, "S3_evidence_matrix_41.tsv"), "w", encoding="utf-8") as fh:
        fh.write(data)

    # ---- S4: literature verification log ----
    with open("results/tables/phase9_literature_verification.tsv", newline="", encoding="utf-8") as fh:
        data = fh.read()
    with open(os.path.join(SUP, "S4_literature_verification.tsv"), "w", encoding="utf-8") as fh:
        fh.write(data)

    # ---- figure package ----
    os.makedirs(FIGDIR, exist_ok=True)
    manifest = [["figure", "filename", "description", "source_script", "source_data", "format", "resolution"]]
    src = "results/figures"
    if os.path.isdir(src):
        for f in sorted(os.listdir(src)):
            if f.endswith((".svg", ".png")):
                shutil.copy(os.path.join(src, f), os.path.join(FIGDIR, f))
    # manifest from figure_mapping (descriptive)
    for num, name in [
        ("Fig 1", "fig1_workflow"), ("Fig 2", "fig2_revel_distribution"),
        ("Fig 3", "fig3_sift_distribution"), ("Fig 4", "fig4_polyphen_distribution"),
        ("Fig 5", "fig5_agreement_matrix"), ("Fig 6", "fig8_revel_by_presence"),
        ("Fig 7", "fig9_log10_af_vs_revel"), ("Fig 8", "fig11_cohort_composition"),
        ("Fig 9", "fig12_predictor_vs_functional"), ("Fig 10", "fig13_concordance"),
        ("Fig 11", "fig14_functional_by_category"), ("Fig 12", "fig16_domain_comparison"),
    ]:
        manifest.append([num, name + ".svg", "", "src/figures.py", "annotated/with-functional", "SVG+PNG", "300 DPI"])
    with open("manuscript/figure_manifest.tsv", "w", newline="", encoding="utf-8") as fh:
        csv.writer(fh, delimiter="\t").writerows(manifest)

    print(f"S1 rows: {len(rows)}")
    print(f"S2 rows: {len(s2_rows)}")
    print(f"S3/S4 copied; figures copied: {len(os.listdir(FIGDIR))} files")


if __name__ == "__main__":
    main()
