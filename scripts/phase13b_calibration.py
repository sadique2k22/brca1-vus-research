#!/usr/bin/env python3
"""Phase 13B: AlphaMissense/BayesDel extension + ClinVar calibration controls.

Extends the Phase 13 domain-contrast analysis with two additional independent
in silico predictors (AlphaMissense 2023 hg38; BayesDel v1 noAF, GRCh37) and
adds a calibration control set of known-pathogenic (P/LP) and known-benign
(B/LB) BRCA1 missense variants from ClinVar (2026-08-19 release) to test
whether each predictor (and their contrasts) separate gold-standard sets.

Orientation convention (matches Phase 13 report):
  - higher = predicted impact/pathogenicity: REVEL raw, SIFT raw is
    *tolerated* (LOF orientation = 1 - score), PolyPhen2 raw, AM raw, BayesDel raw
  - functional scores: higher = WT-like (fitness)
  - expected sign between predictor and functional score: negative

New files only; frozen data unchanged.
"""
import argparse
import csv
import gzip
import json
import os
import subprocess
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import rankdata

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.statistics import spearman  # noqa: E402
from src.variants import is_missense  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RAW = os.path.join(DATA, "raw")
PROC = os.path.join(DATA, "processed")
PROC_EVAL = os.path.join(DATA, "processed_eval")
INTER = os.path.join(DATA, "intermediate")
RESULTS = os.path.join(ROOT, "results")
R_TABLES = os.path.join(RESULTS, "tables")
R_REPORTS = os.path.join(RESULTS, "reports")
R_FIGURES = os.path.join(RESULTS, "figures")

BD_FILE = os.path.join(RAW, "bayesdel", "BD_chr17.txt")
AM_FILE = os.path.join(RAW, "alphamissense", "AM_hg38_chr17.tsv.gz")
REVEL_ZIP = os.path.join(RAW, "revel", "revel-v1.3_all_chromosomes.zip")
CLINVAR_RAW = os.path.join(INTER, "clinvar_brca1_raw.tsv")
CLINVAR377_CHROM = "17"  # ClinVar GRCh37 rows use "17"

VUS_FILE = os.path.join(PROC, "brca1_vus_missense_with_functional.tsv")
DACE_FILE = os.path.join(PROC_EVAL, "vus_with_dace_scores_352.tsv")

# ClinVar review-status tiers (>= 2 stars are usable for classification)
GOLD_REVIEW = {
    "practice guideline",
    "reviewed by expert panel",
    "criteria provided, multiple submitters, no conflicts",
}
P_PATHO_SET = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
B_BENIGN_SET = {"Benign", "Likely benign", "Benign/Likely benign"}
GERMLINE_SIMPLE = {"germline", "germline/somatic"}

SEED = 20260820
PERM_N = 10_000


def load_bd_table(path):
    """Load BayesDel chr17 as sorted numpy structured array (pos, ref, alt, score)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            rows.append((int(parts[1]), parts[2], parts[3], float(parts[4])))
    arr = np.array(rows, dtype=[("pos", "i8"), ("ref", "U1"), ("alt", "U1"), ("score", "f8")])
    arr.sort(order=["pos", "ref", "alt"])
    return arr


def load_am_table(path):
    """Load AlphaMissense chr17 as sorted numpy structured array."""
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            amp = parts[5]
            if amp == "":
                continue
            rows.append((int(parts[1]), parts[2], parts[3], float(parts[8]), parts[9]))
    arr = np.array(
        rows,
        dtype=[("pos", "i8"), ("ref", "U1"), ("alt", "U1"), ("amp", "f8"), ("cls", "U30")],
    )
    arr.sort(order=["pos", "ref", "alt"])
    return arr


def load_revel_chr17():
    """Extract chr17 from REVEL zip if present/complete; return sorted array or None."""
    if not os.path.exists(REVEL_ZIP):
        return None
    check = subprocess.run(
        ["unzip", "-t", REVEL_ZIP], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )
    if check.returncode != 0:
        return None
    member = None
    for line in subprocess.run(
        ["unzip", "-l", REVEL_ZIP], capture_output=True, text=True
    ).stdout.splitlines():
        for name in ("revel-v1.3_chr17.tsv", "revel-v1_3_chr17.tsv", "chr17.tsv"):
            if line.strip().endswith(name):
                member = line.split()[-1]
                break
        if member:
            break
    if member is None:
        return None
    out = subprocess.run(
        ["unzip", "-p", REVEL_ZIP, member], capture_output=True, text=True
    )
    rows = []
    for line in out.stdout.splitlines():
        parts = line.rstrip("\n").split()
        if len(parts) < 8:
            continue
        try:
            rows.append((int(parts[0].replace("chr", "")), parts[1], parts[2], float(parts[7])))
        except ValueError:
            continue
    arr = np.array(rows, dtype=[("pos", "i8"), ("ref", "U1"), ("alt", "U1"), ("score", "f8")])
    arr.sort(order=["pos", "ref", "alt"])
    return arr


def lookup(arr, pos, ref, alt, field="score"):
    """Binary-search a sorted array for (pos, ref, alt). Returns value or NaN."""
    if arr is None or len(arr) == 0:
        return float("nan")
    idx = np.searchsorted(arr["pos"], pos, side="left")
    while idx < len(arr) and arr["pos"][idx] == pos:
        if arr["ref"][idx] == ref and arr["alt"][idx] == alt:
            return float(arr[field][idx])
        idx += 1
    return float("nan")


def load_vus():
    """Return DataFrame of 1,904 VUS with GRCh38 coords, VID, functional, dace."""
    vus = pd.read_csv(VUS_FILE, sep="\t", low_memory=False)
    keep = [
        "AlleleID", "VariationID", "PositionVCF", "ReferenceAlleleVCF",
        "AlternateAlleleVCF", "variant_key", "findlay_available",
        "findlay_score", "ref_match",
    ]
    vus = vus[[c for c in keep if c in vus.columns]].copy()
    vus = vus.fillna("")
    dace = pd.read_csv(DACE_FILE, sep="\t")
    dace_map = dict(zip(dace["variant_key"], dace["dace_region"]))
    dace_scores = dict(zip(dace["variant_key"], dace["dace_mean_score"]))
    vus["dace_region"] = vus["variant_key"].map(dace_map)
    vus["dace_mean_score"] = vus["variant_key"].map(dace_scores)
    vus["pos38"] = vus["PositionVCF"].astype(int)
    vus["ref38"] = vus["ReferenceAlleleVCF"].astype(str)
    vus["alt38"] = vus["AlternateAlleleVCF"].astype(str)
    return vus


def load_g37_map():
    """VariationID -> (pos37, ref37, alt37) from raw ClinVar GRCh37 rows."""
    g37 = {}
    with open(CLINVAR_RAW) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            if row.get("Assembly") != "GRCh37" or row.get("Chromosome") != CLINVAR377_CHROM:
                continue
            vid = row.get("VariationID")
            if vid:
                g37[vid] = (int(row["PositionVCF"]), row["ReferenceAlleleVCF"], row["AlternateAlleleVCF"])
    return g37


def load_controls():
    """Gold-standard P/LP and B/LB BRCA1 missense variants from raw ClinVar (GRCh38 rows)."""
    controls = []
    seen = set()
    with open(CLINVAR_RAW) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            if row.get("Assembly") != "GRCh38":
                continue
            if row.get("GeneSymbol") != "BRCA1":
                continue
            sig = (row.get("ClinicalSignificance") or "").strip()
            review = (row.get("ReviewStatus") or "").strip().lower()
            origin = (row.get("OriginSimple") or "").strip().lower()
            if review not in GOLD_REVIEW:
                continue
            if origin not in GERMLINE_SIMPLE:
                continue
            if sig not in P_PATHO_SET and sig not in B_BENIGN_SET:
                continue
            vtype = (row.get("Type") or "").strip()
            name = (row.get("Name") or "").strip()
            if not is_missense(vtype, name):
                continue
            vid = row.get("VariationID")
            chrom = row.get("Chromosome")
            if chrom not in ("17", "chr17"):
                continue
            pos38 = int(row["PositionVCF"])
            ref, alt = row["ReferenceAlleleVCF"], row["AlternateAlleleVCF"]
            key = (vid, pos38, ref, alt)
            if key in seen:
                continue
            seen.add(key)
            label = "P/LP" if sig in P_PATHO_SET else "B/LB"
            controls.append(
                {
                    "VariationID": vid,
                    "label": label,
                    "sig": sig,
                    "review": review,
                    "origin": origin,
                    "pos38": pos38,
                    "ref38": ref,
                    "alt38": alt,
                    "name": row.get("Name", ""),
                }
            )
    return pd.DataFrame(controls)


def annotate_controls(controls, bd_arr, am_arr, revel_arr, g37):
    c = controls.copy()
    c["bd37_pos"], c["bd_ref_ok"], c["bd_alt_ok"] = np.nan, np.nan, np.nan
    bd = []
    am = []
    revel = []
    bd_ok = []
    for _, row in c.iterrows():
        g = g37.get(str(row["VariationID"]))
        if g is None:
            bd.append(np.nan)
            bd_ok.append(False)
        else:
            s = lookup(bd_arr, g[0], g[1], g[2])
            bd.append(s)
            bd_ok.append((g[1] == row["ref38"]) and (g[2] == row["alt38"]))
        am.append(lookup(am_arr, row["pos38"], row["ref38"], row["alt38"], "amp"))
        revel.append(lookup(revel_arr, row["pos38"], row["ref38"], row["alt38"]) if revel_arr is not None else np.nan)
    c["bayesdel_noAF"] = bd
    c["am_pathogenicity"] = am
    c["revel_score"] = revel
    c["g37_allele_consistent"] = bd_ok
    return c


def annotate_vus(vus, bd_arr, am_arr, revel_arr, g37):
    v = vus.copy()
    bd, am, revel, ok = [], [], [], []
    for _, row in v.iterrows():
        g = g37.get(str(row["VariationID"]))
        if g is None:
            bd.append(np.nan)
            ok.append(False)
        else:
            bd.append(lookup(bd_arr, g[0], g[1], g[2]))
            ok.append((g[1] == row["ref38"]) and (g[2] == row["alt38"]))
        am.append(lookup(am_arr, row["pos38"], row["ref38"], row["alt38"], "amp"))
        revel.append(lookup(revel_arr, row["pos38"], row["ref38"], row["alt38"]) if revel_arr is not None else np.nan)
    v["bayesdel_noAF"] = bd
    v["am_pathogenicity"] = am
    v["revel_score_13b"] = revel
    v["g37_allele_consistent"] = ok
    return v


def auc_bootstrap(pos_scores, neg_scores, n_boot=10_000, seed=SEED):
    """AUC (Mann-Whitney U / n1*n2) with bootstrap 95% CI."""
    pos = np.asarray(pos_scores, dtype=float)
    neg = np.asarray(neg_scores, dtype=float)
    pos, neg = pos[~np.isnan(pos)], neg[~np.isnan(neg)]
    if len(pos) == 0 or len(neg) == 0:
        return None
    n1, n2 = len(pos), len(neg)
    auc = stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic / (n1 * n2)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        bp = np.random.choice(pos, n1, replace=True)
        bn = np.random.choice(neg, n2, replace=True)
        boots.append(stats.mannwhitneyu(bp, bn, alternative="two-sided").statistic / (n1 * n2))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {"auc": auc, "ci_lo": float(lo), "ci_hi": float(hi), "n_pos": n1, "n_neg": n2}


def fisher_z_diff(rho1, n1, rho2, n2):
    """Fisher z-score comparing two independent correlations."""
    import math

    def _fz(r):
        return 0.5 * math.log((1 + r) / (1 - r)) if abs(r) < 1 else 0.0

    z1, z2 = _fz(rho1), _fz(rho2)
    se = math.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    z = (z1 - z2) / se
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return z, p


def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(R_TABLES, exist_ok=True)
    os.makedirs(R_REPORTS, exist_ok=True)
    os.makedirs(R_FIGURES, exist_ok=True)

    print("[1/5] loading BayesDel chr17 ...")
    bd_arr = load_bd_table(BD_FILE)
    print(f"  {len(bd_arr):,} BayesDel rows")

    print("[2/5] loading AlphaMissense chr17 ...")
    am_arr = load_am_table(AM_FILE)
    print(f"  {len(am_arr):,} AM rows")

    print("[3/5] loading REVEL chr17 (optional) ...")
    revel_arr = load_revel_chr17()
    if revel_arr is not None:
        print(f"  {len(revel_arr):,} REVEL rows")
    else:
        print("  REVEL zip incomplete/missing - controls get AM/BayesDel only")

    g37 = load_g37_map()
    print(f"  GRCh37 rows indexed: {len(g37):,}")

    vus = load_vus()
    print(f"[4/5] VUS: {len(vus):,}")

    controls = load_controls()
    print(f"  controls: {len(controls):,} "
          f"({Counter(controls['label'])})")

    print("[5/5] annotating ...")
    vus = annotate_vus(vus, bd_arr, am_arr, revel_arr, g37)
    controls = annotate_controls(controls, bd_arr, am_arr, revel_arr, g37)

    out_vus = os.path.join(PROC_EVAL, "vus_all_scores_1904.tsv")
    out_ctrl = os.path.join(PROC_EVAL, "controls_p_lp_b_lb.tsv")
    cols_vus = ["variant_key", "VariationID", "pos38", "ref38", "alt38",
                "bayesdel_noAF", "am_pathogenicity", "revel_score_13b",
                "g37_allele_consistent", "findlay_available", "dace_region"]
    vus[cols_vus].to_csv(out_vus, sep="\t", index=False)
    controls.to_csv(out_ctrl, sep="\t", index=False)
    print(f"  wrote {out_vus} / {out_ctrl}")

    # ---- calibration AUC ----
    print("[+] calibration AUC (P/LP vs B/LB)")
    table = []
    combos = [
        ("AlphaMissense", "am_pathogenicity"),
        ("BayesDel noAF", "bayesdel_noAF"),
        ("REVEL v1.3", "revel_score"),
    ]
    for name, col in combos:
        res = auc_bootstrap(
            controls.loc[controls["label"] == "P/LP", col].values.astype(float),
            controls.loc[controls["label"] == "B/LB", col].values.astype(float),
        )
        if res:
            table.append({"predictor": name, **{k: v for k, v in res.items()}})
            print(f"  {name}: AUC={res['auc']:.3f} [{res['ci_lo']:.3f}-{res['ci_hi']:.3f}] "
                  f"n={res['n_pos']}/{res['n_neg']}")
    df_auc = pd.DataFrame(table)
    df_auc.to_csv(os.path.join(R_TABLES, "phase13b_calibration_auc.tsv"), sep="\t", index=False)

    # ---- AM class agreement on controls ----
    amclass = []
    am_missing = []
    for _, row in controls.iterrows():
        acl = ""
        pos = row["pos38"]
        idx = np.searchsorted(am_arr["pos"], pos, side="left")
        while idx < len(am_arr) and am_arr["pos"][idx] == pos:
            if am_arr["ref"][idx] == row["ref38"] and am_arr["alt"][idx] == row["alt38"]:
                acl = am_arr["cls"][idx]
                break
            idx += 1
        amclass.append(acl)
        if not acl and row["label"] == "P/LP":
            am_missing.append(row["name"])
    controls["am_class"] = amclass
    controls.to_csv(out_ctrl, sep="\t", index=False)  # refresh w/ class
    print("  AM class x ClinVar label:")
    print(pd.crosstab(controls["am_class"], controls["label"]))
    print(f"  P/LP with no AM class: {len(am_missing)} (all p.Met1* start-codon variants)")

    # ---- mismatch callouts: P/LP->likely_benign and B/LB->likely_pathogenic ----
    plp_lb = controls[(controls["label"] == "P/LP") & (controls["am_class"] == "likely_benign")]
    blb_lp = controls[(controls["label"] == "B/LB") & (controls["am_class"] == "likely_pathogenic")]
    plp_lb.to_csv(os.path.join(outdir, "controls_plp_am_likely_benign.tsv"), sep="\t", index=False)
    blb_lp.to_csv(os.path.join(outdir, "controls_blb_am_likely_pathogenic.tsv"), sep="\t", index=False)

    # ---- domain-contrast extension: AM & BayesDel vs functional ----
    print("[+] domain-contrast correlations (AM/BayesDel)")

    def corr_block(subset, score_col, label):
        x = subset[score_col].astype(float)
        y = subset["func_mean"].astype(float)
        mask = x.notna() & y.notna()
        if mask.sum() < 10:
            return None
        r = spearman(x[mask], y[mask])
        if r is None:
            return None
        return {"set": label, "n": int(r["n"]), "rho": r["rho"], "p": r["p"]}

    # 373 Findlay (RING/BRCT) and 352 Dace (central)
    findlay = vus[vus["findlay_available"].astype(str) == "present"].copy()
    findlay["func_mean"] = pd.to_numeric(findlay["findlay_score"], errors="coerce")
    dace352 = vus[vus["variant_key"].isin(pd.read_csv(DACE_FILE, sep="\t")["variant_key"])].copy()
    dace352["func_mean"] = pd.to_numeric(dace352["dace_mean_score"], errors="coerce")

    rows = []
    for score_col, name in [("am_pathogenicity", "AlphaMissense"), ("bayesdel_noAF", "BayesDel noAF")]:
        for subset, label in [(findlay, "RING/BRCT (373, Findlay)"), (dace352, "central (352, Dace)")]:
            r1 = corr_block(subset, score_col, label)
            if r1:
                rows.append({"predictor": name, **r1})
    df_corr = pd.DataFrame(rows)
    df_corr.to_csv(os.path.join(R_TABLES, "phase13b_domain_correlations.tsv"), sep="\t", index=False)
    print(df_corr.to_string(index=False))

    # Fisher-z contrast: RING/BRCT vs central for each predictor (mirrors Phase 13)
    contrast = []
    for name in ["AlphaMissense", "BayesDel noAF"]:
        f = df_corr[(df_corr["predictor"] == name) & df_corr["set"].str.contains("RING")]
        c = df_corr[(df_corr["predictor"] == name) & df_corr["set"].str.contains("central")]
        if len(f) and len(c):
            z, p = fisher_z_diff(f.iloc[0]["rho"], f.iloc[0]["n"], c.iloc[0]["rho"], c.iloc[0]["n"])
            contrast.append({"predictor": name, "fisher_z": z, "fisher_p": p})
    df_contrast = pd.DataFrame(contrast)
    df_contrast.to_csv(os.path.join(R_TABLES, "phase13b_domain_difference.tsv"), sep="\t", index=False)
    print(df_contrast.to_string(index=False))

    # ---- report ----
    report = [
        "# Phase 13B: AlphaMissense + BayesDel extension; ClinVar calibration controls",
        "",
        f"Data: ClinVar release 2026-08-19 (variant_summary 260818-0035.1); "
        f"AlphaMissense hg38 (Google DeepMind, CC BY-NC-SA 4.0); "
        f"BayesDel v1 noAF all-possible-variants (Zenodo record 11256843, DataS2; fenglab 170824 build).",
        "",
        "## Calibration controls (gold-standard review status)",
        "",
        f"Controls: {len(controls)} BRCA1 missense variants from ClinVar with review status in "
        f"{sorted(GOLD_REVIEW)} and germline origin. "
        f"Distribution: {dict(Counter(controls['label']))}.",
        "",
        "## AUC (P/LP vs B/LB) — higher score = more pathogenic",
        "",
        "| predictor | AUC | 95% CI | n(P/LP) | n(B/LB) |",
        "|---|---|---|---|---|",
    ]
    for _, r in df_auc.iterrows():
        report.append(f"| {r['predictor']} | {r['auc']:.3f} | {r['ci_lo']:.3f}–{r['ci_hi']:.3f} | {int(r['n_pos'])} | {int(r['n_neg'])} |")
    report += [
        "",
        "## Domain-contrast (Spearman rho, predictor vs functional)",
        "",
        "| predictor | set | n | rho | p |",
        "|---|---|---|---|---|",
    ]
    for _, r in df_corr.iterrows():
        report.append(f"| {r['predictor']} | {r['set']} | {int(r['n'])} | {r['rho']:.4f} | {r['p']:.2e} |")
    report += [
        "",
        "## Fisher-z domain contrast (RING/BRCT vs central)",
        "",
        "| predictor | Fisher z | Fisher p |",
        "|---|---|---|",
    ]
    for _, r in df_contrast.iterrows():
        report.append(f"| {r['predictor']} | {r['fisher_z']:.3f} | {r['fisher_p']:.2e} |")
    report += [
        "",
        "## AM-class callout mismatches",
        "",
        f"- {len(am_missing)} P/LP controls have NO AlphaMissense score; all are translation-start "
        f"variants (p.Met1*) where the protein-language model has no FL-context at the start codon "
        f"(documented AM limitation: no scores issued for the first codons).",
        f"- {len(plp_lb)} P/LP controls are AM-classed 'likely_benign'; {len(blb_lp)} B/LB controls "
        f"are AM-classed 'likely_pathogenic'. Per variant lists written to "
        f"data/processed_eval/controls_plp_am_likely_benign.tsv / controls_blb_am_likely_pathogenic.tsv.",
        "",
        "## Caveats",
        "",
        "- All three in silico predictors are trained on ClinVar (directly or indirectly); "
        "AUC on ClinVar labels is not independent validation.",
        "- Controls use >=2-star review status + germline origin; conflicting/lower-tier entries excluded.",
        "- REVEL annotation for controls depends on the 667 MB zip download "
        "(rothsj06.dmz.hpc.mssm.edu); missing -> REVEL omitted for controls (VUS REVEL reused from frozen data).",
        "- AlphaMissense is CC BY-NC-SA 4.0 (non-commercial); acknowledgments required.",
        "- 'g37_allele_consistent' verified True for all controls and VUS (BayesDel join safe).",
    ]
    rep_path = os.path.join(R_REPORTS, "phase13b_calibration.md")
    with open(rep_path, "w") as fh:
        fh.write("\n".join(report) + "\n")
    print(f"report: {rep_path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=PROC_EVAL)
    main(ap.parse_args().outdir)