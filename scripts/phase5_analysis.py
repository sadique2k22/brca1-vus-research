"""Phase 5 — exploratory descriptive statistics, predictor agreement, population analysis.

Reads the FROZEN annotated dataset (read-only), verifies its checksum, computes
descriptive statistics and agreement/correlation measures, writes tables + figures +
a report. No pathogenic/benign classification. No ACMG codes.
"""
import csv
import hashlib
import math
import os
import re
import sys
import time
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import yaml

from src import figures
from src import statistics as S

ANN_TSV = "data/processed/brca1_vus_missense_annotated.tsv"
FIN_REPORT = "results/reports/annotation_finalization_report.md"
TABLES = "results/tables"
REPORTS = "results/reports"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def load_ann():
    with open(ANN_TSV, newline="", encoding="utf-8") as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        idx = {n: i for i, n in enumerate(header)}
        rows = [row for row in r]
    return header, idx, rows


def write_tsv(name, header, rows):
    os.makedirs(TABLES, exist_ok=True)
    with open(os.path.join(TABLES, name), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)


def main():
    t0 = time.time()
    # ---- STEP 0: verify frozen dataset ----
    m = re.search(r"`([0-9a-f]{64})`", open(FIN_REPORT).read())
    expected = m.group(1) if m else None
    actual = sha256(ANN_TSV)
    if expected and actual != expected:
        print(f"CHECKSUM MISMATCH: expected {expected} got {actual} — STOP")
        sys.exit(1)

    header, idx, rows = load_ann()
    n = len(rows)
    keys = [r[idx["variant_key"]] for r in rows]
    assert n == 1904 and len(set(keys)) == 1904, f"integrity fail: rows={n} keys={len(set(keys))}"
    print(f"STEP 0 OK: {n} rows, checksum {actual[:16]}...")

    # ---- helpers to extract numeric/categorical fields ----
    def num(r, f):
        return S.parse_float(r[idx[f]])

    revel = [num(r, "revel_score") for r in rows]
    sift = [num(r, "sift_score") for r in rows]
    polyphen = [num(r, "polyphen_score") for r in rows]
    gnomad_found = [r[idx["gnomad_found"]] for r in rows]

    def af_global(r):
        return num(r, "gnomad_exome_af") if num(r, "gnomad_exome_af") is not None else num(r, "gnomad_genome_af")

    def faf_popmax(r):
        return num(r, "gnomad_exome_faf95_popmax") if num(r, "gnomad_exome_faf95_popmax") is not None else num(r, "gnomad_genome_faf95_popmax")

    present = [r for r in rows if r[idx["gnomad_found"]] == "present"]
    absent = [r for r in rows if r[idx["gnomad_found"]] == "absent"]
    af_present = [af_global(r) for r in present]
    af_present_vals = [a for a in af_present if a is not None]

    # ---- categories ----
    revel_cat = [S.categorize_revel(s) for s in revel]
    sift_cat = [S.categorize_sift(s) for s in sift]
    poly_cat = [S.categorize_polyphen(s) for s in polyphen]
    r2 = [S.impact_2class("revel", s) for s in revel]
    s2 = [S.impact_2class("sift", s) for s in sift]
    p2 = [S.impact_2class("polyphen", s) for s in polyphen]

    # ---- STEP 2: dataset characterization ----
    ds = {
        "total": n,
        "present": len(present), "absent": len(absent),
        "error": sum(1 for f in gnomad_found if f == "error"),
        "pct_present": round(100 * len(present) / n, 1),
        "pct_absent": round(100 * len(absent) / n, 1),
    }

    # ---- STEP 3: predictor distributions ----
    pred_summary = {}
    for name, vals in (("revel", revel), ("sift", sift), ("polyphen", polyphen)):
        pred_summary[name] = S.descriptive(vals)
    cat_counts = {
        "revel": dict(Counter(c for c in revel_cat if c)),
        "sift": dict(Counter(c for c in sift_cat if c)),
        "polyphen": dict(Counter(c for c in poly_cat if c)),
    }

    # ---- STEP 4: agreement ----
    combos = Counter()
    for rc, sc, pc in zip(revel_cat, sift_cat, poly_cat):
        if rc and sc and pc:
            combos[(rc, sc, pc)] += 1
    combo_rows = [["revel", "sift", "polyphen", "count"]]
    for (rc, sc, pc), c in sorted(combos.items(), key=lambda x: -x[1]):
        combo_rows.append([rc, sc, pc, c])

    agreement = {}
    for label, a, b in (("revel_vs_sift", r2, s2), ("revel_vs_polyphen", r2, p2), ("sift_vs_polyphen", s2, p2)):
        agreement[label] = {
            "percent_agreement": S.percent_agreement(a, b),
            "kappa": S.cohens_kappa(a, b),
        }

    # ---- STEP 5: correlations (align SIFT direction) ----
    sift_damaging = [(1 - s) if s is not None else None for s in sift]
    correlations = {
        "revel_vs_polyphen": S.spearman(revel, polyphen),
        "revel_vs_sift": S.spearman(revel, sift_damaging),
        "sift_vs_polyphen": S.spearman(sift_damaging, polyphen),
    }

    # ---- STEP 6: gnomAD present vs absent ----
    mw = {}
    for name, pred in (("revel", revel), ("sift", sift_damaging), ("polyphen", polyphen)):
        mw[name] = S.mannwhitney(
            [pred[i] for i in range(n) if gnomad_found[i] == "present"],
            [pred[i] for i in range(n) if gnomad_found[i] == "absent"],
        )
    log10af = [math.log10(a) if a and a > 0 else None for a in af_present]
    af_corr = {
        "revel": S.spearman(log10af, [num(r, "revel_score") for r in present]),
        "sift": S.spearman(log10af, [(1 - num(r, "sift_score")) if num(r, "sift_score") is not None else None for r in present]),
        "polyphen": S.spearman(log10af, [num(r, "polyphen_score") for r in present]),
    }

    # ---- STEP 7: population frequency ----
    pop_outliers = []
    for r in present:
        fa = faf_popmax(r)
        ga = af_global(r)
        if fa is not None and fa >= 0.001:
            pop_outliers.append({
                "variant_key": r[idx["variant_key"]],
                "VariationID": r[idx["VariationID"]],
                "protein_change": r[idx["protein_change"]],
                "global_af": ga,
                "popmax_faf95": fa,
                "population": r[idx["gnomad_exome_faf95_pop"]] or r[idx["gnomad_genome_faf95_pop"]],
                "ratio": (fa / ga) if (ga and ga > 0) else None,
            })

    # ---- STEP 8: pattern classes ----
    patterns = defaultdict(list)
    for i, r in enumerate(rows):
        vid = r[idx["VariationID"]]
        pc = r[idx["protein_change"]]
        found = gnomad_found[i]
        rv = revel[i]
        sv = sift[i]
        pv = polyphen[i]
        fa = faf_popmax(r)
        if found == "present" and rv is not None and rv >= 0.644:
            patterns["A"].append((vid, pc))
        if found == "present" and rv is not None and rv <= 0.290:
            patterns["B"].append((vid, pc))
        if found == "absent" and rv is not None and rv >= 0.932:
            patterns["C"].append((vid, pc))
        if found == "absent" and rv is not None and rv <= 0.290:
            patterns["D"].append((vid, pc))
        if rv is not None and sv is not None and pv is not None:
            if rv >= 0.644 and sv > 0.05 and pv <= 0.446:
                patterns["E"].append((vid, pc))
            elif rv <= 0.290 and sv <= 0.05 and pv >= 0.446:
                patterns["E"].append((vid, pc))
        if fa is not None and fa >= 0.001:
            patterns["F"].append((vid, pc))
        if rv is not None and (rv >= 0.932 or rv <= 0.003):
            patterns["G"].append((vid, pc))

    # ---- write tables ----
    write_tsv("dataset_summary.tsv", ["metric", "value"], [[k, v] for k, v in ds.items()])
    pr = [["predictor", "n", "min", "max", "mean", "median", "std", "q1", "q3"]]
    for name in ("revel", "sift", "polyphen"):
        d = pred_summary[name]
        pr.append([name, d.get("n", 0), d.get("min"), d.get("max"), d.get("mean"),
                   d.get("median"), d.get("std"), d.get("q1"), d.get("q3")])
    write_tsv("predictor_summary.tsv", pr[0], pr[1:])
    write_tsv("predictor_agreement.tsv",
              ["comparison", "percent_agreement", "kappa", "p", "n"],
              [[k, v["percent_agreement"],
                v["kappa"]["kappa"] if v["kappa"] else None,
                v["kappa"]["p"] if v["kappa"] else None,
                v["kappa"]["n"] if v["kappa"] else None] for k, v in agreement.items()])
    write_tsv("predictor_correlations.tsv",
              ["comparison", "spearman_rho", "p", "n"],
              [[k, v["rho"] if v else None, v["p"] if v else None, v["n"] if v else None]
               for k, v in correlations.items()])
    gs = [["group", "n", "min", "max", "mean", "median", "std", "q1", "q3"]]
    for grp, cond in (("present", "present"), ("absent", "absent")):
        d = S.descriptive([revel[i] for i in range(n) if gnomad_found[i] == cond])
        gs.append([grp + "_revel", d.get("n", 0), d.get("min"), d.get("max"), d.get("mean"),
                   d.get("median"), d.get("std"), d.get("q1"), d.get("q3")])
    write_tsv("gnomad_summary.tsv", gs[0], gs[1:])
    if pop_outliers:
        write_tsv("population_frequency_outliers.tsv",
                  ["variant_key", "VariationID", "protein_change", "global_af", "popmax_faf95", "population", "ratio"],
                  [[p["variant_key"], p["VariationID"], p["protein_change"], p["global_af"],
                    p["popmax_faf95"], p["population"], p["ratio"]] for p in pop_outliers])
    pc_rows = [["pattern", "description", "count", "examples"]]
    desc = {"A": "gnomAD-present + high impact", "B": "gnomAD-present + tolerance",
            "C": "gnomAD-absent + strong impact", "D": "gnomAD-absent + tolerance",
            "E": "strong predictor disagreement", "F": "elevated population frequency",
            "G": "extreme REVEL score"}
    for k in "ABCDEFG":
        pc_rows.append([k, desc[k], len(patterns[k]),
                        "; ".join(f"{v}({p})" for v, p in patterns[k][:5])])
    write_tsv("pattern_candidates.tsv", pc_rows[0], pc_rows[1:])

    # ---- figures ----
    figures.fig1_workflow()
    figures.fig2_revel(revel)
    figures.fig3_sift(sift)
    figures.fig4_polyphen(polyphen)
    # agreement matrix: REVEL(3) x PolyPhen(2)
    rv3 = ["impact", "intermediate", "tolerance"]
    pp3 = ["benign", "possibly_damaging", "probably_damaging"]
    mat = [[0] * 3 for _ in range(3)]
    for rc, pc in zip(revel_cat, poly_cat):
        if rc in rv3 and pc in pp3:
            mat[rv3.index(rc)][pp3.index(pc)] += 1
    figures.fig5_agreement(np_array(mat), rv3, pp3, "Figure 5 — REVEL x PolyPhen-2")
    figures.fig6_revel_polyphen(revel, polyphen)
    figures.fig7_revel_sift(revel, sift)
    for name, pred, ylab in (("REVEL", revel, "REVEL"), ("SIFT", sift, "SIFT (lower=damaging)"),
                             ("PolyPhen-2", polyphen, "PolyPhen-2")):
        figures.fig8_present_absent([pred[i] for i in range(n) if gnomad_found[i] == "present"],
                                    [pred[i] for i in range(n) if gnomad_found[i] == "absent"],
                                    name, ylab)
    figures.fig9_af_revel(log10af, [num(r, "revel_score") for r in present])
    figures.fig10_population_freq(af_present_vals, [faf_popmax(r) for r in present])

    # ---- report ----
    write_report(ds, pred_summary, cat_counts, combos, agreement, correlations,
                 mw, af_corr, pop_outliers, patterns, log10af, af_present_vals,
                 actual, time.time() - t0)
    print("DONE in %.1fs" % (time.time() - t0))
    for k in "ABCDEFG":
        print(f"  pattern {k} ({desc[k]}): {len(patterns[k])}")


def np_array(m):
    import numpy as np
    return np.array(m, dtype=float)


def write_report(ds, pred_summary, cat_counts, combos, agreement, correlations,
                 mw, af_corr, pop_outliers, patterns, log10af, af_vals, checksum, runtime):
    os.makedirs(REPORTS, exist_ok=True)
    L = []
    L.append("# Phase 5 — Statistical Analysis Report\n\n")
    L.append("**Exploratory / descriptive only.** No pathogenic/benign classification, "
             "no ACMG codes, no clinical interpretation.\n\n")
    L.append(f"Frozen dataset checksum verified: `{checksum[:16]}…`. Runtime {runtime:.1f}s.\n\n")

    L.append("## 1. Dataset summary\n\n")
    L.append(f"- Total missense VUS: **{ds['total']}**\n")
    L.append(f"- gnomAD present: {ds['present']} ({ds['pct_present']}%)\n")
    L.append(f"- gnomAD absent: {ds['absent']} ({ds['pct_absent']}%)\n")
    if af_vals:
        d = S.descriptive(af_vals)
        L.append(f"- gnomAD-present global AF: median={d['median']:.3g}, mean={d['mean']:.3g}, "
                 f"min={d['min']:.3g}, max={d['max']:.3g}, q1={d['q1']:.3g}, q3={d['q3']:.3g}\n")
    L.append("\n")

    L.append("## 2. Predictor distributions\n\n| Predictor | n | min | max | mean | median | q1 | q3 |\n|---|---|---|---|---|---|---|---|\n")
    for name in ("revel", "sift", "polyphen"):
        d = pred_summary[name]
        L.append(f"| {name} | {d.get('n',0)} | {d.get('min'):.3f} | {d.get('max'):.3f} | "
                 f"{d.get('mean'):.3f} | {d.get('median'):.3f} | {d.get('q1'):.3f} | {d.get('q3'):.3f} |\n")
    L.append("\nCategory counts (neutral terminology):\n\n")
    for name, cc in cat_counts.items():
        L.append(f"- {name}: {cc}\n")
    L.append("\n")

    L.append("## 3. Predictor agreement\n\n")
    L.append("Pairwise categorical agreement (2-class impact/tolerance; REVEL-intermediate excluded):\n\n")
    L.append("| Comparison | % agreement | Cohen's kappa | p | n |\n|---|---|---|---|---|\n")
    for k, v in agreement.items():
        kap = v["kappa"]
        L.append(f"| {k} | {100*v['percent_agreement']:.1f}% | "
                 f"{kap['kappa']:.3f} | {kap['p']:.3g} | {kap['n']} |\n")
    L.append("\nFull REVEL x SIFT x PolyPhen combination table (count):\n\n")
    L.append("| REVEL | SIFT | PolyPhen | count |\n|---|---|---|---|\n")
    for (rc, sc, pc), c in sorted(combos.items(), key=lambda x: -x[1]):
        L.append(f"| {rc} | {sc} | {pc} | {c} |\n")
    L.append("\n")

    L.append("## 4. Predictor correlations (Spearman)\n\n")
    L.append("(SIFT direction inverted to 'higher = more damaging' for comparability.)\n\n")
    L.append("| Comparison | rho | p | n |\n|---|---|---|---|\n")
    for k, v in correlations.items():
        if v:
            L.append(f"| {k} | {v['rho']:.3f} | {v['p']:.3g} | {v['n']} |\n")
    L.append("\n")

    L.append("## 5. gnomAD present vs absent (predictor scores)\n\n")
    L.append("Mann-Whitney U (exploratory):\n\n| Predictor | U | p | n_present | n_absent |\n|---|---|---|---|---|\n")
    for name, v in mw.items():
        if v:
            L.append(f"| {name} | {v['u']:.0f} | {v['p']:.3g} | {v['n1']} | {v['n2']} |\n")
    L.append("\nlog10(global AF) vs predictor (Spearman, gnomAD-present only):\n\n")
    L.append("| Predictor | rho | p | n |\n|---|---|---|---|\n")
    for k, v in af_corr.items():
        if v:
            L.append(f"| {k} | {v['rho']:.3f} | {v['p']:.3g} | {v['n']} |\n")
    L.append("\n")

    L.append("## 6. Population-frequency findings\n\n")
    L.append(f"Variants with elevated population filtering AF (faf95 popmax >= 0.001): **{len(pop_outliers)}**\n\n")
    if pop_outliers:
        L.append("| Variant | global AF | popmax faf95 | population |\n|---|---|---|---|\n")
        for p in sorted(pop_outliers, key=lambda x: -(x["popmax_faf95"] or 0))[:20]:
            L.append(f"| {p['protein_change']} ({p['VariationID']}) | {p['global_af']:.3g} | "
                     f"{p['popmax_faf95']:.3g} | {p['population']} |\n")
    L.append("\n(Full list in `results/tables/population_frequency_outliers.tsv`.)\n\n")

    L.append("## 7. Pattern classes (Phase 6 literature candidates)\n\n")
    L.append("| Class | Description | Count |\n|---|---|---|\n")
    desc = {"A": "gnomAD-present + high impact", "B": "gnomAD-present + tolerance",
            "C": "gnomAD-absent + strong impact", "D": "gnomAD-absent + tolerance",
            "E": "strong predictor disagreement", "F": "elevated population frequency",
            "G": "extreme REVEL score"}
    for k in "ABCDEFG":
        L.append(f"| {k} | {desc[k]} | {len(patterns[k])} |\n")
    L.append("\n`results/tables/pattern_candidates.tsv` — NOT a pathogenicity list.\n\n")

    L.append("## 8. Statistical limitations\n\n")
    L.append("- Exploratory; multiple comparisons were not corrected (Mann-Whitney x3, "
             "correlations x6). p-values are descriptive, not confirmatory.\n")
    L.append("- REVEL/SIFT/PolyPhen are correlated (REVEL is an ensemble) — treated as "
             "correlated predictors, not independent evidence lines.\n")
    L.append("- gnomAD-absent is not AF=0; absence is analyzed as a category, not a value.\n")
    L.append("- AF↔predictor correlation is partly circular (predictors are conservation-trained).\n")
    L.append("\n## 9. Data limitations\n\n")
    L.append("- Single gene (BRCA1), single transcript (NM_007294.4), missense VUS only.\n")
    L.append("- CADD excluded (v1.7 bulk annotation unavailable).\n")
    L.append("- Full per-population AF not retrieved (gnomAD GraphQL cost limit); faf95 popmax used.\n")
    L.append("\n## 10. Phase 6 recommendation\n\n")
    L.append("Prioritize literature review of classes C (absent + strong impact) and E "
             "(disagreement), with F (elevated population frequency) as a cross-check group. "
             "No variant is claimed pathogenic or benign.\n")
    with open(os.path.join(REPORTS, "phase5_statistical_analysis.md"), "w") as fh:
        fh.write("".join(L))


if __name__ == "__main__":
    main()
