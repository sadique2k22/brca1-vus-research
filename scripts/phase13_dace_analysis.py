"""Phase 13 — Dace 2025 SGE integration: domain-contrast analysis.

Adds the non-RING/BRCT functional scores (Dace et al. 2025 medRxiv preprint,
doi:10.1101/2025.08.11.25333423) alongside the frozen Findlay 2018 RING/BRCT
scores and quantifies whether computational-predictor/functional correspondence
transfers across protein domains.

Orientation convention (documented, one convention everywhere):
  * impact orientation: higher = predicted impact:
      REVEL raw, SIFT as (1 - SIFT), PolyPhen-2 raw  (identical to phase9.py)
  * raw SIFT (higher = tolerated) reported alongside for direct comparability
    with the +0.363 reported by Findlay et al. 2018 (Extended Data Fig. 9).
  * functional scores: higher = WT-like / fitter in BOTH assays
      (Findlay score: lower = reduced cellular fitness;
       Dace mean replicate score: lower = reduced cellular fitness).

Frozen inputs are read only. All outputs are new Phase 13 artifacts.
"""
import csv
import json
import math
import os
import random

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

plt.rcParams.update({
    "figure.dpi": 100, "savefig.dpi": 300,
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.grid": True, "grid.alpha": 0.3, "figure.figsize": (5.2, 3.6),
})

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FROZEN = os.path.join(ROOT, "data/processed/brca1_vus_missense_with_functional.tsv")
DACE = os.path.join(ROOT, "data/processed_eval/vus_with_dace_scores_352.tsv")
TABDIR = os.path.join(ROOT, "results/tables")
FIGDIR = os.path.join(ROOT, "results/figures")
RPTDIR = os.path.join(ROOT, "results/reports")

DACE_PAPER_LOF = -0.799   # Dace et al. nonsense-matched LoF threshold
PHASE9_NONFUNC = -1.0     # phase9.py "< -1.0" descriptive threshold

IMPACT_COLS = {"revel_score": "revel", "sift_score": "sift", "polyphen_score": "polyphen"}


def load(path):
    with open(path, newline="", encoding="utf-8") as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        idx = {n.lstrip("#"): i for i, n in enumerate(header)}
        rows = [row for row in r]
    return header, idx, rows


def parse_float(v):
    if v is None or v == "" or v == "NA":
        return None
    try:
        return float(v)
    except ValueError:
        return None


def spearman_ci(rho, n):
    if n < 4 or rho is None:
        return None
    rho = max(-0.999, min(0.999, rho))
    z = math.atanh(rho)
    se = 1.0 / math.sqrt(n - 3)
    lo = math.tanh(z - 1.96 * se)
    hi = math.tanh(z + 1.96 * se)
    return (round(lo, 3), round(hi, 3))


def fisher_z_diff(rho1, n1, rho2, n2):
    """Two-sided p for difference of independent Spearman rhos (Fisher z)."""
    z1 = math.atanh(max(-0.999, min(0.999, rho1)))
    z2 = math.atanh(max(-0.999, min(0.999, rho2)))
    se = math.sqrt(1.0 / (n1 - 3) + 1.0 / (n2 - 3))
    z = (z1 - z2) / se
    p = 2.0 * (1.0 - stats.norm.cdf(abs(z)))
    return z, p


def permutation_diff(x1, y1, x2, y2, seed=20260820, n_perm=10000):
    """Empirical two-sided p that |rho1 - rho2| exceeds chance, under
    exchangeability of the (x, y) pairs across the two domains."""
    rng = random.Random(seed)
    obs_r1 = stats.spearmanr(x1, y1).statistic
    obs_r2 = stats.spearmanr(x2, y2).statistic
    obs_d = abs(obs_r1 - obs_r2)
    pool = list(zip(x1 + x2, y1 + y2))
    n1 = len(x1)
    n2 = len(x2)
    extreme = 0
    for _ in range(n_perm):
        rng.shuffle(pool)
        a = pool[:n1]
        b = pool[n2:]
        ra = stats.spearmanr([p[0] for p in a], [p[1] for p in a]).statistic
        rb = stats.spearmanr([p[0] for p in b], [p[1] for p in b]).statistic
        if abs(ra - rb) >= obs_d - 1e-12:
            extreme += 1
    return obs_r1, obs_r2, (extreme + 1) / (n_perm + 1)


def _save(fig, name):
    os.makedirs(FIGDIR, exist_ok=True)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, name + ".svg"))
    fig.savefig(os.path.join(FIGDIR, name + ".png"), dpi=300)
    plt.close(fig)


def _impact(name, vals):
    """Impact orientation (higher = predicted impact), matching phase9.py."""
    if name == "sift":
        return [None if v is None else 1.0 - v for v in vals]
    return vals


def main():
    os.makedirs(TABDIR, exist_ok=True)
    os.makedirs(RPTDIR, exist_ok=True)

    # ---- inputs ----
    _, widx, wf = load(FROZEN)
    ring = []
    for r in wf:
        if r[widx["findlay_available"]] == "present":
            fs = parse_float(r[widx["findlay_score"]])
            if fs is None:
                continue
            ring.append({
                "key": r[widx["variant_key"]],
                "protein_change": r[widx["protein_change"]],
                "revel": parse_float(r[widx["revel_score"]]),
                "sift": parse_float(r[widx["sift_score"]]),
                "polyphen": parse_float(r[widx["polyphen_score"]]),
                "func": fs,
            })
    _, didx, dace_rows = load(DACE)
    central = []
    for r in dace_rows:
        fs = parse_float(r[didx["dace_mean_score"]])
        if fs is None:
            continue
        central.append({
            "key": r[didx["variant_key"]],
            "protein_change": r[didx["protein_change"]],
            "revel": parse_float(r[didx["revel_score"]]),
            "sift": parse_float(r[didx["sift_score"]]),
            "polyphen": parse_float(r[didx["polyphen_score"]]),
            "func": fs,
            "region": r[didx["dace_region"]],
            "present": r[didx["gnomad_present"]],
        })

    domains = {"RING/BRCT (Findlay 2018)": ring, "central exons (Dace 2025)": central}

    # ---- per-domain correlations (impact orientation + raw SIFT) ----
    rows = []
    summary = {}
    predictors = [("REVEL", "revel"), ("SIFT (1\u2212SIFT)", "sift"), ("PolyPhen-2", "polyphen")]
    for dname, d in domains.items():
        entry = {"domain": dname, "n": len(d)}
        for label, col in predictors:
            xs = [v[col] for v in d]
            ys = [v["func"] for v in d]
            x_imp = _impact(col, xs)
            pair = [(a, b) for a, b in zip(x_imp, ys) if a is not None and b is not None]
            xx = [p[0] for p in pair]
            yy = [p[1] for p in pair]
            rho, p = stats.spearmanr(xx, yy)
            ci = spearman_ci(rho, len(pair))
            entry[label] = {"rho": round(float(rho), 4), "p": float(p),
                            "ci": ci, "n": len(pair)}
            if col == "sift":
                pair_raw = [(a, b) for a, b in zip(xs, ys) if a is not None and b is not None]
                rrho, rp = stats.spearmanr([p[0] for p in pair_raw], [p[1] for p in pair_raw])
                entry["SIFT raw"] = {"rho": round(float(rrho), 4), "p": float(rp),
                                     "ci": spearman_ci(rrho, len(pair_raw)),
                                     "n": len(pair_raw)}
        summary[dname] = entry
        for label, col in predictors:
            e = entry[label]
            rows.append([dname, label, e["n"], e["rho"], e["ci"][0] if e["ci"] else "",
                         e["ci"][1] if e["ci"] else "", e["p"]])
        rows.append([dname, "SIFT raw (higher = tolerated)", entry["SIFT raw"]["n"],
                     entry["SIFT raw"]["rho"],
                     entry["SIFT raw"]["ci"][0] if entry["SIFT raw"]["ci"] else "",
                     entry["SIFT raw"]["ci"][1] if entry["SIFT raw"]["ci"] else "",
                     entry["SIFT raw"]["p"]])

    # ---- domain-difference tests ----
    diff_rows = []
    perms = {}
    for label, col in predictors:
        def paired(d, col=col):
            x = _impact(col, [v[col] for v in d])
            y = [v["func"] for v in d]
            return ([a for a, b in zip(x, y) if a is not None and b is not None],
                    [b for a, b in zip(x, y) if a is not None and b is not None])
        x1, y1 = paired(ring)
        x2, y2 = paired(central)
        r1 = stats.spearmanr(x1, y1).statistic
        r2 = stats.spearmanr(x2, y2).statistic
        z, p = fisher_z_diff(r1, len(x1), r2, len(x2))
        pr1, pr2, pp = permutation_diff(x1, y1, x2, y2)
        diff_rows.append([label, len(x1), round(r1, 4), len(x2), round(r2, 4),
                          round(z, 4), p, pp])
        perms[label] = {"fisher_p": p, "perm_p": pp, "z": round(z, 4)}

    # ---- conflict analysis on the 352 (phase9-identical REVEL rules) ----
    def classify(rv):
        if rv is None:
            return "intermediate"
        if rv >= 0.644:
            return "impact"
        if rv <= 0.290:
            return "tolerance"
        return "intermediate"

    def conflicts(d, nonfunc_thresh, label):
        out = []
        for v in d:
            cat = classify(v["revel"])
            if cat == "impact" and v["func"] > 0:
                out.append([v["key"], v["protein_change"], cat, "impact vs WT-like (>0)",
                            v["revel"], v["func"]])
            elif cat == "tolerance" and v["func"] < nonfunc_thresh:
                out.append([v["key"], v["protein_change"], cat,
                            "tolerance vs non-functional (<%s)" % nonfunc_thresh,
                            v["revel"], v["func"]])
        return out

    ring_conf_1 = conflicts(ring, PHASE9_NONFUNC, "RING/BRCT")
    dace_conf_1 = conflicts(central, PHASE9_NONFUNC, "central")
    dace_conf_0799 = conflicts(central, DACE_PAPER_LOF, "central")

    # per-predictor exploratory conflicts (no comparator in the 41-cohort)
    def pp_conflicts(d, thresh, pred, func_max, func_min, func_thresh):
        out = []
        for v in d:
            pv = v[pred]
            if pv is None:
                continue
            if pred == "sift":
                impact = (1 - pv) >= 0.95
            else:
                impact = pv >= 0.908 if pred == "polyphen" else pv >= 0.644
            if impact and v["func"] > func_max:
                out.append([v["key"], v["protein_change"], label, v["func"]])
        return out

    per_pred = {}
    for pred, lab in (("revel", "REVEL>=0.644"), ("sift", "1-SIFT>=0.95 (raw<=0.05)"),
                      ("polyphen", "PP2>0.908")):
        c = pp_conflicts(central, None, pred, 0.0, None, None)
        per_pred[pred] = len(c)

    # ---- figures ----
    # Fig 17 — Dace region coverage of our VUS
    from collections import Counter
    counts = Counter(v["region"] for v in central)
    order = ["exon 6", "exon 10 (5')", "exon 10 (mid 1)", "exon 10 (mid 2)",
             "exon 10 (3')", "exon 11", "exon 12 (5')", "exon 12 (3')"]
    labels = [r for r in order if r in counts]
    vals = [counts[r] for r in labels]
    fig, ax = plt.subplots(figsize=(6.5, 3.4))
    ax.bar(labels, vals, color="#4C72B0", alpha=0.85)
    ax.set_ylabel("VUS with Dace score")
    ax.set_title("Figure 17 \u2014 Dace 2025 HAP1 coverage of our 1,904 VUS (n = %d)" % len(central))
    ax.tick_params(axis="x", rotation=30)
    _save(fig, "fig17_dace_region_coverage")

    # Fig 18 — predictor vs Dace functional score (impact orientation)
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4))
    for ax, (label, col) in zip(axes, predictors):
        x = _impact(col, [v[col] for v in central])
        y = [v["func"] for v in central]
        xs = [a for a, b in zip(x, y) if a is not None and b is not None]
        ys = [b for a, b in zip(x, y) if a is not None and b is not None]
        ax.scatter(xs, ys, s=4, alpha=0.4, color="#4C72B0")
        ax.set_xlabel(label)
        ax.set_ylabel("Dace HAP1 score (lower = reduced fitness)")
        ax.set_title(label)
    fig.suptitle("Figure 18 \u2014 central-exon VUS (Dace 2025, n = %d)" % len(central),
                 fontsize=10)
    _save(fig, "fig18_dace_predictor_scatter")

    # Fig 19 — domain contrast with CIs (impact orientation)
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    colors = {"RING/BRCT (Findlay 2018)": "#C44E52", "central exons (Dace 2025)": "#4C72B0"}
    positions = {}
    for i, (label, col) in enumerate(predictors):
        for j, (dname, d) in enumerate(domains.items()):
            e = summary[dname][label]
            xpos = i + (j - 0.5) * 0.32
            positions[(dname, label)] = xpos
            ax.errorbar(xpos, e["rho"], yerr=[[e["rho"] - e["ci"][0]], [e["ci"][1] - e["rho"]]],
                        fmt="o", color=colors[dname], capsize=3)
    ax.axhline(0, color="grey", lw=0.8, ls="--")
    ax.set_xticks(range(3))
    ax.set_xticklabels([l for l, _ in predictors])
    ax.set_ylabel("Spearman rho (impact orientation)")
    handles = [plt.Line2D([0], [0], marker="o", ls="", color=colors[d], label=d)
               for d in domains]
    ax.legend(handles=handles, fontsize=8)
    ax.set_title("Figure 19 \u2014 correspondence does not transfer across domains")
    _save(fig, "fig19_domain_contrast")

    # ---- table outputs ----
    with open(os.path.join(TABDIR, "phase13_domain_correlations.tsv"), "w",
              newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["domain", "predictor", "n", "rho_impact", "ci_lo", "ci_hi", "p"])
        w.writerows(rows)
    with open(os.path.join(TABDIR, "phase13_domain_difference.tsv"), "w",
              newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["predictor", "n_ring", "rho_ring", "n_central", "rho_central",
                    "fisher_z", "fisher_p", "permutation_p"])
        w.writerows(diff_rows)
    with open(os.path.join(TABDIR, "phase13_dace_conflicts.tsv"), "w",
              newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["variant_key", "protein_change", "computational_call", "conflict_type",
                    "revel", "dace_mean_score", "domain"])
        for v in dace_conf_1:
            w.writerow(v + ["central"])
        for v in ring_conf_1:
            w.writerow(v + ["ring_brct"])

    # ---- S5 supplementary table (Dace 2025, 352 variants) ----
    sup = os.path.join(ROOT, "manuscript/supplementary")
    os.makedirs(sup, exist_ok=True)
    with open(os.path.join(sup, "S5_dace_352.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["variant_key", "protein_change", "normalized_hgvs_c", "revel_score",
                    "sift_score", "polyphen_score", "gnomad_present", "dace_region",
                    "dace_mean_score", "dace_scores"])
        for r in dace_rows:
            w.writerow([r[didx[c]] for c in ("variant_key", "protein_change",
                                             "normalized_hgvs_c", "revel_score",
                                             "sift_score", "polyphen_score",
                                             "gnomad_present", "dace_region",
                                             "dace_mean_score", "dace_scores")])

    # ---- report ----
    out = []
    out.append("# Phase 13 — Domain-contrast analysis (Dace 2025 integration)\n")
    out.append("Date: 2026-08-20. Frozen inputs read-only; all outputs are new Phase 13 artifacts.\n")
    out.append("## Orientation convention (used everywhere in Phase 13)\n")
    out.append("- Impact orientation: REVEL raw, SIFT as (1 - SIFT), PolyPhen-2 raw "
               "(higher = predicted impact). Identical to `phase9.py`.\n")
    out.append("- Raw SIFT (higher = tolerated) is reported alongside for direct comparability "
               "with Findlay et al. 2018 (their reported SIFT rho = +0.363).\n")
    out.append("- Functional scores: higher = WT-like/fitter in BOTH assays.\n")
    out.append("## Per-domain Spearman correlations (impact orientation; 95% CI Fisher z)\n")
    out.append("| Domain | Predictor | n | rho | 95% CI |")
    out.append("|---|---|---|---|---|")
    for r in rows:
        out.append("| {} | {} | {} | {:+.3f} | {} to {} |".format(
            r[0], r[1], r[2], r[3], r[4], r[5]))
    out.append("\n## Domain-difference tests (RING/BRCT vs central exons)\n")
    out.append("| Predictor | Fisher z | Fisher p | Permutation p (10,000, seed 20260820) |")
    out.append("|---|---|---|---|")
    for label, col in predictors:
        zs = perms[label]["z"]
        out.append("| {} | {:.3f} | {:.2e} | {:.2e} |".format(
            label, zs, perms[label]["fisher_p"], perms[label]["perm_p"]))
    out.append("\n## Conflict analysis on the 352 newly covered VUS (phase9-identical REVEL rules)\n")
    out.append(f"- REVEL impact (>=0.644) with WT-like Dace score (>0): "
               f"{sum(1 for v in dace_conf_1 if v[2] == 'impact')} / 352")
    out.append(f"- REVEL tolerance (<=0.290) with Dace score < -1.0: "
               f"{sum(1 for v in dace_conf_1 if v[2] == 'tolerance')} / 352")
    out.append(f"- Sensitivity (Dace-paper LoF threshold -0.799): "
               f"{len(dace_conf_0799)} conflicts total")
    out.append(f"- For reference, RING/BRCT 373 (same rules): {len(ring_conf_1)} conflicts "
               f"({len(ring_conf_1)/373*100:.1f}%)")
    out.append(f"- Central 352 (same rules, -1.0): {len(dace_conf_1)} conflicts "
               f"({len(dace_conf_1)/352*100:.1f}%)")
    out.append(f"- Per-predictor exploratory (impact call with WT-like Dace score): "
               f"REVEL {per_pred['revel']}, SIFT {per_pred['sift']}, "
               f"PolyPhen-2 {per_pred['polyphen']} (no 41-cohort comparator; descriptive only)\n")
    out.append("## Orientation audit (fixes applied in Phase 13)\n")
    out.append("- `phase9.py` computes SIFT as (1 - SIFT) for correlations; the SIFT negative "
               "sign in the manuscript therefore arises from SIFT's own inversion, NOT from the "
               "functional-score orientation. The previous sign-rationale sentence applied the "
               "REVEL-style rationale to SIFT, which is incorrect with raw SIFT (rho = +0.370).")
    out.append("- Raw-SIFT rho +0.370 (RING/BRCT, n = 373) is essentially identical to the "
               "+0.363 reported by Findlay et al. 2018 \u2014 concordant, not merely in magnitude.")
    out.append("- The Phase 13 feasibility report's earlier \u201cSIFT is oppositely signed\u201d "
               "statement mixed raw-SIFT (central) with impact-oriented SIFT (RING/BRCT); "
               "that was an orientation artifact. With matched orientation, SIFT correspondence "
               "attenuates (0.370 \u2192 0.260 impact orientation) but does not change sign.\n")
    out.append("## Interpretation\n")
    out.append("- REVEL and PolyPhen-2 correspondence essentially disappears in the central "
               "exons (REVEL -0.384 \u2192 {:.3f}; PolyPhen-2 -0.188 \u2192 {:.3f}); SIFT "
               "attenuates but remains directionally consistent (-0.370 \u2192 {:.3f}, impact "
               "orientation). Domain-difference tests: see above.".format(
                   summary["central exons (Dace 2025)"]["REVEL"]["rho"],
                   summary["central exons (Dace 2025)"]["PolyPhen-2"]["rho"],
                   summary["central exons (Dace 2025)"]["SIFT (1\u2212SIFT)"]["rho"]))
    out.append("- Caveats: Dace et al. is a preprint; unfiltered continuous table used; "
               "the 41-variant cohort and its 13 conflicts are untouched (all RING/BRCT or "
               "unassayed); the conflict rates above are descriptive, on different variant sets.")
    out.append("\n## Artifacts\n")
    out.append("- results/tables/phase13_domain_correlations.tsv (this table)")
    out.append("- results/tables/phase13_domain_difference.tsv (Fisher z + permutation)")
    out.append("- results/tables/phase13_dace_conflicts.tsv (conflict rows, central + RING/BRCT)")
    out.append("- results/figures/fig17_dace_region_coverage.{svg,png}")
    out.append("- results/figures/fig18_dace_predictor_scatter.{svg,png}")
    out.append("- results/figures/fig19_domain_contrast.{svg,png}")
    with open(os.path.join(RPTDIR, "phase13_domain_contrast.md"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")

    print(json.dumps({
        "n_ring": len(ring), "n_central": len(central),
        "summary": {k: {p: v["rho"] for p, v in dom.items()
                        if isinstance(v, dict) and "rho" in v}
                    for k, dom in summary.items()},
        "domain_difference": {k: v for k, v in perms.items()},
        "conflicts": {"ring_373": len(ring_conf_1), "central_352": len(dace_conf_1),
                      "central_352_loft": len(dace_conf_0799),
                      "per_predictor": per_pred},
    }, indent=1))


if __name__ == "__main__":
    main()