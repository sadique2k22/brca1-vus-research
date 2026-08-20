#!/usr/bin/env python3
"""Phase 13B figure: calibration AUC + AM/BayesDel domain contrast.

fig20: 2-panel figure.
  (a) AUC of each in silico predictor separating ClinVar gold-standard
      P/LP vs B/LB BRCA1 missense controls (bootstrap 95% CI).
  (b) Domain-contrast Spearman rho (impact orientation) for AlphaMissense
      and BayesDel noAF vs functional scores (RING/BRCT Findlay 2018 vs
      central-exon Dace 2025), with 95% CI.

Conventions match Phase 13 (scripts/phase13_dace_analysis.py).
"""
import csv
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.phase13b_calibration import auc_bootstrap, fisher_z_diff  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABDIR = os.path.join(ROOT, "results/tables")
FIGDIR = os.path.join(ROOT, "results/figures")

plt.rcParams.update({
    "figure.dpi": 100, "savefig.dpi": 300,
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.grid": True, "grid.alpha": 0.3,
})
COLORS = {"RING/BRCT": "#C44E52", "central": "#4C72B0"}


def _save(fig, name):
    os.makedirs(FIGDIR, exist_ok=True)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, name + ".svg"))
    fig.savefig(os.path.join(FIGDIR, name + ".png"), dpi=300)
    plt.close(fig)


def spearman_ci(rho, n):
    if n < 4 or rho is None:
        return None
    rho = max(-0.999, min(0.999, rho))
    z = math.atanh(rho)
    se = 1.0 / math.sqrt(n - 3)
    lo = math.tanh(z - 1.96 * se)
    hi = math.tanh(z + 1.96 * se)
    return (round(lo, 3), round(hi, 3))


def _load_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        rows = [dict(row) for row in r]
    return rows


def main():
    ctrl = os.path.join(ROOT, "data/processed_eval/controls_p_lp_b_lb.tsv")
    corr = os.path.join(TABDIR, "phase13b_domain_correlations.tsv")

    plp, blb = [], []
    with open(ctrl) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            d = (plp if row["label"] == "P/LP" else blb)
            d.append(row)

    def col(rows, name):
        vals = []
        for row in rows:
            v = row.get(name, "")
            vals.append(float(v) if v not in ("", "nan") else math.nan)
        return np.array(vals)

    aucs = []
    for name, c in [("AlphaMissense", "am_pathogenicity"),
                    ("BayesDel noAF", "bayesdel_noAF"),
                    ("REVEL v1.3", "revel_score")]:
        res = auc_bootstrap(col(plp, c), col(blb, c))
        if res:
            aucs.append((name, res))
    has_revel = any(n == "REVEL v1.3" for n, _ in aucs)

    # panel (a) AUC
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(10.5, 3.8))
    names = [n for n, _ in aucs]
    vals = [r["auc"] for _, r in aucs]
    los = [r["auc"] - r["ci_lo"] for _, r in aucs]
    his = [r["ci_hi"] - r["auc"] for _, r in aucs]
    xpos = np.arange(len(names))
    axa.errorbar(xpos, vals, yerr=[los, his], fmt="o", capsize=4,
                 color="#4C72B0", markersize=7)
    axa.axhline(0.5, color="grey", lw=0.8, ls="--")
    axa.set_xticks(xpos)
    axa.set_xticklabels(names, rotation=15, ha="right")
    axa.set_ylabel("AUC (P/LP vs B/LB, bootstrap 95% CI)")
    axa.set_title("(a) Calibration: 159 P/LP vs 188 B/LB BRCA1 missense")
    axa.set_ylim(0.5, 1.02)

    # panel (b) domain contrast rho
    corr_rows = _load_tsv(corr)
    xlabels = ["AlphaMissense", "BayesDel noAF"]
    for i, pred in enumerate(xlabels):
        for j, dom in enumerate(["RING/BRCT", "central"]):
            m = [r for r in corr_rows
                 if r["predictor"] == pred and dom in r["set"]]
            if not m:
                continue
            r = m[0]
            rho = float(r["rho"])
            n = int(r["n"])
            ci = spearman_ci(rho, n)
            xpos = i + (j - 0.5) * 0.32
            axb.errorbar(xpos, rho, yerr=[[rho - ci[0]], [ci[1] - rho]],
                         fmt="o", color=COLORS[dom], capsize=3)
    axb.axhline(0, color="grey", lw=0.8, ls="--")
    axb.set_xticks(range(2))
    axb.set_xticklabels(xlabels)
    axb.set_ylabel("Spearman rho (impact orientation)")
    axb.set_title("(b) Domain contrast: predictor vs functional score")
    handles = [plt.Line2D([0], [0], marker="o", ls="", color=COLORS[d], label=d)
               for d in ("RING/BRCT", "central")]
    axb.legend(handles=handles, fontsize=8)

    fig.suptitle("Figure 20 \u2014 Phase 13B: in silico calibration + domain contrast",
                 fontsize=10)
    _save(fig, "fig20_phase13b_calibration")


if __name__ == "__main__":
    main()