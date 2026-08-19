"""Phase 5 figure generation (matplotlib, headless Agg). Vector (SVG) + PNG (300 dpi)."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "figure.dpi": 100, "savefig.dpi": 300,
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.grid": True, "grid.alpha": 0.3, "figure.figsize": (5.2, 3.6),
})

OUTDIR = "results/figures"


def _save(fig, name):
    os.makedirs(OUTDIR, exist_ok=True)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, name + ".svg"))
    fig.savefig(os.path.join(OUTDIR, name + ".png"), dpi=300)
    plt.close(fig)


def fig1_workflow():
    fig, ax = plt.subplots(figsize=(7.2, 2.2))
    ax.axis("off")
    steps = ["ClinVar", "BRCA1", "GRCh38", "VUS", "missense", "normalized",
             "annotated", "1,904 analyzed"]
    n = len(steps)
    for i, s in enumerate(steps):
        ax.add_patch(plt.Rectangle((i / n, 0.3), 0.9 / n, 0.4, facecolor="#4C72B0",
                                   edgecolor="k", alpha=0.85))
        ax.text(i / n + 0.45 / n, 0.5, s, ha="center", va="center", fontsize=7, color="white")
        if i < n - 1:
            ax.annotate("", xy=((i + 1) / n, 0.5), xytext=(i / n + 0.9 / n, 0.5),
                        arrowprops=dict(arrowstyle="->", lw=1.2, color="k"))
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0, 1)
    ax.set_title("Figure 1 — Study workflow: filtering to the analysis set")
    _save(fig, "fig1_workflow")


def _hist(ax, values, bins, color, xlabel, title):
    vals = [v for v in values if v is not None]
    ax.hist(vals, bins=bins, color=color, edgecolor="k", alpha=0.85)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    ax.set_title(title)


def fig2_revel(revel):
    fig, ax = plt.subplots()
    _hist(ax, revel, np.linspace(0, 1, 31), "#4C72B0", "REVEL score", "Figure 2 — REVEL distribution")
    for x in (0.290, 0.644, 0.932):
        ax.axvline(x, color="r", ls="--", lw=0.8)
    _save(fig, "fig2_revel_distribution")


def fig3_sift(sift):
    fig, ax = plt.subplots()
    _hist(ax, sift, np.linspace(0, 1, 31), "#DD8452", "SIFT score", "Figure 3 — SIFT distribution")
    ax.axvline(0.05, color="r", ls="--", lw=0.8)
    _save(fig, "fig3_sift_distribution")


def fig4_polyphen(polyphen):
    fig, ax = plt.subplots()
    _hist(ax, polyphen, np.linspace(0, 1, 31), "#55A868", "PolyPhen-2 score", "Figure 4 — PolyPhen-2 distribution")
    for x in (0.446, 0.908):
        ax.axvline(x, color="r", ls="--", lw=0.8)
    _save(fig, "fig4_polyphen_distribution")


def fig5_agreement(matrix, row_labels, col_labels, title):
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    im = ax.imshow(matrix, cmap="Blues", aspect="auto")
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=40, ha="right", fontsize=7)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=7)
    for i in range(len(row_labels)):
        for j in range(len(col_labels)):
            ax.text(j, i, int(matrix[i, j]), ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, label="count")
    ax.set_title(title)
    _save(fig, "fig5_agreement_matrix")


def _scatter(ax, x, y, xlab, ylab, title):
    xs, ys = [], []
    for a, b in zip(x, y):
        if a is not None and b is not None:
            xs.append(a)
            ys.append(b)
    ax.scatter(xs, ys, s=4, alpha=0.4, color="#4C72B0")
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)


def fig6_revel_polyphen(revel, polyphen):
    fig, ax = plt.subplots()
    _scatter(ax, revel, polyphen, "REVEL", "PolyPhen-2", "Figure 6 — REVEL vs PolyPhen-2")
    _save(fig, "fig6_revel_vs_polyphen")


def fig7_revel_sift(revel, sift):
    fig, ax = plt.subplots()
    _scatter(ax, revel, sift, "REVEL", "SIFT (lower = more damaging)", "Figure 7 — REVEL vs SIFT")
    _save(fig, "fig7_revel_vs_sift")


def fig8_present_absent(present_scores, absent_scores, name, ylab):
    fig, ax = plt.subplots()
    data = [present_scores, absent_scores]
    parts = ax.violinplot(data, showmedians=True)
    for pc in parts["bodies"]:
        pc.set_facecolor("#4C72B0")
        pc.set_alpha(0.7)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["gnomAD present", "gnomAD absent"])
    ax.set_ylabel(ylab)
    ax.set_title(f"Figure 8 — {name} by gnomAD presence")
    _save(fig, f"fig8_{name.lower().replace(' ', '_')}_by_presence")


def fig9_af_revel(log10af, revel):
    fig, ax = plt.subplots()
    xs, ys = [], []
    for a, b in zip(log10af, revel):
        if a is not None and b is not None:
            xs.append(a)
            ys.append(b)
    ax.scatter(xs, ys, s=5, alpha=0.5, color="#4C72B0")
    ax.set_xlabel("log10(global AF)")
    ax.set_ylabel("REVEL")
    ax.set_title("Figure 9 — log10(global AF) vs REVEL (gnomAD-present)")
    _save(fig, "fig9_log10_af_vs_revel")


def fig10_population_freq(global_af, popmax_af):
    fig, ax = plt.subplots()
    xs, ys = [], []
    for a, b in zip(global_af, popmax_af):
        if a is not None and b is not None and a > 0 and b > 0:
            xs.append(a)
            ys.append(b)
    ax.scatter(xs, ys, s=5, alpha=0.5, color="#DD8452")
    ax.plot([min(xs + [1e-6]), max(xs)], [min(xs + [1e-6]), max(xs)], "k--", lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("global AF")
    ax.set_ylabel("max population filtering AF (faf95 popmax)")
    ax.set_title("Figure 10 — global vs population frequency")
    _save(fig, "fig10_population_frequency")
