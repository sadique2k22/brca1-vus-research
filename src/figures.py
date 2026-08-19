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


# ---- Phase 7 figures ----

def fig11_cohort(stratum_counts):
    fig, ax = plt.subplots(figsize=(4.5, 3.2))
    order = ["A", "B", "C", "D", "E"]
    vals = [stratum_counts.get(k, 0) for k in order]
    ax.bar(order, vals, color="#4C72B0", edgecolor="k")
    ax.set_xlabel("Stratum")
    ax.set_ylabel("Variants")
    ax.set_title("Figure 11 — final cohort composition")
    _save(fig, "fig11_cohort_composition")


def fig12_predictor_functional(revel, fscore):
    fig, ax = plt.subplots()
    _scatter(ax, revel, fscore, "REVEL", "Findlay function score (lower = LOF)",
             "Figure 12 — REVEL vs Findlay functional score")
    _save(fig, "fig12_predictor_vs_functional")


def fig13_concordance(evidence):
    import numpy as np
    # evidence rows: [.., stratum(2), gnomAD(3), comp(4), findlay(5), ..., conflict(12)]
    cats = ["impact+LOF", "impact+normal", "tol+LOF", "tol+normal", "intermediate"]
    mat = np.zeros((5, 5))
    for e in evidence:
        comp = e[4]
        func = e[5]
        i = {"impact": 0, "tolerance": 1, "intermediate": 2}.get(comp, 2)
        j = {"LOF": 0, "normal": 1, "NA": 2}.get(func, 2)
        mat[i, j] += 1
    fig, ax = plt.subplots(figsize=(4.5, 3.5))
    im = ax.imshow(mat[:3, :2], cmap="Blues")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["LOF", "normal"])
    ax.set_yticks([0, 1, 2]); ax.set_yticklabels(["impact", "tolerance", "intermediate"])
    for i in range(3):
        for j in range(2):
            ax.text(j, i, int(mat[i, j]), ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax)
    ax.set_title("Figure 13 — computational vs functional concordance")
    _save(fig, "fig13_concordance")


def fig14_functional_by_category(scored):
    # scored: list of (fs, rv, sv, pv)
    def cat(rv):
        if rv is None:
            return "intermediate"
        return "impact" if rv >= 0.644 else ("tolerance" if rv <= 0.290 else "intermediate")
    groups = {"impact": [], "tolerance": [], "intermediate": []}
    for fs, rv, sv, pv in scored:
        groups[cat(rv)].append(fs)
    fig, ax = plt.subplots()
    data = [groups["impact"], groups["tolerance"], groups["intermediate"]]
    ax.violinplot(data, showmedians=True)
    ax.set_xticks([1, 2, 3]); ax.set_xticklabels(["impact", "tolerance", "intermediate"])
    ax.set_ylabel("Findlay function score")
    ax.set_title("Figure 14 — functional score by REVEL category")
    _save(fig, "fig14_functional_by_category")


def fig15_evidence_availability(evidence):
    import numpy as np
    # evidence rows: [.., findlay(5), other(6), clinical(7), segregation(8), expert(9), ...]
    n = len(evidence)
    avail = {
        "Findlay": sum(1 for e in evidence if e[5] in ("LOF", "normal")),
        "other functional": sum(1 for e in evidence if e[6] != "none identified"),
        "clinical cases": sum(1 for e in evidence if e[7] != "not assessed"),
        "segregation": sum(1 for e in evidence if e[8] != "not assessed"),
        "expert curation": sum(1 for e in evidence if e[9] not in ("NA", "")),
        "literature": sum(1 for e in evidence if e[11] > 0),
    }
    fig, ax = plt.subplots(figsize=(5.5, 3.2))
    names = list(avail.keys())
    ax.barh(names, [avail[k] for k in names], color="#55A868", edgecolor="k")
    ax.set_xlabel("variants (of %d)" % n)
    ax.set_title("Figure 15 — evidence availability across cohort")
    _save(fig, "fig15_evidence_availability")


def fig16_domain(ring, brct):
    fig, ax = plt.subplots()
    data = [ring, brct]
    ax.violinplot(data, showmedians=True)
    ax.set_xticks([1, 2]); ax.set_xticklabels(["RING", "BRCT"])
    ax.set_ylabel("Findlay function score")
    ax.set_title("Figure 16 — functional score by domain")
    _save(fig, "fig16_domain_comparison")
