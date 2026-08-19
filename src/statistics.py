"""Descriptive statistics, predictor agreement, and correlation helpers (Phase 5).

Pure functions operating on plain values for testability. Directions are handled
explicitly: SIFT is lower-is-more-damaging, REVEL/PolyPhen are higher-is-more-damaging.
"""
import math

import scipy.stats as st


def parse_float(v):
    """Parse a numeric cell; empty string or 'NA' -> None."""
    if v is None or v == "" or v == "NA":
        return None
    try:
        return float(v)
    except (ValueError, TypeError):
        return None


def categorize_revel(score):
    if score is None:
        return None
    if score >= 0.644:
        return "impact"
    if score <= 0.290:
        return "tolerance"
    return "intermediate"


def categorize_sift(score):
    if score is None:
        return None
    return "deleterious" if score <= 0.05 else "tolerated"


def categorize_polyphen(score):
    if score is None:
        return None
    if score >= 0.908:
        return "probably_damaging"
    if score <= 0.446:
        return "benign"
    return "possibly_damaging"


def impact_2class(kind, score):
    """Map a predictor to 2-class 'impact'/'tolerance' (None if indeterminate).

    REVEL intermediate (0.290 < s < 0.644) maps to None (excluded from 2-class agreement).
    """
    if score is None:
        return None
    if kind == "revel":
        return {"impact": "impact", "tolerance": "tolerance"}.get(categorize_revel(score))
    if kind == "sift":
        return "impact" if categorize_sift(score) == "deleterious" else "tolerance"
    if kind == "polyphen":
        return "impact" if categorize_polyphen(score) in ("possibly_damaging", "probably_damaging") else "tolerance"
    return None


def descriptive(values):
    vals = [v for v in values if v is not None]
    if not vals:
        return {"n": 0}
    q1, q2, q3 = st.mstats.mquantiles(vals, [0.25, 0.5, 0.75])
    return {
        "n": len(vals),
        "min": float(min(vals)),
        "max": float(max(vals)),
        "mean": float(st.tmean(vals)),
        "median": float(q2),
        "std": float(st.tstd(vals)),
        "q1": float(q1),
        "q3": float(q3),
    }


def spearman(x, y):
    xs, ys = _paired(x, y)
    if len(xs) < 3:
        return None
    rho, p = st.spearmanr(xs, ys)
    return {"rho": float(rho), "p": float(p), "n": len(xs)}


def cohens_kappa(a, b):
    aa, bb = _paired(a, b)
    if len(aa) < 2:
        return None
    cats = sorted(set(aa) | set(bb))
    n = len(aa)
    table = {}
    for x, y in zip(aa, bb):
        table[(x, y)] = table.get((x, y), 0) + 1
    po = sum(table.get((c, c), 0) for c in cats) / n
    pe = 0.0
    for c in cats:
        row = sum(table.get((c, y), 0) for y in cats) / n
        col = sum(table.get((x, c), 0) for x in cats) / n
        pe += row * col
    if abs(1 - pe) < 1e-12:
        return {"kappa": 1.0, "p": 0.0, "n": n}
    kappa = (po - pe) / (1 - pe)
    se = math.sqrt(po * (1 - po) / (n * (1 - pe) ** 2))
    z = kappa / se if se > 0 else float("inf")
    p = 2.0 * (1.0 - st.norm.cdf(abs(z)))
    return {"kappa": float(kappa), "p": float(p), "n": n}


def percent_agreement(a, b):
    aa, bb = _paired(a, b)
    if not aa:
        return None
    return sum(1 for x, y in zip(aa, bb) if x == y) / len(aa)


def mannwhitney(x, y):
    xs = [v for v in x if v is not None]
    ys = [v for v in y if v is not None]
    if len(xs) < 3 or len(ys) < 3:
        return None
    u, p = st.mannwhitneyu(xs, ys, alternative="two-sided")
    return {"u": float(u), "p": float(p), "n1": len(xs), "n2": len(ys)}


def _paired(a, b):
    aa, bb = [], []
    for x, y in zip(a, b):
        if x is not None and y is not None:
            aa.append(x)
            bb.append(y)
    return aa, bb
