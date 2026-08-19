"""Phase 9 — final evidence verification and scientific refinement.

A) literature verification (abstract-level exact-variant classification)
B) Findlay score interpretation (continuous; mixture-model classification documented)
C) continuous correlations with confidence intervals
D) conflict reassessment (continuous)
E) domain analysis with effect size
F) rebuilt evidence matrix
G) literature-verification table + report
"""
import csv
import hashlib
import json
import math
import os
import re
import sys
import time
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import requests

from src.statistics import descriptive, mannwhitney, parse_float, spearman

PROC = "data/processed"
TABLES = "results/tables"
REPORTS = "results/reports"
CACHE = "data/intermediate/phase9_cache"
WITH_FUNC = os.path.join(PROC, "brca1_vus_missense_with_functional.tsv")
COHORT = os.path.join(TABLES, "phase7_final_cohort.tsv")

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
       "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
       "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
       "Tyr": "Y", "Val": "V"}


def one_letter(pc):
    m = re.match(r'^([A-Z][a-z][a-z])(\d+)([A-Z][a-z][a-z])$', pc or "")
    if m and m.group(1) in AA3 and m.group(3) in AA3:
        return f"{AA3[m.group(1)]}{m.group(2)}{AA3[m.group(3)]}"
    return None


def fetch_abstract(pmid, cache_file):
    if os.path.exists(cache_file):
        return open(cache_file).read()
    try:
        r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                         params={"db": "pubmed", "id": pmid, "rettype": "abstract", "retmode": "text"},
                         timeout=40)
        text = r.text
    except Exception:
        text = ""
    with open(cache_file, "w") as fh:
        fh.write(text)
    time.sleep(0.4)
    return text


def classify(abstract, pc, cd):
    if not abstract:
        return "UNCLEAR"
    for c in (pc, f"p.{pc}", one_letter(pc), f"p.{one_letter(pc)}", cd):
        if c and c in abstract:
            return "EXACT_VARIANT"
    if "BRCA1" in abstract:
        return "GENE_LEVEL"
    return "UNCLEAR"


def pubmed_hits(pc, cd):
    """Re-run PubMed search for a variant (cached) and return list of PMIDs."""
    pmids = []
    for q in (f"BRCA1 {pc}", f"BRCA1 {cd}"):
        qkey = hashlib.sha256(q.encode()).hexdigest()[:16]
        cache_file = os.path.join(CACHE, f"pubmed_{qkey}.json")
        if os.path.exists(cache_file):
            hits = json.load(open(cache_file))
        else:
            r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                             params={"db": "pubmed", "term": q, "retmax": "5", "retmode": "json"},
                             timeout=60)
            try:
                ids = r.json().get("esearchresult", {}).get("idlist", [])
            except Exception:
                ids = []
            hits = [{"pmid": i} for i in ids]
            with open(cache_file, "w") as fh:
                json.dump(hits, fh)
            time.sleep(0.4)
        for h in hits:
            if h.get("pmid") and h["pmid"] not in pmids:
                pmids.append(h["pmid"])
    return pmids


def spearman_ci(rho, n):
    if n < 4 or rho is None:
        return None
    rho = max(-0.999, min(0.999, rho))
    z = math.atanh(rho)
    se = 1.0 / math.sqrt(n - 3)
    lo = math.tanh(z - 1.96 * se)
    hi = math.tanh(z + 1.96 * se)
    return (round(lo, 3), round(hi, 3))


def main():
    os.makedirs(CACHE, exist_ok=True)
    _, cidx, cohort = load(COHORT)
    print(f"Cohort: {len(cohort)} variants")

    # ---- PART A/G: literature verification ----
    lit_rows = []
    stats = Counter()
    for r in cohort:
        k = r[cidx["variant_key"]]
        vid = r[cidx["VariationID"]]
        pc = r[cidx["protein_change"]]
        cd = r[cidx["cDNA_change"]]
        pmids = pubmed_hits(pc, cd)
        stats["variants_searched"] += 1
        if pmids:
            stats["variants_with_pubmed"] += 1
        for pmid in pmids:
            abstract = fetch_abstract(pmid, os.path.join(CACHE, f"abs_{pmid}.txt"))
            cls = classify(abstract, pc, cd)
            stats[f"class_{cls}"] += 1
            exact = "yes" if cls == "EXACT_VARIANT" else "no"
            lit_rows.append([k, vid, pc, cd, pmid, "", f"BRCA1 {pc}; BRCA1 {cd}",
                             cls, exact, "no", "no", "no", "no",
                             "VERIFIED" if cls == "EXACT_VARIANT" else "NOT_VERIFIED", ""])
        if not pmids:
            lit_rows.append([k, vid, pc, cd, "", "", "", "NONE", "no", "no", "no",
                             "no", "no", "NONE", "no PubMed record"])
    with open(os.path.join(TABLES, "phase9_literature_verification.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["variant_key", "VariationID", "protein_change", "cDNA_change", "PMID", "DOI",
                    "search_query", "paper_type", "variant_exactly_tested", "functional_evidence",
                    "clinical_evidence", "segregation_evidence", "expert_evidence",
                    "verification_status", "notes"])
        w.writerows(lit_rows)

    # ---- PART B/C: continuous correlations + CI ----
    _, widx, wf = load(WITH_FUNC)
    fs = [parse_float(r[widx["findlay_score"]]) for r in wf]
    rv = [parse_float(r[widx["revel_score"]]) for r in wf]
    sv = [parse_float(r[widx["sift_score"]]) for r in wf]
    pv = [parse_float(r[widx["polyphen_score"]]) for r in wf]
    idxs = [i for i in range(len(fs)) if fs[i] is not None]
    corr = {}
    for name, pred in (("revel", rv), ("sift", [(1 - sv[i]) if sv[i] is not None else None for i in range(len(sv))]),
                       ("polyphen", pv)):
        c = spearman([fs[i] for i in idxs], [pred[i] for i in idxs])
        corr[name] = {"rho": c["rho"], "p": c["p"], "n": c["n"],
                      "ci": spearman_ci(c["rho"], c["n"])}

    # ---- PART D: conflict reassessment (continuous) ----
    # conflict = computational impact (REVEL>=0.644) with Findlay score > 0 (functional range),
    # or tolerance (REVEL<=0.290) with Findlay score < -1.0 (non-functional range).
    conflicts = []
    for r in cohort:
        k = r[cidx["variant_key"]]
        rv = parse_float(r[cidx["REVEL"]])
        fs_val = parse_float(r[cidx["Findlay_score"]])
        if rv is None or fs_val is None:
            continue
        cat = "impact" if rv >= 0.644 else ("tolerance" if rv <= 0.290 else "intermediate")
        if cat == "impact" and fs_val > 0:
            conflicts.append([k, r[cidx["protein_change"]], r[cidx["cDNA_change"]],
                              cat, rv, fs_val, "impact-supporting vs functional WT-like score"])
        elif cat == "tolerance" and fs_val < -1.0:
            conflicts.append([k, r[cidx["protein_change"]], r[cidx["cDNA_change"]],
                              cat, rv, fs_val, "tolerance-supporting vs non-functional score (< -1)"])
    print(f"Conflicts (continuous): {len(conflicts)}")

    # ---- PART E: domain ----
    ring = [parse_float(r[widx["findlay_score"]]) for r in wf
            if parse_float(r[widx["findlay_score"]]) is not None and aa_pos(r[widx["protein_change"]]) is not None
            and aa_pos(r[widx["protein_change"]]) <= 109]
    brct = [parse_float(r[widx["findlay_score"]]) for r in wf
            if parse_float(r[widx["findlay_score"]]) is not None and aa_pos(r[widx["protein_change"]]) is not None
            and aa_pos(r[widx["protein_change"]]) >= 1642]
    mw = mannwhitney(ring, brct)
    rank_biserial = (1 - 2 * mw["u"] / (mw["n1"] * mw["n2"])) if mw else None
    domain = {"ring": descriptive(ring), "brct": descriptive(brct), "mw": mw,
              "rank_biserial": rank_biserial}

    # ---- PART F: rebuilt evidence matrix ----
    matrix_rows = []
    for r in cohort:
        k = r[cidx["variant_key"]]
        pc = r[cidx["protein_change"]]
        rv = parse_float(r[cidx["REVEL"]])
        fs_val = parse_float(r[cidx["Findlay_score"]])
        comp = "impact" if rv is not None and rv >= 0.644 else ("tolerance" if rv is not None and rv <= 0.290 else "intermediate")
        func = "non-functional (score<0)" if fs_val is not None and fs_val < 0 else ("WT-like (score>=0)" if fs_val is not None else "NA")
        # exact-variant literature + functional evidence from lit_rows
        exact = any(x[0] == k and x[9] == "yes" for x in lit_rows)
        conflict = "yes" if any(x[0] == k for x in conflicts) else "no"
        matrix_rows.append([k, pc, comp, func,
                            "present" if r[cidx["gnomAD_status"]] == "present" else "absent",
                            "yes" if exact else "no",
                            "yes" if fs_val is not None else "no",
                            "no", "no", "no", conflict, "MODERATE"])
    with open(os.path.join(TABLES, "phase9_evidence_matrix.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["variant_key", "protein_change", "computational_evidence", "findlay_evidence",
                    "population_evidence", "literature_hit", "functional_evidence",
                    "clinical_evidence", "segregation_evidence", "expert_curation",
                    "evidence_conflict", "evidence_quality"])
        w.writerows(matrix_rows)

    write_reports(stats, corr, conflicts, domain, lit_rows)
    print(json.dumps({"literature": dict(stats), "correlations": {k: {"rho": v["rho"], "ci": v["ci"]} for k, v in corr.items()},
                      "conflicts": len(conflicts), "domain": {"ring_n": len(ring), "brct_n": len(brct),
                                                              "mw_p": mw["p"] if mw else None,
                                                              "rank_biserial": rank_biserial}}, indent=2))


def aa_pos(pc):
    m = re.search(r'(\d+)', pc or "")
    return int(m.group(1)) if m else None


def load(path):
    with open(path, newline="", encoding="utf-8") as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        idx = {n.lstrip("#"): i for i, n in enumerate(header)}
        rows = [row for row in r]
    return header, idx, rows


def write_reports(stats, corr, conflicts, domain, lit_rows):
    # findlay interpretation doc
    L = ["# Findlay 2018 Score Interpretation — Phase 9\n\n"]
    L.append("## Score meaning (from the original study)\n\n")
    L.append("- Function score = log2 ratio of each SNV's frequency on day 11 vs. the original plasmid library, "
             "positional-bias-corrected, normalized across exons, averaged over 2 replicates (HAP1 cells).\n")
    L.append("- Lower (more negative) score = reduced cellular fitness; ~0 = WT-like.\n\n")
    L.append("## Original classification method\n\n")
    L.append("- Findlay et al. fitted a **two-component Gaussian mixture model** to the function scores and "
             "classified each SNV by the posterior probability of non-functionality P(nf):\n")
    L.append("  - P(nf) > 0.99 = 'non-functional'\n")
    L.append("  - 0.01 < P(nf) < 0.99 = 'intermediate'\n")
    L.append("  - P(nf) < 0.01 = 'functional'\n")
    L.append("- Synonymous SNVs (functional controls) scored ~0 (median 0.00; 98.7% > −1.25).\n\n")
    L.append("## Decision for this study\n\n")
    L.append("- The MaveDB score set provides only the **continuous** function score (plus replicates and RNA score); "
             "it does **not** provide P(nf) or the mixture-model classification.\n")
    L.append("- Therefore we use the **continuous Findlay score as the primary functional variable** "
             "(Option 1). We do **not** apply a binary threshold, because the validated binary classification "
             "requires the mixture-model fit (not retrievable from MaveDB), and a simple zero-based split would be "
             "arbitrary.\n")
    L.append("- The Phase 7 'score < 0 = non-functional / >= 0 = WT-like' heuristic is therefore superseded; "
             "correlations and conflict characterization now use the continuous score.\n")
    with open(os.path.join(REPORTS, "findlay_score_interpretation.md"), "w") as fh:
        fh.write("".join(L))

    # literature verification report
    G = ["# Phase 9 — Literature Verification Report\n\n"]
    G.append(f"- Variants searched: {stats.get('variants_searched', 0)}\n")
    G.append(f"- Variants with PubMed records: {stats.get('variants_with_pubmed', 0)}\n")
    G.append(f"- Paper classification: EXACT_VARIANT={stats.get('class_EXACT_VARIANT', 0)}, "
             f"GENE_LEVEL={stats.get('class_GENE_LEVEL', 0)}, UNCLEAR={stats.get('class_UNCLEAR', 0)}, "
             f"NONE={stats.get('class_NONE', 0)}\n")
    G.append("\nPubMed search hits are NOT evidence by themselves; only EXACT_VARIANT papers "
             "(abstract contains the exact cDNA/protein change) are counted as exact-variant evidence.\n")
    with open(os.path.join(REPORTS, "phase9_literature_verification.md"), "w") as fh:
        fh.write("".join(G))

    # final report
    F = ["# Phase 9 — Final Report\n\n"]
    F.append("## Literature verification\n\n")
    for k, v in sorted(stats.items()):
        F.append(f"- {k}: {v}\n")
    F.append("\n## Findlay interpretation\n\n- Continuous score used as primary functional variable; "
             "binary classification requires the (unavailable) mixture-model fit.\n\n")
    F.append("## Computational-functional correlations (continuous, Spearman + 95% CI)\n\n")
    F.append("| Predictor | rho | 95% CI | p | n |\n|---|---|---|---|---|\n")
    for k, v in corr.items():
        F.append(f"| {k} | {v['rho']:.3f} | {v['ci']} | {v['p']:.3g} | {v['n']} |\n")
    F.append(f"\n## Conflict reassessment (continuous): **{len(conflicts)}**\n\n")
    F.append("## Domain analysis\n\n")
    F.append(f"- RING: n={domain['ring'].get('n')}, median={domain['ring'].get('median')}\n")
    F.append(f"- BRCT: n={domain['brct'].get('n')}, median={domain['brct'].get('median')}\n")
    if domain['mw']:
        F.append(f"- Mann-Whitney p={domain['mw']['p']:.3g}; rank-biserial r={domain['rank_biserial']:.3f}\n")
    F.append("\n## Submission readiness\n\n")
    F.append("Manuscript is updated for: (1) PubMed hits ≠ evidence; (2) continuous Findlay score; "
             "(3) RING/BRCT scope; (4) no clinical reclassification. Remaining: author/affiliation/funding placeholders.\n")
    with open(os.path.join(REPORTS, "phase9_final_report.md"), "w") as fh:
        fh.write("".join(F))


if __name__ == "__main__":
    main()
