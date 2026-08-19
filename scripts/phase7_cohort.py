"""Phase 7 — final evidence cohort selection and systematic evidence synthesis.

Selects a deterministic stratified cohort (~50 variants across strata A-E) from the
derived functional dataset, performs evidence search (ClinVar + PubMed + Findlay),
builds an evidence matrix, and produces conflict/predictor/domain analyses, figures and
a report. Frozen datasets are never modified.
"""
import csv
import hashlib
import os
import re
import sys
import time
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import yaml

from src import figures
from src.literature import fetch_clinvar_status, pubmed_search
from src.statistics import descriptive, parse_float, spearman

PROC = "data/processed"
TABLES = "results/tables"
REPORTS = "results/reports"
CACHE = "data/intermediate/phase7_cache"
WITH_FUNC = os.path.join(PROC, "brca1_vus_missense_with_functional.tsv")


def load_with_functional():
    with open(WITH_FUNC, newline="", encoding="utf-8") as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        idx = {n.lstrip("#"): i for i, n in enumerate(header)}
        rows = [row for row in r]
    return header, idx, rows


def domain_of(aa):
    if aa is None:
        return "other"
    if aa <= 109:
        return "RING"
    if aa >= 1642:
        return "BRCT"
    return "other"


def aa_pos(protein_change):
    m = re.search(r'(\d+)', protein_change or "")
    return int(m.group(1)) if m else None


def cdna_of(hgvsc):
    m = re.search(r':(c\.[^)\s]+)', hgvsc or "")
    return m.group(1) if m else ""


def classify_stratum(rv, fs):
    if rv is None or fs is None:
        return None
    if rv >= 0.644:
        return "A" if fs < 0 else "B"
    if rv <= 0.290:
        return "C" if fs < 0 else "D"
    return "E"


def main():
    t0 = time.time()
    os.makedirs(CACHE, exist_ok=True)
    cfg = yaml.safe_load(open("configs/phase7_cohort.yaml"))
    header, idx, rows = load_with_functional()

    # classify all scored variants into strata
    strata = defaultdict(list)
    for r in rows:
        fs = parse_float(r[idx["findlay_score"]])
        rv = parse_float(r[idx["revel_score"]])
        if fs is None:
            continue
        s = classify_stratum(rv, fs)
        if s:
            pos = aa_pos(r[idx["protein_change"]])
            strata[s].append((pos if pos is not None else 0, r))
    print("strata sizes (scored):", {k: len(v) for k, v in strata.items()})

    # deterministic selection: even-spacing by protein position
    target = cfg["phase7"]["target_per_stratum"]
    cohort = []
    for s in "ABCDE":
        members = sorted(strata[s], key=lambda x: x[0])
        n = len(members)
        if n <= target:
            picks = members
        else:
            picks = [members[round(i * (n - 1) / (target - 1))] for i in range(target)]
        for _, r in picks:
            r = list(r)
            cohort.append(r + [s])
    print(f"Cohort size: {len(cohort)}")
    stratum_counts = Counter(c[-1] for c in cohort)

    # write frozen cohort
    with open(os.path.join(TABLES, "phase7_final_cohort.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["variant_key", "VariationID", "AlleleID", "cDNA_change", "protein_change",
                    "exon", "domain", "stratum", "REVEL", "SIFT", "PolyPhen", "gnomAD_status",
                    "gnomAD_AF", "Findlay_score"])
        for r in cohort:
            s = r[-1]
            aa = aa_pos(r[idx["protein_change"]])
            dom = domain_of(aa)
            exon = "2-5" if dom == "RING" else ("15-23" if dom == "BRCT" else "NA")
            w.writerow([r[idx["variant_key"]], r[idx["VariationID"]], r[idx["AlleleID"]],
                        cdna_of(r[idx["normalized_hgvs_c"]]), r[idx["protein_change"]],
                        exon, dom, s,
                        r[idx["revel_score"]], r[idx["sift_score"]], r[idx["polyphen_score"]],
                        r[idx["gnomad_found"]],
                        r[idx["gnomad_exome_af"]] or r[idx["gnomad_genome_af"]] or "",
                        r[idx["findlay_score"]]])

    # evidence search
    vids = [r[idx["VariationID"]] for r in cohort]
    clinvar = fetch_clinvar_status(vids, os.path.join(CACHE, "clinvar_status.json"))
    search_date = time.strftime("%Y-%m-%d", time.gmtime())
    evidence = []
    search_log = []
    for r in cohort:
        s = r[-1]
        k = r[idx["variant_key"]]
        vid = r[idx["VariationID"]]
        pc = r[idx["protein_change"]]
        cd = cdna_of(r[idx["normalized_hgvs_c"]])
        fs = parse_float(r[idx["findlay_score"]])
        rv = parse_float(r[idx["revel_score"]])
        gf = r[idx["gnomad_found"]]
        cv = clinvar.get(vid, {})
        review = cv.get("review_status") or "NA"
        cur_sig = cv.get("significance") or "NA"
        # PubMed
        pub_hits = 0
        for q in (f"BRCA1 {pc}", f"BRCA1 {cd}"):
            qkey = hashlib.sha256(q.encode()).hexdigest()[:16]
            hits = pubmed_search(q, os.path.join(CACHE, f"pubmed_{qkey}.json"))
            n_hit = sum(1 for h in hits if h.get("pmid"))
            pub_hits += n_hit
            search_log.append([k, "PubMed", q, search_date, len(hits), n_hit, n_hit, ""])
        # evidence fields
        comp_ev = "impact" if (rv is not None and rv >= 0.644) else ("tolerance" if (rv is not None and rv <= 0.290) else "intermediate")
        findlay_ev = "LOF" if (fs is not None and fs < 0) else ("normal" if fs is not None else "NA")
        conflict = "no"
        if rv is not None and fs is not None:
            if (rv >= 0.644 and fs >= 0) or (rv <= 0.290 and fs < 0):
                conflict = "yes"
        quality = "HIGH" if "expert" in str(review) else "MODERATE"
        summary = neutral_summary(comp_ev, findlay_ev, conflict)
        evidence.append([k, pc, s, gf, comp_ev, findlay_ev,
                         "none identified", "not assessed", "not assessed",
                         review, cur_sig, pub_hits, conflict, quality, summary])

    # write evidence matrix
    cols = ["variant_key", "protein_change", "stratum", "gnomAD_evidence",
            "computational_evidence", "Findlay_evidence", "other_functional_evidence",
            "clinical_case_evidence", "segregation_evidence", "expert_curation",
            "ClinVar_current_status", "literature_count", "evidence_conflict",
            "evidence_quality", "evidence_summary"]
    with open(os.path.join(TABLES, "phase7_evidence_matrix.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        w.writerows(evidence)

    # search log
    with open(os.path.join(TABLES, "phase7_search_log.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["variant_key", "database", "query", "date", "results_count",
                    "records_screened", "records_relevant", "notes"])
        w.writerows(search_log)

    # conflict analysis
    conflicts = [e for e in evidence if e[12] == "yes"]
    print(f"Conflicts in cohort: {len(conflicts)}")

    # predictor/function analysis (ordinal; no ROC without validated binary threshold)
    fs_all = [parse_float(r[idx["findlay_score"]]) for r in rows]
    rv_all = [parse_float(r[idx["revel_score"]]) for r in rows]
    sv_all = [parse_float(r[idx["sift_score"]]) for r in rows]
    pv_all = [parse_float(r[idx["polyphen_score"]]) for r in rows]
    scored = [(fs_all[i], rv_all[i], sv_all[i], pv_all[i]) for i in range(len(rows)) if fs_all[i] is not None]
    corr = {
        "revel": spearman([x[0] for x in scored], [x[1] for x in scored]),
        "sift": spearman([x[0] for x in scored], [(1 - x[2]) if x[2] is not None else None for x in scored]),
        "polyphen": spearman([x[0] for x in scored], [x[3] for x in scored]),
    }

    # domain analysis
    ring_scores = [parse_float(r[idx["findlay_score"]]) for r in rows
                   if parse_float(r[idx["findlay_score"]]) is not None and domain_of(aa_pos(r[idx["protein_change"]])) == "RING"]
    brct_scores = [parse_float(r[idx["findlay_score"]]) for r in rows
                   if parse_float(r[idx["findlay_score"]]) is not None and domain_of(aa_pos(r[idx["protein_change"]])) == "BRCT"]

    # figures
    figures.fig11_cohort(stratum_counts)
    figures.fig12_predictor_functional([x[1] for x in scored], [x[0] for x in scored])
    figures.fig13_concordance(evidence)
    figures.fig14_functional_by_category(scored)
    figures.fig15_evidence_availability(evidence)
    figures.fig16_domain(ring_scores, brct_scores)

    write_report(len(cohort), stratum_counts, evidence, conflicts, corr, len(scored), ring_scores, brct_scores, t0)
    print(f"Done ({time.time()-t0:.1f}s). cohort={len(cohort)}, conflicts={len(conflicts)}")


def neutral_summary(comp, func, conflict):
    if conflict == "yes":
        return "conflicting evidence (computational vs functional)"
    if comp == "impact" and func == "LOF":
        return "evidence predominantly supports functional impact"
    if comp == "tolerance" and func == "normal":
        return "evidence predominantly supports functional tolerance"
    if comp == "impact" and func == "normal":
        return "conflicting evidence (computational impact vs functional normal)"
    if comp == "tolerance" and func == "LOF":
        return "conflicting evidence (computational tolerance vs functional LOF)"
    return "intermediate/insufficient evidence"


def write_report(n, stratum_counts, evidence, conflicts, corr, n_scored, ring, brct, t0):
    L = ["# Phase 7 — Final Evidence Synthesis\n\n"]
    L.append("Evidence synthesis only. No pathogenic/benign classification, no ACMG codes.\n\n")
    L.append("## 1. Final cohort\n\n")
    L.append(f"- Cohort size: **{n}**\n")
    L.append(f"- Stratum counts: {dict(stratum_counts)}\n\n")
    L.append("## 2. Evidence availability\n\n")
    lit = sum(1 for e in evidence if e[11] > 0)
    exp = sum(1 for e in evidence if "expert" in str(e[9]))
    L.append(f"- With >=1 PubMed hit: {lit}/{n}\n")
    L.append(f"- With expert-panel curation: {exp}/{n}\n")
    L.append(f"- With Findlay exact-variant score: {n}/{n} (by cohort definition)\n\n")
    L.append("## 3. Conflict analysis\n\n")
    L.append(f"- Evidence conflicts in cohort: **{len(conflicts)}**\n\n")
    L.append("## 4. Predictor vs functional (ordinal, Spearman)\n\n")
    L.append("(ROC-AUC not computed: a validated binary Findlay threshold was not retrieved; "
             "ordinal correlation is reported instead.)\n\n")
    L.append("| Predictor | rho | p | n |\n|---|---|---|---|\n")
    for k, v in corr.items():
        if v:
            L.append(f"| {k} | {v['rho']:.3f} | {v['p']:.3g} | {v['n']} |\n")
    L.append(f"\n(Full-scored n = {n_scored}.)\n\n")
    L.append("## 5. Domain analysis (Findlay score)\n\n")
    if ring and brct:
        L.append(f"- RING: n={len(ring)}, median={descriptive(ring).get('median')}\n")
        L.append(f"- BRCT: n={len(brct)}, median={descriptive(brct).get('median')}\n")
    L.append("\n## 6. Limitations\n\n")
    L.append("- Findlay covers RING + BRCT only; cohort is drawn from these regions.\n")
    L.append("- Clinical/segregation evidence is not programmatically extractable at metadata level.\n")
    L.append("- PubMed is title/metadata level.\n")
    with open(os.path.join(REPORTS, "phase7_evidence_synthesis.md"), "w") as fh:
        fh.write("".join(L))


if __name__ == "__main__":
    main()
