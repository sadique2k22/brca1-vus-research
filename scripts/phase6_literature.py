"""Phase 6 — literature and experimental evidence review (evidence collection only).

Computes the candidate union from the FROZEN annotated dataset (re-deriving Phase 5
pattern classes deterministically), prioritizes into tiers (frozen in
configs/analysis_config.yaml), then retrieves: Findlay 2018 saturation-genome-editing
scores (MaveDB), current ClinVar status (E-utilities), and PubMed records (esearch),
all cached and rate-limited. No pathogenic/benign classification.
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

from src.literature import fetch_clinvar_status, fetch_findlay_scores, pubmed_search
from src.statistics import parse_float

ANN_TSV = "data/processed/brca1_vus_missense_annotated.tsv"
FIN_REPORT = "results/reports/annotation_finalization_report.md"
TABLES = "results/tables"
REPORTS = "results/reports"
CACHE = "data/intermediate/phase6_cache"

CLASS_DESC = {"A": "gnomAD-present + high impact", "B": "gnomAD-present + tolerance",
              "C": "gnomAD-absent + strong impact", "D": "gnomAD-absent + tolerance",
              "E": "strong predictor disagreement", "F": "elevated population frequency",
              "G": "extreme REVEL score"}


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
        idx = {n.lstrip("#"): i for i, n in enumerate(header)}
        rows = [row for row in r]
    return header, idx, rows


def classify_patterns(idx, rows):
    """Return dict variant_key -> set of pattern classes (deterministic, matches Phase 5)."""
    pat = defaultdict(set)
    for r in rows:
        vid_key = r[idx["variant_key"]]
        found = r[idx["gnomad_found"]]
        rv = parse_float(r[idx["revel_score"]])
        sv = parse_float(r[idx["sift_score"]])
        pv = parse_float(r[idx["polyphen_score"]])
        faf = parse_float(r[idx["gnomad_exome_faf95_popmax"]]) or parse_float(r[idx["gnomad_genome_faf95_popmax"]])
        if found == "present" and rv is not None and rv >= 0.644:
            pat[vid_key].add("A")
        if found == "present" and rv is not None and rv <= 0.290:
            pat[vid_key].add("B")
        if found == "absent" and rv is not None and rv >= 0.932:
            pat[vid_key].add("C")
        if found == "absent" and rv is not None and rv <= 0.290:
            pat[vid_key].add("D")
        if rv is not None and sv is not None and pv is not None:
            if (rv >= 0.644 and sv > 0.05 and pv <= 0.446) or (rv <= 0.290 and sv <= 0.05 and pv >= 0.446):
                pat[vid_key].add("E")
        if faf is not None and faf >= 0.001:
            pat[vid_key].add("F")
        if rv is not None and (rv >= 0.932 or rv <= 0.003):
            pat[vid_key].add("G")
    return pat


def cdna_from(name):
    m = re.search(r':(c\.[^)\s]+)', name or "")
    return m.group(1) if m else ""


def main():
    t0 = time.time()
    os.makedirs(CACHE, exist_ok=True)
    cfg = yaml.safe_load(open("configs/analysis_config.yaml"))
    tier_map = {}
    for tier, classes in (("1", cfg["phase6"]["tier1_classes"]),
                          ("2", cfg["phase6"]["tier2_classes"]),
                          ("3", cfg["phase6"]["tier3_classes"])):
        for c in classes:
            tier_map[c] = tier

    # STEP 0: verify frozen dataset unchanged
    m = re.search(r"`([0-9a-f]{64})`", open(FIN_REPORT).read())
    expected = m.group(1) if m else None
    actual = sha256(ANN_TSV)
    if expected and actual != expected:
        print("CHECKSUM MISMATCH — STOP")
        sys.exit(1)
    header, idx, rows = load_ann()

    # STEP 1: candidate union
    patterns = classify_patterns(idx, rows)
    info = {r[idx["variant_key"]]: r for r in rows}
    union_rows = []
    for r in rows:
        k = r[idx["variant_key"]]
        cls = patterns.get(k, set())
        if not cls:
            continue
        union_rows.append(r)
    n_union = len(union_rows)
    print(f"Candidate union: {n_union} unique variants")

    # STEP 2: prioritize
    tiers = {}
    for r in union_rows:
        k = r[idx["variant_key"]]
        cls = patterns[k]
        tier = min(tier_map[c] for c in cls) if cls else "3"
        tiers[k] = tier
    tier_counts = Counter(tiers.values())
    print(f"Tiers: {dict(tier_counts)}")
    prioritized = [r for r in union_rows if tiers[r[idx["variant_key"]]] in ("1", "2")]
    print(f"Prioritized (Tier 1+2): {len(prioritized)}")

    # write union table
    with open(os.path.join(TABLES, "phase6_candidate_union.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["variant_key", "VariationID", "AlleleID", "protein_change", "cDNA_change",
                    "candidate_classes", "REVEL", "SIFT", "PolyPhen", "gnomAD_status", "gnomAD_AF", "tier"])
        for r in sorted(union_rows, key=lambda x: tiers[x[idx["variant_key"]]]):
            k = r[idx["variant_key"]]
            w.writerow([k, r[idx["VariationID"]], r[idx["AlleleID"]], r[idx["protein_change"]],
                        cdna_from(r[idx["normalized_hgvs_c"]]), ",".join(sorted(patterns[k])),
                        r[idx["revel_score"]], r[idx["sift_score"]], r[idx["polyphen_score"]],
                        r[idx["gnomad_found"]],
                        r[idx["gnomad_exome_af"]] or r[idx["gnomad_genome_af"]] or "",
                        tiers[k]])

    # STEP 6/7: Findlay 2018 functional scores (all candidates)
    by_cdna, by_protein = fetch_findlay_scores(os.path.join(CACHE, "findlay_scores.csv"))
    def findlay(r):
        cd = cdna_from(r[idx["normalized_hgvs_c"]])
        if cd in by_cdna:
            return by_cdna[cd]
        pc = r[idx["protein_change"]]
        return by_protein.get(pc)

    # STEP 8: ClinVar current status (all candidates)
    vids = [r[idx["VariationID"]] for r in union_rows]
    clinvar = fetch_clinvar_status(vids, os.path.join(CACHE, "clinvar_status.json"))

    # STEP 3-5: PubMed search (prioritized Tier 1+2)
    evidence_rows = []
    search_date = time.strftime("%Y-%m-%d", time.gmtime())
    for r in union_rows:
        k = r[idx["variant_key"]]
        pc = r[idx["protein_change"]]
        cd = cdna_from(r[idx["normalized_hgvs_c"]])
        vid = r[idx["VariationID"]]
        # Findlay functional row (partial coverage: RING + BRCT domains)
        fs = findlay(r)
        evidence_rows.append([k, vid, pc, cd, search_date, "Findlay SGE (MaveDB)", "MaveDB",
                              "", "", "Findlay et al. 2018 (PMID 30209399)", "saturation genome editing",
                              "saturation genome editing (HDR)", pc, fs, "", "", fs, "HIGH",
                              "partial coverage: RING+BRCT domains only"])
        # ClinVar row
        cv = clinvar.get(vid, {})
        cur_sig = cv.get("significance")
        frozen_sig = r[idx["ClinicalSignificance"]]
        note = "" if cur_sig == frozen_sig else "Current external evidence differs from the frozen ClinVar classification."
        evidence_rows.append([k, vid, pc, cd, search_date, "ClinVar esummary", "ClinVar",
                              "", "", "", "clinical interpretation", "NA", "NA", "NA", "",
                              cur_sig, f"current={cur_sig}; review={cv.get('review_status')}",
                              "HIGH" if cv.get("review_status") and "expert" in str(cv.get("review_status")) else "MODERATE", note])
        # PubMed rows (only prioritized)
        if tiers[k] in ("1", "2"):
            for q in (f"BRCA1 {pc}", f"BRCA1 {cd}"):
                qkey = hashlib.sha256(q.encode()).hexdigest()[:16]
                hits = pubmed_search(q, os.path.join(CACHE, f"pubmed_{qkey}.json"))
                for h in hits:
                    if h.get("error"):
                        evidence_rows.append([k, vid, pc, cd, search_date, q, "PubMed",
                                              "", "", "", "NA", "NA", "NA", "NA", "", "",
                                              "NA", "LOW", "search failed"])
                    elif h.get("pmid"):
                        evidence_rows.append([k, vid, pc, cd, search_date, q, "PubMed",
                                              h.get("pmid"), h.get("doi"), h.get("title"),
                                              "publication", "NA", "NA", "NA", "", "",
                                              "NA", "UNCLEAR", ""])
                else:
                    continue
                if not hits or (len(hits) == 1 and hits[0].get("error")):
                    continue
            # if no PubMed hits at all, record "no relevant publication identified"
            # (handled below in summary)
    write_evidence_table(evidence_rows)

    # STEP 12: conflicts
    conflicts = detect_conflicts(union_rows, idx, patterns, findlay, clinvar)
    write_conflicts(conflicts)

    # STEP 13: representative variants
    write_representatives(union_rows, idx, patterns, findlay, clinvar)

    # summarize
    reviewed = n_union
    with_findlay = sum(1 for r in union_rows if findlay(r) is not None)
    with_clinvar = sum(1 for r in union_rows if clinvar.get(r[idx["VariationID"]], {}).get("significance"))
    # PubMed "no relevant literature" count for prioritized
    no_lit = 0
    for r in prioritized:
        k = r[idx["variant_key"]]
        pc = r[idx["protein_change"]]
        cd = cdna_from(r[idx["normalized_hgvs_c"]])
        any_hit = False
        for q in (f"BRCA1 {pc}", f"BRCA1 {cd}"):
            qkey = hashlib.sha256(q.encode()).hexdigest()[:16]
            hits = pubmed_search(q, os.path.join(CACHE, f"pubmed_{qkey}.json"))
            if any(h.get("pmid") for h in hits):
                any_hit = True
                break
        if not any_hit:
            no_lit += 1

    write_report(n_union, tier_counts, reviewed, with_findlay, with_clinvar, no_lit, len(prioritized), conflicts, t0)


def write_evidence_table(rows):
    cols = ["variant_key", "VariationID", "protein_change", "cDNA_change", "search_date",
            "search_query", "source", "PMID", "DOI", "paper_title", "study_type",
            "functional_assay", "exact_variant_tested", "functional_result", "clinical_cases",
            "expert_curation", "evidence_summary", "evidence_quality", "notes"]
    with open(os.path.join(TABLES, "variant_literature_evidence.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for row in rows:
            w.writerow(row)


def detect_conflicts(union_rows, idx, patterns, findlay, clinvar):
    conflicts = []
    for r in union_rows:
        k = r[idx["variant_key"]]
        rv = parse_float(r[idx["revel_score"]])
        fs = findlay(r)
        vid = r[idx["VariationID"]]
        cv = clinvar.get(vid, {})
        cur = cv.get("significance")
        frozen = r[idx["ClinicalSignificance"]]
        reasons = []
        if rv is not None and fs is not None:
            if rv >= 0.644 and fs > 0:
                reasons.append("computational impact vs functional-normal (Findlay score > 0)")
            if rv <= 0.290 and fs < 0:
                reasons.append("computational tolerance vs functional-LOF (Findlay score < 0)")
        if cur and cur != frozen and cur != "Uncertain significance":
            reasons.append(f"ClinVar changed from frozen '{frozen}' to '{cur}'")
        if reasons:
            conflicts.append([k, r[idx["VariationID"]], r[idx["protein_change"]],
                              rv, fs, "; ".join(reasons)])
    return conflicts


def write_conflicts(conflicts):
    with open(os.path.join(TABLES, "phase6_evidence_conflicts.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["variant_key", "VariationID", "protein_change", "REVEL", "Findlay_score", "conflict"])
        w.writerows(conflicts)


def write_representatives(union_rows, idx, patterns, findlay, clinvar):
    cats = {"damaging_supported": [], "damaging_contradicted": [], "tolerated_evidence": [],
            "disagreement": [], "insufficient": []}
    for r in union_rows:
        k = r[idx["variant_key"]]
        rv = parse_float(r[idx["revel_score"]])
        fs = findlay(r)
        if rv is not None and rv >= 0.644 and fs is not None and fs < 0:
            cats["damaging_supported"].append(r)
        elif rv is not None and rv >= 0.644 and fs is not None and fs > 0:
            cats["damaging_contradicted"].append(r)
        elif rv is not None and rv <= 0.290 and fs is not None:
            cats["tolerated_evidence"].append(r)
        elif "E" in patterns[k]:
            cats["disagreement"].append(r)
        elif fs is None:
            cats["insufficient"].append(r)
    with open(os.path.join(TABLES, "representative_variants.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["category", "variant_key", "VariationID", "protein_change", "REVEL", "Findlay_score"])
        for cat, rr in cats.items():
            for r in rr[:3]:
                w.writerow([cat, r[idx["variant_key"]], r[idx["VariationID"]], r[idx["protein_change"]],
                            r[idx["revel_score"]], findlay(r)])


def write_report(n_union, tier_counts, reviewed, with_findlay, with_clinvar, no_lit, n_prioritized, conflicts, t0):
    L = ["# Phase 6 — Literature & Experimental Evidence Review\n\n"]
    L.append("Evidence-review only. No pathogenic/benign classification, no ACMG codes.\n\n")
    L.append("## 1. Candidate selection\n\n")
    L.append(f"- Unique candidate union: **{n_union}**\n")
    L.append(f"- Tiers: {dict(tier_counts)}\n")
    L.append(f"- Prioritized for literature search (Tier 1+2): {n_prioritized}\n")
    L.append("Methodology: `results/reports/phase6_selection_methodology.md`\n\n")
    L.append("## 2. Evidence coverage\n\n")
    L.append(f"- Variants reviewed: {reviewed}\n")
    L.append(f"- With Findlay 2018 functional score: {with_findlay}\n")
    L.append(f"- With current ClinVar curation: {with_clinvar}\n")
    L.append(f"- Prioritized variants with no PubMed hit identified: {no_lit}\n\n")
    L.append("## 3. Major evidence conflicts\n\n")
    L.append(f"- Conflicts identified: **{len(conflicts)}** "
             "(see `results/tables/phase6_evidence_conflicts.tsv`)\n\n")
    L.append("## 4. Representative variants\n\n")
    L.append("See `results/tables/representative_variants.tsv` (selected only after systematic collection).\n\n")
    L.append("## 5. Limitations\n\n")
    L.append("- PubMed search is title/metadata level (esearch/esummary), not full-text.\n")
    L.append("- 'No PubMed hit' means 'No relevant publication identified using the documented "
             "search strategy' — NOT 'no functional evidence exists'.\n")
    L.append("- Findlay 2018 functional score retrieved from MaveDB covers only the RING and "
             "BRCT domains (partial dataset); whole-gene scores require the Nature supplementary "
             "table (recommended follow-up).\n")
    L.append("- Literature-derived evidence does not modify the frozen computational dataset.\n\n")
    L.append("## 6. Recommendations\n\n")
    L.append("Focus further investigation on conflicts and on absent+strong-impact variants. "
             "No variant is claimed pathogenic or benign.\n")
    with open(os.path.join(REPORTS, "phase6_literature_review.md"), "w") as fh:
        fh.write("".join(L))
    print(f"Report written ({time.time()-t0:.1f}s). Union={n_union}, findlay={with_findlay}, "
          f"clinvar={with_clinvar}, no_pubmed={no_lit}, conflicts={len(conflicts)}")


if __name__ == "__main__":
    main()
