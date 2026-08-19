"""Phase 6 literature/evidence clients: MaveDB (Findlay 2018), ClinVar, PubMed.

All results are returned exactly as retrieved (no fabrication); network results are
cached to disk and rate-limit-aware.
"""
import csv
import io
import json
import os
import re
import time

import requests

FINDLY_URN = "urn:mavedb:00000097-0-2"
MAVEDB_URL = "https://api.mavedb.org/api/v1/score-sets/%s/scores" % FINDLY_URN
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def fetch_findlay_scores(cache_file, pause=1.0):
    """Download the Findlay 2018 BRCA1 SGE normalized scores (CSV) -> lookup dicts.

    Returns (by_cdna, by_protein) where by_cdna maps 'c.4219C>G' -> score and
    by_protein maps 'Leu1407Val' -> score.
    """
    if os.path.exists(cache_file):
        text = open(cache_file).read()
    else:
        r = requests.get(MAVEDB_URL, params={"limit": 100000}, timeout=300)
        r.raise_for_status()
        text = r.text
        with open(cache_file, "w") as fh:
            fh.write(text)
        time.sleep(pause)
    by_cdna, by_protein = {}, {}
    reader = csv.DictReader(io.StringIO(text))
    for row in reader:
        score = row.get("score")
        try:
            score = float(score)
        except (ValueError, TypeError):
            score = None
        nt = (row.get("hgvs_nt") or "")
        m = re.search(r':(c\.[^,]+)', nt)
        if m:
            by_cdna[m.group(1)] = score
        prot = (row.get("hgvs_pro") or "").lstrip("p.")
        if prot and prot != "NA":
            by_protein[prot] = score
    return by_cdna, by_protein


def fetch_clinvar_status(variation_ids, cache_file, chunk=200, pause=0.4):
    """Batch-fetch current ClinVar status for VariationIDs via esummary.

    Returns {vid: {significance, review_status, last_evaluated}}.
    """
    if os.path.exists(cache_file):
        return json.load(open(cache_file))
    out = {}
    ids = [v for v in variation_ids if v]
    for i in range(0, len(ids), chunk):
        batch = ids[i:i + chunk]
        r = requests.get(f"{EUTILS}/esummary.fcgi",
                         params={"db": "clinvar", "id": ",".join(batch), "retmode": "json"},
                         timeout=120)
        r.raise_for_status()
        res = r.json().get("result", {})
        for vid in batch:
            rec = res.get(vid, {})
            if isinstance(rec, dict):
                gc = rec.get("germline_classification") or {}
                out[vid] = {
                    "significance": gc.get("description") if isinstance(gc, dict) else None,
                    "review_status": gc.get("review_status") if isinstance(gc, dict) else None,
                    "last_evaluated": gc.get("last_evaluated") if isinstance(gc, dict) else None,
                }
            else:
                out[vid] = {"significance": None, "review_status": None, "last_evaluated": None}
        time.sleep(pause)
    with open(cache_file, "w") as fh:
        json.dump(out, fh)
    return out


def pubmed_search(term, cache_file, retmax=5, pause=0.4):
    """esearch PubMed for `term`, then esummary top hits. Returns list of {pmid, year, title, journal}."""
    if os.path.exists(cache_file):
        return json.load(open(cache_file))
    hits = []
    try:
        r = requests.get(f"{EUTILS}/esearch.fcgi",
                         params={"db": "pubmed", "term": term, "retmax": str(retmax), "retmode": "json"},
                         timeout=60)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        time.sleep(pause)
        if ids:
            r2 = requests.get(f"{EUTILS}/esummary.fcgi",
                              params={"db": "pubmed", "id": ",".join(ids), "retmode": "json"},
                              timeout=60)
            r2.raise_for_status()
            res = r2.json().get("result", {})
            for pid in ids:
                x = res.get(pid, {})
                if isinstance(x, dict):
                    hits.append({
                        "pmid": pid,
                        "year": (x.get("pubdate") or "")[:4],
                        "title": x.get("title", ""),
                        "journal": x.get("fulljournalname", ""),
                        "doi": None,
                    })
        time.sleep(pause)
    except Exception:
        hits = [{"pmid": None, "year": None, "title": None, "journal": None, "doi": None,
                 "error": "search_failed"}]
    with open(cache_file, "w") as fh:
        json.dump(hits, fh)
    return hits
