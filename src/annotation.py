"""Ensembl VEP REST client for transcript-consequence annotation (cached).

Used in Phase 4A to re-derive HGVS on the approved transcript from genomic coordinates
and to validate the reference allele against GRCh38 (two-pass validation).
"""
import hashlib
import json
import os
import time

import requests

VEP_URL = "https://rest.ensembl.org/vep/human/region"
MANE_ENST = "ENST00000357654"  # = NM_007294.4 (BRCA1 MANE Select)
CHUNK = 200


def _batch_key(batch):
    h = hashlib.sha256("\n".join(batch).encode()).hexdigest()[:16]
    return h


def vep_annotate(variant_strings, cache_dir, chunk=50, pause=0.5, retries=3):
    """Annotate variant strings via VEP REST; cache each batch's raw JSON.

    Returns a dict {input_string: result_dict}. Batches that fail after retries are
    skipped (their variants are absent from the result and reported as unresolved).
    """
    os.makedirs(cache_dir, exist_ok=True)
    out = {}
    failed = []
    for i in range(0, len(variant_strings), chunk):
        batch = variant_strings[i:i + chunk]
        key = _batch_key(batch)
        cache_file = os.path.join(cache_dir, f"batch_{key}.json")
        if os.path.exists(cache_file):
            data = json.load(open(cache_file))
        else:
            data = None
            for attempt in range(retries):
                try:
                    resp = requests.post(
                        VEP_URL + "?hgvs=1",
                        headers={"Content-Type": "application/json", "Accept": "application/json"},
                        data=json.dumps({"variants": batch}),
                        timeout=120,
                    )
                    resp.raise_for_status()
                    data = resp.json()
                    break
                except Exception as exc:  # noqa: BLE001
                    if attempt == retries - 1:
                        failed.append((i, str(exc)))
                    time.sleep(2 ** attempt)
            if data is None:
                continue
            with open(cache_file, "w") as fh:
                json.dump(data, fh)
            time.sleep(pause)
        for d in data:
            out[d.get("input", "")] = d
    return out, failed


def extract_mane_consequence(result):
    """From one VEP result dict, return (hgvsc, hgvsp, amino_acids, consequence, error).

    Picks the transcript consequence for the MANE Select transcript (ENST00000357654).
    """
    error = result.get("error")
    hgvsc = hgvsp = amino_acids = consequence = None
    for tc in result.get("transcript_consequences", []):
        if str(tc.get("transcript_id", "")).startswith(MANE_ENST):
            hgvsc = tc.get("hgvsc")
            hgvsp = tc.get("hgvsp")
            amino_acids = tc.get("amino_acids")
            consequence = ",".join(tc.get("consequence_terms", []))
            break
    return hgvsc, hgvsp, amino_acids, consequence, error
