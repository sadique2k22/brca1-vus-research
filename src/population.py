"""gnomAD population-frequency retrieval via the GraphQL API (cached, batched, resumable).

Uses the public gnomAD GraphQL endpoint (https://gnomad.broadinstitute.org/api), dataset
`gnomad_r4` (v4, GRCh38). One request per unique biological variant (batched via GraphQL
aliases); responses cached to disk keyed by variant key.
"""
import json
import os
import time

import requests

GNOMAD_URL = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r4"
BATCH = 50

_VARIANT_FIELDS = """
  variantId
  genome { ac an af homozygote_count faf95 { popmax popmax_population } populations { id ac an } }
  exome  { ac an af homozygote_count faf95 { popmax popmax_population } populations { id ac an } }
"""


def _pop_af(populations):
    """Compute per-population AF from ac/an; return {id: af} and popmax AF."""
    if not populations:
        return {}, None
    afs = {}
    for p in populations:
        ac, an = p.get("ac"), p.get("an")
        if ac is not None and an:
            afs[p.get("id", "")] = ac / an
    return afs, (max(afs.values()) if afs else None)


def query_gnomad(variant_ids, cache_dir, batch=BATCH, pause=1.0, retries=3):
    """Query gnomAD for a list of variant IDs; return {variant_id: record}. Cached + resumable."""
    os.makedirs(cache_dir, exist_ok=True)
    out = {}
    for i in range(0, len(variant_ids), batch):
        chunk = variant_ids[i:i + batch]
        cache_file = os.path.join(cache_dir, f"gnomad_{i}.json")
        if os.path.exists(cache_file):
            results = json.load(open(cache_file))
        else:
            aliases = "\n".join(
                f'  v{j}: variant(variantId: "{vid}", dataset: {DATASET}) {{ {_VARIANT_FIELDS} }}'
                for j, vid in enumerate(chunk)
            )
            query = "query {" + aliases + "}"
            results = None
            for attempt in range(retries):
                try:
                    r = requests.post(GNOMAD_URL, json={"query": query},
                                      headers={"Content-Type": "application/json"}, timeout=120)
                    r.raise_for_status()
                    data = r.json()
                    results = {vid: (data.get("data", {}).get(f"v{j}") or {})
                               for j, vid in enumerate(chunk)}
                    break
                except Exception:
                    if attempt == retries - 1:
                        results = {vid: {"_error": "request_failed"} for vid in chunk}
                    time.sleep(2 ** attempt)
            with open(cache_file, "w") as fh:
                json.dump(results, fh)
            time.sleep(pause)
        out.update(results)
    return out


def parse_gnomad_record(rec):
    """Flatten a gnomAD GraphQL variant record into a dict of annotation fields.

    Missing values are represented as None (caller maps to 'NA').
    """
    if not rec or rec.get("_error"):
        return {"gnomad_found": "error" if rec.get("_error") else "absent"}
    out = {"gnomad_found": "present"}
    for scope in ("genome", "exome"):
        d = rec.get(scope) or {}
        out[f"gnomad_{scope}_af"] = d.get("af")
        out[f"gnomad_{scope}_ac"] = d.get("ac")
        out[f"gnomad_{scope}_an"] = d.get("an")
        out[f"gnomad_{scope}_hom"] = d.get("homozygote_count")
        faf = d.get("faf95") or {}
        out[f"gnomad_{scope}_faf95_popmax"] = faf.get("popmax")
        out[f"gnomad_{scope}_faf95_pop"] = faf.get("popmax_population")
        afs, popmax = _pop_af(d.get("populations"))
        out[f"gnomad_{scope}_af_popmax"] = popmax
        out[f"gnomad_{scope}_populations"] = json.dumps(afs) if afs else None
    return out
