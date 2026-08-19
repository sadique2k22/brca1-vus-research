"""gnomAD population-frequency retrieval via the GraphQL API (cached, batched, resumable).

Uses the public gnomAD GraphQL endpoint (dataset `gnomad_r4`, v4, GRCh38). The API enforces
a GraphQL *cost limit* (max 25), so we query a compact field set in batches of 25.

Fields retrieved: global AC/AN/AF, homozygote count, and faf95 popmax (filtering AF 95% CI,
population max — used as the population-aware maximum frequency). The full per-population
breakdown is NOT retrieved (would exceed the API cost limit for bulk queries); documented
in the annotation metadata.
"""
import json
import os
import time

import requests

GNOMAD_URL = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r4"
BATCH = 25

_FIELDS = """
  variantId
  genome { ac an af homozygote_count faf95 { popmax popmax_population } }
  exome  { ac an af homozygote_count faf95 { popmax popmax_population } }
"""


def _run_query(aliases, cache_file, retries=3):
    if os.path.exists(cache_file):
        return json.load(open(cache_file))
    query = "query {" + aliases + "}"
    reason = "unknown"
    for attempt in range(retries):
        try:
            r = requests.post(GNOMAD_URL, json={"query": query},
                              headers={"Content-Type": "application/json"}, timeout=120)
            if r.status_code == 429:
                reason = "rate_limited"
                time.sleep(int(r.headers.get("Retry-After", 5)))
                continue
            r.raise_for_status()
            data = r.json()
            msgs = [e.get("message", "") for e in (data.get("errors") or [])]
            if any("too expensive" in m or "cost" in m.lower() for m in msgs):
                reason = "cost_limit"
                time.sleep(5 * (2 ** attempt))
                continue
            with open(cache_file, "w") as fh:
                json.dump(data, fh)
            return data
        except requests.exceptions.Timeout:
            reason = "timeout"
        except requests.exceptions.HTTPError as e:
            reason = f"http_{e.response.status_code if e.response is not None else 'error'}"
        except requests.exceptions.ConnectionError:
            reason = "connection_error"
        except Exception:
            reason = "error"
        time.sleep(5 * (2 ** attempt))
    return {"_failed": True, "_reason": reason}


def query_gnomad(variant_ids, cache_dir, batch=BATCH, pause=2.0):
    """Query gnomAD for variant IDs; return {variant_id: record}. Cached + resumable.

    'Variant not found' (genuine absence) maps to {} -> 'absent'; only real failures
    (cost limit, network) map to {'_error': ...}.
    """
    os.makedirs(cache_dir, exist_ok=True)
    out = {}
    for i in range(0, len(variant_ids), batch):
        chunk = variant_ids[i:i + batch]
        cache_file = os.path.join(cache_dir, f"gnomad_{i}.json")
        aliases = "\n".join(
            f'  v{j}: variant(variantId: "{vid}", dataset: {DATASET}) {{ {_FIELDS} }}'
            for j, vid in enumerate(chunk)
        )
        data = _run_query(aliases, cache_file)
        if data.get("_failed"):
            for vid in chunk:
                out[vid] = {"_error": data.get("_reason", "request_failed")}
        else:
            d = data.get("data") or {}
            for j, vid in enumerate(chunk):
                out[vid] = d.get(f"v{j}") or {}
        time.sleep(pause)
    return out


def parse_gnomad_record(rec):
    """Flatten a gnomAD GraphQL variant record into annotation fields.

    Missing values are represented as None (caller maps to 'NA').
    """
    if not rec:
        return {"gnomad_found": "absent"}
    if rec.get("_error"):
        return {"gnomad_found": "error", "gnomad_error": rec.get("_error", "request_failed")}
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
    return out
