# Annotation Source Audit — Phase 4B (STEP 0)

Date: 2026-08-19. Verified live against each source; versions/behaviour recorded below.

## Summary table

| Source | Official | Version (verified) | GRCh38 | API | Batch | Auth | Response size |
|---|---|---|---|---|---|---|---|
| gnomAD | Broad Institute | v4 (dataset `gnomad_r4`) | ✅ | GraphQL POST | aliases (multi-variant) | none (rate-limited) | ~KB/variant |
| Ensembl VEP | ensembl.org | release **116** | ✅ | REST POST `/vep/human/region` | ≤200 variants/request | none (rate-limited) | ~10–50 KB/variant |
| SIFT | via VEP (precomputed) | VEP 116 | ✅ | VEP REST | ✅ | — | included in VEP |
| PolyPhen-2 | via VEP (HumVar) | VEP 116 | ✅ | VEP REST | ✅ | — | included in VEP |
| CADD | cadd.gs.washington.edu | v1.7 | ✅ | web form (`/score`) | small chunks only | none | — |
| REVEL | sites.google.com/site/revelgenomics | standalone GRCh38 TSV | ✅ | none (file download) | — | none | ~GB file |

## Per-source detail

### 1. gnomAD (population frequency)
- **Official:** https://gnomad.broadinstitute.org — Broad Institute.
- **Version/assembly:** gnomAD **v4** (`gnomad_r4` dataset), **GRCh38**, exomes + genomes.
- **API:** GraphQL POST at `https://gnomad.broadinstitute.org/api`. Variant query
  `variant(variantId: "17-43092412-C-T", dataset: gnomad_r4)`.
- **Verified live:** positive control returned AF/AC/AN for common variants
  (e.g. rs4986852 AF≈0.0143 genome / 0.0198 exome); rare VUS correctly return
  "Variant not found" (i.e. genuinely absent). Dataset enum confirmed by introspection.
- **Fields:** `genome{ac,an,af}`, `exome{ac,an,af}`, plus `populations` (AF per ancestry),
  `homozygote_count`, `faf95`/`faf99` (filtering AF), `joint`.
- **Batch:** GraphQL `aliases` allow many variants per request (size-limited; cap ~100/request).
- **Rate limit:** not officially published → batching + disk cache + exponential backoff.
- **Auth:** none.
- **Cache:** required (per-variant JSON, keyed by variant key).

### 2. Ensembl VEP (consequence + SIFT + PolyPhen)
- **Official:** https://rest.ensembl.org — Ensembl **release 116**, GRCh38.
- **API:** POST `/vep/human/region` (`hgvs=1`), up to **200 variants/request**.
- **Verified live:** returns `consequence_terms`, `hgvsc`, `hgvsp`, `amino_acids`,
  `sift_score`, `sift_prediction`, `polyphen_score`, `polyphen_prediction` (HumVar),
  `gene_symbol`, `hgnc_id`, `codons`, `impact`.
- **Rate limit:** ~15 req/s anonymous (documented).
- **Auth:** none.
- **Cache:** already cached in Phase 4A (`data/intermediate/vep_cache/`); reuse + extract
  SIFT/PolyPhen fields (no re-query needed where cache is present).

### 3. SIFT
- Source: VEP (precomputed SIFT for GRCh38); original ref Ng & Henikoff 2003, PMID 12824425.
- No separate API needed — read `sift_score`/`sift_prediction` from VEP.

### 4. PolyPhen-2
- Source: VEP (HumVar model); original ref Adzhubei 2010, PMID 20354512.
- No separate API needed — read `polyphen_score`/`polyphen_prediction` from VEP.

### 5. CADD
- **Official:** https://cadd.gs.washington.edu — version **v1.7**, GRCh38.
- **API:** web form at `/score` (paste/upload VCF); **no clean bulk REST API**. POST to
  `/score` returned `405`. Precomputed whole-genome scores are ~300 GB (prohibited).
- **For ~1,904 variants:** only the web form in small chunks (slow/fragile), or skip.
  **DECISION POINT.**

### 6. REVEL
- **Official:** https://sites.google.com/site/revelgenomics (reachable); original ref
  Ioannidis 2016, PMID 27666373.
- **Distribution:** standalone GRCh38 TSV via Google Drive (~GB); also a column in
  dbNSFP (jpopgen site reachable; full dbNSFP ≈ 15–30 GB → discouraged under the 50 GB cap).
- **No API.** **For ~1,904 variants:** download the standalone GRCh38 REVEL file (~GB,
  one-time), extract the chr17 BRCA1 region (~43.0–43.2 Mb). **DECISION POINT (download
  approval + size verification).**

## Proposed strategy (mobile-safe, restartable)

1. **VEP REST** (batched 200, cached) → consequence + HGVS + **SIFT + PolyPhen-2**.
2. **gnomAD GraphQL** (`gnomad_r4`) → AF/AC/AN/AF_popmax/hom/faf, batched via aliases, cached.
3. **CADD** — *pending decision* (web form chunks vs skip).
4. **REVEL** — *pending decision* (standalone GRCh38 file vs skip).

All requests: disk-cached keyed by `variant_key`; retries with exponential backoff;
checkpointing so an interrupted run resumes without re-fetching.

## Decision points (need approval before STEP 1)
- **D1 — REVEL:** approve downloading the standalone REVEL GRCh38 file (~GB, verify size
  first)? Or leave REVEL missing?
- **D2 — CADD:** attempt the CADD web form in small chunks (slow, may rate-limit), or skip
  CADD and rely on SIFT/PolyPhen/REVEL?
