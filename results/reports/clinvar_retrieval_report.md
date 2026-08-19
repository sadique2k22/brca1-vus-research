# ClinVar Retrieval Report — Phase 2

Date: 2026-08-19. Protocol v1.0. Gene: BRCA1.

## Source & release

| Item | Value |
|---|---|
| Source URL | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz |
| File | `data/raw/clinvar/variant_summary.txt.gz` (immutable) |
| Release / build | ClinVar build **260818-0035.1** |
| Last modified (server) | 2026-08-17 18:51:44 GMT |
| Retrieval date (UTC) | 2026-08-19T02:59:49Z |
| File size | 442,283,371 bytes (~442 MB) |
| SHA-256 | `d2a8c9c2f038c70325a55a82801b059aa6e9325f958eca8a280625d19ddb27a8` |
| MD5 | `90c969bbf917a9eced1a8768d2cc33f1` (matches official ClinVar `.md5`) ✅ |
| Columns | 43 |

Metadata recorded in `data/raw/clinvar/metadata.json`.

## Extraction summary

| Metric | Count |
|---|---|
| Raw records (data rows, whole file) | **9,044,810** |
| Records for BRCA1 (rows, GeneSymbol match) | **31,542** |
| Unique variants / unique Variation IDs | **16,023** |
| Unique Allele IDs | 16,023 |

Rows per variant reflect ClinVar's per-assembly duplication:
GRCh37 = 15,818; GRCh38 = 15,529; Assembly "na" = 195.
All 15,529 GRCh38-mapped variants are on chromosome 17 ✅ (sanity check).

## Clinical significance distribution (unique Variation IDs)

| Significance | Count |
|---|---|
| Pathogenic | 4,133 |
| Conflicting classifications of pathogenicity | 2,928 |
| Likely benign | 2,902 |
| **Uncertain significance** | **2,583** |
| Benign | 742 |
| Likely pathogenic | 337 |
| Pathogenic/Likely pathogenic | 229 |
| Benign/Likely benign | 95 |
| not provided | 50 |
| no classification for the single variant | 3 |
| no classifications from unflagged records | 1 |
| VUS-high | 1 |
| (no germline classification, "-") | 2,019 |

## Review status distribution (unique Variation IDs)

| Review status | Count |
|---|---|
| criteria provided, single submitter | 4,808 |
| reviewed by expert panel | 3,408 |
| criteria provided, conflicting classifications | 2,926 |
| criteria provided, multiple submitters, no conflicts | 2,324 |
| no assertion criteria provided | 485 |
| no classification provided | 49 |
| no classification for the single variant | 3 |
| no classifications from unflagged records | 1 |
| ("-") | 2,019 |

## Completeness

| Metric | Count (of 16,023 unique VIDs) |
|---|---|
| Missing HGVS c. (no `:c.` in Name) | 516 |
| Missing protein consequence (no `p.` in Name) | 3,819 |

## Data-quality findings (discovered, not fabricated)

1. **Per-assembly row duplication.** Each variant has separate GRCh37 and GRCh38 rows
   (`Assembly` column). 289 variants have a GRCh37 row but no GRCh38 row; 195 rows have
   `Assembly = "na"`. Phase 3 normalization must select the GRCh38 row and deduplicate by
   VariationID.
2. **Deprecated allele columns.** `ReferenceAllele` / `AlternateAllele` are `na`; the actual
   alleles are in `ReferenceAlleleVCF` / `AlternateAlleleVCF`. (This resolves the Phase 0
   audit item C3 — the build/column question was verified empirically, not assumed.)
3. **Aggregate VUS = 2,583, not the Phase 0 API estimate of 7,930.** The Phase 0 `esearch`
   query matched *submission-level* significance, whereas `variant_summary.txt`
   `ClinicalSignificance` is the *aggregate* germline classification. The parsed value
   (2,583) is authoritative for this project (resolves protocol U3).
4. **2,019 variants have no germline classification (`-`)** and will fall outside the VUS
   set; they are retained in the intermediate file for transparency, not silently dropped.
5. **One `VUS-high` value** and a few `no classification…` values are unusual ClinVar
   categories; retained verbatim and reported.
6. **516 variants lack HGVS c. and 3,819 lack protein consequence** — expected for
   non-coding/synonymous/structural records; relevant to Phase 3 missense filtering, where
   the p. notation (or VEP) will be used to define "missense".

## Validation (automated)

`scripts/validate_clinvar.py` → **PASSED**:
raw file exists + size matches; sha256 matches; md5 matches metadata *and* the official
ClinVar `.md5`; a full re-parse reproduces the intermediate dataset byte-for-byte; the raw
file is confirmed unmodified.

## Artifacts produced

- `data/raw/clinvar/variant_summary.txt.gz` (+ `.md5`, `metadata.json`)
- `data/intermediate/clinvar_brca1_raw.tsv` (31,542 BRCA1 rows, verbatim 43 columns)
- `data/intermediate/clinvar_brca1_summary.json`
- `src/clinvar.py`, `scripts/validate_clinvar.py`

**Stopped here per instructions — no gnomAD annotation, no predictors, no filtering yet.**
