# Data Provenance

## RAW PUBLIC DATA (not redistributed; fetched and checksum-verified in CI)

| Resource | Source | Version / release | Access |
|---|---|---|---|
| ClinVar | `ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz` | build 260818-0035.1 | 2026-08-19 (sha256 `d2a8c9c2…`) |
| gnomAD | `gnomad.broadinstitute.org` GraphQL API | v4 (`gnomad_r4`) | 2026-08-19 |
| Ensembl VEP | `rest.ensembl.org` | release 116 | 2026-08-19 |
| REVEL | `rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip` | v1.3 (GRCh38) | 2026-08-19 |
| Findlay SGE scores | MaveDB `urn:mavedb:00000097-0-2` | "BRCA1 SGE Normalized Scores" | 2026-08-19 |
| PubMed | NCBI E-utilities (`esearch`/`esummary`/`efetch`) | — | 2026-08-19 |

## DERIVED DATA (committed)

- `data/processed/brca1_vus_missense_annotated.tsv` — 1,904 annotated VUS (FROZEN; sha256 `7afe54db…`).
- `data/processed/brca1_vus_missense_with_functional.tsv` — annotated + Findlay fields.
- `data/intermediate/*` — normalized, biological-variant map, unique-variant annotations, summaries.

## ANALYSIS OUTPUTS (committed)

- `results/tables/` — 18 generated tables (dataset, predictor, agreement, correlation, cohort,
  evidence, literature, conflicts).
- `results/figures/` — fig1–fig16 (SVG + 300-dpi PNG).
- `results/reports/` — 34 phase + audit reports.

## Distinction

- **Raw public data** are external and re-fetched (not redistributed).
- **Derived data** are our processed datasets (frozen, checksummed).
- **Analysis outputs** are regenerable tables/figures/reports produced by `scripts/`.

All access dates and checksums are taken from the project's existing `data/raw/*/metadata.json`
and report files; none are invented.
