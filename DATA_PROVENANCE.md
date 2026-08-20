# Data Provenance

## RAW PUBLIC DATA (not redistributed; fetched and checksum-verified in CI)

| Resource | Source | Version / release | Access |
|---|---|---|---|
| ClinVar | `ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz` | build 260818-0035.1 | 2026-08-19 (sha256 `d2a8c9c2…`) |
| gnomAD | `gnomad.broadinstitute.org` GraphQL API | v4 (`gnomad_r4`) | 2026-08-19 |
| Ensembl VEP | `rest.ensembl.org` | release 116 | 2026-08-19 |
| REVEL | `rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip` | v1.3 (GRCh38) | 2026-08-19 |
| Findlay SGE scores | MaveDB `urn:mavedb:00000097-0-2` | "BRCA1 SGE Normalized Scores" | 2026-08-19 |
| Dace 2025 SGE HAP1 | `github.com/phoebedace/BRCA1_SGE_HAP1_HMEC` `data_for_plots/HAP1_variants_unfiltered.csv` | preprint 2025-08-16 (medRxiv doi:10.1101/2025.08.11.25333423; CC-BY) | 2026-08-20 |
| PubMed | NCBI E-utilities (`esearch`/`esummary`/`efetch`) | — | 2026-08-19 |
| AlphaMissense (Phase 13B) | `storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz` | hg38 (Google DeepMind 2023; CC BY-NC-SA 4.0) | 2026-08-20 (sha256 `0516cfd7…`) |
| BayesDel noAF (Phase 13B) | Zenodo record 11256843, `DataS2_BayesDel_all_possible_variants.zip` (fenglab `BayesDel_170824_noAF`) | BayesDel v1 noAF (GRCh37) | 2026-08-20 (sha256 `1e47aff0…`) |
| REVEL (Phase 13B controls) | `rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip` | v1.3 (GRCh38) | 2026-08-20 (downloading; chr17 only required) |

## DERIVED DATA (committed)

- `data/processed/brca1_vus_missense_annotated.tsv` — 1,904 annotated VUS (FROZEN; sha256 `7afe54db…`).
- `data/processed/brca1_vus_missense_with_functional.tsv` — annotated + Findlay fields.
- `data/processed_eval/vus_with_dace_scores_352.tsv` — PHASE 13 (new; NOT frozen): 1,904-set-derived
  VUS + Dace 2025 HAP1 replicate scores + region, generated from the Dace preprint table.
- `data/intermediate/*` — normalized, biological-variant map, unique-variant annotations, summaries.

## ANALYSIS OUTPUTS (committed)

- `results/tables/` — 21 generated tables (dataset, predictor, agreement, correlation, cohort,
  evidence, literature, conflicts; + 3 Phase 13 domain/conflict tables).
- `results/figures/` — fig1–fig19 (SVG + 300-dpi PNG; fig17–19 from `scripts/phase13_dace_analysis.py`).
- `results/reports/` — 36 phase + audit reports (incl. `phase13_feasibility_dace2025.md`,
  `phase13_domain_contrast.md`).
- `manuscript/supplementary/S5_dace_352.tsv` — Phase 13 Dace table (352 rows).

## PHASE 13B (derived data, new; NOT frozen)

- `data/processed_eval/vus_all_scores_1904.tsv` — 1,904 VUS + AlphaMissense (`am_pathogenicity`)
  and BayesDel noAF (`bayesdel_noAF`) scores; `g37_allele_consistent` verified True (0/1904
  mismatched GRCh37 ref/alt alleles vs GRCh38 — BayesDel join is safe).
- `data/processed_eval/controls_p_lp_b_lb.tsv` — 347 gold-standard ClinVar BRCA1 missense
  controls (≥2-star review status, germline origin; 159 P/LP, 188 B/LB) with AM/BayesDel
  scores + AM class.
- `data/processed_eval/controls_plp_am_likely_benign.tsv` / `controls_blb_am_likely_pathogenic.tsv`
  — AM-class callout mismatch lists.
- `results/tables/phase13b_{calibration_auc,domain_correlations,domain_difference}.tsv`
- `results/reports/phase13b_calibration.md` — PHASE 13B report (new; NOT frozen).

## Distinction

- **Raw public data** are external and re-fetched (not redistributed).
- **Derived data** are our processed datasets (frozen, checksummed).
- **Analysis outputs** are regenerable tables/figures/reports produced by `scripts/`.

All access dates and checksums are taken from the project's existing `data/raw/*/metadata.json`
and report files; none are invented.
