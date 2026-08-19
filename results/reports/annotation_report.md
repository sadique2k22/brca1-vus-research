# Annotation Report — Phase 4B (QC + resources)

## Quality control

- total_records: 1904
- unique_variants: 1904
- gnomad_present: 0
- gnomad_absent: 4
- missing_cadd: 1904
- missing_revel: 0
- missing_sift: 0
- missing_polyphen: 0

## Resource usage

- Runtime: 401.6s
- gnomad_cache storage: 0.09 MB
- vep_cache storage: 52.95 MB

## Annotation metadata (provenance)

| Field | Source | Version | Assembly | Date | Method |
|---|---|---|---|---|---|
| gnomad_* | gnomAD | v4 (`gnomad_r4`) | GRCh38 | 2026-08-19 | GraphQL API |
| vep_*, sift_*, polyphen_* | Ensembl VEP | release 116 | GRCh38 | 2026-08-19 | REST /vep/human/region |
| revel_score | REVEL | v1.3 | GRCh38 | 2021-05-03 | standalone file (chr17 region) |
| cadd_phred | CADD | v1.6 (v1.7 degraded) | GRCh38 | 2026-08-19 | web service (best-effort) |

Missing values are recorded as empty (NA); no field is estimated.
