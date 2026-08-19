# Annotation Report — Phase 4B (QC + resources)

## Quality control

- total_records: 1904
- unique_variants: 1904
- gnomad_present: 424
- gnomad_absent: 1480
- missing_cadd: 1904
- missing_revel: 0
- missing_sift: 0
- missing_polyphen: 0

## Resource usage

- Runtime: 161.1s
- gnomad_cache storage: 0.16 MB
- vep_cache storage: 52.95 MB

## Annotation metadata (provenance)

| Field | Source | Version | Assembly | Date | Method |
|---|---|---|---|---|---|
| gnomad_* | gnomAD | v4 (`gnomad_r4`) | GRCh38 | 2026-08-19 | GraphQL API |
| vep_*, sift_*, polyphen_* | Ensembl VEP | release 116 | GRCh38 | 2026-08-19 | REST /vep/human/region |
| revel_score | REVEL | v1.3 | GRCh38 | 2021-05-03 | standalone file (chr17 region) |
| cadd_phred | CADD | v1.6 (v1.7 degraded) | GRCh38 | 2026-08-19 | web service (best-effort) |

**CADD exclusion:** CADD was not included in the primary analysis because the available v1.7 service did not provide a reliable bulk annotation pathway compatible with the project's resource constraints. No older CADD version was substituted. Missing CADD is not treated as a negative result and CADD is excluded from subsequent statistical analyses unless later added through a validated source.

Missing values are recorded as empty (NA); no field is estimated.
