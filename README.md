# BRCA1 VUS Evidence Integration

A reproducible, provenance-preserving framework for integrating population, computational, and
experimental functional evidence around *BRCA1* missense variants of uncertain significance (VUS).

## Overview

We analyzed **1,904 *BRCA1* missense VUS** (ClinVar, GRCh38, transcript NM_007294.4) by annotating
each variant with gnomAD v4 population frequency, three computational predictors (REVEL, SIFT,
PolyPhen-2), and — where available — Findlay *et al.* (2018) saturation genome-editing
**HAP1 cellular-fitness** scores (**373 variants**, RING and BRCT domains only). A **41-variant
descriptive evidence cohort** was selected deterministically from computational/functional strata
and assessed against expert curation and abstract-level literature verification. Computational
predictors showed only moderate correspondence with the functional measurement
(Spearman ρ = −0.38 to −0.19), with substantial variant-level discordance (13/41 conflicts).
**No variant was clinically reclassified, and no pathogenic/benign label was assigned.**

## Research question

How consistently do population-frequency evidence and computational missense-variant predictors
correspond with available experimental functional evidence for *BRCA1* missense VUS?

## Important scope

- This is **not a clinical diagnostic tool**; no variant is classified pathogenic or benign.
- Findlay functional measurements cover the **RING and BRCT domains only** (not the DNA-binding domain).
- Literature verification was **abstract/metadata-level**; exact-variant evidence in full text was out of scope.
- **AlphaMissense was not included** (documented as future work; see Limitations).

## Repository structure

```
config/ , configs/     pipeline + analysis + cohort configuration (frozen)
data/                 raw (provenance) / intermediate / processed (frozen datasets)
src/                  pipeline modules
scripts/              phase orchestrators + validators
tests/                62 unit tests
results/tables/       18 generated tables
results/reports/      34 phase + audit reports
results/figures/      fig1–fig16 (SVG + 300-dpi PNG)
manuscript/           manuscript.md (+ historical versions) + supplementary/
```

## Reproducibility

- Python 3.12 (CI); dependencies in `requirements.txt`.
- All compute runs on GitHub Actions (`.github/workflows/pipeline.yml`); results are auto-committed.
- Frozen annotated dataset SHA-256: `7afe54db14718bcc…`.
- `python -m unittest discover -s tests` → 62/62.
- Full instructions: `REPRODUCIBILITY.md`.

## Results (validated highlights)

- 1,904 VUS → 424 gnomAD-present (22.3%), 1,480 absent (77.7%); present variants ultra-rare (median AF 6.8×10⁻⁷).
- Predictor agreement moderate (κ = 0.27–0.45); predictor–functional correspondence moderate for
  REVEL (ρ = −0.384) / SIFT (ρ = −0.370) and weak for PolyPhen-2 (ρ = −0.188), within the 373 RING/BRCT subset.
- 13/41 descriptive-cohort conflicts between computational and functional evidence.

## Manuscript

Canonical manuscript: `manuscript/manuscript.md` (Phase 11 revision; historical versions preserved).

## Supplementary material

- **S1** — all 1,904 variants (`manuscript/supplementary/S1_all_1904_variants.tsv`/`.xlsx`)
- **S2** — 373 variants with Findlay measurements (`S2_findlay_373.tsv`)
- **S3** — 41-variant evidence matrix (`S3_evidence_matrix_41.tsv`)
- **S4** — literature verification log (`S4_literature_verification.tsv`)
- Supplementary Methods: `manuscript/supplementary/supplementary_methods.md`

## Limitations

ClinVar ascertainment; VUS instability; gnomAD ancestry imbalance; predictor dependence/circularity;
exclusion of CADD and AlphaMissense; Findlay HAP1 single-assay and RING/BRCT-only scope; abstract-level
literature review; exploratory statistics (no multiple-testing correction); descriptive (non-powered)
41-variant cohort; no clinical validation.

## Citation

See `CITATION.cff` (author list pending) and `DATA_PROVENANCE.md`.
