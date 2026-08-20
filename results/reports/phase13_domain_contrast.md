# Phase 13 — Domain-contrast analysis (Dace 2025 integration)

Date: 2026-08-20. Frozen inputs read-only; all outputs are new Phase 13 artifacts.

## Orientation convention (used everywhere in Phase 13)

- Impact orientation: REVEL raw, SIFT as (1 - SIFT), PolyPhen-2 raw (higher = predicted impact). Identical to `phase9.py`.

- Raw SIFT (higher = tolerated) is reported alongside for direct comparability with Findlay et al. 2018 (their reported SIFT rho = +0.363).

- Functional scores: higher = WT-like/fitter in BOTH assays.

## Per-domain Spearman correlations (impact orientation; 95% CI Fisher z)

| Domain | Predictor | n | rho | 95% CI |
|---|---|---|---|---|
| RING/BRCT (Findlay 2018) | REVEL | 373 | -0.384 | -0.468 to -0.294 |
| RING/BRCT (Findlay 2018) | SIFT (1−SIFT) | 373 | -0.369 | -0.454 to -0.278 |
| RING/BRCT (Findlay 2018) | PolyPhen-2 | 373 | -0.188 | -0.284 to -0.088 |
| RING/BRCT (Findlay 2018) | SIFT raw (higher = tolerated) | 373 | +0.369 | 0.278 to 0.454 |
| central exons (Dace 2025) | REVEL | 352 | -0.049 | -0.153 to 0.055 |
| central exons (Dace 2025) | SIFT (1−SIFT) | 352 | -0.260 | -0.355 to -0.16 |
| central exons (Dace 2025) | PolyPhen-2 | 352 | +0.047 | -0.057 to 0.151 |
| central exons (Dace 2025) | SIFT raw (higher = tolerated) | 352 | +0.260 | 0.16 to 0.355 |

## Domain-difference tests (RING/BRCT vs central exons)

| Predictor | Fisher z | Fisher p | Permutation p (10,000, seed 20260820) |
|---|---|---|---|
| REVEL | -4.767 | 1.87e-06 | 1.00e-04 |
| SIFT (1−SIFT) | -1.631 | 1.03e-01 | 7.29e-02 |
| PolyPhen-2 | -3.181 | 1.47e-03 | 1.10e-03 |

## Conflict analysis on the 352 newly covered VUS (phase9-identical REVEL rules)

- REVEL impact (>=0.644) with WT-like Dace score (>0): 14 / 352
- REVEL tolerance (<=0.290) with Dace score < -1.0: 14 / 352
- Sensitivity (Dace-paper LoF threshold -0.799): 30 conflicts total
- For reference, RING/BRCT 373 (same rules): 23 conflicts (6.2%)
- Central 352 (same rules, -1.0): 28 conflicts (8.0%)
- Per-predictor exploratory (impact call with WT-like Dace score): REVEL 14, SIFT 71, PolyPhen-2 11 (no 41-cohort comparator; descriptive only)

## Orientation audit (fixes applied in Phase 13)

- `phase9.py` computes SIFT as (1 - SIFT) for correlations; the SIFT negative sign in the manuscript therefore arises from SIFT's own inversion, NOT from the functional-score orientation. The previous sign-rationale sentence applied the REVEL-style rationale to SIFT, which is incorrect with raw SIFT (rho = +0.370).
- Raw-SIFT rho +0.370 (RING/BRCT, n = 373) is essentially identical to the +0.363 reported by Findlay et al. 2018 — concordant, not merely in magnitude.
- The Phase 13 feasibility report's earlier “SIFT is oppositely signed” statement mixed raw-SIFT (central) with impact-oriented SIFT (RING/BRCT); that was an orientation artifact. With matched orientation, SIFT correspondence attenuates (0.370 → 0.260 impact orientation) but does not change sign.

## Interpretation

- REVEL and PolyPhen-2 correspondence essentially disappears in the central exons (REVEL -0.384 → -0.049; PolyPhen-2 -0.188 → 0.047); SIFT attenuates but remains directionally consistent (-0.370 → -0.260, impact orientation). Domain-difference tests: see above.
- Caveats: Dace et al. is a preprint; unfiltered continuous table used; the 41-variant cohort and its 13 conflicts are untouched (all RING/BRCT or unassayed); the conflict rates above are descriptive, on different variant sets.

## Artifacts

- results/tables/phase13_domain_correlations.tsv (this table)
- results/tables/phase13_domain_difference.tsv (Fisher z + permutation)
- results/tables/phase13_dace_conflicts.tsv (conflict rows, central + RING/BRCT)
- results/figures/fig17_dace_region_coverage.{svg,png}
- results/figures/fig18_dace_predictor_scatter.{svg,png}
- results/figures/fig19_domain_contrast.{svg,png}
