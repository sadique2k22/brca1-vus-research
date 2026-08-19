# Phase 6 Revision — After Full Findlay Integration

## What was wrong / incomplete in Phase 6

- Phase 6 mapped Findlay scores onto only the **436 candidates** (reporting 43 with scores), not the full 1,904 VUS.
- The 'partial (RING + BRCT)' description was correct, but the coverage number (43) understated the available functional evidence.

## What the full integration adds

- Full-VUS Findlay coverage: **373/1904** variants gained a functional score.
- Findlay score distribution (n=373): median=-0.330071693871008, min=-4.591698941323781, max=0.844991185551513, q1=-0.9977668305910576, q3=-0.026930242173993858.

## Computational vs functional comparison

- intermediate: 201
- computational_impact + functional_LOF (agreement): 141
- computational_impact + functional_normal: 20
- computational_tolerance + functional_LOF: 9
- computational_tolerance + functional_normal (agreement): 2

## Findlay score vs predictors (Spearman)

| Comparison | rho | p | n |
|---|---|---|---|
| findlay_vs_revel | -0.384 | 1.4e-14 | 373 |
| findlay_vs_sift | -0.370 | 1.63e-13 | 373 |
| findlay_vs_polyphen | -0.188 | 0.000268 | 373 |

## Remaining limitations

- Findlay 2018 covers only RING + BRCT (13 exons); the DNA-binding domain (exons 6–14) has no functional score.
- Single-assay readout (HAP1 viability); does not capture all mechanisms.
- Computational-vs-functional comparison is descriptive; no ACMG PS3/BS3 applied.
