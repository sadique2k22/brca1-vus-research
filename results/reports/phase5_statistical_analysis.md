# Phase 5 — Statistical Analysis Report

**Exploratory / descriptive only.** No pathogenic/benign classification, no ACMG codes, no clinical interpretation.

Frozen dataset checksum verified: `7afe54db14718bcc…`. Runtime 3.3s.

## 1. Dataset summary

- Total missense VUS: **1904**
- gnomAD present: 424 (22.3%)
- gnomAD absent: 1480 (77.7%)
- gnomAD-present global AF: median=6.84e-07, mean=1.69e-06, min=0, max=2.63e-05, q1=6.84e-07, q3=1.37e-06

## 2. Predictor distributions

| Predictor | n | min | max | mean | median | q1 | q3 |
|---|---|---|---|---|---|---|---|
| revel | 1904 | 0.001 | 0.956 | 0.475 | 0.522 | 0.371 | 0.600 |
| sift | 1904 | 0.000 | 1.000 | 0.108 | 0.040 | 0.000 | 0.140 |
| polyphen | 1904 | 0.000 | 1.000 | 0.303 | 0.129 | 0.019 | 0.568 |

Category counts (neutral terminology):

- revel: {'tolerance': 367, 'intermediate': 1246, 'impact': 291}
- sift: {'deleterious': 1071, 'tolerated': 833}
- polyphen: {'benign': 1328, 'possibly_damaging': 412, 'probably_damaging': 164}

## 3. Predictor agreement

Pairwise categorical agreement (2-class impact/tolerance; REVEL-intermediate excluded):

| Comparison | % agreement | Cohen's kappa | p | n |
|---|---|---|---|---|
| revel_vs_sift | 71.1% | 0.448 | 0 | 658 |
| revel_vs_polyphen | 72.9% | 0.441 | 0 | 658 |
| sift_vs_polyphen | 61.8% | 0.272 | 0 | 1904 |

Full REVEL x SIFT x PolyPhen combination table (count):

| REVEL | SIFT | PolyPhen | count |
|---|---|---|---|
| intermediate | tolerated | benign | 542 |
| intermediate | deleterious | benign | 363 |
| intermediate | deleterious | possibly_damaging | 175 |
| tolerance | tolerated | benign | 169 |
| tolerance | deleterious | benign | 137 |
| impact | deleterious | benign | 111 |
| impact | deleterious | possibly_damaging | 93 |
| intermediate | tolerated | possibly_damaging | 88 |
| impact | deleterious | probably_damaging | 77 |
| intermediate | deleterious | probably_damaging | 72 |
| tolerance | deleterious | possibly_damaging | 36 |
| tolerance | tolerated | possibly_damaging | 17 |
| tolerance | deleterious | probably_damaging | 7 |
| impact | tolerated | benign | 6 |
| intermediate | tolerated | probably_damaging | 6 |
| impact | tolerated | possibly_damaging | 3 |
| impact | tolerated | probably_damaging | 1 |
| tolerance | tolerated | probably_damaging | 1 |

## 4. Predictor correlations (Spearman)

(SIFT direction inverted to 'higher = more damaging' for comparability.)

| Comparison | rho | p | n |
|---|---|---|---|
| revel_vs_polyphen | 0.289 | 8.05e-38 | 1904 |
| revel_vs_sift | 0.329 | 3.55e-49 | 1904 |
| sift_vs_polyphen | 0.408 | 2.29e-77 | 1904 |

## 5. gnomAD present vs absent (predictor scores)

Mann-Whitney U (exploratory):

| Predictor | U | p | n_present | n_absent |
|---|---|---|---|---|
| revel | 315102 | 0.893 | 424 | 1480 |
| sift | 306074 | 0.436 | 424 | 1480 |
| polyphen | 289283 | 0.0142 | 424 | 1480 |

log10(global AF) vs predictor (Spearman, gnomAD-present only):

| Predictor | rho | p | n |
|---|---|---|---|
| revel | 0.011 | 0.832 | 406 |
| sift | 0.027 | 0.587 | 406 |
| polyphen | -0.040 | 0.422 | 406 |

## 6. Population-frequency findings

Variants with elevated population filtering AF (faf95 popmax >= 0.001): **0**


(Full list in `results/tables/population_frequency_outliers.tsv`.)

## 7. Pattern classes (Phase 6 literature candidates)

| Class | Description | Count |
|---|---|---|
| A | gnomAD-present + high impact | 61 |
| B | gnomAD-present + tolerance | 84 |
| C | gnomAD-absent + strong impact | 4 |
| D | gnomAD-absent + tolerance | 283 |
| E | strong predictor disagreement | 49 |
| F | elevated population frequency | 0 |
| G | extreme REVEL score | 7 |

`results/tables/pattern_candidates.tsv` — NOT a pathogenicity list.

## 8. Statistical limitations

- Exploratory; multiple comparisons were not corrected (Mann-Whitney x3, correlations x6). p-values are descriptive, not confirmatory.
- REVEL/SIFT/PolyPhen are correlated (REVEL is an ensemble) — treated as correlated predictors, not independent evidence lines.
- gnomAD-absent is not AF=0; absence is analyzed as a category, not a value.
- AF↔predictor correlation is partly circular (predictors are conservation-trained).

## 9. Data limitations

- Single gene (BRCA1), single transcript (NM_007294.4), missense VUS only.
- CADD excluded (v1.7 bulk annotation unavailable).
- Full per-population AF not retrieved (gnomAD GraphQL cost limit); faf95 popmax used.

## 10. Phase 6 recommendation

Prioritize literature review of classes C (absent + strong impact) and E (disagreement), with F (elevated population frequency) as a cross-check group. No variant is claimed pathogenic or benign.
