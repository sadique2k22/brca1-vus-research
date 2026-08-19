# Phase 7 — Final Evidence Synthesis

Evidence synthesis only. No pathogenic/benign classification, no ACMG codes.

## 1. Final cohort

- Cohort size: **41**
- Stratum counts: {'A': 10, 'B': 10, 'C': 9, 'D': 2, 'E': 10}

## 2. Evidence availability

- With >=1 PubMed hit: 36/41
- With expert-panel curation: 0/41
- With Findlay exact-variant score: 41/41 (by cohort definition)

## 3. Conflict analysis

- Evidence conflicts in cohort: **19**

## 4. Predictor vs functional (ordinal, Spearman)

(ROC-AUC not computed: a validated binary Findlay threshold was not retrieved; ordinal correlation is reported instead.)

| Predictor | rho | p | n |
|---|---|---|---|
| revel | -0.384 | 1.4e-14 | 373 |
| sift | -0.370 | 1.63e-13 | 373 |
| polyphen | -0.188 | 0.000268 | 373 |

(Full-scored n = 373.)

## 5. Domain analysis (Findlay score)

- RING: n=124, median=-0.52312978576273
- BRCT: n=236, median=-0.29921714020573553

## 6. Limitations

- Findlay covers RING + BRCT only; cohort is drawn from these regions.
- Clinical/segregation evidence is not programmatically extractable at metadata level.
- PubMed is title/metadata level.
