# Predictor Interpretation Thresholds — Phase 5

Frozen 2026-08-19. Locked BEFORE inspecting variant-level results. These thresholds are
used only for **neutral descriptive categorization** (impact-supporting / intermediate /
tolerance-supporting). They are NOT used to classify variants as pathogenic or benign.

## REVEL

| Threshold | Value | Interpretation |
|---|---|---|
| tolerance-supporting (BP4) | ≤ 0.290 | computational-tolerance-supporting |
| indeterminate | 0.290 < s < 0.644 | intermediate/indeterminate |
| impact-supporting (PP3) | ≥ 0.644 | computational-impact-supporting |
| impact-moderate (PP3_Moderate) | ≥ 0.773 | (reported for reference) |
| impact-strong (PP3_Strong) | ≥ 0.932 | strong impact signal |
| tolerance-strong (very-strong benign) | ≤ 0.003 | (reported for reference) |

**Source:** Pejaver V, Byrne AB, Feng BJ, et al. "Calibration of computational tools for
missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4
criteria." *Am J Hum Genet* 2022;109(12):2163-2177. PMID 36413997 (Table 2).
DOI 10.1016/j.ajhg.2022.10.013.

**BRCA1-specific note:** the ENIGMA/ClinGen BRCA1/BRCA2 VCEP (Parsons et al. 2024,
PMID 39142283) applies its *own* bioinformatic code using **BayesDel** (BRCA1 PP3 ≥ 0.28),
not REVEL. No BRCA1-specific REVEL threshold exists; therefore the generic ClinGen-calibrated
REVEL thresholds above are used for descriptive purposes only.

## SIFT

| Threshold | Value | Interpretation |
|---|---|---|
| deleterious | ≤ 0.05 | predicted deleterious (impact) |
| tolerated | > 0.05 | predicted tolerated |

**Source:** Ng PC, Henikoff S. "SIFT: predicting amino acid changes that affect protein
function." *Nucleic Acids Res* 2003;31(13):3812-3814. PMID 12824425 (developer-recommended
threshold).

## PolyPhen-2 (HumVar)

| Threshold | Value | Interpretation |
|---|---|---|
| benign | ≤ 0.446 | predicted benign (tolerance) |
| possibly damaging | 0.446 < s ≤ 0.908 | predicted possibly damaging |
| probably damaging | > 0.908 | predicted probably damaging (impact) |

**Source:** Adzhubei IA, et al. "A method and server for predicting damaging missense
mutations." *Nat Methods* 2010;7(4):248-249. PMID 20354512.

## gnomAD

- `absent` is **never** encoded as AF=0; it is a category, kept separate from numerical AF.
- Reference point for "elevated population-specific frequency" reporting: filtering AF
  (faf95 popmax) ≥ 0.001 (BRCA1 BS1 threshold; reported, not applied as a classification).
  Source: ENIGMA/ClinGen BRCA1/BRCA2 VCEP, PMID 39142283.

## Important caveats

- REVEL is an **ensemble** meta-predictor; SIFT/PolyPhen are constituents of several
  meta-predictors. The three are correlated and are **not** treated as independent
  evidence lines in this project.
- These categories describe *computational evidence patterns only*.
