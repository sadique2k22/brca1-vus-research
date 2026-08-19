# Phase 11 Manuscript Numerical Audit

Every numerical statement in `manuscript/manuscript_phase11.md` is traced to a frozen or
Phase 9-validated output. Status: VERIFIED.

| # | Claim (value) | Source | Verified |
|---|---|---|---|
| 1 | 1,904 BRCA1 missense VUS | data/processed/brca1_vus_missense_annotated.tsv | ✅ |
| 2 | 424 present (22.3%); 1,480 absent (77.7%) | results/reports/phase5_statistical_analysis.md | ✅ |
| 3 | median AF 6.8×10⁻⁷; max 2.6×10⁻⁵ | results/reports/phase5_statistical_analysis.md | ✅ |
| 4 | REVEL 367/1246/291 | results/reports/phase5_statistical_analysis.md | ✅ |
| 5 | SIFT 1071/833 | results/reports/phase5_statistical_analysis.md | ✅ |
| 6 | PolyPhen 1328/412/164 | results/reports/phase5_statistical_analysis.md | ✅ |
| 7 | κ 0.448 / 0.441 / 0.272 | results/tables/predictor_agreement.tsv | ✅ |
| 8 | predictor–predictor ρ 0.29–0.41 | results/tables/predictor_correlations.tsv | ✅ |
| 9 | Mann-Whitney REVEL 0.89 / SIFT 0.44 / PolyPhen 0.014 | results/reports/phase5_statistical_analysis.md | ✅ |
| 10 | 373 scored (19.6%); 124 RING / 236 BRCT | results/reports/findlay_coverage_report.md | ✅ |
| 11 | REVEL ρ=−0.384 (CI −0.468,−0.294) | results/reports/phase9_final_report.md | ✅ |
| 12 | SIFT ρ=−0.370 (CI −0.454,−0.278) | results/reports/phase9_final_report.md | ✅ |
| 13 | PolyPhen ρ=−0.188 (CI −0.284,−0.088) | results/reports/phase9_final_report.md | ✅ |
| 14 | 13 conflicts (32%) | results/reports/phase9_final_report.md | ✅ |
| 15 | 36/41 with PubMed record; 295 records; 0 exact-variant | results/reports/phase9_literature_verification.md | ✅ |
| 16 | RING −0.52 / BRCT −0.30; p=0.059; r=0.12 | results/reports/phase9_final_report.md | ✅ |

**Result: all 16 numerical claims VERIFIED.**
