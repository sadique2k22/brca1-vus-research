# Phase 9 Manuscript Numerical Audit

Every numerical statement in the revised `manuscript/manuscript.md` is traced to a frozen
or Phase 9-validated output. Status: VERIFIED.

| # | Claim (value) | Source file | Verified |
|---|---|---|---|
| 1 | 1,904 BRCA1 missense VUS | data/processed/brca1_vus_missense_annotated.tsv | ✅ |
| 2 | 424 gnomAD-present (22.3%) | results/reports/phase5_statistical_analysis.md | ✅ |
| 3 | 1,480 gnomAD-absent (77.7%) | results/reports/phase5_statistical_analysis.md | ✅ |
| 4 | median AF 6.8×10⁻⁷; max 2.6×10⁻⁵ | results/reports/phase5_statistical_analysis.md | ✅ |
| 5 | REVEL 367/1246/291 | results/reports/phase5_statistical_analysis.md | ✅ |
| 6 | SIFT 1071/833 | results/reports/phase5_statistical_analysis.md | ✅ |
| 7 | PolyPhen 1328/412/164 | results/reports/phase5_statistical_analysis.md | ✅ |
| 8 | κ 0.448 / 0.441 / 0.272 | results/tables/predictor_agreement.tsv | ✅ |
| 9 | 373 Findlay-scored (19.6%) | results/reports/findlay_coverage_report.md | ✅ |
| 10 | REVEL vs Findlay ρ=−0.384 (95% CI −0.468,−0.294) | results/reports/phase9_final_report.md | ✅ |
| 11 | SIFT vs Findlay ρ=−0.370 (95% CI −0.454,−0.278) | results/reports/phase9_final_report.md | ✅ |
| 12 | PolyPhen vs Findlay ρ=−0.188 (95% CI −0.284,−0.088) | results/reports/phase9_final_report.md | ✅ |
| 13 | 13 conflicts in cohort (32%) | results/reports/phase9_final_report.md | ✅ |
| 14 | 36/41 with PubMed record; 0 exact-variant in abstract | results/reports/phase9_literature_verification.md | ✅ |
| 15 | RING median −0.52 (n=124); BRCT median −0.30 (n=236) | results/reports/phase9_final_report.md | ✅ |
| 16 | Mann-Whitney p=0.059; rank-biserial r=0.12 | results/reports/phase9_final_report.md | ✅ |
| 17 | 0/41 expert curation | results/reports/phase7_evidence_synthesis.md | ✅ |

**Result: all revised numerical claims VERIFIED.**
