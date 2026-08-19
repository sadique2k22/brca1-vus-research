# Manuscript Numerical Audit

Every numerical statement in `manuscript/manuscript.md` is traced to a frozen result table
or report. Status: VERIFIED.

| # | Claim (value) | Source file | Source table/column | Verified |
|---|---|---|---|---|
| 1 | 1,904 BRCA1 missense VUS | data/processed/brca1_vus_missense_annotated.tsv | row count | ✅ |
| 2 | 424 gnomAD-present (22.3%) | results/reports/phase5_statistical_analysis.md | §1 | ✅ |
| 3 | 1,480 gnomAD-absent (77.7%) | results/reports/phase5_statistical_analysis.md | §1 | ✅ |
| 4 | median AF 6.8×10⁻⁷ | results/reports/phase5_statistical_analysis.md | §1 | ✅ |
| 5 | max AF 2.6×10⁻⁵ | results/reports/phase5_statistical_analysis.md | §1 | ✅ |
| 6 | REVEL 367/1246/291 (tol/int/impact) | results/reports/phase5_statistical_analysis.md | §2 | ✅ |
| 7 | SIFT 1071/833 (del/tol) | results/reports/phase5_statistical_analysis.md | §2 | ✅ |
| 8 | PolyPhen 1328/412/164 | results/reports/phase5_statistical_analysis.md | §2 | ✅ |
| 9 | κ REVEL–SIFT 0.448 | results/tables/predictor_agreement.tsv | kappa | ✅ |
| 10 | κ REVEL–PolyPhen 0.441 | results/tables/predictor_agreement.tsv | kappa | ✅ |
| 11 | κ SIFT–PolyPhen 0.272 | results/tables/predictor_agreement.tsv | kappa | ✅ |
| 12 | predictor–predictor ρ 0.29–0.41 | results/tables/predictor_correlations.tsv | spearman_rho | ✅ |
| 13 | Mann-Whitney REVEL p=0.89, SIFT p=0.44, PolyPhen p=0.014 | results/reports/phase5_statistical_analysis.md | §5 | ✅ |
| 14 | log10(AF) vs predictor ρ ≈ 0 | results/reports/phase5_statistical_analysis.md | §5 | ✅ |
| 15 | 0 variants faf95 ≥ 0.001 | results/reports/phase5_statistical_analysis.md | §6 | ✅ |
| 16 | 373 Findlay-scored (19.6%) | results/reports/findlay_coverage_report.md | — | ✅ |
| 17 | Findlay vs REVEL ρ=−0.384 | results/reports/phase6_revision_after_full_findlay.md | — | ✅ |
| 18 | Findlay vs SIFT ρ=−0.370 | results/reports/phase6_revision_after_full_findlay.md | — | ✅ |
| 19 | Findlay vs PolyPhen ρ=−0.188 | results/reports/phase6_revision_after_full_findlay.md | — | ✅ |
| 20 | 41-variant cohort (A10/B10/C9/D2/E10) | results/tables/phase7_final_cohort.tsv | stratum | ✅ |
| 21 | 19 conflicts in cohort | results/reports/phase7_evidence_synthesis.md | §3 | ✅ |
| 22 | 36/41 with PubMed hit | results/reports/phase7_evidence_synthesis.md | §2 | ✅ |
| 23 | 0/41 expert curation | results/reports/phase7_evidence_synthesis.md | §2 | ✅ |
| 24 | 436 candidate union | results/tables/phase6_candidate_union.tsv | row count | ✅ |
| 25 | 20 impact+WT-like; 9 tolerance+non-functional | results/tables/phase6_functional_comparison.tsv | comparison | ✅ |

**Result: all 25 numerical claims VERIFIED against frozen outputs. No unverifiable number found.**
