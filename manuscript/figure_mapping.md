# Manuscript Figure Mapping

| Manuscript figure | Source figure | Description | Data source | Generating script |
|---|---|---|---|---|
| Figure 1 | fig1_workflow | Study workflow | — | `src/figures.py` |
| Figure 2 | fig2_revel_distribution | REVEL score distribution | annotated dataset | `src/figures.py` |
| Figure 3 | fig3_sift_distribution | SIFT score distribution | annotated dataset | `src/figures.py` |
| Figure 4 | fig4_polyphen_distribution | PolyPhen-2 distribution | annotated dataset | `src/figures.py` |
| Figure 5 | fig5_agreement_matrix | Predictor agreement | annotated dataset | `src/figures.py` |
| Figure 6 | fig6_revel_vs_polyphen | REVEL vs PolyPhen-2 | annotated dataset | `src/figures.py` |
| Figure 7 | fig7_revel_vs_sift | REVEL vs SIFT | annotated dataset | `src/figures.py` |
| Figure 8 | fig8_revel_by_presence | Predictor distributions by gnomAD presence | annotated dataset | `src/figures.py` |
| Figure 9 | fig9_log10_af_vs_revel | log10(AF) vs REVEL | annotated dataset | `src/figures.py` |
| Figure 10 | fig10_population_frequency | gnomAD population frequency | annotated dataset | `src/figures.py` |
| Figure 11 | fig11_cohort_composition | Final cohort composition | cohort | `src/figures.py` |
| Figure 12 | fig12_predictor_vs_functional | REVEL vs Findlay score | with-functional dataset | `src/figures.py` |
| Figure 13 | fig13_concordance | Computational-functional concordance | cohort | `src/figures.py` |
| Figure 14 | fig14_functional_by_category | Functional score by REVEL category | with-functional dataset | `src/figures.py` |
| Figure 15 | fig15_evidence_availability | Evidence availability | cohort | `src/figures.py` |
| Figure 16 | fig16_domain_comparison | Functional score by domain (RING vs BRCT) | with-functional dataset | `src/figures.py` |
| Figure 17 | fig17_dace_region_coverage | Dace 2025 HAP1 coverage of our VUS by region | Phase 13 derived | `scripts/phase13_dace_analysis.py` |
| Figure 18 | fig18_dace_predictor_scatter | Predictor vs Dace score, central exons | Phase 13 derived | `scripts/phase13_dace_analysis.py` |
| Figure 19 | fig19_domain_contrast | Correspondence contrast across domains | Phase 13 derived | `scripts/phase13_dace_analysis.py` |

All figures are generated as SVG + 300-dpi PNG. Scientific figures are not modified by hand.