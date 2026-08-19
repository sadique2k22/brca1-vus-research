# Manuscript Figure Mapping

| Manuscript figure | Source figure | Description | Data source | Generating script |
|---|---|---|---|---|
| Figure 1 | fig1_workflow | Study workflow | — | `src/figures.py` |
| Figure 2 | fig2_revel_distribution | REVEL score distribution | annotated dataset | `src/figures.py` |
| Figure 3 | fig3_sift_distribution | SIFT score distribution | annotated dataset | `src/figures.py` |
| Figure 4 | fig4_polyphen_distribution | PolyPhen-2 distribution | annotated dataset | `src/figures.py` |
| Figure 5 | fig5_agreement_matrix | Predictor agreement | annotated dataset | `src/figures.py` |
| Figure 6 | fig6_revel_vs_polyphen | REVEL vs PolyPhen-2 | annotated dataset | `src/figures.py` |
| Figure 7 | fig8_*_by_presence | Predictor distributions by gnomAD presence | annotated dataset | `src/figures.py` |
| Figure 8 | fig9_log10_af_vs_revel | log10(AF) vs REVEL | annotated dataset | `src/figures.py` |
| Figure 9 | fig11_cohort_composition | Final cohort composition | cohort | `src/figures.py` |
| Figure 10 | fig12_predictor_vs_functional | REVEL vs Findlay score | with-functional dataset | `src/figures.py` |
| Figure 11 | fig13_concordance | Computational-functional concordance | cohort | `src/figures.py` |
| Figure 12 | fig14_functional_by_category | Functional score by REVEL category | with-functional dataset | `src/figures.py` |
| Figure 13 | fig15_evidence_availability | Evidence availability | cohort | `src/figures.py` |
| Figure 14 | fig16_domain_comparison | Functional score by domain (RING vs BRCT) | with-functional dataset | `src/figures.py` |

All figures are generated as SVG + 300-dpi PNG. Scientific figures are not modified by hand.
