# Phase 13B: AlphaMissense + BayesDel extension; ClinVar calibration controls

Data: ClinVar release 2026-08-19 (variant_summary 260818-0035.1); AlphaMissense hg38 (Google DeepMind, CC BY-NC-SA 4.0); BayesDel v1 noAF all-possible-variants (Zenodo record 11256843, DataS2; fenglab 170824 build).

## Calibration controls (gold-standard review status)

Controls: 347 BRCA1 missense variants from ClinVar with review status in ['criteria provided, multiple submitters, no conflicts', 'practice guideline', 'reviewed by expert panel'] and germline origin. Distribution: {'P/LP': 159, 'B/LB': 188}.

## AUC (P/LP vs B/LB) — higher score = more pathogenic

| predictor | AUC | 95% CI | n(P/LP) | n(B/LB) |
|---|---|---|---|---|
| AlphaMissense | 0.961 | 0.937–0.981 | 154 | 188 |
| BayesDel noAF | 0.938 | 0.911–0.962 | 159 | 188 |

## Domain-contrast (Spearman rho, predictor vs functional)

| predictor | set | n | rho | p |
|---|---|---|---|---|
| AlphaMissense | RING/BRCT (373, Findlay) | 373 | -0.4918 | 4.07e-24 |
| AlphaMissense | central (352, Dace) | 352 | -0.2910 | 2.69e-08 |
| BayesDel noAF | RING/BRCT (373, Findlay) | 373 | -0.3926 | 3.39e-15 |
| BayesDel noAF | central (352, Dace) | 352 | -0.0895 | 9.35e-02 |

## Fisher-z domain contrast (RING/BRCT vs central)

| predictor | Fisher z | Fisher p |
|---|---|---|
| AlphaMissense | -3.201 | 1.37e-03 |
| BayesDel noAF | -4.357 | 1.32e-05 |

## AM-class callout mismatches

- 5 P/LP controls have NO AlphaMissense score; all are BRCT/start-codon variants (p.Met1*) where the protein-language model has no FL-annotated context at the translation start (documented AM limitation for residues 1-codons).
- 13 P/LP controls are AM-classed 'likely_benign'; 4 B/LB controls are AM-classed 'likely_pathogenic'. Per variant lists written to data/processed_eval/controls_plp_am_likely_benign.tsv / controls_blb_am_likely_pathogenic.tsv.

## Caveats

- All three in silico predictors are trained on ClinVar (directly or indirectly); AUC on ClinVar labels is not independent validation.
- Controls use >=2-star review status + germline origin; conflicting/lower-tier entries excluded.
- REVEL annotation for controls depends on the 667 MB zip download (rothsj06.dmz.hpc.mssm.edu); missing -> REVEL omitted for controls (VUS REVEL reused from frozen data).
- AlphaMissense is CC BY-NC-SA 4.0 (non-commercial); acknowledgments required.
- 'g37_allele_consistent' verified True for all controls and VUS (BayesDel join safe).
