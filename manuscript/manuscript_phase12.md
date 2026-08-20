# A reproducible framework for integrating population, computational, and functional evidence for BRCA1 missense variants of uncertain significance

*Authors: [to be completed]*

*Affiliations: [to be completed]*

---

## Abstract

**Background.** Missense variants of uncertain significance (VUS) in *BRCA1* are a persistent challenge for genetic counseling. Computational pathogenicity predictors and population-frequency data are widely used to prioritize such variants, but their correspondence with experimental functional evidence is incompletely characterized.

**Methods.** We performed a retrospective computational and evidence-integration study of 1,904 *BRCA1* missense VUS from ClinVar (GRCh38; transcript NM_007294.4). Each variant was annotated with gnomAD v4 population frequency and REVEL, SIFT, and PolyPhen-2 scores, and with Findlay *et al.* (2018) saturation genome-editing (HAP1 cellular-fitness) functional scores where available. A 41-variant descriptive evidence cohort was selected deterministically from predefined computational/functional strata.

**Results.** Of 1,904 variants, 424 (22.3%) were present in gnomAD and 1,480 (77.7%) were absent; present variants were themselves extremely rare (median allele frequency 6.8×10⁻⁷). REVEL, SIFT, and PolyPhen-2 showed only moderate pairwise agreement (Cohen's κ = 0.27–0.45). Among the **373 variants with available Findlay HAP1 functional measurements (RING and BRCT domains only)**, computational predictors showed moderate correspondence with the functional score: Spearman ρ = −0.384 (95% CI −0.468 to −0.294) for REVEL, −0.370 (−0.454 to −0.278) for SIFT, and −0.188 (−0.284 to −0.088) for PolyPhen-2. Within the 41-variant descriptive cohort, 13 variants showed conflict between computational and functional evidence.

**Conclusions.** Computational predictors provide useful prioritization information but correspond only moderately with this experimental functional measurement, with substantial variant-level discordance. The primary contribution of this work is a reproducible, provenance-preserving framework for integrating these evidence types, not a clinical reclassification of *BRCA1* variants.

**Keywords.** BRCA1; variants of uncertain significance; missense variant; computational prediction; population frequency; saturation genome editing; functional assay; reproducibility.

---

## 1. Introduction

Germline loss-of-function variants in *BRCA1* predispose carriers to breast and ovarian cancer, and *BRCA1* is among the most frequently sequenced cancer-predisposition genes. A substantial fraction of observed *BRCA1* variants are missense variants of uncertain significance (VUS), for which the associated risk cannot be definitively assigned. The resulting clinical uncertainty limits the utility of genetic testing.

Several classes of evidence are used to interpret missense VUS. *Population-frequency* evidence (e.g., from gnomAD) can indicate that a variant is too common to be a highly penetrant pathogenic variant, but rarity is not itself evidence of pathogenicity. *Computational* missense predictors (e.g., REVEL, SIFT, PolyPhen-2) estimate the likely impact of an amino-acid substitution, but they are imperfect proxies and are themselves correlated (REVEL is an ensemble meta-predictor). *Experimental functional* assays provide a more direct readout of variant effect; for *BRCA1*, saturation genome editing (Findlay *et al.*, 2018) produced cellular-fitness scores for thousands of variants in CRISPR-edited HAP1 cells. Computational prediction, population frequency, and experimental functional evidence are distinct and should not be conflated.

The central question addressed here is: **how consistently do population-frequency evidence and computational missense-variant predictors correspond with available experimental functional evidence for *BRCA1* missense VUS?**

**Positioning and novelty.** The comparison of computational predictors against multiplexed functional data is established prior work (e.g., *In silico* prediction tools for variant curation in cancer genes; ACMG/AMP interpretation of *BRCA1* missense variants with structure-informed scores; RING-domain *BRCA1* VUS assessment). The present study does **not** claim novelty for the predictors, the Findlay assay, or the general predictor–functional comparison. Its contribution is **methodological and framework-level**: a reproducible, provenance-preserving, deterministic, end-to-end pipeline that traces a defined public VUS cohort from ClinVar through normalization, population annotation, computational prediction, functional-data integration, deterministic candidate selection, and evidence-source verification, **without performing clinical reclassification**. The findings are descriptive and are not intended to reclassify any variant.

---

## 2. Methods

### Study design

This is a retrospective, observational computational and evidence-integration study using exclusively public data. No primary clinical or experimental work was performed, and no variant was reclassified.

### ClinVar dataset

*BRCA1* variants were obtained from the NCBI ClinVar `variant_summary.txt.gz` (release build 260818-0035.1, retrieved 2026-08-19; SHA-256 `d2a8c9c2…`). The **aggregate** germline clinical significance was used rather than submission-level labels, because a single ClinVar Variation record can carry multiple submissions and multiple conditions; submission-level "uncertain significance" counts therefore overestimate the number of variants with an aggregate VUS classification. Variants were included if the aggregate germline significance was exactly "Uncertain significance" and the molecular consequence was a single-nucleotide amino-acid substitution (missense), identified from the ClinVar protein-change notation together with variant type "single nucleotide variant". Variants with "Conflicting interpretations of pathogenicity" were tracked separately and excluded from the VUS set.

### Variant normalization

All variants were mapped to GRCh38. HGVS cDNA/protein consequences were re-derived from genomic coordinates using Ensembl VEP (release 116) and cross-validated against ClinVar annotations (two-pass validation; 0 discrepancies). The analysis transcript was **NM_007294.4** (MANE Select; ENST00000357654), which is also the ENIGMA/ClinGen BRCA1 reference transcript. NM_007294.3 (used by the Findlay dataset) and NM_007294.4 have identical coding sequence (CDS), so cDNA numbering is directly comparable.

### gnomAD

Population-frequency data were retrieved from the gnomAD v4 GraphQL API (dataset `gnomad_r4`, GRCh38), recording global allele frequency, allele count/number, homozygote count, and filtering allele frequency (faf95) population maximum. Exome and genome datasets were recorded separately. **Absence from gnomAD was treated as a category (absent), never as an allele frequency of zero.**

### Computational predictors

REVEL, SIFT, and PolyPhen-2 (HumVar) scores were obtained from Ensembl VEP / REVEL, using documented thresholds for neutral descriptive categorization: REVEL ≤ 0.290 (tolerance-supporting) and ≥ 0.644 (impact-supporting) per Pejaver *et al.* (2022); SIFT ≤ 0.05 (deleterious) per Ng & Henikoff (2003); PolyPhen-2 ≤ 0.446 (benign) and > 0.908 (probably damaging) per Adzhubei *et al.* (2010). REVEL is an ensemble predictor and the three predictors are correlated; they were therefore treated as correlated predictors, **not** as independent lines of evidence. The ENIGMA/ClinGen BRCA1/BRCA2 VCEP applies its own bioinformatic code using BayesDel; the generic REVEL thresholds were used descriptively only.

**Predictor-panel limitations.** CADD was excluded because its v1.7 bulk-annotation pathway was unavailable during this study. **AlphaMissense** (Cheng *et al.* 2023, PMID 37733863) was **not** included: integrating it would require re-annotating and re-running the frozen analysis, and its non-commercial license constrains redistribution; it is noted as a planned extension (see Limitations).

### Findlay functional dataset

The Findlay *et al.* (2018) *BRCA1* saturation genome-editing functional scores were retrieved from MaveDB (URN `urn:mavedb:00000097-0-2`, "BRCA1 SGE Normalized Scores"). This assay measures **HAP1 cellular fitness** after CRISPR-based editing of the endogenous *BRCA1* locus; lower (more negative) scores indicate reduced cellular fitness. The dataset comprises 3,893 SNVs across **13 exons** encoding the RING and BRCT domains; **it does not cover the DNA-binding domain (exons 6–14)**. Scores were mapped to our variants by cDNA identity (never by protein change alone). Findlay *et al.* classified variants using a two-component Gaussian mixture model with posterior-probability thresholds (P(non-functional) > 0.99 / < 0.01); because the MaveDB score set provides only the **continuous** score (not the mixture-model classification), we use the continuous score as the primary functional variable and do **not** apply a binary threshold. A negative score is interpreted as *Findlay HAP1 cellular-fitness evidence consistent with reduced cellular fitness*, and is **not** a direct clinical pathogenicity measurement.

### Statistical analysis

Descriptive statistics; pairwise categorical agreement (Cohen's κ); Spearman rank correlation with 95% confidence intervals; Mann-Whitney U with rank-biserial effect size for two-group comparisons. All analyses are **exploratory**; no multiple-testing correction was applied, and p-values are reported descriptively and accompanied by effect sizes, not presented as standalone findings. ROC-AUC was **not** computed because a validated binary functional threshold was not available.

### Candidate selection

The analysis proceeded as: 1,904 total VUS → Phase 5 computational/population pattern analysis → 436-variant candidate union → 41-variant final descriptive evidence cohort. The final cohort was selected **deterministically** from five strata (computational impact/tolerance × functional non-functional/WT-like, plus intermediate), with target 10 per stratum and even-spacing by protein position to preserve regional diversity. Selection used only pre-existing computational and functional information (frozen in `configs/phase7_cohort.yaml`); **literature availability was not a selection criterion**. The 41-variant cohort is **descriptive and not statistically powered for population-level inference**.

### Literature search and verification

PubMed was searched for each prioritized variant using `BRCA1 <protein change>` and `BRCA1 <cDNA change>` queries (E-utilities), recording PMID/title/year/journal. In Phase 9, each retrieved record was **verified at the abstract level** against the exact variant identifiers (cDNA and protein change), classifying records as exact-variant, gene-level, or unclear. Searches are abstract-level, **not full-text**; exact-variant evidence that appears only in full text or supplementary tables was therefore not captured. All PMIDs/DOIs originate from E-utilities responses; no citation was fabricated.

---

## 3. Results

### Dataset characteristics

The analysis set comprised **1,904 *BRCA1* missense VUS** (GRCh38; NM_007294.4). Of these, **424 (22.3%) were present** in gnomAD v4 and **1,480 (77.7%) were absent**. Present variants were themselves extremely rare (median global allele frequency 6.8×10⁻⁷; maximum 2.6×10⁻⁵); the present-vs-absent comparison is therefore effectively an ultra-rare-versus-absent comparison. No VUS had a filtering allele frequency (faf95) population maximum ≥ 0.001.

### Computational predictions

REVEL categorized 367 variants as tolerance-supporting, 1,246 as intermediate, and 291 as impact-supporting. SIFT categorized 1,071 as deleterious and 833 as tolerated. PolyPhen-2 categorized 1,328 as benign, 412 as possibly damaging, and 164 as probably damaging. Pairwise categorical agreement was moderate: Cohen's κ = 0.448 (REVEL–SIFT), 0.441 (REVEL–PolyPhen-2), and 0.272 (SIFT–PolyPhen-2). Continuous-score correlations among predictors were weak-to-moderate (Spearman ρ = 0.29–0.41).

### Population versus computational evidence

Predictor scores did not materially differ between gnomAD-present and gnomAD-absent variants (Mann-Whitney: REVEL p = 0.89; SIFT p = 0.44; PolyPhen-2 p = 0.014). The PolyPhen-2 difference is **uncorrected exploratory evidence** and should not be interpreted as a standalone finding given the multiple tests performed and the small effect. Among gnomAD-present variants, allele frequency did not correlate with predictor scores (Spearman ρ ≈ 0).

### Functional evidence

Findlay HAP1 cellular-fitness scores were available for **373 of 1,904 variants (19.6%)**, reflecting the RING and BRCT scope of the assay (124 RING and 236 BRCT scored; the remaining variants lie in regions not assayed, including the DNA-binding domain).

### Computational–functional correspondence (RING/BRCT subset)

Among the **373 RING/BRCT variants with Findlay HAP1 functional measurements**, computational predictors showed **moderate** correspondence with the functional score: Spearman ρ = −0.384 (95% CI −0.468 to −0.294) for REVEL; −0.370 (−0.454 to −0.278) for SIFT; and −0.188 (−0.284 to −0.088) for PolyPhen-2. The negative sign reflects the opposite orientation of the Findlay score (lower = reduced cellular fitness) relative to the predictors (higher = predicted impact). These results apply only to the measured RING/BRCT subset and cannot be generalized to all 1,904 variants or to the unmeasured DNA-binding region.

### Evidence conflicts

Within the 41-variant descriptive cohort, **13 variants (32%)** showed conflict between computational and functional evidence when assessed using the continuous functional score (computational impact-supporting together with a WT-like functional score, or computational tolerance-supporting together with a clearly negative functional score). This is a descriptive count, not an inferential estimate.

### Evidence-source verification

Of the 41 cohort variants, 36 returned at least one PubMed record (295 records total); however, **abstract-level verification identified no record whose abstract explicitly discusses the exact variant** (all retrieved records were gene-level *BRCA1* papers). No exact-variant publication was identified under the documented abstract/metadata-level verification procedure; full-text/supplementary verification was outside the automated scope. None of the 41 variants had ClinVar expert-panel curation at the time of analysis (all remained aggregate VUS).

---

## 4. Discussion

Three descriptive findings emerge. First, the three computational predictors agree only moderately with one another on *BRCA1* missense VUS (κ = 0.27–0.45), so a variant's categorization is highly tool-dependent. Second, computational predictors correspond only moderately with the Findlay HAP1 functional measurement within the measured RING/BRCT subset (|ρ| ≈ 0.19–0.38), and weakly for PolyPhen-2. Third, a substantial minority of the descriptive evidence cohort (13/41) exhibited discordance between computational and functional evidence.

These results do **not** imply that computational predictors are "inaccurate." Rather, their correspondence with this particular functional assay is moderate and heterogeneous, and their mutual dependence (REVEL is an ensemble; conservation-trained predictors share a common signal, introducing circularity with population/evolutionary data) limits their value as independent evidence.

Population rarity was largely uninformative in this VUS set: most variants were ultra-rare or absent, and no variant crossed the 0.001 filtering-AF reference point.

The observed functional-score distributions differed numerically between the sampled RING and BRCT variants (median −0.52 vs. −0.30), but this difference did not reach statistical significance (Mann–Whitney p = 0.059; rank-biserial r = 0.12, a small effect) and is confounded by variant composition and unequal coverage. It is reported descriptively and is not interpreted as a difference in domain "pathogenicity."

A further caveat concerns the literature: although most cohort variants returned PubMed records, abstract-level verification found no record explicitly discussing the exact variant, so those records are gene-level rather than variant-specific. This does not imply that no variant-specific evidence exists — only that it would reside in full text or supplementary data beyond the scope of this review.

The Findlay assay is a powerful but specific readout: it measures HAP1 cellular fitness across the RING and BRCT domains, does not cover the DNA-binding domain, and cannot capture every mechanism (e.g., splicing or tissue-specific effects). The 13 conflicts highlight variants where computational and functional evidence diverge and where neither source should be treated as decisive.

### Relationship to prior and concurrent analyses

This study is a re-analysis of public data rather than a new biological experiment, so its findings must be read alongside the numbers already reported by the producers of the functional data. Findlay *et al.* (2018; ref 7) themselves reported, in their supplementary analyses (Extended Data Fig. 9), that missense SGE function scores correlated with SIFT (Spearman ρ = 0.363) and PolyPhen-2 (ρ = −0.277), with P < 1 × 10⁻³⁷ for all correlations, using their own annotation pipeline. Our independently computed estimates for the RING/BRCT subset (SIFT ρ = −0.370; PolyPhen-2 ρ = −0.188) are concordant in magnitude with those values once score orientation is normalized (the sign difference reflects the opposite orientation of the functional score, lower = reduced cellular fitness, relative to the predictors, higher = predicted impact). We therefore do not claim the correspondence itself as new. What this work adds is: (i) a fully reproducible pipeline from ClinVar through to the final descriptive cohort, with provenance, checksums, and tests; (ii) the clinically relevant VUS subset with allele-frequency context from gnomAD v4, rather than all assayed variants; (iii) REVEL, which Findlay *et al.* did not evaluate; (iv) confidence intervals and rank-biserial effect sizes for every comparison; (v) a deterministic 41-variant descriptive evidence cohort and its 13 computational–functional conflicts; and (vi) an explicit decision to make no clinical classification claims. Concurrently, a successor from the same laboratory has extended HAP1 SGE to 11 regions outside the RING and BRCT domains (Dace *et al.*, 2025; 4,113 newly assayed variants) and introduced a mammary-epithelial (HMEC) assay; we cite this work as the immediate extension path and did not integrate it here because it is a preprint whose integration would require re-running the frozen analysis (see Limitations).

**Contribution.** The primary contribution of this work is not a new variant-classification algorithm or a new biological classification of *BRCA1* variants. It is a reproducible, auditable, provenance-preserving, deterministic, end-to-end framework for tracing a defined public VUS cohort through normalization, population annotation, computational prediction, functional-data integration, deterministic candidate selection, and evidence-source verification, without performing clinical reclassification.

---

## 5. Limitations

This study has important limitations: ClinVar ascertainment bias (submitted variants are not a random sample); instability of VUS classification over time; gnomAD population representation and ancestry imbalance (and the present-variant group being itself ultra-rare); dependence and training-set circularity among computational predictors; the exclusion of CADD and AlphaMissense; the Findlay assay's HAP1 single-readout design and RING/BRCT-only scope: at the time of analysis the only peer-reviewed saturation-scale functional resource for *BRCA1* covered the RING and BRCT domains, with no measurements in the central/DNA-binding exons; a 2025 preprint has since expanded HAP1 SGE coverage to exons 6, 10, 11 and 12 (Dace *et al.*, ref 13), which is the natural extension of this work but was not integrated here because it is not yet peer-reviewed and integration would require re-running the frozen analysis. Consequently the functional-correlation results apply only to the measured RING/BRCT subset and cannot be generalized to all 1,904 variants or to the DNA-binding region; abstract-level literature search (no full-text review, no systematic clinical phenotype or segregation extraction); the absence of expert-panel curation among the cohort; exploratory statistics without multiple-testing correction; the descriptive, non-powered nature of the 41-variant cohort; and the absence of any clinical validation or ACMG reclassification. Accordingly, no variant is reported as pathogenic or benign, and the 41-variant cohort is not used to estimate predictor performance for all *BRCA1* VUS.

---

## 6. Conclusion

This study demonstrates a reproducible framework for integrating population, computational, functional, and literature evidence around *BRCA1* missense VUS. Computational predictors showed only moderate correspondence with available experimental functional measurements and substantial variant-level discordance. These findings support using computational predictions as prioritization tools rather than as standalone determinants of variant interpretation. The study does not reclassify any variant and does not establish pathogenicity or benignity.

---

## 7. Data and Code Availability

All raw data originate from public resources (NCBI ClinVar, gnomAD, Ensembl VEP, MaveDB, PubMed). Derived datasets, code, configuration, and generated results are available in the project repository (`github.com/sadique2k22/brca1-vus-research`); exact versions, retrieval dates, and SHA-256 checksums are documented in `data/raw/*/metadata.json` and the reports directory. The full pipeline is executed on GitHub Actions (`.github/workflows/pipeline.yml`) and is reproducible from the committed code and configuration. No external repository deposition has been performed.

---

## 8. Author Contributions

[placeholder — to be completed]

## 9. Conflict of Interest

[placeholder — to be completed]

## 10. Funding

[placeholder — to be completed]

---

## 11. References

1. Landrum MJ, Lee JM, Benson M, et al. ClinVar: improving access to variant interpretations and supporting evidence. Nucleic Acids Res. 2018;46(D1):D1062–D1067. PMID 29165669.
2. Chen S, Francioli LC, Goodrich JK, et al. A genomic mutational constraint map using variation in 76,156 human genomes. Nature. 2024;625(7993):92–100. PMID 38057664.
3. Karczewski KJ, Francioli LC, Tiao G, et al. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature. 2020;581(7809):434–443. PMID 32461654.
4. Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants. Genet Med. 2015;17(5):405–424. PMID 25741868.
5. Pejaver V, Byrne AB, Feng BJ, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. Am J Hum Genet. 2022;109(12):2163–2177. PMID 36413997.
6. Parsons MT, de la Hoya M, Richardson ME, et al. Evidence-based recommendations for gene-specific ACMG/AMP variant classification from the ClinGen ENIGMA BRCA1 and BRCA2 Variant Curation Expert Panel. Am J Hum Genet. 2024;111(9):2044–2058. PMID 39142283.
7. Findlay GM, Daza RM, Martin B, et al. Accurate classification of BRCA1 variants with saturation genome editing. Nature. 2018;562(7726):217–222. PMID 30209399.
8. Ioannidis NM, Rothstein JH, Pejaver V, et al. REVEL: an ensemble method for predicting the pathogenicity of rare missense variants. Am J Hum Genet. 2016;99(4):877–885. PMID 27666373.
9. Kircher M, Witten DM, Jain P, O'Roak BJ, Cooper GM, Shendure J. A general framework for estimating the relative pathogenicity of human genetic variants. Nat Genet. 2014;46(3):310–315. PMID 24487276.
10. Ng PC, Henikoff S. SIFT: predicting amino acid changes that affect protein function. Nucleic Acids Res. 2003;31(13):3812–3814. PMID 12824425.
11. Adzhubei IA, Schmidt S, Peshkin L, et al. A method and server for predicting damaging missense mutations. Nat Methods. 2010;7(4):248–249. PMID 20354512.
12. Cheng J, Novati G, Pan J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. Science. 2023;381(6664):eadg7492. PMID 37733863.
13. Dace P, Forrester NM, Zanti M, et al. Saturation genome editing of BRCA1 across cell types accurately resolves cancer risk. medRxiv preprint. 2025. doi:10.1101/2025.08.11.25333423. (Preprint, not peer-reviewed.)

*Database/tool resources (no PMID):* Ensembl VEP release 116; gnomAD v4; MaveDB (URN `urn:mavedb:00000097-0-2`).

---

## 12. Figure Legends

See `manuscript/figure_mapping.md`.

## 13. Table Legends

- **Table 1.** Dataset characteristics (n = 1,904 *BRCA1* missense VUS).
- **Table 2.** Computational predictor distributions and pairwise agreement.
- **Table 3.** Final 41-variant descriptive evidence cohort.
- **Table 4.** Representative computational–functional conflicts.
