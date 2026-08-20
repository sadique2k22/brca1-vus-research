# Phase 13 Feasibility — Dace 2025 SGE data for non-RING/BRCT regions

Date: 2026-08-20. Status: **FEASIBLE — 352 additional VUS gain functional scores (18.5%).**

## Sources

- Dace P, Forrester NM, Zanti M, et al. *Saturation genome editing of BRCA1 across cell types
  accurately resolves cancer risk.* medRxiv 2026.08.11.25333423 (posted 2026-08-16; **preprint,
  not peer-reviewed**; CC-BY 4.0).
- Data/code: `github.com/phoebedace/BRCA1_SGE_HAP1_HMEC` (MIT/CC-BY per repository).
- File used: `data_for_plots/HAP1_variants_unfiltered.csv` (GRCh38 positions; r1/r2 function scores).

## Method

Our 1,904 GRCh38-bi-allelic VUS (`S1_all_1904_variants.tsv`) were intersected with the Dace HAP1
variant table by exact (chrom, pos, ref, alt). Duplicate rows (variants at the exon 12 5'/3' region
boundary) were deduplicated by taking the mean of replicate scores.

## Results

| Metric | Value |
|---|---|
| Our VUS with a Dace score (NEW regions only) | **352 / 1,904 (18.5%)** |
| Regions represented | exon 6 (81), exon 12 5' (67), exon 12 3' (67), exon 11 (58), exon 10 5' (29), exon 10 mid1 (23), exon 10 mid2 (26), exon 10 3' (14) |
| Findlay 2018 (RING/BRCT) | 373 / 1,904 (19.6%) — unchanged |
| Combined functional coverage | **725 / 1,904 (38.1%)** (was 19.6%) |
| Coverage of the 13 conflict cohort | unchanged (all remain RING/BRCT or unassayed; none gain new data) |

### Correspondence in the NEW regions (mean replicate Dace score)

| Predictor | Spearman rho (95% CI) | n |
|---|---|---|
| REVEL | −0.049 (−0.153 to 0.055) | 352 |
| SIFT | +0.260 (+0.160 to +0.355) | 352 |
| PolyPhen-2 | +0.047 (−0.057 to +0.151) | 352 |

For comparison, in the RING/BRCT subset (Findlay): REVEL −0.384, SIFT −0.370, PolyPhen-2 −0.188.
The correspondence **does not transfer** to the newly assayed central/coiled-coil-adjacent exons:
REVEL and PolyPhen-2 lose essentially all association, and SIFT is oppositely signed.

## Interpretation

- This is a **structured, potentially novel finding**: predictor-versus-function correspondence is
  domain-specific — weak-to-moderate in RING/BRCT (where most known pathogenic missense cluster) and
  near-zero for REVEL/PolyPhen-2 in the newly assayed exons, where the majority of VUS actually lie.
- It directly answers the "predictor domain bias" question and converts the project's biggest named
  limitation into a differentiator.
- The correlation contrast is the analytic core that would be new; the 13-conflict cohort and the
  41-cohort remain untouched (they are RING/BRCT or unassayed).

## Caveats / prerequisites for Phase 13

1. Preprint data — not peer-reviewed; use final filtered score tables (S1/S2) for binary classes and
   this unfiltered continuous table for correlations; document provenance + license terms.
2. Confirm GRCh38 coordinate consistency against our VEP annotations for all 352 (spot-checked OK).
3. Adding BayesDel and/or AlphaMissense is a separate decision (re-annotation of the panel).
4. Any integration changes the frozen dataset (`7afe54db…`) → run as a **new phase, not an edit**.

## Artifacts

- `data/processed_eval/vus_with_dace_scores_352.tsv` — our 352 VUS + Dace replicate scores
- This report.