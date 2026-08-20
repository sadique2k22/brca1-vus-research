# Submission Checklist

## Scientific
- [x] scientific audit GREEN (`pre_manuscript_scientific_audit.md`)
- [x] peer-review issues addressed (`phase11_peer_review_response.md`)
- [x] **Findlay-engagement fix**: manuscript now cites Findlay's Extended Data Fig. 9 correlations (SIFT ρ=0.363, PolyPhen-2 ρ=-0.277) and states what this analysis adds
- [x] **Dace 2025 preprint integrated as Phase 13**: 352/1,904 VUS gain central-exon functional scores; domain-contrast analysis (REVEL/PolyPhen-2 correspondence does not transfer; SIFT attenuates); preprint-grade caveats documented
- [x] **orientation audit**: SIFT reported on a single documented impact orientation (1−SIFT) with raw-SIFT (+0.370 vs Findlay +0.363) alongside; sign rationale corrected in Results and §4; Methods document the inversion
- [x] no unsupported claims (numerical + citation audits pass; "11.1%" stat verified absent from repo; not carried into any manuscript)
- [x] frozen dataset verified (`7afe54db…`)
- [x] 41-variant cohort verified (unchanged)

## Manuscript
- [x] title
- [x] abstract
- [x] methods
- [x] results
- [x] discussion
- [x] limitations
- [x] references (13, verified)
- [x] figures (fig1–fig19, SVG + 300-dpi PNG)
- [x] tables (4, with legends)

## Supplementary
- [x] S1 (1,904 variants)
- [x] S2 (373 Findlay)
- [x] S3 (41 evidence matrix)
- [x] S4 (literature verification)
- [x] S5 (352 Dace 2025, Phase 13)
- [x] supplementary methods

## Repository
- [x] README (finished-project)
- [x] reproducibility guide (`REPRODUCIBILITY.md`)
- [x] provenance (`DATA_PROVENANCE.md`)
- [x] citation metadata (`CITATION.cff`)
- [x] license decision (documented; requires author choice)
- [x] CI (path-triggered, cached, idempotent)
- [x] tests (62/62)
- [x] no secrets / API keys / credentials

## Manual information still required
- [ ] final author names
- [ ] affiliations
- [ ] author contributions
- [ ] funding statement
- [ ] conflict-of-interest statement
- [ ] corresponding author
- [ ] license selection
- [ ] journal selection
- [ ] journal-specific formatting
