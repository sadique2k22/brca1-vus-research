# Changelog

All notable changes to this project are documented here.
Format: date — summary (author/agent).

## 2026-08-19 — Phase 0
- Environment inspected (see Phase 0 report).
- Project directory structure created.
- `README.md`, DRAFT `protocol.md`, `.gitignore`, `CHANGELOG.md`, `config/config.yaml`
  template, and `src/` module stubs created.
- Candidate genes researched via NCBI ClinVar (build 260818-0035.1) and PubMed.
- Gene comparison table produced; recommendation pending user approval.
- **No variants analyzed yet.**

## 2026-08-19 — Phase 0.1 (self-audit + fixes)
- Self-audit performed (see `results/reports/phase0_audit.md`).
- **CRITICAL fix:** corrected mis-cited PMID (Giacomelli 2018 → PMID 30224644, was 30061595).
- Protocol revised: VCEP-specified transcript (exact version); "Conflicting" vs "Uncertain"
  distinction; multi-consequence missense rule; ClinVar build verification; gene-specific
  threshold override; anti-fabrication/provenance section.

## 2026-08-19 — Phase 0.2 (full re-verification of citations & thresholds)
- All 10 cited PMIDs re-verified against PubMed E-utilities (full bibliographic detail).
- REVEL PP3 (≥0.932) / BP4 (≤0.290) confirmed from Pejaver 2022 (PMID 36413997) Table 1.
- BA1 "above 5% in ESP/1000G/ExAC" confirmed from Richards 2015 (PMID 25741868) Table 4;
  protocol wording corrected to ">5%" and BS1-vs-BA1 nuance documented.

## 2026-08-19 — Phase 0.3 (mobile-environment audit)
- Confirmed Android/Termux/proot-distro host (Cortex-X4 + 7× A720, f2fs, ~2.3 GB RAM free).
- Created `results/reports/environment_report.md`; README "Computational Environment"
  section added.
- Established resource rules: API-first, streaming/chunked, caching, 1–4 workers, no Docker,
  `/tmp` avoided (RAM-backed), dbNSFP/VEP-cache/whole-genome-gnomAD prohibited.
- Resource estimate: workflow SAFE under 2 GB RAM / 5 GB temp / 10 GB storage if streamed.
- No software installed; no large downloads.

## 2026-08-19 — Phase 0.4 (storage cap)
- User set a **50 GB total project storage cap**; recorded in `config/config.yaml`
  (`resources.storage_cap_gb: 50`) and updated `environment_report.md` + README.
- Revised dataset limits: gnomAD whole-genome / full CADD still prohibited (exceed cap);
  VEP cache & dbNSFP now "fits but discouraged"; standalone REVEL feasible within budget.
- Baseline workflow stays < 1 GB (~2% of cap).

## 2026-08-19 — Phase 1 (protocol freeze)
- Gene approved: **BRCA1**; predictors approved: **REVEL, CADD, SIFT, PolyPhen-2**.
- `protocol.md` frozen at **v1.0** (decision register, evidence-category definitions,
  statistical plan, sources table).
- Verified from primary sources: BRCA1 VCEP transcript NM_007294.4 (= MANE Select);
  BRCA1 BS1 = faf ≥ 0.001; gene-specific bioinformatic code = BayesDel ≥ 0.28.
- **Correction:** REVEL PP3 supporting threshold is **≥ 0.644** (0.932 = strong), per
  Pejaver 2022 Table 2; BP4 ≤ 0.290 confirmed.
- `config/config.yaml` frozen (no secrets); `scripts/validate_config.py` + consistency
  check report added; README "How the analysis will proceed" section added.
- No datasets downloaded.

## 2026-08-19 — Phase 2 (ClinVar retrieval & basic parsing)
- Downloaded ClinVar `variant_summary.txt.gz` (442 MB, build 260818-0035.1); md5 verified
  against official `.md5`; sha256 recorded; `metadata.json` written.
- `src/clinvar.py` streaming parser: 9,044,810 raw rows → 31,542 BRCA1 rows → 16,023 unique
  Variation IDs (GRCh37 15,818 / GRCh38 15,529 / na 195).
- Aggregate "Uncertain significance" = 2,583 (not Phase 0 API estimate 7,930 — API counted
  submission-level significance; resolves U3).
- Discovered per-assembly row duplication + deprecated `ReferenceAllele`/`AlternateAllele`
  columns (use `*VCF` columns).
- `scripts/validate_clinvar.py` (checksum + reproducibility checks) PASSED.
- Report: `results/reports/clinvar_retrieval_report.md`. No annotation/predictors run.

## 2026-08-19 — Phase 3 (ClinVar audit + VUS/missense filtering)
- `src/variants.py`: protein-change classifier (3-letter AA codes) + `filter_vus_missense`;
  `tests/test_variants.py` (12 unit tests, pass).
- `scripts/clinvar_audit.py`: 13-point audit; applies FINAL inclusion criteria.
- Final candidate set: **1,904 missense VUS** (GRCh38) after dedup (31,542→16,023), GRCh38
  (→15,529), VUS (→2,504), missense (→1,904).
- Outputs: `data/processed/clinvar_vus_missense.tsv`, `results/reports/clinvar_audit.md`.
- **Fixed .gitignore** (was silently excluding reports/metadata/processed data from git).
- No gnomAD annotation / predictors run (stopped per instructions).

## 2026-08-19 — Phase 4A (transcript & variant normalization; CI pivot)
- Transcript investigation (NCBI RefSeq, Ensembl): NM_007294.3 vs .4 = identical CDS
  (5592 nt, NP_009225.1); .3→.4 is 5' UTR-only. NM_001408514.1 = distinct isoform
  (NP_001395443.1). All 1,904 missense VUS are already on NM_007294.4 (approved).
- `src/annotation.py` (VEP REST client, cached, retry), `src/variants.py` normalization
  helpers, `scripts/transcript_normalize.py`, `tests/test_normalize.py` (21 tests pass).
- **Pivot to GitHub Actions**: `requirements.txt` + `.github/workflows/pipeline.yml`
  (download→parse→filter→normalize→commit-results + artifacts). Compute runs in CI, not
  on the phone.
- Public repo `sadique2k22/brca1-vus-research`; auto-commit results + artifacts.

## 2026-08-19 — Phase 4B (population + computational annotation)
- **STEP 0 source audit** (`results/reports/annotation_source_audit.md`): gnomAD v4 GraphQL
  (gnomad_r4) + VEP release 116 (SIFT/PolyPhen included) verified live; REVEL = standalone
  ~636 MB zip (direct URL); **CADD v1.7 online scoring degraded** (hardware failure, long
  queue, no clean bulk API) → best-effort v1.6 web only.
- New modules: `src/population.py` (gnomAD), `src/predictors.py` (REVEL + CADD),
  `src/annotation.py` (SIFT/PolyPhen extraction), `src/variants.py` (variant_key).
- `scripts/prepare_revel.py` (download + extract chr17 region),
  `scripts/annotate_variants.py` (dedup → annotate → map-back → QC/reports).
- Outputs: `biological_variant_map.tsv`, `annotation_unique_variants.tsv`,
  `brca1_vus_missense_annotated.tsv`; reports `duplicate_variant_report.md`,
  `annotation_report.md`. Tests: 30 pass. Workflow extended (REVEL cache + annotate).

## 2026-08-20 — Phase 13 (Dace 2025 integration + domain-contrast analysis)
- Integrated Dace et al. 2025 medRxiv preprint HAP1 SGE scores (central exons 6/10/11/12) for 352/1,904 VUS (18.5%); combined functional coverage 725/1,904 (38.1%). NEW files only; frozen Phase 7-12 inputs byte-identical.
- Domain-contrast analysis (`scripts/phase13_dace_analysis.py`): predictor-functional correspondence does NOT transfer to central exons — REVEL -0.384→-0.049 (Fisher p=1.9e-6, perm p=1e-4), PolyPhen-2 -0.188→+0.047 (p=1.5e-3/1.1e-3), SIFT -0.370→-0.260 (p=0.10/0.073, attenuates but directional).
- ORIENTATION AUDIT (fix): phase9.py computes SIFT as (1-SIFT); manuscript's sign rationale was wrong for SIFT (raw SIFT ρ=+0.370 matches Findlay's +0.363). Impact orientation now documented in Methods; Results and §4 corrected; feasibility report orientation mix corrected.
- Conflict analysis on the 352: 28 (8.0%) vs 23/373 (6.2%) RING/BRCT (phase9-identical REVEL rules); 30 with Dace-paper LoF threshold (-0.799); 41-cohort/13 conflicts untouched.
- Artifacts: S5_dace_352.tsv, fig17-19 (SVG+PNG), phase13_domain_contrast.md, 3 new tables; checksums/frozen dataset re-verified.
