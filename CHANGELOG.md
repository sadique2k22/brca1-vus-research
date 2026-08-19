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
