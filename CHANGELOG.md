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
