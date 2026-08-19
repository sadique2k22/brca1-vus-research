# Protocol Consistency Check (v1.0, frozen)

Date: 2026-08-19. Rule: every parameter used by the planned pipeline is either defined in
`protocol.md` (§N) or marked N/A / deferred (§1.4 U1–U5).

| Parameter (config.yaml) | Defined in protocol.md | Status |
|---|---|---|
| protocol.version / date_frozen | §0 | ✅ |
| study.gene | §3 (D1) | ✅ |
| study.disease | §3 | ✅ |
| study.inheritance | §3 | ✅ |
| study.genome_build | §4 | ✅ |
| study.transcript / transcript_ensembl / genomic_ref | §4 | ✅ |
| clinvar.source_url / build / file | §5.1 | ✅ |
| clinvar.access_date | §5.1 | ✅ |
| clinvar.include_significance / include_consequence | §6, §7 | ✅ |
| clinvar.exclude_conflicting | §6, §7 | ✅ |
| gnomad.version | §5.2 (U1) | ✅ deferred-pin |
| gnomad.datasets / fields | §5.2 | ✅ |
| predictors[].name/source/range/higher_is_damaging | §5.3 | ✅ |
| thresholds.ba1_af | §9 | ✅ |
| thresholds.brca1_bs1_faf | §9 | ✅ |
| thresholds.revel_bp4 / revel_pp3 / revel_pp3_strong | §9 | ✅ |
| thresholds.enigma_bayesdel_pp3 | §9 | ✅ (reference) |
| evidence.min_predictors / frequent_faf | §7.1 (A1, A2) | ✅ |
| resources.* (storage/ram/temp/workers) | protocol §15 + environment report | ✅ |
| output_dirs.* | README layout + §15 | ✅ |

## Deferred parameters (not N/A — resolved at execution)
| ID | Parameter | Resolution rule | protocol.md |
|---|---|---|---|
| U1 | gnomAD exact release | verify + pin at access | §1.4, §5.2 |
| U2 | predictor score source/version | pin at annotation | §1.4, §5.3 |
| U3 | exact BRCA1 missense VUS count | parsed table, not API | §1.4, §5.1 |
| U4 | CADD/SIFT/PolyPhen calibrated bins | transcribe Pejaver Table 2 | §1.4, §9 |
| U5 | OMIM IDs | verify at literature step | §1.4, §3 |

## Not applicable (explicitly out of scope)
- API keys / secrets → excluded by policy (`config.yaml` note).
- Whole-genome gnomAD / dbNSFP / VEP cache / full CADD → prohibited (§5.2, environment report).
- Clinical classification verdicts → out of scope (§2, §10, §16).

## Automated check
`scripts/validate_config.py` (exit 0 = valid) validates the config schema above.
