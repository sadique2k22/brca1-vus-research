# Phase 0 Self-Audit — Skeptical Reviewer Report

Date: 2026-08-19. Reviewer role: skeptical computational genetics reviewer.

Findings are classified CRITICAL (must fix before Phase 1), IMPORTANT (must fix before the
corresponding analysis step), MINOR (cosmetic / wording). Fixes already applied are noted.

## Findings by category

### 1. Unsupported assumptions
| # | Severity | Finding | Resolution |
|---|---|---|---|
| A1 | **CRITICAL** | Transcript assumed "MANE Select by default"; VCEP-specified transcript + exact RefSeq version can differ and changes HGVS p. numbering. | Fixed — §3 now defers to VCEP transcript with exact version. |
| A2 | IMPORTANT | Assumed gnomAD "v4.1" is current without verifying. | Fixed — §10 verifies/pins version at access. |
| A3 | IMPORTANT | Assumed "≥1 star" review-status filter is appropriate. | Fixed — §6 records status, does not filter by default. |
| A4 | **CRITICAL** | Did not separate "Conflicting interpretations" from "Uncertain significance". | Fixed — §5/§8 treat them as distinct. |
| A5 | MINOR | "gold standard" overstates Findlay 2018 (single-assay readout). | Fixed in report wording. |

### 2. Outdated database/tool information
| # | Severity | Finding | Resolution |
|---|---|---|---|
| B1 | IMPORTANT | gnomAD version not pinned (see A2). | Fixed. |
| B2 | MINOR | VCEP status inferred from published specs, not verified on clingen.org roster. | Documented as inference; verify in Phase 1. |
| B3 | none | ClinVar build/date recorded correctly. | — |

### 3. Incorrect interpretation of ClinVar
| # | Severity | Finding | Resolution |
|---|---|---|---|
| C1 | IMPORTANT | One-Variant-many-conditions structure flagged but needed explicit rule. | §5 reinforced. |
| C2 | **CRITICAL** | Molecular consequence can be multi-valued ("missense" + "splice"); esearch counts therefore over-count plain missense. | Fixed — §9 multi-consequence rule. |
| C3 | IMPORTANT | "ClinVar reports GRCh38" stated without verifying column build at download. | Fixed — §4 verify, don't assume. |
| C4 | IMPORTANT | API counts presented as precise and not saved to disk. | Fixed — report caveat + §22 saves queries/responses. |

### 4. Incorrect interpretation of gnomAD
| # | Severity | Finding | Resolution |
|---|---|---|---|
| D1 | IMPORTANT | Generic BA1=5% not appropriate for dominant highly-penetrant genes (BS1 is the relevant code; gene-specific cutoffs). | Fixed — §13 override rule. |
| D2 | IMPORTANT | Population-AF confidence and ancestry imbalance understated. | Fixed — §10 ancestry caveats. |
| D3 | MINOR | "faf95" described loosely. | Fixed — §11 defined. |

### 5. Variant normalization
| # | Severity | Finding | Resolution |
|---|---|---|---|
| E1 | IMPORTANT | HGVS re-derivation must use a robust library (VEP REST / hgvs), not string parsing. | Fixed — §17. |
| E2 | IMPORTANT | ref/alt validation needs a named reference FASTA. | Fixed — §18. |
| E3 | MINOR | Multi-nucleotide variant representation not addressed. | Fixed — §18. |

### 6. Transcript selection
| # | Severity | Finding | Resolution |
|---|---|---|---|
| F1 | **CRITICAL** | Same root cause as A1 (MANE vs VCEP, version mismatch). | Fixed — §3/§17. |

### 7. Inappropriate ACMG/AMP usage
| # | Severity | Finding | Resolution |
|---|---|---|---|
| G1 | IMPORTANT | Generic Pejaver REVEL thresholds / BA1 may not match gene-specific spec. | Fixed — §13 defer to gene-specific spec. |
| G2 | IMPORTANT | Risk of summing criteria into a clinical class. | Fixed — §13 explicit "no verdict" rule. |

### 8. Population-bias issues
| # | Severity | Finding | Resolution |
|---|---|---|---|
| H1 | IMPORTANT | gnomAD ancestry imbalance; rare-variant AF wide CIs. | Fixed — §10. |
| H2 | IMPORTANT | AF↔predictor-score correlation is partly circular (conservation-trained tools). | Fixed — §19 caveat. |
| H3 | MINOR | ClinVar/founder-testing ascertainment bias in VUS set. | Fixed — §21. |

### 9. Reproducibility problems
| # | Severity | Finding | Resolution |
|---|---|---|---|
| I1 | IMPORTANT | No `requirements.txt` / venv pinning yet. | Scheduled for Phase 1. |
| I2 | IMPORTANT | Phase 0 count queries/responses not saved. | Fixed — §22. |
| I3 | IMPORTANT | Must freeze raw files + checksum at Phase 1 start. | §22 + download log ready. |
| I4 | MINOR | Notebooks non-reproducible if authoritative. | §22 scripts authoritative. |

### 10. Fabrication / inference risk
| # | Severity | Finding | Resolution |
|---|---|---|---|
| J1 | **CRITICAL** | **Incorrect PMID** `30061595` cited for Giacomelli 2018 (that PMID is an IBD paper; correct is `30224644`). | **Fixed** in report; verified via E-utilities. |
| J2 | IMPORTANT | Predictor scores must never be model-generated. | Fixed — §12 hard rule. |
| J3 | IMPORTANT | Functional dataset values must be downloaded, not approximated. | Fixed — §23. |
| J4 | IMPORTANT | All PMIDs/DOIs resolved via E-utilities + logged. | Fixed — §14/§23. |

## Summary

- **CRITICAL (4):** A1/F1 transcript selection; A4 conflicting-vs-uncertain; C2 multi-consequence
  missense; J1 incorrect PMID. **All fixed.**
- **IMPORTANT:** 16 items; fixed in the revised protocol where protocol-level, remaining are
  scheduled for their Phase 1 step (requirements.txt, clingen.org roster verification).
- **MINOR:** 5 items; wording/documentation corrections applied or noted.

Phase 1 remains blocked pending gene-selection approval, but the protocol and Phase 0 report
are now defensible.
