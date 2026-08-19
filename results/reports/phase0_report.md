# Phase 0 Report — Environment Inventory & Gene Selection

Date: 2026-08-19

## 1. Environment inventory

| Item | Status |
|---|---|
| Operating system | Ubuntu 25.10 (aarch64 / ARM64), kernel 6.17.0 |
| Python | 3.13.7 (active venv at `/root/frida-venv`) |
| CPU / RAM / Disk | 6 cores / 10 GiB RAM (~2.3 GiB free) + 11 GiB swap / 177 GiB free |
| git | 2.51.0 ✅ |
| curl / wget | 8.14.1 ✅ / 1.25.0 ✅ |
| Docker | ❌ not installed |
| VEP | ❌ not installed |
| bcftools / tabix / samtools / bgzip | ❌ not installed |
| Internet (NCBI eutils, ClinVar FTP, gnomAD API, Ensembl REST) | ✅ all reachable |

Key Python packages (in active venv): pandas 3.0.5, numpy 2.5.2, scipy 1.18.0,
matplotlib 3.11.1, biopython (Bio) 1.88, requests, openpyxl, PyYAML, tabulate, jupyter.
Not present: seaborn, statsmodels, requests-cache, scikit-learn, pysam.

**Note:** the active venv (`frida-venv`) is a general-purpose environment (frida, telethon,
whisper, video tooling). A dedicated project venv is recommended for reproducibility.

## 2. Proposed minimal environment (to be confirmed — not installed yet)

- `python3 -m venv .venv` in the project root.
- Core: `pandas numpy scipy matplotlib seaborn requests requests-cache openpyxl PyYAML`.
- Optional (Phase 2+): `statsmodels` (regression), `scikit-learn` (ROC-AUC),
  `pyensembl` or Ensembl REST client.
- VEP: **not installed locally.** Options (propose before acting): (a) Dockerized VEP
  (Docker unavailable → no), (b) conda/micromamba VEP, (c) Ensembl **VEP REST** API or
  precomputed CADD/gnomAD VCFs consumed via `tabix` (tabix not installed). We will prefer
  precomputed scores + VEP REST over a heavy local VEP cache unless justified.

## 3. Candidate gene comparison

ClinVar counts are **Variation-record** counts retrieved 2026-08-19 via NCBI eutils
(build 260818-0035.1). "VUS" = aggregate germline "Uncertain significance"; "missense"
= molecular consequence "missense variant". Note: total < sum of significance categories
because one Variation record can be linked to multiple conditions/significances.

> **Caveat:** these are *indexed search-result counts* (approximate, for gene selection
> only) — they are **not** a curated parse of `variant_summary.txt`. They may be affected
> by multi-consequence indexing (e.g. "missense + splice") and by "Conflicting
> interpretations" records. They will be superseded by a parsed, validated table in Phase 1.

| Gene | Disease association | VUS (all) | Missense VUS | Gene-specific ACMG/AMP spec | Saturation genome editing (SGE) functional map | gnomAD |
|---|---|---|---|---|---|---|
| **ATM** | Ataxia-telangiectasia (biallelic); monoallelic breast-cancer risk (magnitude debated) | 12,118 | 8,498 | Yes — ClinGen ATM VCEP, PMID 39317201 (2024) | No whole-gene SGE | Yes |
| **TP53** | Li-Fraumeni syndrome (AD) | 2,267 | 1,343 | Yes — TP53 VCEP, PMID 33300245 (2021) | Yes — Kotler 2018 (PMID 29979965); Giacomelli 2018 (PMID 30224644) | Yes |
| **BRCA1** | Hereditary breast/ovarian cancer | 7,930 | 5,027 | Yes — ENIGMA/ClinGen, PMID 39142283 (2024) | Yes — Findlay 2018 (PMID 30209399; whole-gene SGE, single-assay readout) | Yes |
| **BRCA2** | Hereditary breast/ovarian cancer | 13,463 | 9,398 | Yes — ENIGMA/ClinGen, PMID 39142283 (2024) | Yes — two SGE studies (PMID 42110911, 2026) | Yes |
| **MLH1** | Lynch syndrome | 3,690 | 2,251 | InSiGHT 5-tier, PMID 24362816 (2014) | No whole-gene SGE | Yes |
| **MSH2** | Lynch syndrome | 4,965 | 3,105 | InSiGHT 5-tier, PMID 24362816 (2014) | No whole-gene SGE | Yes |

## 4. Recommendation

**Primary recommendation: BRCA1.** Runner-up: TP53.

Rationale (not "many VUS"):
1. **Functional ground truth** — the project's objective #8/13/17 is to compare
   computational predictions against *experimental* evidence. BRCA1 has the best-established,
   freely downloadable functional map (Findlay 2018 saturation genome editing, PMID 30209399),
   enabling predictor benchmarking and ROC-AUC.
2. **Modern gene-specific ACMG/AMP specifications** (ENIGMA/ClinGen, PMID 39142283) permit a
   faithful evidence-framework discussion (PP3/BP4, BA1/BS1) without inventing rules.
3. **Adequate sample size** (5,027 missense VUS) for stable statistics, without being
   unwieldy.
4. **Extensive VUS-reclassification literature** allows the "compare against prior
   reclassifications" analysis.

TP53 is a strong alternative (excellent multi-assay functional data and a clean VCEP spec,
but fewer missense VUS, ~1,343, and more nuanced functional readouts — dominant-negative /
gain-of-function — which raises methodological complexity for a first pass).

ATM/MLH1/MSH2 were set aside because they lack whole-gene saturation functional maps, which
weakens the "predictions vs. experimental evidence" comparison that is central to this
project.

**Decision required:** approve BRCA1 (or choose TP53 / another) before Phase 1 begins.

## 5. Files created in Phase 0

- `README.md`, DRAFT `protocol.md`, `CHANGELOG.md`, `.gitignore`
- `config/config.yaml` (parameters), `data/raw/README.md` (download log)
- `src/` module stubs (8 modules + `__init__.py`)
- this report

**No variants have been downloaded or analyzed yet.**
