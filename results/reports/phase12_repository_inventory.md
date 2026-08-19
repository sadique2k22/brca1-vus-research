# Phase 12 Repository Inventory

Date: 2026-08-19. Total tracked files: 152 (before Phase 12 additions).

## Categorization

| Category | Path | Count | Notes |
|---|---|---|---|
| SOURCE | `src/*.py` | 11 | pipeline modules |
| SCRIPT | `scripts/*.py` | 14 | phase orchestrators + validators |
| CONFIG | `config/`, `configs/` | 3 | pipeline + analysis + cohort configs |
| TEST | `tests/*.py` | 9 | unit tests (62 cases) |
| DATA-RAW | `data/raw/` | 3 | README (download log), ClinVar metadata.json, .md5 |
| DATA-PROCESSED | `data/processed/` | 4 | annotated (frozen) + with-functional + .gitkeep |
| DATA-INTERMEDIATE | `data/intermediate/` | 6 | normalized, map, unique, summary, unresolved, .gitkeep |
| RESULT-TABLE | `results/tables/` | 18 | generated tables |
| REPORT | `results/reports/` | 34 | phase + audit reports |
| FIGURE | `results/figures/` | 37 | fig1–fig16 (SVG + 300-dpi PNG) |
| MANUSCRIPT | `manuscript/` | 6+ | manuscript.md + versioned + mapping (+ Phase 12 additions) |
| DOCUMENTATION | root `.md`, `requirements.txt` | 4 | README, protocol, CHANGELOG, requirements |
| CI | `.github/workflows/` | 1 | pipeline.yml |

## Findings

- **Stale / obsolete files:** none identified. `protocol.md` is the Phase 1 frozen protocol (intentional historical artifact).
- **Duplicate files:** none.
- **Temporary files:** none tracked (all `*_cache/` and raw `.gz`/`.zip` are gitignored).
- **Accidentally tracked caches:** none.
- **Missing files (referenced but absent):** none — all scripts/reports referenced by the README and manifests exist.
- **Undocumented files:** none — all tracked files are under documented directories.

## Git-ignored (correctly untracked)

`data/raw/*.gz`, `*.zip`, `data/intermediate/*` (regenerable), all `*_cache/`, `__pycache__/`,
`.venv/`, secrets. The 442 MB ClinVar gzip and 636 MB REVEL zip are **not** committed (provenance
in `metadata.json`); they are re-fetched and checksum-verified in CI.
