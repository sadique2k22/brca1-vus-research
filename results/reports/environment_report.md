# Environment Report — Resource-Constrained Mobile Linux (Android / Termux / PRoot)

Date: 2026-08-19

## 1. Context

The project runs entirely on an Android phone:
`Android → Termux → proot-distro (Ubuntu) → OpenCode`.
The Ubuntu rootfs lives under
`/data/data/com.termux/files/usr/var/lib/proot-distro/containers/ubuntu/rootfs`
(confirmed via `/proc/self/root`), i.e. on Android internal flash storage.

**PRoot caveats (important for interpretation):** some `df`/`/proc` values are misleading
under PRoot. The *mount table* (`mount`) is more authoritative than `df` for `/tmp`.
`/proc` and `/sys` reflect the host; no cgroup memory limit was readable.

## 2. Hardware / system

| Item | Value | Notes |
|---|---|---|
| CPU architecture | aarch64 (ARM64) | |
| CPU cores | 8 physical: 1× Cortex-X4 + 7× Cortex-A720 | big.LITTLE |
| CPUs usable to us | **6** (`nproc`=6; cpuinfo lists 8) | affinity-limited |
| CPU max freq | ~3.01 GHz | |
| Throttling | **already active** (~46–55% scaling MHz) | overheat risk confirmed |
| Kernel | 6.17.0-PRoot-Distro (proot@termux) | |
| Total RAM | 11.5 GB (MemTotal 11,510,480 kB) | |
| **Available RAM** | **~2.3 GB** (MemAvailable 2,374,984 kB) | most RAM held by Android/other apps |
| Free RAM (instant) | ~84 MB | tight |
| Swap | 12 GB total / 7.3 GB free | Android zram |
| Storage total | 479 GB (f2fs) | |
| **Storage free** | **~172 GB** | ample |
| Filesystem (project) | **f2fs** (block size 4096) | Flash-Friendly FS, on `/dev/block/dm-56` |

## 3. Software / tools

| Tool | Status | Detail |
|---|---|---|
| Python | ✅ 3.13.7 | in `/root/frida-venv` |
| pip | ✅ 26.2 | |
| git | ✅ 2.51.0 | |
| curl / wget | ✅ 8.14.1 / 1.25.0 | curl has zstd + brotli |
| Java | ✅ OpenJDK 17.0.11 (Temurin) | via sdkman (`/root/.sdkman/...`) |
| gzip / bzip2 / xz | ✅ | single-threaded (no `pigz`) |
| zstd / 7z / zip / unzip / tar | ✅ | |
| zcat | ✅ | |
| pigz / lz4 / lzip | ❌ | |
| bcftools / tabix / samtools / bgzip | ❌ | not installed |
| VEP | ❌ | not installed |
| Docker | ❌ | not installed (expected under PRoot) |
| conda / micromamba | ❌ | not installed |
| Internet (NCBI, ClinVar FTP, gnomAD API, Ensembl REST) | ✅ | all reachable |

## 4. Current project size

~527 KB total (398 KB is `.git`). Negligible.

## 5. Resource-budget estimate vs. proposed workflow

### 5.1 Will we exceed 2 GB RAM?
- **Only risk:** naively loading the full ClinVar `variant_summary.txt`
  (~1.5–2 GB uncompressed) into pandas → **4–6 GB peak**, i.e. **would exceed 2 GB and
  risk OOM** (only ~2.3 GB free).
- **Mitigation:** stream-filter first (`zcat … | grep/awk` → small gene-specific file), then
  load. Peak memory then stays **< 500 MB**.
- gnomAD / VEP / predictors are API-per-variant (cached), tiny per request.
- **Verdict: SAFE (< 500 MB) *only if* we stream-filter ClinVar. Never load the full table
  into memory.**

### 5.2 Will we exceed 5 GB temporary storage?
- Full decompression of `variant_summary.txt.gz` to disk ≈ **~2 GB single file** (would fit
  but wasteful; streaming avoids it).
- gnomAD cache ≈ 10–50 MB; VEP cache ≈ 25–125 MB.
- **Verdict: SAFE (< 1 GB)** with streaming. **Do NOT use `/tmp` for large temp files** —
  `/tmp` is a RAM-backed tmpfs (mount size ≈ 5.6 GB) sharing the already-tight RAM; use
  `data/intermediate/` on f2fs instead.

### 5.3 Will we exceed the 50 GB total project storage cap?
- **Cap: 50 GB** (user-specified hard ceiling for this project).
- Raw: `variant_summary.txt.gz` ≈ 0.2 GB (+ optional gene VCF ≈ 0.05 GB).
- Caches + processed tables/figures/reports ≈ < 0.5 GB.
- **Verdict: SAFE.** Baseline workflow ≈ < 1 GB, i.e. ~2% of the cap. Headroom is generous.

### 5.4 Dataset/tool storage limits (vs the 50 GB cap)
| Dataset/tool | Est. size | Status under 50 GB cap |
|---|---|---|
| Whole-genome gnomAD v4 VCF | ~100+ GB | ❌ prohibited (exceeds cap) |
| Full CADD genome scores | ~300 GB | ❌ prohibited (exceeds cap) |
| Full Ensembl VEP cache (GRCh38) | ~20–40 GB | ⚠️ fits but discouraged; needs explicit approval |
| dbNSFP (REVEL+CADD+SIFT+PolyPhen bundle) | ~15–30 GB | ⚠️ fits but discouraged; needs explicit approval |
| Standalone REVEL GRCh38 file | ~GB-scale (verify) | ✅ feasible within budget (verify size first) |

## 6. Workflow adjustments for the mobile environment

1. **API-first annotation.** Use Ensembl VEP REST (SIFT, PolyPhen-2) + CADD web service +
   gnomAD GraphQL API. No local VEP cache, no dbNSFP.
2. **Streaming/chunked I/O.** `zcat | grep` to build the gene-specific ClinVar subset before
   any pandas load; process in chunks; never hold the full table in memory.
3. **Caching.** Cache every API response to disk (per-variant JSON) so retries don't
   re-fetch; respect rate limits.
4. **No duplicate large copies.** One raw copy, derive new files, delete intermediates when
   superseded.
5. **Concurrency.** 1–4 workers max (throttling already active; phone can overheat). API
   calls serial or lightly concurrent (2–4), rate-limit-aware.
6. **No Docker.** PRoot + no Docker; use pure-Python/streaming tools.
7. **No sudo/systemd assumption.** The environment is already root inside PRoot, but no
   system services.
8. **Temp on f2fs, not `/tmp`.** `/tmp` is RAM-backed.

## 7. Predictor feasibility (Phase 1 decision point — not decided here)

| Predictor | Mobile-feasible route | Notes |
|---|---|---|
| SIFT | Ensembl VEP REST | ✅ lightweight |
| PolyPhen-2 | Ensembl VEP REST | ✅ lightweight |
| CADD | CADD web service (chunked) | ✅ feasible for ~5k variants |
| REVEL | standalone REVEL file (~GB, within 50 GB budget) or dbNSFP (discouraged) | ⚠️ needs decision + size check |
| BayesDel | dbNSFP (discouraged) | ⚠️ needs decision |

Recommendation to confirm in Phase 1: run SIFT + PolyPhen-2 (VEP REST) + CADD (web) as the
core predictor set; treat REVEL as *optional* and only add it with explicit approval + size
verification of the per-chromosome REVEL file.

## 8. Remaining uncertainty

- Exact current sizes of `variant_summary.txt.gz`, standalone REVEL, and CADD files (verify
  via `HEAD` before any download).
- Whether the gnomAD public GraphQL API will comfortably serve ~5k variant queries within
  rate limits (test with a small batch first).
- No cgroup memory cap readable — assume the ~2.3 GB available figure is the real ceiling.

## 9. Files created/changed

- This report: `results/reports/environment_report.md`.
- `README.md` updated with a "Computational Environment" section.
- `CHANGELOG.md` updated.

**No software installed. No large datasets downloaded. STOP — awaiting approval.**
