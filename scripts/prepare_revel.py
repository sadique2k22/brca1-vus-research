"""Phase 4B — download the REVEL GRCh38 file and extract the BRCA1 (chr17) region.

REVEL is distributed as a single ~636 MB zip (ZIP64, one entry 'revel_with_transcript_ids',
uncompressed >4 GB). This script downloads it (skippable if cached), streams it line-by-line,
detects columns from the header, and writes only the chr17 BRCA1 region to a small TSV.
"""
import os
import sys
import zipfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import requests

REVEL_URL = "https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip"
ZIP_PATH = "data/raw/revel/revel-v1.3_all_chromosomes.zip"
OUT_TSV = "data/intermediate/revel_brca1.tsv"
CHROM = "17"
REGION_LO = 43_000_000
REGION_HI = 43_200_000


def download(url, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    print(f"Downloading {url} ...")
    with requests.get(url, stream=True, timeout=600) as r:
        r.raise_for_status()
        total = 0
        with open(path, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
                total += len(chunk)
    print(f"  downloaded {total} bytes")


def extract(zip_path, out_tsv, chrom=CHROM, lo=REGION_LO, hi=REGION_HI):
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
        # pick the single data file (skip directories)
        entry = next((n for n in names if not n.endswith("/")), None)
        print(f"Zip entries: {names[:5]}; using '{entry}'")
        with zf.open(entry) as raw:
            header = raw.readline().decode("utf-8", "replace").rstrip("\n").split("\t")
            cols = {c: i for i, c in enumerate(header)}
            pos_col = next((c for c in ("grch38_pos", "pos", "hg19_pos") if c in cols), None)
            if pos_col is None:
                raise SystemExit(f"REVEL position column not found in header: {header}")
            chrom_col = cols.get("chr", cols.get("chrom", 0))
            os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
            n_kept = 0
            with open(out_tsv, "w") as out:
                out.write("\t".join(header) + "\n")
                for line in raw:
                    parts = line.decode("utf-8", "replace").rstrip("\n").split("\t")
                    if len(parts) <= pos_col:
                        continue
                    if parts[chrom_col] != chrom:
                        continue
                    pos = parts[pos_col]
                    if pos.isdigit() and lo <= int(pos) <= hi:
                        out.write("\t".join(parts) + "\n")
                        n_kept += 1
            print(f"Header: {header}")
            print(f"Extracted {n_kept} REVEL rows for chr{chrom}:{lo}-{hi} -> {out_tsv}")


if __name__ == "__main__":
    if not os.path.exists(ZIP_PATH):
        download(REVEL_URL, ZIP_PATH)
    extract(ZIP_PATH, OUT_TSV)
