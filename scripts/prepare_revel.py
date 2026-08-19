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
        entry = next((n for n in names if not n.endswith("/")), None)
        print(f"Zip entries: {names[:5]}; using '{entry}'")
        with zf.open(entry) as raw:
            header_line = raw.readline().decode("utf-8", "replace").rstrip("\n")
            delim = "," if "," in header_line else "\t"
            header = [h.strip() for h in header_line.split(delim)]
            cols = {c: i for i, c in enumerate(header)}
            print(f"REVEL header ({delim}-separated): {header}")
            if "grch38_pos" not in cols:
                raise SystemExit(f"grch38_pos column not found in header: {header}")
            pos_col = cols["grch38_pos"]
            chrom_col = cols.get("chr", 0)
            ref_col = cols["ref"]
            alt_col = cols["alt"]
            score_col = cols["REVEL"]
            os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
            n_kept = 0
            with open(out_tsv, "w") as out:
                out.write("chr\tgrch38_pos\tref\talt\tREVEL\n")
                for line in raw:
                    parts = line.decode("utf-8", "replace").rstrip("\n").split(delim)
                    if len(parts) <= max(chrom_col, pos_col, ref_col, alt_col, score_col):
                        continue
                    if parts[chrom_col].lstrip("chr") != chrom:
                        continue
                    pos = parts[pos_col]
                    if pos.isdigit() and lo <= int(pos) <= hi:
                        out.write(f"{chrom}\t{pos}\t{parts[ref_col]}\t{parts[alt_col]}\t{parts[score_col]}\n")
                        n_kept += 1
            print(f"Extracted {n_kept} REVEL rows for chr{chrom}:{lo}-{hi} -> {out_tsv}")


if __name__ == "__main__":
    if not os.path.exists(ZIP_PATH):
        download(REVEL_URL, ZIP_PATH)
    extract(ZIP_PATH, OUT_TSV)
