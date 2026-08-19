"""ClinVar variant_summary.txt retrieval/parsing for a single gene (Phase 2).

Streams the raw .gz (never fully in memory, never modifies the raw file), extracts all
rows whose GeneSymbol equals the target gene, and writes a verbatim gene-level subset plus
a summary JSON. No biological filtering is applied at this stage.

Column notes (verified against build 260818-0035.1, 2026-08-19):
- one row per (AlleleID x Assembly); each VariationID appears as a GRCh37 row AND a GRCh38 row.
- Chromosome/Start/Stop and PositionVCF are on the row's Assembly.
- ReferenceAllele/AlternateAllele are deprecated ('na'); use ReferenceAlleleVCF/
  AlternateAlleleVCF (GRCh38-consistent) for actual alleles.
- Name carries HGVS ("NM_007294.4(BRCA1):c.190T>G (p.Cys64Gly)").
- PhenotypeList is pipe-delimited (multi-condition).
"""
import gzip
import json
import os

import yaml

DEFAULT_CONFIG = "config/config.yaml"
GENE_COL = 4          # GeneSymbol
TYPE_COL = 1
NAME_COL = 2
CLNSIG_COL = 6
ASSEMBLY_COL = 16
CHROM_COL = 18
REVIEW_COL = 24
VID_COL = 30
N_COLUMNS = 43


def load_config(path=DEFAULT_CONFIG):
    with open(path) as fh:
        return yaml.safe_load(fh)


def parse_gene_subset(gz_path, gene_symbol, out_tsv, out_summary):
    """Stream gz_path, keep rows for gene_symbol, write out_tsv + out_summary.

    Returns the stats dict.
    """
    header = None
    total_rows = 0
    gene_rows = 0
    vids = set()
    allele_ids = set()
    vid_sig = {}
    vid_review = {}
    vid_name = {}
    assembly_rows = {}
    chrom_grch38 = {}

    with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as fh, \
            open(out_tsv, "w", encoding="utf-8") as out:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                out.write("\t".join(fields) + "\n")
                continue
            total_rows += 1
            if len(fields) < N_COLUMNS:
                continue  # malformed short row — counted as skipped (see stats note)
            if fields[GENE_COL] == gene_symbol:
                gene_rows += 1
                out.write("\t".join(fields) + "\n")

                vid = fields[VID_COL]
                aid = fields[0]
                assembly = fields[ASSEMBLY_COL]
                if vid:
                    vids.add(vid)
                if aid not in ("", "-", "na"):
                    allele_ids.add(aid)
                assembly_rows[assembly] = assembly_rows.get(assembly, 0) + 1
                if vid and vid not in vid_sig:
                    vid_sig[vid] = fields[CLNSIG_COL]
                    vid_review[vid] = fields[REVIEW_COL]
                    vid_name[vid] = fields[NAME_COL]
                if assembly == "GRCh38" and vid:
                    chrom_grch38[vid] = fields[CHROM_COL]

    stats = {
        "gene": gene_symbol,
        "n_columns": len(header) if header else None,
        "total_rows_in_file": total_rows,
        "gene_rows": gene_rows,
        "unique_variation_ids": len(vids),
        "unique_allele_ids": len(allele_ids),
        "assembly_rows": assembly_rows,
        "chromosome_of_unique_variants_grch38": _counter(chrom_grch38.values()),
        "clinical_significance": _counter(vid_sig.values()),
        "review_status": _counter(vid_review.values()),
        "missing_hgvs_c": sum(1 for n in vid_name.values() if ":c." not in n),
        "missing_protein": sum(1 for n in vid_name.values() if "p." not in n),
    }
    os.makedirs(os.path.dirname(out_summary), exist_ok=True)
    with open(out_summary, "w") as fh:
        json.dump(stats, fh, indent=2)
    return stats


def _counter(iterable):
    out = {}
    for x in iterable:
        out[x] = out.get(x, 0) + 1
    return out


def main():
    cfg = load_config()
    gene = cfg["study"]["gene"]
    raw_dir = cfg["output_dirs"]["raw"]
    int_dir = cfg["output_dirs"]["intermediate"]
    gz_path = os.path.join(raw_dir, "clinvar", cfg["clinvar"]["file"])
    out_tsv = os.path.join(int_dir, f"clinvar_{gene.lower()}_raw.tsv")
    out_summary = os.path.join(int_dir, f"clinvar_{gene.lower()}_summary.json")
    if os.path.exists(out_tsv):
        print(f"Output exists ({out_tsv}); skipping parse (idempotent).")
        return
    stats = parse_gene_subset(gz_path, gene, out_tsv, out_summary)
    print(json.dumps(stats, indent=2))
    print(f"\nWrote {out_tsv}")
    print(f"Wrote {out_summary}")


if __name__ == "__main__":
    main()
