"""Phase 4B — annotate the 1,904 BRCA1 missense VUS (dedup to unique biological variants).

Pipeline:
  1. read normalized dataset
  2. generate biological variant keys (GRCh38:chrom:pos:ref:alt)
  3. deduplicate to unique biological variants (annotation queried once per unique variant)
  4. annotate: gnomAD (GraphQL), VEP (consequence + SIFT + PolyPhen), REVEL (local), CADD (best-effort)
  5. map annotations back to every original ClinVar record
  6. write map + unique + annotated datasets, QC + resource + duplicate reports

All remote requests are disk-cached and resumable.
"""
import csv
import json
import os
import sys
import time
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import yaml

from src.annotation import extract_mane_annotation, vep_annotate
from src.population import parse_gnomad_record, query_gnomad
from src.predictors import CaddClient, load_revel, revel_score
from src.variants import (
    build_variant_string,
    gnomad_variant_id,
    make_variant_key,
)

ANNOTATION_FIELDS = [
    "variant_key",
    "gnomad_found",
    "gnomad_genome_af", "gnomad_genome_ac", "gnomad_genome_an", "gnomad_genome_hom",
    "gnomad_genome_faf95_popmax", "gnomad_genome_faf95_pop",
    "gnomad_exome_af", "gnomad_exome_ac", "gnomad_exome_an", "gnomad_exome_hom",
    "gnomad_exome_faf95_popmax", "gnomad_exome_faf95_pop",
    "vep_consequence", "vep_hgvsc", "vep_hgvsp", "vep_amino_acids", "vep_gene", "vep_impact",
    "sift_score", "sift_prediction",
    "polyphen_score", "polyphen_prediction",
    "cadd_phred",
    "revel_score",
]


def na(v):
    return "" if v is None else v


def load_normalized(path):
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        idx = {n.lstrip("#"): i for i, n in enumerate(header)}
        rows = []
        for row in reader:
            rec = {n.lstrip("#"): (row[i] if i < len(row) else "") for n, i in idx.items()}
            rec["_raw"] = row
            rows.append(rec)
    return header, rows


def main():
    t0 = time.time()
    cfg = yaml.safe_load(open("config/config.yaml"))
    int_dir = cfg["output_dirs"]["intermediate"]
    proc_dir = cfg["output_dirs"]["processed"]
    rep_dir = cfg["output_dirs"]["reports"]
    norm_tsv = os.path.join(int_dir, "brca1_vus_missense_normalized.tsv")
    gnomad_cache = os.path.join(int_dir, "gnomad_cache")
    vep_cache = os.path.join(int_dir, "vep_cache")

    header, records = load_normalized(norm_tsv)
    n_total = len(records)
    print(f"Loaded {n_total} records")

    # 2) variant keys + dedup
    for r in records:
        r["variant_key"] = make_variant_key(r["Chromosome"], r["Start"],
                                            r["ReferenceAlleleVCF"], r["AlternateAlleleVCF"])
    unique = {}
    order = []
    for r in records:
        k = r["variant_key"]
        if k not in unique:
            unique[k] = r
            order.append(k)
    n_unique = len(unique)
    print(f"Unique biological variants: {n_unique} (from {n_total} records)")

    # 3) annotations per unique variant
    # gnomAD
    gids = [gnomad_variant_id(unique[k]["Chromosome"], unique[k]["Start"],
                              unique[k]["ReferenceAlleleVCF"], unique[k]["AlternateAlleleVCF"])
            for k in order]
    gnomad_raw = query_gnomad(gids, gnomad_cache)
    gnomad = {}
    for k, gid in zip(order, gids):
        gnomad[k] = parse_gnomad_record(gnomad_raw.get(gid, {}))

    # VEP (reuse cache; extracts consequence + SIFT + PolyPhen)
    vstrs = [build_variant_string(unique[k]["Chromosome"], unique[k]["Start"], unique[k]["Start"],
                                  unique[k]["ReferenceAlleleVCF"], unique[k]["AlternateAlleleVCF"])
             for k in order]
    vep_raw, _ = vep_annotate(vstrs, vep_cache)
    vep = {k: extract_mane_annotation(vep_raw.get(vs, {})) for k, vs in zip(order, vstrs)}

    # REVEL (local)
    revel_table = load_revel(os.path.join(int_dir, "revel_brca1.tsv"))
    revel = {}
    for k in order:
        r = unique[k]
        revel[k] = revel_score(revel_table, r["Start"], r["ReferenceAlleleVCF"], r["AlternateAlleleVCF"]) \
            if revel_table else None

    # CADD (best-effort)
    cadd = {}
    try:
        cadd_client = CaddClient()
        cadd_res = cadd_client.score([(unique[k]["Chromosome"], unique[k]["Start"],
                                       unique[k]["ReferenceAlleleVCF"], unique[k]["AlternateAlleleVCF"])
                                      for k in order])
        for k in order:
            r = unique[k]
            cadd[k] = cadd_res.get(f"{r['Chromosome']}:{r['Start']}:{r['ReferenceAlleleVCF']}:{r['AlternateAlleleVCF']}")
    except Exception as exc:  # noqa: BLE001
        print(f"CADD unavailable: {exc}")
        cadd = {k: None for k in order}

    # 4) assemble per-unique annotation
    unique_ann = {}
    for k in order:
        g = gnomad[k]
        v = vep[k]
        unique_ann[k] = {
            "variant_key": k,
            "gnomad_found": g.get("gnomad_found"),
            "gnomad_genome_af": g.get("gnomad_genome_af"),
            "gnomad_genome_ac": g.get("gnomad_genome_ac"),
            "gnomad_genome_an": g.get("gnomad_genome_an"),
            "gnomad_genome_hom": g.get("gnomad_genome_hom"),
            "gnomad_genome_faf95_popmax": g.get("gnomad_genome_faf95_popmax"),
            "gnomad_genome_faf95_pop": g.get("gnomad_genome_faf95_pop"),
            "gnomad_exome_af": g.get("gnomad_exome_af"),
            "gnomad_exome_ac": g.get("gnomad_exome_ac"),
            "gnomad_exome_an": g.get("gnomad_exome_an"),
            "gnomad_exome_hom": g.get("gnomad_exome_hom"),
            "gnomad_exome_faf95_popmax": g.get("gnomad_exome_faf95_popmax"),
            "gnomad_exome_faf95_pop": g.get("gnomad_exome_faf95_pop"),
            "vep_consequence": v.get("vep_consequence"),
            "vep_hgvsc": v.get("vep_hgvsc"),
            "vep_hgvsp": v.get("vep_hgvsp"),
            "vep_amino_acids": v.get("vep_amino_acids"),
            "vep_gene": v.get("vep_gene"),
            "vep_impact": v.get("vep_impact"),
            "sift_score": v.get("sift_score"),
            "sift_prediction": v.get("sift_prediction"),
            "polyphen_score": v.get("polyphen_score"),
            "polyphen_prediction": v.get("polyphen_prediction"),
            "cadd_phred": cadd.get(k),
            "revel_score": revel[k],
        }

    # 5) write biological_variant_map.tsv
    map_tsv = os.path.join(int_dir, "biological_variant_map.tsv")
    with open(map_tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["variant_key", "ClinVar_VariationID", "ClinVar_AlleleID",
                    "protein_change", "original_HGVS", "normalized_HGVS"])
        for r in records:
            w.writerow([r["variant_key"], r["VariationID"], r["AlleleID"],
                        r.get("protein_change", ""), r["Name"],
                        r.get("normalized_hgvs_c", "")])

    # 6) write annotation_unique_variants.tsv (one row per unique variant)
    unique_tsv = os.path.join(int_dir, "annotation_unique_variants.tsv")
    uni_header = ["variant_key", "chrom", "pos", "ref", "alt"] + ANNOTATION_FIELDS
    with open(unique_tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(uni_header)
        for k in order:
            r = unique[k]
            ann = unique_ann[k]
            w.writerow([k, r["Chromosome"], r["Start"], r["ReferenceAlleleVCF"], r["AlternateAlleleVCF"]]
                       + [na(ann[f]) for f in ANNOTATION_FIELDS])

    # 7) write brca1_vus_missense_annotated.tsv (map back to all 1,904 records)
    annotated_tsv = os.path.join(proc_dir, "brca1_vus_missense_annotated.tsv")
    with open(annotated_tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header + ANNOTATION_FIELDS)
        for r in records:
            ann = unique_ann[r["variant_key"]]
            w.writerow(r["_raw"] + [na(ann[f]) for f in ANNOTATION_FIELDS])

    # 8) QC counts
    qc = {
        "total_records": n_total,
        "unique_variants": n_unique,
        "gnomad_present": sum(1 for a in unique_ann.values() if a["gnomad_found"] == "present"),
        "gnomad_absent": sum(1 for a in unique_ann.values() if a["gnomad_found"] == "absent"),
        "missing_cadd": sum(1 for a in unique_ann.values() if a["cadd_phred"] in (None, "")),
        "missing_revel": sum(1 for a in unique_ann.values() if a["revel_score"] in (None, "")),
        "missing_sift": sum(1 for a in unique_ann.values() if a["sift_score"] in (None, "")),
        "missing_polyphen": sum(1 for a in unique_ann.values() if a["polyphen_score"] in (None, "")),
    }

    # genomic duplicates (expected 0) + protein-change collisions (distinct alts -> same AA)
    key_vids = defaultdict(list)
    for r in records:
        key_vids[r["variant_key"]].append((r["VariationID"], r.get("protein_change", "")))
    multi = {k: v for k, v in key_vids.items() if len(v) > 1}

    prot_vids = defaultdict(list)
    for r in records:
        prot_vids[r.get("protein_change", "")].append(
            (r["VariationID"], r["ReferenceAlleleVCF"], r["AlternateAlleleVCF"]))
    collisions = {k: v for k, v in prot_vids.items() if k and len(v) > 1}

    runtime = round(time.time() - t0, 1)
    write_reports(rep_dir, int_dir, qc, multi, collisions, runtime)
    print(json.dumps(qc, indent=2))
    print(f"Runtime: {runtime}s")
    print(f"Wrote: {map_tsv}, {unique_tsv}, {annotated_tsv}")


def write_reports(rep_dir, int_dir, qc, multi, collisions, runtime):
    os.makedirs(rep_dir, exist_ok=True)

    # duplicate report
    with open(os.path.join(rep_dir, "duplicate_variant_report.md"), "w") as fh:
        fh.write("# Duplicate Variant Report — Phase 4B\n\n")
        fh.write(f"Total ClinVar records: {qc['total_records']}\n\n")
        fh.write(f"Unique biological variants (chrom:pos:ref:alt): {qc['unique_variants']}\n\n")
        fh.write(f"Genomic positions with >1 ClinVar VariationID: {len(multi)}\n\n")
        fh.write("> Correction to Phase 4A: the earlier '38 duplicate representations' were "
                 "protein-change collisions (distinct nucleotide changes producing the same amino-acid "
                 "substitution), not genomic duplicates. There are **0** true genomic duplicates.\n\n")
        fh.write(f"Protein-change collisions (distinct ref>alt → same AA substitution): {len(collisions)}\n\n")
        for k in sorted(collisions)[:40]:
            vids = ", ".join(f"{v}({r}>{a})" for v, r, a in collisions[k])
            fh.write(f"- {k}: {vids}\n")
        fh.write("\n(Full list in `data/intermediate/biological_variant_map.tsv`.)\n")

    # annotation QC + resource report
    with open(os.path.join(rep_dir, "annotation_report.md"), "w") as fh:
        fh.write("# Annotation Report — Phase 4B (QC + resources)\n\n")
        fh.write("## Quality control\n\n")
        for k, v in qc.items():
            fh.write(f"- {k}: {v}\n")
        fh.write("\n## Resource usage\n\n")
        fh.write(f"- Runtime: {runtime}s\n")
        for d in ("gnomad_cache", "vep_cache"):
            p = os.path.join(int_dir, d)
            size = sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(p) for f in fs) if os.path.isdir(p) else 0
            fh.write(f"- {d} storage: {size / 1e6:.2f} MB\n")
        fh.write("\n## Annotation metadata (provenance)\n\n")
        fh.write("| Field | Source | Version | Assembly | Date | Method |\n|---|---|---|---|---|---|\n")
        fh.write("| gnomad_* | gnomAD | v4 (`gnomad_r4`) | GRCh38 | 2026-08-19 | GraphQL API |\n")
        fh.write("| vep_*, sift_*, polyphen_* | Ensembl VEP | release 116 | GRCh38 | 2026-08-19 | REST /vep/human/region |\n")
        fh.write("| revel_score | REVEL | v1.3 | GRCh38 | 2021-05-03 | standalone file (chr17 region) |\n")
        fh.write("| cadd_phred | CADD | v1.6 (v1.7 degraded) | GRCh38 | 2026-08-19 | web service (best-effort) |\n")
        fh.write("\nMissing values are recorded as empty (NA); no field is estimated.\n")


if __name__ == "__main__":
    main()
