"""Computational pathogenicity predictors: REVEL (local file) and CADD (best-effort web)."""
import os

REVEL_DEFAULT_PATH = "data/intermediate/revel_brca1.tsv"


def load_revel(path=REVEL_DEFAULT_PATH):
    """Load the extracted BRCA1-region REVEL file into {pos: {f'{ref}>{alt}': score}}.

    Expects a TSV with a header. Detects the position column (grch38_pos/pos/hg19_pos),
    ref/alt columns, and the REVEL score column.
    """
    if not os.path.exists(path):
        return None
    table = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cols = {c: i for i, c in enumerate(header)}
        pos_col = next((c for c in ("grch38_pos", "pos", "hg19_pos", "grch37_pos") if c in cols), None)
        ref_col = next((c for c in ("ref", "REF") if c in cols), None)
        alt_col = next((c for c in ("alt", "ALT") if c in cols), None)
        score_col = next((c for c in ("REVEL", "revel", "score") if c in cols), None)
        if None in (pos_col, ref_col, alt_col, score_col):
            raise ValueError(f"REVEL columns not detected: {header}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(pos_col, ref_col, alt_col, score_col):
                continue
            pos, ref, alt, score = f[pos_col], f[ref_col], f[alt_col], f[score_col]
            try:
                table.setdefault(pos, {})[f"{ref}>{alt}"] = float(score)
            except ValueError:
                continue
    return table


def revel_score(table, pos, ref, alt):
    """Look up a REVEL score; return None if absent."""
    if not table:
        return None
    return table.get(str(pos), {}).get(f"{ref}>{alt}")


class CaddClient:
    """Best-effort CADD web-service client.

    NOTE: CADD v1.7 online scoring is degraded (hardware failure, long queue) and there is
    no clean bulk REST API. This client posts a VCF to the /upload form (v1.6, which is not
    delayed) and attempts to parse a result. If unusable, scores are left missing and the
    limitation is documented.
    """

    URL = "https://cadd.gs.washington.edu/upload"
    VERSION = "v1.6"

    def score(self, variants):
        """variants: list of (chrom, pos, ref, alt). Returns dict keyed by 'chrom:pos:ref:alt'."""
        import requests

        vcf = ["##fileformat=VCFv4.2",
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
        keys = []
        for chrom, pos, ref, alt in variants:
            vcf.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.")
            keys.append(f"{chrom}:{pos}:{ref}:{alt}")
        body = "\n".join(vcf) + "\n"
        try:
            r = requests.post(self.URL,
                              data={"version": self.VERSION, "email": ""},
                              files={"file": ("variants.vcf", body)},
                              timeout=300)
            r.raise_for_status()
            return self._parse(r.text, keys)
        except Exception:
            return {k: None for k in keys}

    @staticmethod
    def _parse(text, keys):
        out = {k: None for k in keys}
        for line in text.splitlines():
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 6 and parts[0].isdigit():
                key = f"{parts[0]}:{parts[1]}:{parts[3]}:{parts[4]}"
                if key in out:
                    try:
                        out[key] = float(parts[5])
                    except (ValueError, TypeError):
                        out[key] = None
        return out
