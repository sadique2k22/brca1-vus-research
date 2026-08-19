"""Variant representation, classification, and VUS/missense filtering.

Pure functions so the filtering logic is unit-testable. Classifies ClinVar protein-change
notation (3-letter amino-acid codes, as used in variant_summary.txt "Name" field) and
implements the FINAL inclusion criteria from protocol.md v1.0 section 6/7.
"""
import re

AA3 = ("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
       "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")
_AA3 = "(" + "|".join(AA3) + ")"


def extract_protein_change(name):
    """Return the p. token from a ClinVar Name (e.g. 'Cys64Gly') or None."""
    if not name:
        return None
    m = re.search(r'p\.([^)\s]+)', name)
    return m.group(1) if m else None


def classify_protein_change(name):
    """Classify a ClinVar protein change into a consequence class.

    Returns one of: missense, synonymous, nonsense, frameshift, inframe_indel,
    unknown, no_protein, other.
    """
    ch = extract_protein_change(name)
    if ch is None:
        return "no_protein"
    if "=" in ch:
        return "synonymous"
    if "?" in ch:
        return "unknown"
    if "fs" in ch:
        return "frameshift"
    if "del" in ch or "ins" in ch or "dup" in ch:
        return "inframe_indel"
    if ch.startswith("Ter") or "*" in ch:
        return "nonsense"
    m = re.match(r'^%s([0-9]+)%s$' % (_AA3, _AA3), ch)
    if m:
        return "synonymous" if m.group(1) == m.group(3) else "missense"
    return "other"


def is_vus(clinical_significance):
    """True iff aggregate germline clinical significance is exactly VUS."""
    return clinical_significance == "Uncertain significance"


def is_missense(vtype, name):
    """Protocol missense: single-nucleotide substitution that changes one amino acid."""
    return (vtype == "single nucleotide variant"
            and classify_protein_change(name) == "missense")


def filter_vus_missense(records):
    """Apply FINAL inclusion criteria to a list of record dicts.

    Each record dict must contain at least: VariationID, Assembly,
    ClinicalSignificance, Type, Name.

    Returns (kept, steps) where `steps` is a list of
    (step_name, before, removed, remaining, reason).
    """
    steps = []

    # Step 1: deduplicate by VariationID, preferring the GRCh38 row.
    n0 = len(records)
    by_vid = {}
    for r in records:
        key = r.get("VariationID") or r.get("AlleleID") or id(r)
        cur = by_vid.get(key)
        if cur is None or (r.get("Assembly") == "GRCh38" and cur.get("Assembly") != "GRCh38"):
            by_vid[key] = r
    deduped = list(by_vid.values())
    steps.append(("deduplicate by VariationID (prefer GRCh38 row)",
                  n0, n0 - len(deduped), len(deduped),
                  "remove duplicate per-assembly rows"))

    # Step 2: require GRCh38 coordinates.
    grch38 = [r for r in deduped if r.get("Assembly") == "GRCh38"]
    steps.append(("require GRCh38 coordinates",
                  len(deduped), len(deduped) - len(grch38), len(grch38),
                  "no GRCh38 row (GRCh37-only or unmapped)"))

    # Step 3: require VUS (exact aggregate 'Uncertain significance').
    vus = [r for r in grch38 if is_vus(r.get("ClinicalSignificance"))]
    steps.append(("require ClinVar VUS (Uncertain significance)",
                  len(grch38), len(grch38) - len(vus), len(vus),
                  "aggregate significance not 'Uncertain significance'"))

    # Step 4: require missense.
    miss = [r for r in vus if is_missense(r.get("Type"), r.get("Name"))]
    steps.append(("require missense (SNP + amino-acid substitution)",
                  len(vus), len(vus) - len(miss), len(miss),
                  "not a single-nucleotide amino-acid substitution"))

    return miss, steps


# ---- Normalization helpers (Phase 4A) ----

VALID_BASES = {"A", "C", "G", "T"}

# BRCA1 MANE Select / ENIGMA-ClinGen VCEP transcript (protocol v1.0 section 4)
NORMALIZED_TRANSCRIPT = "NM_007294.4"
MANE_ENST = "ENST00000357654"


def validate_alleles(ref, alt):
    """Return (ok, reason). ok iff ref/alt are single valid bases and differ."""
    if not ref or not alt:
        return False, "empty allele"
    if ref == alt:
        return False, "ref == alt"
    if ref not in VALID_BASES or alt not in VALID_BASES:
        return False, "non-single-base or invalid allele"
    return True, ""


def build_variant_string(chrom, start, stop, ref, alt):
    """VEP REST region-input string, e.g. '17 43082542 43082542 G/C 1'."""
    return f"{chrom} {start} {stop} {ref}/{alt} 1"


def extract_protein_substitution(hgvsp):
    """Extract the amino-acid substitution token (e.g. 'Leu1407Val') from a VEP hgvsp
    string like 'ENSP00000350283.3:p.Leu1407Val' (tolerant of parentheses/prefixes)."""
    if not hgvsp:
        return None
    idx = hgvsp.rfind("p.")
    tok = hgvsp[idx + 2:] if idx >= 0 else hgvsp
    return tok.strip("()[] \t")


def canonical_representation(chrom, pos, ref, alt):
    return f"chr{chrom}:{pos}:{ref}>{alt}"


def extract_c_change(s):
    """Extract the cDNA change token (e.g. '4219C>G') from an HGVS c. string."""
    if not s:
        return None
    m = re.search(r':c\.([^)\s]+)', s) or re.search(r'c\.([^)\s]+)', s)
    return m.group(1) if m else None

