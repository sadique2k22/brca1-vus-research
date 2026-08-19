"""Unit tests for transcript/variant normalization logic (Phase 4A)."""
import unittest

from src.variants import (
    NORMALIZED_TRANSCRIPT,
    build_variant_string,
    canonical_representation,
    extract_c_change,
    extract_protein_substitution,
    validate_alleles,
)


class TestAlleleValidation(unittest.TestCase):
    def test_valid_snp(self):
        self.assertEqual(validate_alleles("G", "C"), (True, ""))

    def test_empty(self):
        self.assertFalse(validate_alleles("", "C")[0])

    def test_same_allele(self):
        ok, reason = validate_alleles("A", "A")
        self.assertFalse(ok)
        self.assertIn("==", reason)

    def test_non_base(self):
        self.assertFalse(validate_alleles("N", "C")[0])
        self.assertFalse(validate_alleles("GG", "C")[0])


class TestRepresentation(unittest.TestCase):
    def test_variant_string(self):
        self.assertEqual(build_variant_string("17", 43082542, 43082542, "G", "C"),
                         "17 43082542 43082542 G/C 1")

    def test_canonical(self):
        self.assertEqual(canonical_representation("17", 43082542, "G", "C"),
                         "chr17:43082542:G>C")


class TestHGVSParsing(unittest.TestCase):
    def test_c_change(self):
        self.assertEqual(extract_c_change("NM_007294.4(BRCA1):c.4219C>G (p.Leu1407Val)"), "4219C>G")
        self.assertEqual(extract_c_change("ENST00000357654.9:c.4219C>G"), "4219C>G")

    def test_protein_substitution(self):
        self.assertEqual(extract_protein_substitution("ENSP00000350283.3:p.Leu1407Val"), "Leu1407Val")
        self.assertEqual(extract_protein_substitution("NP_009225.1:p.Leu1407Val"), "Leu1407Val")

    def test_transcript_constants(self):
        self.assertEqual(NORMALIZED_TRANSCRIPT, "NM_007294.4")


if __name__ == "__main__":
    unittest.main()
