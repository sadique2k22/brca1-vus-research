"""Phase 9 tests: exact-variant matching, one-letter conversion, classification."""
import unittest

from scripts.phase9 import classify, one_letter


class TestOneLetter(unittest.TestCase):
    def test_convert(self):
        self.assertEqual(one_letter("Leu1407Val"), "L1407V")
        self.assertEqual(one_letter("Asp96Glu"), "D96E")
        self.assertIsNone(one_letter(""))  # no protein change


class TestClassify(unittest.TestCase):
    def test_exact_variant(self):
        abs = "We characterized BRCA1 Leu1407Val and found reduced function."
        self.assertEqual(classify(abs, "Leu1407Val", "c.4219C>G"), "EXACT_VARIANT")

    def test_exact_variant_cdna(self):
        abs = "BRCA1 c.4219C>G was tested."
        self.assertEqual(classify(abs, "Leu1407Val", "c.4219C>G"), "EXACT_VARIANT")

    def test_exact_variant_one_letter(self):
        abs = "The L1407V variant..."
        self.assertEqual(classify(abs, "Leu1407Val", "c.4219C>G"), "EXACT_VARIANT")

    def test_gene_level(self):
        abs = "BRCA1 is a tumor suppressor gene involved in DNA repair."
        self.assertEqual(classify(abs, "Leu1407Val", "c.4219C>G"), "GENE_LEVEL")

    def test_unclear_no_abstract(self):
        self.assertEqual(classify("", "Leu1407Val", "c.4219C>G"), "UNCLEAR")


if __name__ == "__main__":
    unittest.main()
