"""Unit tests for Phase 5 statistics/threshold logic."""
import unittest

from src.statistics import (
    categorize_polyphen,
    categorize_revel,
    categorize_sift,
    cohens_kappa,
    descriptive,
    impact_2class,
    parse_float,
    percent_agreement,
    spearman,
)


class TestParseFloat(unittest.TestCase):
    def test_absent_not_zero(self):
        # gnomAD-absent must NOT be parsed as 0.0
        self.assertIsNone(parse_float(""))
        self.assertIsNone(parse_float("NA"))
        self.assertIsNone(parse_float(None))

    def test_valid(self):
        self.assertEqual(parse_float("0.123"), 0.123)


class TestCategories(unittest.TestCase):
    def test_revel(self):
        self.assertEqual(categorize_revel(0.1), "tolerance")
        self.assertEqual(categorize_revel(0.5), "intermediate")
        self.assertEqual(categorize_revel(0.7), "impact")
        self.assertIsNone(categorize_revel(None))

    def test_sift_direction(self):
        # lower SIFT = deleterious (damaging)
        self.assertEqual(categorize_sift(0.01), "deleterious")
        self.assertEqual(categorize_sift(0.5), "tolerated")

    def test_polyphen(self):
        self.assertEqual(categorize_polyphen(0.1), "benign")
        self.assertEqual(categorize_polyphen(0.6), "possibly_damaging")
        self.assertEqual(categorize_polyphen(0.95), "probably_damaging")

    def test_impact_2class(self):
        self.assertEqual(impact_2class("revel", 0.7), "impact")
        self.assertEqual(impact_2class("revel", 0.1), "tolerance")
        self.assertIsNone(impact_2class("revel", 0.5))  # intermediate
        self.assertEqual(impact_2class("sift", 0.01), "impact")
        self.assertEqual(impact_2class("polyphen", 0.95), "impact")
        self.assertEqual(impact_2class("polyphen", 0.1), "tolerance")


class TestAgreement(unittest.TestCase):
    def test_kappa_perfect(self):
        k = cohens_kappa(["impact", "tolerance", "impact"], ["impact", "tolerance", "impact"])
        self.assertAlmostEqual(k["kappa"], 1.0)

    def test_percent_agreement(self):
        self.assertAlmostEqual(percent_agreement(["a", "b", "a"], ["a", "a", "a"]), 2 / 3)

    def test_missing_pairs_dropped(self):
        # None values are dropped from agreement
        pa = percent_agreement(["impact", None, "impact"], ["impact", "tolerance", "impact"])
        self.assertEqual(pa, 1.0)


class TestDescriptiveAndCorrelation(unittest.TestCase):
    def test_descriptive(self):
        d = descriptive([1, 2, 3, 4, 5])
        self.assertEqual(d["n"], 5)
        self.assertEqual(d["min"], 1)
        self.assertEqual(d["max"], 5)
        self.assertAlmostEqual(d["median"], 3)

    def test_spearman_handles_none(self):
        r = spearman([0.1, 0.2, None, 0.4], [0.2, 0.3, 0.9, 0.5])
        self.assertIsNotNone(r)
        self.assertEqual(r["n"], 3)


if __name__ == "__main__":
    unittest.main()
