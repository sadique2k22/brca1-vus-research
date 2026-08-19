"""Phase 7 tests: deterministic cohort selection, strata, frozen datasets, identifiers."""
import unittest

from scripts.phase7_cohort import aa_pos, cdna_of, classify_stratum, domain_of, neutral_summary


class TestStratification(unittest.TestCase):
    def test_classify(self):
        self.assertEqual(classify_stratum(0.7, -1.0), "A")
        self.assertEqual(classify_stratum(0.7, 0.5), "B")
        self.assertEqual(classify_stratum(0.1, -1.0), "C")
        self.assertEqual(classify_stratum(0.1, 0.5), "D")
        self.assertEqual(classify_stratum(0.5, -1.0), "E")
        self.assertIsNone(classify_stratum(None, 0.5))
        self.assertIsNone(classify_stratum(0.7, None))

    def test_domain(self):
        self.assertEqual(domain_of(50), "RING")
        self.assertEqual(domain_of(1700), "BRCT")
        self.assertEqual(domain_of(500), "other")


class TestHelpers(unittest.TestCase):
    def test_aa_pos(self):
        self.assertEqual(aa_pos("Leu1407Val"), 1407)
        self.assertIsNone(aa_pos(""))

    def test_cdna_strips_transcript(self):
        self.assertEqual(cdna_of("NM_007294.4:c.4219C>G"), "c.4219C>G")
        self.assertEqual(cdna_of("NM_007294.3:c.4219C>G"), "c.4219C>G")

    def test_no_pathogenic_language(self):
        for s in [neutral_summary("impact", "LOF", "no"), neutral_summary("impact", "normal", "yes"),
                  neutral_summary("tolerance", "normal", "no")]:
            self.assertNotIn("pathogenic", s)
            self.assertNotIn("benign", s)


if __name__ == "__main__":
    unittest.main()
