"""Phase 6.5 tests: Findlay dataset mapping, cDNA matching, score ranges, preservation."""
import unittest

from scripts.phase6_5_findlay import aa_pos, cdna_of, classify_comparison


class TestCdnaMapping(unittest.TestCase):
    def test_transcript_stripped(self):
        # NM_007294.3 and NM_007294.4 map to the SAME c. token (CDS identical)
        self.assertEqual(cdna_of("NM_007294.3:c.4219C>G"), "c.4219C>G")
        self.assertEqual(cdna_of("ENST00000357654.9:c.4219C>G"), "c.4219C>G")
        self.assertEqual(cdna_of("NM_007294.4:c.4219C>G (p.Leu1407Val)"), "c.4219C>G")

    def test_aa_pos(self):
        self.assertEqual(aa_pos("Leu1407Val"), 1407)
        self.assertEqual(aa_pos("Asp96Glu"), 96)
        self.assertIsNone(aa_pos(""))


class TestComparisonLabels(unittest.TestCase):
    def test_no_pathogenic_terms(self):
        labels = set()
        for rv, sv, pv, fs in [(0.1, 0.01, 0.9, -1.5), (0.7, 0.5, 0.1, 0.3),
                                (0.7, 0.01, 0.9, -1.0), (0.1, 0.5, 0.1, 0.5)]:
            labels.add(classify_comparison(rv, sv, pv, fs))
        for l in labels:
            self.assertNotIn("pathogenic", l)
            self.assertNotIn("benign", l)

    def test_disagreement_classes(self):
        self.assertEqual(classify_comparison(0.1, 0.5, 0.1, -1.5),
                         "computational_tolerance + functional_LOF")
        self.assertEqual(classify_comparison(0.7, 0.5, 0.1, 0.3),
                         "computational_impact + functional_normal")
        self.assertEqual(classify_comparison(0.7, 0.01, 0.9, -1.0),
                         "computational_impact + functional_LOF (agreement)")


class TestScoreHandling(unittest.TestCase):
    def test_missing_is_empty(self):
        # score fields for unscored variants are empty strings, not 0
        self.assertIsNone(_to_float(""))
        self.assertIsNone(_to_float("NA"))


def _to_float(v):
    try:
        return float(v)
    except (ValueError, TypeError):
        return None


if __name__ == "__main__":
    unittest.main()
