"""Unit tests for the ClinVar VUS/missense filtering logic (src/variants.py)."""
import unittest

from src.variants import (
    classify_protein_change,
    extract_protein_change,
    filter_vus_missense,
    is_missense,
    is_vus,
)


class TestProteinChange(unittest.TestCase):
    def test_extract(self):
        self.assertEqual(extract_protein_change("NM_007294.4(BRCA1):c.4219C>G (p.Leu1407Val)"), "Leu1407Val")
        self.assertEqual(extract_protein_change("NM_007294.4(BRCA1):c.4096+1G>A"), None)

    def test_missense(self):
        self.assertEqual(classify_protein_change("(p.Leu1407Val)"), "missense")
        self.assertEqual(classify_protein_change("(p.Asn1468Asp)"), "missense")
        self.assertEqual(classify_protein_change("(p.Ser1486Cys)"), "missense")

    def test_synonymous(self):
        self.assertEqual(classify_protein_change("(p.=)"), "synonymous")
        self.assertEqual(classify_protein_change("(p.Leu1407=)"), "synonymous")

    def test_nonsense(self):
        self.assertEqual(classify_protein_change("(p.Glu23*)"), "nonsense")
        self.assertEqual(classify_protein_change("(p.Ter23Gln)"), "nonsense")

    def test_frameshift(self):
        self.assertEqual(classify_protein_change("(p.Glu23fs)"), "frameshift")
        self.assertEqual(classify_protein_change("(p.Glu23Valfs*17)"), "frameshift")

    def test_inframe_indel(self):
        self.assertEqual(classify_protein_change("(p.Pro1603_Val1607del)"), "inframe_indel")

    def test_no_protein(self):
        self.assertEqual(classify_protein_change("NM_007294.4(BRCA1):c.4096+1G>A"), "no_protein")

    def test_unknown(self):
        self.assertEqual(classify_protein_change("(p.?)"), "unknown")

    def test_is_missense_requires_snp(self):
        self.assertTrue(is_missense("single nucleotide variant", "(p.Leu1407Val)"))
        self.assertFalse(is_missense("Deletion", "(p.Leu1407Val)"))
        self.assertFalse(is_missense("single nucleotide variant", "(p.Leu1407=)"))
        self.assertFalse(is_missense("single nucleotide variant", "c.4096+1G>A"))

    def test_is_vus(self):
        self.assertTrue(is_vus("Uncertain significance"))
        self.assertFalse(is_vus("Conflicting classifications of pathogenicity"))


class TestFilter(unittest.TestCase):
    def _rec(self, vid, assembly, sig, vtype, name):
        return {"VariationID": vid, "Assembly": assembly,
                "ClinicalSignificance": sig, "Type": vtype, "Name": name}

    def test_filter_end_to_end(self):
        records = [
            # duplicate GRCh37/GRCh38 rows for same VID
            self._rec("1", "GRCh37", "Uncertain significance", "single nucleotide variant", "(p.Leu1407Val)"),
            self._rec("1", "GRCh38", "Uncertain significance", "single nucleotide variant", "(p.Leu1407Val)"),
            # missense VUS, kept
            self._rec("2", "GRCh38", "Uncertain significance", "single nucleotide variant", "(p.Asn1468Asp)"),
            # VUS but splice SNV (no protein) -> removed
            self._rec("3", "GRCh38", "Uncertain significance", "single nucleotide variant", "c.4096+1G>A"),
            # VUS but not missense (frameshift) -> removed
            self._rec("4", "GRCh38", "Uncertain significance", "single nucleotide variant", "(p.Glu23fs)"),
            # missense but Pathogenic -> removed
            self._rec("5", "GRCh38", "Pathogenic", "single nucleotide variant", "(p.Cys64Gly)"),
            # missense VUS but GRCh37 only -> removed at GRCh38 step
            self._rec("6", "GRCh37", "Uncertain significance", "single nucleotide variant", "(p.Ser1486Cys)"),
        ]
        kept, steps = filter_vus_missense(records)
        kept_vids = {r["VariationID"] for r in kept}
        self.assertEqual(kept_vids, {"1", "2"})
        # verify step remaining counts: dedup 7->6, GRCh38 6->5, VUS 5->4, missense 4->2
        self.assertEqual([s[3] for s in steps], [6, 5, 4, 2])

    def test_dedup_prefers_grch38(self):
        records = [
            self._rec("1", "GRCh37", "Uncertain significance", "single nucleotide variant", "(p.Leu1407Val)"),
            self._rec("1", "GRCh38", "Uncertain significance", "single nucleotide variant", "(p.Leu1407Val)"),
        ]
        kept, _ = filter_vus_missense(records)
        self.assertEqual(len(kept), 1)
        self.assertEqual(kept[0]["Assembly"], "GRCh38")


if __name__ == "__main__":
    unittest.main()
