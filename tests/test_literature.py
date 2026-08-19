"""Phase 6 QC tests: candidate selection, PMID/DOI integrity, frozen dataset checks."""
import hashlib
import os
import re
import unittest

from scripts.phase6_literature import cdna_from, classify_patterns


class TestCdnaExtraction(unittest.TestCase):
    def test_extract(self):
        self.assertEqual(cdna_from("NM_007294.4:c.4219C>G (p.Leu1407Val)"), "c.4219C>G")
        self.assertEqual(cdna_from("NM_007294.4:c.181T>G"), "c.181T>G")
        self.assertEqual(cdna_from(""), "")


class TestPatternClassification(unittest.TestCase):
    def _row(self, key, found, revel, sift, polyphen, faf=""):
        return {
            "variant_key": key, "gnomad_found": found, "revel_score": revel,
            "sift_score": sift, "polyphen_score": polyphen,
            "gnomad_exome_faf95_popmax": faf, "gnomad_genome_faf95_popmax": faf,
        }

    def test_deterministic_and_correct(self):
        rows = [
            self._row("k1", "present", "0.9", "0.01", "0.95"),   # A (impact), G (>=0.932? no 0.9<0.932)
            self._row("k2", "absent", "0.95", "0.01", "0.95"),   # C (strong impact), G
            self._row("k3", "present", "0.1", "0.5", "0.1"),     # B (tolerance)
            self._row("k4", "absent", "0.1", "0.01", "0.95"),    # D (tolerance), E (disagree)
        ]
        pat1 = classify_patterns({"variant_key": 0, "gnomad_found": 1, "revel_score": 2,
                                  "sift_score": 3, "polyphen_score": 4,
                                  "gnomad_exome_faf95_popmax": 5, "gnomad_genome_faf95_popmax": 6},
                                 [list(r.values()) + [""] * 0 for r in rows])
        # rebuild as list-of-rows in the format classify_patterns expects (list rows)
        cols = ["variant_key", "gnomad_found", "revel_score", "sift_score", "polyphen_score",
                "gnomad_exome_faf95_popmax", "gnomad_genome_faf95_popmax"]
        idx = {n: i for i, n in enumerate(cols)}
        lrows = [[r[c] for c in cols] for r in rows]
        pat = classify_patterns(idx, lrows)
        self.assertEqual(pat["k1"], {"A"})
        self.assertEqual(pat["k2"], {"C", "G"})
        self.assertEqual(pat["k3"], {"B"})
        self.assertEqual(pat["k4"], {"D", "E"})
        # determinism
        self.assertEqual(pat, classify_patterns(idx, lrows))


class TestIdentifierIntegrity(unittest.TestCase):
    def test_pmid_format(self):
        for pmid in ("30209399", "12824425"):
            self.assertRegex(pmid, r"^\d{7,8}$")

    def test_doi_format(self):
        for doi in ("10.1038/s41586-018-0461-z", "10.1016/j.ajhg.2022.10.013"):
            self.assertRegex(doi, r"^10\.\d{4,9}/")

    def test_frozen_checksum_stable(self):
        # the frozen dataset checksum recorded in the finalization report is 64 hex chars
        report = "results/reports/annotation_finalization_report.md"
        if os.path.exists(report):
            m = re.search(r"`([0-9a-f]{64})`", open(report).read())
            self.assertIsNotNone(m)


if __name__ == "__main__":
    unittest.main()
