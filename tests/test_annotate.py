"""Unit tests for Phase 4B annotation logic (variant keys, dedup, gnomAD parse, ranges)."""
import unittest

from src.population import parse_gnomad_record
from src.predictors import CaddClient, revel_score
from src.variants import gnomad_variant_id, make_variant_key


class TestVariantKey(unittest.TestCase):
    def test_deterministic(self):
        k1 = make_variant_key("17", "43082542", "G", "C")
        k2 = make_variant_key("17", "43082542", "G", "C")
        self.assertEqual(k1, k2)
        self.assertEqual(k1, "GRCh38:17:43082542:G:C")

    def test_distinct(self):
        self.assertNotEqual(make_variant_key("17", "1", "G", "C"),
                            make_variant_key("17", "2", "G", "C"))

    def test_gnomad_id(self):
        self.assertEqual(gnomad_variant_id("17", "43082542", "G", "C"), "17-43082542-G-C")


class TestGnomadParse(unittest.TestCase):
    def test_absent(self):
        self.assertEqual(parse_gnomad_record({}), {"gnomad_found": "absent"})

    def test_error(self):
        self.assertEqual(parse_gnomad_record({"_error": "x"})["gnomad_found"], "error")

    def test_present_and_range(self):
        rec = {"genome": {"af": 0.014, "ac": 2177, "an": 152184, "homozygote_count": 20,
                          "faf95": {"popmax": 0.02, "popmax_population": "nfe"},
                          "populations": [{"id": "nfe", "ac": 100, "an": 10000}]},
               "exome": {"af": 0.019, "ac": 28896, "an": 1461662, "homozygote_count": 318,
                         "faf95": {"popmax": 0.023, "popmax_population": "nfe"},
                         "populations": [{"id": "nfe", "ac": 200, "an": 10000}]}}
        p = parse_gnomad_record(rec)
        self.assertEqual(p["gnomad_found"], "present")
        self.assertAlmostEqual(p["gnomad_genome_af"], 0.014)
        self.assertEqual(p["gnomad_genome_hom"], 20)
        self.assertAlmostEqual(p["gnomad_genome_faf95_popmax"], 0.02)
        for f in ("gnomad_genome_af", "gnomad_exome_af"):
            self.assertIsNotNone(p[f])
            self.assertGreaterEqual(p[f], 0.0)
            self.assertLessEqual(p[f], 1.0)


class TestRevel(unittest.TestCase):
    def test_lookup(self):
        table = {"43082542": {"G>C": 0.5, "G>A": 0.1}}
        self.assertEqual(revel_score(table, "43082542", "G", "C"), 0.5)
        self.assertIsNone(revel_score(table, "43082542", "G", "T"))
        self.assertIsNone(revel_score(table, "999", "G", "C"))


class TestCaddParse(unittest.TestCase):
    def test_parse(self):
        text = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n17\t43082542\t.\tG\tC\t.\t.\t.\n"
        out = CaddClient._parse(text, ["17:43082542:G:C"])
        self.assertEqual(out["17:43082542:G:C"], None)  # QUAL column, not a score -> None

    def test_parse_no_crash(self):
        self.assertEqual(CaddClient._parse("", ["a"]), {"a": None})


if __name__ == "__main__":
    unittest.main()
