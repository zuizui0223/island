"""Exact Jaccard beta identities and fixed matched-cohort integration."""
import csv
import gzip
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from island_v2.floral_beta_partition import partition, analyze, SCOPES, REGIMES


class BetaPartitionTests(unittest.TestCase):
    def test_identical_disjoint_nested_mixed(self):
        x = partition(2, 2, 2)
        self.assertEqual((x["beta_jaccard"], x["beta_turnover"], x["beta_nestedness"]),
                         (0, 0, 0))
        self.assertIsNone(x["turnover_share"])
        self.assertEqual(x["richness_normalized_jaccard"], 1)
        x = partition(2, 2, 0)
        self.assertEqual((x["beta_jaccard"], x["beta_turnover"], x["beta_nestedness"]),
                         (1, 1, 0))
        x = partition(2, 4, 2)
        self.assertEqual((x["beta_jaccard"], x["beta_turnover"], x["beta_nestedness"]),
                         (.5, 0, .5))
        self.assertEqual(x["richness_normalized_jaccard"], 1)
        x = partition(3, 5, 2)
        self.assertAlmostEqual(x["beta_turnover"], .5)
        self.assertAlmostEqual(x["beta_nestedness"], 1/6)
        self.assertAlmostEqual(x["richness_normalized_jaccard"], 5/9)

    def test_exact_additivity_symmetry_and_bounds(self):
        for a in range(1, 10):
            for b in range(1, 10):
                for shared in range(min(a, b) + 1):
                    x, y = partition(a, b, shared), partition(b, a, shared)
                    self.assertEqual(x, y)
                    self.assertAlmostEqual(x["beta_jaccard"],
                                           x["beta_turnover"] + x["beta_nestedness"])
                    for number in x.values():
                        if number is not None:
                            self.assertGreaterEqual(number, -1e-12)
                            self.assertLessEqual(number, 1+1e-12)

    def test_rejects_bad_counts(self):
        for a,b,c in ((0,1,0), (1,0,0), (1,1,2), (1,1,-1), (1.0,1,1), (True,1,1)):
            with self.assertRaises(ValueError):
                partition(a,b,c)

    @staticmethod
    def fixture(root):
        input_path = root / "input.csv.gz"
        columns = ["evidence_scope", "analysis_regime", "island_a", "island_b",
                   "n_profiles_a", "n_profiles_b", "profile_intersection",
                   "profile_union", "profile_jaccard", "n_shared_taxa",
                   "n_unique_a", "n_unique_b"]
        with gzip.open(input_path, "wt", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=columns)
            writer.writeheader()
            for scope in SCOPES:
                for regime in REGIMES:
                    for a,b,ka,kb,overlap in (("A","B",3,5,2), ("A","C",2,4,2)):
                        writer.writerow(dict(evidence_scope=scope, analysis_regime=regime,
                                             island_a=a, island_b=b,
                                             n_profiles_a=ka, n_profiles_b=kb,
                                             profile_intersection=overlap,
                                             profile_union=ka+kb-overlap,
                                             profile_jaccard=overlap/(ka+kb-overlap),
                                             n_shared_taxa=1, n_unique_a=9, n_unique_b=9))
        return input_path

    def test_16_groups_aggregate_only(self):
        with TemporaryDirectory() as d:
            root = Path(d)
            source = self.fixture(root)
            summaries = analyze(source, root/"out")
            self.assertEqual(len(summaries), 16)
            self.assertEqual(sum(s["n_pairs"] for s in summaries), 32)
            self.assertEqual(summaries[0]["n_islands"], 3)
            with gzip.open(root/"out"/"pairwise_beta_partition.csv.gz", "rt") as f:
                rows = list(csv.DictReader(f))
            self.assertEqual(len(rows), 32)
            self.assertNotIn("accepted_species", rows[0])
            self.assertNotIn("p_two_sided", rows[0])
            self.assertIn("descriptive only", (root/"out"/"manifest.json").read_text())

    def test_bad_union_rejected(self):
        with TemporaryDirectory() as d:
            root = Path(d)
            source = self.fixture(root)
            with gzip.open(source, "rt") as f:
                rows = list(csv.DictReader(f))
            rows[0]["profile_union"] = 100
            with gzip.open(source, "wt", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(ValueError, "Inconsistent"):
                analyze(source, root/"out")
            self.assertFalse((root/"out").exists())

    def test_unmatched_all_direct_rejected(self):
        with TemporaryDirectory() as d:
            root = Path(d)
            source = self.fixture(root)
            with gzip.open(source, "rt") as f:
                rows = list(csv.DictReader(f))
            for row in rows:
                if (row["evidence_scope"] == "direct"
                        and row["analysis_regime"] == REGIMES[0]
                        and row["island_b"] == "C"):
                    row["island_b"] = "D"
            with gzip.open(source, "wt", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(ValueError, "not matched"):
                analyze(source, root/"out")


if __name__ == "__main__":
    unittest.main()
