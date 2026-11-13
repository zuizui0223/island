"""Tests for exact margin preservation and correct separation of estimands."""
import random
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from island_v2.floral_fixed_margins import (margins, curveball_trade, jaccard, _adjust,
                                           run_comparison, run_suite, write_csv, read_rows)

class FixedMarginTests(unittest.TestCase):
    def test_margin_invariance_after_thousands_of_trades(self):
        original = [{0, 1, 4}, {0, 2, 3}, {1, 3, 4}, {2, 3, 4}]
        rows = [set(x) for x in original]
        rng = random.Random(22)
        expected = margins(rows)
        self.assertTrue(any(curveball_trade(rows, rng) for _ in range(100)))
        for _ in range(1000):
            curveball_trade(rows, rng)
        self.assertEqual(margins(rows), expected)

    def test_reproducible_and_no_mutation(self):
        rows = [{0, 1}, {0, 2}, {1, 2}, {0, 3}]
        original = [set(row) for row in rows]
        a, draws_a = run_comparison(rows, [(0, 1), (0, 2), (1, 3)], seed=6, per_chain=13,
                                    burn_multiplier=2, gap_multiplier=2)
        b, draws_b = run_comparison(rows, [(0, 1), (0, 2), (1, 3)], seed=6, per_chain=13,
                                    burn_multiplier=2, gap_multiplier=2)
        self.assertEqual(a, b)
        self.assertEqual(draws_a, draws_b)
        self.assertEqual(rows, original)
        self.assertEqual(len(draws_a), 26)
        self.assertEqual(a['observed_median_jaccard'], 1/3)

    def test_known_jaccard_and_pair_validation(self):
        self.assertEqual(jaccard({0, 1}, {1, 2}), 1/3)
        with self.assertRaisesRegex(ValueError, "empty union"):
            jaccard(set(), set())
        for p in [[(0, 0)], [(0, 11)], [(-1, 1)]]:
            with self.assertRaisesRegex(ValueError, "Invalid pair"):
                run_comparison([{1}, {2}], p, seed=4, per_chain=3)

    def test_adjustment_fixed_family_and_missing_values(self):
        p = [0.01, 0.03, None, 0.70]
        out = _adjust(p)
        self.assertEqual(out[0], 0.04)
        self.assertEqual(out[1], 0.06)
        self.assertIsNone(out[2])
        self.assertGreaterEqual(_adjust(p, by=True)[0], out[0])

    def test_mismatched_all_direct_pairs_rejected(self):
        with TemporaryDirectory() as d:
            d = Path(d)
            p, root, out = d / 'prior', d / 'root', d / 'out'
            p.mkdir()
            root.mkdir()
            fields = ['island_a', 'island_b', 'evidence_scope', 'analysis_regime', 'n_shared_taxa']
            rows = []
            for scope in ['all', 'direct', 'all_single_state', 'direct_single_state']:
                for regime in ['northern_midlatitude', 'northern_high_latitude', 'tropical',
                               'southern_extratropical']:
                    rows.append(dict(zip(fields, ['I', 'J' if scope != 'direct' else 'K',
                                                  scope, regime, 0])))
            write_csv(p / 'pairwise_observed_and_null.csv.gz', fields, rows)
            with self.assertRaisesRegex(ValueError, "not matched"):
                run_suite(p, root, out, per_chain=3)
            self.assertFalse(out.exists())

    def test_deterministic_gzip_and_roundtrip(self):
        with TemporaryDirectory() as d:
            a, b = Path(d) / 'a.gz', Path(d) / 'b.gz'
            rows = [{'x': 'A', 'y': '1'}, {'x': 'B', 'y': '2'}]
            for p in [a, b]:
                write_csv(p, ['x', 'y'], rows)
            self.assertEqual(a.read_bytes(), b.read_bytes())
            self.assertEqual(a.read_bytes()[4:8], b'\0\0\0\0')
            self.assertEqual(list(read_rows(a)), rows)

if __name__ == '__main__':
    unittest.main()
