"""Synthetic checks that GBIF effort is an observation proxy, not absence."""
from __future__ import annotations

import csv
import numpy as np
import pandas as pd
import pytest

from fd.core import ProfileSpace
from fd.effort_audit import (
    read_effort, pair_effort_label, eligible_pairs, score_subset,
)

def pair(i):
    return {"continental_id": f"c{i}", "oceanic_id": f"o{i}",
            "region": "northern_midlatitude"}

def effort(n, datasets, uncertainty):
    return {"n_unique_gbif_ids": n, "n_datasets": datasets,
            "median_coordinate_uncertainty_m": uncertainty}

def test_no_record_in_effort_does_not_count_as_true_absence():
    pairs = [pair(0), pair(1)]
    known = {"c0": effort(100, 4, 1000), "o0": effort(100, 4, 1000)}
    label, meta = pair_effort_label(pairs[1], known)
    assert label == "no_effort_metadata"
    assert not meta["has_both"]
    included, audit = eligible_pairs(pairs, known, "both_effort_available")
    assert len(included) == 1 and audit["n_pairs_retained"] == 1
    original, _ = eligible_pairs(pairs, known, "all_original_pairs")
    assert len(original) == 2

def test_calipers_and_uncertainty_reject_unequal_data_process():
    pairs = [pair(i) for i in range(4)]
    records = {
        "c0": effort(300, 4, 200), "o0": effort(280, 3, 1000),
        "c1": effort(30, 3, 2000), "o1": effort(3000, 3, 150),
        "c2": effort(150, 4, 500), "o2": effort(100, 4, 20000),
        "c3": effort(5, 1, 100), "o3": effort(500, 5, 1000),
    }
    a, audit = eligible_pairs(pairs, records, "min30records_2datasets")
    assert len(a) == 3
    b, audit = eligible_pairs(pairs, records,
                              "min30_balanced4x_uncertainty10km")
    assert len(b) == 1 and b[0]["continental_id"] == "c0"
    assert audit["n_minimum_effort_both"] == 3

def test_tie_aware_flora_difference_score_gated_by_pair_support():
    profiles = {
        0: (frozenset({"red"}), frozenset({"bell"})),
        1: (frozenset({"white"}), frozenset({"bell"})),
        2: (frozenset({"white"}), frozenset({"tubular"})),
    }
    space = ProfileSpace(profiles)
    rows = [pair(i) for i in range(8)]
    carriers = {
        **{f"c{i}": {0: 5, 1: 3} for i in range(8)},
        **{f"o{i}": {1: 5, 2: 3} for i in range(8)},
    }
    a = score_subset(rows, carriers, space, draws=999, seed=41)
    b = score_subset(rows, carriers, space, draws=999, seed=41)
    assert a == b
    assert len(a) == 2
    assert all(0 < x["p_swap_upper"] <= 1 for x in a)
    assert score_subset(rows[:7], carriers, space, draws=999, seed=41) == []

def test_read_effort_requires_exact_unique_island_rows(tmp_path):
    p = tmp_path / "effort.csv"
    cols = ("island_id", "n_unique_gbif_ids", "n_datasets", "n_species",
            "n_preserved_specimen", "n_human_observation", "n_uncertainty_reported",
            "median_coordinate_uncertainty_m", "n_coordinate_uncertainty_gt_1000m")
    rows = [dict(zip(cols, ["a", "100", "4", "30", "15", "85", "60", "900", "8"]))]
    with p.open("w", newline="") as stream:
        w = csv.DictWriter(stream, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    assert len(read_effort(p)) == 1
    with p.open("a", newline="") as stream:
        w = csv.DictWriter(stream, fieldnames=cols)
        w.writerows(rows)
    with pytest.raises(ValueError, match="repeated"):
        read_effort(p)
