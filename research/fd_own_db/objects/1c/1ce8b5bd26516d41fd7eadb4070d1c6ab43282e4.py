"""Tests for matched, distinct origin cohorts and repertoire statistical contracts."""
from __future__ import annotations

import numpy as np
import pytest

from fd.repertoire import (
    PROFILE_FAMILY, _binary_profile_tests, _paired_swap_statistic,
    _state_contrasts, _two_sided, _upper_tail, match_origin_islands,
)


def island(island_id: str, region: str = "northern_midlatitude",
           area: float = 100.0, isolation: float = 45.0,
           names: int = 120, coverage: float = 0.40):
    return {
        "island_id": island_id,
        "analysis_regime": region,
        "area_km2": str(area),
        "distance_to_continent_km": str(isolation),
        "n_names": str(names),
        "joint_coverage": str(coverage),
    }


def geography(origin: str, *, climate3: float = 0.0):
    return {"geological_origin": origin,
            "climate_pc1": "0.1", "climate_pc2": "-0.2",
            "climate_pc3": str(climate3), "climate_pc4": "0.3"}


def test_geological_matching_is_1_to_1_and_floral_outcome_blind():
    rows = [island(f"c{i}", area=95+i, isolation=40+i)
            for i in range(12)]
    rows += [island(f"o{i}", area=97+i, isolation=42+i)
             for i in range(12)]
    crosswalk = {f"c{i}": geography("continental") for i in range(12)}
    crosswalk.update({f"o{i}": geography("oceanic") for i in range(12)})
    # Direct is nested inside All; profiles / Rao Q do not enter matching.
    a, diag = match_origin_islands(rows, rows, crosswalk)
    assert len(a) == 12 and diag["n_main_pairs"] == 12
    assert len({p["continental_id"] for p in a}) == 12
    assert len({p["oceanic_id"] for p in a}) == 12
    assert all(p["region"] == "northern_midlatitude" for p in a)


def test_crosswalk_requires_calipers_and_region_overlap():
    rows = [island(f"c{i}") for i in range(9)]
    rows += [island(f"o{i}", area=60000) for i in range(9)]
    crosswalk = {f"c{i}": geography("continental") for i in range(9)}
    crosswalk.update({f"o{i}": geography("oceanic") for i in range(9)})
    with pytest.raises(ValueError, match="No geographic regime"):
        match_origin_islands(rows, rows, crosswalk)


def test_paired_signed_permutation_and_two_sided_ties():
    differences = np.array([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
    signs = np.array([[1,1,1], [-1,1,1], [1,-1,1], [-1,-1,-1]], dtype=float)
    observed, null = _paired_swap_statistic(differences, np.eye(2), signs)
    assert observed == pytest.approx(1.0)
    assert null[0] == pytest.approx(1.0)
    assert null[1] == pytest.approx(1/9)
    assert 0 < _upper_tail(observed, null) <= 1
    assert _two_sided(0.0, np.zeros(30)) == 1.0


def test_profile_family_does_not_shrink_to_favorable_results():
    labels = [(i, "red", "bell") for i in range(PROFILE_FAMILY)]
    oceanic = np.zeros((8, PROFILE_FAMILY), dtype=bool)
    continental = np.zeros_like(oceanic)
    oceanic[:, 1] = True
    continental[:, 2] = True
    records = _binary_profile_tests(oceanic, continental, labels, "all")
    assert len(records) == PROFILE_FAMILY
    assert records[1]["n_oceanic_only"] == 8
    assert records[1]["p_exact_pair_swap"] == pytest.approx(2/256)
    assert records[1]["q_bh_496"] >= records[1]["p_exact_pair_swap"]
    assert records[0]["p_exact_pair_swap"] == 1.0


def test_state_records_derive_from_recorded_state_sets():
    profiles = {
        0: (frozenset({"red", "pink"}), frozenset({"bell"})),
        1: (frozenset({"white"}), frozenset({"tubular"})),
    }
    oc = np.array([[True, False], [False, True]], dtype=bool)
    co = np.array([[False, True], [False, True]], dtype=bool)
    result = _state_contrasts(oc, co, profiles, [0, 1], "all")
    red = next(r for r in result if r["dimension"] == "color" and r["state"] == "red")
    assert red["oceanic_fraction"] == 0.5
    assert red["continental_fraction"] == 0
    assert len(result) == 5


def test_geographic_source_neighborhood_rejects_faraway_matched_pairs():
    rows = [island(f"c{i}") for i in range(10)]
    rows += [island(f"o{i}") for i in range(10)]
    labels = {
        **{f"c{i}": {**geography("continental"),
                     "island_latitude": "15.0",
                     "island_longitude": str(20 + i * 0.02)}
           for i in range(10)},
        **{f"o{i}": {**geography("oceanic"),
                     "island_latitude": "15.0",
                     "island_longitude": str(20 + i * 0.02 if i < 9 else -140)}
           for i in range(10)},
    }
    unrestricted, _ = match_origin_islands(rows, rows, labels)
    nearby, audit = match_origin_islands(
        rows, rows, labels, max_pair_geographic_km=500)
    assert len(unrestricted) == 10
    assert len(nearby) == 9
    assert all(x["oceanic_id"] != "o9" for x in nearby)
    assert audit["distance_calipers"]["max_pair_geographic_km"] == 500
