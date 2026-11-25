"""Geological-origin and LGM connectivity inference safeguards."""
import math
import numpy as np
import pytest

from fd.origin import crosswalk, xyz
from fd.origin_analysis import _ols_origin_effect, descriptive


def _island(island_id, lat, lon, area=100.):
    return {"id": island_id, "lat": lat, "lon": lon, "area": area}


def test_distance_and_area_screen_prevent_nearest_only_false_matches():
    source = [_island("exact", 10.01, 20.01, 102.),
              _island("near_wrong_area", 10.011, 20.011, 1200.),
              _island("other", -30., -40., 80.)]
    target = [_island("A", 10., 20., 100.),
              _island("B", 5., 5., 100.)]
    result, audit = crosswalk(target, source)
    assert set(result) == {"A"}
    assert result["A"]["source"]["id"] == "exact"
    assert result["A"]["area_ratio"] < 1.1
    assert audit["no_candidate"] == 1


def test_ambiguous_and_reused_geography_is_not_forced_into_origin():
    sources = [_island("1", 10., 20.0), _island("2", 10., 20.002)]
    m, report = crosswalk([_island("X", 10., 20.001)], sources)
    assert not m
    assert report["ambiguous"] == 1
    m, report = crosswalk([_island("A", 10., 20.0), _island("B", 10.001, 20.0)],
                           [_island("1", 10., 20.0)])
    assert not m and report["source_reused"] == 2


def test_xyz_angles_and_origin_descriptive_never_auto_assign():
    q = xyz([0, 0, 90], [0, 90, 0])
    assert np.allclose(q @ q.T, np.eye(3), atol=1e-12)
    assert math.isfinite(q[0, 0])
    rows = [
        {"island_id":"a","analysis_regime":"northern_midlatitude","geological_origin":"oceanic",
         "n_names":"100","joint_coverage":"0.6","rao_q":"0.4","rao_q_profile_equal":"0.5",
         "kde90_area_pc12_equal_h035":"13", "distance_to_continent_km":"70",
         "spatial_block":"x"},
        {"island_id":"b","analysis_regime":"northern_midlatitude","geological_origin":"continental",
         "n_names":"100","joint_coverage":"0.7","rao_q":"0.5","rao_q_profile_equal":"0.6",
         "kde90_area_pc12_equal_h035":"12", "distance_to_continent_km":"10",
         "spatial_block":"x"},
        {"island_id":"c","analysis_regime":"northern_midlatitude","geological_origin":"unresolved",
         "n_names":"100","joint_coverage":"0.6","rao_q":"0.6","rao_q_profile_equal":"0.7",
         "kde90_area_pc12_equal_h035":"11", "distance_to_continent_km":"50",
         "spatial_block":"y"},
    ]
    a, b = descriptive(rows, [], key="geological_origin", first="continental",
                       second="oceanic", scope="all")
    assert next(r["n_islands"] for r in a if r["regime"] == "all"
                and r["group"] == "unresolved") == 1
    assert all(r["n_pairs"] == 0 for r in b)


def test_adjusted_cluster_bootstrap_requires_two_origin_classes():
    rows = []
    for i in range(100):
        region = "northern_midlatitude" if i % 2 else "tropical"
        rows.append({"geological_origin": "oceanic" if i < 50 else "continental",
                     "rao_q_profile_equal": 0.4 + 0.08 * int(i < 50) + 0.0004 * i,
                     "area_km2": 100 + (i % 10),
                     "distance_to_continent_km": 10 + (i % 13),
                     "n_names": 100 + (i % 7), "joint_coverage": 0.7,
                     "analysis_regime": region, "spatial_block": str(i % 12),
                     **{f"climate_pc{j}": ((i * (j+1)) % 17) / 17
                        for j in range(1, 5)}})
    fit = _ols_origin_effect(rows, "geological_origin", "continental",
                             "oceanic", "rao_q_profile_equal", seed=11)
    assert fit["n_first"] == 50 and fit["n_second"] == 50
    assert fit["status"] in ("exploratory_adjusted_wild_cluster_bootstrap",)
    assert fit["coefficient"] > 0
    rows[:60] = [dict(r, geological_origin="unresolved") for r in rows[:60]]
    incomplete = _ols_origin_effect(rows, "geological_origin", "continental",
                                    "oceanic", "rao_q_profile_equal", seed=11)
    assert incomplete["status"] == "insufficient_matched_support"
