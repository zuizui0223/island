"""Synthetic-only nearest geographic floral support/reference safeguards."""
import numpy as np
import pytest

from fd.local_geography_ses import (
    _support_summary, great_circle_distance_km,
    local_pool, nearest_regional_donors, paired_reference_summary,
)
from fd.distance_ses import sample_null


def test_great_circle_haversine_and_coordinate_validation():
    dist = great_circle_distance_km(
        0, 0, np.array([0., 0., 30.]), np.array([0., 90., 0.]))
    assert dist[0] == pytest.approx(0.)
    assert dist[1] == pytest.approx(10007.557, rel=1e-4)
    assert dist[2] == pytest.approx(3335.852, rel=1e-4)
    with pytest.raises(ValueError):
        great_circle_distance_km(95, 0, np.array([0.]), np.array([0.]))
    with pytest.raises(ValueError):
        great_circle_distance_km(0, 0, np.array([91.]), np.array([0.]))


def test_nearest_donor_excludes_focal_and_other_regime_even_when_closer():
    m = {
        key: {"analysis_regime": regime}
        for key, regime in (
            ("focal", "tropical"), ("a", "tropical"),
            ("b", "tropical"), ("c", "tropical"),
            ("d", "northern_midlatitude"),
        )
    }
    geo = {
        "focal": {"island_latitude": "0", "island_longitude": "0"},
        "a": {"island_latitude": "0", "island_longitude": "1"},
        "b": {"island_latitude": "0", "island_longitude": "-1"},
        "c": {"island_latitude": "0", "island_longitude": "3"},
        "d": {"island_latitude": "0", "island_longitude": "0"},
    }
    out, radius = nearest_regional_donors("focal", "tropical", m, geo, k=2)
    assert out == ["a", "b"]  # equal distance, tie resolved by ID
    assert radius == pytest.approx(111.195, abs=.01)
    assert "focal" not in out and "d" not in out
    assert nearest_regional_donors("focal","tropical",m,geo,k=4)[0] == []


def test_local_pool_weights_donor_island_incidence_not_carrier_number():
    donors = ["a", "b", "c"]
    carriers = {
        "a": {0: 100, 1: 1}, "b": {1: 50, 2: 1},
        "c": {0: 1, 2: 1},
    }
    pool, weight = local_pool(donors, carriers, {0: 0, 1: 1, 2: 2})
    assert pool.tolist() == [0, 1, 2]
    assert weight.tolist() == [2, 2, 2]
    # A weighted null cannot manufacture new profiles outside the local pool.
    d=np.array([[0, .2, .8], [.2,0,.5], [.8,.5,0]])
    scores=sample_null(d,pool,3,draws=40,seed=11,weights=weight)
    assert np.unique(scores[:,0]).size == 1
    with pytest.raises(ValueError):
        local_pool(["a", "missing"],carriers,{0:0,1:1,2:2})


def test_small_geographic_support_cell_does_not_release_radius():
    common = dict(scope="all", geological_origin="oceanic",
                  has_neighbors=True, has_full=False, has_equal30=True,
                  n_donor_profiles=70, radius_km=150.)
    rows = [{**common, "island_id": f"private_{i}"} for i in range(4)]
    small = _support_summary(rows)
    assert len(small) == 1 and small[0]["status"] == "small_cell_suppressed"
    assert "median_10th_neighbor_distance_km" not in small[0]
    rows.append({**common, "island_id": "private_4"})
    aggregate = _support_summary(rows)[0]
    assert aggregate["median_10th_neighbor_distance_km"] == pytest.approx(150.)
    assert aggregate["n_n1c_equal30_estimable"] == 5
    assert "island_id" not in aggregate


def test_matched_n1b_n1c_summary_uses_only_same_focal_and_unchanged_observed():
    rows=[]
    for i in range(6):
        for null,ses,expected in (
            ("N1b_regional_islandfreq_LOO",-.8,.31),
            ("N1c_nearest10_geo_LOO",-.1,.26)):
            rows.append({
                "scope":"all", "mode":"equal30",
                "endpoint":"mean_nearest_profile_distance",
                "region":"tropical","geological_origin":"oceanic",
                "island_id":f"secret_{i}","n_profiles":38,
                "null":null,"ses":ses,"observed":.25,
                "null_mean":expected, "delta":.25-expected,
            })
    rows.append({
        "scope":"all", "mode":"equal30",
        "endpoint":"mean_nearest_profile_distance",
        "region":"tropical","geological_origin":"continental",
        "island_id":"unpaired","n_profiles":38,
        "null":"N1b_regional_islandfreq_LOO",
        "ses":-.5,"observed":.25,"null_mean":.3,"delta":-.05,
    })
    out=paired_reference_summary(rows)
    assert len(out)==2  # tropical plus pooled all, unpaired island omitted
    assert all(x["n_islands_both_nulls"]==6 for x in out)
    assert all(x["median_N1c_minus_N1b_SES"]==pytest.approx(.7) for x in out)
    assert all(x["median_N1c_minus_N1b_null_mean"]==pytest.approx(-.05) for x in out)
    assert all("island_id" not in x for x in out)
    bad=[dict(x) for x in rows]
    bad[1]["observed"]=.27
    with pytest.raises(AssertionError):
        paired_reference_summary(bad)


def test_unpaired_and_small_cells_emit_no_effect_estimates():
    samples=[]
    for i in range(4):
        for null in ("N1b_regional_islandfreq_LOO","N1c_nearest10_geo_LOO"):
            samples.append(dict(scope="direct",mode="equal30",
                endpoint="mean_pairwise_distance",region="northern_midlatitude",
                geological_origin="continental",island_id=f"private_{i}",
                null=null,observed=.3,ses=.1,delta=.01,null_mean=.29))
    out=paired_reference_summary(samples)
    assert len(out)==2
    assert all(x["status"]=="small_cell_suppressed" for x in out)
    assert all("median_N1c_SES" not in x for x in out)
