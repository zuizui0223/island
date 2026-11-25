"""Synthetic-only validations: archipelago holdout is geographic, not taxon absence."""
from __future__ import annotations

import math

import numpy as np
import pytest

from fd.archipelago_cv import (
    _cluster_test, _kernel_score, donor_selection,
    haversine_km, missing_archipelago, nearest_geographic_donors,
)
from fd.core import ProfileSpace


def test_na_labels_are_never_treated_as_common_archipelago():
    for value in ("NA", " n/a ", "", "unresolved", "unknown", None, "-"):
        assert missing_archipelago(value)
    assert not missing_archipelago("Juan Fernandez")
    assert not missing_archipelago("Ryukyu")


def test_geographic_distances_and_source_pool_exclusion():
    lat = np.array([0.,0.,0.,0.,0.,0.,0.,0.])
    lon = np.array([0.,0.5,1.,1.5,0.6,1.1,1.55,2.5])
    labels=np.array(["arch1"]*4 + ["arch2"]*4)
    distance=haversine_km(lat,lon)
    assert distance.shape==(8,8)
    assert np.allclose(distance,distance.T,atol=1e-8)
    assert np.allclose(np.diag(distance),0,atol=1e-5)
    within, outside, maxwithin, gap=donor_selection(
        0,labels,distance,max_distance_ratio=3,within_radius_km=1500)
    assert len(within)==3
    assert set(within)=={1,2,3}
    assert outside is not None and len(outside)==3
    assert all(labels[j]!=labels[0] for j in outside)
    assert maxwithin>0 and gap<=math.log(3)
    n_geo=nearest_geographic_donors(0,distance)
    assert len(n_geo)==3 and 0 not in n_geo
    assert set(n_geo).issubset(set(range(1,8)))
    a,b,*_=donor_selection(0,labels,distance,within_radius_km=10)
    assert b is None


def test_no_geographically_matched_other_archipelago_if_distant():
    lat=np.zeros(8)
    lon=np.array([0,0.01,0.02,0.03,110,110.01,110.02,110.03])
    label=np.array(["same"]*4+["other"]*4)
    d=haversine_km(lat,lon)
    same, matched, radius, gap=donor_selection(
        0,label,d,max_distance_ratio=3,within_radius_km=1500)
    assert len(same)==3 and matched is None
    assert gap>math.log(3)
    with pytest.raises(ValueError,match="Insufficient"):
        donor_selection(0,np.array(["a","a","a","a","b","b","b","b"])[:3],
                        d[:3,:3])


def test_original_full_category_similarity_kernel_differs_from_exact_id():
    profiles={
        0:(frozenset({"white"}),frozenset({"bell"})),
        1:(frozenset({"red"}),frozenset({"bell"})),
        2:(frozenset({"white"}),frozenset({"tubular"})),
        3:(frozenset({"red"}),frozenset({"tubular"})),
    }
    space=ProfileSpace(profiles)
    obs=np.array([1.,0.,0.,0.])
    prediction=np.array([0.1,.4,.4,.1])
    l2=_kernel_score(obs,prediction,np.eye(4))
    original=_kernel_score(obs,prediction,1-space.matrix)
    assert l2>0
    assert original>0
    assert l2 != pytest.approx(original)
    assert _kernel_score(obs,obs,np.eye(4))==0.


def test_archipelago_is_the_only_signflip_and_bootstrap_unit():
    values=np.array([1.,2.,0.5,1.5,0.8,1.2,0.9,0.6])
    a=_cluster_test(values,permutations=1999,bootstrap=1999,seed=11)
    b=_cluster_test(values,permutations=1999,bootstrap=1999,seed=11)
    assert a==b
    assert a["mean_group_improvement"]==pytest.approx(values.mean())
    assert a["cluster_bootstrap_q025"]>0
    assert 0<a["p_cluster_signflip_one_sided"]<=1
    with pytest.raises(ValueError):
        _cluster_test(values[:7],permutations=999,bootstrap=999,seed=11)


def test_geography_only_donors_ignore_archipelago_and_may_match_exactly():
    # The three closest geographic islands can already be exactly the three
    # closest named archipelago donors; in that case added group identity
    # must have zero predictive incremental value.
    lat=np.zeros(8)
    lon=np.array([0.,.2,.4,.6,20,21,22,23])
    distance=haversine_km(lat,lon)
    labels=np.array(["A"]*4+["B"]*4)
    within, matched,*_=donor_selection(
        0,labels,distance,within_radius_km=1500)
    nearest=nearest_geographic_donors(0,distance)
    assert set(nearest)==set(within)
    with pytest.raises(ValueError,match="Too few"):
        nearest_geographic_donors(0,distance[:3,:3])
