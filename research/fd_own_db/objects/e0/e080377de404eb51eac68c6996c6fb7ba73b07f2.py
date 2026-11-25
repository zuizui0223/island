"""Synthetic matching and exact sign-flip checks, no private island identifiers."""
from __future__ import annotations

import numpy as np
import pytest

from fd.geography_independent import (
    _geographic_resources, _missing_label, _exact_pair_null,
    independent_pair_matching, _build_edges,
)


def edge(c, o, cb, ob, ca="NA", oa="NA", cost=1.0, region="tropical"):
    return {
        "continental_id": c, "oceanic_id": o, "region": region,
        "c_spatial_block": cb, "o_spatial_block": ob,
        "c_archipelago": ca, "o_archipelago": oa,
        "cost": cost,
    }


def test_source_labels_unknown_do_not_merge_distinct_pairs():
    e1 = edge("c1","o1","b1","b2")
    e2 = edge("c2","o2","b3","b4")
    assert _missing_label("NA")
    assert _missing_label(None)
    assert _missing_label("unresolved")
    assert not (_geographic_resources(e1,"both") &
                _geographic_resources(e2,"both"))
    e3=edge("c3","o3","b9","b10","SameArchipelago","Other")
    e4=edge("c4","o4","b11","b12","SameArchipelago","Another")
    assert _geographic_resources(e3,"both") & _geographic_resources(e4,"both")
    with pytest.raises(ValueError):
        _geographic_resources(e1,"unrecognized")


def test_milp_maximizes_independent_pairs_before_cost():
    # This crafted graph requires choosing the two disjoint edges, even
    # if one single edge would have lower cost.
    data=[
        edge("c1","o1","b1","b2", "Arch1","Arch2",cost=.1),
        edge("c1","o2","b1","b3", "Arch1","Arch3",cost=4.1),
        edge("c2","o1","b4","b2", "Arch4","Arch2",cost=3.7),
    ]
    for grouping in ("spatial_block","weigelt_archipelago","both","both_verified"):
        chosen,audit=independent_pair_matching(data,grouping)
        assert len(chosen)==2
        assert audit["solver_status"]==0
        assert sorted((x["continental_id"],x["oceanic_id"]) for x in chosen)==[
            ("c1","o2"),("c2","o1")
        ]


def test_exact_enumeration_has_correct_2_over_2n_resolution():
    n=8
    d=np.zeros((n,3),dtype=float)
    d[:,0]=1.0
    k=np.eye(3)
    result=_exact_pair_null(d,k)
    assert result["n_exact_sign_configurations"]==2**n
    assert result["observed"]==pytest.approx(1)
    assert result["p_exact_upper"]==pytest.approx(2/(2**n))
    assert result["minimum_attainable_p"]==pytest.approx(2/(2**n))
    with pytest.raises(ValueError):
        _exact_pair_null(d[:7],k)


def test_region_restricted_eligibility_uses_original_calipers():
    geography={}
    all_g0=[]
    for i in range(12):
        for origin in ("continental","oceanic"):
            island=f"{origin}_{i}"
            geography[island]={
                "geological_origin":origin,
                "spatial_block":f"block_{island}",
                "weigelt_archipelago":f"arch_{island}",
                "climate_pc1":0.1,
                "climate_pc2":-0.2,
                "climate_pc3":0.3,
                "climate_pc4":0.4,
            }
            all_g0.append({
                "island_id":island,
                "analysis_regime":"tropical",
                "area_km2":"80",
                "distance_to_continent_km":"30",
                "n_names":"100",
                "joint_coverage":"0.5",
            })
    edges,audit=_build_edges(all_g0,all_g0,geography)
    assert len(edges)==144
    matched,result=independent_pair_matching(edges,"both")
    assert len(matched)==12
    assert result["n_tropical"]==12
    assert result["n_distinct_spatial_blocks_used"]==24


def test_verified_archipelago_eligibility_discards_unknown_pairs():
    data=[
        edge("c1","o1","b1","b2","Arch1","Arch2"),
        edge("c2","o2","b3","b4","Arch3","Arch4"),
        edge("c3","o3","b5","b6","NA","Arch6"),
    ]
    selected,audit=independent_pair_matching(data,"both_verified")
    assert len(selected)==2
    assert audit["n_candidate_edges"]==2
    assert audit["n_unknown_archipelago_endpoints"]==0
    assert all(x["continental_id"]!="c3" for x in selected)
    with pytest.raises(ValueError,match="Cannot verify"):
        _geographic_resources(data[-1],"both_verified")


def test_verified_archipelagos_exclude_unknown_not_as_independent_source():
    samples = [
        edge("c1", "o1", "b1", "b2", "Japan", "Izu", 0.4),
        edge("c2", "o2", "b3", "b4", "NA", "Other", 0.3),
        edge("c3", "o3", "b5", "b6", "Hawaii", "Madeira", 0.5),
    ]
    nonverified, audit = independent_pair_matching(samples, "both")
    assert len(nonverified) == 3
    verified, audit = independent_pair_matching(samples, "both_verified")
    assert len(verified) == 2
    assert audit["n_candidate_edges"] == 2
    assert audit["n_unknown_archipelago_endpoints"] == 0
    assert all(x["continental_id"] != "c2" for x in verified)
