"""Synthetic checks: block-shared origin pairs are not independent swaps."""
from __future__ import annotations

import numpy as np
import pytest

from fd.spatial_clusters import (
    connected_block_components, cluster_permutation,
    delete_component_sensitivity,
)


def _pairs(n: int):
    return [{"continental_id":f"c{i}","oceanic_id":f"o{i}",
             "region":"northern_midlatitude"} for i in range(n)]


def test_shared_spatial_blocks_transitively_connect_physical_pairs():
    pairs=_pairs(8)
    geography={}
    for i in range(8):
        geography[f"c{i}"]={"spatial_block":f"continent_{i}"}
        geography[f"o{i}"]={"spatial_block":f"oceanic_{i}"}
    # 0~1 through continental side, 1~2 through oceanic side, 4~5
    geography["c1"]["spatial_block"]="continent_0"
    geography["o2"]["spatial_block"]="oceanic_1"
    geography["o5"]["spatial_block"]="oceanic_4"
    comps=connected_block_components(pairs,geography)
    assert sorted([len(c) for c in comps],reverse=True)==[3,2,1,1,1]
    assert sum(len(c) for c in comps)==8
    del geography["c7"]
    with pytest.raises(ValueError,match="Missing geographic"):
        connected_block_components(pairs,geography)


def test_exhaustive_cluster_signflip_has_sign_symmetry_minimum():
    rng=np.random.default_rng(9)
    pair_delta=np.vstack([
        np.array([2.,0.,0.]) for _ in range(8)
    ])
    groups=[[0,1,2],[3,4],[5],[6],[7]]
    res=cluster_permutation(pair_delta,groups,np.eye(3))
    assert res["n_exact_component_sign_configurations"]==32
    assert res["exact_component_upper_p"]==pytest.approx(2/32)
    assert res["minimum_attainable_p_due_sign_symmetry"]==pytest.approx(2/32)
    assert res["observed"]==pytest.approx(4.)
    with pytest.raises(ValueError,match="partition"):
        cluster_permutation(pair_delta,[[0,1,2],[3,4],[5],[6],[6]],
                            np.eye(3))


def test_leave_one_cluster_out_does_not_reclassify_remaining_pairs():
    rng=np.random.default_rng(79)
    differences=rng.normal(size=(37,5))
    comps=[list(range(23)),list(range(23,34)),[34],[35],[36]]
    rows=delete_component_sensitivity(
        differences,comps,np.eye(5),
        scope="all",subset="all_37",endpoint="exact_profile_L2")
    assert len(rows)==5
    rem_largest=next(r for r in rows if r["component_size_removed"]==23)
    assert rem_largest["n_pairs_remaining"]==14
    assert rem_largest["n_components_remaining"]==4
    assert all(r["ratio_norm_remaining_to_original"]>=0 for r in rows)
    assert all(-1.0000001<=r["cosine_with_original_mean_direction"]<=1.0000001
               for r in rows)


def test_weigelt_archipelago_independent_grouping_does_not_join_unresolved():
    pairs=_pairs(8)
    geo={}
    for i in range(8):
        geo[f"c{i}"]={"spatial_block":f"cb{i}",
                        "weigelt_archipelago":"Japan" if i<3 else ""}
        geo[f"o{i}"]={"spatial_block":f"ob{i}",
                        "weigelt_archipelago":("West Indies" if 3<=i<6 else "")}
    comps=connected_block_components(pairs,geo,field="weigelt_archipelago")
    assert sorted(len(c) for c in comps)==[1,1,3,3]
    # Published Weigelt table uses the literal string NA (not a blank) in
    # some missing source-archipelago cells. These cannot be a common block.
    geo["c6"]["weigelt_archipelago"]="NA"
    geo["o7"]["weigelt_archipelago"]="NA"
    assert sorted(len(c) for c in connected_block_components(
        pairs,geo,field="weigelt_archipelago"))==[1,1,3,3]
    with pytest.raises(ValueError,match="Invalid"):
        connected_block_components(pairs,geo,field="fake_geography")
