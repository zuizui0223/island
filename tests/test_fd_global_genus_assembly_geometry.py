"""Pure synthetic validation of 13 real Direct endpoints and EXACT signed beta algebra."""
from __future__ import annotations

import numpy as np
import pytest
from fd_global_genus_assembly_geometry import (
    Aggregate,exact_components,source_endpoints)
from fd_island_functional_beta import TRAITS


def test_exact_signed_interaction_may_be_negative_and_never_suppressed():
    p=np.array([[1.,0.],[0.,1.],[.5,.5]])
    q=np.array([[.4,.6],[.6,.4],[.5,.5]])
    obs,genus,within,cross=exact_components(p,q)
    assert np.allclose(obs,genus+within+cross,atol=1e-12)
    assert obs[0,1]==pytest.approx(2)
    assert genus[0,1]==pytest.approx(.08)
    assert within[0,1]==pytest.approx(2.88)
    assert cross[0,1]==pytest.approx(-.96)
    assert cross[1,0]==pytest.approx(-.96)
    assert np.allclose(cross,np.transpose(cross))


def test_pair_island_block_additivity_and_correct_one_block_deletion():
    blocks=np.array([0,0,1,2,3,3])
    a=Aggregate(blocks)
    i=np.array([0,1,2,2])
    j=np.array([2,3,3,4])
    # Identity includes negative interaction; every dyad different blocks.
    d=np.array([[.2,.1,.15,-.05],
                [.3,.5,.2,-.4],
                [.6,.2,.25,.15],
                [.7,.45,.35,-.1]])
    a.add(i,j,d)
    result=a.summary(min_pairs=1,min_islands=2,min_blocks=2)
    assert result["n_species_disjoint_cross_block_pairs"]==4
    assert result["n_distinct_pair_end_islands"]==5
    assert result["n_distinct_physical_10degree_blocks"]==4
    for prefix in ("pair_equal_","island_equal_","physical_block_equal_"):
        assert result[prefix+"observed_D2"]==pytest.approx(sum(
            result[prefix+key] for key in Aggregate.KEYS[1:]))
    assert result["pair_equal_observed_D2"]==pytest.approx(.45)
    assert result["delete_one_physical_block_pair_mean_observed_D2_min"]==pytest.approx(.3)
    assert result["delete_one_physical_block_pair_mean_observed_D2_max"]==pytest.approx(.65)


def test_full_original_species_and_all_formal_flower_endpoints_no_triple_gate():
    states={t:{} for t in TRAITS}
    for t in TRAITS:
        states[t]["GenA first"]="x"
        states[t]["GenB second"]="y"
    del states["floral_symmetry"]["GenA first"]
    z=source_endpoints(states)
    assert len(z)==13
    assert z["flower_primary_color__floral_form"]=={
        "GenA first":("x","x"),"GenB second":("y","y")}
    assert z["floral_form__floral_symmetry"]=={
        "GenB second":("y","y")}
    assert len(z["floral_form"])==2
    assert len(z["floral_symmetry"])==1


def test_bad_category_algebra_and_same_block_are_rejected():
    p=np.array([[1.,0.],[0.,1.]])
    with pytest.raises(ValueError):
        exact_components(p,np.array([[1.,0.]]))
    m=Aggregate(np.array([0,0,1]))
    with pytest.raises(ValueError):
        m.add(np.array([0]),np.array([1]),np.array([[1.,1.,0.,0.]]))
    with pytest.raises(AssertionError):
        m.add(np.array([0]),np.array([2]),np.array([[1.,.4,.4,.1]]))


def test_source_region_required_to_have_all_original_ten_traits():
    with pytest.raises(ValueError):
        source_endpoints({"floral_form":{"A a":"tubular"}})
