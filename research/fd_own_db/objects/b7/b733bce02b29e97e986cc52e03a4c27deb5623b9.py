"""Synthetic dated-tree true-tip flora matching: no protected input required."""
import numpy as np
import pytest

from fd.true_tip_native_pairing import (
    _matched_rows, _original_pool_matrices,
    sample_true_tip_pairing,sample_same_focal_distance_bands
)
from fd.regional_native_flora_fd import METRICS


def make_source():
    a=[f"Local{i} species" for i in range(8)]
    near=[f"Near{i} donor" for i in range(8)]
    mid=[f"Mid{i} donor" for i in range(8)]
    source=near+mid
    names=a+source
    families={name:f"Family{i % 8}"
              for i, name in enumerate(a+near+mid)}
    triples={}
    for i,name in enumerate(names):
        triples[name]=(
            ("white","yellow","violet")[i%3],
            "open_radial" if i%2 else "tubular",
            "actinomorphic" if i%2 else "zygomorphic")
    D=np.full((24,24),190.,dtype=float)
    np.fill_diagonal(D,0)
    for i in range(8):
        D[i,8+i]=D[8+i,i]=35.
        D[i,16+i]=D[16+i,i]=75.
        D[8+i,16+i]=D[16+i,8+i]=40.
    return a,source,families,triples,names,D


def test_source_matching_controls_and_dated_tree_distance_bands():
    left,right,family,triple,names,D=make_source()
    a,b,dist,sf,sg=_original_pool_matrices(
        set(left),set(right),families=family,
        tip_index={n:i for i,n in enumerate(names)},D=D)
    assert a==left
    assert b==sorted(right)
    assert sf.sum()==16
    assert sg.sum()==0
    x=sample_true_tip_pairing(a,b,dist,sf,sg,triple,threshold=50.,
                              repeats=24,seed=20261114)
    y=sample_true_tip_pairing(a,b,dist,sf,sg,triple,threshold=50.,
                              repeats=24,seed=20261114)
    assert x==y
    assert x["estimable"]
    assert 0<=x["mean_original_bounded_tree_distance"]<=50
    assert all(np.isfinite(x["bounded_"+m]) for m in METRICS)
    x100=sample_true_tip_pairing(a,b,dist,sf,sg,triple,threshold=100.,
                              repeats=24,seed=20261114)
    assert x100["estimable"]
    z=sample_same_focal_distance_bands(a,b,dist,sf,triple,
                                       repeats=24,seed=20261114)
    assert z["estimable"] and z["n_completed_draws"]==24
    assert z["mean_near_tree_distance"]<=50
    assert 50<z["mean_mid_tree_distance"]<=100
    for m in METRICS:
        assert np.isclose(z["focal_minus_near_"+m]-
                          z["focal_minus_mid_"+m],
                          z["near_minus_mid_source_sensitivity_"+m])


def test_source_matching_prohibits_name_overlap_and_false_genus():
    a,b,family,traits,names,D=make_source()
    with pytest.raises(ValueError):
        _original_pool_matrices(set(a),set(b)|{a[0]},
            families=family,tip_index={x:i for i,x in enumerate(names)},D=D)
    local,source,dist,sf,sg=_original_pool_matrices(
        set(a),set(b),families=family,
        tip_index={x:i for i,x in enumerate(names)},D=D)
    z=sample_true_tip_pairing(local,source,dist,sf,sg,traits,
          threshold=100.,require_genus=True,repeats=24,seed=20261114)
    assert not z["estimable"]
    assert _matched_rows(np.zeros((8,2),bool)).tolist()==[-1]*8


def test_reject_nonfinite_and_unpublished_tree_paths():
    a,b,family,traits,names,D=make_source()
    D[0,8]=np.nan
    with pytest.raises(ValueError):
        _original_pool_matrices(set(a),set(b),
            families=family,tip_index={x:i for i,x in enumerate(names)},D=D)
