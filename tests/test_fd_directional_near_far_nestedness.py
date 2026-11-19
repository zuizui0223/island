"""Synthetic unique original category masks and proper directional nestedness null gates."""
import numpy as np
import pytest
from collections import Counter
from fd_directional_near_far_nestedness import (
    finite_all_allowed_probability,elementary_symmetric,uniform_fixed_k,
    weighted_fixed_k,genus_subset_expectation,outcome_blind_matching,
    report_matching,
)


def test_equal_k_prohibits_strict_set_subset_except_exact_equality():
    far={1,2}
    near={1,2}
    assert far<=near and not far<near
    assert uniform_fixed_k(4,2,2)==pytest.approx(1/6)
    assert uniform_fixed_k(4,3,2)==pytest.approx(.5)
    assert uniform_fixed_k(4,1,2)==0.
    with pytest.raises(ValueError):uniform_fixed_k(4,5,2)


def test_weighted_k_exact_combination_not_n_species_replacement():
    weights=(.25,.25,.25,.25)
    near=(1<<0)|(1<<1)|(1<<2)
    assert weighted_fixed_k(near,2,weights)==pytest.approx(.5)
    weights=(.7,.1,.1,.1)
    assert weighted_fixed_k(near,2,weights)==pytest.approx(.625)
    assert elementary_symmetric((.5,.5),1)==pytest.approx(1.)
    assert elementary_symmetric((.5,.5),2)==pytest.approx(.25)
    assert elementary_symmetric((.1,),2)==0
    with pytest.raises(ValueError):weighted_fixed_k(1<<5,1,weights)


def test_genus_conditional_all_pair_profile_subset_exact_without_replacement():
    near=(1<<0)
    reference={
        "A":Counter({("tube","zygo"):1,("open","radial"):1}),
        "B":Counter({("tube","zygo"):3}),
    }
    code={("tube","zygo"):0,("open","radial"):1}
    assert genus_subset_expectation(
        near,Counter({"A":1,"B":1}),reference,code)==pytest.approx(.5)
    assert genus_subset_expectation(
        near,Counter({"A":2,"B":1}),reference,code)==0
    assert finite_all_allowed_probability(5,3,2)==pytest.approx(.3)
    with pytest.raises(ValueError):
        finite_all_allowed_probability(3,4,2)


def test_identical_outcome_blind_spatial_matching_not_an_independent_history():
    rows=[]
    for i,(a,b) in enumerate(((0,1),(0,2),(1,3),(2,3),(3,4),(4,5))):
        rows.append(dict(block_near=a,block_far=b,
            kfar=1,knear=2,
            subset=1 if i%2 else 0,
            strict_subset=1 if i%2 else 0,
            uniform=.5,weighted=.6,
            genus_null=.3))
    source=outcome_blind_matching(rows,draws=5)
    assert len(source["selections"])==5
    assert source["graph"]["n_candidate_blocks"]==6
    assert all(len(set(x))==len(x) for x in source["selections"])
    assert all(0<=i<len(rows) for i in source["unique"])
    assert len(source["unique"])>0
    out=report_matching(rows,source)
    assert out["n_outcome_blind_block_pair_reselections"]==5
    assert out["n_disjoint_block_pairs_minimum"]>=1
    assert out["n_disjoint_block_pairs_maximum"]<=3
    assert out["genus_conditional_all_pair_expected_subset_probability_median"]==pytest.approx(.3)
    assert out["exact_uniform_fixed_k_reference_median"]==pytest.approx(.5)
    assert out["exact_weighted_fixed_k_reference_median"]==pytest.approx(.6)


def test_matching_is_outcome_blind_when_corresponding_score_swapped():
    rows=[dict(block_near=a,block_far=b,
               kfar=1,knear=2,subset=0,strict_subset=0,
               uniform=.2,weighted=.2,genus_null=.3)
          for a,b in ((0,1),(0,2),(2,3),(3,4))]
    a=outcome_blind_matching(rows,draws=7)
    for x in rows:
        x["subset"]=1
        x["strict_subset"]=1
        x["genus_null"]=.8
    b=outcome_blind_matching(rows,draws=7)
    assert all(np.array_equal(x,y) for x,y in zip(a["selections"],b["selections"]))
