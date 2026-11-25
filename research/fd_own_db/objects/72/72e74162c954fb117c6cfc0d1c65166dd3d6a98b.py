import importlib.util
from pathlib import Path
import numpy as np
import pytest

P=Path(__file__).parents[1]/'src/fd/island_floral_pair_structure.py'
def module():
    assert P.exists(),'Same-mask pair-structure control not implemented'
    spec=importlib.util.spec_from_file_location('fd.pair_structure',P)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def test_identical_marginals_opposite_combinations():
    m=module();cats=[('white','open'),('white','tube'),('purple','open'),('purple','tube')]
    p=np.array([[.5,0,0,.5],[0,.5,.5,0]])
    r=m.summarize_structure(p,cats,np.array([0,1]),np.array([1,0]),minimum_islands=0)
    assert r['same_mask_first_trait_affinity']==pytest.approx(1)
    assert r['same_mask_second_trait_affinity']==pytest.approx(1)
    assert r['observed_joint_affinity']==pytest.approx(0)
    assert r['marginal_product_reference_affinity']==pytest.approx(1)
    assert r['joint_minus_product_reference']==pytest.approx(-1)

def test_joint_minus_product_is_signed_not_always_negative():
    m=module();p=np.array([[.9,.1],[.1,.9]])
    r=m.summarize_structure(p,[('w','open'),('b','tube')],np.array([0,1]),np.array([1,0]),minimum_islands=0)
    assert r['observed_joint_affinity']==pytest.approx(.6)
    assert r['marginal_product_reference_affinity']==pytest.approx(.36)
    assert r['joint_minus_product_reference']==pytest.approx(.24)
    assert sum(v['affinity_contribution'] for v in r['leading_observed_profile_contributions'])+r['unlisted_contribution']==pytest.approx(.6)

def test_independent_island_compositions_equal_reference():
    m=module();p=np.array([[.42,.18,.28,.12],[.1,.1,.4,.4]])
    r=m.summarize_structure(p,[('w','o'),('w','t'),('b','o'),('b','t')],np.array([0,1]),np.array([1,0]),minimum_islands=0)
    assert r['joint_minus_product_reference']==pytest.approx(0,abs=1e-12)

def test_no_pairs_not_biological_zero():
    m=module();p=np.array([[.5,.5]])
    r=m.summarize_structure(p,[('a','b'),('b','a')],np.array([],dtype=int),np.array([],dtype=int))
    assert r['n_directed_comparisons']==0
    assert 'observed_joint_affinity' not in r
