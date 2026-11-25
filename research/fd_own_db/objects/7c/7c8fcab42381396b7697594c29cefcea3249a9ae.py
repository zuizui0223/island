"""Synthetic k-preserving rare/common categorical support tests; no user taxa."""
from itertools import product

import numpy as np
import pandas as pd
import pytest

from fd.rare_profile_decomposition import (
    BIN_NAMES,OUTCOMES,RarityReference,_taxonomy_eligibility
)
from fd.discrete_functional_volume import TRAITS


def _synthetic_realistic_profile_frequency():
    profiles=list(product(
         ("white","purple","yellow","red","blue","pink","green","orange","brown","cream"),
         ("open_radial","tubular","zygomorphic","composite","brush"),
         ("actinomorphic","zygomorphic")))
    assert len(profiles)==100
    # Fixed reference exactly 2234 distinct species; profile counts range
    # 1..45 with a predeclared frequency ordering and tie thresholds.
    weights=np.arange(1,101,dtype=float)**1.2
    weights=weights/weights.sum()
    baseline=np.ones(100,dtype=int)
    values=np.floor(weights*(2234-100)).astype(int)
    baseline+=values
    for j in range(2234-int(baseline.sum())):
        baseline[j%100]+=1
    assert baseline.sum()==2234
    names=[f"Genus{i%80} sp_{i}" for i in range(2234)]
    p=[]
    for x,num in zip(profiles,baseline):
        p.extend([x]*int(num))
    return pd.DataFrame(p,index=names,columns=TRAITS)


def test_rarity_bins_include_source_ties_and_are_disjoint():
    triple=_synthetic_realistic_profile_frequency()
    r=RarityReference.from_original(triple,null_draws=199,seed=17)
    assert sum(r.n_bins.values())==100
    assert min(r.n_bins.values())>=10
    assert r.low_quartile_species_count<r.high_quartile_species_count
    assert all(r.expected_uniform[k].shape==(3,) for k in range(1,31))
    assert all(r.expected_speciesfreq[k].shape==(3,) for k in range(1,31))
    assert all(r.bins[k]&r.bins[j]==0 for i,k in enumerate(BIN_NAMES) for j in BIN_NAMES[i+1:])
    assert r.profile_n_species.sum()==2234


def test_same_k_with_changed_profile_placement_has_different_rarity_coverage():
    triple=_synthetic_realistic_profile_frequency()
    r=RarityReference.from_original(triple,null_draws=199,seed=7)
    low_indices=np.where(r.profile_n_species==r.profile_n_species.min())[0]
    rare_indices=np.argsort(r.profile_n_species)[:5]
    common_indices=np.argsort(r.profile_n_species)[-5:]
    rare=np.asarray([r.volume.profiles[int(i)] for i in rare_indices],dtype=str)
    common=np.asarray([r.volume.profiles[int(i)] for i in common_indices],dtype=str)
    a,ka=r.measure(rare)
    b,kb=r.measure(common)
    assert ka==kb==5
    assert a.shape==(3,) and b.shape==(3,)
    assert not np.allclose(a,b,atol=1e-12)


def test_k_fixed_but_weighted_null_uses_different_reference_probabilities():
    triple=_synthetic_realistic_profile_frequency()
    r=RarityReference.from_original(triple,null_draws=199,seed=51)
    assert not np.allclose(r.expected_uniform[10],r.expected_speciesfreq[10])
    arr=triple.to_numpy(dtype=str)
    score=r.sample(arr,n=30,repeats=20,seed=52)
    assert set(score)==set(OUTCOMES)
    assert 0<=score["covered_rare"]<=1
    assert 0<=score["covered_common"]<=1
    assert score==r.sample(arr,n=30,repeats=20,seed=52)


def test_source_genus_matching_is_a_data_support_gate_not_causal_null():
    local={"X a","X b","Y c","Z d","Q e"}
    source={"X f","X g","X h","Y i","A j","B k"}
    assert _taxonomy_eligibility(local,source)==3
    assert _taxonomy_eligibility(local,source-set(["X f","X g","X h"]))==1
