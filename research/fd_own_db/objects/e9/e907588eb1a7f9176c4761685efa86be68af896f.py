"""Rights-safe synthetic tests for direct reproductive assurance and floral pairing."""
import numpy as np
import pandas as pd
import pytest

from fd.reproductive_assurance_functional_diversity import (
   direct_binary_assurance,_joint,_sampling
)


def test_unobserved_reproductive_states_are_missing_not_zero():
    d=pd.DataFrame([
       ("A","self_incompatibility","SC","species_direct","high"),
       ("B","self_incompatibility","SI","species_direct","high"),
       ("C","self_incompatibility","mixed_or_variable","species_direct","high"),
       ("D","self_incompatibility","SC","genus_consensus","high"),
       ("A","autonomous_selfing_capacity","autonomous","species_direct","high"),
       ("B","autonomous_selfing_capacity","absent","species_direct","medium"),
       ("C","autonomous_selfing_capacity","facilitated","species_direct","high"),
       ("C","mating_system","predominantly_outcrossing","species_direct","high"),
       ("A","mating_system","predominantly_outcrossing","species_direct","high"),
       ("B","mating_system","mixed_mating","species_direct","high"),
    ],columns=["accepted_species","trait_name","normalized_value","evidence_scope","evidence_quality"])
    x=direct_binary_assurance(d)
    assert x.loc["A","self_compatibility"]==1
    assert x.loc["B","self_compatibility"]==0
    assert pd.isna(x.loc["C","self_compatibility"])
    assert "D" not in x.index
    assert x.loc["B","autonomous_selfing"]==0
    assert pd.isna(x.loc["C","autonomous_selfing"])
    assert x.loc["A","selfing_mating_system"]==0
    assert pd.isna(x.loc["B","selfing_mating_system"])
    with pytest.raises(ValueError):
        direct_binary_assurance(pd.concat([d,d.iloc[[0]]],ignore_index=True))


def test_same_carrier_morphology_and_binary_assurance():
    y=np.array([1.]*10+[0.]*10)
    morph=np.column_stack((["red"]*10+["blue"]*10,
        ["open_radial"]*10+["tubular"]*10,
        ["actinomorphic"]*10+["zygomorphic"]*10))
    a=_sampling(y,morph,n=20,repeats=10,rng=np.random.default_rng(1))
    assert a[0]==pytest.approx(.5)
    assert a[1]==pytest.approx(.5)
    assert a[2]==pytest.approx(.5)
    assert _joint(morph)==pytest.approx(.5)
    assert _joint(morph)!=_joint(np.column_stack((
        ["red"]*5+["blue"]*5+["red"]*5+["blue"]*5,morph[:,1],morph[:,2])))
