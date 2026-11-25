"""Synthetic tests of source-only lineage opportunity -> flower-profile mixing."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from fd.two_stage_assembly import (
    genus_flower_model, predict_island_flora,
    score_flora, source_genus_weights,
)
from fd.core import ProfileSpace


def training_fixture():
    donor = pd.DataFrame({
        "accepted_species": ["A one","A two","B one","B two","C one","D one"],
        "genus": ["A","A","B","B","C","D"],
        "profile_id": [0,0,3,3,2,1],
    })
    return donor


def source_fixture():
    source_maps = {
        "X": {"A":3,"B":1,"C":1,"D":1,"E":1},
        "Y": {"A":2,"B":3,"C":1,"D":1,"E":2},
        "Z": {"A":1,"B":3,"C":1,"D":2,"E":1},
    }
    sizes = {k:sum(v.values()) for k,v in source_maps.items()}
    global_genus = {g:sum(g in v for v in source_maps.values())
                    for g in source_maps["X"]}
    near = [("X", 30.0), ("Y", 75.0), ("Z", 125.0)]
    return source_maps,sizes,global_genus,near


def test_source_genus_weights_are_valid_and_invariant_to_island_genus():
    source,sizes,prevalence,near=source_fixture()
    genus,weights=source_genus_weights(near,source,sizes,prevalence)
    assert len(genus)==5
    assert np.isclose(weights["source_genus_density"].sum(),1)
    assert np.isclose(weights["source_genus_prevalence"].sum(),1)
    with pytest.raises(ValueError):
        source_genus_weights([("X",10.0)],source,sizes,prevalence)


def test_two_stage_source_prediction_is_blind_to_island_genus_names():
    source,sizes,prevalence,near=source_fixture()
    prior,gpdf=genus_flower_model(training_fixture(),4,shrinkage=8)
    a=pd.DataFrame([
        {"accepted_species":"A heldout1","genus":"A","profile_id":0},
        {"accepted_species":"A heldout2","genus":"A","profile_id":1},
    ])
    b=pd.DataFrame([
        {"accepted_species":"B heldout1","genus":"B","profile_id":0},
        {"accepted_species":"B heldout2","genus":"B","profile_id":1},
    ])
    p1,_=predict_island_flora(a,near,source,sizes,prevalence,prior,gpdf)
    p2,_=predict_island_flora(b,near,source,sizes,prevalence,prior,gpdf)
    for model in ("global_prior","source_genus_density",
                  "source_genus_prevalence"):
        assert np.allclose(p1[model][1],p2[model][1])
    # The post-hoc known-genus benchmark can use observed genus identity.
    assert not np.allclose(p1["known_genus_oracle"][1],
                           p2["known_genus_oracle"][1])
    assert np.allclose([p[1].sum() for p in p1.values()],1)


def test_scoring_positive_probabilities_and_original_distance():
    profiles={
        0:(frozenset({"red"}),frozenset({"bell"})),
        1:(frozenset({"blue"}),frozenset({"bell"})),
        2:(frozenset({"red"}),frozenset({"tubular"})),
        3:(frozenset({"blue"}),frozenset({"tubular"})),
    }
    dist=ProfileSpace(profiles).matrix
    obs=np.array([1,0,0,0],dtype=float)
    p=np.array([.7,.1,.1,.1])
    metrics=score_flora(obs,p,dist)
    assert metrics["name_log_loss"]==pytest.approx(-np.log(.7))
    assert metrics["frequency_L2"]==pytest.approx(np.sum((obs-p)**2))
    assert metrics["frequency_kernel"]>=0
    assert metrics["name_trait_distance"]>=0
    with pytest.raises(ValueError):
        score_flora(obs,np.array([1,0,0,0]),dist)


def test_genus_model_requires_source_donors_and_shrinkage():
    with pytest.raises(ValueError):
        genus_flower_model(training_fixture(),4,shrinkage=0)
    with pytest.raises(ValueError):
        genus_flower_model(training_fixture().iloc[:0],4,shrinkage=8)
    p0,pergenus=genus_flower_model(training_fixture(),4,shrinkage=8)
    assert abs(p0.sum()-1)<1e-10
    assert all(abs(x.sum()-1)<1e-10 for x in pergenus.values())
