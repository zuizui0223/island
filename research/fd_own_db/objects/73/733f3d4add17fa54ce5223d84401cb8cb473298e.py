"""Synthetic GIFT source-pool data validation; no provider data in unit tests."""
from __future__ import annotations

import pandas as pd
import numpy as np
import pytest

from fd.mainland_source_residual import (
    SOURCE_MODES, _normalize_profile_set, _source_and_observed,
    load_candidate_sources,
)


def test_distinct_source_regional_profiles_and_coverage(tmp_path):
    islands={"I1","I2"}
    assignments = pd.DataFrame([
        {"island_id":island, "source_mode":mode,
         "entity_ID":("E1" if island=="I1" else "E2"), "source_rank":1}
        for island in sorted(islands) for mode in SOURCE_MODES
    ])
    assignments_path=tmp_path/"source_assignments.csv.gz"
    assignments.to_csv(assignments_path,index=False,compression="gzip")
    source_path=tmp_path/"gift_native.csv.gz"
    pd.DataFrame({
        "entity_ID":["E1","E1","E2","E2","EXTRA"],
        "work_species":["A alpha","B beta","C gamma","D delta","E epsilon"],
    }).to_csv(source_path,index=False,compression="gzip")
    assigned,native,audit=load_candidate_sources(
        assignments_path,source_path,islands)
    assert audit["assigned_source_entities"]==2
    assert len(assigned)==8
    assert native.accepted_species.nunique()==4

    private=pd.DataFrame({
        "accepted_species":["A alpha","B beta","C gamma","D delta"],
        "genus":["A","B","C","D"],"profile_id":[0,1,2,3],
    })
    flora_path=tmp_path/"island_flora.csv.gz"
    pd.DataFrame({
        "island_id":["I1","I1","I2","I2"],
        "accepted_species":["A alpha","B beta","C gamma","D delta"],
    }).to_csv(flora_path,index=False,compression="gzip")
    observed,mainland,coverage=_source_and_observed(
        flora_path,private,native,assigned,islands,n_reference=4,
        min_source_profiles=2)
    assert len(observed)==2
    assert len(mainland)==8
    assert all(value.sum()==pytest.approx(1.0) for value in observed.values())
    assert coverage["I1","geo_k5"]["source_overlap"]==pytest.approx(1.0)
    assert coverage["I2","geo_k5"]["source_profiles"]==2


def test_source_profiles_are_normalized_and_empty_is_rejected():
    x=_normalize_profile_set({0,2},4)
    assert np.array_equal(x,np.array([0.5,0,0.5,0]))
    with pytest.raises(ValueError,match="Empty"):
        _normalize_profile_set(set(),4)
    with pytest.raises(ValueError,match="Unknown"):
        _normalize_profile_set({90},4)


def test_missing_source_mode_is_not_silently_imputed(tmp_path):
    assign=tmp_path/"ass.csv"
    pd.DataFrame([{"island_id":"I1","source_mode":"geo_k5",
                   "entity_ID":"E1","source_rank":1}]).to_csv(assign,index=False)
    native=tmp_path/"native.csv"
    pd.DataFrame([{"entity_ID":"E1","work_species":"A alpha"}]).to_csv(native,index=False)
    with pytest.raises(ValueError,match="Required fixed source-matching modes"):
        load_candidate_sources(assign,native,{"I1"})
