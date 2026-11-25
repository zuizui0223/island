"""No original sites or species; exact-study/year and no fabricated mainland checks."""
import json

import pandas as pd
import pytest

from fd.euppollnet_island_mainland import (
    haversine_km, original_same_plant_island_mainland_opportunity
)
from fd.euppollnet_island_overlap import SOURCE_COLUMNS


def synthetic():
    r=[]
    for study,net,date,lat,lon,plant in (
        ("s","island","2019-06-01",36.0,20.0,"Plant A"),
        ("s","mainland","2019-06-01",36.5,20.0,"Plant A"),
        ("s","island","2019-06-02",36.0,20.0,"Plant B"),
        ("another","mainland","2019-06-02",36.5,20.0,"Plant B"),
        ("s","island","2018-06-01",36.0,20.0,"Plant C"),
        ("s","mainland","2020-06-01",36.5,20.0,"Plant C"),
    ):
        r.append((study,net,date,lat,lon,plant,"SPECIES","Bombus a",
                  "Hymenoptera","Bombus","Apidae"))
    visits=pd.DataFrame(r,columns=SOURCE_COLUMNS)
    sites=pd.DataFrame([
        ("s","island",36.0,20.0,"fd_island","i1","tropical","block1"),
        ("s","mainland",36.5,20.0,"mainland",None,None,None),
        ("another","mainland",36.5,20.0,"mainland",None,None,None)
    ],columns=["Study_id","Network_id","Latitude","Longitude",
               "habitat_class","island_id","analysis_regime","spatial_block"])
    flora=pd.DataFrame([
        ("i1","Plant A","native"),("i1","Plant B","unresolved"),
        ("i1","Plant C","native")
    ],columns=["island_id","accepted_species","origin_status"])
    meta=pd.DataFrame([
        ("s",2018,"Focal_observation","Hymenoptera, Diptera"),
        ("s",2019,"Focal_observation","Hymenoptera, Diptera"),
        ("s",2020,"Focal_observation","Hymenoptera, Diptera"),
        ("another",2019,"Focal_observation","Hymenoptera, Diptera")
    ],columns=["Study_id","Year","Sampling_method","Taxa_recorded"])
    return visits,sites,flora,meta


def test_same_study_year_exact_plant_mainland_500km():
    vis,loc,flora,meta=synthetic()
    r=original_same_plant_island_mainland_opportunity(vis,loc,flora,meta)
    assert r["n_same_study_year_plant_candidate_groups_any_distance"]==1
    assert r["n_pair_groups_nearest_mainland_at_most_500_km"]==1
    assert r["n_pair_groups_nearest_mainland_at_most_250_km"]==1
    assert r["n_eligible_strict_native_pair_groups_500km"]==1
    assert r["n_eligible_islands_500km"]==1
    assert r["n_eligible_spatial_blocks_500km"]==1
    assert r["status"]=="SAME_PLANT_ISLAND_MAINLAND_COMPARISON_SUPPORT_INSUFFICIENT"
    for s in ("island_id","accepted_species","Plant_accepted_name","Study_id",
              "Latitude","Longitude","study_id","plant"):
        assert '"'+s+'"' not in json.dumps(r)
    assert r["actual_pollinator_replacement_measured"] is False


def test_missing_verified_mainland_not_negative_control():
    v,loc,f,m=synthetic()
    loc.loc[loc.habitat_class.eq("mainland"),"habitat_class"]="unresolved"
    r=original_same_plant_island_mainland_opportunity(v,loc,f,m)
    assert r["n_matching_same_study_year_species_island_mainland"]==0
    assert r["status"]=="NO_SAME_STUDY_ISLAND_MAINLAND_ACTUAL_VISIT_SUPPORT"


def test_different_source_study_year_and_species_not_paired():
    v,loc,f,m=synthetic()
    v=v.loc[v.Plant_accepted_name.ne("Plant A")].copy()
    r=original_same_plant_island_mainland_opportunity(v,loc,f,m)
    assert r["n_matching_same_study_year_species_island_mainland"]==0


def test_great_circle_and_site_ambiguity_fail():
    v,loc,f,m=synthetic()
    assert 50<haversine_km(36.,20.,pd.Series([36.5]).to_numpy(),pd.Series([20.]).to_numpy())[0]<60
    loc=pd.concat((loc,loc.iloc[[0]]),ignore_index=True)
    with pytest.raises(ValueError,match="ambiguous"):
        original_same_plant_island_mainland_opportunity(v,loc,f,m)
