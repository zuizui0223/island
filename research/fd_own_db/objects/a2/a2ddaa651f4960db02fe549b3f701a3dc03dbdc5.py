"""Synthetic source-rights and original-network geography checks only."""
import json
import pandas as pd
import pytest

from fd.euppollnet_island_overlap import (
    SOURCE_COLUMNS, validate_original_visits, audit_original_island_euppollnet
)

def synthetic():
    records=[
        ("s1","n1","2016-06-05",35.0,-5.0,"Plant A","SPECIES",
         "Bombus example","Hymenoptera","Bombus","Apidae"),
        ("s1","n1","2018-06-06",35.0,-5.0,"Plant A","SPECIES",
         "Syrphus example","Diptera","Syrphus","Syrphidae"),
        ("s2","n2","2019-06-06",34.0,-5.0,"Plant B","SPECIES",
         "Apis example","Hymenoptera","Apis","Apidae"),
    ]
    d=pd.DataFrame(records,columns=SOURCE_COLUMNS)
    sites=pd.DataFrame([
        ("s1","n1",35.0,-5.0,"isla"),
        ("s2","n2",34.0,-5.0,"islb")
    ],columns=["Study_id","Network_id","Latitude","Longitude","island_id"])
    flora=pd.DataFrame([
        ("isla","Plant A","native"),("islb","Plant C","native")
    ],columns=["island_id","accepted_species","origin_status"])
    geo=pd.DataFrame([
        ("isla","tropical","b1"),("islb","tropical","b2")
    ],columns=["island_id","analysis_regime","spatial_block"])
    return d,sites,flora,geo


def test_actual_guilds_native_exact_species_and_repeated_years():
    data,sites,flora,geo=synthetic()
    z=audit_original_island_euppollnet(data,sites,flora,geo)
    assert z["n_original_928_flora_islands_with_exact_network_site"]==2
    assert z["n_direct_species_joined_interaction_rows"]==2
    assert z["n_strict_native_direct_species_interactions"]==2
    assert z["n_study_network_coordinate_sites_with_two_years"]==1
    assert z["n_regions_meeting_10_islands_5blocks_3studies"]==0
    assert z["regions"][0]["n_documented_bombus_interaction_rows"]==1
    assert z["regions"][0]["n_documented_diptera_interaction_rows"]==1
    assert z["source_island_pollinator_shift_or_loss_identified"] is False
    exported=json.dumps(z)
    for forbidden in ('"island_id"', '"Study_id"', '"Plant_accepted_name"',
                      '"Latitude"','"Longitude"','"Pollinator_accepted_name"'):
        assert forbidden not in exported


def test_no_site_overlap_is_explicit_stop_not_absence():
    d,sites,flora,geo=synthetic()
    z=audit_original_island_euppollnet(
      d,sites.assign(island_id=["different","elsewhere"]),flora,geo)
    assert z["status"]=="NO_SOURCE_VERIFIED_FD_ISLAND_NETWORK_INTERACTIONS"
    assert z["source_island_pollinator_shift_or_loss_identified"] is False


def test_original_data_schema_and_identification_fail_closed():
    d,sites,flora,geo=synthetic()
    with pytest.raises(ValueError, match="header"):
        validate_original_visits(d.drop(columns=["Date"]))
    with pytest.raises(ValueError, match="site"):
        audit_original_island_euppollnet(
            d,pd.concat([sites,sites.iloc[[0]]],ignore_index=True),flora,geo)
