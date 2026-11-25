"""Synthetic counterfactual network guild-removal, no biological loss claims."""
import json
import pandas as pd
from fd.euppollnet_island_overlap import SOURCE_COLUMNS, validate_original_visits
from fd.euppollnet_guild_loss_simulation import (
    _group_remove, original_island_network_guild_ablation
)

def _data():
    rows=[]
    for i in range(13):
        plant=f"Test Plant {i}"
        # Bee link to 10 plants; sampled Diptera also visit 11 plants,
        # of which 8 overlap; source sampling scope supports both.
        if i<10:
            rows.append(("study","network","2019-07-01",42.,12.,plant,"SPECIES",
                  "Bombus demo","Hymenoptera","Bombus","Apidae"))
        if i<8 or i>=10:
            rows.append(("study","network","2019-07-01",42.,12.,plant,"SPECIES",
                  "Syrphus demo","Diptera","Syrphus","Syrphidae"))
        if i<6:
            rows.append(("study","network","2019-07-01",42.,12.,plant,"SPECIES",
                  "Andrena demo","Hymenoptera","Andrena","Andrenidae"))
    d=pd.DataFrame(rows,columns=SOURCE_COLUMNS)
    s=pd.DataFrame([("study","network",42.,12.,"isla")],columns=[
        "Study_id","Network_id","Latitude","Longitude","island_id"])
    flora=pd.DataFrame([("isla",f"Test Plant {i}","native") for i in range(13)],
       columns=["island_id","accepted_species","origin_status"])
    direct=pd.DataFrame([
        (f"Test Plant {i}","white" if i%2 else "red",
         "open_radial" if i%3 else "tubular",
         "actinomorphic" if i%3 else "zygomorphic")
        for i in range(13)
    ],columns=["accepted_species","flower_primary_color","floral_form","floral_symmetry"]).set_index("accepted_species")
    geo=pd.DataFrame([("isla","northern_midlatitude","b1")],
        columns=["island_id","analysis_regime","spatial_block"])
    meta=pd.DataFrame([("study",2019,"Hymenoptera, Diptera, Lepidoptera and Coleoptera")],
        columns=["Study_id","Year","Taxa_recorded"])
    return d,s,flora,direct,geo,meta


def test_simulated_loss_measures_remaining_plant_visits_not_seedset():
    d,s,f,t,g,m=_data()
    z=original_island_network_guild_ablation(d,s,f,t,g,m)
    assert z["n_studies_with_source_documented_bee_and_nonbee_sampling_scope"]==1
    assert z["observed_pollinator_loss_shift_or_effective_pollen_delivery"] is False
    assert z["scenario_aggregates"]
    bee=[x for x in z["scenario_aggregates"] if x["scenario"]=="bombus"][0]
    assert bee["n_original_islands"]==1
    assert bee["n_spatial_blocks"]==1
    assert bee["n_networks_with_8_matched_Direct_morph_plants_after_loss"]==1
    assert bee["functional_diversity_claim_supported"] is False
    assert bee["n_networks_with_equal_removed_visitor_species_null"]==1
    assert 0<bee["median_alternative_visitor_link_coverage_among_exposed_plants"]<1
    out=json.dumps(z)
    for x in ('"island_id"','"Study_id"','"Network_id"',
              '"Plant_accepted_name"','"Pollinator_accepted_name"'):
        assert x not in out


def test_bee_only_original_study_does_not_create_false_compensation_zero():
    d,s,f,t,g,m=_data()
    m.loc[0,"Taxa_recorded"]="Bumblebees"
    z=original_island_network_guild_ablation(d,s,f,t,g,m)
    assert z["scenario_aggregates"]==[]
    assert z["n_studies_with_source_documented_bee_and_nonbee_sampling_scope"]==0
