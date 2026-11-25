"""Synthetic source full flower choice set never manufactures animal absences."""
import pandas as pd
import pytest

from fd.euppollnet_choice_set_audit import audited_nonvisited_flower_census_opportunity


def _data():
    flora=pd.DataFrame([
        dict(island_id=f"I{k}",accepted_species="Dummy alba")
        for k in range(4505)])
    flora=pd.concat([flora,pd.DataFrame([
        dict(island_id="I0",accepted_species="Alpha alba"),
        dict(island_id="I0",accepted_species="Beta rubra"),
    ])],ignore_index=True)
    visit=pd.DataFrame([
        dict(Study_id="S1",Network_id="N1",Site_id="A",
             Date="2020-06-18",Units="Flowers",island_id="I0",
             plant="Alpha alba",Flower_data_merger="K1"),
    ])
    source=pd.DataFrame([
        dict(Study_id="S1",Network_id="N1",
             Latitude="5",Longitude="10"),
        dict(Study_id="S1",Network_id="N1",
             Latitude="5",Longitude="10"),
    ])
    flower=pd.DataFrame([
        dict(Study_id="S1",Network_id="N1",Site_id="A",
             Day="18",Month="6",Year="2020",Units="Flowers",
             Plant_species="Alpha alba",Flower_count="12",
             Flower_data_merger="K1"),
        dict(Study_id="S1",Network_id="N1",Site_id="A",
             Day="18",Month="6",Year="2020",Units="Flowers",
             Plant_species="Beta rubra",Flower_count="20",
             Flower_data_merger="K2"),
        dict(Study_id="S1",Network_id="N1",Site_id="A",
             Day="18",Month="6",Year="2020",Units="Flowers",
             Plant_species="Other alba",Flower_count="22",
             Flower_data_merger="K3"),
        dict(Study_id="S1",Network_id="N1",Site_id="B",
             Day="18",Month="6",Year="2020",Units="Flowers",
             Plant_species="Gamma alba",Flower_count="1000",
             Flower_data_merger="K4"),
        dict(Study_id="S1",Network_id="N1",Site_id="A",
             Day="18",Month="6",Year="2020",Units="Flower heads",
             Plant_species="Alpha alba",Flower_count="1000",
             Flower_data_merger="K5"),
    ])
    traits=dict(
        floral_symmetry={"Alpha alba":"actinomorphic",
                          "Beta rubra":"zygomorphic"},
        tube_depth_class={"Alpha alba":"deep","Beta rubra":"deep"},
        floral_form={"Alpha alba":"tubular",
                    "Beta rubra":"bell_campanulate"},
        flower_primary_color={"Alpha alba":"white",
                             "Beta rubra":"red_pink"},
    )
    return visit,flower,source,flora,traits


def test_real_census_adds_all_accepted_flowers_not_visited_as_biological_absences():
    v,f,source,flora,traits=_data()
    q=audited_nonvisited_flower_census_opportunity(v,f,source,flora,traits)
    assert q["status"]=="GIFT_4505_EUPPOLLNET_FULL_SITE_DAY_FLOWER_CENSUS_CHOICESET_PREFLIGHT"
    a=q["source_waterfall"]
    assert a["n_original_distinct_positive_anchor_site_day_unit_contexts"]==1
    assert a["n_original_all_source_flower_rows_on_anchored_study_site_day_unit_before_GIFT_check"]==3
    assert a["n_source_census_flower_rows_with_original_same_island_accepted_species"]==2
    assert a["n_source_census_original_flower_species_contexts_without_positive_reference"]==1
    assert a["n_original_source_studies_with_flower_census_on_GIFT_island"]==1
    assert a["n_original_GIFT_islands_with_flower_census"]==1
    x=next(x for x in q[
        "true_same_species_botanical_morphology_source_census_opportunities"]
        if x["true_same_species_direct_motif"]=="actinomorphic_deep")
    assert x["n_original_site_day_unit_source_flora_contexts_with_both_scored_morphology_states"]==1
    assert q["limits"]["surveyed_flower_without_recorded_animal_link_not_biological_nondetection"]


def test_reused_study_network_at_another_GPS_is_not_a_single_island_choice_set():
    v,f,source,flora,traits=_data()
    source=pd.concat([source,pd.DataFrame([dict(
        Study_id="S1",Network_id="N1",
        Latitude="5.5",Longitude="10")])],ignore_index=True)
    q=audited_nonvisited_flower_census_opportunity(v,f,source,flora,traits)
    assert q["status"]=="GIFT_4505_EUPPOLLNET_COMPLETE_SITE_FLOWER_CHOICESET_UNIDENTIFIABLE"
    a=q["source_waterfall"]
    assert a["n_source_study_network_keys_with_multiple_original_GPS"]==1
    assert a["n_positive_rows_on_single_GPS_source_study_network"]==0
    assert q["absent_source_choice_set_is_NOT_no_animals"]


def test_unresolved_source_site_plant_and_study_coverage_rejected():
    v,f,source,flora,traits=_data()
    with pytest.raises(ValueError,match="4505"):
        audited_nonvisited_flower_census_opportunity(
            v,f,source,flora,traits,parent_botanical_floras=927)
    with pytest.raises(ValueError,match="schema missing"):
        audited_nonvisited_flower_census_opportunity(
            v,f.drop(columns=["Plant_species"]),source,flora,traits)
