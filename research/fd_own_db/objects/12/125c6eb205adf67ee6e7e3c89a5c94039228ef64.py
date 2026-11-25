"""Synthetic guards: a flower count for another plant is NOT a denominator."""
import copy
import pandas as pd
import pytest

from fd.euppollnet_species_identity_gate import (
    qualified_rows_and_identity_receipt,
    original_same_plant_exposure_audit,
)


def _frames():
    v=pd.DataFrame([
      dict(Study_id="S1",Network_id="N1",Date="2020-06-18",
           Flower_data_merger="K1",island_id="GIFT1",plant="Alpha alba",
           origin_status="native",endemic_status="not_endemic",
           Pollinator_order="Hymenoptera",Pollinator_family="Apidae",
           Pollinator_accepted_name="Bombus terrestris"),
      dict(Study_id="S1",Network_id="N1",Date="2020-06-18",
           Flower_data_merger="K2",island_id="GIFT1",plant="Beta rubra",
           origin_status="native",endemic_status="not_endemic",
           Pollinator_order="Lepidoptera",Pollinator_family="Papilionidae",
           Pollinator_accepted_name="Papilio machaon"),
      dict(Study_id="S1",Network_id="N1",Date="2020-06-18",
           Flower_data_merger="K3",island_id="GIFT1",plant="Alpha alba",
           origin_status="native",endemic_status="not_endemic",
           Pollinator_order="Hymenoptera",Pollinator_family="Apidae",
           Pollinator_accepted_name="Apis mellifera"),
      dict(Study_id="S1",Network_id="N1",Date="2020-06-18",
           Flower_data_merger="K4",island_id="GIFT1",plant="Alpha alba",
           origin_status="native",endemic_status="not_endemic",
           Pollinator_order="Hymenoptera",Pollinator_family="Apidae",
           Pollinator_accepted_name="Apis mellifera"),
      dict(Study_id="S1",Network_id="N1",Date="2020-06-18",
           Flower_data_merger="K5",island_id="GIFT1",plant="Beta rubra",
           origin_status="native",endemic_status="not_endemic",
           Pollinator_order="Lepidoptera",Pollinator_family="Papilionidae",
           Pollinator_accepted_name="Papilio machaon"),
      dict(Study_id="S1",Network_id="N1",Date="2020-06-18",
           Flower_data_merger="K6",island_id="GIFT1",plant="Alpha alba",
           origin_status="native",endemic_status="not_endemic",
           Pollinator_order="Hymenoptera",Pollinator_family="Apidae",
           Pollinator_accepted_name="Bombus terrestris"),
      dict(Study_id="S1",Network_id="N1",Date="2020-06-18",
           Flower_data_merger="K7",island_id="GIFT1",plant="Alpha alba",
           origin_status="native",endemic_status="not_endemic",
           Pollinator_order="Hymenoptera",Pollinator_family="Apidae",
           Pollinator_accepted_name="Bombus terrestris"),
      dict(Study_id="S1",Network_id="N1",Date="2020-06-18",
           Flower_data_merger="K8",island_id="GIFT1",plant="Alpha alba",
           origin_status="native",endemic_status="not_endemic",
           Pollinator_order="Hymenoptera",Pollinator_family="Apidae",
           Pollinator_accepted_name="Bombus terrestris"),
    ])
    f=pd.DataFrame([
      dict(Flower_data_merger="K1",Study_id="S1",Network_id="N1",
           Day="18",Month="6",Year="2020",Plant_species="Alpha alba",
           Flower_count="12",Units="Flowers"),
      dict(Flower_data_merger="K2",Study_id="S1",Network_id="N1",
           Day="18",Month="6",Year="2020",Plant_species="Beta rubra",
           Flower_count="20",Units="Flowers"),
      dict(Flower_data_merger="K3",Study_id="S1",Network_id="N1",
           Day="18",Month="6",Year="2020",Plant_species="Another plant",
           Flower_count="40",Units="Flowers"),
      dict(Flower_data_merger="K4",Study_id="S1",Network_id="N1",
           Day="18",Month="6",Year="2020",Plant_species="Alpha alba var. minor",
           Flower_count="10",Units="Flowers"),
      dict(Flower_data_merger="K5",Study_id="S1",Network_id="N1",
           Day="19",Month="6",Year="2020",Plant_species="Beta rubra",
           Flower_count="17",Units="Flowers"),
      dict(Flower_data_merger="K6",Study_id="S1",Network_id="N1",
           Day="18",Month="6",Year="2020",Plant_species="Alpha alba",
           Flower_count="0",Units="Flowers"),
      dict(Flower_data_merger="K7",Study_id="S1",Network_id="N1",
           Day="18",Month="6",Year="2020",Plant_species="Alpha alba",
           Flower_count="8",Units="Flowers"),
      dict(Flower_data_merger="K7",Study_id="S1",Network_id="N1",
           Day="18",Month="6",Year="2020",Plant_species="Alpha alba",
           Flower_count="8",Units="Flowers"),
      dict(Flower_data_merger="K8",Study_id="S2",Network_id="N1",
           Day="18",Month="6",Year="2020",Plant_species="Alpha alba",
           Flower_count="8",Units="Flowers"),
    ])
    traits=dict(
      floral_symmetry={"Alpha alba":"actinomorphic","Beta rubra":"zygomorphic"},
      tube_depth_class={"Alpha alba":"deep","Beta rubra":"deep"},
      floral_form={"Alpha alba":"tubular","Beta rubra":"bell_campanulate"},
      flower_primary_color={"Alpha alba":"white","Beta rubra":"red_pink"},
    )
    return v,f,traits


def test_same_network_day_unique_key_but_different_flower_species_cannot_be_exposure():
    visits,flowers,_=_frames()
    linked,r=qualified_rows_and_identity_receipt(visits,flowers)
    assert r["n_original_same_island_positive_visitor_rows"]==8
    assert r["n_original_flower_count_keys_repeated"]==1
    assert r["n_positive_rows_matching_unique_original_flower_key"]==7
    assert r["n_positive_rows_same_published_study_and_network"]==6
    assert r["n_positive_rows_same_published_study_network_calendar_day"]==5
    assert r["n_positive_rows_with_source_nonzero_flower_count_and_unit"]==4
    assert r["n_positive_rows_with_unresolvable_original_flower_plant"]==1
    assert r["n_positive_rows_with_different_source_flower_plant"]==1
    assert r["n_positive_rows_with_exact_same_flower_and_visitor_plant"]==2
    assert set(linked.plant)=={"Alpha alba","Beta rubra"}
    assert r["n_distinct_exact_flower_survey_keys"]==2


def test_two_traits_and_two_guilds_only_in_actually_shared_flower_supply_context():
    visits,flowers,traits=_frames()
    a=original_same_plant_exposure_audit(visits,flowers,traits)
    q=a["source_same_study_site_date_and_unit_comparability"]
    assert q["n_same_source_study_network_date_island_unit_contexts"]==1
    assert q["n_contexts_with_two_botanical_flower_species"]==1
    assert q["n_contexts_with_two_documented_visitor_guilds"]==1
    assert q["n_contexts_with_two_botanical_genera"]==1
    x=next(x for x in q["source_motif_comparison_opportunities"]
        if x["real_same_species_direct_motif"]=="actinomorphic_deep")
    assert x["n_source_contexts_with_true_and_false_scored_flower_species"]==1
    assert x["n_source_contexts_with_flower_states_two_guilds_and_two_genera"]==1
    assert x["n_source_studies_with_comparable_flower_states"]==1
    assert x["n_original_GIFT_islands_with_comparable_flower_states"]==1
    assert a["no_visitor_preference_effect_or_pollen_fitness_coefficient"]


def test_original_4505_parent_must_be_preserved_and_missing_plant_col_rejected():
    visits,flowers,traits=_frames()
    with pytest.raises(ValueError,match="parent floras changed"):
        original_same_plant_exposure_audit(visits,flowers,traits,original_parent_n=927)
    with pytest.raises(ValueError,match="fields missing"):
        original_same_plant_exposure_audit(visits,flowers.drop(columns=["Plant_species"]),traits)
