"""Verified GBIF family/order from original Mangal source taxon key only."""
import pandas as pd
import pytest

from fd.gbif_mangal_visitor_taxonomy import (
    gbif_taxonomy_visitor_guild,
    target_original_morphology_partner_keys,
    summarize_gbif_resolved_source,
)


def original_source():
    traits={
        "floral_symmetry":{
            "Ipomoea batatas":"actinomorphic",
            "Impatiens balsamina":"zygomorphic",
        },
        "floral_form":{
            "Ipomoea batatas":"tubular",
            "Impatiens balsamina":"tubular",
        },
        "tube_depth_class":{
            "Ipomoea batatas":"deep",
            "Impatiens balsamina":"shallow",
        },
    }
    lines=pd.DataFrame([
        dict(network_id="network1",island_id="island1",
             plant="Ipomoea batatas",original_partner_gbif_key=123,
             original_partner_scientific_name="Hyles livornica",
             interaction_id="int1"),
        dict(network_id="network1",island_id="island1",
             plant="Ipomoea batatas",original_partner_gbif_key=124,
             original_partner_scientific_name="Apis mellifera",
             interaction_id="int2"),
        dict(network_id="network1",island_id="island1",
             plant="Impatiens balsamina",original_partner_gbif_key=125,
             original_partner_scientific_name="Hyles gallii",
             interaction_id="int3"),
        dict(network_id="network2",island_id="island2",
             plant="Ipomoea batatas",original_partner_gbif_key=None,
             original_partner_scientific_name="unidentified fly",
             interaction_id="int4"),
    ])
    return lines,traits


def test_verified_original_gbif_taxon_key_and_animal_kingdom_required():
    assert gbif_taxonomy_visitor_guild({
        "key":123,"kingdom":"Animalia","order":"Lepidoptera",
        "family":"Sphingidae",
    },123)=="sphingidae_moth"
    assert gbif_taxonomy_visitor_guild({
        "key":123,"kingdom":"Animalia","order":"Lepidoptera",
        "family":"Nymphalidae",
    },124)=="unresolved_visitor_taxonomy"
    assert gbif_taxonomy_visitor_guild({
        "key":123,"kingdom":"Plantae","order":"Lepidoptera",
        "family":"Sphingidae",
    },123)=="unresolved_visitor_taxonomy"
    assert gbif_taxonomy_visitor_guild({},123)=="unresolved_visitor_taxonomy"


def test_only_original_true_same_species_both_direct_morphological_axes():
    src,trait=original_source()
    v=target_original_morphology_partner_keys(src,trait)
    assert len(v)==3
    assert set(v.plant)=={"Ipomoea batatas"}
    assert v.actinomorphic_deep_true.eq(True).all()
    assert v.zygomorphic_deep_true.eq(False).all()
    assert v.partner_gbif_key.dropna().nunique()==2


def test_unresolved_gbif_key_does_not_erase_original_positive_visit_record():
    src,trait=original_source()
    original=target_original_morphology_partner_keys(src,trait)
    original["origin_status"]="native"
    original["endemic_status"]="nonendemic"
    gbif={
        123:dict(key=123,kingdom="Animalia",
                 family="Sphingidae",order="Lepidoptera"),
        124:dict(key=124,kingdom="Animalia",
                 family="Apidae",order="Hymenoptera"),
    }
    result=summarize_gbif_resolved_source(original,gbif)
    assert result["n_original_GIFT_botanical_islands"]==4505
    assert result["n_original_direct_flower_pair_partner_source_interactions"]==3
    assert result["n_GBIF_keys_retrieved_from_public_backbone"]==2
    assert result["n_source_guild_support_cells"]==4*3*7
    a=next(z for z in result["support"] if
           z["same_accepted_plant_true_flower_pair"]=="actinomorphic_deep"
           and z["GBIF_verified_partner_guild"]=="sphingidae_moth")
    assert a["n_source_positive_network_island_plant_partner_units"]==1
    x=next(z for z in result["support"] if
           z["same_accepted_plant_true_flower_pair"]=="actinomorphic_deep"
           and z["GBIF_verified_partner_guild"]=="unresolved_visitor_taxonomy")
    assert x["n_source_positive_network_island_plant_partner_units"]==1


def test_no_fabricated_full_GIFT_study_population():
    src,trait=original_source()
    r=target_original_morphology_partner_keys(src,trait)
    r["origin_status"]="native"
    r["endemic_status"]="nonendemic"
    with pytest.raises(ValueError,match="full botanical"):
        summarize_gbif_resolved_source(r,{},total_original_islands=40)
