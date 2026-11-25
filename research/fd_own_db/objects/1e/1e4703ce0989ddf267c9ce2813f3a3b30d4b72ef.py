"""No invented island pollinators: public network metadata and GIFT joins."""
import numpy as np
import pandas as pd
import pytest

from fd.independent_pollination_networks import (
    original_binomial,visitor_guild,network_location_metadata,
    mangal_location_metadata,identify_actual_plant_edges,
    source_bounded_site_plant_audit,real_same_species_morphologies,
)


def test_strict_taxonomy_never_fills_unidentified_pollinator_or_plant():
    assert original_binomial("Impatiens balsamina L.")=="Impatiens balsamina"
    assert original_binomial("Ipomoea_pes-caprae")=="Ipomoea pes-caprae"
    assert original_binomial("Canavalia sp.") is None
    assert original_binomial("unidentified flower") is None
    assert visitor_guild("Lepidoptera","Sphingidae")=="sphingidae_moth"
    assert visitor_guild("Lepidoptera","Nymphalidae")=="butterfly_family_supported"
    assert visitor_guild("Lepidoptera","Unknown")=="other_lepidoptera"
    assert visitor_guild("Hymenoptera","Apidae")=="bees_family_supported"
    assert visitor_guild("Diptera","Syrphidae")=="flies_order_supported"
    assert visitor_guild("Passeriformes","Nectariniidae")=="bird_family_supported"
    assert visitor_guild("Hymenoptera","")=="unresolved_visitor_taxonomy"


def test_web_of_life_original_site_catalogue_requires_true_coordinates():
    n=pd.DataFrame([
        {"network":"M_PL_044","interaction_type":"Pollination","latitude":28.377248,
         "longitude":129.493741,"reference":"Kato 2000"},
        {"network":"M_PL_045","interaction_type":"Pollination","latitude":"NA",
         "longitude":140.0,"reference":"source"},
        {"network":"M_FR_000","interaction_type":"Herbivory","latitude":30,
         "longitude":120,"reference":"irrelevant"},
    ])
    y,q=network_location_metadata(n)
    assert list(y.network_id)==["M_PL_044"]
    assert q["n_declared_pollination_networks"]==2
    assert q["n_pollination_networks_with_coordinates"]==1
    assert q["n_pollination_networks_without_supported_coordinates"]==1
    with pytest.raises(ValueError,match="conflicting"):
        network_location_metadata(pd.concat([n,n.iloc[:1]]))


def test_mangal_source_point_is_search_candidate_not_proven_visit():
    networks=[
        {"id":123,"dataset_id":4,"geom":{"type":"Point","coordinates":[130.,28.]}},
        {"id":124,"dataset_id":5,"geom":{"type":"Point","coordinates":[13.,51.]}},
        {"id":125,"dataset_id":4,"geom":{"type":"Polygon","coordinates":[]}},
    ]
    dataset=[
        {"id":4,"name":"Bee pollination networks","ref_id":77},
        {"id":5,"name":"Fish trophic networks","ref_id":99},
    ]
    r,a=mangal_location_metadata(networks,dataset)
    assert list(r.network_id)==["123"]
    assert a["n_mangal_keyword_candidate_networks"]==2
    assert a["n_mangal_candidate_georeferenced_points"]==1
    assert a["source_interaction_types_and_nodes_not_yet_verified"]


def test_real_plant_edge_role_not_default_matrix_orientation():
    edges=pd.DataFrame([
        dict(row_species="Impatiens balsamina",column_species="Macroglossum stellatarum",
             row_role="plant",column_role="pollinator"),
        dict(row_species="Bombus terrestris",column_species="Impatiens balsamina",
             row_role="pollinator",column_role="plant"),
        dict(row_species="Other plant sp.",column_species="unknown animal",
             row_role="row",column_role="column"),
    ])
    sp=pd.DataFrame([
        dict(species="Macroglossum stellatarum",order="Lepidoptera",family="Sphingidae"),
        dict(species="Bombus terrestris",order="Hymenoptera",family="Apidae"),
    ])
    out,meta=identify_actual_plant_edges(edges,sp)
    assert len(out)==2
    assert set(out.guild)=={"bees_family_supported","sphingidae_moth"}
    assert meta["n_unresolved_plant_role_or_name"]==1
    with pytest.raises(ValueError,match="edge role"):
        identify_actual_plant_edges(edges.drop(columns="row_role"),sp)


def test_original_same_accepted_species_two_axes_not_cartesian_product():
    traits={
        "floral_symmetry":{
            "Impatiens balsamina":"zygomorphic","Camellia japonica":"actinomorphic",
            "Ipomoea pes-caprae":"actinomorphic"
        },
        "floral_form":{
            "Impatiens balsamina":"tubular","Camellia japonica":"open_radial",
        },
        "tube_depth_class":{
            "Impatiens balsamina":"deep","Camellia japonica":"shallow"
        },
    }
    assert real_same_species_morphologies("Camellia japonica",traits)[
        "actinomorphic_deep"] is False
    assert real_same_species_morphologies("Impatiens balsamina",traits)[
        "zygomorphic_deep"] is True
    assert real_same_species_morphologies("Ipomoea pes-caprae",traits) is None


def test_real_network_plant_guild_audit_preserves_full_original_population():
    traits={
        "floral_symmetry":{"Impatiens balsamina":"zygomorphic",
                           "Camellia japonica":"actinomorphic"},
        "floral_form":{"Impatiens balsamina":"tubular",
                       "Camellia japonica":"open_radial"},
        "tube_depth_class":{"Impatiens balsamina":"deep",
                            "Camellia japonica":"shallow"},
    }
    flora=pd.DataFrame([
        dict(island_id="strict_geo_original",accepted_species="Impatiens balsamina",
             origin_status="native",endemic_status="endemic"),
        dict(island_id="strict_geo_original",accepted_species="Camellia japonica",
             origin_status="unresolved",endemic_status="unresolved"),
    ])
    edges=pd.DataFrame([
        dict(network_id="M_PL_044",island_id="strict_geo_original",
             plant="Impatiens balsamina",guild="bees_family_supported"),
        dict(network_id="M_PL_044",island_id="strict_geo_original",
             plant="Camellia japonica",guild="bird_family_supported"),
        dict(network_id="M_PL_044",island_id="other island",
             plant="Impatiens balsamina",guild="sphingidae_moth"),
    ])
    q=source_bounded_site_plant_audit(edges,flora,traits)
    assert q["n_original_global_GIFT_botanical_islands"]==4505
    assert q["n_original_species_network_plant_guild_units"]==2
    a=next(z for z in q["cells"] if z["botanical_scope"]=="native_label"
           and z["verified_real_joint_morphology"]=="zygomorphic_deep"
           and z["source_network_visitor_taxonomy"]=="bees_family_supported")
    assert a["n_plant_guild_network_documented_units"]==1
    assert a["effective_pollen_transfer_not_measured"]
    assert a["no_network_edge_not_a_biological_absence"]
    b=next(z for z in q["cells"] if z["botanical_scope"]=="native_label"
           and z["verified_real_joint_morphology"]=="actinomorphic_deep"
           and z["source_network_visitor_taxonomy"]=="bird_family_supported")
    assert b["n_plant_guild_network_documented_units"]==0
    endemic=next(z for z in q["cells"] if
        z["botanical_scope"]=="endemic_label" and
        z["verified_real_joint_morphology"]=="zygomorphic_deep" and
        z["source_network_visitor_taxonomy"]=="bees_family_supported")
    assert endemic["n_plant_guild_network_documented_units"]==1


def test_no_flora_genus_imputation_or_fabricated_island_link():
    flora=pd.DataFrame([dict(island_id="island1",accepted_species="Impatiens balsamina",
                             origin_status="unresolved",endemic_status="unresolved")])
    with pytest.raises(ValueError,match="frame reduced"):
        source_bounded_site_plant_audit(
            pd.DataFrame(columns=["network_id","island_id","plant","guild"]),
            flora,{},global_original_islands=300)
    with pytest.raises(ValueError,match="duplicated"):
        source_bounded_site_plant_audit(
            pd.DataFrame(columns=["network_id","island_id","plant","guild"]),
            pd.concat([flora,flora]),{},global_original_islands=4505)
