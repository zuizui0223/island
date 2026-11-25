"""Synthetic original GIFT 4,505 botanical accepted floral community × European field-site audit."""
import numpy as np
import pandas as pd
import pytest
# The default fd CI installs only the core numerical package. Geospatial
# source validation is optional and fully tested in the dedicated source job.
gpd=pytest.importorskip("geopandas")
box=pytest.importorskip("shapely.geometry").box
from fd.euppollnet_gift_exact_site_bridge import (
    original_guild,strict_original_accepted_plant,
    exact_source_network_site_union,evaluate_original_european_island_visits,
)


def _frame():
    ids=[f"I{i:04d}" for i in range(4505)]
    flora=pd.DataFrame([dict(
        island_id=i,accepted_species="Dummy alba",
        origin_status="unresolved",endemic_status="not_endemic")
        for i in ids])
    flora=pd.concat([flora,pd.DataFrame([
        dict(island_id="I0000",accepted_species="Alpha alba",
             origin_status="native",endemic_status="endemic"),
        dict(island_id="I0001",accepted_species="Beta rubra",
             origin_status="native",endemic_status="not_endemic"),
    ])],ignore_index=True)
    polygons=gpd.GeoDataFrame(dict(island_id=["I0000","I0001"]),
        geometry=[box(9.9,4.9,10.1,5.1),box(10.9,5.9,11.1,6.1)],
        crs="EPSG:4326")
    source=pd.DataFrame([
        dict(Study_id="S1",Network_id="N1",Latitude="5.0",
             Longitude="10.0",Date="2020-06-18",
             Plant_accepted_name="Alpha alba",Pollinator_accepted_name="Bombus terrestris",
             Pollinator_order="Hymenoptera",Pollinator_family="Apidae"),
        dict(Study_id="S1",Network_id="N1",Latitude="5.0",
             Longitude="10.0",Date="2020-06-18",
             Plant_accepted_name="Alpha alba",Pollinator_accepted_name="Papilio machaon",
             Pollinator_order="Lepidoptera",Pollinator_family="Papilionidae"),
        dict(Study_id="S2",Network_id="N2",Latitude="6.0",
             Longitude="11.0",Date="2021-07-03",
             Plant_accepted_name="Beta rubra",Pollinator_accepted_name="Macroglossum stellatarum",
             Pollinator_order="Lepidoptera",Pollinator_family="Sphingidae"),
        dict(Study_id="S2",Network_id="N2",Latitude="6.0",
             Longitude="11.0",Date="2021-07-03",
             Plant_accepted_name="Another species",Pollinator_accepted_name="Papilio machaon",
             Pollinator_order="Lepidoptera",Pollinator_family="Papilionidae"),
        dict(Study_id="S3",Network_id="N3",Latitude="45.0",
             Longitude="10.0",Date="2020-06-18",
             Plant_accepted_name="Alpha alba",Pollinator_accepted_name="Bombus terrestris",
             Pollinator_order="Hymenoptera",Pollinator_family="Apidae"),
    ])
    traits=dict(floral_symmetry={"Alpha alba":"actinomorphic",
            "Beta rubra":"zygomorphic"},
        tube_depth_class={"Alpha alba":"deep","Beta rubra":"deep"},
        floral_form={"Alpha alba":"tubular","Beta rubra":"bell_campanulate"},
        flower_primary_color={"Alpha alba":"white","Beta rubra":"red_pink"})
    regime={"I0000":"northern_midlatitude","I0001":"northern_midlatitude"}
    return source,polygons,flora,traits,regime


def test_animal_taxonomy_requires_actual_family():
    assert original_guild("Hymenoptera","Apidae","Bombus terrestris")=="bombus"
    assert original_guild("Hymenoptera","Apidae","Apis mellifera")=="other_bees"
    assert original_guild("Lepidoptera","Sphingidae","Macroglossum stellatarum")=="sphingidae_moth"
    assert original_guild("Lepidoptera","Papilionidae","Papilio machaon")=="butterfly_family"
    assert original_guild("Lepidoptera","","Unresolved")=="other_lepidoptera"
    assert strict_original_accepted_plant("Alpha alba L.")=="Alpha alba"
    assert strict_original_accepted_plant("Alpha alba var. alba") is None
    assert strict_original_accepted_plant("Alpha sp.") is None
    assert strict_original_accepted_plant("Alpha cf. alba") is None


def test_original_site_and_study_id_alignment_with_csv_text_lat_lon():
    source,polygons,_,_,_= _frame()
    site,stats=exact_source_network_site_union(source,polygons)
    assert stats["n_source_original_study_ids_with_exact_GIFT_polygon"]==2
    assert stats["n_original_source_study_network_site_units_on_verified_GIFT_polygon"]==2
    assert stats["n_original_distinct_source_coordinate_pairs"]==3
    assert stats["n_original_source_coordinate_pairs_inside_unique_GIFT_polygons"]==2
    assert len(site)==2


def test_exact_gift_4505_source_accepted_plant_and_animal_records():
    source,polygons,flora,traits,regime=_frame()
    q=evaluate_original_european_island_visits(source,polygons,flora,traits,regime)
    assert q["status"]=="SOURCE_EUPPOLLNET_EUROPEAN_FIELD_SITE_EXACT_GIFT_FLORA_OPPORTUNITY"
    assert q["n_original_global_GIFT_botanical_floras"]==4505
    assert q["n_original_source_positive_recorded_rows_matching_same_island_GIFT_plant"]==3
    assert q["n_original_source_distinct_islands_with_matched_flora_plant"]==2
    assert q["n_original_source_study_ids_with_matched_flora_plant"]==2
    assert q["n_original_unique_study_network_site_year_plant_guild_positive_units"]==3
    assert len(q["true_source_direct_flower_and_animal_guild_documentation"])==6*3*5*8
    cells=q["true_source_direct_flower_and_animal_guild_documentation"]
    radial_butterfly=next(x for x in cells if
       x["region"]=="northern_midlatitude" and
       x["botanical_status"]=="original_native" and
       x["same_accepted_species_direct_flower_pair"]=="actinomorphic_deep" and
       x["original_european_source_guild"]=="butterfly_family")
    assert radial_butterfly["n_source_study_network_site_year_positive_units"]==1
    assert radial_butterfly["n_source_independently_verified_GIFT_islands"]==1
    hawkmoth=next(x for x in cells if
       x["region"]=="GLOBAL_ALL" and x["botanical_status"]=="all_recorded" and
       x["same_accepted_species_direct_flower_pair"]=="zygomorphic_deep" and
       x["original_european_source_guild"]=="sphingidae_moth")
    assert hawkmoth["n_source_studies_with_positive_motif_visits"]==1
    assert q["no_source_interaction_edge_implicit_negative"]
    assert q["source_original_study_overlap_with_GloBI_WOL_DoPI_not_checked"]


def test_unverified_site_is_not_a_negative_visit_and_origin_mask_is_not_inferred():
    source,polygons,flora,traits,regime=_frame()
    source=source[source.Study_id.eq("S3")].copy()
    q=evaluate_original_european_island_visits(source,polygons,flora,traits,regime)
    assert q["status"]=="SOURCE_FIELD_NETWORK_ISLAND_OPPORTUNITY_UNIDENTIFIABLE"
    assert q["n_original_same_island_accepted_flower_visitor_units"]==0
    assert q["positive_source_absence_not_biological_no_visitation"]
