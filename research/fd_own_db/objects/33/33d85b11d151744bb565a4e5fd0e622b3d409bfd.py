"""Synthetic original GIFT GloBI full-flora morphology × positive source tests."""
import json
import pandas as pd
import pytest
from fd.fullframe_globi_local_support import (
    genuine_direct_motifs,original_source_guild_morphology_availability,
    strict_source_geography,
)

def _traits():
    return dict(floral_symmetry={
        "Alpha alba":"actinomorphic","Beta rubra":"zygomorphic"},
      floral_form={"Alpha alba":"tubular","Beta rubra":"bell_campanulate"},
      tube_depth_class={"Alpha alba":"deep","Beta rubra":"deep"},
      flower_primary_color={"Alpha alba":"white","Beta rubra":"red_pink"})

def _flora():
    return pd.DataFrame([
      dict(island_id="I1",accepted_species="Alpha alba",
           origin_status="native",endemic_status="endemic"),
      dict(island_id="I1",accepted_species="Beta rubra",
           origin_status="introduced",endemic_status="not_endemic")])

def _links():
    return pd.DataFrame([
      dict(island_id="I1",accepted_species="Alpha alba",origin_status="native",
           endemic_status="endemic",guild="flower_visiting_birds",
           reference_id="paper1",provider_id="archive_A",evidence_class="flower_visit",
           latitude=5.,longitude=10.),
      dict(island_id="I1",accepted_species="Beta rubra",origin_status="introduced",
           endemic_status="not_endemic",guild="lepidoptera",
           reference_id="paper2",provider_id="archive_B",evidence_class="pollination",
           latitude=5.,longitude=10.),
      dict(island_id="I1",accepted_species="Alpha alba",origin_status="native",
           endemic_status="endemic",guild="non_bombus_bees",
           reference_id="paper3",provider_id="archive_B",evidence_class="flower_visit",
           latitude=5.,longitude=10.),
    ])

def test_source_true_pairs_require_same_accepted_species_direct_axes():
    z=genuine_direct_motifs(_traits())
    assert z["actinomorphic_deep"]["Alpha alba"]
    assert not z["actinomorphic_deep"]["Beta rubra"]
    assert z["zygomorphic_deep"]["Beta rubra"]
    assert z["white_tube_bell"]["Alpha alba"]
    assert z["red_pink_tube_bell"]["Beta rubra"]
    # A documented colour and an independently documented shape on
    # DIFFERENT species cannot create a false joint flower phenotype.
    partial=_traits()
    del partial["flower_primary_color"]["Alpha alba"]
    assert "Alpha alba" not in genuine_direct_motifs(partial)["white_tube_bell"]

def test_original_source_omission_zeros_documentation_not_ecological_absence():
    out=original_source_guild_morphology_availability(
        _links(),_flora(),_traits(),"archive_A",{"I1":"tropical"})
    assert len(out["cells"])==2*5*3*5*5
    x=[q for q in out["cells"] if
       q["source_mode"]=="all_cited" and
       q["region"]=="tropical" and
       q["botanical_status"]=="original_native_endemic" and
       q["real_species_direct_pair"]=="actinomorphic_deep" and
       q["original_source_visitor_guild"]=="flower_visiting_birds"][0]
    assert x["n_positive_cited_relation_units"]==1
    assert x["n_documented_reference_labels"]==1
    after=[q for q in out["cells"] if
       q["source_mode"]=="largest_archive_excluded" and
       q["region"]=="tropical" and
       q["botanical_status"]=="original_native_endemic" and
       q["real_species_direct_pair"]=="actinomorphic_deep" and
       q["original_source_visitor_guild"]=="flower_visiting_birds"][0]
    assert after["n_positive_cited_relation_units"]==0
    assert after["no_recorded_positive_is_NOT_ecological_absence"]
    assert out["n_original_source_bibliographic_keys"]==3
    raw=json.dumps(out)
    assert not any('"'+k+'"' in raw for k in (
        "island_id","accepted_species","reference_id",
        "provider_id","latitude","longitude","source_genus"))

def test_positive_geopoint_must_match_accepted_plant_and_verified_polygon():
    gpd=pytest.importorskip("geopandas")
    shp=pytest.importorskip("shapely.geometry")
    poly=gpd.GeoDataFrame(
        dict(island_id=["I1"]),
        geometry=[shp.box(9.9,4.9,10.1,5.1)],crs="EPSG:4326")
    rows=[
      ("Alpha alba","flower_visiting_birds","paper1",
       "archive_A","flower_visit",5.,10.),
      ("Another sp","non_bombus_bees","paper2",
       "archive_B","flower_visit",5.,10.),
      ("Beta rubra","lepidoptera","paper2",
       "archive_B","pollination",5.,10.),
      ("Alpha alba","flower_visiting_birds","paper3",
       "archive_A","flower_visit",50.,10.),
    ]
    linked,audit=strict_source_geography(rows,poly,_flora())
    assert len(linked)==2
    assert audit["n_original_candidate_points"]==4
    assert audit["n_candidate_points_in_exact_island_polygons"]==3
    assert audit["n_original_accepted_island_plant_positive_source_rows"]==2
