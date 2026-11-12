import json

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import box

from island_v2.fd_globi_species_shift_spatial_v09 import (
    PIN, file_sha256, original_candidate_island_locality_gate,
)


def fixture(tmp_path):
    islands=[f"ISL{k}" for k in range(13)]
    names=[f"Plant{k} species{k}" for k in range(54)]
    # Exactly the same original cohort dimensions as the 54-species gate.
    pair=[]
    for i,island in enumerate(islands):
        for j in range(15):
            pair.append((island,names[(i*15+j)%54]))
    pair.extend([(islands[0],names[25]),(islands[0],names[26])])
    pairs=pd.DataFrame(pair,columns=["island_id","accepted_species"])
    assert len(pairs)==197 and not pairs.duplicated().any()
    assert pairs.accepted_species.nunique()==54 and pairs.island_id.nunique()==13
    geometry=gpd.GeoDataFrame({"island_id":islands},
        geometry=[box(i*3,0,i*3+1,1) for i in range(13)],crs="EPSG:4326")
    p=tmp_path/"original_islands.gpkg"
    geometry.to_file(p,layer="islands",driver="GPKG")
    receipt={
        "status":"FOCAL_21_POLYGONS_IDENTICAL_BUT_FULL_UNIVERSE_DRIFT",
        "GSHHG_archive_SHA256":PIN,
        "n_original_focal_islands":21,
        "n_focal_islands_with_exact_historical_polygon":21,
        "geometry_file_SHA256":file_sha256(p),
        "global_full_universe_identity_verified":False,
    }
    j=tmp_path/"receipt.json"
    j.write_text(json.dumps(receipt))
    return p,j,pairs


def test_same_native_plant_same_original_island_point_required(tmp_path):
    gpkg,receipt,pairs=fixture(tmp_path)
    plant=pairs[pairs.island_id.eq("ISL0")].accepted_species.iloc[0]
    visitor=[(plant,"bombus","refA","provA","flower_visit",.5,.5),
             (plant,"lepidoptera","refB","provB","flower_visit",.5,.5),
             (plant,"bombus","refC","provC","pollination",.5,.5)]
    result=original_candidate_island_locality_gate(
        set(visitor),{"refA":{"provA"},"refB":{"provB"},"refC":{"provC"}},
        pairs,gpkg,receipt)
    assert result["n_source_georeferenced_rows_inside_any_original_13_island_polygon"]==3
    assert result["n_source_georeferenced_rows_matching_native_species_on_exact_same_island"]==3
    assert result["n_distinct_focal_islands_with_record_matching_native_species"]==1
    assert result["n_same_island_plant_species_with_2plus_documented_guilds"]==1
    assert result["n_original_explicit_pollination_relations_inside_island_with_native_species"]==1
    assert result["n_confirmed_same_species_mainland_to_island_pollinator_shifts"]==0
    assert not {"plant_species","accepted_species","island_id","reference_id","provider_id"} & set(result)


def test_global_species_match_elsewhere_cannot_make_local_plant_visit(tmp_path):
    gpkg,receipt,pairs=fixture(tmp_path)
    plant=pairs[pairs.island_id.eq("ISL0")].accepted_species.iloc[0]
    p=pairs.loc[~(pairs.island_id.eq("ISL1") & pairs.accepted_species.eq(plant))]
    # If plant occurs on I0 but no corresponding record on I1,
    # a point on I1 is not evidence of that plant's island-local guild.
    z=original_candidate_island_locality_gate(
        {(plant,"bombus","refA","provA","flower_visit",.5,3.5)},
        {"refA":{"provA"}},p,gpkg,receipt)
    assert z["n_source_georeferenced_rows_inside_any_original_13_island_polygon"]==1
    assert z["n_source_georeferenced_rows_matching_native_species_on_exact_same_island"]==0


def test_ambiguous_or_mismatched_original_geometry_fails_closed(tmp_path):
    gpkg,receipt,pairs=fixture(tmp_path)
    j=json.loads(receipt.read_text())
    j["GSHHG_archive_SHA256"]="bad"
    receipt.write_text(json.dumps(j))
    with pytest.raises(ValueError,match="source-hash verification failed"):
        original_candidate_island_locality_gate(
            set(),{},pairs,gpkg,receipt)


def test_shared_reference_across_provider_labels_is_not_spatial_evidence(tmp_path):
    gpkg,receipt,pairs=fixture(tmp_path)
    plant=pairs[pairs.island_id.eq("ISL0")].accepted_species.iloc[0]
    result=original_candidate_island_locality_gate(
        {(plant,"lepidoptera","shared","providerA","flower_visit",.5,.5)},
        {"shared":{"providerA","providerB"}},pairs,gpkg,receipt)
    assert result["n_source_georeferenced_rows_matching_native_species_on_exact_same_island"]==0
