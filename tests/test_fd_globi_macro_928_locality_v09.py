import json

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import box

from island_v2.fd_globi_macro_928_locality_v09 import (
    strict_direct_928_flora, macro_exact_928_locality_gate,
)
from island_v2.fd_globi_species_shift_spatial_v09 import file_sha256, PIN


def _data(tmp_path):
    ledger=[]
    for j in range(30):
        s=f"Plant{j%3} species{j}"
        for trait,val in (
            ("flower_primary_color","white"),("floral_form","open_radial"),
            ("floral_symmetry","actinomorphic")):
            ledger.append(dict(accepted_species=s,trait_name=trait,
                normalized_value=val,evidence_scope="species_direct",
                evidence_quality="high"))
    flora=[]
    for i in range(928):
        for j in range(30):
            flora.append(dict(island_id=f"ISL{i}",accepted_species=f"Plant{j%3} species{j}",
                              origin_status="native" if j%2==0 else "unresolved"))
    pairs,summary=strict_direct_928_flora(pd.DataFrame(ledger),pd.DataFrame(flora))
    assert summary["n_original_eligible_islands"]==928
    assert len(pairs)==928*30
    geometry=gpd.GeoDataFrame({"island_id":["ISL0","ISL1","ISL2"]},
        geometry=[box(0,0,1,1),box(3,0,4,1),box(6,0,7,1)],crs="EPSG:4326")
    p=tmp_path/"source_gshhg.gpkg"
    geometry.to_file(p,layer="islands",driver="GPKG")
    receipt={
        "status":"FOCAL_21_POLYGONS_IDENTICAL_BUT_FULL_UNIVERSE_DRIFT",
        "GSHHG_archive_SHA256":PIN,
        "n_original_focal_islands":21,
        "n_focal_islands_with_exact_historical_polygon":21,
        "geometry_file_SHA256":file_sha256(p),
    }
    j=tmp_path/"geometry_receipt.json"
    j.write_text(json.dumps(receipt))
    return pairs,p,j


def test_macro_locality_positive_only_records_not_nonobservations(tmp_path):
    plants,geometry,receipt=_data(tmp_path)
    rows={
      ("Plant0 species0","flower_visiting_birds","rA","archiveA","flower_visit",.3,.4),
      ("Plant0 species0","non_bombus_bees","rB","archiveB","pollination",.3,.4),
      ("Plant1 species1","lepidoptera","rC","archiveC","flower_visit",.3,3.4),
      ("Plant2 species2","bombus","rD","archiveD","flower_visit",.3,100.0),
    }
    providers={k:{p} for *_,k,p in []}
    providers={"rA":{"archiveA"},"rB":{"archiveB"},"rC":{"archiveC"},
               "rD":{"archiveD"}}
    a=macro_exact_928_locality_gate(rows,providers,plants,geometry,receipt)
    assert a["n_original_Direct_triple_complete_islands"]==928
    assert a["n_islands_matching_original_GSHHG_geometric_id"]==3
    assert a["n_islands_unmatchable_to_original_geometry"]==925
    assert a["n_georeferenced_relations_inside_evaluable_polygon"]==3
    assert a["n_georeferenced_relations_on_exact_same_island_flora_species"]==3
    assert a["n_distinct_islands_with_species_locality_match"]==2
    assert a["n_island_species_pairs_with_2plus_documented_visitor_guilds"]==1
    assert a["n_strict_native_flora_island_species_pairs_with_visitor_relation"]==1
    assert a["n_origin_unresolved_flora_island_species_pairs_with_visitor_relation"]==1
    assert a["n_biological_pollinator_shifts_demonstrated"]==0
    assert not {"island_id","accepted_species","provider_id","reference_id"} & set(a)


def test_crossprovider_shared_bibliography_is_excluded(tmp_path):
    plants,geometry,receipt=_data(tmp_path)
    a=macro_exact_928_locality_gate(
        {("Plant0 species0","bombus","shared","A","flower_visit",.5,.5)},
        {"shared":{"A","B"}},plants,geometry,receipt)
    assert a["n_georeferenced_relations_on_exact_same_island_flora_species"]==0
    assert a["n_distinct_islands_with_species_locality_match"]==0


def test_original_geometry_SHA_mismatch_fails_closed(tmp_path):
    plants,geometry,receipt=_data(tmp_path)
    d=json.loads(receipt.read_text())
    d["GSHHG_archive_SHA256"]="wrong"
    receipt.write_text(json.dumps(d))
    with pytest.raises(ValueError,match="Geometry is not source-pin"):
        macro_exact_928_locality_gate(set(),{},plants,geometry,receipt)


def test_trait_genus_consensus_never_imputed_as_strict_species_direct(tmp_path):
    d=[]
    for i in range(30):
        for trait in ("flower_primary_color","floral_form","floral_symmetry"):
            d.append(dict(accepted_species=f"GenusA species{i}",trait_name=trait,
                normalized_value="white",evidence_scope="genus_consensus",
                evidence_quality="high"))
    flora=pd.DataFrame([
        {"island_id":f"ISL{island}",
         "accepted_species":f"GenusA species{i}","origin_status":"native"}
        for island in range(928) for i in range(30)])
    with pytest.raises(ValueError,match="All 3 original"):
        strict_direct_928_flora(pd.DataFrame(d),flora)
