"""Source-protected 4505 GIFT × EuPPollNet actual flower-count support.

The original European public source Flower_counts.csv and positive interaction
rows are Git-blob-pinned, every original GIFT polygon/accepted plant ID is
SHA-authenticated; same study, network, published calendar day and unique
flower-survey key required. No plant/study/island/site identity exported.
This is availability support only, NOT preference/pollen success estimation.
"""
from __future__ import annotations

from io import BytesIO
from pathlib import Path
from urllib.request import Request,urlopen
import argparse
import json
import pandas as pd
import geopandas as gpd

from fd_island_functional_beta import strict_traits
from island_v2.fd_globi_macro_928_locality_v09 import strict_direct_928_flora
from island_v2.fd_globi_species_shift_spatial_v09 import file_sha256
from fd.euppollnet_public_preflight import (
    SOURCE_REPO,SOURCE_COMMIT,qualify_original_blob
)
from fd.euppollnet_flower_availability import (
    source_original_flower_receipt,
    qualified_gift_flower_supply_on_positive_visits,
)
from fd.euppollnet_gift_exact_site_bridge import (
    exact_source_network_site_union,strict_original_accepted_plant
)
from fd.fullframe_wol_network_bridge import full_botanical_exact_geometry_gate

FLORA="6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a"
TRAITS="9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745"

def public_bytes(name,limit):
    if name not in ("Flower_counts.csv","Interaction_data.csv.gz"):
        raise ValueError("Unexpected public source file")
    u=("https://raw.githubusercontent.com/"+SOURCE_REPO+"/"+
       SOURCE_COMMIT+"/Data/3_Final_data/"+name)
    with urlopen(Request(u,headers={
        "User-Agent":"GIFT-4505-European-flower-supply-source-replay"}),
        timeout=150) as response:
        blob=response.read(limit+1)
    if not blob or len(blob)>limit:
        raise ValueError("Source original public file truncated or over cap")
    return blob

def run(original_root:Path,polygon_path:Path,receipt_path:Path):
    traits_path=original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora_path=original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    if file_sha256(traits_path)!=TRAITS or file_sha256(flora_path)!=FLORA:
        raise ValueError("GIFT original accepted plants or Direct species traits changed")
    t=pd.read_csv(traits_path,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"])
    f=pd.read_csv(flora_path,dtype=str,usecols=[
        "island_id","accepted_species","origin_status","endemic_status"])
    if f.island_id.nunique()!=4505:
        raise ValueError("Source global GIFT flora sample drifted")
    old,_=strict_direct_928_flora(t,f)
    if old.island_id.nunique()!=928:
        raise ValueError("Original selected triple-Direct 928 flora count drifted")
    meta=json.loads(receipt_path.read_text())
    if file_sha256(polygon_path)!=meta.get("geometry_file_SHA256"):
        raise ValueError("Original GSHHG geopackage receipt failed")
    original_geo=gpd.read_file(polygon_path,layer="islands")
    if original_geo.crs is None or original_geo.crs.to_epsg()!=4326:
        raise ValueError("Source GSHHG geometry not WGS84")
    polygons,geo=full_botanical_exact_geometry_gate(
        original_geo,f,set(old.island_id),meta)
    if geo["n_original_island_floras_with_exact_hash_rebuilt_polygon"]!=4499:
        raise ValueError("Source 4499 GIFT exact polygon baseline drifted")

    original_visits=public_bytes("Interaction_data.csv.gz",5000000)
    original_flowers=public_bytes("Flower_counts.csv",25000000)
    original_receipts=dict(
        visits=qualify_original_blob("Interaction_data.csv.gz",original_visits),
        flowers=source_original_flower_receipt(original_flowers),
    )
    source=pd.read_csv(BytesIO(original_visits),compression="gzip",
        dtype=str,low_memory=False,encoding="utf-8",encoding_errors="replace")
    flowers=pd.read_csv(BytesIO(original_flowers),dtype=str,
        low_memory=False,encoding="utf-8",encoding_errors="replace")
    if len(source)!=623476 or len(flowers)!=199692:
        raise ValueError("Original published positive-visit/flower count rows drifted")
    for fld in ("Study_id","Network_id","Latitude","Longitude",
                 "Date","Plant_accepted_name","Pollinator_accepted_name",
                 "Flower_data_merger"):
        if fld not in source:
            raise ValueError("Original source study/site/plant/flower link missing")
    source["Latitude"]=pd.to_numeric(source.Latitude,errors="coerce")
    source["Longitude"]=pd.to_numeric(source.Longitude,errors="coerce")
    valid=(source.Latitude.between(-90,90)&source.Longitude.between(-180,180))
    if not valid.all():
        raise ValueError("Source EuPPollNet originally valid 623476 visit coordinates drifted")
    locations,site_counts=exact_source_network_site_union(
        source,polygons[["island_id","geometry"]])
    selected=source.merge(locations,on=[
        "Study_id","Network_id","Latitude","Longitude"],how="inner",
        validate="many_to_one")
    unique=set(selected.Plant_accepted_name.dropna())
    canonical={k:strict_original_accepted_plant(k) for k in unique}
    selected["plant"]=selected.Plant_accepted_name.map(canonical)
    selected=selected.loc[selected.plant.notna()].copy()
    site=selected.merge(f[[
        "island_id","accepted_species","origin_status","endemic_status"]]
        .rename(columns={"accepted_species":"plant"}),
        on=["island_id","plant"],how="inner",validate="many_to_one")
    if len(site)!=36384 or site.island_id.nunique()!=13:
        raise ValueError("Earlier original 13-island/36384 positive EuPPollNet source link drifted")
    supply=qualified_gift_flower_supply_on_positive_visits(
        site,flowers,strict_traits(t))
    if supply["parent_original_GIFT_botanical_floras"]!=4505:
        raise ValueError("Original GIFT sample no longer full botanical 4505")
    return dict(
        status="GIFT_4505_ORIGINAL_EUPPOLLNET_FIELD_FLOWER_SUPPLY_FEASIBILITY",
        original_source_GIFT_flora_SHA256=FLORA,
        original_source_direct_traits_SHA256=TRAITS,
        original_EuPPollNet_public_Git_blob_receipts=original_receipts,
        source_original_GIFT_verified_polygon_geography=geo,
        source_original_published_site_polygons=site_counts,
        n_original_raw_european_visits=623476,
        n_original_raw_european_flower_surveys=199692,
        n_replayed_original_source_positive_GIFT_island_plant_events=len(site),
        n_replayed_GIFT_island_plant_positive_islands=13,
        source_actual_original_flower_availability_on_GIFT_positive_visits=supply,
        source_positive_is_not_effective_pollination_or_true_animal_preference=True,
        unknown_source_study_overlap_with_GloBI_WOL_DoPI=True,
        no_source_study_island_plant_site_flower_key_ids_exported=True)

def main():
    p=argparse.ArgumentParser()
    for arg in ("original-root","original-gshhg-geometry",
                "original-geometry-receipt","output"):
        p.add_argument("--"+arg,type=Path,required=True)
    a=p.parse_args()
    report=run(a.original_root,a.original_gshhg_geometry,
               a.original_geometry_receipt)
    body=json.dumps(report,sort_keys=True,indent=2,allow_nan=False)
    for forbidden in ("accepted_species","island_id","source_signature",
        "Study_id","Network_id","site_id","latitude","longitude",
        "reference_id","provider_id","source_region_id","work_species"):
        if '"'+forbidden+'"' in body:
            raise ValueError("Source GIFT island/flower/plant/study key leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(body+"\n",encoding="utf-8")
    r=report["source_actual_original_flower_availability_on_GIFT_positive_visits"]
    print("ANONYMOUS_GIFT_EUPPOLLNET_FLOWER_SUPPLY="+json.dumps({
        "status":report["status"],"original_parent_floras":4505,
        "original_replayed_islands":13,
        "original_site_positive_visit_rows":r.get(
            "n_original_source_positive_same_island_botanical_visits",0),
        "rows_unique_flower_key":r.get(
            "n_original_visitor_rows_with_unique_original_flower_key",0),
        "rows_matching_day_and_study":r.get(
            "n_original_visitors_with_same_published_study_network_calendar_day",0),
        "rows_positive_numeric_flower_exposure":r.get(
            "n_original_positive_visitors_with_nonzero_numeric_flower_count_and_unit",0),
        "n_original_GIFT_islands_with_flower_supply":r.get(
            "n_original_GIFT_islands_with_qualified_flower_exposure",0),
        "n_original_studies_with_flower_supply":r.get(
            "n_original_source_studies_with_qualified_flower_exposure",0),
        "n_original_botanical_plants_with_flower_supply":r.get(
            "n_original_accepted_botanical_plants_with_qualified_flower_exposure",0),
        "source_original_count_units":r.get(
            "original_qualified_source_flower_unit_labels_counts",{}),
        "no_trait_guild_preference_coefficient_estimated":True,
    },sort_keys=True))
    if r.get("status")=="ORIGINAL_GIFT_EUPPOLLNET_EXACT_POSITIVE_FLOWER_SUPPLY_LINK_AUDIT":
        for x in r["true_source_same_species_flower_x_animal_positive_exposure_cells"]:
            if x["botanical_scope"]=="all_recorded" and x[
                "n_original_qualified_flower_count_linked_visitor_rows"]>0:
                print("ANONYMOUS_GIFT_EUPPOLLNET_FLOWER_MOTIF="+json.dumps(
                    x,sort_keys=True))
if __name__=="__main__":
    main()
