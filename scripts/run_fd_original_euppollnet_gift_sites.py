"""Source-original GIFT 4505 island accepted flora × EuPPollNet sites/visitors.

Source botanicals and exact GSHHG polygons stay inside protected island CI.
EuPPollNet is a Europe-only field network dataset, not a tropical guild
experiment and not source-study-independent of GloBI/WOL/DoPI without
primary paper matching. NEVER export original taxon, island, network,
study, or exact site coordinates; all summaries are anonymous.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from urllib.request import Request,urlopen
import geopandas as gpd
import pandas as pd

from fd_island_functional_beta import strict_traits
from island_v2.fd_globi_macro_928_locality_v09 import strict_direct_928_flora
from island_v2.fd_globi_species_shift_spatial_v09 import file_sha256
from fd.fullframe_wol_network_bridge import full_botanical_exact_geometry_gate
from fd.euppollnet_public_preflight import (
    SOURCE_REPO,SOURCE_COMMIT,SOURCE_FILES,
    qualify_original_blob,available_columns
)
from fd.euppollnet_gift_exact_site_bridge import (
    evaluate_original_european_island_visits,
)

PLANTS_SHA="6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a"
DIRECT_SHA="9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745"
N_EXPECTED_SOURCE_EVENTS=623476
N_EXPECTED_SOURCE_STUDIES=54

def source_bytes(name:str)->tuple[bytes,dict]:
    if name not in SOURCE_FILES:
        raise ValueError("Unexpected unpinned original source filename")
    url=("https://raw.githubusercontent.com/"+SOURCE_REPO+"/"+
         SOURCE_COMMIT+"/Data/3_Final_data/"+name)
    req=Request(url,headers={"User-Agent":"island-source-euppollnet-field-guild-1.0"})
    with urlopen(req,timeout=180) as stream:
        data=stream.read(8000000)
    return data,qualify_original_blob(name,data)

def study(root:Path,original_geo:Path,receipt_path:Path,regime_path:Path)->dict:
    ledger=root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora_file=root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    if file_sha256(ledger)!=DIRECT_SHA or file_sha256(flora_file)!=PLANTS_SHA:
        raise ValueError("Original accepted GIFT 4505 botanical/traits source SHA drift")
    direct=pd.read_csv(ledger,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"])
    flora=pd.read_csv(flora_file,dtype=str,usecols=[
        "island_id","accepted_species","origin_status",
        "endemic_status","floristic_status"])
    if flora.island_id.nunique()!=4505:
        raise ValueError("Original botanical global island sample no longer 4505")
    old,_=strict_direct_928_flora(direct,flora)
    if old.island_id.nunique()!=928:
        raise ValueError("Historical source GIFT 928 trait complete cohort drift")
    source_receipt=json.loads(receipt_path.read_text())
    if file_sha256(original_geo)!=source_receipt.get("geometry_file_SHA256"):
        raise ValueError("Original GSHHG SHA-verified geopackage not source-authenticated")
    gpkg=gpd.read_file(original_geo,layer="islands")
    if gpkg.crs is None or gpkg.crs.to_epsg()!=4326:
        raise ValueError("Original GSHHG CRS not EPSG 4326")
    polygons,cover=full_botanical_exact_geometry_gate(
        gpkg,flora,set(old.island_id),source_receipt)
    if cover["n_original_island_floras_with_exact_hash_rebuilt_polygon"]!=4499:
        raise ValueError("Already-verified 4499 original source GIFT polygons drifted")

    relations,pin_rel=source_bytes("Interaction_data.csv.gz")
    metadata,pin_meta=source_bytes("Metadata.csv")
    from io import BytesIO
    meta=pd.read_csv(BytesIO(metadata),encoding="utf-8",
                     encoding_errors="replace",dtype=str)
    interactions=pd.read_csv(BytesIO(relations),compression="gzip",
        encoding="utf-8",encoding_errors="replace",dtype=str,low_memory=False)
    if len(interactions)!=N_EXPECTED_SOURCE_EVENTS:
        raise ValueError("Pinned original European positive interaction count drifted")
    usable={
       "Study_id","Network_id","Latitude","Longitude","Date",
       "Plant_accepted_name","Pollinator_accepted_name",
       "Pollinator_order","Pollinator_family",
    }
    if not usable<=set(interactions):
        raise ValueError("Source European original field network schema changed")
    if interactions.Study_id.nunique()!=N_EXPECTED_SOURCE_STUDIES:
        raise ValueError("Source European original study count changed")
    # Original published site/network/study/plant/visitor name encoding must
    # be usable exactly. Unknown or replacement-decoded names are quarantined.
    invalid=interactions[[
       "Study_id","Network_id","Plant_accepted_name",
       "Pollinator_accepted_name"]].fillna("").apply(
         lambda c:c.astype(str).str.contains("\ufffd",regex=False)).any(axis=1)
    dropped=int(invalid.sum())
    clean=interactions.loc[~invalid].copy()
    if not len(clean):
        raise ValueError("No source original identified positive interactions")
    regions=pd.read_csv(regime_path,dtype=str,usecols=[
        "island_id","analysis_regime"])
    if regions.island_id.duplicated().any():
        raise ValueError("Independent GIFT original geography duplicates islands")
    audit=evaluate_original_european_island_visits(
        clean,polygons[["island_id","geometry"]],flora,
        strict_traits(direct),regions.set_index("island_id").analysis_regime.to_dict(),
    )
    if audit["n_original_global_GIFT_botanical_floras"]!=4505:
        raise ValueError("GIFT original parent lost during European locality test")

    return dict(
        status="ORIGINAL_GIFT_4505_EUPPOLLNET_EUROPE_FIELD_NETWORK_CROSSWALK",
        parent_botanical_source_floras=4505,
        original_direct_ledger_SHA256=DIRECT_SHA,
        original_GIFT_accepted_flora_SHA256=PLANTS_SHA,
        original_euppollnet_source_pins=dict(interactions=pin_rel,metadata=pin_meta),
        original_euppollnet_metadata_source_rows=len(meta),
        original_source_rows_without_replacement_names=len(clean),
        original_unusable_identifier_encoded_rows_quarantined=dropped,
        original_GIFT_source_geometry_coverage=cover,
        original_european_source_interaction_count_before_quality_gate=len(interactions),
        source_original_european_network_to_GIFT_true_flora=audit,
        no_cross_provider_duplicate_study_independence_inferred=True,
        source_published_positive_visitor_record_not_effective_pollination=True,
        no_original_taxon_study_site_location_or_island_identifier_exported=True,
    )

def main():
    p=argparse.ArgumentParser()
    for key in ("original-root","original-gshhg-geometry",
                "original-geometry-receipt","original-geography","output"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    r=study(a.original_root,a.original_gshhg_geometry,
            a.original_geometry_receipt,a.original_geography)
    body=json.dumps(r,indent=2,sort_keys=True,allow_nan=False)
    for key in ("island_id","accepted_species","study_id","site_id",
                "network_id","latitude","longitude","reference_id",
                "source_region_id","source_signature","visitor_species",
                "work_species","entity_ID"):
        if '"'+key+'"' in body:
            raise ValueError("Original island/study/species/geographic ID leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(body+"\n",encoding="utf-8")
    q=r["source_original_european_network_to_GIFT_true_flora"]
    summary=dict(
        original_parent_GIFT_floras=4505,
        original_pinned_source_rows=r["original_european_source_interaction_count_before_quality_gate"],
        n_original_GIFT_geographically_supported_islands=4499,
        source_status=q["status"],
        source_network_exact_site_opportunity=q["source_verified_site_opportunity"],
        original_source_qualified_island_plant_event_rows=q.get(
            "n_original_source_positive_recorded_rows_matching_same_island_GIFT_plant",0),
        original_source_qualified_islands=q.get(
            "n_original_source_distinct_islands_with_matched_flora_plant",0),
        original_source_qualified_studies=q.get(
            "n_original_source_study_ids_with_matched_flora_plant",0),
        no_effective_pollen_or_evolution_claim=True,
    )
    print("ANONYMOUS_ORIGINAL_EUPPOLLNET_GIFT="+json.dumps(summary,sort_keys=True))
    if q["status"]=="SOURCE_EUPPOLLNET_EUROPEAN_FIELD_SITE_EXACT_GIFT_FLORA_OPPORTUNITY":
        selected=[x for x in q["true_source_direct_flower_and_animal_guild_documentation"]
            if x["region"] in ("GLOBAL_ALL","northern_midlatitude") and
            x["botanical_status"] in ("all_recorded","original_native") and
            x["n_source_study_network_site_year_positive_units"]>0]
        print("ANONYMOUS_ORIGINAL_EUPPOLLNET_MOTIFS="+json.dumps(
            selected,sort_keys=True))
if __name__=="__main__":
    main()
