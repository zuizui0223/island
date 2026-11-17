"""Independent Web of Life + Mangal geographically verified network preflight.

Run ONLY in original protected island GitHub Actions, from downloaded public
network publication archives and exact source-hash GIFT flora and historical
GSHHG polygons. Do not publicly release GIFT accepted names, IDs, coordinates
or source matches. 4,505-island original macro population is NEVER replaced
by the 927-polygon observationally covered subset.

For Web of Life, require documented plant/animal roles in each OBSERVED edge
and original species names in the exact island flora. Mangal contributes ONLY
a potentially relevant point/site metadata inventory until independent network
nodes and positive interaction type are checked. No implicit pseudoabsences.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import io
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from fd_island_functional_beta import strict_traits
from fd_direct_three_stage_floristic_comparison import status_groups
from island_v2.fd_globi_macro_928_locality_v09 import strict_direct_928_flora
from island_v2.fd_globi_species_shift_spatial_v09 import PIN as GSHHG_PIN, file_sha256
from fd.independent_pollination_networks import (
    network_location_metadata,mangal_location_metadata,
    identify_actual_plant_edges,source_bounded_site_plant_audit,
)

ORIGINAL_SHA={
    "direct_only_trait_ledger.csv.gz":
      "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
    "chapter1_status_flora.csv.gz":
      "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
}
WOL_BASE="https://www.web-of-life.es"
MANGAL_BASE="https://mangal.io/api/v2"
MAX_WOL_NETWORKS=150
MAX_NETWORK_CSV_BYTES=14000000


def digest(path:Path)->str:
    return file_sha256(path)


def get_bytes(session:requests.Session,url:str,*,max_bytes:int=MAX_NETWORK_CSV_BYTES):
    response=session.get(url,timeout=(15,90),headers={
        "Accept":"text/csv, application/json, */*",
        "User-Agent":"fd-independent-network-source-preflight/1.0 (source-reproducibility)"
    })
    response.raise_for_status()
    blob=response.content
    if not blob or len(blob)>max_bytes:
        raise ValueError("independent source empty or over explicit size cap")
    return blob


def get_csv(session,url):
    blob=get_bytes(session,url)
    df=pd.read_csv(io.BytesIO(blob),low_memory=False)
    return df,sha256(blob).hexdigest()


def get_all_mangal(session,section:str,max_pages=4)->tuple[list,dict]:
    all_rows=[]
    hashes=[]
    for page in range(max_pages):
        blob=get_bytes(
            session,f"{MANGAL_BASE}/{section}?count=1000&page={page}",
            max_bytes=18000000
        )
        hashes.append(sha256(blob).hexdigest())
        body=json.loads(blob)
        if isinstance(body,list):rows=body
        elif isinstance(body,dict):
            rows=body.get("data",body.get("results"))
            if not isinstance(rows,list):
                raise ValueError("Mangal API data not a list")
        else:raise ValueError("Mangal source schema changed")
        all_rows.extend(rows)
        if len(rows)<1000:break
    return all_rows,dict(
        original_page_SHA256=hashes,
        n_original_rows_fetched=len(all_rows),
        is_exhaustive_with_page_cap=len(rows)<1000,
        page_cap=max_pages,
    )


def strict_source_geometry(geo_path:Path,receipt_path:Path,
                           direct:pd.DataFrame,flora:pd.DataFrame):
    import geopandas as gpd
    if not geo_path.is_file() or not receipt_path.is_file():
        raise ValueError("original verified GSHHG polygons unavailable")
    receipt=json.loads(receipt_path.read_text())
    if (receipt.get("GSHHG_archive_SHA256")!=GSHHG_PIN or
        receipt.get("geometry_file_SHA256")!=file_sha256(geo_path) or
        receipt.get("n_focal_islands_with_exact_historical_polygon")!=21):
        raise ValueError("independent original historical polygon SHA or source receipt drift")
    full=status_groups(flora)
    if len(full)!=4505:raise ValueError("original GIFT 4505 botanical island population changed")
    observed,_=strict_direct_928_flora(direct,flora)
    if observed.island_id.nunique()!=928:
        raise ValueError("original 928 Direct-flora documentary support subset drift")
    geo=gpd.read_file(geo_path,layer="islands")
    if (geo.crs is None or geo.crs.to_epsg()!=4326 or
        geo.island_id.isna().any() or geo.island_id.duplicated().any()):
        raise ValueError("original island polygon identity / WGS84 changed")
    polygon=geo.loc[
        geo.island_id.isin(set(observed.island_id))
    ][["island_id","geometry"]].copy()
    if polygon.island_id.nunique()!=927:
        raise ValueError("original GIFT 927 authenticated polygon subset drift")
    if not set(polygon.island_id)<=set(full):
        raise ValueError("polygon is not one of original botanical GIFT islands")
    return polygon,full


def same_original_island_network_points(frame:pd.DataFrame,polygon)->pd.DataFrame:
    import geopandas as gpd
    if frame.empty:
        return pd.DataFrame(columns=["network_id","island_id"])
    if not frame.network_id.is_unique:
        raise ValueError("network ID duplication means location unclear")
    points=gpd.GeoDataFrame(
        frame[["network_id"]].copy(),
        geometry=gpd.points_from_xy(frame.lon,frame.lat),
        crs="EPSG:4326"
    )
    located=gpd.sjoin(points,polygon,how="inner",predicate="within")
    n=located.groupby("network_id").island_id.nunique()
    unique_ids=set(n.index[n.eq(1)])
    matched=located.loc[located.network_id.isin(unique_ids),
                        ["network_id","island_id"]].drop_duplicates()
    if matched.network_id.duplicated().any():
        raise ValueError("ambiguous original island for a public network site")
    return matched.reset_index(drop=True)


def run(original_root:Path,geo:Path,receipt:Path)->dict:
    ledger=original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    florapath=original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    for path in (ledger,florapath):
        if not path.exists() or digest(path)!=ORIGINAL_SHA[path.name]:
            raise ValueError("original botanical Direct GIFT immutable source SHA failure")
    traits=pd.read_csv(ledger,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"
    ])
    flora=pd.read_csv(florapath,dtype=str,usecols=[
        "island_id","accepted_species","origin_status"
    ])
    polygon,full=strict_source_geometry(geo,receipt,traits,flora)
    original_direct=strict_traits(traits)

    session=requests.Session()
    catalogue_errors=[]
    WOL=None;wol_sha=None;wolmeta={}
    for url in (WOL_BASE+"/api/v1/networks.csv",
                WOL_BASE+"/get_network_info.php"):
        try:
            raw,sh=get_csv(session,url)
            candidate,meta=network_location_metadata(raw)
            WOL=candidate;wol_sha=sh;wolmeta=meta
            wolmeta["original_network_source_url"]=url
            break
        except (requests.RequestException,ValueError,KeyError) as exc:
            catalogue_errors.append(type(exc).__name__+":"+str(exc)[:140])

    dataset=[];networks=[];mangalmeta={}
    try:
        dataset,ds_meta=get_all_mangal(session,"dataset")
        networks,net_meta=get_all_mangal(session,"network")
        m,mangalmeta=mangal_location_metadata(networks,dataset)
        mangalmeta["original_data_page_integrity"]={
            "datasets":ds_meta,"networks":net_meta
        }
        mangalmeta["source_retrieval_status"]="RETRIEVED_METADATA_NOT_VALIDATED_POLLINATION_EDGES"
        matchm=same_original_island_network_points(m,polygon)
        mangalmeta["n_keyword_candidates_in_original_927_polygons"]=int(
            len(matchm))
        mangalmeta["n_original_islands_with_keyword_candidates"]=int(
            matchm.island_id.nunique())
    except (requests.RequestException,ValueError,KeyError) as exc:
        mangalmeta=dict(source_retrieval_status="SOURCE_UNAVAILABLE_OR_UNSUPPORTED_SCHEMA",
                        original_error_class=type(exc).__name__,
                        n_source_confirmed_visitor_edges=None)

    edges_list=[]
    wol_node_summary=dict(
        source_retrieval_status="SOURCE_UNAVAILABLE_OR_UNSUPPORTED_SCHEMA"
        if WOL is None else "RETRIEVED_SITE_METADATA",
        n_original_network_site_polygon_matches=0,
        n_original_island_network_sites=0,
        n_verified_web_of_life_plant_role_edges=0,
        n_public_networks_with_valid_plant_roles=0,
        n_networks_with_unsupported_edge_schema=0,
        n_networks_downloaded=0,
        n_networks_skipped_by_explicit_download_cap=0,
    )
    if WOL is not None:
        match=same_original_island_network_points(WOL,polygon)
        wol_node_summary["n_original_network_site_polygon_matches"]=int(len(match))
        wol_node_summary["n_original_island_network_sites"]=int(
            match.island_id.nunique())
        wol_node_summary["original_wol_network_catalogue_sha256"]=wol_sha
        wol_node_summary.update(wolmeta)
        to_load=match.sort_values("network_id")
        wol_node_summary["n_networks_skipped_by_explicit_download_cap"]=int(
            max(0,len(to_load)-MAX_WOL_NETWORKS))
        schemas={}
        errors=[]
        for row in to_load.head(MAX_WOL_NETWORKS).itertuples(index=False):
            prefix=f"{WOL_BASE}/api/v1/networks/{row.network_id}"
            try:
                edges,es=get_csv(session,prefix+"/edges.csv")
                names,ns=get_csv(session,prefix+"/species.csv")
                good,quality=identify_actual_plant_edges(edges,names)
                schemas[row.network_id]={
                    "edge_headers":list(edges.columns),
                    "species_headers":list(names.columns),
                    "edge_source_sha256":es,
                    "species_source_sha256":ns,
                }
            except (requests.RequestException,ValueError,KeyError) as exc:
                wol_node_summary["n_networks_with_unsupported_edge_schema"]+=1
                errors.append(type(exc).__name__+":"+str(exc)[:160])
                continue
            wol_node_summary["n_networks_downloaded"]+=1
            wol_node_summary["n_verified_web_of_life_plant_role_edges"]+=len(good)
            if len(good):
                wol_node_summary["n_public_networks_with_valid_plant_roles"]+=1
                good["network_id"]=row.network_id
                good["island_id"]=row.island_id
                edges_list.append(good[["network_id","island_id","plant","guild"]])
        wol_node_summary["source_errors_number"]=len(errors)
        wol_node_summary["source_error_categories_sample"]=errors[:5]
        wol_node_summary["n_source_networks_with_verified_header_snapshot"]=len(schemas)
        # The headers of a public data provider are a useful machine-readable
        # protocol, but source network IDs/plant IDs never join the public
        # output because these are used to infer private GIFT island matches.
        wol_node_summary["original_edge_header_sets"]=[
            sorted(set(x["edge_headers"])) for x in list(schemas.values())[:3]]
        wol_node_summary["original_species_header_sets"]=[
            sorted(set(x["species_headers"])) for x in list(schemas.values())[:3]]
        if wol_node_summary["n_networks_downloaded"]==0 and len(match):
            wol_node_summary["source_retrieval_status"]="MATCHED_NETWORK_CSV_UNAVAILABLE"
        else:
            wol_node_summary["source_retrieval_status"]="VERIFIED_ONLY_POSITIVE_EDGES_AND_SOURCE_SITE"

    matched_edges=pd.concat(edges_list,ignore_index=True) if edges_list else pd.DataFrame(
        columns=["network_id","island_id","plant","guild"])
    flora=flora.loc[flora.island_id.isin(set(polygon.island_id))].copy()
    original_audit=source_bounded_site_plant_audit(
        matched_edges,flora,original_direct,global_original_islands=4505
    )
    return dict(
        status="ORIGINAL_4505_INDEPENDENT_EXTERNAL_POLLINATION_NETWORK_PREFLIGHT",
        n_original_botanical_GIFT_island_floras=4505,
        n_original_geometrically_evaluable_Direct_observation_subset=927,
        n_original_direct_plant_phenotype_observation_floras=928,
        WOL_network_metadata=wol_node_summary,
        WOL_catalogue_retrieval_errors=catalogue_errors,
        Mangal_metadata=mangalmeta,
        exact_site_and_original_accepted_plant_WOL_evidence=original_audit,
        source_GIFT_flora_sha256=digest(florapath),
        source_Direct_species_trait_ledger_sha256=digest(ledger),
        methodological_ceiling=dict(
            source_web_of_life_and_mangal_different_aggregators_not_proven_new_original_studies=True,
            source_network_site_gps_precision_not_certified_by_original_metadata=True,
            source_network_observed_edges_not_effective_pollen_transfer=True,
            only_927_of_original_4505_geometrically_evaluable_under_source_contract=True,
            Web_of_Life_species_plant_role_and_accepted_binomial_match_required=True,
            Mangal_keyword_site_candidates_are_NOT_verified_plant_pollinator_interactions=True,
            no_unobserved_visitors_treated_as_negative_records=True,
            original_citation_study_level_GloBI_overlap_not_independently_audited=True,
            source_web_API_snapshot_mutable_not_historical_original_research_lock=True,
            no_island_trait_evolution_or_guild_causal_shift_inferred=True,
        )
    )


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--original-gshhg-geometry",type=Path,required=True)
    p.add_argument("--original-geometry-receipt",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    x=p.parse_args()
    result=run(x.original_root,x.original_gshhg_geometry,
               x.original_geometry_receipt)
    body=json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False,
                    allow_nan=False)
    for forbidden in ("island_id","accepted_species","work_species",
        "entity_ID","spatial_block","source_signature","source_region_id",
        "reference_id","provider_id","latitude","longitude","genus",
        "source_genus","visitor_species"):
        if '"'+forbidden+'"' in body:
            raise ValueError("protected botanical source/site taxon IDs leaked")
    x.output.parent.mkdir(parents=True,exist_ok=True)
    x.output.write_text(body+"\n",encoding="utf-8")
    print(json.dumps(dict(
        status=result["status"],original_botanical_islands=4505,
        WOL_source_status=result["WOL_network_metadata"]["source_retrieval_status"],
        WOL_georeferenced_networks_inside_authenticated_polygons=result[
            "WOL_network_metadata"]["n_original_network_site_polygon_matches"],
        WOL_real_plant_guild_units=result[
            "exact_site_and_original_accepted_plant_WOL_evidence"][
                "n_original_species_network_plant_guild_units"],
        mangal_candidate_site_count=result["Mangal_metadata"].get(
            "n_keyword_candidates_in_original_927_polygons"),
        no_effective_pollination_inferred=True,
    ),sort_keys=True))


if __name__=="__main__":
    main()
