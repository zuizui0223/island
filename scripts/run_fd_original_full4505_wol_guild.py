"""Full botanical GIFT 4,505 × Web of Life real positive network-site audit.

Original source-private accepted plant/island/network/site IDs never leave this
run. Historical 927-of-928 source polygon WOL result must be reproduced.
Additional 4,505 botanical flora polygons are admitted ONLY when the exact
original GSHHG geometry hash and separate signed GeoPackage receipt agree.
This is documentary site/edge opportunity, NOT local pollination effectiveness.
"""
from __future__ import annotations
from hashlib import sha256
from pathlib import Path
from collections import Counter
import argparse
import json
import pandas as pd
import requests
import geopandas as gpd

from fd_island_functional_beta import strict_traits
from island_v2.fd_globi_macro_928_locality_v09 import strict_direct_928_flora
from island_v2.fd_globi_species_shift_spatial_v09 import file_sha256
from run_fd_independent_network_site_preflight import (
    ORIGINAL_SHA, WOL_BASE, get_csv, same_original_island_network_points,
)
from fd.fullframe_wol_network_bridge import (
    full_botanical_exact_geometry_gate,
    anonymous_network_extension_counts,
    source_scoped_region_network_audit,
)
from fd.independent_pollination_networks import (
    network_location_metadata,identify_actual_plant_edges,
    source_bounded_site_plant_audit,
)

WOL_CATALOGUE_SHA="f0109a68510bf347006aa97582dd448d94daf29b91b9061e7d1c2a88850856f3"
MAX_COORDINATED_NETWORKS=200

def verify(p:Path,key:str)->None:
    if not p.is_file() or file_sha256(p)!=ORIGINAL_SHA[key]:
        raise ValueError("Immutable source GIFT botanical input changed: "+key)

def original_full_study(root:Path,geometry:Path,receipt_path:Path,
                        geography_path:Path)->dict:
    ledger=root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora_path=root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    verify(ledger,"direct_only_trait_ledger.csv.gz")
    verify(flora_path,"chapter1_status_flora.csv.gz")
    original_direct=pd.read_csv(ledger,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"])
    flora=pd.read_csv(flora_path,dtype=str,usecols=[
        "island_id","accepted_species","origin_status",
        "endemic_status","floristic_status"])
    historical,_=strict_direct_928_flora(original_direct,flora)
    original_old_ids=set(historical.island_id)
    if flora.island_id.nunique()!=4505 or len(original_old_ids)!=928:
        raise ValueError("Parent 4505 botanical or historical 928 Direct source changed")
    receipt=json.loads(receipt_path.read_text())
    if file_sha256(geometry)!=receipt.get("geometry_file_SHA256"):
        raise ValueError("Independent original polygon GeoPackage SHA256 drift")
    raw=gpd.read_file(geometry,layer="islands")
    if raw.crs is None or raw.crs.to_epsg()!=4326:
        raise ValueError("Original verified GSHHG geometry not WGS84")
    verified,coverage=full_botanical_exact_geometry_gate(
        raw,flora,original_old_ids,receipt)
    if not set(verified.island_id).issubset(set(flora.island_id)):
        raise ValueError("Original island geometry cannot create new botanical study units")
    regions=pd.read_csv(geography_path,dtype=str,usecols=[
        "island_id","analysis_regime"])
    if regions.island_id.duplicated().any():
        raise ValueError("Independent published GSHHG geographical mapping duplicated")
    mapping=regions.set_index("island_id").analysis_regime.to_dict()

    session=requests.Session()
    catalogue,catalogue_sha=get_csv(session,WOL_BASE+"/api/v1/networks.csv")
    if catalogue_sha!=WOL_CATALOGUE_SHA:
        raise ValueError("Web of Life original 2026-11-17 catalogue differs: "
                         "new snapshot requires separate original source audit")
    locations,site_census=network_location_metadata(catalogue)
    if (site_census["n_all_web_of_life_network_records"]!=330
            or site_census["n_declared_pollination_networks"]!=179
            or site_census["n_pollination_networks_with_coordinates"]!=175):
        raise ValueError("Historical Web of Life source catalogue count changed")
    matched=same_original_island_network_points(
        locations,verified[["island_id","geometry"]])
    site_counts=anonymous_network_extension_counts(
        matched,original_old_ids)
    if len(matched)>MAX_COORDINATED_NETWORKS:
        raise ValueError("Source WOL network list exceeds frozen full-retrieval cap")
    old_site_ids=set(matched.loc[
        matched.island_id.isin(original_old_ids),"network_id"])
    edges=[]
    errors=Counter()
    n_parsed_source_edges=0
    n_taxon_classified_source_edges=0
    raw_net_sha=[]
    for item in matched.sort_values("network_id").itertuples(index=False):
        prefix=WOL_BASE+"/api/v1/networks/"+str(item.network_id)
        try:
            record,rs=get_csv(session,prefix+"/edges.csv")
            taxonomy,ts=get_csv(session,prefix+"/species.csv")
            positive,quality=identify_actual_plant_edges(
                record,taxonomy)
            n_parsed_source_edges+=quality["n_explicit_plant_role_binomial_edge_records"]
            n_taxon_classified_source_edges+=quality[
                "n_taxonomy_classified_guild_edge_records"]
            raw_net_sha.append(sha256((rs+"|"+ts).encode()).hexdigest())
            if len(positive):
                positive["network_id"]=str(item.network_id)
                positive["island_id"]=str(item.island_id)
                edges.append(positive[[
                    "network_id","island_id","plant","guild"]])
        except (requests.RequestException,ValueError,KeyError) as exc:
            errors[type(exc).__name__]+=1
    if errors:
        raise ValueError("Full original GIFT GSHHG-site Web of Life network "
                         "retrieval not complete: "+str(dict(errors)))
    source_edges=pd.concat(edges,ignore_index=True) if edges else pd.DataFrame(
        columns=["network_id","island_id","plant","guild"])
    # The historical old 927 polygon cohort must reproduce the old source
    # EXACT network-site botanical accepted plant and guild count (2469).
    old_edges=source_edges[source_edges.network_id.isin(old_site_ids)]
    old_flora=flora[flora.island_id.isin(original_old_ids)]
    historical_audit=source_bounded_site_plant_audit(
        old_edges,old_flora,strict_traits(original_direct),
        global_original_islands=4505)
    if (historical_audit["n_original_external_networks_with_plant_guild_overlap"]<=0
            or historical_audit["n_original_species_network_plant_guild_units"]!=2469):
        raise ValueError("Previous source-original 927-island Web of Life "
                         "accepted plant×guild union failed to reproduce")
    supported=flora[flora.island_id.isin(set(verified.island_id))]
    full=source_bounded_site_plant_audit(
        source_edges,supported,strict_traits(original_direct),
        global_original_islands=4505)
    if full["n_original_species_network_plant_guild_units"] < 2469:
        raise AssertionError("Full GIFT frame cannot lose historical visitor support")
    regional=source_scoped_region_network_audit(
        source_edges,supported,strict_traits(original_direct),
        mapping,source_audit=site_counts)
    return dict(
        status="ORIGINAL_GIFT_4505_EXACT_HASH_FULLFRAME_WOL_LOCAL_POLLINATOR_SUPPORT",
        parent_botanical_floras_original=4505,
        source_web_of_life_catalogue_sha256=catalogue_sha,
        source_gshhg_GIFT_geometry_file_sha256=file_sha256(geometry),
        original_flora_source_sha256=file_sha256(flora_path),
        original_species_Direct_source_sha256=file_sha256(ledger),
        original_real_botanical_geometry_coverage=coverage,
        verified_published_network_site_coverage=site_counts,
        source_WOL_site_index=site_census,
        n_network_source_edges_with_original_role_and_binomial=int(n_parsed_source_edges),
        n_source_edges_with_curated_animal_taxonomy=int(n_taxon_classified_source_edges),
        n_network_source_files_retrieved=int(len(matched)),
        n_network_api_errors=int(sum(errors.values())),
        original_network_csv_snapshot_pair_count=len(raw_net_sha),
        historical_927_source_pair_reconciliation=dict(
            historical_927_networks=75,historical_original_island_sites=16,
            original_accepted_plant_network_guild_units=2469,
            exact_replay_passed=True,
        ),
        full_original_botanical_flora_island_visitor_evidence=full,
        original_source_island_regional_guild_and_real_flower_pairs=regional,
        no_effective_pollination_or_plant_local_phenotype_assumption=True,
        strict_native_status_not_inferred_from_unresolved=True,
        original_species_network_positive_only_no_negative_visits=True,
        source_WOL_study_independence_from_GloBI_Mangal_NOT_established=True,
        original_global_GSHHG_nine_unreproduced_IDs_not_imputed=True,
    )

def main():
    a=argparse.ArgumentParser()
    for key in ("original-root","original-gshhg-geometry",
                "original-geometry-receipt","original-geography","output"):
        a.add_argument("--"+key,type=Path,required=True)
    args=a.parse_args()
    result=original_full_study(args.original_root,args.original_gshhg_geometry,
        args.original_geometry_receipt,args.original_geography)
    body=json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False,
                    allow_nan=False)
    for key in ("island_id","accepted_species","network_id","site_id",
                "latitude","longitude","source_region_id","spatial_block",
                "reference_id","provider_id","source_genus","visitor_species",
                "taxon_id","work_species","entity_ID"):
        if '"'+key+'"' in body:
            raise ValueError("Private original botanical/network source identifier leaked")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+"\n",encoding="utf-8")
    print("ANONYMOUS_FULL4505_WOL="+json.dumps({
        "geometry":result["original_real_botanical_geometry_coverage"],
        "sites":result["verified_published_network_site_coverage"],
        "full_global":{
           "n_documented_networks":result["full_original_botanical_flora_island_visitor_evidence"][
               "n_original_external_networks_with_plant_guild_overlap"],
           "n_documented_islands":result["full_original_botanical_flora_island_visitor_evidence"][
               "n_original_external_network_islands_with_plant_guild_overlap"],
           "n_plant_guild_units":result["full_original_botanical_flora_island_visitor_evidence"][
               "n_original_species_network_plant_guild_units"],
        },
        "source_edges":result["n_network_source_edges_with_original_role_and_binomial"],
        "regional": [{
           "region":r["region"],
           "n_islands_with_botanical_positive_networks":r[
               "n_original_islands_with_source_positive_plant_guild_documentation"],
           "n_architecture_guild_units":sum(x["n_plant_guild_network_documented_units"]
              for x in r["true_direct_paired_flower_guild_documentation_cells"])
        } for r in result["original_source_island_regional_guild_and_real_flower_pairs"][
                "regional_original_source_opportunity"]],
    },sort_keys=True))
    print(json.dumps(dict(status=result["status"],
        original_parent_islands=4505,
        no_animal_absence_or_fitness_conclusion=True)))
if __name__=="__main__":
    main()
