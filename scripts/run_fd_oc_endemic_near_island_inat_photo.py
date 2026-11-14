"""Private source-exact endemic public photos near geological island centroids.

Original island native/endemic names, GIFT exact mainland source IDs,
geological archipelago and Direct plant colour source data stay in CI.
A bounded iNaturalist official REST lookup queries at most 13 matched
genera (16 target accepted endemic species); NO pixels are downloaded
or classified; no observation IDs, photo IDs or endemic taxon names
leave source CI. Missing provider data is NOT absence of wild photos.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import pandas as pd
from fd.oceanic_endemic_inat_photo_pilot import (
    MAX_PAIRED_GENERA,select_original_endemic_colour_pairs)
from fd.oceanic_endemic_near_island_photo import original_near_island_photo_provider_pilot


def run(root:Path,geology:Path,*,max_pairs:int=MAX_PAIRED_GENERA)->dict:
    ledger=pd.read_csv(
        root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        dtype=str,usecols=["accepted_species","trait_name","normalized_value",
                           "evidence_scope","evidence_quality"])
    status=pd.read_csv(
        root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        dtype=str,usecols=["island_id","accepted_species","origin_status",
                           "endemic_status","floristic_status"])
    geo=pd.read_csv(geology,dtype=str,
        usecols=["island_id","geological_origin","analysis_regime",
                 "spatial_block","weigelt_archipelago",
                 "island_latitude","island_longitude","area_km2"])
    typed=geo.loc[geo.geological_origin.isin(("oceanic","continental"))]
    if typed.island_id.nunique()!=434:
        raise ValueError("Source verified 434 geological island universe drift")
    if typed.geological_origin.value_counts().get("oceanic")!=258 or (
        typed.geological_origin.value_counts().get("continental")!=176):
        raise ValueError("Independently sourced oceanic/continental geology labels drift")
    source=pd.read_csv(
        root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    native=pd.read_csv(
        root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"])
    selected,selection=select_original_endemic_colour_pairs(
        ledger,status,geo,source,native,max_pairs=max_pairs)
    if not len(selected)==2*max_pairs or (
        selection["n_original_distinct_source_pair_eligible_genera"]<max_pairs):
        raise ValueError("Original source pairs insufficient for outcome-blind taxon pilot")
    metadata=geo.set_index("island_id")
    if metadata.index.has_duplicates:
        raise ValueError("Original GSHHG geologic centroids not unique")
    for row in selected:
        key=row["island_id"]
        if key not in metadata.index:
            raise ValueError("Original selected island has no audited GSHHG centroid")
        row["island_latitude"]=float(metadata.at[key,"island_latitude"])
        row["island_longitude"]=float(metadata.at[key,"island_longitude"])
        row["island_area_km2"]=float(metadata.at[key,"area_km2"])
    output=original_near_island_photo_provider_pilot(selected,selection)
    output["n_original_GSHHG_geologically_typed_islands"]=434
    output["n_original_tropical_endemic_candidates_chosen_before_API"]=2*max_pairs
    output["live_public_provider_not_the_frozen_source_of_taxonomic_endemism"]=True
    return output


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--geology-crosswalk",type=Path,required=True)
    p.add_argument("--max-pairs",type=int,default=13)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    result=run(a.original_root,a.geology_crosswalk,max_pairs=a.max_pairs)
    body=json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)
    forbidden=(
       "accepted_species","island_id","entity_ID","work_species",
       "source_signature","spatial_block","source_genus",
       "photo_id","observation_id","latitude","longitude","inat_taxon_id",
       "photo_url","taxon_name","genus","species")
    if any('"'+name+'"' in body for name in forbidden):
        raise ValueError("Private source biological or photo identity leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(body+"\n",encoding="utf-8")
    print(json.dumps({
        "status":result["status"],
        "n_source_paired_genera":result[
            "original_selection_before_photo_query"]["n_selected_original_distinct_genera"],
        "n_public_near_centroid_photo_taxa_by_geology":{
            t:z["n_public_near_island_centroid_qualified_photo_taxa"]
            for t,z in result["original_public_geo_photo_provider_by_island_type"].items()},
        "n_public_near_centroid_photo_matched_genera":result[
            "n_matched_genera_with_centroid_near_public_photo_in_both_geologies"],
        "source_geographic_stage":"GSHHG centroid radius NOT original exact polygon",
        "no_photo_pixels_or_colour_roi_measured":True,
    },ensure_ascii=False,sort_keys=True))


if __name__=="__main__":
    main()
