"""Protected original island-endemic vs mainland-GIFT native source congener photos.

This is an executed SOURCE/photographic coverage test, not a test of
evolutionary adaptation. A source GIFT region is not an ancestral population
and public iNaturalist photographed plants need not be in that region.
The protected source taxon names/island IDs/metadata never leave CI.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import pandas as pd

from fd.oceanic_endemic_inat_photo_pilot import (
    MAX_PAIRED_GENERA,select_original_endemic_colour_pairs)
from fd.oceanic_endemic_mainland_congener_photo import (
    source_independent_mainland_controls,
    original_mainland_congener_public_photo_gate)


def run(root:Path,geo_path:Path)->dict:
    ledger=pd.read_csv(
        root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        dtype=str,usecols=["accepted_species","trait_name","normalized_value",
                           "evidence_scope","evidence_quality"])
    status=pd.read_csv(
        root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        dtype=str,usecols=["island_id","accepted_species",
                           "origin_status","endemic_status","floristic_status"])
    geology=pd.read_csv(
        geo_path,dtype=str,
        usecols=["island_id","geological_origin","analysis_regime",
                 "spatial_block","weigelt_archipelago"])
    typed=geology.loc[geology.geological_origin.isin(("oceanic","continental"))]
    if (len(typed)!=434 or
        typed.geological_origin.value_counts().get("oceanic")!=258 or
        typed.geological_origin.value_counts().get("continental")!=176):
        raise ValueError("Original Roeble/Zizka 434 geology source frame changed")
    assignments=pd.read_csv(
        root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(
        root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"])
    selected,selection=select_original_endemic_colour_pairs(
        ledger,status,geology,assignments,mainland,max_pairs=MAX_PAIRED_GENERA)
    if len(selected)!=26 or selection["n_selected_original_distinct_genera"]!=13:
        raise ValueError("Source-original paired 13 genus endemic target sample changed")
    matched,audit=source_independent_mainland_controls(
        selected,ledger,status,assignments,mainland)
    out=original_mainland_congener_public_photo_gate(matched,audit)
    out["original_qualified_island_endemic_13_genus_source"]=selection
    out["n_original_geologically_typed_islands"]=434
    out["n_frozen_endemic_to_GIFT_native_mainland_congener_pairs"]=26
    out["live_photo_provider_not_original_native_or_endemic_taxonomic_authority"]=True
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--geology-crosswalk",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    z=p.parse_args()
    result=run(z.original_root,z.geology_crosswalk)
    text=json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)
    confidential=(
        "accepted_species","island_id","entity_ID","work_species",
        "source_signature","spatial_block","source_genus","source_family",
        "photo_id","observation_id","latitude","longitude",
        "inat_taxon_id","photo_url","taxon_name","species","genus")
    if any('"'+key+'"' in text for key in confidential):
        raise ValueError("Private original plant/GIFT photo source identity escaped")
    z.output.parent.mkdir(parents=True,exist_ok=True)
    z.output.write_text(text+"\n",encoding="utf-8")
    print(json.dumps(dict(
        status=result["status"],
        n_source_original_pairs=result[
            "n_frozen_endemic_to_GIFT_native_mainland_congener_pairs"],
        n_independent_GIFT_mainland_congener_names=result[
            "n_unique_source_mainland_accepted_congener_names_queried"],
        photo_ready_by_geological_type={
            typ:r["n_source_genus_pairs_with_endemic_and_mainland_global_photo"]
            for typ,r in result["public_photo_opportunity_by_geological_island_type"].items()
        },
        no_images_measured=True,
    ),sort_keys=True))


if __name__=="__main__":main()
