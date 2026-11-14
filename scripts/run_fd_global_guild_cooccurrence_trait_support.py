"""Original 4,505 botanical islands × GBIF pollinator guild co-observation preflight.

A qualified non-detection, detected other guild and Direct plant trait
documentation on SAME island constitutes a study OPPORTUNITY only. It is
NOT an observed same-plant pollinator switch or functional replacement.
Original source island, plant, GIFT, GBIF animal and geolocation identities
stay inside protected original source workflow; aggregate JSON only.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import pandas as pd

from fd.global_guild_cooccurrence_trait_support import guild_cooccurrence_global


def file_sha(path:Path)->str:
    d=sha256()
    with path.open("rb") as f:
        for buf in iter(lambda:f.read(1<<20),b""):
            d.update(buf)
    return d.hexdigest()


def run(root:Path,geo_path:Path,channels_path:Path)->dict:
    botan=root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    ledger_path=root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    if file_sha(botan)!=(
        "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a"):
        raise ValueError("Original full flora source 4505 changed")
    if file_sha(ledger_path)!=(
        "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745"):
        raise ValueError("Original source Direct trait ledger changed")
    flora=pd.read_csv(botan,dtype=str,
        usecols=["island_id","accepted_species","origin_status"])
    if flora.island_id.nunique()!=4505:
        raise ValueError("Global original GIFT botanical island membership no longer 4505")
    ledger=pd.read_csv(ledger_path,dtype=str,
        usecols=["accepted_species","trait_name","normalized_value",
                 "evidence_scope","evidence_quality"])
    geo=pd.read_csv(geo_path,dtype=str,
        usecols=["island_id","geological_origin",
                 "analysis_regime","spatial_block"])
    if (geo.island_id.nunique()<8200 or
        geo.geological_origin.eq("oceanic").sum()!=258 or
        geo.geological_origin.eq("continental").sum()!=176):
        raise ValueError("Original independent Roeble/Zizka geology source changed")
    channels=pd.read_csv(channels_path,dtype=str,
        usecols=["island_id","channel_id","channel_state","source_state"])
    output=guild_cooccurrence_global(flora,ledger,geo,channels)
    if output["n_original_global_plant_flora_islands"]!=4505:
        raise AssertionError("Original all flora denominator changed")
    output["original_GBIF_channel_source_SHA256"]=file_sha(channels_path)
    output["original_GIFT_4505_flora_SHA256"]=file_sha(botan)
    output["source_original_data_not_exported"]=True
    return output


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--geology-crosswalk",type=Path,required=True)
    p.add_argument("--channel-evidence-csv",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    result=run(args.original_root,args.geology_crosswalk,
               args.channel_evidence_csv)
    body=json.dumps(result,sort_keys=True,indent=2,allow_nan=False)
    forbidden=(
        "island_id","accepted_species","entity_ID","work_species",
        "source_signature","spatial_block","source_region_id",
        "longitude","latitude","taxon_id","GBIF_id","observer_id",
        "observation_id","pollinator_species","plant_species",
        "original_plant_source_id","source_family")
    if any('"'+x+'"' in body for x in forbidden):
        raise ValueError("Original rights-safe source island/plant/guild ID escaped")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+"\n")
    global_info=next(x for x in result[
        "original_global_guild_co_observation_by_region_and_geology"]
        if x["region"]=="GLOBAL_ALL" and x["geological_origin"]=="all")
    print(json.dumps(dict(
        status=result["status"],
        n_global_original_flora=result["n_original_global_plant_flora_islands"],
        n_islands_5_guilds_all_evaluable=result[
            "original_GBIF_evaluable_all_five_guilds"],
        n_islands_with_two_recorded_guilds=global_info[
            "n_at_least_two_recorded_guilds"],
        n_islands_with_any_qualified_guild_nondetection=global_info[
            "n_any_source_qualified_guild_nondetection"],
        no_historical_guild_shift_inferred=True,
    ),sort_keys=True))


if __name__=="__main__":
    main()
