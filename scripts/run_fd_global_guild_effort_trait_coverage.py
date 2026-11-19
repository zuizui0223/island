"""Original full 4,505 GIFT flora and NEE GBIF BACKGROUND EFFORT documentation.

Routinely repeatable source provenance audit. None of the original
island/plant/observer/GIFT source/GBIF taxon identities leave original CI.
Missing GBIF background rows are unknown, NOT zero guild members or
confirmed adequate non-detection. This study estimates documentation
covariation, not true floral visitor or historical guild loss.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import pandas as pd
from fd.global_guild_effort_trait_coverage import (
    original_global_flora_guild_effort_audit
)


def digest(path:Path)->str:
    q=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            q.update(chunk)
    return q.hexdigest()


def run(original:Path,geology:Path,channels:Path)->dict:
    ledger=original/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora=original/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    if digest(ledger)!="9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745":
        raise ValueError("Original strictly Direct source provenance mismatch")
    if digest(flora)!="6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a":
        raise ValueError("Original full 4505 GIFT flora source provenance mismatch")
    fl=pd.read_csv(flora,dtype=str,
        usecols=["island_id","accepted_species","origin_status"])
    if fl.island_id.nunique()!=4505:
        raise ValueError("Full botanical flora original population drift")
    tt=pd.read_csv(ledger,dtype=str,
        usecols=["accepted_species","trait_name","normalized_value",
                 "evidence_scope","evidence_quality"])
    geo=pd.read_csv(geology,dtype=str,
        usecols=["island_id","analysis_regime",
                 "geological_origin","spatial_block"])
    if (geo.island_id.nunique()<8200 or
        geo.geological_origin.eq("continental").sum()!=176 or
        geo.geological_origin.eq("oceanic").sum()!=258):
        raise ValueError("Independent Roeble/Zizka original geological classification changed")
    co=pd.read_csv(channels,dtype=str,
        usecols=["island_id","channel_id","channel_state","source_state",
                 "background_record_count","background_spatial_units",
                 "background_temporal_units","distinct_dataset_count"])
    out=original_global_flora_guild_effort_audit(fl,tt,geo,co)
    out["original_source_frozen_plant_ledger_SHA256"]=digest(ledger)
    out["original_source_qualified_GBIF_background_SHA256"]=digest(channels)
    out["original_source_rights_protected_plant_and_animal_rows_not_exported"]=True
    return out


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--original-root",required=True,type=Path)
    a.add_argument("--geology-crosswalk",required=True,type=Path)
    a.add_argument("--channel-evidence-csv",required=True,type=Path)
    a.add_argument("--output",required=True,type=Path)
    q=a.parse_args()
    obj=run(q.original_root,q.geology_crosswalk,q.channel_evidence_csv)
    body=json.dumps(obj,sort_keys=True,indent=2,allow_nan=False)
    forbidden=(
        "island_id","accepted_species","entity_ID","work_species",
        "source_signature","source_region_id","spatial_block",
        "latitude","longitude","gbif_id","taxon_id",
        "observer_id","source_genus","source_family",
        "evidence_id","photo_id")
    if any('"'+key+'"' in body for key in forbidden):
        raise ValueError("Private original GIFT/GBIF source identifers leaked")
    q.output.parent.mkdir(parents=True,exist_ok=True)
    q.output.write_text(body+"\n")
    x=next(v for v in obj[
        "original_plant_flora_and_GBIF_background_effort_by_group"]
        if v["region"]=="GLOBAL_ALL" and v["geological_origin"]=="all"
        and v["guild"]=="bombus")
    print(json.dumps({
        "status":obj["status"],
        "n_original_flora_islands":obj["n_original_botanical_islands"],
        "n_full_original_animal_state_rows":obj["n_original_GBIF_island_guild_cells"],
        "n_Bombus_detected":x["n_GBIF_guild_recorded"],
        "n_Bombus_adequate_nondetected":x["n_guild_source_qualified_NONdetection"],
        "n_Bombus_background_record_effort_bins_with_both_states":x[
            "n_background_record_effort_bins_with_both_qualified_states"],
        "no_plant_or_guild_effectiveness_claim":True,
    },sort_keys=True))


if __name__=="__main__": main()
