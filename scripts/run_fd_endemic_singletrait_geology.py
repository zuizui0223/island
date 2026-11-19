"""Original protected geology-native endemic floral single-trait evidence preflight.

This reads immutable original source-Direct trait ledger and separately
source-backed GIFT island endemic-status, not a synthetic floral profile.
No original accepted taxa, island IDs, endemic region label, coordinates or
donor source entities are published. Only cumulative geographic-origin,
region, strict endemic status and direct-trait evidence counts.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import pandas as pd
from fd.geological_endemic_trait_opportunity import original_endemic_trait_opportunity


def run(original:Path,crosswalk:Path):
    direct=pd.read_csv(
        original/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        dtype=str,usecols=["accepted_species","trait_name","normalized_value",
                           "evidence_scope","evidence_quality"])
    flora=pd.read_csv(
        original/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        dtype=str,usecols=["island_id","accepted_species","origin_status",
                           "endemic_status","floristic_status"])
    geo=pd.read_csv(crosswalk,dtype=str,usecols=[
        "island_id","geological_origin","analysis_regime","spatial_block"])
    typed=geo.loc[geo.geological_origin.isin(("continental","oceanic"))]
    if (not 430<=len(typed)<=450 or
        typed.geological_origin.value_counts().get("continental",0)<150 or
        typed.geological_origin.value_counts().get("oceanic",0)<200):
        raise ValueError("Independent original Roeble/Zizka geologically typed source drift")
    nearest=pd.read_csv(
        original/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    donors=pd.read_csv(
        original/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"])
    r=original_endemic_trait_opportunity(direct,flora,geo,nearest,donors)
    if r["n_original_geological_typed_islands"]!=len(typed):
        raise ValueError("Original geological/source island label drift")
    if r["n_geologically_typed_islands_with_source_flora"]<300:
        raise ValueError("Original published geology and flora overlap too sparse")
    return r


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--geology-crosswalk",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    data=run(a.original_root,a.geology_crosswalk)
    text=json.dumps(data,sort_keys=True,indent=2,allow_nan=False)
    private=(
       "island_id","accepted_species","origin_status","endemic_status",
       "floristic_status","entity_ID","work_species","source_signature",
       "spatial_block","source_family","source_genus","geological_record_id",
       "native_l3_codes","tdwg_l3_code","longitude","latitude")
    if any('"'+w+'"' in text for w in private):
        raise ValueError("Original taxon/island/endemic/source identity escaped")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(text+"\n")
    print(json.dumps({
        "status":data["status"],
        "n_geologically_typed_original_islands":data["n_original_geological_typed_islands"],
        "n_typed_source_flora_islands":data["n_geologically_typed_islands_with_source_flora"],
        "n_trait_scope_groups":len(data["strict_original_status_counts_by_region_geology_trait"]),
        "n_direct_species_by_trait":data["original_global_single_trait_Direct_species_counts"],
    },sort_keys=True))


if __name__=="__main__":
    main()
