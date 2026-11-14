"""Original FCP photographs × candidate CoL XR synonym corroboration audit.

This restores frozen 4,475-species genuine measured-photo FCP source and joins
with protected original island/GIFT native status + Roeble/Zizka geology.
Only geographic/regional opportunity counts leave original protected CI.
FCP white/D/W are NOT local-island flower observations or adaptive phenotypes.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import pandas as pd

from fd.fcp_endemic_colxr_synonym_gate import original_endemic_fcp_taxonomy_sensitivity

FCP_LEDGER_SHA256="23d3daadc9799c2e7b5867d54ccc1ed98369ba79a81b72dc7f831e5197e63c49"
COLXR_SHA256="079d9c27b0ca3aa51f0c7310f43ad4901b07a319d6063eab5a7d62b291bdd9e3"


def validate_bytes(p:Path,sha:str):
    h=hashlib.sha256()
    with p.open("rb") as r:
        for chunk in iter(lambda:r.read(1<<20),b""):
            h.update(chunk)
    if h.hexdigest()!=sha:
        raise ValueError("Original FCP 20261021 species photographic archive changed")


def run(root:Path,geology:Path,fcp_path:Path,colxr:Path)->dict:
    validate_bytes(fcp_path,FCP_LEDGER_SHA256)
    validate_bytes(colxr,COLXR_SHA256)
    fcp=pd.read_csv(fcp_path,low_memory=False)
    if fcp.species.isna().any() or fcp.species.duplicated().any() or len(fcp)!=4475:
        raise ValueError("FCP original 4475 accepted name roster changed")
    crosswalk=pd.read_csv(colxr,dtype=str).fillna("")
    status=pd.read_csv(
        root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        dtype=str,usecols=["island_id","accepted_species","origin_status",
                           "endemic_status","floristic_status"])
    geo=pd.read_csv(geology,dtype=str,usecols=[
        "island_id","geological_origin","analysis_regime","spatial_block"])
    typed=geo.loc[geo.geological_origin.isin(("continental","oceanic"))]
    if (typed.island_id.nunique()!=434 or
        typed.geological_origin.value_counts().get("continental",0)!=176 or
        typed.geological_origin.value_counts().get("oceanic",0)!=258):
        raise ValueError("Independently classified 176/258 original geology sample drift")
    result=original_endemic_fcp_taxonomy_sensitivity(
        set(fcp.species),status,geo,crosswalk)
    result["FCP_original_ledger_verified_sha256"]=FCP_LEDGER_SHA256
    result["CoL_XR_original_candidate_ledger_verified_sha256"]=COLXR_SHA256
    if result["n_original_geologically_typed_islands"]!=434:
        raise ValueError("Original ecology-geology intersection no longer stable")
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--geology-crosswalk",type=Path,required=True)
    p.add_argument("--fcp-ledger",type=Path,required=True)
    p.add_argument("--colxr-crosswalk",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    data=run(a.original_root,a.geology_crosswalk,a.fcp_ledger,a.colxr_crosswalk)
    text=json.dumps(data,sort_keys=True,indent=2,allow_nan=False)
    forbidden=(
        "accepted_species","island_id","endemic_status","origin_status",
        "floristic_status","entity_ID","work_species","source_signature",
        "spatial_block","latitude","longitude","inat_taxon_id",
        "submitted_name","candidate_accepted_name","candidate_accepted_key")
    if any('"'+x+'"' in text for x in forbidden):
        raise ValueError("Private original taxa, photo identifiers or island leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(text+"\n")
    print(json.dumps({
        "status":data["status"],
        "n_original_typed_islands":data["n_original_geologically_typed_islands"],
        "n_source_FCP_4475":data["original_taxonomy_candidate_audit"]["n_FCP_source_species"],
        "n_synonym_candidate_source_total":data["original_taxonomy_candidate_audit"]["n_one_to_one_original_photo_synonym_candidates"],
        "n_original_regional_status_panels":len(data["geological_native_endemic_candidate_overlap"]),
        "n_unreviewed_tropical_oceanic_endemic_candidates":next((
            q["n_one_to_one_CoL_synonym_candidate_FCP_species"]
            for q in data["geological_native_endemic_candidate_overlap"]
            if q["region"]=="tropical" and q["geological_origin"]=="oceanic"
            and q["verified_plant_origin_status"]=="endemic"),None)
    },sort_keys=True))

if __name__=="__main__":main()
