"""Original rights-safe FCP species-wide photo evidence × oceanic island endemics.

This restores frozen 4,475-species genuine measured-photo FCP source and joins
with protected original island/GIFT native status + Roeble/Zizka geology.
Only geographic/regional opportunity counts leave original protected CI.
FCP white/D/W are NOT local-island flower observations or adaptive phenotypes.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import pandas as pd

from fd.fcp_geological_endemic_opportunity import audit_fcp_endemic_geology_opportunity

FCP_LEDGER_SHA256="23d3daadc9799c2e7b5867d54ccc1ed98369ba79a81b72dc7f831e5197e63c49"


def validate_bytes(p:Path,sha:str):
    h=hashlib.sha256()
    with p.open("rb") as r:
        for chunk in iter(lambda:r.read(1<<20),b""):
            h.update(chunk)
    if h.hexdigest()!=sha:
        raise ValueError("Original FCP 20261021 species photographic archive changed")


def run(root:Path,geology:Path,fcp_path:Path)->dict:
    validate_bytes(fcp_path,FCP_LEDGER_SHA256)
    fcp=pd.read_csv(fcp_path,low_memory=False)
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
    src=pd.read_csv(
        root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    donor=pd.read_csv(
        root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"])
    result=audit_fcp_endemic_geology_opportunity(fcp,status,geo,src,donor)
    result["FCP_original_ledger_verified_sha256"]=FCP_LEDGER_SHA256
    if result["n_original_geologically_typed_islands"]!=434:
        raise ValueError("Original ecology-geology intersection no longer stable")
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--geology-crosswalk",type=Path,required=True)
    p.add_argument("--fcp-ledger",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    data=run(a.original_root,a.geology_crosswalk,a.fcp_ledger)
    text=json.dumps(data,sort_keys=True,indent=2,allow_nan=False)
    forbidden=(
        "accepted_species","island_id","endemic_status","origin_status",
        "floristic_status","entity_ID","work_species","source_signature",
        "spatial_block","latitude","longitude","species","inat_taxon_id")
    if any('"'+x+'"' in text for x in forbidden):
        raise ValueError("Private original taxa, photo identifiers or island leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(text+"\n")
    print(json.dumps({
        "status":data["status"],
        "n_original_typed_islands":data["n_original_geologically_typed_islands"],
        "n_photo_species_original":data["n_FCP_species_original"],
        "n_cross_geo_endemic_photo_panels":len(data["original_regional_photo_opportunity"]),
        "n_fcp_endemic_distinct_species_oceanic_tropical":next((
            q["n_original_exact_FCP_distinct_species"]
            for q in data["original_regional_photo_opportunity"]
            if q["region"]=="tropical" and q["geological_origin"]=="oceanic"
            and q["native_endemic_status"]=="endemic"),None)
    },sort_keys=True))

if __name__=="__main__":main()
