"""Run fd source-provider coverage audit from original pinned GloBI v0.9 TSV.

Input target genera are read from the PRE-EXISTING outcome-independent source
genus breadth artifact. Only a privacy-safe JSON receipt is published.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd

from island_v2.fd_globi_provider_v09 import analyze_v09_gzip, VERSION, PIN
from island_v2.fd_globi_provider_family_transfer import family_transfer

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--interactions",required=True,type=Path)
    p.add_argument("--frozen-source-genus-csv",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--original-direct-ledger",type=Path)
    p.add_argument("--original-island-taxonomy",type=Path)
    p.add_argument("--chunksize",type=int,default=200_000)
    a=p.parse_args()
    d=pd.read_csv(a.frozen_source_genus_csv)
    needed={"genus","n_reference_channel_units","n_independent_references"}
    if not needed.issubset(d):raise ValueError("Missing frozen source-genus breadth contract")
    if d.genus.isna().any() or d.genus.duplicated().any():
        raise ValueError("Unresolved or duplicate focal plant genus")
    original, after, result=analyze_v09_gzip(
        a.interactions, set(d.genus.astype(str)), chunksize=a.chunksize
    )
    if bool(a.original_direct_ledger) != bool(a.original_island_taxonomy):
        raise ValueError("Need BOTH the original direct ledger and frozen GBIF taxonomy")
    if a.original_direct_ledger:
        for pth in (a.original_direct_ledger, a.original_island_taxonomy):
            if not pth.is_file():
                raise ValueError("Missing original source-backed morphology or taxonomy")
        original_ledger=pd.read_csv(a.original_direct_ledger,
           usecols=["accepted_species","trait_name","normalized_value",
                    "evidence_scope","evidence_quality"])
        original_taxonomy=pd.read_csv(a.original_island_taxonomy,
           usecols=["genus","family"])
        result["provider_omitted_family_transfer"]=family_transfer(
            original,after,original_ledger,d,original_taxonomy)

    result.update({
        "n_outcome_blind_frozen_original_genus_list":int(len(d)),
        "v09_full_TSV_source_version":VERSION,
        "v09_full_TSV_SHA256":PIN,
        "focal_island_traits_used":False,
        # Source universe was fixed before looking at morphology. The
        # prediction DOES use Direct morphology in its feature matrix;
        # avoid implying that no phenotype data entered the prediction.
        "source_universe_selected_without_flower_traits":True,
        "direct_flower_traits_used_for_family_prediction":
            bool(a.original_direct_ledger),
        "provider_family_prediction_rerun_completed":
            bool(result.get("provider_omitted_family_transfer",{}).get(
                "heldout_original_provider_prediction_completed",False)),
        "provider_heldout_effect_on_flower_morphology_NOT_TESTED":
            not result.get("provider_omitted_family_transfer",{}).get(
                "heldout_original_provider_prediction_completed",False),
        "original_archive_source_granularity_unresolved":True,
    })
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
