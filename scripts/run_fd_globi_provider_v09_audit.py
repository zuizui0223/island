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

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--interactions",required=True,type=Path)
    p.add_argument("--frozen-source-genus-csv",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--chunksize",type=int,default=200_000)
    a=p.parse_args()
    d=pd.read_csv(a.frozen_source_genus_csv)
    needed={"genus","n_reference_channel_units","n_independent_references"}
    if not needed.issubset(d):raise ValueError("Missing frozen source-genus breadth contract")
    if d.genus.isna().any() or d.genus.duplicated().any():
        raise ValueError("Unresolved or duplicate focal plant genus")
    _,_,result=analyze_v09_gzip(a.interactions,set(d.genus.astype(str)),chunksize=a.chunksize)
    result.update({
        "n_outcome_blind_frozen_original_genus_list":int(len(d)),
        "v09_full_TSV_source_version":VERSION,
        "v09_full_TSV_SHA256":PIN,
        "focal_island_traits_used":False,
        "provider_heldout_effect_on_flower_morphology_NOT_TESTED":True
    })
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
