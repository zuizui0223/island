"""Original-source paired morphology comparison, anonymized regional output.

Reuses source-only congener prediction rows; this is NOT out-of-archipelago
or provider heldout and cannot establish a causal pollinator mechanism.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from run_fd_original_source_only_joint_congener_prediction import (
    original_source_only_prediction_frame, digest, SHA,
)
from fd.origin import run_crosswalk

SECTORS=("actinomorphic_deep","tubular_bell_deep")
MODES=("geo5","climate10")
SCOPES=("all_recorded","strict_native","strict_native_nonendemic")
GATES=("any_congener","at_least_two_congeners")

def paired_aggregate(frame:pd.DataFrame,*,draws:int=999)->dict:
    if frame.island_id.nunique()!=4505 or len(frame)!=4505*3*3*2:
        raise ValueError("original botanical population reduced")
    results=[]
    for region in sorted(frame.regime.unique()):
        for scope in SCOPES:
            for gate in GATES:
                a=frame.loc[frame.regime.eq(region)&frame.scope.eq(scope)&frame.donor_gate.eq(gate)&frame.sector.isin(SECTORS)].copy()
                for mode in MODES:
                    n=a[mode+"__n_target_direct"].to_numpy(float)
                    donor=a[mode+"__n_source_matched_direct"].to_numpy(float)
                    b=a[mode+"__brier_pooled"].to_numpy(float)
                    g=a[mode+"__brier_genus"].to_numpy(float)
                    a[mode+"__gain"]=b-g
                    a[mode+"__supported"]=(n>=5)&(donor>=10)&np.isfinite(b)&np.isfinite(g)
                for mode in MODES:
                    valid=a.loc[a[mode+"__supported"]].copy()
                    p=valid.pivot(index="island_id",columns="sector",values=mode+"__gain").dropna()
                    original=a.drop_duplicates("island_id").set_index("island_id")
                    z=original.loc[p.index,["spatial_block","source_signature"]].copy()
                    z["delta"]=p[SECTORS[0]]-p[SECTORS[1]]
                    z=z.dropna(subset=["spatial_block","source_signature"])
                    block=z.groupby("spatial_block",sort=True).delta.mean()
                    source=z.groupby("source_signature",sort=True).delta.mean()
                    cell={"region":region,"scope":scope,"donor_gate":gate,"source_mode":mode,
                          "n_original_botanical_islands_in_region":int(a.island_id.nunique()),
                          "n_common_source_supported_islands":len(z),
                          "n_geographical_blocks":len(block),
                          "n_first_source_entities":len(source)}
                    if len(z)<30 or len(block)<10 or len(source)<10:
                        cell["status"]="NOT_IDENTIFIABLE_COMMON_ISLAND_SUPPORT"
                    else:
                        cell["status"]="DESCRIPTIVE_SOURCE_ONLY_COMMON_ISLAND_PAIRED"
                        for key,v in (("block",block),("source",source)):
                            arr=v.to_numpy(float)
                            rng=np.random.default_rng(20261117+(0 if key=="block" else 1))
                            boot=arr[rng.integers(len(arr),size=(draws,len(arr)))].mean(1)
                            cell[key+"_equal_delta_gain_radial_minus_tubular"]=float(arr.mean())
                            cell[key+"_bootstrap95"]=[float(x) for x in np.quantile(boot,[.025,.975])]
                    results.append(cell)
    return {"status":"ORIGINAL_4505_COMMON_ISLAND_PAIRED_CONGENER_GAIN",
            "n_original_botanical_islands":4505,
            "n_panels":len(results),"panels":results,
            "limits":{"source_only_distinct_accepted_species":True,
                      "archipelago_and_provider_holdout_not_done":True,
                      "same_island_but_not_necessarily_same_target_species":True,
                      "source_candidate_not_arrival":True,
                      "not_causal_guild_or_adaptation":True}}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--geology-crosswalk",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    x=p.parse_args()
    root=x.original_root
    ledger=root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora=root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    for path in (ledger,flora):
        if digest(path)!=SHA[path.name]:raise ValueError("original botanical SHA drift")
    direct=pd.read_csv(ledger,dtype=str,usecols=["accepted_species","trait_name","normalized_value","evidence_scope","evidence_quality"])
    status=pd.read_csv(flora,dtype=str,usecols=["island_id","accepted_species","origin_status","endemic_status","floristic_status"])
    geo=pd.read_csv(x.geology_crosswalk,dtype=str,usecols=["island_id","analysis_regime","spatial_block"])
    sources=pd.read_csv(root/"fixed/source/result/island_source_assignments.csv.gz",dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(root/"fixed/source/source/gift_native_mainland_flora.csv.gz",dtype=str,usecols=["entity_ID","work_species"])
    rows=original_source_only_prediction_frame(direct,status,geo,sources,mainland)
    out=paired_aggregate(rows)
    raw=json.dumps(out,sort_keys=True,indent=2,allow_nan=False)
    for name in ("island_id","accepted_species","work_species","entity_ID","source_signature","spatial_block","latitude","longitude"):
        if '"'+name+'"' in raw:raise ValueError("source identifiers must remain private")
    x.output.parent.mkdir(parents=True,exist_ok=True)
    x.output.write_text(raw+"\n")
    print(json.dumps({"status":out["status"],"n_panels":out["n_panels"],
                      "supported_panels":sum(p["status"].startswith("DESCRIPTIVE") for p in out["panels"])}))

if __name__=="__main__":main()
