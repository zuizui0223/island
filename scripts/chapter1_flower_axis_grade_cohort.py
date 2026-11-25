#!/usr/bin/env python3
"""Outcome-motivated coverage/grade decomposition of the tropical H1 extension."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from chapter1_flower_axis_decoupling import (LOCK, digest, trait_pairs, make_frame, rates, design, estimate, METRICS, bootstrap, SEED)

def main():
    p=argparse.ArgumentParser()
    for k in ("flora","ledger","geography","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    for k,h in LOCK.items():
        if digest(getattr(a,k))!=h:raise ValueError("source hash mismatch "+k)
    f=pd.read_csv(a.flora,dtype=str,usecols=["island_id","accepted_species","origin_status"])
    ledger=pd.read_csv(a.ledger,dtype=str)
    geo=pd.read_csv(a.geography,dtype={"island_id":str,"spatial_block":str})
    def setup(high):
        v=trait_pairs(ledger,high)
        d,M,C,g=make_frame(f,v,geo,"high_only" if high else "all_recorded")
        d=d.loc[d.analysis_regime.eq("tropical")&d.carrier_n.ge(10)&d.deep_n.ge(5)&d.shallow_n.ge(5)].sort_values("island_id")
        return d,M,C,g
    all_d,all_m,all_c,all_g=setup(False)
    hi_d,hi_m,hi_c,hi_g=setup(True)
    assert len(all_d)==257 and len(hi_d)==62
    assert set(hi_d.island_id)<=set(all_d.island_id)
    common=all_d.set_index("island_id").loc[hi_d.island_id].reset_index()
    assert common.island_id.to_list()==hi_d.island_id.to_list()
    m=METRICS.index("conditional_difference")
    y0=rates(np.asarray(all_m[all_d.mi.to_numpy(int)]@all_c))[:,m]
    y1=rates(np.asarray(all_m[common.mi.to_numpy(int)]@all_c))[:,m]
    y2=rates(np.asarray(hi_m[hi_d.mi.to_numpy(int)]@hi_c))[:,m]
    x0,j=design(all_d);x1,j1=design(common);x2,j2=design(hi_d)
    if any(x is None for x in (x0,x1,x2)):raise ValueError("constant design")
    b0=all_d.spatial_block.to_numpy(str);b1=common.spatial_block.to_numpy(str);b2=hi_d.spatial_block.to_numpy(str)
    assert np.array_equal(b1,b2)
    v0=estimate(x0,y0[:,None],b0,j)
    v1=estimate(x1,np.column_stack([y1,y2,y2-y1]),b1,j1)
    v2=estimate(x2,y2[:,None],b2,j2)
    if any(v is None for v in (v0,v1,v2)):raise ValueError("matched-cohort fit not supported")
    raw=v0["0"]["beta"];same=v1["0"]["beta"];high_common=v1["1"]["beta"];high_original=v2["0"]["beta"]
    assert abs(raw-0.02762020672544175)<1e-10 and abs(high_original-(-0.0007223441697008538))<1e-10
    effects={"island_selection_and_design_shift":same-raw,"carrier_change_same_islands_same_covariates":high_common-same,"coverage_covariate_shift":high_original-high_common}
    assert np.isclose(sum(effects.values()),high_original-raw)
    result={"scope":"tropical conditional symmetry difference","full_n":len(all_d),"matched_n":len(common),"matched_blocks":len(set(b1)),
      "all_full":v0["0"],"all_common_islands":v1["0"],"high_common_islands_all_covariates":v1["1"],
      "high_minus_all_matched_paired":v1["2"],"high_original_covariates":v2["0"],"signed_change_decomposition":effects,
      "warning":"Post-outcome matched-cohort diagnostic. High-only vs High/Medium compares different documented species, not repeat measurements; composition, documentation, and geography can differ. All intervals pointwise conditional; no causal measurement-bias conclusion."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("FLORAL_DECOUPLING_MATCHED_GRADE "+json.dumps(result,allow_nan=False),flush=True)
if __name__=="__main__":main()
