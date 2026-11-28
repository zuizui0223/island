#!/usr/bin/env python3
"""Outcome-motivated H2 link for the tropical paired-axis H1 signal.
No new source data. No causal mediation. Anonymous aggregates only.
"""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from chapter1_flower_axis_decoupling import LOCK,digest,trait_pairs,make_frame,rates,design,estimate,METRICS,holm,SEED

AUDIT_SHA="4890da7df53b2e240c4b1493ef93c0814a68435a5e1670e3c0e7c3e6717b90e7"
REPRO={
 "self_incompatibility":({"SC"},{"SI"}),
 "mating_system":({"predominantly_selfing","obligate_selfing"},{"predominantly_outcrossing"}),
 "autonomous_selfing_capacity":({"autonomous","delayed"},{"absent"}),
}
REGION="tropical"

def build_scores(audit):
    # Exact H2 selfing_core scoring rule (equal weights, min2/complete3).
    q=audit.loc[audit.resolved_for_primary.fillna("").astype(str).str.lower().isin({"true","1","yes"})&audit.trait_name.isin(REPRO)].copy()
    if q.duplicated(["accepted_species","trait_name"]).any():raise ValueError("duplicate state audit key")
    def signed(r):
        tok={x.strip() for x in str(r.canonical_signature).split("|") if x.strip()}
        plus,minus=REPRO[r.trait_name]
        return 1. if tok and tok<=plus else -1. if tok and tok<=minus else np.nan
    q["score"]=[signed(r) for r in q.itertuples()]
    w=q.pivot(index="accepted_species",columns="trait_name",values="score").reindex(columns=list(REPRO))
    n=w.notna().sum(axis=1)
    w["min2"]=w.mean(axis=1).where(n>=2)
    w["complete3"]=w[list(REPRO)].mean(axis=1).where(n==3)
    return w

def aggregate(f,s,key):
    x=f[["island_id","accepted_species"]].merge(s[[key]].dropna(),left_on="accepted_species",right_index=True,validate="many_to_one")
    return x.groupby("island_id")[key].agg(["mean","count"]).rename(columns={"mean":"score","count":"score_n"})

def one_panel(d,y,score,key,roster,scope,nmin):
    z=d.copy()
    z["y"]=y
    z=z.merge(score,left_on="island_id",right_index=True,how="left",validate="one_to_one")
    z=z.loc[z.score_n.ge(nmin)].dropna(subset=["score","y"])
    b=z.spatial_block.to_numpy(str)
    out={"scope":scope,"score":key,"roster":roster,"min_scored":nmin,"n":len(z),"blocks":len(set(b)),"status":"insufficient_support","p":1.}
    if len(z)<50 or len(set(b))<10:return out
    X,j=design(z)
    if X is None or z.score.std()<=1e-12:return out
    s=z.score.to_numpy(float);Y=z.y.to_numpy(float)[:,None]
    X2=np.column_stack([X,s-s.mean()])
    if np.linalg.matrix_rank(X2)<X2.shape[1]:return out|{"status":"rank_failure"}
    a=estimate(X,Y,b,j);c=estimate(X2,Y,b,j)
    if a is None or c is None:return out
    # Direct conditional difference on identical island mask.
    f=estimate(X2,Y,b,X2.shape[1]-1)
    return out|{"status":c["0"]["status"],"p":c["0"]["p"],
        "unadjusted_same_islands":a["0"],"adjusted_same_islands":c["0"],"reproduction_coefficient":f["0"] if f else None,
        "median_scored_species":float(z.score_n.median()),"signed_slope_change":c["0"]["beta"]-a["0"]["beta"]}

def self_test():
    a=pd.DataFrame([{"accepted_species":"A alba","trait_name":"self_incompatibility","canonical_signature":"SC","resolved_for_primary":True},
      {"accepted_species":"A alba","trait_name":"mating_system","canonical_signature":"predominantly_outcrossing","resolved_for_primary":True},
      {"accepted_species":"B alba","trait_name":"self_incompatibility","canonical_signature":"SC|SI","resolved_for_primary":True}])
    w=build_scores(a);assert w.loc["A alba","min2"]==0 and np.isnan(w.loc["A alba","complete3"]) and np.isnan(w.loc["B alba","min2"])
    print("H2_DECOUPLING_SCORE_TEST_PASS")

def main():
    p=argparse.ArgumentParser()
    for k in ("flora","ledger","geography","state-audit","out"):p.add_argument("--"+k,type=Path)
    p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    for k,h in LOCK.items():
        if digest(getattr(a,k))!=h:raise ValueError("source hash mismatch "+k)
    if digest(a.state_audit)!=AUDIT_SHA:raise ValueError("frozen H2 state audit mismatch")
    f=pd.read_csv(a.flora,dtype=str,usecols=["island_id","accepted_species","origin_status"])
    ledger=pd.read_csv(a.ledger,dtype=str)
    geo=pd.read_csv(a.geography,dtype={"island_id":str,"spatial_block":str})
    audit=pd.read_csv(a.state_audit,dtype=str)
    scores=build_scores(audit)
    # SC is a distinct atomic reproduction trait; do not infer autonomy.
    scores["SC"]=scores.self_incompatibility
    scores["autonomy"]=scores.autonomous_selfing_capacity
    rows=[];supports=[]
    for scope in ("all_recorded","high_only","literal_native"):
        v=trait_pairs(ledger,scope=="high_only")
        d,M,C,g=make_frame(f,v,geo,scope)
        d=d.loc[d.analysis_regime.eq(REGION)&d.carrier_n.ge(10)&d.deep_n.ge(5)&d.shallow_n.ge(5)].sort_values("island_id")
        y=rates(np.asarray(M[d.mi.to_numpy(int)]@C))[:,METRICS.index("conditional_difference")]
        pool=f if scope!="literal_native" else f.loc[f.origin_status.eq("native")]
        pair=pool.loc[pool.accepted_species.isin(v.index)]
        supports.append({"scope":scope,"morphology_islands":len(d),"morphology_pair_species":len(v),"min2_species":int(scores.min2.notna().sum()),"complete3_species":int(scores.complete3.notna().sum()),"same_pair_min2_species":int(scores.loc[scores.index.intersection(v.index),"min2"].notna().sum())})
        for roster,ff in (("whole_flora",pool),("pair_carriers",pair)):
            for key in ("min2","complete3","SC","autonomy"):
                # Fixed minimum 10 score carriers for every panel; same-pair
                # score uses a subset, not all pair-carriers.
                rows.append(one_panel(d,y,aggregate(ff,scores,key),key,roster,scope,10))
    assert len(rows)==24
    for scope in ("all_recorded","high_only","literal_native"):holm([r for r in rows if r["scope"]==scope])
    result={"owner":"island","status":"exploratory_H2_extension_not_main_promotion","state_audit_sha256":AUDIT_SHA,"support":supports,"panels":rows,
      "primary":"all_recorded / whole_flora / min2, original H2 selfing score; Holm over 8 fixed rows per scope",
      "warning":"Outcome-motivated H2 conditional diagnostic on H1 paired morphology outcome. Same-island before/after estimates; pair-carrier score is only known subset. High-only refers to morphology pair grade, not reproduction score grade (all use frozen Direct state audit). Not causal mediation, not independent evidence, not calibrated grade robustness. Unsupported rows retained."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("H2_DECOUPLING_SUPPORT "+json.dumps(supports,allow_nan=False),flush=True)
    for row in rows:print("H2_DECOUPLING_PANEL "+json.dumps(row,allow_nan=False),flush=True)
    print("H2_DECOUPLING_COMPLETE 24",flush=True)

if __name__=="__main__":main()
