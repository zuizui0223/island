#!/usr/bin/env python3
"""Post-hoc exact atomic support sensitivity for broad-flora H1 and H2.

Reconstruct the frozen axis-Direct baseline first. Retain original informative
state signatures only if the same accepted species/trait has the exact state
set in the original atomic ledger, resolved High/Medium. No alternate state,
quality upgrade, source-policy change, refitting WCVP, or outcome-selected
subset. All original model functions/configuration retained. Aggregate output.
"""
from __future__ import annotations
import argparse, hashlib, io, json, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t
import yaml
from island_v2.chapter1_h1_traitwise import run_scope, fit_trait
from island_v2.chapter1_all_data_probability import build_broad_counts, _prepare
from island_v2.chapter1_h2_complete_selfing_sensitivity import (
    build_species_score, aggregate_island_score, fit_models,
)
from island_v2.chapter1_v14_h2_decomposition import run_decomposition, _bh

ROOT=Path(__file__).resolve().parents[3]
PINS={
    "analysis":"b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999",
    "ledger":"5e9c887e093c42eb720f3abd2dc85701ce123fb7fe888c916dc4098ac898ccdd",
    "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
    "audit":"4890da7df53b2e240c4b1493ef93c0814a68435a5e1670e3c0e7c3e6717b90e7",
    "atomic":"7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577",
    "geography":"08218617b5880153986012e79817de65a74dd4e44b5a65e8a537dab8cb2cd0a7",
    "h1_config_historical_crlf":"b20a672f9a23463f21d1be6d713c5550ced69cc303d620f49fb1f54871d1e9a7",
}
def sha(b):return hashlib.sha256(b).hexdigest()
def assert_hash(b,k):
    if sha(b)!=PINS[k]:raise ValueError("source digest mismatch: "+k)
def read_zip(z,p,k,**kwargs):
    b=z.read(p);assert_hash(b,k)
    return pd.read_csv(io.BytesIO(b),compression="gzip",**kwargs)
def truthy(s):return s.astype(str).str.lower().isin(["true","1","yes"])
def strict_mask(a,d):
    key=["accepted_species","trait_name"]
    assert not a.duplicated(key).any() and not d.duplicated(key).any()
    good=d[d.quality.isin(["high","medium"]) & d.resolution_status.eq("resolved")].copy()
    good["atomic_signature"]=good.state_set.map(lambda x:"|".join(sorted(set(json.loads(x)))))
    x=a.merge(good[key+["atomic_signature"]],on=key,how="left",validate="one_to_one",sort=False)
    # Use exactly the original state-audit signature; preserve all original
    # columns and row order, including multistate signatures if eligible.
    keep=truthy(x.resolved_for_primary)&x.canonical_signature.eq(x.atomic_signature)
    assert len(x)==len(a)
    return a.loc[keep.to_numpy()].copy(),x
def finite_t(frame,prefix="distance"):
    beta=pd.to_numeric(frame[prefix+"_estimate"],errors="coerce")
    se=pd.to_numeric(frame[prefix+"_se"],errors="coerce")
    df=pd.to_numeric(frame.n_clusters,errors="coerce")-1
    return pd.Series(2*t.sf(abs(beta/se),df),index=frame.index).where((df>0)&(se>0)&frame.status.eq("fit"))
def bh_fixed(p):
    # Retain all planned family slots, count untestable as 1 internally only.
    return _bh(p.fillna(1)).where(p.notna())
def h2(f,a,c,h1cfg,sc,h2cfg,scope):
    min_n=int(sc["score_definition"]["minimum_informative_traits"])
    scores=[];sizes={}
    for name in ("selfing_core","generalized_accessible"):
        s=build_species_score(a,sc["syndromes"][name],minimum_informative_traits=min_n,require_all=False)
        sizes[name]=int(s.score.notna().sum())
        w=aggregate_island_score(f,s,score_name="syndrome_score")
        w["stratum"]="all_observed";w["syndrome"]=name;scores.append(w)
    cnt=build_broad_counts(f,a,h1cfg)
    out,_=run_decomposition(cnt,pd.concat(scores,ignore_index=True),c,h2cfg,evidence_scope=scope)
    out["distance_p_finite_cluster_t"]=finite_t(out)
    primary=out.response.isin(["plain_colour","generalized_accessible"])
    assert int(primary.sum())==8
    out["primary_H2b_q_finite_cluster_t"]=np.nan
    out.loc[primary,"primary_H2b_q_finite_cluster_t"]=bh_fixed(out.loc[primary,"distance_p_finite_cluster_t"])
    return out,sizes
def comp3(f,a,c,sc,h2cfg):
    out,s=fit_models(f,a,c,sc,h2cfg)
    # Original complete3 table uses normal-reference p and BH4. Add finite t
    # BH4 as a distinct stricter sensitivity, not silently replace the old rule.
    out["distance_p_finite_cluster_t"]=finite_t(out)
    out["distance_q_finite_cluster_t"]=np.nan
    for v,g in out.groupby("variant"):
        assert len(g)==4
        out.loc[g.index,"distance_q_finite_cluster_t"]=bh_fixed(g.distance_p_finite_cluster_t)
    return out,s
def gate(new,old,keys,cols,tols,exact):
    m=new.merge(old,on=keys,suffixes=("_new","_old"),validate="one_to_one")
    if len(m)!=len(new) or len(m)!=len(old):raise ValueError("baseline row roster differs")
    checks={}
    for k in exact:
        checks[k]=bool(m[k+"_new"].eq(m[k+"_old"]).all())
    for k,tol in zip(cols,tols,strict=True):
        x=pd.to_numeric(m[k+"_new"],errors="coerce")
        y=pd.to_numeric(m[k+"_old"],errors="coerce")
        same_missing=x.isna().eq(y.isna()).all()
        delta=(x-y).abs().max()
        checks[k]={"max_abs_diff":float(delta) if pd.notna(delta) else None,"tolerance":tol,
                   "pass":bool(same_missing and (pd.isna(delta) or delta<=tol))}
    return {"rows":len(m),"checks":checks,"pass":all(v if isinstance(v,bool) else v["pass"] for v in checks.values())}
def safe_rows(frame):
    return json.loads(frame.to_json(orient="records",double_precision=15))
def emit(path,obj):
    path.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n")
def run(analysis,ledger,out):
    for p,k in [(analysis,"analysis"),(ledger,"ledger")]:assert_hash(p.read_bytes(),k)
    with zipfile.ZipFile(analysis) as z,zipfile.ZipFile(ledger) as d:
        f=read_zip(z,"fixed/canonical/input/chapter1_status_flora.csv.gz","flora",dtype={"island_id":str})
        a=read_zip(z,"snapshot/chapter1_trait_state_audit_direct.csv.gz","audit",dtype=str,keep_default_na=False)
        atomic=read_zip(d,"direct_species_trait_ledger.csv.gz","atomic",dtype=str,keep_default_na=False)
    covpath=ROOT/"results/geography_20260924/corrected_geography_covariates.csv"
    assert_hash(covpath.read_bytes(),"geography");cov=pd.read_csv(covpath,dtype={"island_id":str})
    cpath=ROOT/"config/chapter1_h1_final_traitwise_t_20261004.yml"
    raw=cpath.read_bytes().replace(b"\r\n",b"\n")
    assert_hash(raw.replace(b"\n",b"\r\n"),"h1_config_historical_crlf")
    cfg=yaml.safe_load(raw)
    scpath=ROOT/"config/chapter1_pr138_pollination_syndromes.yml"
    h2path=ROOT/"config/chapter1_v14_h2_decomposition.yml"
    sc=yaml.safe_load(scpath.read_bytes());hcfg=yaml.safe_load(h2path.read_bytes())
    b,x=strict_mask(a,atomic)
    out.mkdir(parents=True,exist_ok=True)
    print("ATOMIC_SCOPE_ROSTER "+json.dumps({"original_audit_rows":len(a),"atomic_exact_resolved_rows":len(b),
          "kept_by_trait":{k:int(v) for k,v in b.trait_name.value_counts().items()}}),flush=True)
    # Gate original baseline BEFORE interpreting sensitivity.
    h1old=pd.read_csv(ROOT/"results/h1_final_traitwise_t_20261004/direct_broad.csv")
    h2old=pd.read_csv(ROOT/"results/h1_final_directional_20261003/h2_direct_finite_cluster_audit.csv")
    c3old=pd.read_csv(ROOT/"results/h2_complete_selfing_corrected_20261003/h2_complete_selfing_sensitivity.csv")
    h1b=run_scope(f,a,cov,cfg,"all_observed")
    h2b,n0=h2(f,a,cov,cfg,sc,hcfg,"axis_direct")
    c3b,s0=comp3(f,a,cov,sc,hcfg)
    h1b.to_csv(out/"h1_axis_replay.csv",index=False)
    h2b.to_csv(out/"h2_axis_replay.csv",index=False)
    c3b.to_csv(out/"h2_complete_axis_replay.csv",index=False)
    # Pinned original Python/scipy version differs across platforms. Small
    # predeclared tolerance for numerical optimizer; exact row/island/trial
    # and cluster counts. OLS comparisons are separately nearly exact.
    g1=gate(h1b,h1old,["context","outcome"],["estimate","se","p_two_sided"],[.003,.003,.02],["status","n_islands","n_species_trials","n_clusters"])
    g2=gate(h2b,h2old,["context","response"],["distance_estimate","distance_se","distance_p_finite_cluster_t","primary_H2b_q_finite_cluster_t"],[.003,.003,.02,.03],["status","n_islands","n_clusters"])
    g2ols=gate(h2b[h2b.model_family.eq("clustered_ols")],h2old[h2old.model_family.eq("clustered_ols")],["context","response"],["distance_estimate","distance_se"],[1e-9,1e-9],["status","n_islands","n_clusters"])
    g3=gate(c3b,c3old,["context","variant"],["distance_estimate","distance_se","distance_p","distance_q_within_variant"],[1e-9,1e-9,1e-8,1e-8],["status","n_islands","n_clusters"])
    gates={"h1":g1,"h2":g2,"h2_ols":g2ols,"complete3":g3}
    emit(out/"BASELINE_REPLAY_GATES.json",gates)
    print("ATOMIC_BASELINE_GATES "+json.dumps(gates,allow_nan=False),flush=True)
    if not all(g["pass"] for g in gates.values()):raise SystemExit("Frozen baseline replay failed; do not interpret strict model")
    h1s=run_scope(f,b,cov,cfg,"all_observed")
    h2s,n1=h2(f,b,cov,cfg,sc,hcfg,"atomic_exact_direct")
    c3s,s1=comp3(f,b,cov,sc,hcfg)
    h1s.to_csv(out/"h1_atomic.csv",index=False);h2s.to_csv(out/"h2_atomic.csv",index=False);c3s.to_csv(out/"h2_complete_atomic.csv",index=False)
    # Same-island H1 comparison separates carrier/trait change from island
    # inclusion change. Both regressions restandardize within common mask.
    p0=_prepare(build_broad_counts(f,a,cfg),cov,cfg);p1=_prepare(build_broad_counts(f,b,cfg),cov,cfg)
    common=[]
    for context in cfg["contexts"]:
        for outcome in cfg["model_outcomes"]:
            q0=p0[p0[cfg["context_column"]].eq(context)&p0.outcome.eq(outcome)]
            q1=p1[p1[cfg["context_column"]].eq(context)&p1.outcome.eq(outcome)]
            ids=set(q0.island_id)&set(q1.island_id)
            assert set(q1.island_id)<=set(q0.island_id)
            r0=fit_trait(q0[q0.island_id.isin(ids)],outcome,cfg)
            r1=fit_trait(q1[q1.island_id.isin(ids)],outcome,cfg)
            common.append({"context":context,"outcome":outcome,"common_islands":len(ids),
                "axis_on_common_status":r0["status"],"atomic_on_common_status":r1["status"],
                "axis_on_common_estimate":r0.get("estimate"),"atomic_on_common_estimate":r1.get("estimate"),
                "axis_on_common_p":r0.get("p_two_sided"),"atomic_on_common_p":r1.get("p_two_sided")})
    pd.DataFrame(common).to_csv(out/"h1_common_islands.csv",index=False)
    result={"source_sha256":PINS,"h2_config_sha256":sha(h2path.read_bytes()),"syndrome_config_sha256":sha(scpath.read_bytes()),
            "baseline_replay_pass":True,"h1":safe_rows(h1s),"h2":safe_rows(h2s),"complete3":safe_rows(c3s),
            "h1_common_islands":common,"score_counts_axis":n0,"score_counts_atomic":n1,
            "complete3_axis_summary":s0,"complete3_atomic_summary":s1,
            "scope":"post-hoc exact state-set overlap with selected High/Medium atomic ledger; includes same-sign multistate signatures per unchanged model. Broad flora only. Not independent replication or proof of trait correctness/permissions.",
            "frozen_selector_changed":False,"wcvp_replayed":False,"H4_refitted":False}
    emit(out/"ATOMIC_SENSITIVITY_RESULT.json",result)
    # concise log rows, exact numeric result JSON also retained as aggregate artifact
    print("ATOMIC_SENSITIVITY_RESULT "+json.dumps(result,allow_nan=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--analysis",type=Path,required=True);p.add_argument("--ledger",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();run(a.analysis,a.ledger,a.out)
