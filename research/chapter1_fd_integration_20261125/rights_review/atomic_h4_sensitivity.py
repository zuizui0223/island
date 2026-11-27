#!/usr/bin/env python3
"""Post-hoc strict atomic H4: preserve original experiment parents, weights and models."""
from __future__ import annotations
import argparse, hashlib, json, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
from scipy.stats import t
from atomic_direct_sensitivity import ROOT, PINS, assert_hash, read_zip, strict_mask, gate, safe_rows, emit
from island_v2.chapter1_h2_complete_selfing_sensitivity import build_species_score
from island_v2.chapter1_v14_h4_exact_h2_score_bridge import run_exact_h2_bridge, analyse_exact_family, exact_h2_scores

GIT_PINS={
    "results/geography_20260924/h4_ra_corrected_effect_rows.csv.gz":"fef8346705bdd2168c81fa2af783f7e5ad0e66a5",
    "results/geography_20260924/h4_arch_corrected_effect_rows.csv.gz":"0b10e077ecb798abed59a375cecfa22f3a4757f7",
    "results/geography_20260924/h4_exact_corrected.csv":"3ece77d515dc36821b56d7530df90fbd4486f8dd",
    "config/chapter1_v14_h4_exact_h2_score_bridge.yml":"6a6c0c9c42a8ca4e0df3456db8a475f4cc3d89ec",
}
def git_pins():
    out={}
    for p,h in GIT_PINS.items():
        b=(ROOT/p).read_bytes()
        got=hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
        if got!=h:raise ValueError("frozen Git blob mismatch: "+p)
        out[p]={"git_blob":got,"sha256":hashlib.sha256(b).hexdigest()}
    return out
def scores(a,cfg):
    parts=[]
    for name in ("selfing_core","generalized_accessible"):
        z=build_species_score(a,cfg["syndromes"][name],minimum_informative_traits=int(cfg["score_definition"]["minimum_informative_traits"]),require_all=False).dropna(subset=["score"]).copy()
        z["syndrome"]=name
        z["soft_membership"]=(z.score+1)/2
        assert z.soft_membership.between(0,1).all()
        parts.append(z[["accepted_species","syndrome","soft_membership"]])
    return pd.concat(parts,ignore_index=True)
def finite_t(frame):
    x=frame.copy()
    if "estimate" not in x:return x
    g=pd.to_numeric(x.n_publications,errors="coerce")-1
    z=pd.to_numeric(x.estimate,errors="coerce")/pd.to_numeric(x.se,errors="coerce")
    valid=x.evaluable & (g>0) & (x.se>0)
    x["finite_publication_t_two_sided_p"]=np.where(valid,2*t.sf(abs(z),g),np.nan)
    x["finite_publication_t_one_sided_negative_p"]=np.where(valid,t.cdf(z,g),np.nan)
    q=t.ppf(.975,g)
    x["ci_low"]=x.estimate-q*x.se
    x["ci_high"]=x.estimate+q*x.se
    return x
def run(analysis,ledger,out):
    out.mkdir(parents=True,exist_ok=True)
    for p,k in [(analysis,"analysis"),(ledger,"ledger")]:assert_hash(p.read_bytes(),k)
    pins=git_pins()
    with zipfile.ZipFile(analysis) as z,zipfile.ZipFile(ledger) as d:
        a=read_zip(z,"snapshot/chapter1_trait_state_audit_direct.csv.gz","audit",dtype=str,keep_default_na=False)
        atomic=read_zip(d,"direct_species_trait_ledger.csv.gz","atomic",dtype=str,keep_default_na=False)
    b,_=strict_mask(a,atomic)
    synd=yaml.safe_load((ROOT/"config/chapter1_pr138_pollination_syndromes.yml").read_bytes())
    cfg=yaml.safe_load((ROOT/"config/chapter1_v14_h4_exact_h2_score_bridge.yml").read_bytes())
    ra=pd.read_csv(ROOT/"results/geography_20260924/h4_ra_corrected_effect_rows.csv.gz")
    arch=pd.read_csv(ROOT/"results/geography_20260924/h4_arch_corrected_effect_rows.csv.gz")
    old=pd.read_csv(ROOT/"results/geography_20260924/h4_exact_corrected.csv")
    s0=scores(a,synd);s1=scores(b,synd)
    assert s0.groupby("syndrome").size().to_dict()=={"generalized_accessible":3389,"selfing_core":2651}
    assert s1.groupby("syndrome").size().to_dict()=={"generalized_accessible":3385,"selfing_core":2298}
    baseline,bs=run_exact_h2_bridge(s0,ra,arch,cfg)
    baseline.to_csv(out/"H4_axis_replay.csv",index=False)
    g=gate(baseline,old,["family","analysis"],["estimate","se","two_sided_p"],[1e-9,1e-9,1e-9],["evaluable","n_cells","n_publications","n_sites","n_species"])
    emit(out/"H4_BASELINE_GATE.json",g)
    print("H4_ATOMIC_BASELINE_GATE "+json.dumps(g),flush=True)
    if not g["pass"]:raise SystemExit("Frozen H4 baseline mismatch; no strict inference")
    strict,ss=run_exact_h2_bridge(s1,ra,arch,cfg)
    strict=finite_t(strict)
    strict.to_csv(out/"H4_atomic.csv",index=False)
    parents={"reproductive_assurance":ra,"floral_architecture":arch}
    outrows=[];cohort=[]
    for fam,fc in cfg["families"].items():
        syn=str(fc["H2_syndrome"]);sn=str(fc["score_name"])
        x=exact_h2_scores(s0,syndrome=syn,score_name=sn)
        y=exact_h2_scores(s1,syndrome=syn,score_name=sn)
        ids=set(x.species_key)&set(y.species_key)
        assert set(y.species_key)<=set(x.species_key)
        xy=x.merge(y,on="species_key",suffixes=("_old","_strict"),validate="one_to_one")
        cohort.append({"family":fam,"axis_scored_species":len(x),"atomic_scored_species":len(y),"shared_scored_species":len(ids),"changed_scores_on_shared_species":int((abs(xy[sn+"_old"]-xy[sn+"_strict"])>1e-12).sum())})
        for label,frame in [("axis_on_atomic_roster",x),("atomic_on_atomic_roster",y)]:
            frame=frame.loc[frame.species_key.isin(ids)].rename(columns={"species_key":"accepted_species",sn:"soft_membership"}).copy()
            frame["syndrome"]=syn
            rows,_=analyse_exact_family(frame,parents[fc["parent"]],family=fam,family_config=fc,publication_total_weight=float(cfg["analysis"]["publication_total_weight"]))
            outrows.extend({"scope":label,**r} for r in rows)
    common=finite_t(pd.DataFrame(outrows))
    for _,group in common.groupby(["family","analysis"]):
        assert len(group)==2
        if group.evaluable.all():
            for k in ["n_cells","n_publications","n_sites","n_species"]:assert group[k].nunique()==1,k
    common.to_csv(out/"H4_common_species.csv",index=False)
    result={"source_archive_pins":{k:PINS[k] for k in ["analysis","ledger","audit","atomic"]},"git_inputs":pins,
        "baseline_gate":g,"axis_support":bs,"atomic_support":ss,"cohort":cohort,"strict":safe_rows(strict),"common_species":safe_rows(common),
        "scope":"posthoc exact atomic state support of original H2-to-H4 score associations; same corrected parent experiments; no independent replication, causal mediation or rights promotion",
        "H3_changed":False,"original_H4_changed":False}
    emit(out/"H4_ATOMIC_RESULT.json",result)
    print("H4_ATOMIC_RESULT "+json.dumps(result,allow_nan=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--analysis",type=Path,required=True)
    p.add_argument("--ledger",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.analysis,a.ledger,a.out)
