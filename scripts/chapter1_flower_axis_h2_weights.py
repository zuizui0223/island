#!/usr/bin/env python3
"""Coherent genus and spatial sensitivity of the frozen tropical H2 extension.
Weights shared across repeated accepted-name genera in both morphological and
reproductive aggregates; fixed carrier masks/covariates. Anonymous results.
"""
import argparse,json,re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from chapter1_flower_axis_decoupling import LOCK,digest,trait_pairs,make_frame,rates,design,estimate,METRICS,SEED,DRAWS
from chapter1_flower_axis_h2 import AUDIT_SHA,build_scores,aggregate,one_panel

def genus_label(name):
    m=re.match(r"^([A-Z][A-Za-z-]+)\s+[a-z][A-Za-z-]+",name)
    # Unparsed accepted names receive their own weight, not a single pooled genus.
    return m.group(1) if m else "unparsed:"+name

def get_repro_matrix(z,pool,scores,key):
    sc=scores[key].dropna().sort_index()
    j=pool.loc[pool.island_id.isin(z.island_id),["island_id","accepted_species"]].merge(sc.rename("score"),left_on="accepted_species",right_index=True,validate="many_to_one")
    names=sorted(j.accepted_species.unique())
    ix=pd.Series(np.arange(len(z)),index=z.island_id)
    si=pd.Series(np.arange(len(names)),index=names)
    R=sparse.csr_matrix((np.ones(len(j)),(j.island_id.map(ix).to_numpy(),j.accepted_species.map(si).to_numpy())),shape=(len(z),len(names)))
    assert R.nnz==len(j) and np.all(R.data==1)
    y=sc.loc[names].to_numpy(float)
    assert np.all(np.asarray(R.sum(1)).ravel()>=10)
    return R,names,y

def weighted_score(R,s,w):
    d=np.asarray(R@w).ravel()
    if np.any(d<=0):raise ValueError("zero reweighted denominator")
    return np.asarray(R@(w*s)).ravel()/d

def slopes(X,j,y,sc,bw):
    A=X*np.sqrt(bw[:,None])
    Z=np.column_stack([X,sc-np.mean(sc)])*np.sqrt(bw[:,None])
    if np.linalg.matrix_rank(A)<A.shape[1] or np.linalg.matrix_rank(Z)<Z.shape[1]:return None
    t=y*np.sqrt(bw)
    raw=np.linalg.lstsq(A,t,rcond=None)[0][j]
    fit=np.linalg.lstsq(Z,t,rcond=None)[0]
    return [float(raw),float(fit[j]),float(fit[j]-raw),float(fit[-1])]

def quantiles(a):
    return {"valid_draws":len(a),"range_025_975":np.quantile(np.asarray(a),[.025,.975],axis=0).T.tolist() if len(a)>=.9*DRAWS else None}

def sensitivity(X,j,M,C,R,s,m_names,r_names,blocks):
    labels=sorted({genus_label(x) for x in [*m_names,*r_names]})
    g={v:i for i,v in enumerate(labels)}
    mg=np.array([g[genus_label(x)] for x in m_names])
    rg=np.array([g[genus_label(x)] for x in r_names])
    y0=rates(np.asarray(M@C))[:,METRICS.index("conditional_difference")]
    sc0=weighted_score(R,s,np.ones(len(s)))
    _,bi=np.unique(blocks,return_inverse=True);nb=int(bi.max())+1
    rng=np.random.default_rng(SEED);out={k:[] for k in ("geography_only","genus_only","joint")}
    for _ in range(DRAWS):
        w=rng.exponential(1,len(labels))
        y=rates(np.asarray(M@(C*w[mg,None])))[:,METRICS.index("conditional_difference")]
        sc=weighted_score(R,s,w[rg])
        bw=np.bincount(rng.integers(nb,size=nb),minlength=nb)[bi]
        for k,yy,ss,ww in (("geography_only",y0,sc0,bw),("genus_only",y,sc,np.ones(len(y))),("joint",y,sc,bw)):
            v=slopes(X,j,yy,ss,ww)
            if v is not None and np.all(np.isfinite(v)):out[k].append(v)
    return {"metric_order":["unadjusted_isolation","conditioned_isolation","signed_slope_change","reproduction_score"],
        "weight_groups":len(labels),"unparsed_score_species":sum(genus_label(x).startswith("unparsed:") for x in r_names),
        "observed":slopes(X,j,y0,sc0,np.ones(len(y0))),
        **{k:quantiles(v) for k,v in out.items()},"not_calibrated_ci":True}

def self_test():
    R=sparse.csr_matrix([[1,1,0],[0,1,1]])
    s=np.array([1.,-1.,1.])
    assert np.allclose(weighted_score(R,s,np.ones(3)),0)
    assert np.allclose(weighted_score(R,s,np.array([2.,1.,2.])),1/3)
    assert genus_label("Abies alba")==genus_label("Abies firma")
    assert genus_label("unparsed")=="unparsed:unparsed"
    rng=np.random.default_rng(1);a=rng.normal(size=80);b=rng.normal(size=80);X=np.column_stack([np.ones(80),a])
    y=.2*a+.4*b
    assert np.isclose(slopes(X,1,y,b,np.ones(80))[1],.2)
    print("H2_COHERENT_WEIGHT_SELF_TEST_PASS")

def main():
    p=argparse.ArgumentParser()
    for k in ("flora","ledger","geography","state-audit","parent","out"):p.add_argument("--"+k,type=Path)
    p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    for k,h in LOCK.items():
        if digest(getattr(a,k))!=h:raise ValueError("source hash mismatch "+k)
    if digest(a.state_audit)!=AUDIT_SHA:raise ValueError("state audit drift")
    expected=json.loads(a.parent.read_text())
    assert expected["run"]==45584985080 and len(expected["panels"])==24
    old={(x["scope"],x["roster"],x["score"]):x for x in expected["panels"]}
    f=pd.read_csv(a.flora,dtype=str,usecols=["island_id","accepted_species","origin_status"])
    ledger=pd.read_csv(a.ledger,dtype=str);geo=pd.read_csv(a.geography,dtype={"island_id":str,"spatial_block":str})
    scores=build_scores(pd.read_csv(a.state_audit,dtype=str));scores["SC"]=scores.self_incompatibility;scores["autonomy"]=scores.autonomous_selfing_capacity
    panels=[]
    for scope in ("all_recorded","high_only","literal_native"):
        v=trait_pairs(ledger,scope=="high_only")
        d,M,C,_=make_frame(f,v,geo,scope)
        d=d.loc[d.analysis_regime.eq("tropical")&d.carrier_n.ge(10)&d.deep_n.ge(5)&d.shallow_n.ge(5)].sort_values("island_id")
        y=rates(np.asarray(M[d.mi.to_numpy(int)]@C))[:,METRICS.index("conditional_difference")]
        pool=f if scope!="literal_native" else f.loc[f.origin_status.eq("native")]
        pair=pool.loc[pool.accepted_species.isin(v.index)]
        for roster,ff in (("whole_flora",pool),("pair_carriers",pair)):
            for key in ("min2","complete3","SC","autonomy"):
                sc=aggregate(ff,scores,key)
                row=one_panel(d,y,sc,key,roster,scope,10);prior=old[(scope,roster,key)]
                assert row["n"]==prior["n"] and row["blocks"]==prior["blocks"] and row["status"]==prior["status"]
                out={"scope":scope,"score":key,"roster":roster,"n":row["n"],"blocks":row["blocks"],"status":row["status"],"parent_holm_p":prior["holm_p"]}
                if row["status"]=="fit":
                    for k in ("adjusted_same_islands","unadjusted_same_islands"):
                        assert np.isclose(row[k]["beta"],prior[k]["beta"],atol=1e-10,rtol=0)
                    z=d.merge(sc,left_on="island_id",right_index=True,how="left",validate="one_to_one").loc[lambda t:t.score_n.ge(10)].dropna(subset=["score"])
                    assert len(z)==row["n"]
                    X,j=design(z);blocks=z.spatial_block.to_numpy(str)
                    R,r_names,s=get_repro_matrix(z,ff,scores,key)
                    assert np.allclose(weighted_score(R,s,np.ones(len(s))),z.score)
                    mm=M[z.mi.to_numpy(int)];cols=np.flatnonzero(np.asarray(mm.sum(0)).ravel()>0)
                    mm=mm[:,cols];cc=C[cols];names=v.index[cols].to_list()
                    sens=sensitivity(X,j,mm,cc,R,s,names,r_names,blocks)
                    assert np.isclose(sens["observed"][1],prior["adjusted_same_islands"]["beta"],atol=1e-10,rtol=0)
                    out["sensitivity"]=sens
                    ran=sens["joint"]["range_025_975"]
                    out["positive_parent_and_joint_sensitivity"]=bool(prior["holm_p"]<.05 and ran is not None and ran[1][0]>0)
                panels.append(out)
    assert len(panels)==24
    result={"owner":"island","parent_run":45584985080,"draws":DRAWS,"panels":panels,
        "warning":"Post-outcome taxon/geography weighting sensitivity, not calibrated confidence intervals or new p-values. All fitted H2 estimates reconstructed; same genus weight for morphology and reproductive score, fixed cohort and original nuisance covariates. Not causal mediation, independent replication, or evidence-quality invariance."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    for row in panels:print("H2_COHERENT_WEIGHT_PANEL "+json.dumps(row,allow_nan=False),flush=True)
    print("H2_COHERENT_WEIGHT_COMPLETE 24",flush=True)

if __name__=="__main__":main()
