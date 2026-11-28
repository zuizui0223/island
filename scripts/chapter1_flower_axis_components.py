#!/usr/bin/env python3
"""Decompose the observed tropical conditional symmetry contrast into its two depth-specific slopes.
Post-outcome interpretation, identical frozen H1/H2 masks. No new endpoint search or functional causal claim.
"""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from chapter1_flower_axis_decoupling import LOCK,digest,trait_pairs,make_frame,design,estimate,SEED,DRAWS,holm
from chapter1_flower_axis_h2 import AUDIT_SHA,build_scores,aggregate
from chapter1_flower_axis_h2_weights import get_repro_matrix,weighted_score,genus_label

METRICS=("radial_given_deep","radial_given_shallow_open","deep_minus_shallow")
def components(c):
    deep=c[:,0]+c[:,1];shallow=c[:,2]+c[:,3]
    if np.any(deep<=0) or np.any(shallow<=0):raise ValueError("empty depth group")
    a=c[:,0]/deep;b=c[:,2]/shallow
    return np.column_stack([a,b,a-b])

def weighted_fit(X,j,y,score,bw):
    Z=X if score is None else np.column_stack([X,score-np.mean(score)])
    A=Z*np.sqrt(bw[:,None])
    if np.linalg.matrix_rank(A)<A.shape[1]:return None
    v=np.linalg.lstsq(A,y*np.sqrt(bw[:,None]),rcond=None)[0][j,:]
    assert abs(v[2]-v[0]+v[1])<1e-9
    return v

def stress(X,j,M,C,mnames,blocks,R=None,rs=None,rnames=None):
    rnames=[] if rnames is None else rnames
    gs=sorted({genus_label(x) for x in [*mnames,*rnames]});gm={x:i for i,x in enumerate(gs)}
    mg=np.array([gm[genus_label(x)] for x in mnames]);rg=np.array([gm[genus_label(x)] for x in rnames],dtype=int)
    _,bi=np.unique(blocks,return_inverse=True);nb=bi.max()+1
    Y=components(np.asarray(M@C));s=None if R is None else weighted_score(R,rs,np.ones(len(rs)))
    rng=np.random.default_rng(SEED);out={k:[] for k in ("geography_only","genus_only","joint")}
    for _ in range(DRAWS):
        w=rng.exponential(1,len(gs));yy=components(np.asarray(M@(C*w[mg,None])))
        ss=None if R is None else weighted_score(R,rs,w[rg])
        bw=np.bincount(rng.integers(nb,size=nb),minlength=nb)[bi]
        for k,y,q,b in (("geography_only",Y,s,bw),("genus_only",yy,ss,np.ones(len(Y))),("joint",yy,ss,bw)):
            z=weighted_fit(X,j,y,q,b)
            if z is not None:out[k].append(z.tolist())
    return {k:{"valid_draws":len(v),"range_025_975":np.quantile(v,[.025,.975],axis=0).T.tolist() if len(v)>=.9*DRAWS else None} for k,v in out.items()}

def self_test():
    x=components(np.array([[3.,1.,1.,3.],[1.,3.,3.,1.]]))
    assert np.allclose(x,[[.75,.25,.5],[.25,.75,-.5]])
    t=np.arange(80,dtype=float);X=np.column_stack([np.ones(80),t])
    y=np.column_stack([.1+.02*t,.3-.01*t,-.2+.03*t])
    assert np.allclose(weighted_fit(X,1,y,None,np.ones(80)),[.02,-.01,.03])
    print("FLORAL_COMPONENT_SELF_TEST_PASS")

def main():
    p=argparse.ArgumentParser()
    for k in ("flora","ledger","geography","state-audit","h1-parent","h2-parent","out"):p.add_argument("--"+k,type=Path)
    p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    for k,h in LOCK.items():
        if digest(getattr(a,k))!=h:raise ValueError("source hash mismatch "+k)
    if digest(a.state_audit)!=AUDIT_SHA:raise ValueError("audit hash drift")
    h1=json.loads(a.h1_parent.read_text());h2=json.loads(a.h2_parent.read_text())
    assert h1["run"]==45579499069 and h2["run"]==45584985080
    h1p=next(x for x in h1["panels"] if x["scope"]=="all_recorded" and x["region"]=="tropical" and x["endpoint"]=="conditional_difference" and x["nmin"]==10)
    old={(x["scope"],x["roster"],x["score"]):x for x in h2["panels"]}
    f=pd.read_csv(a.flora,dtype=str,usecols=["island_id","accepted_species","origin_status"])
    v=trait_pairs(pd.read_csv(a.ledger,dtype=str))
    geo=pd.read_csv(a.geography,dtype={"island_id":str,"spatial_block":str})
    s=build_scores(pd.read_csv(a.state_audit,dtype=str));s["SC"]=s.self_incompatibility;s["autonomy"]=s.autonomous_selfing_capacity
    d,M,C,_=make_frame(f,v,geo,"all_recorded")
    d=d.loc[d.analysis_regime.eq("tropical")&d.carrier_n.ge(10)&d.deep_n.ge(5)&d.shallow_n.ge(5)].sort_values("island_id")
    pair=f.loc[f.accepted_species.isin(v.index)]
    panels=[]
    for roster,score in [("none","none"),("whole_flora","min2"),("pair_carriers","SC"),("pair_carriers","autonomy")]:
        prior=h1p if score=="none" else old[("all_recorded",roster,score)]
        if score=="none":z=d.copy();R=None;rs=None;rnames=None
        else:
            pool=f if roster=="whole_flora" else pair
            z=d.merge(aggregate(pool,s,score),left_on="island_id",right_index=True,how="left",validate="one_to_one").loc[lambda x:x.score_n.ge(10)].dropna(subset=["score"])
            R,rnames,rs=get_repro_matrix(z,pool,s,score)
            assert np.allclose(weighted_score(R,rs,np.ones(len(rs))),z.score)
        assert len(z)==prior["n"]
        X,j=design(z);b=z.spatial_block.to_numpy(str)
        mm=M[z.mi.to_numpy(int)];ix=np.flatnonzero(np.asarray(mm.sum(0)).ravel()>0);mm=mm[:,ix];cc=C[ix];mnames=v.index[ix].to_list()
        Y=components(np.asarray(mm@cc))
        X2=X if score=="none" else np.column_stack([X,z.score.to_numpy(float)-z.score.mean()])
        fit=estimate(X2,Y,b,j)
        assert fit is not None
        expected=prior["beta"] if score=="none" else prior["adjusted_same_islands"]["beta"]
        assert np.isclose(fit["2"]["beta"],expected,atol=1e-10,rtol=0)
        assert np.isclose(fit["2"]["beta"],fit["0"]["beta"]-fit["1"]["beta"],atol=1e-10,rtol=0)
        metrics={k:fit[str(i)] for i,k in enumerate(METRICS)}
        # Only two component signs in unadjusted H1 are additional formal tests.
        # H2 components and the already-tested difference are descriptive.
        if score=="none":holm([metrics[k] for k in METRICS[:2]])
        r={"roster":roster,"score":score,"n":len(z),"blocks":len(set(b)),"species":len(mnames),"metrics":metrics,"mean_radial_given_depth":Y[:,:2].mean(0).tolist(),"sensitivity":stress(X,j,mm,cc,mnames,b,R,rs,rnames)}
        panels.append(r)
        print("FLORAL_COMPONENT_PANEL "+json.dumps(r,allow_nan=False),flush=True)
    result={"owner":"island","panels":panels,"metric_order":list(METRICS),"draws":DRAWS,
      "warning":"Post-hoc decomposition of prior conditional contrast on identical masks, not independent evidence. H1 two component slopes Holm2; H2 component estimates descriptive. Coherent genus and geographic weights are sensitivity ranges, not calibrated confidence intervals. No direct functional or causal claims."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("FLORAL_COMPONENT_COMPLETE 4",flush=True)

if __name__=="__main__":main()
