#!/usr/bin/env python3
"""FD post-result common-carrier isolation gradients; protected source runner."""
import argparse,functools,hashlib,io,json,pathlib
import numpy as np
import pandas as pd
import fd_pairwise_flower_space as pairspace
from fd_pairwise_flower_space import PIN,read_pin,REGIONS,SEED,incidence
from fd_common_carrier_partner_contrast import C,F,S,triple_eval

REG=("north","tropical","tropical_minus_north")
METRICS=("total","marginal","coupling")
ROLES=("target","donor","delta")
COLS=["distance","area","log_n","log_m","coverage","donor_coverage","abs_lat","climate_pc1","climate_pc2","climate_pc3","climate_pc4"]
PRIMARY="all_recorded__geo_k5__disjoint"

def design(d):
    x=d[COLS].to_numpy(float);n=len(d);p=x.shape[1]
    blocks=sorted(set(d.block));bid={b:i for i,b in enumerate(blocks)}
    groups={}
    for region in REG[:2]:
        rr=[]
        for b,v in d[d.region.eq(region)].groupby("block",sort=True):
            ids=v.index.to_numpy()
            if len(ids)<2 or np.ptp(x[ids,0])<=1e-12:continue
            z=x[ids]-x[ids].mean(0);z/=np.sqrt(len(ids))
            rr.append((bid[b],ids,z))
        groups[region]=rr
    def row(region,w):
        rr=[(b,ids,z) for b,ids,z in groups[region] if w[b]>0]
        if not rr or sum(len(ids)-1 for b,ids,z in rr)-p<4:return None
        X=np.vstack([z*np.sqrt(w[b]) for b,ids,z in rr])
        scale=np.sqrt(np.sum(X**2,axis=0))
        if np.any(scale<=1e-12):return None
        X=X/scale
        if np.linalg.matrix_rank(X)<p or np.linalg.cond(X)>30:return None
        co=np.linalg.pinv(X)[0]/scale[0]
        H=np.zeros(n);pos=0
        for b,ids,z in rr:
            h=co[pos:pos+len(ids)]*np.sqrt(w[b]/len(ids));H[ids]=h-h.mean();pos+=len(ids)
        return H
    def rows(w):
        a=row("north",w);b=row("tropical",w)
        return [a,b,None if a is None or b is None else b-a]
    return blocks,groups,rows

def audit(d,tar,src,tab,perms=999,boots=999):
    d=d[d.region.isin(REG[:2])].replace([np.inf,-np.inf],np.nan).dropna(subset=COLS+["block"]).reset_index(drop=True)
    if not len(d):return {"status":"empty"}
    blocks,groups,rows=design(d);K=len(blocks)
    H0=rows(np.ones(K))
    if any(h is None for h in H0):return {"status":"not_estimable","complete_islands":len(d),"informative_blocks":{r:len(v) for r,v in groups.items()}}
    H=np.vstack(H0)
    ids=d.island_id.tolist();I=len(ids);rr=[tar[i] for i in ids];ss=[src[i] for i in ids]
    names=sorted(set().union(*rr,*ss));ix={s:i for i,s in enumerate(names)}
    A,n=incidence(rr+ss,ix);assert min(n)>=5
    t=tab.loc[names];v=list(t[[C,F,S]].itertuples(index=False,name=None));prof=sorted(set(v));idx={z:i for i,z in enumerate(prof)}
    code=np.array([idx[z] for z in v]);fn=triple_eval(A,n,prof)
    def outcomes(c):
        y=fn(c);return np.column_stack([y[:I],y[I:],y[:I]-y[I:]])
    Y=outcomes(code);point=H@Y
    codes=t.groupby(["genus","grade"],sort=True).ngroup().to_numpy();order=np.argsort(codes,kind="stable")
    movable=np.zeros(len(t),bool)
    for g in t.groupby(["genus","grade"],sort=True).indices.values():
        if len(g)>1 and len(np.unique(code[g]))>1:movable[g]=True
    rng=np.random.default_rng(SEED);draw=np.empty((perms,3,9));E=np.zeros_like(Y)
    for k in range(perms):
        u=np.empty_like(code);u[order]=code[np.lexsort((rng.random(len(code)),codes))]
        yy=outcomes(u);E+=yy/perms;draw[k]=H@yy
    expected=H@E;res=point-expected
    assert np.allclose(point[:,6:],point[:,:3]-point[:,3:6])
    assert np.allclose(point[:,[0,3,6]],point[:,[1,4,7]]+point[:,[2,5,8]])
    rb=[[] for _ in REG];eb=[[] for _ in REG];leave=[[] for _ in REG]
    for k in range(boots):
        for j,h in enumerate(rows(np.bincount(rng.integers(K,size=K),minlength=K))):
            if h is not None:rb[j].append(h@Y);eb[j].append(h@(Y-E))
    for k in range(K):
        w=np.ones(K);w[k]=0
        for j,h in enumerate(rows(w)):
            if h is not None:leave[j].append(h@(Y-E))
    counts={}
    for r in REG[:2]:
        rid=np.concatenate([ids for b,ids,z in groups[r]]) if groups[r] else np.array([],int)
        present=np.asarray(A[np.r_[rid,rid+I]].sum(0)).ravel()>0
        counts[r]={"informative_islands":len(rid),"informative_blocks":len(groups[r]),"movable_species":int((present&movable).sum()),"movable_genera":int(t.loc[present&movable,"genus"].nunique())}
    result={}
    for j,r in enumerate(REG):
        rel=[r] if r in counts else REG[:2]
        gate=all(counts[x]["informative_islands"]>=200 and counts[x]["informative_blocks"]>=30 and counts[x]["movable_species"]>=30 and counts[x]["movable_genera"]>=10 for x in rel) and len(rb[j])>=.9*boots and len(leave[j])==K
        out={"gate":bool(gate),"bootstrap_valid":len(rb[j]),"leave_valid":len(leave[j]),"metrics":{}}
        for k,name in enumerate([a+"__"+b for a in ROLES for b in METRICS]):
            v=draw[:,j,k];o=point[j,k]
            p=min(1.,2*min((1+int(np.sum(v<=o+1e-12)))/(perms+1),(1+int(np.sum(v>=o-1e-12)))/(perms+1)))
            def ci(x):return np.quantile(np.asarray(x)[:,k],[.025,.975]).tolist() if len(x)>=.9*boots else None
            lr=np.asarray(leave[j])[:,k] if len(leave[j]) else None
            out["metrics"][name]={"observed":float(o),"observed_ci":ci(rb[j]),"null_mean":float(expected[j,k]),"residual":float(res[j,k]),"residual_ci":ci(eb[j]),"coherent_p":float(p),"residual_leave_range":[float(lr.min()),float(lr.max())] if lr is not None else None}
        result[r]=out
    return {"status":"estimated","complete_islands":I,"blocks":K,"support":counts,"regions":result,"permutations":perms}

def self_test():
    rng=np.random.default_rng(SEED);n=240
    d=pd.DataFrame(rng.normal(size=(n,len(COLS))),columns=COLS)
    d["region"]=np.repeat(REG[:2],n//2);d["block"]=np.repeat(np.arange(40),6)
    y=d.distance.to_numpy()*np.where(d.region.eq("north"),.2,-.4)+rng.normal(size=40).repeat(6)
    blocks,g,rows=design(d);H=np.vstack(rows(np.ones(len(blocks))))
    assert np.allclose(H@y,[.2,-.4,-.6],atol=1e-10)
    print("FD_PARTNER_ISOLATION_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--self-test",action="store_true");p.add_argument("--original-root",type=pathlib.Path);p.add_argument("--geography",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);a=p.parse_args()
    if a.self_test:self_test();return
    r=a.original_root
    f=pd.read_csv(io.BytesIO(read_pin(r/"fixed/canonical/input/chapter1_status_flora.csv.gz","flora")),compression="gzip",dtype=str,usecols=["island_id","accepted_species","origin_status"])
    t=pd.read_csv(io.BytesIO(read_pin(r/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz","traits")),compression="gzip",dtype=str)
    ass=pd.read_csv(io.BytesIO(read_pin(r/"fixed/source/result/island_source_assignments.csv.gz","assign")),compression="gzip",dtype=str)
    raw=(r/"fixed/source/source/gift_native_mainland_flora.csv.gz").read_bytes()
    mainland=pd.read_csv(io.BytesIO(raw),compression="gzip",dtype=str,usecols=["entity_ID","work_species"])
    geo=pd.read_csv(io.BytesIO(read_pin(a.geography,"geography")),dtype={"island_id":str})
    assert f.island_id.nunique()==4505 and not f.duplicated(["island_id","accepted_species"]).any()
    assert mainland.entity_ID.nunique()==748 and geo.island_id.is_unique and len(geo)==8264
    full=f.groupby("island_id").accepted_species.agg(set).to_dict();pools=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    ass=ass[ass.source_mode.isin(["geo_k5","geo50_climate10"])].copy();ass.source_rank=pd.to_numeric(ass.source_rank)
    assert not ass.duplicated(["island_id","source_mode","source_rank"]).any()
    keys=ass.sort_values("source_rank").groupby(["island_id","source_mode"]).entity_ID.agg(tuple).to_dict()
    @functools.lru_cache(maxsize=128)
    def pool(ids):return frozenset().union(*(pools[i] for i in ids))
    g=geo[geo.island_id.isin(full)].dropna(subset=["island_latitude","island_longitude"]).copy()
    g["region"]=np.where(g.analysis_regime.eq("tropical"),"tropical",np.where(g.island_latitude.gt(23.5),"north",np.where(g.island_latitude.lt(-23.5),"south","unclassified")))
    g=g[g.region.isin(REGIONS)].copy();g["block"]=np.floor(g.island_latitude/10).astype(int).astype(str)+"_"+np.floor(g.island_longitude/10).astype(int).astype(str)
    pairspace.AXES=(C,F,S);base=pairspace.states(t);panels={}
    for scope in ["all_recorded","original_native","high_quality_only"]:
        ff=f[f.origin_status.eq("native")] if scope=="original_native" else f
        tab=base[base.grade.eq("high_high_high")] if scope=="high_quality_only" else base
        scored=set(tab.index);total=ff.groupby("island_id").accepted_species.nunique()
        tar=ff[ff.accepted_species.isin(scored)].groupby("island_id").accepted_species.agg(set).to_dict()
        for mode,nr in [("geo_k5",5),("geo50_climate10",10)]:
            src={};denom={}
            for i in g.island_id:
                if len(tar.get(i,set()))<5:continue
                k=keys.get((i,mode),())
                if len(k)!=nr or any(x not in pools for x in k):continue
                s=pool(k)-full[i];v=s&scored
                if len(v)>=5:src[i]=v;denom[i]=len(s)
            d=g[g.island_id.isin(src)].copy()
            d["distance"]=np.log2(1+d.distance_to_continent_km.where(d.distance_to_continent_km.ge(0)))
            area=next((x for x in ["island_area_km2","area_km2","area"] if x in d),None)
            if area is None:raise ValueError("area")
            d["area"]=np.log2(d[area].where(d[area].gt(0)))
            n=d.island_id.map(lambda i:len(tar[i]));m=d.island_id.map(lambda i:len(src[i]))
            d["log_n"]=np.log2(n);d["log_m"]=np.log2(m)
            d["coverage"]=n/d.island_id.map(total);d["donor_coverage"]=m/d.island_id.map(denom)
            d["abs_lat"]=d.island_latitude.abs()
            k=scope+"__"+mode+"__disjoint"
            panels[k]=audit(d,tar,src,tab,perms=999 if k==PRIMARY else 199)
    p=panels[PRIMARY].get("regions",{}).get("tropical_minus_north",{})
    row=p.get("metrics",{}).get("delta__total",{})
    ci=row.get("residual_ci");lr=row.get("residual_leave_range")
    row["primary_supported"]=bool(p.get("gate") and row.get("coherent_p",1)<.05 and ci and lr and ((ci[0]>0 and lr[0]>0) or (ci[1]<0 and lr[1]<0)))
    out={"botanical_parent":4505,"geography_parent":8264,"source_pins":PIN,"mainland_sha256":hashlib.sha256(raw).hexdigest(),"primary":PRIMARY+":tropical_minus_north:delta__total","controls":COLS[1:],"panels":panels,"warning":"Post-result exploratory isolation association on common triple species. Distance log2(1+km); separate region within-10-degree block equal weights, controls fixed, standardized design condition <=30. Common species-coherent whole triple exchange across target/donor union, paired block bootstrap conditional on fixed observed rosters. One nominated primary; other outcomes and sensitivities descriptive. No causal or pollinator inference, no raw identities exported."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False,indent=2))
    for k,v in panels.items():print("FD_PARTNER_ISOLATION_PANEL",k,json.dumps(v,allow_nan=False),flush=True)
    print("FD_PARTNER_ISOLATION_COMPLETE",json.dumps({"panels":len(panels)}),flush=True)
if __name__=="__main__":main()
