#!/usr/bin/env python3
"""FD post-result audit: color/symmetry versus form/symmetry on identical carriers."""
import argparse,functools,hashlib,io,json,pathlib
from collections import Counter
import numpy as np
import pandas as pd
from scipy import sparse
import fd_pairwise_flower_space as pairspace
from fd_pairwise_flower_space import PIN,read_pin,REGIONS,REG,SEED,incidence
from fd_global_source_architecture import region_design

C,F,S=pairspace.ALL_AXES[:3]
PAIRS=((C,S),(F,S))
STAGES=("pair_specific","matched_islands","common_carriers","island_selection","target_selection","donor_selection","marginal_component","coupling_component")

def q(seq):
    n=len(seq)
    assert n>=2
    return 1-sum(v*(v-1) for v in Counter(seq).values())/(n*(n-1))

def values(tab):
    axes=list(tab.columns)
    return {i:tuple(x) for i,x in zip(tab.index,tab[axes].itertuples(index=False,name=None))}

def rawpair(roster,v):
    return q([v[x] for x in roster])

def triple_eval(A,n,profiles):
    N=A.shape[1];P=len(profiles)
    maps=[]
    for axes in [(0,2),(1,2),(0,),(1,),(2,)]:
        cats=[tuple(v[i] for i in axes) for v in profiles]
        u={x:i for i,x in enumerate(sorted(set(cats)))}
        codes=np.array([u[x] for x in cats])
        maps.append(sparse.csr_matrix((np.ones(P),(np.arange(P),codes)),shape=(P,len(u))))
    den=n*(n-1)
    def qmat(c):
        c=c.tocoo()
        return 1-np.bincount(c.row,weights=c.data*(c.data-1),minlength=len(n))/den
    def fn(code):
        Z=sparse.csr_matrix((np.ones(N),(np.arange(N),code)),shape=(N,P))
        H=A@Z
        qc,qf,c,f,s=[qmat(H@M) for M in maps]
        diff=qc-qf
        marginal=(c-f)*(1-s)
        return np.column_stack([diff,marginal,diff-marginal])
    return fn

def boot_raw(d,common,boots=999):
    keys=sorted(set(d.block));ix={k:i for i,k in enumerate(keys)};K=len(keys)
    def table(x,cols):
        v=x.groupby(["region","block"],sort=True)[cols].mean().reset_index()
        return v[cols].to_numpy(float),np.array([REGIONS.index(z) for z in v.region]),np.array([ix[z] for z in v.block])
    full=table(d,["cs","fs"])
    mat=table(common,["cs","fs","common","target_change","donor_change","marginal","coupling"])
    def aggregate(tab,w):
        y,r,b=tab;out=np.full((3,y.shape[1]),np.nan)
        for j in range(3):
            for c in range(y.shape[1]):
                k=(r==j)&np.isfinite(y[:,c]);ww=w[b[k]]
                if ww.sum()>0:out[j,c]=np.average(y[k,c],weights=ww)
        return np.vstack([out,out[1]-out[0]])
    def calc(w):
        a=aggregate(full,w);b=aggregate(mat,w)
        z0=a[:,0]-a[:,1];z1=b[:,0]-b[:,1];z2=b[:,2]
        ok=np.isfinite(z2)
        assert np.allclose((z2-z1)[ok],(b[:,3]+b[:,4])[ok],atol=1e-10)
        assert np.allclose(z2[ok],(b[:,5]+b[:,6])[ok],atol=1e-10)
        return np.column_stack([z0,z1,z2,z1-z0,b[:,3:7]])
    obs=calc(np.ones(K));rng=np.random.default_rng(SEED)
    draws=np.array([calc(np.bincount(rng.integers(K,size=K),minlength=K)) for _ in range(boots)])
    out={}
    for j,r in enumerate(REG):
        required=[r] if r in REGIONS else ["north","tropical"]
        support=all(len(common[common.region.eq(x)])>=200 and common[common.region.eq(x)].block.nunique()>=30 for x in required)
        row={"common_islands":int(common.region.isin(required).sum()),"common_blocks":int(common.loc[common.region.isin(required),"block"].nunique()),"geographic_gate":bool(support),"stages":{}}
        for i,k in enumerate(STAGES):
            if not np.isfinite(obs[j,i]):continue
            v=draws[:,j,i];v=v[np.isfinite(v)]
            row["stages"][k]={"estimate":float(obs[j,i]),"ci":np.quantile(v,[.025,.975]).tolist() if len(v)>=.9*boots else None}
        out[r]=row
    return out

def coherent(d,tar,src,tab,perms):
    ids=d.island_id.tolist();I=len(ids)
    rr=[tar[i] for i in ids];ss=[src[i] for i in ids]
    species=sorted(set().union(*rr,*ss));ix={s:i for i,s in enumerate(species)}
    A,n=incidence(rr+ss,ix);assert min(n)>=5
    t=tab.loc[species]
    v=list(t[[C,F,S]].itertuples(index=False,name=None))
    profiles=sorted(set(v));idx={z:i for i,z in enumerate(profiles)}
    code=np.array([idx[z] for z in v])
    fn=triple_eval(A,n,profiles)
    G,reg,bid,blocks=region_design(d)
    L=np.zeros((3,len(reg)))
    for j in range(3):
        k=reg==j
        if k.any():L[j,k]=1/k.sum()
    W=sparse.csr_matrix(np.vstack([L,L[1]-L[0]]))@G
    def eval(c):
        z=fn(c);return np.asarray(W@(z[:I]-z[I:]))
    obs=eval(code);rng=np.random.default_rng(SEED)
    codes=t.groupby(["genus","grade"],sort=True).ngroup().to_numpy();order=np.argsort(codes,kind="stable")
    groups=t.groupby(["genus","grade"],sort=True).indices
    movable=np.zeros(len(t),bool)
    for g in groups.values():
        if len(g)>1 and len(np.unique(code[g]))>1:movable[g]=True
    draw=np.empty((perms,4,3))
    for k in range(perms):
        u=np.empty_like(code);u[order]=code[np.lexsort((rng.random(len(code)),codes))]
        draw[k]=eval(u)
    coverage={}
    for r in REGIONS:
        ids=np.flatnonzero(d.region.eq(r).to_numpy())
        present=np.asarray(A[np.r_[ids,ids+I]].sum(0)).ravel()>0
        coverage[r]={"islands":len(ids),"blocks":int(d.loc[d.region.eq(r),"block"].nunique()),"movable_species":int((present&movable).sum()),"movable_genera":int(t.loc[present&movable,"genus"].nunique())}
    out={}
    for j,r in enumerate(REG):
        needed=[r] if r in REGIONS else ["north","tropical"]
        if any(coverage[x]["islands"]==0 for x in needed):continue
        gate=all(coverage[x]["islands"]>=200 and coverage[x]["blocks"]>=30 and coverage[x]["movable_species"]>=30 and coverage[x]["movable_genera"]>=10 for x in needed)
        out[r]={"gate":bool(gate),"metrics":{}}
        for i,name in enumerate(["common_carriers","marginal_component","coupling_component"]):
            x=draw[:,j,i];o=obs[j,i];e=x.mean()
            p=min(1.,2*min((1+int(np.sum(x<=o+1e-12)))/(perms+1),(1+int(np.sum(x>=o-1e-12)))/(perms+1)))
            out[r]["metrics"][name]={"observed":float(o),"null_mean":float(e),"residual":float(o-e),"p":float(p)}
    return {"coverage":coverage,"permutations":perms,"regions":out}

def self_test():
    prof=[("a","a","x"),("a","b","y"),("b","a","x"),("b","a","y")]
    A=sparse.csr_matrix(np.ones((1,8)))
    code=np.repeat(np.arange(4),2)
    got=triple_eval(A,np.array([8]),prof)(code)[0]
    z=[prof[i] for i in code]
    expected=q([(a,c) for a,b,c in z])-q([(b,c) for a,b,c in z])
    assert np.isclose(got[0],expected) and np.isclose(got[0],got[1]+got[2])
    d=pd.DataFrame({"region":["north","tropical"],"block":["a","b"],"cs":[.1,.2],"fs":[0,.1]})
    com=d.assign(common=[.2,.2],target_change=.08,donor_change=.02,marginal=.15,coupling=.05)
    assert np.isclose(boot_raw(d,com,boots=10)["north"]["stages"]["common_carriers"]["estimate"],.2)
    print("FD_PARTNER_MATCH_SELF_TEST_PASS",flush=True)
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
    tabs=[]
    for axes in [*PAIRS,(C,F,S)]:
        pairspace.AXES=axes;tabs.append(pairspace.states(t))
    assert set(tabs[2].index)==set(tabs[0].index)&set(tabs[1].index)
    panels={}
    for scope in ["all_recorded","original_native","high_quality_only"]:
        ff=f[f.origin_status.eq("native")] if scope=="original_native" else f
        ts=[z[z.grade.eq("_".join(["high"]*n))] if scope=="high_quality_only" else z for z,n in zip(tabs,[2,2,3])]
        v=[values(z[list(axes)]) for z,axes in zip(ts,[*PAIRS,(C,F,S)])]
        rt=[ff[ff.accepted_species.isin(z.index)].groupby("island_id").accepted_species.agg(set).to_dict() for z in ts]
        for mode,nr in [("geo_k5",5),("geo50_climate10",10)]:
            rows=[];matched=[];tr={};sr={}
            for row in g[["island_id","region","block"]].itertuples(index=False):
                i=row.island_id;k=keys.get((i,mode),())
                if len(k)!=nr or any(x not in pools for x in k):continue
                s=pool(k)-full[i]
                roster=[x.get(i,set()) for x in rt]
                donor=[s&set(z.index) for z in ts]
                t0=[rawpair(roster[j],v[j]) if len(roster[j])>=5 else np.nan for j in range(2)]
                s0=[rawpair(donor[j],v[j]) if len(donor[j])>=5 else np.nan for j in range(2)]
                delta=np.array(t0)-s0
                if not np.isfinite(delta).any():continue
                base={"region":row.region,"block":row.block,"cs":delta[0],"fs":delta[1]}
                rows.append(base)
                if min(len(roster[2]),len(donor[2]))<5:continue
                assert np.isfinite(delta).all()
                tr[i]=roster[2];sr[i]=donor[2]
                def calc(x):
                    a=[v[2][z] for z in x]
                    zc=q([(c,s) for c,f,s in a]);zf=q([(f,s) for c,f,s in a])
                    marginal=(q([c for c,f,s in a])-q([f for c,f,s in a]))*(1-q([s for c,f,s in a]))
                    return np.array([zc-zf,marginal,zc-zf-marginal])
                x=calc(roster[2]);y=calc(donor[2]);m=x-y
                matched.append(base|{"common":m[0],"target_change":x[0]-(t0[0]-t0[1]),"donor_change":-(y[0]-(s0[0]-s0[1])),"marginal":m[1],"coupling":m[2]})
            key=scope+"__"+mode+"__disjoint"
            d=pd.DataFrame(rows,columns=["region","block","cs","fs"]);dd=pd.DataFrame(matched)
            if not len(dd):panels[key]={"status":"empty"};continue
            raw=boot_raw(d,dd)
            null=coherent(g[g.island_id.isin(tr)].copy(),tr,sr,ts[2],999 if key=="all_recorded__geo_k5__disjoint" else 199)
            for r,z in null["regions"].items():
                assert np.isclose(z["metrics"]["common_carriers"]["observed"],raw[r]["stages"]["common_carriers"]["estimate"],atol=1e-10)
            panels[key]={"status":"estimated","raw":raw,"coherent":null}
    primary=panels["all_recorded__geo_k5__disjoint"].get("coherent",{}).get("regions",{}).get("tropical_minus_north",{})
    test=primary.get("metrics",{}).get("common_carriers",{})
    test["primary_supported"]=bool(primary.get("gate") and test.get("p",1)<.05)
    out={"botanical_parent":4505,"geography_parent":8264,"source_pins":PIN,"panels":panels,"primary":"all_recorded__geo_k5__disjoint:tropical_minus_north:common_carriers","warning":"Post-result audit, one nominated primary null test before rerun; not independent confirmation. Joint Rao Q contrast color-symmetry minus form-symmetry, on common triple carriers. Pointwise paired block intervals conditional on observed rosters. Carrier changes and marginal/coupling decomposition are accounting identities, not causal. No raw IDs exported."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False,indent=2))
    for k,v in panels.items():print("FD_PARTNER_MATCH_PANEL",k,json.dumps(v,allow_nan=False),flush=True)
    print("FD_PARTNER_MATCH_COMPLETE",json.dumps({"panels":len(panels)}),flush=True)
if __name__=="__main__":main()
