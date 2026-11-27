#!/usr/bin/env python3
"""FD source-relative beta diversity at fixed repertoire sample size.
Post-pilot protected analysis: functional cell deficits vs balanced mismatch,
not historical loss/colonization or measured pollination.
"""
import argparse,functools,hashlib,io,json,pathlib
import numpy as np
import pandas as pd
from scipy import sparse
import fd_pairwise_flower_space as pairspace
from fd_pairwise_flower_space import PIN,read_pin,REGIONS,REG,SEED,incidence
PAIRS=(("flower_primary_color","floral_form"),("floral_symmetry","tube_depth_class"))
METRICS=("shared_cells_per_draw","target_only_cells_per_draw","source_only_cells_per_draw","source_minus_target_repertoire","balanced_exclusive_cells_per_draw","unbalanced_exclusive_cells_per_draw","target_repertoire","source_repertoire")
PRIMARY="all_recorded__geo_k5__disjoint__n5"

def occupancy(C,n,k):
    miss=np.ones(C.shape,float)
    for j in range(k):
        miss*=np.maximum(n[:,None]-C-j,0)/(n[:,None]-j)
    return 1-miss

def calc(C,n,k,I):
    p=occupancy(C,n,k)
    a=p[:I];b=p[I:]
    shared=np.sum(a*b,axis=1)/k
    to=np.sum(a*(1-b),axis=1)/k
    so=np.sum(b*(1-a),axis=1)/k
    dif=so-to;bal=2*np.minimum(so,to)
    # Identities concern expected numbers of occupied cells at equal sample size,
    # not expected samplewise min/max or expected Sørensen/Simpson ratios.
    out=np.column_stack([shared,to,so,dif,bal,np.abs(dif),shared+to,shared+so])
    assert np.allclose(to+so,bal+np.abs(dif))
    assert np.allclose(dif,out[:,7]-out[:,6])
    return out

def design(d,tag):
    keys=sorted(set(d[tag]));ki={v:i for i,v in enumerate(keys)}
    cells=sorted(set(zip(d.region,d[tag])));ix={v:i for i,v in enumerate(cells)}
    c=np.array([ix[v] for v in zip(d.region,d[tag])]);n=np.bincount(c)
    G=sparse.csr_matrix((1/n[c],(c,np.arange(len(d)))),shape=(len(cells),len(d)))
    region=np.array([REGIONS.index(x[0]) for x in cells]);ids=np.array([ki[x[1]] for x in cells])
    def agg(B,w):
        out=np.full((len(REGIONS),B.shape[1]),np.nan)
        for j in range(len(REGIONS)):
            a=region==j;ww=w[ids[a]]
            if ww.sum()>0:out[j]=np.average(B[a],axis=0,weights=ww)
        return np.vstack([out,out[1]-out[0]])
    return G,agg,len(keys)

def summarize(d,tar,src,tab,axes,k,perms,boots=999):
    if not len(d):return {"status":"empty","pairs":0}
    d=d.reset_index(drop=True);I=len(d);ids=d.island_id.tolist()
    rr=[tar[i] for i in ids];ss=[src[i] for i in ids]
    names=sorted(set().union(*rr,*ss));ix={s:i for i,s in enumerate(names)}
    A,n=incidence(rr+ss,ix);assert min(n)>=k
    t=tab.loc[names];v=list(t[list(axes)].itertuples(index=False,name=None))
    profiles=sorted(set(v));idx={z:i for i,z in enumerate(profiles)};code=np.array([idx[z] for z in v]);P=len(profiles)
    def fn(c):
        Z=sparse.csr_matrix((np.ones(len(names)),(np.arange(len(names)),c)),shape=(len(names),P))
        return calc((A@Z).toarray(),n,k,I)
    gc=t.groupby(["genus","grade"],sort=True).ngroup().to_numpy();order=np.argsort(gc,kind="stable")
    movable=np.zeros(len(t),bool)
    for g in t.groupby(["genus","grade"],sort=True).indices.values():
        if len(np.unique(code[g]))>1:movable[g]=True
    support={}
    for r in REGIONS:
        ii=np.flatnonzero(d.region.eq(r).to_numpy());present=np.asarray(A[np.r_[ii,ii+I]].sum(0)).ravel()>0
        support[r]={"islands":len(ii),"blocks":int(d.loc[ii,"block"].nunique()),"source_groups":int(d.loc[ii,"source_group"].nunique()),"movable_species":int((present&movable).sum()),"movable_genera":int(t.loc[present&movable,"genus"].nunique()),"median_target":float(np.median(n[ii])) if len(ii) else None,"median_donor":float(np.median(n[ii+I])) if len(ii) else None}
    W={tag:design(d,tag) for tag in ("block","source_group")};rng=np.random.default_rng(SEED)
    Y=fn(code);E=np.zeros_like(Y)
    draws={tag:np.empty((perms,len(REG),len(METRICS))) for tag in W}
    for j in range(perms):
        u=np.empty_like(code);u[order]=code[np.lexsort((rng.random(len(code)),gc))]
        z=fn(u);E+=z/perms
        for tag,(G,agg,N) in W.items():draws[tag][j]=agg(G@z,np.ones(N))
    groups={}
    for tag,(G,agg,N) in W.items():
        B=G@Y;D=G@E;point=agg(B,np.ones(N));ex=agg(D,np.ones(N));res=point-ex
        raw=[];eb=[];leave=[]
        for j in range(boots):
            w=np.bincount(rng.integers(N,size=N),minlength=N);raw.append(agg(B,w));eb.append(agg(B-D,w))
        for j in range(N):
            w=np.ones(N);w[j]=0;leave.append(agg(B-D,w))
        raw=np.array(raw);eb=np.array(eb);leave=np.array(leave);out={}
        for j,r in enumerate(REG):
            req=[r] if r in REGIONS else ["north","tropical"]
            if any(support[q]["islands"]==0 for q in req):continue
            gate=all(support[q]["islands"]>=200 and support[q]["blocks"]>=30 and support[q]["source_groups"]>=30 and support[q]["movable_species"]>=30 and support[q]["movable_genera"]>=10 for q in req)
            out[r]={"gate":bool(gate),"metrics":{}}
            for i,m in enumerate(METRICS):
                z=draws[tag][:,j,i];o=point[j,i]
                p=min(1.,2*min((1+int(np.sum(z<=o+1e-12)))/(perms+1),(1+int(np.sum(z>=o-1e-12)))/(perms+1)))
                def ci(v):
                    x=v[:,j,i];x=x[np.isfinite(x)]
                    return np.quantile(x,[.025,.975]).tolist() if len(x)>=.9*boots else None
                lv=leave[:,j,i]
                out[r]["metrics"][m]={"observed":float(o),"observed_ci":ci(raw),"null_mean":float(ex[j,i]),"residual":float(res[j,i]),"residual_ci":ci(eb),"coherent_p":float(p),"residual_leave_range":[float(lv.min()),float(lv.max())] if np.isfinite(lv).all() else None}
        groups[tag]=out
    return {"status":"estimated","pairs":I,"species":len(names),"profiles":P,"k":k,"support":support,"group_weightings":groups}

def self_test():
    from itertools import combinations
    x=[0,0,1,2,3,3];y=[0,1,1,1,4,4];k=3
    C=np.array([np.bincount(z,minlength=5) for z in [x,y]])
    a=calc(C,np.array([6,6]),k,1)[0]
    v=[]
    for p in combinations(x,k):
        for q in combinations(y,k):
            u,w=set(p),set(q);v.append([len(u&w)/k,len(u-w)/k,len(w-u)/k])
    assert np.allclose(a[:3],np.mean(v,axis=0)) and np.isclose(a[3],a[2]-a[1])
    assert np.isclose(a[4]+a[5],a[1]+a[2])
    z=calc(np.array([[5,0],[0,5]]),np.array([5,5]),5,1)[0]
    assert np.allclose(z[:6],[0,.2,.2,0,.4,0])
    print("FD_SOURCE_REPERTOIRE_SELF_TEST_PASS",flush=True)

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
    panels={}
    specs=[("all_recorded","geo_k5",5,"disjoint",5),("all_recorded","geo_k5",5,"inclusive",5),("all_recorded","geo50_climate10",10,"disjoint",5),("original_native","geo_k5",5,"disjoint",5),("high_quality_only","geo_k5",5,"disjoint",5),("all_recorded","geo_k5",5,"disjoint",10)]
    for axes in PAIRS:
        pairspace.AXES=axes;base=pairspace.states(t);pk="__".join(axes)
        for scope,mode,nr,over,k in specs:
            ff=f[f.origin_status.eq("native")] if scope=="original_native" else f
            tab=base[base.grade.eq("high_high")] if scope=="high_quality_only" else base
            scored=set(tab.index);tar=ff[ff.accepted_species.isin(scored)].groupby("island_id").accepted_species.agg(set).to_dict()
            src={};group={}
            for i in g.island_id:
                if len(tar.get(i,set()))<k:continue
                ids=keys.get((i,mode),())
                if len(ids)!=nr or any(x not in pools for x in ids):continue
                z=pool(ids)&scored;disj=z-full[i]
                if len(disj)<k:continue
                src[i]=disj if over=="disjoint" else z;group[i]=ids[0]
            d=g[g.island_id.isin(src)].copy();d["source_group"]=d.island_id.map(group)
            kk="__".join([scope,mode,over,"n"+str(k)]);key=pk+"__"+kk
            panels[key]=summarize(d,tar,src,tab,axes,k,999 if kk==PRIMARY else 199)
            print("FD_SOURCE_REPERTOIRE_PROGRESS",key,json.dumps({"pairs":len(d)},allow_nan=False),flush=True)
    tests=[]
    for axes in PAIRS:
        v=panels["__".join(axes)+"__"+PRIMARY].get("group_weightings",{})
        B=v.get("block",{}).get("tropical_minus_north",{});D=v.get("source_group",{}).get("tropical_minus_north",{})
        for m in ("source_minus_target_repertoire","balanced_exclusive_cells_per_draw"):
            x=B.get("metrics",{}).get(m,{})
            y=D.get("metrics",{}).get(m,{})
            tests.append((x.get("coherent_p",1.) if B.get("gate") else 1.,x,y,bool(B.get("gate") and D.get("gate"))))
    last=0
    for j,(p,x,y,gate) in enumerate(sorted(tests,key=lambda z:z[0])):
        last=max(last,min(1.,p*(len(tests)-j)));x["primary_holm_p"]=last
        def sign(e):
            c=e.get("residual_ci");l=e.get("residual_leave_range")
            if not c or not l:return 0
            return 1 if c[0]>0 and l[0]>0 else -1 if c[1]<0 and l[1]<0 else 0
        x["primary_supported"]=bool(gate and last<.05 and y.get("coherent_p",1)<.05 and sign(x)!=0 and sign(x)==sign(y))
    out={"botanical_parent":4505,"geography_parent":8264,"source_pins":PIN,"mainland_sha256":hashlib.sha256(raw).hexdigest(),"primary":PRIMARY,"primary_family_size":4,"panels":panels,"warning":"Post-pilot source-relative expected functional-cell occupancy at equal independent species draw size; expected unique cell counts divided by k. Balanced exclusive = 2 min(expected target-only, expected source-only); unbalanced = absolute difference, not expectation of random-sample min or historical turnover/nestedness. Signed deficit positive means fewer expected target cells. Within-genus exact evidence-grade whole-pair exchange, geographic and first-source grouping; overlapping panels not independent. Not ecological loss, colonization or pollinator replacement. No protected identities exported."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    for k,v in panels.items():print("FD_SOURCE_REPERTOIRE_PANEL",k,json.dumps(v,allow_nan=False),flush=True)
    print("FD_SOURCE_REPERTOIRE_COMPLETE",json.dumps({"panels":len(panels),"primary_family_size":4}),flush=True)
if __name__=="__main__":main()
