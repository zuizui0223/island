#!/usr/bin/env python3
"""FD-owned protected replay: regional source-relative floral architecture.
Post-pilot; no animal data; candidate mainland pools are not ancestors.
Species-coherent label permutation across source and target union.
Only anonymous aggregate output.
"""
import argparse, functools, hashlib, io, json, pathlib
import numpy as np
import pandas as pd
from scipy import sparse
from fd_original_tropical_pollinator_bridge import PIN, read_pin, pair_states
from fd_original_source_relative_pollinator import membership

REGIONS=("north","tropical","south")
ENDS=("radial_deep","bilateral_deep","deep_any","radial_minus_bilateral")
SEED=20261125

def region_design(d):
    d=d.reset_index(drop=True)
    key=list(zip(d.region,d.block))
    ks=sorted(set(key)); idx={k:i for i,k in enumerate(ks)}
    ci=np.array([idx[k] for k in key]);n=np.bincount(ci)
    G=sparse.csr_matrix((1/n[ci],(ci,np.arange(len(d)))),shape=(len(ks),len(d)))
    reg=np.array([REGIONS.index(k[0]) for k in ks])
    blocks=sorted({k[1] for k in ks});bi={k:i for i,k in enumerate(blocks)}
    bid=np.array([bi[k[1]] for k in ks])
    return G,reg,bid,blocks

def summarize(d,tar,src,tab,perms=999,boots=999,matched_high=None):
    if not len(d):return {"status":"empty"}
    ids=d.island_id.tolist()
    rr=[tar[i] for i in ids];ss=[src[i] for i in ids]
    species=sorted(set().union(*rr,*ss));ix={s:i for i,s in enumerate(species)}
    A,n=membership(rr,ix);M,m=membership(ss,ix)
    assert min(n)>=5 and min(m)>=5
    t=tab.loc[species].copy()
    vals=np.array([[a=="actinomorphic" and b=="deep",a=="zygomorphic" and b=="deep",b=="deep"] for a,b in t.pair],float)
    V=np.column_stack([vals,vals[:,0]-vals[:,1]])
    codes=t.groupby(["genus","grade"],sort=True).ngroup().to_numpy()
    order=np.argsort(codes,kind="stable")
    # Conditional expectation of whole pair labels, preserving their coupling.
    _,inv=np.unique(codes,return_inverse=True)
    E=np.zeros((int(inv.max())+1,V.shape[1]));np.add.at(E,inv,V);E/=np.bincount(inv)[:,None];E=E[inv]
    G,reg,bid,blocks=region_design(d)
    W=(G@(A-M)).tocsr()
    extra={}
    if matched_high is not None:
        ht,hs=matched_high
        hr=[ht[i] for i in ids]; hsr=[hs[i] for i in ids]
        assert all(h<=l for h,l in zip(hr,rr)) and all(h<=l for h,l in zip(hsr,ss))
        AH,nh=membership(hr,ix);MH,mh=membership(hsr,ix)
        assert min(nh)>=5 and min(mh)>=5
        WH=(G@(AH-MH)).tocsr()
        W=WH-W
        extra={"estimand":"matched high minus all source contrast; common genus-grade reference over all-grade union","median_high_target":float(np.median(nh)),"median_high_source":float(np.median(mh))}
    Y=W@V; EY=W@E
    K=len(blocks);R=len(REGIONS)
    def means(z,w):
        ww=w[bid];out=np.full((R,z.shape[1]),np.nan)
        for j in range(R):
            ok=reg==j
            if ww[ok].sum()>0:out[j]=np.average(z[ok],axis=0,weights=ww[ok])
        return np.vstack([out,out[1]-out[0]])
    point=means(Y,np.ones(K));expected=means(EY,np.ones(K));res=point-expected
    rng=np.random.default_rng(SEED)
    weights=[np.bincount(rng.integers(K,size=K),minlength=K) for _ in range(boots)]
    boot=np.array([means(Y-EY,w) for w in weights])
    rawboot=np.array([means(Y,w) for w in weights])
    expboot=np.array([means(EY,w) for w in weights])
    leave=np.array([means(Y-EY,1-np.eye(K,dtype=int)[k]) for k in range(K)])
    # Collapse to regional linear operators so permutations retain all repeated taxa.
    L=np.zeros((R,len(reg)))
    for j in range(R):
        z=reg==j
        if z.any():L[j,z]=1/z.sum()
    H=sparse.csr_matrix(L)@W
    H=sparse.vstack([H,H[1]-H[0]]).tocsr()
    draws=np.empty((perms,R+1,V.shape[1]))
    for k in range(perms):
        shuffled=np.lexsort((rng.random(len(V)),codes))
        U=np.empty_like(V);U[order]=V[shuffled]
        draws[k]=H@(U-E)
    movable=np.any(np.abs(V-E)>1e-12,axis=1)
    assert np.allclose(point[:R][np.isfinite(point[:R])],(H@V)[:R][np.isfinite(point[:R])])
    counts={j:{"islands":int((d.region==j).sum()),"blocks":int((reg==i).sum()),"movable_taxa":int(np.sum(movable & (np.asarray(np.abs(H[i]).sum(0)).ravel()>1e-14)))} for i,j in enumerate(REGIONS)}
    gc=pd.factorize(t.genus,sort=True)[0]
    def interval(v):
        v=v[np.isfinite(v)]
        return np.quantile(v,[.025,.975]).tolist() if len(v)>=.9*boots else None
    def influence(j,i):
        h=H[j].toarray().ravel()
        weights=h*(V[:,i]-E[:,i])
        contributions=np.bincount(gc,weights=weights)
        assert np.isclose(contributions.sum(),res[j,i],atol=1e-12)
        active=contributions[np.abs(contributions)>1e-14]
        ordered=active[np.argsort(-np.abs(active))]
        total=float(np.abs(ordered).sum())
        movable=(np.abs(h)>1e-14)&(np.abs(V[:,i]-E[:,i])>1e-12)
        return {"nonzero_genera":int(len(active)),"endpoint_weighted_movable_taxa":int(movable.sum()),"endpoint_weighted_movable_genera":int(len(np.unique(gc[movable]))),"top1_absolute_share":float(abs(ordered[0])/total) if total else None,"top5_absolute_share":float(np.abs(ordered[:5]).sum()/total) if total else None,"fixed_design_drop_top":{"1":float(res[j,i]-ordered[:1].sum()),"3":float(res[j,i]-ordered[:3].sum()),"5":float(res[j,i]-ordered[:5].sum()),"10":float(res[j,i]-ordered[:10].sum())},"single_genus_drop_range":[float(np.min(res[j,i]-active)),float(np.max(res[j,i]-active))] if len(active) else None}
    out={}
    for j,name in enumerate((*REGIONS,"tropical_minus_north")):
        c=counts[name] if name in counts else None
        gate=all(counts[z]["islands"]>=200 and counts[z]["blocks"]>=30 and counts[z]["movable_taxa"]>=30 for z in ([name] if c else ["north","tropical"]))
        out[name]={"support":c,"gate":bool(gate),"endpoints":{}}
        for i,k in enumerate(ENDS):
            if not np.isfinite(res[j,i]):continue
            v=boot[:,j,i];v=v[np.isfinite(v)]
            q=leave[:,j,i];q=q[np.isfinite(q)]
            ci=np.quantile(v,[.025,.975]).tolist() if len(v)>=.9*boots else None
            lr=[float(q.min()),float(q.max())] if len(q)==K else None
            sim=draws[:,j,i]
            p=float((1+(np.abs(sim)>=abs(res[j,i])-1e-12).sum())/(perms+1))
            out[name]["endpoints"][k]={"observed_source_delta":float(point[j,i]),"genus_grade_expected_delta":float(expected[j,i]),"residual":float(res[j,i]),"residual_ci":ci,"observed_ci":interval(rawboot[:,j,i]),"expected_ci":interval(expboot[:,j,i]),"genus_influence":influence(j,i),"leave_range":lr,"coherent_p":p,"valid_bootstrap":len(v),"valid_leave":len(q)}
    return {"status":"estimated","islands":len(d),"blocks":K,"species":len(species),"median_target":float(np.median(n)),"median_source":float(np.median(m)),"regions":out,**extra}

def main():
    p=argparse.ArgumentParser();p.add_argument("--original-root",required=True,type=pathlib.Path);p.add_argument("--geography",required=True,type=pathlib.Path);p.add_argument("--out",required=True,type=pathlib.Path);p.add_argument("--permutations",type=int,default=999);a=p.parse_args();r=a.original_root
    f=pd.read_csv(io.BytesIO(read_pin(r/"fixed/canonical/input/chapter1_status_flora.csv.gz","flora")),compression="gzip",dtype=str,usecols=["island_id","accepted_species","origin_status"])
    t=pd.read_csv(io.BytesIO(read_pin(r/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz","traits")),compression="gzip",dtype=str)
    ass=pd.read_csv(io.BytesIO(read_pin(r/"fixed/source/result/island_source_assignments.csv.gz","assign")),compression="gzip",dtype=str)
    raw=(r/"fixed/source/source/gift_native_mainland_flora.csv.gz").read_bytes()
    mainland=pd.read_csv(io.BytesIO(raw),compression="gzip",dtype=str,usecols=["entity_ID","work_species"])
    geo=pd.read_csv(io.BytesIO(read_pin(a.geography,"geography")),dtype={"island_id":str})
    assert f.island_id.nunique()==4505 and not f.duplicated(["island_id","accepted_species"]).any()
    assert mainland.entity_ID.nunique()==748 and geo.island_id.is_unique and len(geo)==8264
    full=f.groupby("island_id").accepted_species.agg(set).to_dict()
    pools=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    ass=ass[ass.source_mode.isin(["geo_k5","geo50_climate10"])].copy();ass.source_rank=pd.to_numeric(ass.source_rank)
    keys=ass.sort_values("source_rank").groupby(["island_id","source_mode"]).entity_ID.agg(tuple).to_dict()
    @functools.lru_cache(maxsize=128)
    def pool(ids):return frozenset().union(*(pools[i] for i in ids))
    g=geo[geo.island_id.isin(full)].dropna(subset=["island_latitude","island_longitude"]).copy()
    g["region"]=np.where(g.analysis_regime.eq("tropical"),"tropical",np.where(g.island_latitude.gt(23.5),"north",np.where(g.island_latitude.lt(-23.5),"south","unclassified")))
    g=g[g.region.isin(REGIONS)].copy()
    g["block"]=np.floor(g.island_latitude/10).astype(int).astype(str)+"_"+np.floor(g.island_longitude/10).astype(int).astype(str)
    pair=pair_states(t);panels={};matched={};all_inputs={}
    for scope in ["all_recorded","original_native","high_quality_only"]:
        ff=f[f.origin_status.eq("native")] if scope=="original_native" else f
        tab=pair[pair.grade.eq("high")] if scope=="high_quality_only" else pair
        scored=set(tab.index)
        tar=ff[ff.accepted_species.isin(scored)].groupby("island_id").accepted_species.agg(set).to_dict()
        for mode,nr in [("geo_k5",5),("geo50_climate10",10)]:
            incl={};disj={}
            for i in g.island_id:
                if len(tar.get(i,set()))<5:continue
                k=keys.get((i,mode),())
                if len(k)!=nr or any(x not in pools for x in k):continue
                s=pool(k)&scored;v=s-full[i]
                if len(v)<5:continue
                incl[i]=s;disj[i]=v
                assert not (v&tar[i])
            d=g[g.island_id.isin(disj)].copy()
            for label,src in [("inclusive",incl),("disjoint",disj)]:
                key="__".join([scope,mode,label])
                panels[key]=summarize(d,tar,src,tab,perms=a.permutations)
                if scope=="all_recorded":all_inputs[(mode,label)]=(tar,src)
                if scope=="high_quality_only":
                    at,asrc=all_inputs[(mode,label)]
                    assert set(d.island_id)<=set(asrc)
                    matched[mode+"__"+label]=summarize(d,at,asrc,pair,perms=a.permutations,matched_high=(tar,src))
    tests=[e for v in panels.values() for rr in v.get("regions",{}).values() if rr["gate"] for e in rr["endpoints"].values()]
    last=0.
    for j,e in enumerate(sorted(tests,key=lambda q:q["coherent_p"])):
        last=max(last,min(1.,e["coherent_p"]*(len(tests)-j)));e["holm_p"]=last
    qt=[e for v in matched.values() for rr in v.get("regions",{}).values() if rr["gate"] for e in rr["endpoints"].values()]
    last=0.
    for j,e in enumerate(sorted(qt,key=lambda q:q["coherent_p"])):
        last=max(last,min(1.,e["coherent_p"]*(len(qt)-j)));e["holm_p"]=last
    out={"quality_matched_panels":matched,"quality_matched_holm_family_size":len(qt),"botanical_parent":4505,"geography_parent":8264,"source_pins":PIN,"mainland_sha256":hashlib.sha256(raw).hexdigest(),"panels":panels,"holm_family_size":len(tests),"warning":"Post-pilot exploratory source-pool contrasts, not causal or isolation slopes. Source assignment is not ancestral identity. Geography-block equal weighting, complete pair labels exchanged coherently in global genus/combined evidence stratum. No direct pollinator evidence. Native ascertainment incomplete. No raw identities exported."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    for k,v in panels.items():print("FD_SOURCE_ARCH_PANEL",k,json.dumps(v,allow_nan=False),flush=True)
    for k,v in matched.items():print("FD_SOURCE_ARCH_MATCHED",k,json.dumps(v,allow_nan=False),flush=True)
    print("FD_SOURCE_ARCH_COMPLETE",json.dumps({"panels":len(panels),"holm_family_size":len(tests),"matched":len(matched),"matched_holm_family_size":len(qt)}),flush=True)

if __name__=="__main__":main()
