#!/usr/bin/env python3
"""FD four-axis categorical space: source-relative diversity, not continuous volume."""
import argparse, functools, hashlib, io, json, pathlib
import numpy as np
import pandas as pd
from scipy import sparse
from fd_original_tropical_pollinator_bridge import PIN,read_pin,accepted_genus,pair_states
from fd_global_source_architecture import REGIONS,region_design
AXES=("flower_primary_color","floral_form","floral_symmetry","tube_depth_class")
METRICS=("joint_Q","additive_Q","joint_minus_fixed_margin_Q","rarefied_profiles5_per_draw","radial_deep_minus_bilateral_deep")
REG=(*REGIONS,"tropical_minus_north")
PRIMARY="all_recorded__geo_k5__disjoint"
SEED=20261125

def states(t):
    z=t[t.trait_name.isin(AXES)&t.evidence_scope.eq("species_direct")&t.evidence_quality.isin(["high","medium"])&t.normalized_value.notna()].copy()
    z["normalized_value"]=z.normalized_value.str.strip()
    z=z[z.normalized_value.ne("")&~z.normalized_value.str.contains("|",regex=False)&~z.normalized_value.str.lower().isin(["none","nan","null","unknown","unresolved"])]
    assert not z.duplicated(["accepted_species","trait_name"]).any()
    v=z.pivot(index="accepted_species",columns="trait_name",values="normalized_value").dropna(subset=list(AXES))
    q=z.pivot(index="accepted_species",columns="trait_name",values="evidence_quality").loc[v.index,list(AXES)]
    v["grade"]=q.agg("_".join,axis=1)
    v["genus"]=[accepted_genus(s) for s in v.index]
    return v.dropna(subset=["genus"])

def incidence(rosters,idx):
    rr=[];cc=[];n=np.array(list(map(len,rosters)),int)
    for i,s in enumerate(rosters):
        rr.extend([i]*len(s));cc.extend(idx[z] for z in sorted(s))
    return sparse.csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(rosters),len(idx))),n

def evaluator(A,n,profiles):
    P=len(profiles);N=A.shape[1]
    maps=[]
    for i in range(len(AXES)):
        _,c=np.unique([v[i] for v in profiles],return_inverse=True)
        maps.append(sparse.csr_matrix((np.ones(P),(np.arange(P),c)),shape=(P,int(c.max())+1)))
    direction=np.array([(v[2]=="actinomorphic" and v[3]=="deep")-(v[2]=="zygomorphic" and v[3]=="deep") for v in profiles],float)
    den=n*(n-1)
    def row_sum(c,fun):
        c=c.tocoo()
        return np.bincount(c.row,weights=fun(c.data,c.row),minlength=len(n))
    def q(c):return 1-row_sum(c,lambda x,r:x*(x-1))/den
    def metrics(code):
        Z=sparse.csr_matrix((np.ones(N),(np.arange(N),code)),shape=(N,P))
        C=(A@Z).tocsr()
        j=q(C)
        marg=np.column_stack([q(C@M) for M in maps])
        # Exact fixed-marginal expected joint Q under independent axis assignments.
        expected=1-np.prod(1-marg,axis=1)
        def rare(x,r):
            miss=np.ones(len(x))
            for k in range(5):miss*=np.maximum(n[r]-x-k,0)/(n[r]-k)
            return (1-miss)/5
        return np.column_stack([j,marg.mean(1),j-expected,row_sum(C,rare),np.asarray(C@direction).ravel()/n])
    return metrics

def summarize(d,tar,src,tab,perms=399,boots=999):
    if not len(d):return {"status":"empty"}
    d=d.reset_index(drop=True)
    ids=d.island_id.tolist();rr=[tar[i] for i in ids];ss=[src[i] for i in ids]
    names=sorted(set().union(*rr,*ss));idx={s:i for i,s in enumerate(names)}
    A,n=incidence(rr+ss,idx);assert min(n)>=5
    t=tab.loc[names];vectors=list(t[list(AXES)].itertuples(index=False,name=None))
    profiles=sorted(set(vectors));ix={v:i for i,v in enumerate(profiles)}
    code=np.array([ix[v] for v in vectors])
    fn=evaluator(A,n,profiles);I=len(d)
    G,reg,bid,blocks=region_design(d);K=len(blocks);R=len(REGIONS)
    def regional(z,w):
        ww=w[bid];out=np.full((R,len(METRICS)),np.nan)
        for j in range(R):
            ok=reg==j
            if ww[ok].sum()>0:out[j]=np.average(z[ok],axis=0,weights=ww[ok])
        return np.vstack([out,out[1]-out[0]])
    def blockdelta(c):
        z=fn(c);return G@(z[:I]-z[I:])
    B=blockdelta(code);point=regional(B,np.ones(K))
    groups=t.groupby(["genus","grade"],sort=True).indices
    groups=[np.array(v,int) for v in groups.values() if len(v)>1 and len(np.unique(code[v]))>1]
    movable=np.zeros(len(t),bool)
    for g in groups:movable[g]=True
    rng=np.random.default_rng(SEED);null=np.empty((perms,R+1,len(METRICS)));expected=np.zeros_like(B)
    for k in range(perms):
        c=code.copy()
        for g in groups:c[g]=rng.permutation(c[g])
        z=blockdelta(c);expected+=z/perms;null[k]=regional(z,np.ones(K))
    E=regional(expected,np.ones(K));res=point-E
    W=[np.bincount(rng.integers(K,size=K),minlength=K) for _ in range(boots)]
    rawboot=np.array([regional(B,w) for w in W]);resboot=np.array([regional(B-expected,w) for w in W])
    leave=np.array([regional(B-expected,1-np.eye(K,dtype=int)[k]) for k in range(K)])
    counts={}
    for j,r in enumerate(REGIONS):
        ids=np.flatnonzero(d.region.eq(r).to_numpy())
        present=np.asarray(A[np.r_[ids,ids+I]].sum(0)).ravel()>0
        counts[r]={"islands":len(ids),"blocks":int((reg==j).sum()),"movable_species_present":int(np.sum(present&movable)),"movable_genera_present":int(t.loc[present&movable,"genus"].nunique()),"median_target":float(np.median(n[ids])) if len(ids) else None,"median_donor":float(np.median(n[ids+I])) if len(ids) else None}
    def ci(v):
        v=v[np.isfinite(v)]
        return np.quantile(v,[.025,.975]).tolist() if len(v)>=.9*boots else None
    out={}
    for j,r in enumerate(REG):
        relevant=[r] if r in counts else ["north","tropical"]
        gate=all(counts[z]["islands"]>=200 and counts[z]["blocks"]>=30 and counts[z]["movable_species_present"]>=30 and counts[z]["movable_genera_present"]>=10 for z in relevant)
        out[r]={"support":counts.get(r),"gate":bool(gate),"metrics":{}}
        for i,k in enumerate(METRICS):
            if not np.isfinite(point[j,i]):continue
            z=null[:,j,i]
            # Conservative two-tail rank p, valid for asymmetric null distributions.
            p=min(1.,2*min((1+int(np.sum(z<=point[j,i]+1e-12)))/(perms+1),(1+int(np.sum(z>=point[j,i]-1e-12)))/(perms+1)))
            lv=leave[:,j,i];lv=lv[np.isfinite(lv)]
            out[r]["metrics"][k]={"observed_delta":float(point[j,i]),"observed_ci":ci(rawboot[:,j,i]),"genus_grade_null_mean":float(E[j,i]),"residual":float(res[j,i]),"residual_ci":ci(resboot[:,j,i]),"coherent_p":float(p),"leave_range":[float(lv.min()),float(lv.max())] if len(lv)==K else None}
    return {"status":"estimated","islands":I,"distinct_species":len(names),"distinct_profiles":len(profiles),"permutations":perms,"regions":out}

def self_test():
    p=[("a","a","actinomorphic","deep"),("b","b","zygomorphic","deep")]
    A=sparse.csr_matrix(np.ones((1,6)));n=np.array([6])
    q=evaluator(A,n,p)(np.array([0,0,0,1,1,1]))[0]
    assert np.isclose(q[0],.6) and np.isclose(q[1],.45)
    assert np.isclose(q[2],.6-(1-.4**3))
    assert np.isclose(q[3],.4) and np.isclose(q[4],0)
    # Equal profile counts with different pairing: additive Q unchanged, joint Q changes.
    p=[("a","a","actinomorphic","deep"),("a","b","actinomorphic","deep"),("b","a","actinomorphic","deep"),("b","b","actinomorphic","deep")]
    A=sparse.csr_matrix(np.ones((1,8)));n=np.array([8]);fn=evaluator(A,n,p)
    x=fn(np.array([0]*4+[3]*4))[0];y=fn(np.repeat(np.arange(4),2))[0]
    assert np.isclose(x[1],y[1]) and y[0]>x[0] and y[2]>x[2]
    print("FD_FOUR_AXIS_SELF_TEST_PASS",flush=True)

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
    four=states(t);pair=pair_states(t)
    assert set(four.index)<=set(pair.index)
    coverage={}
    for scope in ["all_recorded","original_native","high_quality_only"]:
        ff=f[f.origin_status.eq("native")] if scope=="original_native" else f
        tab=four[four.grade.eq("high_high_high_high")] if scope=="high_quality_only" else four
        pp=pair[pair.grade.eq("high")] if scope=="high_quality_only" else pair
        z={}
        for name,tt in [("pair",pp),("four",tab)]:
            counts=ff[ff.accepted_species.isin(tt.index)].groupby("island_id").accepted_species.nunique()
            c=g.assign(n=g.island_id.map(counts).fillna(0))
            z[name]={r:{"n5":int((v.n>=5).sum()),"n30":int((v.n>=30).sum())} for r,v in c.groupby("region")}
        coverage[scope]=z
    print("FD_FOUR_AXIS_COVERAGE",json.dumps(coverage),flush=True)
    panels={}
    for scope in ["all_recorded","original_native","high_quality_only"]:
        ff=f[f.origin_status.eq("native")] if scope=="original_native" else f
        tab=four[four.grade.eq("high_high_high_high")] if scope=="high_quality_only" else four
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
            d=g[g.island_id.isin(disj)].copy()
            for label,src in [("inclusive",incl),("disjoint",disj)]:
                key="__".join([scope,mode,label])
                panels[key]=summarize(d,tar,src,tab,perms=999 if key==PRIMARY else 199)
                print("FD_FOUR_AXIS_PANEL",key,json.dumps(panels[key],allow_nan=False),flush=True)
    # Primary family fixed before run: five outcomes x north/tropics/contrast.
    # Ineligible slots count as p=1; sensitivities cannot rescue primary inference.
    rows=[]
    for r in ["north","tropical","tropical_minus_north"]:
        z=panels[PRIMARY].get("regions",{}).get(r,{})
        for k in METRICS:
            row=z.get("metrics",{}).get(k,{})
            rows.append((row.get("coherent_p",1.) if z.get("gate") else 1.,row))
    last=0.
    for j,(pv,row) in enumerate(sorted(rows,key=lambda v:v[0])):
        last=max(last,min(1.,pv*(len(rows)-j)));row["primary_holm_p"]=last
    out={"botanical_parent":4505,"geography_parent":8264,"coverage":coverage,"source_pins":PIN,"mainland_sha256":hashlib.sha256(raw).hexdigest(),"primary":PRIMARY,"primary_family_size":15,"panels":panels,"warning":"Post-pilot exploratory categorical four-axis space; not continuous hypervolume, realized pollination or ancestral pool. Pair vs four coverage audit; complete-case target and donor rosters. Geography intervals conditional on fixed species and donor roster; coherent whole-profile permutation across union within genus x exact four-trait evidence pattern. Sensitivities nominal only. Rarefied repertoire at 5 draws is an expected distinct profile count divided by five. No raw identifiers exported."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False,indent=2))
    for k,v in panels.items():print("FD_FOUR_AXIS_FINAL",k,json.dumps(v,allow_nan=False),flush=True)
    print("FD_FOUR_AXIS_COMPLETE",json.dumps({"panels":len(panels),"primary_family_size":15}),flush=True)
if __name__=="__main__":main()
