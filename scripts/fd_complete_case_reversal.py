#!/usr/bin/env python3
"""Post-result audit of complete-case reversal; no new ecological mechanism test."""
import argparse,functools,hashlib,io,json,pathlib
import numpy as np
import pandas as pd
from fd_four_axis_space import PIN,read_pin,states,pair_states,REGIONS,AXES
SEED=20261125
REG=(*REGIONS,"tropical_minus_north")
TIERS=("color","form","four")
STAGES=("pair_full","pair_matched","target_restricted","donor_restricted","both_restricted","island_selection","target_selection","donor_selection","total_change")
def scores(pair):
    return {s:float((a=="actinomorphic" and b=="deep")-(a=="zygomorphic" and b=="deep")) for s,(a,b) in pair.pair.items()}
def tier_sets(t,pair,scope):
    z=t[t.trait_name.isin(AXES[:2])&t.evidence_scope.eq("species_direct")&t.evidence_quality.isin(["high","medium"])&t.normalized_value.notna()].copy()
    z["normalized_value"]=z.normalized_value.str.strip()
    z=z[z.normalized_value.ne("")&~z.normalized_value.str.contains("|",regex=False)&~z.normalized_value.str.lower().isin(["none","nan","null","unknown","unresolved"])]
    if scope=="high_quality_only":z=z[z.evidence_quality.eq("high")]
    base=set(pair.index)
    color=base&set(z.loc[z.trait_name.eq(AXES[0]),"accepted_species"])
    form=base&set(z.loc[z.trait_name.eq(AXES[1]),"accepted_species"])
    return {"color":color,"form":form,"four":color&form}

def paired_audit(d,df,boots=999):
    # d: pair-eligible islands with pair target/donor scores. df: subset
    # with restricted target and donor scores, identical geography.
    keys=sorted(set(d.block));ix={b:i for i,b in enumerate(keys)};K=len(keys)
    def blocktable(x,cols):
        z=x.groupby(["region","block"],sort=True)[cols].mean().reset_index()
        return z[cols].to_numpy(float),np.array([REGIONS.index(k) for k in z.region]),np.array([ix[k] for k in z.block])
    all_y,all_r,all_b=blocktable(d,["pair"])
    mat_y,mat_r,mat_b=blocktable(df,["pair","target_only","donor_only","both"])
    def aggregate(y,r,b,w):
        out=np.full((len(REGIONS),y.shape[1]),np.nan)
        for j in range(len(REGIONS)):
            k=r==j;ww=w[b[k]]
            if ww.sum()>0:out[j]=np.average(y[k],axis=0,weights=ww)
        return np.vstack([out,out[1]-out[0]])
    def calc(w):
        a=aggregate(all_y,all_r,all_b,w)[:,0];m=aggregate(mat_y,mat_r,mat_b,w)
        p,t,s,b=m.T
        # Order-independent Shapley accounting of two roster restrictions.
        tsel=((t-p)+(b-s))/2;ssel=((s-p)+(b-t))/2
        total=b-a
        out=np.column_stack([a,p,t,s,b,p-a,tsel,ssel,total])
        ok=np.isfinite(total)
        assert np.allclose(out[ok,5:8].sum(1),total[ok],atol=1e-12)
        return out
    point=calc(np.ones(K));rng=np.random.default_rng(SEED)
    draw=np.array([calc(np.bincount(rng.integers(K,size=K),minlength=K)) for _ in range(boots)])
    result={}
    for j,r in enumerate(REG):
        a=d[d.region.eq(r)] if r in REGIONS else d[d.region.isin(["north","tropical"])]
        b=df[df.region.eq(r)] if r in REGIONS else df[df.region.isin(["north","tropical"])]
        relevant=[r] if r in REGIONS else ["north","tropical"]
        support=all(len(df[df.region.eq(z)])>=200 and df[df.region.eq(z)].block.nunique()>=30 for z in relevant)
        q={"pair_islands":len(a),"matched_islands":len(b),"matched_blocks":int(b.block.nunique()),"geographic_support":bool(support),"stages":{}}
        for i,k in enumerate(STAGES):
            if not np.isfinite(point[j,i]):continue
            v=draw[:,j,i];v=v[np.isfinite(v)]
            q["stages"][k]={"value":float(point[j,i]),"ci":np.quantile(v,[.025,.975]).tolist() if len(v)>=.9*boots else None}
        result[r]=q
    return result
def self_test():
    d=pd.DataFrame({"region":["north","tropical"],"block":["a","b"],"pair":[.1,.2]})
    f=d.assign(target_only=[.2,.3],donor_only=[.3,.4],both=[.4,.5])
    z=paired_audit(d,f,boots=10)
    q=z["tropical"]["stages"]
    assert np.isclose(q["target_selection"]["value"],.1) and np.isclose(q["donor_selection"]["value"],.2)
    print("FD_CARRIER_AUDIT_SELF_TEST_PASS",flush=True)
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
    four=states(t);pair=pair_states(t);panels={}
    for scope in ["all_recorded","original_native","high_quality_only"]:
        ff=f[f.origin_status.eq("native")] if scope=="original_native" else f
        pp=pair[pair.grade.eq("high")] if scope=="high_quality_only" else pair
        scored=set(pp.index);sc=scores(pp)
        tiers=tier_sets(t,pp,scope)
        expected=set(four[four.grade.eq("high_high_high_high")].index) if scope=="high_quality_only" else set(four.index)
        assert tiers["four"]==expected
        tar=ff[ff.accepted_species.isin(scored)].groupby("island_id").accepted_species.agg(set).to_dict()
        def mean(s):return sum(sc[x] for x in s)/len(s)
        for mode,nr in [("geo_k5",5),("geo50_climate10",10)]:
            for label in ["inclusive","disjoint"]:
                rows=[];restricted={tier:[] for tier in TIERS}
                for row in g[["island_id","region","block"]].itertuples(index=False):
                    i=row.island_id
                    if len(tar.get(i,set()))<5:continue
                    k=keys.get((i,mode),())
                    if len(k)!=nr or any(x not in pools for x in k):continue
                    source=pool(k)&scored;disjoint=source-full[i]
                    if len(disjoint)<5:continue
                    s=source if label=="inclusive" else disjoint
                    a0=tar[i];ta=mean(a0);sa=mean(s)
                    base={"region":row.region,"block":row.block,"pair":ta-sa}
                    rows.append(base)
                    for tier,eligible in tiers.items():
                        a1=a0&eligible;s1=s&eligible;sd=disjoint&eligible
                        if min(len(a1),len(sd))<5:continue
                        at=mean(a1);st=mean(s1)
                        restricted[tier].append(base|{"target_only":at-sa,"donor_only":ta-st,"both":at-st})
                d=pd.DataFrame(rows,columns=["region","block","pair"])
                for tier in TIERS:
                    df=pd.DataFrame(restricted[tier],columns=["region","block","pair","target_only","donor_only","both"])
                    if not len(d) or not len(df):continue
                    key="__".join([scope,mode,label,tier])
                    panels[key]=paired_audit(d,df)
    out={"botanical_parent":4505,"geography_parent":8264,"source_pins":PIN,"mainland_sha256":hashlib.sha256(raw).hexdigest(),"panels":panels,"warning":"Post-result complete-case selection accounting. Regional block averages, same resampled geographic blocks for all stages; fixed observed species and donor pools. Shapley split is descriptive, not causal. Pointwise intervals, no multiple-test claims. Restricted source may change composition despite unchanged candidate geography; no protected identities exported."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False,indent=2))
    for k,v in panels.items():print("FD_CARRIER_PANEL",k,json.dumps(v,allow_nan=False),flush=True)
    print("FD_CARRIER_COMPLETE",json.dumps({"panels":len(panels)}),flush=True)
if __name__=="__main__":main()
