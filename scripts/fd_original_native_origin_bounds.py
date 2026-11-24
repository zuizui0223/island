#!/usr/bin/env python3
"""Exact conditional partial-identification bounds for missing native status on frozen designs."""
import argparse, functools, hashlib, io, json, math, pathlib
import numpy as np
import pandas as pd
from fd_original_tropical_pollinator_bridge import PIN,read_pin
from fd_original_source_relative_pollinator import membership
from fd_original_tubular_depth_diet_split import cached_bird_animals,states,MODES

BIRD_PIN="6ae6aecf0a2850f91147cea78f0f3e5b6c348f4e1104fa78bf12a43372d7feab"
SCOPES=(("all_recorded",5),("all_recorded",30),("original_native",5),("high_quality_only",5))
ENDS={"form":("tubular",),"form_depth":("tubular","deep","tubular_and_deep")}
def bird_input(path):
    raw=path.read_bytes()
    assert hashlib.sha256(raw).hexdigest()==BIRD_PIN
    d=json.loads(raw)
    assert d["source_artifact"]==14901918629 and d["source_bird_count"]==377 and d["source_island_count"]==38
    lab=np.asarray(d["diet_codes"],int);fam=np.asarray(d["family_blocks"],int);gen=np.asarray(d["genus_blocks"],int)
    assert len(lab)==len(fam)==len(gen)==377 and set(lab)=={0,1,2}
    inc={str(k):np.asarray(v,int) for k,v in d["positive_incidence"]}
    assert len(inc)==38 and all(len(v)==len(set(v)) and (v>=0).all() and (v<377).all() for v in inc.values())
    return lab,fam,gen,inc
def weights(d,counts,blocks):
    # Equal original bird archipelago weights; only within-group differences.
    x=np.column_stack([np.log1p(d.bee),np.log1p(counts[:,2]),np.log1p(counts[:,1]),np.log1p(counts[:,0]),
                      d.log_area,d.log_distance,d.log_n,d.log_source_n])
    if not np.isfinite(x).all():return None
    w=np.zeros((len(d),len(d)))
    for k in blocks:
        n=len(k)
        w[np.ix_(k,k)]=(np.eye(n)-np.ones((n,n))/n)/math.sqrt(n)
    X=w@x
    if np.linalg.matrix_rank(X)<X.shape[1] or sum(len(v)-1 for v in blocks)-X.shape[1]<4:return None
    C=np.linalg.pinv(X)
    return (C[1]-C[2])@w
def interval(known,unknown,min_count,k_lo=0,k_hi=None):
    """Exact extreme mean for fixed known natives and arbitrary subset of unresolved.
    Extremize each endpoint separately: for a given k choose k smallest/largest.
    Empty native rosters are not treated as zeros.
    """
    known=np.asarray(known,float);unknown=np.asarray(unknown,float)
    assert known.ndim==unknown.ndim==2 and known.shape[1]==unknown.shape[1]
    n,u=len(known),len(unknown)
    lo=max(0,min_count-n,k_lo);hi=u if k_hi is None else min(u,k_hi)
    if lo>hi:return None
    order=np.sort(unknown,axis=0)
    zero=np.zeros((1,known.shape[1]))
    ks=np.arange(lo,hi+1)
    small=np.vstack([zero,np.cumsum(order,axis=0)])[ks]
    large=np.vstack([zero,np.cumsum(order[::-1],axis=0)])[ks]
    denom=(n+ks)[:,None]
    assert (denom>0).all()
    return ((known.sum(0)+small)/denom).min(0),((known.sum(0)+large)/denom).max(0)
def brute_guard():
    # Exhaustive synthetic check, including empty known set, variable k and signed residual.
    import itertools
    rng=np.random.default_rng(314159)
    for n in (0,1,5):
        for u in (0,1,5):
            a=rng.normal(size=(n,3));b=rng.normal(size=(u,3))
            for minc in (1,5):
                for cap in (0,u//2,u):
                    vals=[]
                    for k in range(cap+1):
                        if n+k<minc:continue
                        for ix in itertools.combinations(range(u),k):
                            vals.append(np.vstack([a,b[list(ix)]]).mean(0))
                    q=interval(a,b,minc,k_hi=cap)
                    if not vals:assert q is None
                    else:assert np.allclose(q[0],np.min(vals,axis=0)) and np.allclose(q[1],np.max(vals,axis=0))
def fit_origin(d,tab,all_rosters,native_rosters,intro_rosters,unknown_rosters,src,ends,minimum):
    d=d.reset_index(drop=True)
    y=[all_rosters[i] for i in d.island_id];z=[src[i] for i in d.island_id]
    na=[native_rosters.get(i,set()) for i in d.island_id]
    intr=[intro_rosters.get(i,set()) for i in d.island_id]
    unk=[unknown_rosters.get(i,set()) for i in d.island_id]
    for a,n,j,k in zip(y,na,intr,unk):assert not (n&j or n&k or j&k) and a==n|j|k
    sp=sorted(set().union(*y,*z));idx={v:i for i,v in enumerate(sp)}
    t=tab.loc[sp].reset_index(drop=True);V=t[list(ends)].to_numpy(float)
    E=t.groupby(["genus","grade"],sort=False)[list(ends)].transform("mean").to_numpy(float)
    F=np.column_stack([V,V-E])
    def rows(v):return F[[idx[x] for x in sorted(v)]]
    n=np.array(list(map(len,y)));m=np.array(list(map(len,z)));nn=np.array(list(map(len,na)));nu=np.array(list(map(len,unk)))
    Ys=membership(z,idx)[0]@F
    Ya=membership(y,idx)[0]@F
    out={}
    for cohort,sel in [("known_native_min5",nn>=5),("potential_native_min5",nn+nu>=5)]:
        q=d.loc[sel].reset_index(drop=True).copy()
        for ctrl in ("all_recorded_count","known_native_count"):
            tag=cohort+"__"+ctrl
            # Fixed design: no unknown-dependent native richness adjustment.
            # For potential cohort, known-native count can be zero; skip rather than add a pseudocount.
            if ctrl=="known_native_count" and (nn[sel]<5).any():continue
            q["log_n"]=np.log2(n[sel] if ctrl=="all_recorded_count" else nn[sel])
            q["log_source_n"]=np.log2(m[sel])
            groups=[v.index.to_numpy() for _,v in q.groupby("group") if len(v)>=2]
            res={"n":len(q),"informative_groups":len(groups),"known_native_min5_islands":int((nn[sel]>=5).sum()),"median_unresolved_fraction":float(np.median(nu[sel]/n[sel])) if len(q) else None}
            if len(q)<8:
                out[tag]=res|{"status":"insufficient"};continue
            cnt=q[["other","partial","specialist"]].to_numpy(float);w=weights(q,cnt,groups)
            if w is None:
                out[tag]=res|{"status":"insufficient_rank"};continue
            all_slope=w@(Ya[sel]-Ys[sel])
            source=Ys[sel]
            known=[rows(v) for v,k in zip(na,sel) if k]
            unknown=[rows(v) for v,k in zip(unk,sel) if k]
            res["status"]="exploratory";res["max_abs_design_weight"]=float(np.max(np.abs(w)))
            res["effective_weight_islands"]=float(np.square(np.abs(w).sum())/np.square(w).sum())
            scenarios={}
            for cap in (.25,.5,.75,1.):
                bounds=[interval(a,b,5,k_hi=math.floor(cap*len(b)+1e-10)) for a,b in zip(known,unknown)]
                if any(v is None for v in bounds):
                    scenarios[str(cap)]={"status":"infeasible_min5","infeasible_islands":sum(v is None for v in bounds)};continue
                lo=np.stack([v[0] for v in bounds])-source;hi=np.stack([v[1] for v in bounds])-source
                low=np.sum(w[:,None]*np.where(w[:,None]>=0,lo,hi),axis=0)
                high=np.sum(w[:,None]*np.where(w[:,None]>=0,hi,lo),axis=0)
                assert (low<=high+1e-10).all()
                scenarios[str(cap)]={"status":"bounded","endpoints":{}}
                for typ,off in (("observed",0),("within_genus_grade_residual",len(ends))):
                    scenarios[str(cap)]["endpoints"][typ]={e:{"lower":float(low[off+j]),"upper":float(high[off+j]),"sign_identified":bool(low[off+j]>0 or high[off+j]<0)} for j,e in enumerate(ends)}
            # Two explicit status assignments; no inference of actual origin.
            potential=np.stack([np.vstack([a,b]).mean(0) for a,b in zip(known,unknown)])
            pot_slope=w@(potential-source)
            known_slope=w@(np.stack([v.mean(0) for v in known])-source) if all(len(v)>=5 for v in known) else None
            res["scenario_coefficients"]={}
            for typ,off in (("observed",0),("within_genus_grade_residual",len(ends))):
                res["scenario_coefficients"][typ]={e:{"all_recorded":float(all_slope[off+j]),"all_unknown_native":float(pot_slope[off+j]),"literal_native":float(known_slope[off+j]) if known_slope is not None else None} for j,e in enumerate(ends)}
            # Every feasible extreme must contain all-unknown-native assignment at cap=1.
            if scenarios["1.0"]["status"]=="bounded":
                for typ,off in (("observed",0),("within_genus_grade_residual",len(ends))):
                    for j,e in enumerate(ends):
                        v=scenarios["1.0"]["endpoints"][typ][e]
                        assert v["lower"]-1e-9<=pot_slope[off+j]<=v["upper"]+1e-9
                        if known_slope is not None:assert v["lower"]-1e-9<=known_slope[off+j]<=v["upper"]+1e-9
            out[tag]=res|{"scenarios":scenarios}
    return {"n_full":len(d),"cohorts":out}
def main():
    brute_guard()
    p=argparse.ArgumentParser();p.add_argument("--original-root",required=True,type=pathlib.Path);p.add_argument("--geography",required=True,type=pathlib.Path);p.add_argument("--bridge",required=True,type=pathlib.Path);p.add_argument("--bird-incidence",required=True,type=pathlib.Path);p.add_argument("--out",required=True,type=pathlib.Path);a=p.parse_args()
    r=a.original_root
    f=pd.read_csv(io.BytesIO(read_pin(r/"fixed/canonical/input/chapter1_status_flora.csv.gz","flora")),compression="gzip",dtype=str,usecols=["island_id","accepted_species","origin_status"])
    t=pd.read_csv(io.BytesIO(read_pin(r/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz","traits")),compression="gzip",dtype=str,usecols=["accepted_species","trait_name","normalized_value","evidence_scope","evidence_quality"])
    ass=pd.read_csv(io.BytesIO(read_pin(r/"fixed/source/result/island_source_assignments.csv.gz","assign")),compression="gzip",dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    mraw=(r/"fixed/source/source/gift_native_mainland_flora.csv.gz").read_bytes()
    m=pd.read_csv(io.BytesIO(mraw),compression="gzip",dtype=str,usecols=["entity_ID","work_species"])
    geo=pd.read_csv(io.BytesIO(read_pin(a.geography,"geography")),dtype={"island_id":str})
    assert f.island_id.nunique()==4505 and not f.duplicated(["island_id","accepted_species"]).any() and geo.island_id.nunique()==8264 and geo.island_id.is_unique
    assert m.entity_ID.nunique()==748 and not m.duplicated(["entity_ID","work_species"]).any()
    full=f.groupby("island_id").accepted_species.agg(set).to_dict();regions=m.groupby("entity_ID").work_species.agg(set).to_dict()
    ass=ass[ass.source_mode.isin([x for x,n in MODES])].copy();ass.source_rank=pd.to_numeric(ass.source_rank,errors="raise").astype(int)
    assert ass.source_mode.value_counts().to_dict()=={"geo50_climate10":82190,"geo_k5":41325} and not ass.duplicated(["island_id","source_mode","source_rank"]).any()
    keys=ass.sort_values("source_rank").groupby(["island_id","source_mode"]).entity_ID.agg(tuple).to_dict()
    @functools.lru_cache(maxsize=64)
    def pool(ids):return frozenset().union(*(regions.get(i,set()) for i in ids))
    trop=set(geo.loc[geo.analysis_regime.eq("tropical"),"island_id"])
    assert len(set(full)&trop)==1558
    b=cached_bird_animals(a.bridge);bird=bird_input(a.bird_incidence)
    b=b[b.island_id.isin(trop)&b.bee_id.isin(bird[3])].copy()
    assert b.bee_id.is_unique
    assert f.origin_status.isin(["native","introduced","unresolved"]).all()
    b["specialist"]=b.nectar_ge50;b["partial"]=b.nectar_gt0-b.nectar_ge50;b["other"]=b.all_birds-b.nectar_gt0
    assert (b[["specialist","partial","other"]]>=0).all().all()
    assert np.array_equal(b[["specialist","partial","other"]].sum(1),b.all_birds)
    labels,_,_,inc=bird
    for row in b.itertuples():
        cnt=np.bincount(labels[inc[row.bee_id]],minlength=3)
        assert np.array_equal(cnt,[row.other,row.partial,row.specialist])
    specs={"form_depth":("tubular","deep","tubular_and_deep"),"depth":("deep",)}
    st={x:states(t,x) for x in specs};panels={}
    for spec,ends in specs.items():
        for grade in ("high_medium","high_only"):
            tab=st[spec] if grade=="high_medium" else st[spec][st[spec].grade.eq("high")]
            scored=set(tab.index)
            f0=f[f.island_id.isin(b.island_id)&f.accepted_species.isin(scored)]
            all_rosters=f0.groupby("island_id").accepted_species.agg(set).to_dict()
            rosters={k:f0[f0.origin_status.eq(k)].groupby("island_id").accepted_species.agg(set).to_dict()
                     for k in ("native","introduced","unresolved")}
            for mode,nr in MODES:
                src={}
                for i in all_rosters:
                    ids=keys.get((i,mode),())
                    if ids:
                        assert len(ids)==nr
                        src[i]=(set(pool(ids))-full[i])&scored
                for minimum in (5,30):
                    sel={i for i,v in all_rosters.items() if len(v)>=minimum and len(src.get(i,set()))>=5}
                    d=b[b.island_id.isin(sel)].copy()
                    k=f"{spec}_{grade}_n{minimum}_{mode}"
                    panels[k]=fit_origin(d,tab,all_rosters,rosters["native"],rosters["introduced"],rosters["unresolved"],src,ends,minimum)
    assert len(panels)==16
    out={"status":"ORIGINAL_GIFT4505_NATIVE_ORIGIN_PARTIAL_IDENTIFICATION",
         "original_botanical_parent":4505,"original_tropical_parent":1558,"original_geography_parent":8264,"mainland_regions":748,
         "mainland_sha256":hashlib.sha256(mraw).hexdigest(),"source_pins":PIN|{"bird_incidence":BIRD_PIN},"panel_count":len(panels),"panels":panels,
         "warning":"Exact arithmetic decomposition on a fixed taxon pool, not causal mediation, native trait evolution or introduction effect. Literal origin statuses are source annotations; unresolved is never treated as introduced. All-to-native slope change = island selection + scored-count control change + roster composition, latter split introduced/unresolved. Genus-grade means frozen across cohorts. Low independent group support and incomplete origin records remain; no new significance claim. No protected identifiers."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    print("NATIVE_BOUND_SUMMARY",json.dumps({k:v for k,v in out.items() if k!="panels"}),flush=True)
    for k,v in panels.items():print("NATIVE_BOUND_PANEL",k,json.dumps(v),flush=True)
if __name__=="__main__":main()
