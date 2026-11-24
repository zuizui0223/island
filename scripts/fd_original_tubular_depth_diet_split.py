#!/usr/bin/env python3
"""Protected original Direct tubular form versus true depth and bird diet split.
Keep the full botanical parent; release only aggregate model diagnostics.
Predefined endpoints, matched cohorts, joint form/depth sensitivity, group holdouts.
"""
import argparse, functools, hashlib, io, json, math, pathlib
import numpy as np
import pandas as pd
from fd_original_tropical_pollinator_bridge import PIN,read_pin,public_animals,accepted_genus,fit_block,null_p
from fd_original_source_relative_pollinator import membership

SPECS={"form":("floral_form",),"depth":("tube_depth_class",),"form_depth":("floral_form","tube_depth_class")}
OUTCOMES={"form":("tubular",),"depth":("deep",),"form_depth":("tubular","deep","tubular_and_deep")}
FEATURES=["log_bee","log_specialist","log_partial","log_other"]
BASELINE=["log_bee","log_broad","log_other"]
SCOPES=(("all_recorded",5,"all"),("all_recorded",30,"all"),("original_native",5,"all"),("high_quality_only",5,"all"),("all_recorded",5,"tight"))
MODES=(("geo_k5",5),("geo50_climate10",10))

def cached_bird_animals(p):
    # All fields are public original native-bee model and bird checklist summaries,
    # from fd's signed candidate artifact. No botanical records in this input.
    bridge=json.loads(read_pin(p,"bridge"))
    g=pd.DataFrame(bridge["geometry"],columns=bridge["public_bridge_columns"]).rename(columns={"gshhg_geometry_hash_id":"island_id"})
    assert len(g)==234 and g.bee_id.is_unique and g.island_id.is_unique
    assert (g[["bee_overlap","gshhg_overlap"]]>=.9-1e-8).all().all()
    raw=p.with_name("fd_public_bird_model_covariates_20261124.json").read_bytes()
    assert hashlib.sha256(raw).hexdigest()=="872956d06661bc0fcc8e555beed07c910077423c59bb0d6da58d66570144892c"
    x=json.loads(raw);assert x["source_artifact"]==14901918629
    d=pd.DataFrame(x["rows"],columns=x["columns"])
    assert len(d)==52 and d.bee_id.is_unique
    b=pd.DataFrame(bridge["bird_candidate"])
    check=d.merge(b,on="bee_id",validate="one_to_one",suffixes=("_cache","_bridge"))
    assert len(check)==52
    for k in ("group","nectar_gt0","nectar_ge20","nectar_ge50","all_birds","area_ratio","match_type"):
        assert (check[k+"_cache"]==check[k+"_bridge"]).all()
    d=d.rename(columns={"native_bees":"bee","bee_area_km2":"area","dist_km":"distance","cov":"bee_cov"})
    d=g.merge(d,on="bee_id",validate="one_to_one")
    for k in ("bee","area","distance","nectar_gt0","nectar_ge50","all_birds"):
        d[k]=pd.to_numeric(d[k],errors="raise")
    assert (d.bee>=0).all() and (d.area>0).all() and (d.distance>=0).all()
    d["log_bee"]=np.log1p(d.bee);d["log_area"]=np.log2(d.area);d["log_distance"]=np.log2(1+d.distance)
    return d

def states(t,spec):
    names=SPECS[spec]
    a=t[t.trait_name.isin(names)&t.evidence_scope.eq("species_direct")&t.evidence_quality.isin(["high","medium"])&t.normalized_value.notna()].copy()
    a["normalized_value"]=a.normalized_value.str.strip()
    a=a[a.normalized_value.ne("")&~a.normalized_value.str.contains("|",regex=False)&~a.normalized_value.str.lower().isin(["none","nan","null","unknown","unresolved"])]
    assert not a.duplicated(["accepted_species","trait_name"]).any() and a.accepted_species.notna().all()
    s=a.pivot(index="accepted_species",columns="trait_name",values="normalized_value").dropna(subset=list(names))
    q=a.pivot(index="accepted_species",columns="trait_name",values="evidence_quality").reindex(s.index)
    s["grade"]=np.where(q[list(names)].eq("high").all(axis=1),"high","medium")
    s["genus"]=[accepted_genus(i) for i in s.index]
    s=s[s.genus.notna()].copy()
    s["state"]=list(map(tuple,s[list(names)].to_numpy()))
    if "floral_form" in names:s["tubular"]=s.floral_form.eq("tubular").astype(float)
    if "tube_depth_class" in names:s["deep"]=s.tube_depth_class.eq("deep").astype(float)
    if spec=="form_depth":s["tubular_and_deep"]=s.tubular*s.deep
    return s

def design(d,features,controls):
    z=fit_block(d,features+controls)
    if z is None:return None
    B,groups,dof,ninfo=z
    return B,groups,dof,ninfo

def heldout(d,Y,controls):
    """Predict WITHIN held-out group contrasts, not held-out absolute means.
    Test set centering defines the evaluation contrast; neither model fits test outcomes.
    Same eligible train/test groups for both models, no hyperparameter selection.
    """
    loss=[];count=0
    for g in sorted(d.group.unique()):
        test=d.group.eq(g).to_numpy();train=~test
        if test.sum()<2:continue
        train_df=d.loc[train].reset_index(drop=True)
        models=[]
        for fs in (BASELINE,FEATURES):
            z=design(train_df,fs,controls)
            if z is None:break
            B,gr,df,nn=z
            coef=B@Y[train]
            x=d.loc[test,fs+controls].to_numpy(float)
            models.append((x-x.mean(0))@coef)
        if len(models)!=2:continue
        yt=Y[test]-Y[test].mean(0)
        errs=[np.square(yt-m).mean(0) for m in models]
        loss.append(np.stack(errs));count+=int(test.sum())
    if not loss:return None
    arr=np.stack(loss)
    return {"groups":len(loss),"islands":count,"baseline_mse":arr[:,0].mean(0),"split_mse":arr[:,1].mean(0),"delta_mse":(arr[:,0]-arr[:,1]).mean(0),"groups_improved":(arr[:,0]>arr[:,1]).sum(0)}
def analyze(d,tab,tar,src,ends,controls,nperm):
    cols=list(dict.fromkeys(FEATURES+BASELINE+controls))
    d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=cols+["group"]).reset_index(drop=True)
    result={"n":len(d),"archipelagos":int(d.group.nunique())}
    if len(d)<8:return result|{"status":"insufficient"}
    z=design(d,FEATURES,controls)
    if z is None:return result|{"status":"insufficient_rank"}
    B,groups,dof,ninfo=z
    labels=FEATURES+["specialist_minus_partial"]
    # Same design for all outcomes, and a direct coefficient difference.
    C=np.vstack([B[:4],B[1]-B[2]])
    target=[tar[i] for i in d.island_id];source=[src[i] for i in d.island_id]
    assert all(len(s)>=5 for s in source) and all(not a&b for a,b in zip(target,source))
    species=sorted(set().union(*target,*source));idx={s:i for i,s in enumerate(species)}
    A,n=membership(target,idx);M,m=membership(source,idx)
    t=tab.loc[species].reset_index(drop=True);V=t[list(ends)].to_numpy(float)
    yt=A@V;ys=M@V;delta=yt-ys
    obs_t=C@yt;obs_s=C@ys;obs_d=C@delta
    assert np.allclose(obs_t-obs_s,obs_d)
    # Precompute linear coefficient weights on the shared distinct-species pool.
    ca=np.asarray(C@A);cm=np.asarray(C@M)
    assert ca.shape==(len(labels),len(species))
    blocks=[np.array(v,int) for v in t.groupby(["genus","grade"],sort=True).indices.values() if len(v)>1 and len(set(t.state.iloc[v]))>1]
    if not blocks:return result|{"status":"no_movable_taxa"}
    rng=np.random.default_rng(20261124)
    null_t=np.empty((nperm,len(labels),len(ends)));null_d=np.empty_like(null_t)
    for i in range(nperm):
        x=V.copy()
        for k in blocks:x[k]=V[rng.permutation(k)]
        null_t[i]=ca@x;null_d[i]=(ca-cm)@x
    omit_t=[];omit_d=[]
    for g in groups:
        keep=d.group.ne(g).to_numpy();other=design(d.loc[keep].reset_index(drop=True),FEATURES,controls)
        if other is None:continue
        b=other[0];c=np.vstack([b[:4],b[1]-b[2]])
        omit_t.append(c@yt[keep]);omit_d.append(c@delta[keep])
    ht=heldout(d,yt,controls);hd=heldout(d,delta,controls)
    rows={}
    for j,x in enumerate(labels):
        rows[x]={}
        for k,end in enumerate(ends):
            row={"target_slope":float(obs_t[j,k]),"source_slope":float(obs_s[j,k]),"delta_slope":float(obs_d[j,k])}
            for name,obs,sim,om in [("target",obs_t,null_t,omit_t),("delta",obs_d,null_d,omit_d)]:
                v=sim[:,j,k];b=float(obs[j,k]);vals=[float(q[j,k]) for q in om]
                row[name]={"null_mean":float(v.mean()),"observed_minus_null":float(b-v.mean()),"p_two_sided":null_p(v,b),"null_q025_q975":np.quantile(v,[.025,.975]).tolist(),"leave_group_range":[min(vals),max(vals)] if vals else None,"leave_group_n":len(vals)}
            rows[x][end]=row
    def holdout_result(v):
        if v is None:return None
        return {k:{end:float(x[j]) if k!="groups_improved" else int(x[j]) for j,end in enumerate(ends)} if isinstance(x,np.ndarray) else x for k,x in v.items()}
    return result|{"status":"exploratory","informative_groups":len(groups),"informative_islands":ninfo,"dof":dof,"primary_support":len(groups)>=10 and ninfo>=30,
        "distinct_species":len(species),"movable_species":sum(map(len,blocks)),
        "median_target_scored":float(np.median(n)),"median_source_scored":float(np.median(m)),
        "mean_source_relative":dict(zip(ends,delta.mean(0).tolist())),
        "heldout_within_group":{"target":holdout_result(ht),"delta":holdout_result(hd)},
        "outcomes":rows}

def main():
    p=argparse.ArgumentParser();p.add_argument("--original-root",required=True,type=pathlib.Path);p.add_argument("--geography",required=True,type=pathlib.Path);p.add_argument("--bridge",required=True,type=pathlib.Path);p.add_argument("--out",required=True,type=pathlib.Path);p.add_argument("--permutations",type=int,default=399);a=p.parse_args()
    r=a.original_root
    f=pd.read_csv(io.BytesIO(read_pin(r/"fixed/canonical/input/chapter1_status_flora.csv.gz","flora")),compression="gzip",dtype=str,usecols=["island_id","accepted_species","origin_status"])
    t=pd.read_csv(io.BytesIO(read_pin(r/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz","traits")),compression="gzip",dtype=str,usecols=["accepted_species","trait_name","normalized_value","evidence_scope","evidence_quality"])
    ass=pd.read_csv(io.BytesIO(read_pin(r/"fixed/source/result/island_source_assignments.csv.gz","assign")),compression="gzip",dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    raw=(r/"fixed/source/source/gift_native_mainland_flora.csv.gz").read_bytes();mainland_sha=hashlib.sha256(raw).hexdigest()
    m=pd.read_csv(io.BytesIO(raw),compression="gzip",dtype=str,usecols=["entity_ID","work_species"])
    geo=pd.read_csv(io.BytesIO(read_pin(a.geography,"geography")),dtype={"island_id":str})
    assert f.island_id.nunique()==4505 and not f.duplicated(["island_id","accepted_species"]).any() and geo.island_id.nunique()==8264 and geo.island_id.is_unique
    assert m.entity_ID.nunique()==748 and not m.duplicated(["entity_ID","work_species"]).any()
    full=f.groupby("island_id").accepted_species.agg(set).to_dict();rg=m.groupby("entity_ID").work_species.agg(set).to_dict()
    ass=ass[ass.source_mode.isin([i for i,n in MODES])].copy()
    ass.source_rank=pd.to_numeric(ass.source_rank,errors="raise").astype(int)
    assert not ass.duplicated(["island_id","source_mode","source_rank"]).any()
    assert ass.source_mode.value_counts().to_dict()=={"geo50_climate10":82190,"geo_k5":41325}
    key=ass.sort_values("source_rank").groupby(["island_id","source_mode"]).entity_ID.agg(tuple).to_dict()
    @functools.lru_cache(maxsize=64)
    def pool(ids):return frozenset().union(*(rg.get(i,set()) for i in ids))
    tropical=set(geo.loc[geo.analysis_regime.eq("tropical"),"island_id"])
    assert len(set(full)&tropical)==1558
    b=cached_bird_animals(a.bridge)
    b=b[b.island_id.isin(tropical)&b.nectar_gt0.notna()].copy()
    b["specialist"]=b.nectar_ge50;b["partial"]=b.nectar_gt0-b.nectar_ge50;b["other"]=b.all_birds-b.nectar_gt0
    assert (b[["specialist","partial","other"]]>=0).all().all()
    assert np.allclose(b[["specialist","partial","other"]].sum(axis=1),b.all_birds)
    for k in ("specialist","partial","other"):b["log_"+k]=np.log1p(b[k])
    b["log_broad"]=np.log1p(b.nectar_gt0)
    b["area_ratio"]=pd.to_numeric(b.area_ratio,errors="coerce")
    panels={}
    for spec in SPECS:
        st=states(t,spec)
        for scope,nmin,cohort in SCOPES:
            a0=b if cohort=="all" else b[b.area_ratio.between(.9,1.1)&b.match_type.eq("normalized_name")]
            fs=f[f.island_id.isin(a0.island_id)]
            if scope=="original_native":fs=fs[fs.origin_status.eq("native")]
            tab=st if scope!="high_quality_only" else st[st.grade.eq("high")]
            scored=set(tab.index)
            tar=fs[fs.accepted_species.isin(scored)].groupby("island_id").accepted_species.agg(set).to_dict()
            for mode,nr in MODES:
                src={}
                for i in tar:
                    ids=key.get((i,mode),())
                    if ids:
                        assert len(ids)==nr
                        src[i]=(set(pool(ids))-full[i])&scored
                sel=[i for i,v in tar.items() if len(v)>=nmin and len(src.get(i,set()))>=5]
                d=a0[a0.island_id.isin(sel)].copy()
                d["log_n"]=np.log2(d.island_id.map(lambda x:len(tar[x])))
                d["log_source_n"]=np.log2(d.island_id.map(lambda x:len(src[x])))
                controls=["log_area","log_distance","log_n","log_source_n"]
                name=f"{spec}_{scope}_n{nmin}_{cohort}_{mode}"
                panels[name]=analyze(d,tab,tar,src,OUTCOMES[spec],controls,a.permutations)
    tests=[row[k] for p in panels.values() for ys in p.get("outcomes",{}).values() for row in ys.values() for k in ("target","delta")]
    order=sorted(tests,key=lambda x:x["p_two_sided"]);last=1.
    for i in range(len(order)-1,-1,-1):
        last=min(last,order[i]["p_two_sided"]*len(order)/(i+1));order[i]["BH_q_exploratory"]=float(last)
    assert len(panels)==30
    out={"status":"ORIGINAL_GIFT4505_DIRECT_TUBULAR_DEPTH_DIET_SPLIT_EXPLORATORY","original_botanical_parent":4505,"original_tropical_parent":1558,"original_geography_parent":8264,"mainland_regions":748,"mainland_sha256":mainland_sha,"source_pins":PIN|{"public_bird_model":"872956d06661bc0fcc8e555beed07c910077423c59bb0d6da58d66570144892c"},
         "panel_count":len(panels),"test_count":len(tests),"panels":panels,"permutations":a.permutations,
         "warning":"Single-state Direct form and depth separate; common-species form-depth is sensitivity not local measured flowers. Target and source-relative contrast; candidate source roster is not historical ancestor. Fixed genus×grade randomization preserves incidence, taxonomy and documentation, not causal identification. Held-out prediction evaluates centered within-group differences, not absolute unseen-island outcomes. Sparse groups; no confirmatory inference. Bird crosswalk provisional. No individual botanical/geographic identifiers exported."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    for k,v in panels.items():print("TUBULAR_DEPTH_PANEL",k,json.dumps(v),flush=True)
    print("TUBULAR_DEPTH_COMPLETE",json.dumps({"panels":len(panels),"tests":len(tests),"bh_pass":sum(x["BH_q_exploratory"]<.05 for x in tests)}),flush=True)
if __name__=="__main__":main()
