#!/usr/bin/env python3
"""Original GIFT4505 matched-cohort, two-sided-of-mutualism taxonomy nulls.

Plant state labels: shared distinct plant species, genus x evidence quality.
Bird diet labels: shared distinct bird species, family or genus.
Keep both incidence matrices fixed. Separately test each conditional null and
joint independent shuffles; intersection-union gate requires all three.
No botanical/source IDs or individual rows are exported.
"""
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
def shuffle_blocks(key,labels):
    out=[]
    for k in np.unique(key):
        ids=np.flatnonzero(key==k)
        if len(ids)>1 and len(set(labels[ids]))>1:out.append(ids)
    return out
def pvalue(draw,obs):
    if len(draw)==0:return None
    return (1+int(np.sum(np.asarray(draw)>=obs-1e-12)))/(len(draw)+1)
def fit_panel(d,tab,tar,src,ends,bird,nperm,seed):
    d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=["group","bee","log_area","log_distance","log_n","log_source_n"]).reset_index(drop=True)
    groups=[v.index.to_numpy() for _,v in d.groupby("group") if len(v)>=2]
    ninfo=sum(map(len,groups));res={"n":len(d),"informative_groups":len(groups),"informative_islands":ninfo,"support_gate":len(groups)>=10 and ninfo>=30}
    if not groups or ninfo<8:return res|{"status":"insufficient"}
    y=[tar[i] for i in d.island_id];z=[src[i] for i in d.island_id]
    assert all(len(s)>=5 for s in z) and all(not a&b for a,b in zip(y,z))
    sp=sorted(set().union(*y,*z));idx={s:i for i,s in enumerate(sp)}
    A,n=membership(y,idx);M,m=membership(z,idx)
    t=tab.loc[sp].reset_index(drop=True);V=t[list(ends)].to_numpy(float)
    plants=[np.asarray(v,int) for v in t.groupby(["genus","grade"],sort=True).indices.values() if len(v)>1 and len(set(t.state.iloc[v]))>1]
    labels,fam,gen,inc=bird
    # Freeze exactly the species pool observed in this panel, not the 38-island union.
    included=np.unique(np.concatenate([inc[i] for i in d.bee_id]))
    remap={i:j for j,i in enumerate(included)}
    lab=labels[included];family=fam[included];genus=gen[included]
    W=np.zeros((len(d),len(included)))
    for j,i in enumerate(d.bee_id):W[j,[remap[k] for k in inc[i]]]=1
    eye=np.eye(3);counts=W@eye[lab]
    assert np.array_equal(counts[:,2],d.nectar_ge50) and np.array_equal(counts[:,2]+counts[:,1],d.nectar_gt0) and np.array_equal(counts.sum(1),d.all_birds)
    b=weights(d,counts,groups)
    if b is None:return res|{"status":"insufficient_rank"}
    target=A@V;source=M@V;delta=target-source
    obs=np.stack([b@target,b@delta])
    res.update(status="exploratory",plant_species=len(sp),plant_movable=sum(map(len,plants)),bird_species=len(included),contrast=obs.tolist())
    # Plant-only null draws; shared state assignment across target and source.
    rng=np.random.default_rng(seed)
    plant_draws=np.empty((nperm,2,len(ends)))
    # One draw per bird taxonomic label and per joint null, no optimization of labels.
    birds={k:shuffle_blocks(key,lab) for k,key in (("family",family),("genus",genus))}
    values={k:{"bird":[],"joint":[]} for k in birds}
    for i in range(nperm):
        v=V.copy()
        for p in plants:v[p]=V[rng.permutation(p)]
        yt=A@v;yd=yt-M@v
        plant_draws[i]=np.stack([b@yt,b@yd])
        for k,bb in birds.items():
            u=lab.copy()
            for p in bb:u[p]=lab[rng.permutation(p)]
            c=weights(d,W@eye[u],groups)
            if c is None:continue
            values[k]["bird"].append(np.stack([c@target,c@delta]))
            values[k]["joint"].append(np.stack([c@yt,c@yd]))
    out={}
    for j,endpoint in enumerate(ends):
        out[endpoint]={}
        for i,scope in enumerate(("target","source_relative")):
            actual=float(obs[i,j]);p0=pvalue(plant_draws[:,i,j],actual)
            row={"observed":actual,"plant_null":{"movable_species":sum(map(len,plants)),"draws":nperm,"p":p0,"q025_q50_q975":np.quantile(plant_draws[:,i,j],[.025,.5,.975]).tolist()}}
            for k,bb in birds.items():
                q=values[k];nv=len(q["bird"])
                v={"movable_birds":sum(map(len,bb)),"taxonomic_blocks":len(bb),"valid_draws":nv,"rank_valid_fraction":nv/nperm}
                if nv>=math.ceil(.9*nperm) and len(bb) and len(plants):
                    a=np.asarray(q["bird"])[:,i,j];z=np.asarray(q["joint"])[:,i,j]
                    p1=pvalue(a,actual);p2=pvalue(z,actual)
                    v.update(bird_only_p=p1,joint_p=p2,iut_p=max(p0,p1,p2),
                             bird_null_q025_q50_q975=np.quantile(a,[.025,.5,.975]).tolist(),
                             joint_null_q025_q50_q975=np.quantile(z,[.025,.5,.975]).tolist())
                else:v["status"]="insufficient_movable_or_rank_support"
                row[k]=v
            out[endpoint][scope]=row
    res["outcomes"]=out
    return res
def main():
    p=argparse.ArgumentParser();p.add_argument("--original-root",required=True,type=pathlib.Path);p.add_argument("--geography",required=True,type=pathlib.Path);p.add_argument("--bridge",required=True,type=pathlib.Path);p.add_argument("--bird-incidence",required=True,type=pathlib.Path);p.add_argument("--out",required=True,type=pathlib.Path);p.add_argument("--permutations",type=int,default=999);a=p.parse_args()
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
    st={x:states(t,x) for x in ENDS}
    panels={}
    for scope,nmin in SCOPES:
        f0=f[f.island_id.isin(b.island_id)]
        if scope=="original_native":f0=f0[f0.origin_status.eq("native")]
        tables={x:v if scope!="high_quality_only" else v[v.grade.eq("high")] for x,v in st.items()}
        tar={x:f0[f0.accepted_species.isin(v.index)].groupby("island_id").accepted_species.agg(set).to_dict() for x,v in tables.items()}
        for mode,nr in MODES:
            src={}
            for x,v in tables.items():
                scored=set(v.index);src[x]={}
                for i in tar[x]:
                    ids=keys.get((i,mode),())
                    if ids:
                        assert len(ids)==nr
                        src[x][i]=(set(pool(ids))-full[i])&scored
            selected={x:{i for i,s in tar[x].items() if len(s)>=nmin and len(src[x].get(i,set()))>=5} for x in tables}
            common=selected["form"]&selected["form_depth"]
            for tag,x,ids in (("form_full","form",selected["form"]),("form_common_islands","form",common),("form_depth_common_species","form_depth",common)):
                # Common-island comparisons use identical controls derived from the joint scored roster.
                ctrl="form_depth" if tag!="form_full" else "form"
                d=b[b.island_id.isin(ids)].copy()
                d["log_n"]=np.log2(d.island_id.map(lambda i:len(tar[ctrl][i])).astype(float))
                d["log_source_n"]=np.log2(d.island_id.map(lambda i:len(src[ctrl][i])).astype(float))
                k=f"{scope}_n{nmin}_{mode}_{tag}"
                seed=int.from_bytes(hashlib.sha256(k.encode()).digest()[:8],"big")
                panels[k]=fit_panel(d,tables[x],tar[x],src[x],ENDS[x],bird,a.permutations,seed)
    # Holm across all qualified nested tests for each intersection-union requirement,
    # not just favourable family or form panels. Dependence is allowed by Holm.
    tests=[v[k] for p in panels.values() for e in p.get("outcomes",{}).values() for v in e.values() for k in ("family","genus") if "iut_p" in v[k]]
    order=sorted(tests,key=lambda z:z["iut_p"]);last=0.
    for i,v in enumerate(order):
        last=max(last,min(1.,v["iut_p"]*(len(order)-i)));v["holm_iut_p"]=float(last)
    assert len(panels)==24
    out={"status":"ORIGINAL_GIFT4505_DUAL_PLANT_BIRD_TAXONOMY_NULL_EXPLORATORY","original_botanical_parent":4505,"original_tropical_parent":1558,"original_geography_parent":8264,"mainland_regions":748,
         "mainland_sha256":hashlib.sha256(mraw).hexdigest(),"source_pins":PIN|{"bird_incidence":BIRD_PIN},"permutations":a.permutations,"panel_count":len(panels),"qualified_iut_tests":len(tests),"holm_pass":sum(v["holm_iut_p"]<.05 for v in tests),"panels":panels,
         "warning":"Targeted exploratory positive contrast, not independent replication. Conditional nulls on both sides; intersection-union max p plus Holm, rank-qualified permutations. Fixed bird/plant incidence, provisional bird geography, distinct taxonomic labels; no direct pollination, causation or local phenotypic evolution. Source rosters are candidate sources, not ancestors. No individual protected identifiers."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    print("DUAL_NULL_SUMMARY",json.dumps({k:v for k,v in out.items() if k!="panels"}),flush=True)
    for k,v in panels.items():print("DUAL_NULL_PANEL",k,json.dumps(v),flush=True)
if __name__=="__main__":main()
