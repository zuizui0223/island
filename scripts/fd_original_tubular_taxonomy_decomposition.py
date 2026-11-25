#!/usr/bin/env python3
"""Exact plant taxonomy decomposition and name-disjoint mainland donor audit on the original dual-null cohorts."""
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
def decompose(b,observed,expected,ends):
    """Exact linear coefficient decomposition; no inference of evolutionary mediation."""
    o=b@observed;e=b@expected;r=b@(observed-expected)
    assert np.allclose(o,e+r,rtol=1e-9,atol=1e-10)
    return {end:{"observed":float(o[i]),"taxonomy_expected":float(e[i]),"within_block_residual":float(r[i])}
            for i,end in enumerate(ends)}

def fit_panel(d,tab,tar,src,ends,bird,nperm,seed,full):
    d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=["group","bee","log_area","log_distance","log_n","log_source_n"]).reset_index(drop=True)
    groups=[v.index.to_numpy() for _,v in d.groupby("group") if len(v)>=2]
    ninfo=sum(map(len,groups))
    res={"n":len(d),"informative_groups":len(groups),"informative_islands":ninfo,"support_gate":len(groups)>=10 and ninfo>=30}
    if not groups or ninfo<8:return res|{"status":"insufficient"}
    target=[tar[i] for i in d.island_id];source=[src[i] for i in d.island_id]
    assert all(len(s)>=5 for s in source) and all(not a&b for a,b in zip(target,source))
    sp=sorted(set().union(*target,*source));idx={s:i for i,s in enumerate(sp)}
    A,n=membership(target,idx);M,m=membership(source,idx)
    t=tab.loc[sp].reset_index(drop=True);V=t[list(ends)].to_numpy(float)
    labels,fam,gen,inc=bird
    included=np.unique(np.concatenate([inc[i] for i in d.bee_id]))
    remap={i:j for j,i in enumerate(included)}
    lab=labels[included];W=np.zeros((len(d),len(included)))
    for j,i in enumerate(d.bee_id):W[j,[remap[k] for k in inc[i]]]=1
    counts=W@np.eye(3)[lab]
    assert np.array_equal(counts[:,2],d.nectar_ge50) and np.array_equal(counts[:,1]+counts[:,2],d.nectar_gt0) and np.array_equal(counts.sum(1),d.all_birds)
    b=weights(d,counts,groups)
    if b is None:return res|{"status":"insufficient_rank"}
    yt=A@V;ys=M@V;delta=yt-ys
    # The exact expectation of a uniform within-block shuffle is each block mean,
    # with the observed distinct-species source+target roster frozen.
    out={}
    for tag,cols in (("genus",["genus"]),("genus_grade",["genus","grade"]),("grade",["grade"])):
        expected=t.groupby(cols,sort=False)[list(ends)].transform("mean").to_numpy(float)
        et=A@expected;es=M@expected
        out[tag]={"target":decompose(b,yt,et,ends),"source":decompose(b,ys,es,ends),
                  "source_relative":decompose(b,delta,et-es,ends)}
        if tag=="genus_grade":
            # Verify Monte Carlo mean against analytic expectation (same fixed pool);
            # do not run a new direction-selected test.
            blocks=[np.array(v,int) for v in t.groupby(cols,sort=True).indices.values() if len(v)>1 and len(set(t.state.iloc[v]))>1]
            rng=np.random.default_rng(seed);draws=[]
            wA=np.asarray(b@A).ravel();wM=np.asarray(b@M).ravel()
            for i in range(nperm):
                v=V.copy()
                for g in blocks:v[g]=V[rng.permutation(g)]
                draws.append(np.stack([wA@v,(wA-wM)@v]))
            sims=np.stack(draws)
            expected_coef=np.stack([b@et,b@(et-es)])
            mean=sims.mean(0);se=sims.std(0,ddof=1)/math.sqrt(nperm)
            assert np.all(np.abs(mean-expected_coef)<=6*se+1e-9)
            out[tag]["mc_replay"]={"draws":nperm,"max_abs_mean_error":float(np.max(np.abs(mean-expected_coef))),
                "movable_species":sum(map(len,blocks))}
    # Independent donor check: exclude EVERY recorded plant from EVERY matched
    # target island, even unscored species or other origin statuses.
    # Predict only target species with >=2 distinct candidate-mainland donors
    # in the same genus(/grade). This is NOT a prediction of unscored traits.
    banned=set().union(*(full[i] for i in d.island_id))
    donor_set=set().union(*source)-banned
    assert not donor_set&banned
    donor=tab.loc[sorted(donor_set)].copy()
    donor_results={}
    for tag,cols in (("genus",["genus"]),("genus_grade",["genus","grade"])):
        if donor.empty:
            donor_results[tag]={"status":"no_name_disjoint_donors"};continue
        g=donor.groupby(cols,sort=False);mu=g[list(ends)].mean();mu=mu.loc[g.size().ge(2)]
        target_table=tab.loc[sp].reset_index(drop=True)
        keys=pd.Index(target_table[cols[0]]) if len(cols)==1 else pd.MultiIndex.from_frame(target_table[cols])
        external=mu.reindex(keys).to_numpy(float)
        matched=np.isfinite(external).all(1)
        valid_target=[{s for s in v if matched[idx[s]]} for v in target]
        q=np.array([len(v) for v in valid_target]);coverage=q/n
        eligible=(q>=5)&(coverage>=.5)
        rr={"donor_distinct_species":len(donor_set),"donor_blocks_ge2":len(mu),"median_target_coverage":float(np.median(coverage)),
            "n_eligible":int(eligible.sum()),"n_any_coverage":int((q>0).sum())}
        dd=d.loc[eligible].reset_index(drop=True)
        if len(dd)<8:
            donor_results[tag]=rr|{"status":"insufficient_coverage"};continue
        gg=[v.index.to_numpy() for _,v in dd.groupby("group") if len(v)>=2]
        be=weights(dd,counts[eligible],gg)
        if be is None:
            donor_results[tag]=rr|{"status":"insufficient_rank","informative_groups":len(gg)};continue
        Q,_=membership([v for v,k in zip(valid_target,eligible) if k],idx)
        predicted=np.where(np.isfinite(external),external,0.)
        # No source outcome is evaluated against the donor pool used to train means.
        donor_results[tag]=rr|{"status":"exploratory","informative_groups":len(gg),
            "target_only_common_donor_support":decompose(be,Q@V,Q@predicted,ends)}
    return res|{"status":"exploratory","distinct_scored_species":len(sp),"exact_decomposition":out,"external_name_disjoint_donor":donor_results}
def main():
    p=argparse.ArgumentParser();p.add_argument("--original-root",required=True,type=pathlib.Path);p.add_argument("--geography",required=True,type=pathlib.Path);p.add_argument("--bridge",required=True,type=pathlib.Path);p.add_argument("--bird-incidence",required=True,type=pathlib.Path);p.add_argument("--out",required=True,type=pathlib.Path);p.add_argument("--permutations",type=int,default=199);a=p.parse_args()
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
                panels[k]=fit_panel(d,tables[x],tar[x],src[x],ENDS[x],bird,a.permutations,seed,full)
    assert len(panels)==24
    out={"status":"ORIGINAL_GIFT4505_FLORAL_DIET_EXACT_TAXONOMY_DECOMPOSITION","original_botanical_parent":4505,"original_tropical_parent":1558,"original_geography_parent":8264,"mainland_regions":748,
         "mainland_sha256":hashlib.sha256(mraw).hexdigest(),"source_pins":PIN|{"bird_incidence":BIRD_PIN},"panel_count":len(panels),"permutations_for_analytic_replay":a.permutations,
         "panels":panels,"warning":"Exact conditional arithmetic decomposition of observed cross-sectional diet slopes into genus/grade expectation and residual, not causal mediation or evolutionary change. External mainland donor means strictly exclude every matched island's whole flora; only scored target taxa with donor support evaluated, no imputation or independent prediction of original source outcomes. Common island controls fixed; sparse groups, no new significance claim. No protected identifiers."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    print("TAXONOMY_DECOMP_SUMMARY",json.dumps({k:v for k,v in out.items() if k!="panels"}),flush=True)
    for k,v in panels.items():print("TAXONOMY_DECOMP_PANEL",k,json.dumps(v),flush=True)
if __name__=="__main__":main()
