#!/usr/bin/env python3
"""Exact common-island native-status coefficient decomposition, including unknown-origin annotation."""
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
def coefficient(d,Y):
    """Same equal-archipelago design as the dual null; report only supported slopes."""
    d=d.reset_index(drop=True)
    groups=[v.index.to_numpy() for _,v in d.groupby("group") if len(v)>=2]
    cnt=d[["other","partial","specialist"]].to_numpy(float)
    b=weights(d,cnt,groups)
    if b is None:return None
    return b@Y

def means(rosters,index,V):
    a,n=membership(rosters,index)
    return a@V,n

def fit_origin(d,tab,all_rosters,native_rosters,intro_rosters,unknown_rosters,src,ends,minimum):
    # Pool and taxonomy expectation frozen across full and native subsets.
    d=d.reset_index(drop=True)
    y=[all_rosters[i] for i in d.island_id];z=[src[i] for i in d.island_id]
    na=[native_rosters.get(i,set()) for i in d.island_id]
    intr=[intro_rosters.get(i,set()) for i in d.island_id]
    unk=[unknown_rosters.get(i,set()) for i in d.island_id]
    for a,n,j,k in zip(y,na,intr,unk):
        assert not (n&j or n&k or j&k) and a==n|j|k
    sp=sorted(set().union(*y,*z));idx={v:i for i,v in enumerate(sp)}
    t=tab.loc[sp].reset_index(drop=True);V=t[list(ends)].to_numpy(float)
    E=t.groupby(["genus","grade"],sort=False)[list(ends)].transform("mean").to_numpy(float)
    # Joint vector: observed and exact within-genus-grade residual.
    F=np.column_stack([V,V-E])
    Ya,n=means(y,idx,F);Ys,m=means(z,idx,F)
    assert (n>=minimum).all() and (m>=5).all()
    native_n=np.array([len(v) for v in na])
    intro_n=np.array([len(v) for v in intr]);unknown_n=np.array([len(v) for v in unk])
    assert np.array_equal(n,native_n+intro_n+unknown_n)
    selected=native_n>=5
    d["log_n"]=np.log2(n);d["log_source_n"]=np.log2(m)
    full_beta=coefficient(d,Ya-Ys)
    n0=int(selected.sum())
    res={"full_islands":len(d),"common_native_islands":n0,"native_any_scored":int((native_n>0).sum()),
         "introduced_any_scored":int((intro_n>0).sum()),"unresolved_any_scored":int((unknown_n>0).sum()),
         "full_median_scored_native_fraction":float(np.median(native_n/n)),
         "full_median_scored_introduced_fraction":float(np.median(intro_n/n)),
         "full_median_scored_unresolved_fraction":float(np.median(unknown_n/n)),
         "native_minimum":5,"full_minimum":minimum}
    if full_beta is None:return res|{"status":"insufficient_full_rank"}
    q=d.loc[selected].reset_index(drop=True).copy()
    groups=[v.index.to_numpy() for _,v in q.groupby("group") if len(v)>=2]
    res["common_informative_groups"]=len(groups)
    res["common_informative_islands"]=sum(map(len,groups))
    if n0<8:return res|{"status":"insufficient_common_islands"}
    y0=Ya[selected];s0=Ys[selected];nn=n[selected]
    yn,_=means([v for v,k in zip(na,selected) if k],idx,F)
    # Exact per-island mixture terms. Empty status strata contribute zero weight,
    # not an imputed floral mean or a claim of ecological absence.
    def mix_component(sets):
        out=np.zeros_like(yn)
        for j,v in enumerate([v for v,k in zip(sets,selected) if k]):
            if v:
                out[j]=len(v)/nn[j]*(F[[idx[x] for x in sorted(v)]].mean(0)-yn[j])
        return out
    ci=mix_component(intr);cu=mix_component(unk)
    assert np.allclose(y0-yn,ci+cu,atol=1e-12)
    common_all=coefficient(q,y0-s0)
    qn=q.copy();qn["log_n"]=np.log2(native_n[selected])
    native_design_all=coefficient(qn,y0-s0)
    native_beta=coefficient(qn,yn-s0)
    intro_effect=coefficient(qn,ci)
    unresolved_effect=coefficient(qn,cu)
    if any(x is None for x in (common_all,native_design_all,native_beta,intro_effect,unresolved_effect)):
        return res|{"status":"insufficient_common_rank"}
    sample=common_all-full_beta
    control=native_design_all-common_all
    roster=native_beta-native_design_all
    assert np.allclose(roster,-intro_effect-unresolved_effect,atol=1e-10)
    assert np.allclose(native_beta-full_beta,sample+control+roster,atol=1e-10)
    # Common population with all-plant control retained: isolated roster effect.
    native_fixed_all_design=coefficient(q,yn-s0)
    rows={}
    for typ,offset in (("observed",0),("within_genus_grade_residual",len(ends))):
        rows[typ]={}
        for i,end in enumerate(ends):
            j=offset+i
            rows[typ][end]={k:float(v[j]) for k,v in {
                "full_all_slope":full_beta,"common_all_slope":common_all,
                "common_all_native_controls_slope":native_design_all,"common_native_slope":native_beta,
                "common_native_all_controls_slope":native_fixed_all_design,
                "island_selection_change":sample,"denominator_control_change":control,
                "roster_change":roster,"introduced_removal_contribution":-intro_effect,
                "unresolved_origin_removal_contribution":-unresolved_effect}.items()}
    # Audit whether origin annotation completeness itself tracks the bird contrast.
    # Only observed, scored name rows; no extrapolation of native status.
    cov=np.column_stack([native_n[selected]/nn,intro_n[selected]/nn,unknown_n[selected]/nn])
    covb=coefficient(q,cov)
    return res|{"status":"exploratory","common_median_native_fraction":float(np.median(cov[:,0])),
               "common_median_introduced_fraction":float(np.median(cov[:,1])),
               "common_median_unresolved_fraction":float(np.median(cov[:,2])),
               "origin_fraction_slopes":dict(zip(["native","introduced","unresolved"],covb.tolist())),
               "outcomes":rows}
def main():
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
    out={"status":"ORIGINAL_GIFT4505_COMMON_ISLAND_ORIGIN_STATUS_DECOMPOSITION",
         "original_botanical_parent":4505,"original_tropical_parent":1558,"original_geography_parent":8264,"mainland_regions":748,
         "mainland_sha256":hashlib.sha256(mraw).hexdigest(),"source_pins":PIN|{"bird_incidence":BIRD_PIN},"panel_count":len(panels),"panels":panels,
         "warning":"Exact arithmetic decomposition on a fixed taxon pool, not causal mediation, native trait evolution or introduction effect. Literal origin statuses are source annotations; unresolved is never treated as introduced. All-to-native slope change = island selection + scored-count control change + roster composition, latter split introduced/unresolved. Genus-grade means frozen across cohorts. Low independent group support and incomplete origin records remain; no new significance claim. No protected identifiers."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    print("ORIGIN_DECOMP_SUMMARY",json.dumps({k:v for k,v in out.items() if k!="panels"}),flush=True)
    for k,v in panels.items():print("ORIGIN_DECOMP_PANEL",k,json.dumps(v),flush=True)
if __name__=="__main__":main()
