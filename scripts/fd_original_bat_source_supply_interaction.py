#!/usr/bin/env python3
"""Post-pilot source supply x bat opportunity interaction, protected GIFT4505.

The donor is a contemporary candidate mainland pool, not an ancestor.
Source and target trait labels are permuted coherently on their union;
the interaction design is rebuilt at every draw. Shared species are never
independently permuted. Full and island-name-excluded donor pools use the
same eligible islands and the same target rosters. No individual IDs exit.
"""
import argparse, functools, hashlib, io, json, math, pathlib, tempfile
import numpy as np
import pandas as pd
from fd_original_bat_genus_null import (PIN,read_pin,states,animal_overlay,TRAITS,SEED,BASE_CONTROLS)
from fd_original_source_relative_pollinator import membership

MODES=(("geo_k5",5),("geo50_climate10",10))
ENDPOINTS=("tubular","deep","white")
SCOPES=("all_recorded","original_native","high_quality_only")
CONTROLS=BASE_CONTROLS+["log_unknown","log_source_n","source_coverage"]
PRIMARY="tubular_all_recorded_geo_k5_disjoint"

def fit_design(fixed,supply,groups,Y,weights=None):
    # Supply is in ten-percentage-point units, with both diet class main
    # effects AND their interactions, so no constant-total substitution claim.
    x=np.column_stack([fixed[:,0],fixed[:,1],supply,
                       fixed[:,0]*supply,fixed[:,1]*supply,fixed[:,2:]])
    if weights is None:weights=np.ones(len(groups),int)
    xs=[];ys=[];n=0;k=0
    for ids,m in zip(groups,weights):
        if not m:continue
        z=x[ids];y=Y[ids];w=math.sqrt(m/len(ids))
        xs.append((z-z.mean(0))*w);ys.append((y-y.mean(0))*w)
        n+=len(ids);k+=1
    p=x.shape[1]
    if not xs or n-k-p<4:return None
    X=np.vstack(xs);yy=np.vstack(ys)
    if np.linalg.matrix_rank(X)<p:return None
    b=np.linalg.lstsq(X,yy,rcond=None)[0]
    return b[3]-b[4]

def estimate(d,tar,src,tab,perms=999,boots=1999):
    cols=["log_focal","log_other",*CONTROLS]
    d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=cols+["group"]).reset_index(drop=True)
    info=[g for g,v in d.groupby("group",sort=True) if len(v)>=2 and v.log_focal.nunique()>1]
    result={"eligible_islands":len(d),"eligible_groups":int(d.group.nunique()),
            "informative_groups_diet":len(info),"controls":CONTROLS}
    d=d[d.group.isin(info)].reset_index(drop=True)
    result["informative_islands_diet"]=len(d)
    if len(d)<12 or tab.y.nunique()<2:return result|{"status":"not_estimable"}
    rs=[tar[i] for i in d.island_id];ms=[src[i] for i in d.island_id]
    tab=tab.sort_index();idx={s:i for i,s in enumerate(tab.index)}
    A,n=membership(rs,idx);M,m=membership(ms,idx)
    assert (n>=5).all() and (m>=5).all()
    V=tab.y.to_numpy(float)
    E=tab.groupby(["genus","grade"],sort=True).y.transform("mean").to_numpy(float)
    codes=tab.groupby(["genus","grade"],sort=True).ngroup().to_numpy()
    movable=tab.groupby(["genus","grade"],sort=True).y.transform("nunique").gt(1).to_numpy()
    fixed=d[cols].to_numpy(float)
    groups=[v.index.to_numpy() for _,v in d.groupby("group",sort=True)]
    K=len(groups)
    exp_tar=np.asarray(A@E).ravel();exp_src=np.asarray(M@E).ravel()
    def response(values):
        y=np.asarray(A@values).ravel()
        s=np.asarray(M@values).ravel()
        # raw target and target-minus-source interactions are algebraically
        # identical when source share is included in the design.
        # residual delta additionally removes both genus-grade expectations.
        raw=y-s;expected=exp_tar-exp_src
        return 10*s,np.column_stack([y,raw,expected,raw-expected])
    supply,Y=response(V)
    point=fit_design(fixed,supply,groups,Y)
    joint_groups=[ids for ids in groups if np.ptp(supply[ids])>1e-9]
    result|={"informative_groups_diet_and_supply":len(joint_groups),
        "informative_islands_diet_and_supply":sum(map(len,joint_groups)),
        "source_share_mean":float(supply.mean()/10),"source_share_sd":float(supply.std()/10),
        "target_share_mean":float(Y[:,0].mean()),"target_minus_source_mean":float(Y[:,1].mean()),
        "median_target_scored":float(np.median(n)),"median_source_scored":float(np.median(m)),
        "shared_target_source_names_median":float(np.median([len(a&b) for a,b in zip(rs,ms)])),
        "movable_target_taxa":int((movable & (np.asarray(A.sum(0)).ravel()>0)).sum()),
        "movable_source_taxa":int((movable & (np.asarray(M.sum(0)).ravel()>0)).sum())}
    if point is None:return result|{"status":"not_estimable"}
    assert np.isclose(point[0],point[1],atol=1e-9) and np.isclose(point[1]-point[2],point[3],atol=1e-9)
    rng=np.random.default_rng(SEED)
    base=np.argsort(codes,kind="stable")
    null=[]
    for _ in range(perms):
        order=np.lexsort((rng.random(len(V)),codes))
        perm=np.empty_like(V);perm[base]=V[order]
        ss,yy=response(perm)
        q=fit_design(fixed,ss,groups,yy)
        if q is not None:null.append(q)
    boot=[];leave=[]
    for _ in range(boots):
        q=fit_design(fixed,supply,groups,Y,np.bincount(rng.integers(0,K,K),minlength=K))
        if q is not None:boot.append(q)
    for k in range(K):
        w=np.ones(K,int);w[k]=0
        q=fit_design(fixed,supply,groups,Y,w)
        if q is not None:leave.append(q)
    names=("target","source_relative","genus_grade_expected_delta","beyond_genus_grade_delta")
    intervals=np.quantile(np.asarray(boot),[.025,.975],axis=0).T.tolist() if len(boot)>=.9*boots else [None]*4
    lrange=np.column_stack([np.min(leave,0),np.max(leave,0)]).tolist() if leave else [None]*4
    tests={}
    for j,name in enumerate(names):
        sim=np.asarray(null)[:,j] if len(null) else np.array([])
        p=(1+int((np.abs(sim)>=abs(point[j])-1e-12).sum()))/(len(sim)+1) if len(sim)>=.9*perms else None
        tests[name]={"interaction_contrast":float(point[j]),"bootstrap_interval":intervals[j],
                     "leave_group_range":lrange[j],"joint_null_two_sided_p":p,
                     "joint_null_quantiles":np.quantile(sim,[.025,.5,.975]).tolist() if len(sim) else None}
    gate=(len(joint_groups)>=30 and sum(map(len,joint_groups))>=200 and
          len(boot)>=.9*boots and len(leave)==K and len(null)>=.9*perms and
          result["movable_target_taxa"]>=30 and result["movable_source_taxa"]>=30)
    return result|{"status":"exploratory_estimable" if gate else "insufficient_support",
      "valid_bootstrap":len(boot),"valid_leave_group":len(leave),"valid_joint_permutations":len(null),
      "results":tests}

def holm(panels):
    def p(k):
        return panels[k].get("results",{}).get("beyond_genus_grade_delta",{}).get("joint_null_two_sided_p") or 1.
    last=0.
    for j,k in enumerate(sorted(panels,key=p)):
        v=panels[k];last=max(last,min(1.,p(k)*(len(panels)-j)))
        v["holm_p_joint_residual"]=float(last)
        a=v.get("results",{}).get("beyond_genus_grade_delta",{})
        ci=a.get("bootstrap_interval");lr=a.get("leave_group_range")
        v["supported_positive_interaction"]=bool(v["status"]=="exploratory_estimable" and last<.05 and ci and ci[0]>0 and lr and lr[0]>0)

def prepare(flora,t,geo,animal,assign,mainland):
    full=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    region=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    keys=assign.sort_values("source_rank").groupby(["island_id","source_mode"]).entity_ID.agg(tuple).to_dict()
    @functools.lru_cache(maxsize=256)
    def pool(ids):return frozenset().union(*(region[k] for k in ids))
    trop=geo[geo.analysis_regime.eq("tropical")].copy()
    f=flora[flora.island_id.isin(trop.island_id)].copy()
    assert f.island_id.nunique()==1558
    # All candidate source regions linked to the tropical botanical parent,
    # independently of animal eligibility and trait outcomes.
    tropical_regions=set(assign.loc[assign.island_id.isin(full.keys()&set(trop.island_id)),"entity_ID"])
    donor_names=set().union(*(region[k] for k in tropical_regions if k in region))
    source_census={'tropical_candidate_regions':len(tropical_regions),'candidate_regions_without_native_flora':len(tropical_regions-set(region))}
    print('ORIGINAL_BAT_SUPPLY_DONOR_CENSUS',json.dumps(source_census),flush=True)
    d=animal.merge(trop,on="island_id",validate="one_to_one")
    d=d[d.total.gt(0)&d.range_area.ge(20)&d.distance_to_continent_km.ge(0)&d.match_share.ge(.5)].copy()
    d=d.dropna(subset=["island_latitude","island_longitude"])
    d["group"]=(np.floor(d.island_latitude/10).astype("Int64").astype(str)+"_"+np.floor(d.island_longitude/10).astype("Int64").astype(str))
    d["log_area"]=np.log2(d.range_area);d["log_distance"]=np.log2(1+d.distance_to_continent_km);d["abs_lat"]=d.island_latitude.abs()
    d["log_focal"]=np.log1p(d.nectar);d["log_other"]=np.log1p(d.other);d["log_unknown"]=np.log1p(d.unknown)
    panels={}
    for end in ENDPOINTS:
        for scope in SCOPES:
            fs=f[f.origin_status.eq("native")] if scope=="original_native" else f
            tab=states(t,*TRAITS[end],scope)
            tab=tab[tab.index.isin(set(fs.accepted_species)|donor_names)].copy()
            scored=set(tab.index)
            rosters=fs[fs.accepted_species.isin(scored)].groupby("island_id").accepted_species.agg(set).to_dict()
            total=fs.groupby("island_id").accepted_species.nunique()
            for mode,nr in MODES:
                incl={};excl={};ni={};ne={};missing_assign=missing_donors=0
                for i in set(d.island_id)&set(rosters):
                    ids=keys.get((i,mode),())
                    if not ids:
                        missing_assign+=1;continue
                    assert len(ids)==nr
                    if any(k not in region for k in ids):
                        missing_donors+=1;continue
                    source=pool(ids)
                    separate=source-full[i]
                    incl[i]=source&scored;excl[i]=separate&scored
                    ni[i]=len(source);ne[i]=len(separate)
                    assert not (rosters[i]&excl[i])
                eligible=[i for i in incl if len(rosters[i])>=5 and len(excl[i])>=5]
                q=d[d.island_id.isin(eligible)].copy()
                q["log_n"]=np.log2(q.island_id.map(lambda i:len(rosters[i])))
                q["direct_coverage"]=q.island_id.map(lambda i:len(rosters[i]))/q.island_id.map(total)
                # Same complete covariates and islands in both donor variants.
                q=q.replace([np.inf,-np.inf],np.nan).dropna(subset=[x for x in CONTROLS if x not in ("log_source_n","source_coverage")]+["log_focal","log_other","group"])
                for tag,src,denom in (("inclusive",incl,ni),("disjoint",excl,ne)):
                    z=q.copy()
                    z["log_source_n"]=np.log2(z.island_id.map(lambda i:len(src[i])))
                    z["source_coverage"]=z.island_id.map(lambda i:len(src[i]))/z.island_id.map(denom)
                    k=f"{end}_{scope}_{mode}_{tag}"
                    panels[k]=estimate(z,rosters,src,tab)
                    panels[k]['animal_and_trait_islands_without_source_assignment']=missing_assign
                    panels[k]['animal_and_trait_islands_with_missing_candidate_flora']=missing_donors
                x=panels[f"{end}_{scope}_{mode}_inclusive"];y=panels[f"{end}_{scope}_{mode}_disjoint"]
                assert x["eligible_islands"]==y["eligible_islands"] and x["informative_islands_diet"]==y["informative_islands_diet"]
    assert len(panels)==36
    holm(panels)
    return panels,source_census

def self_test():
    # Verify scale, interaction recovery and source-relative identity exactly.
    rng=np.random.default_rng(SEED);n=240
    x=rng.normal(size=(n,3));s=rng.uniform(0,10,n)
    g=[np.arange(i,i+6) for i in range(0,n,6)]
    y=.3*x[:,0]*s-.1*x[:,1]*s+.2*s+.4*x[:,2]
    v=fit_design(x,s,g,np.column_stack([y,y-s/10]))
    assert np.allclose(v,[.4,.4],atol=1e-10)
    # Species-coherent null and common donor construction tested on synthetic pool.
    ids=[f"Taxon {i}" for i in range(120)]
    tab=pd.DataFrame({"genus":np.repeat(np.arange(30),4),"grade":"high","y":np.tile([0.,0.,1.,1.],30)},index=ids)
    d=pd.DataFrame({"group":np.repeat(np.arange(40),6),"island_id":np.arange(n),"log_focal":x[:,0],"log_other":x[:,1]})
    for c in CONTROLS:d[c]=rng.normal(size=n)
    rr={};ss={}
    for i in range(n):
        a=rng.choice(ids,40,replace=False);rr[i]=set(a[:20]);ss[i]=set(a[20:])
    q=estimate(d,rr,ss,tab,perms=99,boots=99)
    assert q["status"]=="exploratory_estimable" and q["valid_joint_permutations"]>=90
    assert abs(q["results"]["target"]["interaction_contrast"]-q["results"]["source_relative"]["interaction_contrast"])<1e-8
    print("SOURCE_SUPPLY_INTERACTION_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--self-test",action="store_true")
    p.add_argument("--original-root",type=pathlib.Path);p.add_argument("--geography",type=pathlib.Path)
    p.add_argument("--fd-root",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);a=p.parse_args()
    if a.self_test:self_test();return
    r=a.original_root
    f=pd.read_csv(io.BytesIO(read_pin(r/"fixed/canonical/input/chapter1_status_flora.csv.gz","flora")),compression="gzip",dtype=str,usecols=["island_id","accepted_species","origin_status"])
    t=pd.read_csv(io.BytesIO(read_pin(r/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz","traits")),compression="gzip",dtype=str,usecols=["accepted_species","trait_name","normalized_value","evidence_scope","evidence_quality"])
    geo=pd.read_csv(io.BytesIO(read_pin(a.geography,"geography")),dtype={"island_id":str})
    ass=pd.read_csv(io.BytesIO(read_pin(r/"fixed/source/result/island_source_assignments.csv.gz","assign")),compression="gzip",dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    raw=(r/"fixed/source/source/gift_native_mainland_flora.csv.gz").read_bytes()
    mainland_sha=hashlib.sha256(raw).hexdigest()
    mainland=pd.read_csv(io.BytesIO(raw),compression="gzip",dtype=str,usecols=["entity_ID","work_species"])
    assert f.island_id.nunique()==4505 and not f.duplicated(["island_id","accepted_species"]).any()
    assert geo.island_id.nunique()==8264 and geo.island_id.is_unique
    assert mainland.entity_ID.nunique()==748 and not mainland.duplicated(["entity_ID","work_species"]).any()
    ass=ass[ass.source_mode.isin([k for k,n in MODES])].copy()
    ass.source_rank=pd.to_numeric(ass.source_rank,errors="raise").astype(int)
    assert not ass.duplicated(["island_id","source_mode","source_rank"]).any()
    assert ass.source_mode.value_counts().to_dict()=={"geo50_climate10":82190,"geo_k5":41325}
    with tempfile.TemporaryDirectory() as td:
        exact,_,source=animal_overlay(a.fd_root,geo,pathlib.Path(td))
    panels,source_census=prepare(f,t,geo,exact,ass,mainland)
    out={"original_botanical_parent":4505,"original_geography_parent":8264,"original_tropical_parent":1558,
         "mainland_candidate_regions":748,"mainland_sha256":mainland_sha,"pinned_sources":PIN,
         "animal_source":source,"source_census":source_census,"primary_post_pilot":PRIMARY,"panels":panels,
         "warning":"All analyses post-pilot exploratory. Source trait availability is from contemporary geographic/climatic candidate mainland species, not ancestors or verified colonists. Disjoint pool depends on focal flora membership; inclusive variant retains shared names on identical cohort. Model is conditional association, not selective loss, historical adaptation or causal pollination. Interaction compares log-count slopes per ten percentage points of source trait supply, not a constant-total guild replacement. Joint null recomputes source and target states coherently across species; no plant individual or geography identifiers exported."}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    for k,v in panels.items():print("ORIGINAL_BAT_SUPPLY_PANEL",k,json.dumps(v,allow_nan=False),flush=True)
    print("ORIGINAL_BAT_SUPPLY_COMPLETE",json.dumps({"panels":len(panels),"support_qualified":sum(v["status"]=="exploratory_estimable" for v in panels.values()),"positive_residual":sum(v["supported_positive_interaction"] for v in panels.values())}),flush=True)
if __name__=="__main__":main()
