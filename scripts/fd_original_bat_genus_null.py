#!/usr/bin/env python3
"""Original protected GIFT4505 bat-opportunity / plant genus-grade test.
Only aggregate counts, coefficients, intervals and null statistics are exported.
Public pilot results were inspected before this design; all tests are exploratory.
"""
import argparse,hashlib,io,json,math,pathlib,sys,tempfile
import numpy as np
import pandas as pd
from fd_original_tropical_pollinator_bridge import PIN,read_pin,accepted_genus
from fd_original_source_relative_pollinator import membership

SEED=20261125
TRAITS={"white":("flower_primary_color","white"),
        "tubular":("floral_form","tubular"),
        "deep":("tube_depth_class","deep")}
SCOPES=(("all_recorded",5,.5),("all_recorded",5,.8),("all_recorded",30,.5),
        ("original_native",5,.5),("high_quality_only",5,.5))
FAMILIES=("pteropodidae","phyllostomidae")
BASE_CONTROLS=["log_area","log_distance","log_n","abs_lat",*[f"climate_pc{i}" for i in range(1,5)],"direct_coverage"]

def states(t,trait,target,scope):
    a=t[t.trait_name.eq(trait)&t.evidence_scope.eq("species_direct")&
        t.evidence_quality.isin(["high","medium"])&t.normalized_value.notna()].copy()
    a["normalized_value"]=a.normalized_value.str.strip()
    a=a[a.normalized_value.ne("")&~a.normalized_value.str.contains("|",regex=False)&
        ~a.normalized_value.str.lower().isin(["none","nan","null","unknown","unresolved"])]
    assert a.accepted_species.notna().all() and not a.accepted_species.duplicated().any()
    if scope=="high_quality_only":a=a[a.evidence_quality.eq("high")]
    a["genus"]=a.accepted_species.map(accepted_genus)
    a=a[a.genus.notna()].copy()
    a["y"]=a.normalized_value.eq(target).astype(float)
    a["grade"]=a.evidence_quality
    return a.set_index("accepted_species")[["genus","grade","y"]]

def summarize(d,rosters,tab,controls,perms=999,boots=1999):
    cols=["log_focal","log_other",*controls]
    d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=cols+["group"]).reset_index(drop=True)
    result={"eligible_islands":len(d),"eligible_groups":int(d.group.nunique()),"controls":controls,
            "regional_pool_species":len(tab),"regional_pool_genera":int(tab.genus.nunique())}
    info=[g for g,v in d.groupby("group",sort=True) if len(v)>=2 and v.log_focal.nunique()>1]
    d=d[d.group.isin(info)].reset_index(drop=True)
    result.update(informative_groups=len(info),informative_islands=len(d))
    if len(d)<8 or not len(tab) or tab.y.nunique()<2:return result|{"status":"not_estimable"}
    # The regional botanical pool is frozen independently of animal data eligibility.
    tab=tab.sort_index()
    index={s:i for i,s in enumerate(tab.index)}
    rs=[rosters[i] for i in d.island_id]
    A,n=membership(rs,index);assert (n>=5).all()
    V=tab.y.to_numpy(float)
    codes=tab.groupby(["genus","grade"],sort=True).ngroup().to_numpy()
    E=tab.groupby(["genus","grade"],sort=True).y.transform("mean").to_numpy(float)
    observed=np.asarray(A@V).ravel()
    expected=np.asarray(A@E).ravel()
    Y=np.column_stack([observed,expected,observed-expected])
    x=d[cols].to_numpy(float)
    groups=[v.index.to_numpy() for _,v in d.groupby("group",sort=True)]
    Xs=[];Ys=[]
    for ids in groups:
        z=x[ids];y=Y[ids];w=1/math.sqrt(len(ids))
        Xs.append((z-z.mean(0))*w);Ys.append((y-y.mean(0))*w)
    K=len(groups);p=len(cols)
    result["scored_distinct_species"]=len(set().union(*rs))
    result["median_scored_species"]=float(np.median(n))
    def fit(mult):
        active=np.flatnonzero(mult)
        if not len(active) or sum(len(groups[k])-1 for k in active)-p<4:return None
        X=np.vstack([Xs[k]*math.sqrt(mult[k]) for k in active])
        yy=np.vstack([Ys[k]*math.sqrt(mult[k]) for k in active])
        if np.linalg.matrix_rank(X)<p:return None
        coef=np.linalg.lstsq(X,yy,rcond=None)[0]
        return coef[0]-coef[1]
    point=fit(np.ones(K,int))
    if point is None:return result|{"status":"not_estimable"}
    assert np.isclose(point[0]-point[1],point[2],atol=1e-10)
    # Exact contrast weight on each distinct species, keeping repeated
    # occurrences of its trait synchronized across every island.
    X=np.vstack(Xs);C=np.linalg.pinv(X);w=np.zeros(len(d));pos=0
    for ids in groups:
        a=(C[0,pos:pos+len(ids)]-C[1,pos:pos+len(ids)])/math.sqrt(len(ids))
        w[ids]=a-a.mean()
        pos+=len(ids)
    W=np.asarray(A.T@w).ravel()
    assert np.allclose([W@V,W@E,W@(V-E)],point,rtol=1e-8,atol=1e-10)
    movable=tab.groupby(["genus","grade"],sort=True).y.transform("nunique").gt(1).to_numpy()
    result["movable_regional_pool_species"]=int(movable.sum())
    result["movable_scored_species"]=int((movable & (np.asarray(A.sum(0)).ravel()>0)).sum())
    rng=np.random.default_rng(SEED);base=np.argsort(codes,kind="stable")
    vg=V[base];wg=W[base]
    draws=np.empty(perms)
    # Shuffle globally once per unique plant species, within genus x grade.
    # Sorted block positions preserve the block's number of trait states.
    for k in range(perms):
        order=np.lexsort((rng.random(len(V)),codes))
        draws[k]=wg@V[order]
    null_mean=float(W@E)
    residual=point[2]
    p2=(1+int((np.abs(draws-null_mean)>=abs(residual)-1e-12).sum()))/(perms+1)
    pp=(1+int((draws>=point[0]-1e-12).sum()))/(perms+1)
    bb=[];leave=[]
    for _ in range(boots):
        q=fit(np.bincount(rng.integers(0,K,K),minlength=K))
        if q is not None:bb.append(q)
    for i in range(K):
        m=np.ones(K,int);m[i]=0;q=fit(m)
        if q is not None:leave.append(q)
    ci=np.quantile(np.asarray(bb),[.025,.975],axis=0).T.tolist() if len(bb)>=.9*boots else [None]*3
    lr=np.column_stack([np.min(leave,0),np.max(leave,0)]).tolist() if leave else [None]*3
    gate=K>=30 and len(d)>=200 and len(bb)>=.9*boots and len(leave)==K and result["movable_scored_species"]>=30
    return result|{"status":"exploratory_estimable" if gate else "insufficient_support",
         "contrast_observed":float(point[0]),"contrast_genus_grade_expected":float(point[1]),
         "contrast_beyond_genus_grade":float(residual),
         "group_bootstrap_intervals":dict(zip(["observed","expected","residual"],ci)),
         "leave_group_ranges":dict(zip(["observed","expected","residual"],lr)),
         "valid_bootstrap":len(bb),"valid_leave_group":len(leave),
         "null_analytic_mean":null_mean,"null_simulated_mean":float(draws.mean()),
         "null_quantiles":np.quantile(draws,[.025,.975]).tolist(),
         "p_two_sided_genus_grade":float(p2),"p_positive_genus_grade":float(pp),
         "permutations":perms}

def animal_overlay(fd,geo,tmp):
    sys.path.insert(0,str(fd/"scripts"))
    from nectar_bat_public_flora import diet,ranges,gshhg,incidence,GSHHG_SHA
    from bat_mdd_synonym_audit import build
    d,dm=diet();p,pm=ranges(tmp)
    assert dm["sha256"]=="51acbe6215b3c1f742463ff3fb37ef607150c052e0f7a9dc0d97fcaca61d3f76"
    assert pm["rds_sha256"]=="56daae8dfb79a7821b5a48200f18c28625f4b14d66e9a02efba76b98b36b7e2f"
    units,gh=gshhg(tmp,set(geo.island_id));assert gh==GSHHG_SHA and len(units)==8256
    exact,ec=incidence(units[units.area_km2.ge(20)],p,d)
    aug,audit,family_map,family_audit=build(set(p.name),d,return_taxonomy=True)
    assert audit["reciprocal_singleton_added"]==73 and family_audit["assigned"]==1382
    new,nc=incidence(units[units.area_km2.ge(20)],p,aug,strict_family_map=family_map,exact_diet=d)
    a=exact[.5];b=new[.5]
    assert list(a.island_id)==list(b.island_id) and len(a)==3692
    assert np.array_equal(a.total,b.total) and np.array_equal(a.unknown,b.local_exact_global_unknown)
    assert ec["0.5"]["potential_bat_positive_islands"]==2600 and ec["0.5"]["complete_bat_diet_positive_islands"]==1283
    return a,b,{"range":pm,"diet":dm,"exact_coverage":ec,"mdd_coverage":nc,"family_taxonomy_counts":family_audit}

def panel(flora,t,geo,animals,endpoint,scope,min_n,coverage,family=None):
    trait,target=TRAITS[endpoint]
    s=states(t,trait,target,scope)
    f=flora if scope!="original_native" else flora[flora.origin_status.eq("native")]
    regional=set(f.accepted_species)
    s=s[s.index.isin(regional)].copy()
    full=f.groupby("island_id").accepted_species.nunique()
    # Source pool contains unique tropical Direct species regardless of animal linkage.
    rr=f[f.accepted_species.isin(s.index)].groupby("island_id").accepted_species.agg(set).to_dict()
    d=animals.merge(geo,on="island_id",validate="one_to_one")
    d=d[d.island_id.isin(rr)&d.range_area.ge(20)&d.distance_to_continent_km.ge(0)&d.total.gt(0)].copy()
    d["n"]=d.island_id.map(lambda i:len(rr[i]))
    d=d[d.n.ge(min_n)].copy()
    d["direct_coverage"]=d.n/d.island_id.map(full)
    d["log_area"]=np.log2(d.range_area)
    d["log_distance"]=np.log2(1+d.distance_to_continent_km)
    d["log_n"]=np.log2(d.n)
    d["abs_lat"]=d.island_latitude.abs()
    d=d.dropna(subset=["island_latitude","island_longitude"])
    d["group"]=(np.floor(d.island_latitude/10).astype("Int64").astype(str)+"_"+np.floor(d.island_longitude/10).astype("Int64").astype(str))
    if family is None:
        d=d[d.match_share.ge(coverage)].copy()
        d["log_focal"]=np.log1p(d.nectar);d["log_other"]=np.log1p(d.other)
        d["log_unknown"]=np.log1p(d.unknown)
        controls=BASE_CONTROLS+["log_unknown"]
    else:
        prefix="local_mdd_"+family
        count=d[prefix+"_nectar"]+d[prefix+"_other"]+d[prefix+"_unknown"]
        d=d[count.gt(0)&d.local_mdd_unknown_family.eq(0)&d[prefix+"_unknown"].eq(0)].copy()
        count=count.loc[d.index]
        d["log_focal"]=np.log1p(d[prefix+"_nectar"]);d["log_other"]=np.log1p(d[prefix+"_other"])
        d["log_nonfamily"]=np.log1p(d.total-count)
        controls=BASE_CONTROLS+["log_nonfamily"]
    return summarize(d,rr,s,controls)

def holm(panels):
    # Preserve unsupported panels in full multiplicity family with p=1.
    order=sorted(panels,key=lambda k:panels[k].get("p_two_sided_genus_grade",1.))
    mx=0.
    for i,k in enumerate(order):
        a=panels[k];p=a.get("p_two_sided_genus_grade",1.)
        mx=max(mx,min(1.,(len(order)-i)*p));a["holm_p_two_sided"]=float(mx)
        ci=a.get("group_bootstrap_intervals",{}).get("residual")
        lr=a.get("leave_group_ranges",{}).get("residual")
        a["supported_positive_residual"]=bool(a["status"]=="exploratory_estimable" and mx<.05 and ci and ci[0]>0 and lr and lr[0]>0)

def self_test():
    rng=np.random.default_rng(SEED)
    # 40 groups; a distinct plant trait per accepted name, never per occurrence.
    N=240;P=120
    ids=[f"Species {i}" for i in range(P)]
    tab=pd.DataFrame({"genus":np.repeat(np.arange(30),4),"grade":"high","y":np.tile([0.,0.,1.,1.],30)},index=ids)
    d=pd.DataFrame({"group":np.repeat(np.arange(40),6),"island_id":np.arange(N),
                    "log_focal":rng.normal(size=N),"log_other":rng.normal(size=N),"c":rng.normal(size=N)})
    rr={i:set(rng.choice(ids,20,replace=False)) for i in range(N)}
    r=summarize(d,rr,tab,["c"],perms=99,boots=99)
    assert r["status"]=="exploratory_estimable" and r["informative_islands"]==N
    assert np.isclose(r["contrast_genus_grade_expected"],0.,atol=1e-10)
    assert np.isclose(r["contrast_observed"],r["contrast_beyond_genus_grade"])
    d["log_focal"]=0.
    assert summarize(d,rr,tab,["c"],perms=99,boots=99)["status"]=="not_estimable"
    print("ORIGINAL_BAT_GENUS_NULL_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--self-test",action="store_true")
    p.add_argument("--original-root",type=pathlib.Path);p.add_argument("--geography",type=pathlib.Path)
    p.add_argument("--fd-root",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path)
    a=p.parse_args()
    if a.self_test:self_test();return
    root=a.original_root
    flora=pd.read_csv(io.BytesIO(read_pin(root/"fixed/canonical/input/chapter1_status_flora.csv.gz","flora")),compression="gzip",dtype=str,usecols=["island_id","accepted_species","origin_status"])
    t=pd.read_csv(io.BytesIO(read_pin(root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz","traits")),compression="gzip",dtype=str,usecols=["accepted_species","trait_name","normalized_value","evidence_scope","evidence_quality"])
    geo=pd.read_csv(io.BytesIO(read_pin(a.geography,"geography")),dtype={"island_id":str})
    assert flora.island_id.nunique()==4505 and not flora.duplicated(["island_id","accepted_species"]).any()
    assert geo.island_id.nunique()==8264 and geo.island_id.is_unique
    trop=geo[geo.analysis_regime.eq("tropical")].copy()
    f=flora[flora.island_id.isin(trop.island_id)].copy()
    assert f.island_id.nunique()==1558
    with tempfile.TemporaryDirectory() as td:
        exact,mdd,source=animal_overlay(a.fd_root,geo,pathlib.Path(td))
    print("ORIGINAL_BAT_OVERLAY_VERIFIED",json.dumps(source,allow_nan=False),flush=True)
    panels={}
    for endpoint in TRAITS:
        for scope,n,cov in SCOPES:
            k=f"{endpoint}_{scope}_n{n}_cov{cov}"
            panels[k]=panel(f,t,trop,exact,endpoint,scope,n,cov)
        for family in FAMILIES:
            k=f"{endpoint}_family_{family}_mdd"
            panels[k]=panel(f,t,trop,mdd,endpoint,"all_recorded",5,.5,family)
    assert len(panels)==21
    holm(panels)
    out={"original_botanical_parent":4505,"original_geography_parent":8264,"original_tropical_parent":1558,
         "source":source,"panels":panels,"primary_post_pilot":"white_all_recorded_n5_cov0.5",
         "warning":"Post-pilot exploratory tests; original species-wide Direct traits, not island-specific phenotypes. Potential range overlap and recorded nectar/pollen/exudate diet do not establish realized pollination. Plant genus-grade expectation from entire eligible tropical plant pool, not ancestral source flora; fixed-incidence species label permutation is not phylogenetic independent contrast. Geographic bootstrap is not taxon-independent. Only anonymous aggregate outputs."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    for k,v in panels.items():print("ORIGINAL_BAT_GENUS_PANEL",k,json.dumps(v,allow_nan=False),flush=True)
if __name__=="__main__":main()
