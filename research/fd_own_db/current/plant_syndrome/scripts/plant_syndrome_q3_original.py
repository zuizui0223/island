#!/usr/bin/env python3
"""Frozen Q3: within-genus/grade morphology-proxy H x reproductive assurance.
Same accepted-name complete Direct records; no missing trait imputation.
Contemporary source candidates are not ancestors. Anonymous outputs only.
"""
import argparse,functools,json,math,pathlib,sys
import numpy as np
import pandas as pd
from plant_syndrome_core import CORE,encode
from plant_syndrome_q2_original import read_archive,trait_table,bootstrap,MODES,ARCHIVE_SHA,PIN,GEO_SHA,SEED

PRIMARY="autonomous_all_recorded_geo_k5_inclusive"
OUTCOMES={"autonomous_selfing_capacity":("autonomous","absent","autonomous"),
          "self_incompatibility":("SC","SI","SC")}
VARIANTS=("all_recorded","original_native","high_quality_only")

def state_table(t,axis,variant):
    pos,neg,_=OUTCOMES[axis]
    tab=trait_table(t,variant)
    a=t[t.trait_name.eq(axis)&t.evidence_scope.eq("species_direct")&
        t.evidence_quality.isin(["high","medium"])&t.normalized_value.isin([pos,neg])].copy()
    assert a.accepted_species.notna().all() and not a.accepted_species.duplicated().any()
    if variant=="high_quality_only":a=a[a.evidence_quality.eq("high")]
    a=a.set_index("accepted_species")
    tab=tab.join(a[["normalized_value","evidence_quality"]],how="inner")
    tab["assurance"]=tab.normalized_value.eq(pos).astype(float)
    tab["grade_pattern"]=tab.grade_pattern+"|"+tab.evidence_quality
    return tab

def covar(x,y,group):
    n=len(x);assert n>1
    scale=n/(n-1);raw=np.mean(x*y)-x.mean()*y.mean()
    _,g=np.unique(group,return_inverse=True)
    k=np.bincount(g);mx=np.bincount(g,weights=x)/k;my=np.bincount(g,weights=y)/k
    b=np.sum(k/n*mx*my)-x.mean()*y.mean()
    result=scale*np.array([raw,raw-b,b])
    assert np.isclose(result[0],result[1]+result[2],atol=1e-12)
    informative=int(sum(len(np.unique(y[g==j]))>1 and np.ptp(x[g==j])>1e-12 for j in range(len(k))))
    return result,informative

def stats(tab):
    _,_,_,_,s,_=encode(tab,CORE)
    H=s[:,:4].mean(1)-s[:,4];y=tab.assurance.to_numpy(float)
    g=pd.factorize(pd.MultiIndex.from_frame(tab[["genus","grade_pattern"]]),sort=True)[0]
    gen=pd.factorize(tab.genus,sort=True)[0]
    def run(roster):
        ix=np.array(sorted(roster),int)
        c,k=covar(H[ix],y[ix],g[ix]);c2,k2=covar(H[ix],y[ix],gen[ix])
        return np.array([c[1],c[0],c[2],c2[1],c2[2]]),k,k2,float(y[ix].mean())
    return run

def summary(rows):
    if not rows:return {"status":"not_estimable","eligible_islands":0}
    d=np.stack([v["delta"] for v in rows])
    g=bootstrap(d,[v["block"] for v in rows]);s=bootstrap(d,[v["first_source"] for v in rows])
    support=all(v["groups"]>=30 and v["islands"]>=200 and v["valid_draws"]>=1800 and v["nondegenerate"] for v in (g,s))
    # Positive-direction tail. Recompute by applying the existing negative
    # tail function to -delta, not 1-p (zero draws count in both tails).
    gp=bootstrap(-d,[v["block"] for v in rows]);sp=bootstrap(-d,[v["first_source"] for v in rows])
    p=max(gp["tail_negative"][0],sp["tail_negative"][0]) if support else None
    return {"status":"exploratory_estimable" if support else "insufficient_support",
        "eligible_islands":len(rows),"metric_order":["within_genus_grade","total","between_genus_grade","within_genus","between_genus"],
        "spatial_block":g,"first_source_group":s,
        "worst_directional_bootstrap_tail":p,
        "interim_four_question_p_upper_bound":min(1.,4*p) if p is not None else None,
        "interim_directional_candidate":bool(support and 4*p<.05 and g["interval"][0][0]>0 and s["interval"][0][0]>0),
        "target_informative_genus_grade_strata_median":float(np.median([v["kt"] for v in rows])),
        "source_informative_genus_grade_strata_median":float(np.median([v["ks"] for v in rows])),
        "island_equal_assurance_share_delta":float(np.mean([v["assurance_delta"] for v in rows]))}

def analyze(f,t,a,m,geo):
    full=f.groupby("island_id").accepted_species.agg(set).to_dict()
    region=m.groupby("entity_ID").work_species.agg(set).to_dict()
    keys=a.sort_values("source_rank").groupby(["island_id","source_mode"]).entity_ID.agg(tuple).to_dict()
    geo=geo.set_index("island_id")
    blocks={}
    for i in sorted(full.keys()&set(geo.index)):
        z=geo.loc[i];lat=pd.to_numeric(z.island_latitude,errors="coerce");lon=pd.to_numeric(z.island_longitude,errors="coerce")
        if np.isfinite([lat,lon]).all():blocks[i]=f"{math.floor(lat/10)}_{math.floor(lon/10)}"
    panels={}
    for axis,(_,__,name) in OUTCOMES.items():
        for variant in VARIANTS:
            frame=state_table(t,axis,variant);idx={s:i for i,s in enumerate(frame.index)}
            calc=stats(frame)
            f0=f[f.origin_status.eq("native")] if variant=="original_native" else f
            tar={i:{idx[s] for s in v if s in idx} for i,v in f0.groupby("island_id").accepted_species.agg(set).items()}
            all_idx={i:{idx[s] for s in v if s in idx} for i,v in full.items()}
            @functools.lru_cache(maxsize=512)
            def pool(ids):return frozenset(idx[s] for e in ids for s in region[e] if s in idx)
            @functools.lru_cache(maxsize=512)
            def cached(ids):return calc(pool(ids))
            for mode,nr in MODES.items():
                rows={"inclusive":[],"disjoint":[]}
                census={"complete_trait_species":len(frame),"target_geographic_eligible":0,"paired_size_eligible":0,
                        "target_within_stratum_informative":0,"both_source_variants_informative":0,"missing_donor":0}
                for i,block in blocks.items():
                    tt=tar.get(i,set())
                    if len(tt)<10:continue
                    census["target_geographic_eligible"]+=1
                    ids=keys.get((i,mode),())
                    if len(ids)!=nr or any(e not in region for e in ids):
                        census["missing_donor"]+=1;continue
                    inc=pool(ids);exc=inc-all_idx[i]
                    if len(exc)<10:continue
                    census["paired_size_eligible"]+=1
                    u,kt,_,py=calc(tt)
                    vi,ki,_,pi=cached(ids);ve,ke,_,pe=calc(exc)
                    # Preserve zero within-genus covariance for monomorphic
                    # or singleton strata in primary. Report identifiability
                    # counts explicitly; never treat these as independent
                    # evidence of no biological association.
                    census["target_within_stratum_informative"]+=int(kt>0)
                    census["both_source_variants_informative"]+=int(kt>0 and ki>0 and ke>0)
                    for tag,v,k,pp in (("inclusive",vi,ki,pi),("disjoint",ve,ke,pe)):
                        rows[tag].append({"delta":u-v,"block":block,"first_source":ids[0],
                            "kt":kt,"ks":k,"assurance_delta":py-pp})
                for tag,rr in rows.items():
                    key=f"{name}_{variant}_{mode}_{tag}"
                    out=summary(rr);out["support_census"]=census
                    # In addition to frozen geographic support, prevent
                    # a degenerate zero result from being called supported.
                    if census["both_source_variants_informative"]<200:
                        out["within_stratum_information_gate"]=False
                        out["interim_directional_candidate"]=False
                        if out["status"]=="exploratory_estimable":out["status"]="insufficient_within_stratum_information"
                    else:out["within_stratum_information_gate"]=True
                    panels[key]=out
                    print("PLANT_Q3_PANEL",key,json.dumps(out,allow_nan=False),flush=True)
            pool.cache_clear();cached.cache_clear()
    assert len(panels)==24
    return panels

def self_test():
    x=np.array([0,1,0,1.]);y=np.array([0,1,0,1.]);g=np.array([0,0,1,1])
    a,k=covar(x,y,g);assert np.allclose(a,[1/3,1/3,0]) and k==2
    a,k=covar(x,y,np.arange(4));assert np.allclose(a,[1/3,0,1/3]) and k==0
    print("PLANT_Q3_COVARIANCE_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--archive",type=pathlib.Path);p.add_argument("--geography",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    try:
        f,t,s,m,g,md=read_archive(a.archive,a.geography)
        panels=analyze(f,t,s,m,g)
        out={"original_botanical_parent":4505,"original_geography_parent":8264,"mainland_candidate_regions":748,
             "archive_sha256":ARCHIVE_SHA,"mainland_sha256":md,"source_file_pins":PIN|{"geography":GEO_SHA},
             "primary_q3":PRIMARY,"panels":panels,"warning":"Genuine same-species Direct morphology and reproductive category. Autonomous versus absent separate from SC versus SI. Zero within-stratum covariance can reflect lack of variation rather than biological absence. Additional conservative within-stratum information gate added before Q3 output. No individual records or causal evolutionary claims."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    except Exception as e:
        print("PLANT_Q3_ORIGINAL_FAILED",type(e).__name__,flush=True);sys.exit(1)
if __name__=="__main__":main()
