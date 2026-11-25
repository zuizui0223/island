#!/usr/bin/env python3
"""Q4: source-roster inclusion associated with fixed floral morphology.
Native listed vs entirely unlisted in inclusive mainland candidate pools.
Nonlisted is not true absence, colonization failure or adaptation.
Only anonymous aggregates leave the original-source analysis.
"""
import argparse,functools,json,math,pathlib,sys
import numpy as np
import pandas as pd
from plant_syndrome_core import CORE,TEMPLATES,encode
from plant_syndrome_q2_original import read_archive,trait_table,bootstrap,MODES,ARCHIVE_SHA,PIN,GEO_SHA,SEED

PRIMARY="native_geo_k5_genus_grade"
VARIANTS=("native","all_recorded","native_high_quality","native_strict_red")
ENDS=("H",*TEMPLATES)

def matched_difference(present,absent,group,genus,V):
    # Equal weights across matched genus-grade strata; no comparisons
    # between different genera or documentation patterns in the primary.
    # A minimum of two distinct genera and five distinct names per class
    # AFTER restriction to common strata is fixed before observing Q4.
    p=np.array(sorted(present),int);a=np.array(sorted(absent),int)
    if len(p)<5 or len(a)<5:return None
    assert not present&absent
    pg={};ag={}
    for i in p:pg.setdefault(int(group[i]),[]).append(i)
    for i in a:ag.setdefault(int(group[i]),[]).append(i)
    common=sorted(pg.keys()&ag.keys())
    if len(common)<2:return None
    lp=[i for k in common for i in pg[k]]
    la=[i for k in common for i in ag[k]]
    ng=len(set(genus[lp]))
    if min(len(lp),len(la))<5 or ng<2:return None
    delta=np.stack([V[pg[k]].mean(0)-V[ag[k]].mean(0) for k in common])
    balance=np.array([min(len(pg[k]),len(ag[k])) for k in common])
    return {"delta":delta.mean(0),"overlap_weighted_delta":np.average(delta,axis=0,weights=balance),
            "matched_present":len(lp),"matched_unlisted":len(la),"matched_genera":ng,
            "matched_strata":len(common),"present_idx":set(lp),"absent_idx":set(la)}

def summarize(rows):
    out={"eligible_islands":len(rows),"endpoint_order":list(ENDS)}
    if not rows:return out|{"status":"not_estimable"}
    x=np.stack([v["delta"] for v in rows])
    g=bootstrap(x,[v["block"] for v in rows])
    s=bootstrap(x,[v["first_source"] for v in rows])
    support=all(v["groups"]>=30 and v["islands"]>=200 and v["valid_draws"]>=1800 and v["nondegenerate"] for v in (g,s))
    p=max(g["tail_negative"][0],s["tail_negative"][0]) if support else None
    w=np.stack([v["overlap_weighted_delta"] for v in rows])
    present=set().union(*(v["present_idx"] for v in rows))
    absent=set().union(*(v["absent_idx"] for v in rows))
    out|={"status":"exploratory_estimable" if support else "insufficient_support",
          "spatial_block":g,"first_source_group":s,
          "worst_directional_bootstrap_tail":p,"interim_four_question_p_upper_bound":min(1.,4*p) if p is not None else None,
          "interim_directional_candidate":bool(support and p*4<.05 and g["interval"][0][1]<0 and s["interval"][0][1]<0),
          "distinct_matched_listed_taxa":len(present),"distinct_matched_unlisted_taxa":len(absent),
          "distinct_matched_union_taxa":len(present|absent),
          "median_matched_listed":float(np.median([v["matched_present"] for v in rows])),
          "median_matched_unlisted":float(np.median([v["matched_unlisted"] for v in rows])),
          "median_matched_genera":float(np.median([v["matched_genera"] for v in rows])),
          "overlap_weighted_sensitivity":{"spatial_block":bootstrap(w,[v["block"] for v in rows]),
                                           "first_source_group":bootstrap(w,[v["first_source"] for v in rows])},
          "regional_descriptive":{}}
    for reg in sorted(set(v["regime"] for v in rows)):
        rr=[v for v in rows if v["regime"]==reg]
        out["regional_descriptive"][reg]=bootstrap(np.stack([v["delta"] for v in rr]),[v["block"] for v in rr])
    return out

def analyze(f,t,a,m,g):
    full=f.groupby("island_id").accepted_species.agg(set).to_dict()
    native=f[f.origin_status.eq("native")].groupby("island_id").accepted_species.agg(set).to_dict()
    region=m.groupby("entity_ID").work_species.agg(set).to_dict()
    keys=a.sort_values("source_rank").groupby(["island_id","source_mode"]).entity_ID.agg(tuple).to_dict()
    geo=g.set_index("island_id")
    geoc={}
    for i in sorted(full.keys()&set(geo.index)):
        z=geo.loc[i];lat=pd.to_numeric(z.island_latitude,errors="coerce");lon=pd.to_numeric(z.island_longitude,errors="coerce")
        if not np.isfinite([lat,lon]).all() or pd.isna(z.analysis_regime):continue
        geoc[i]=(f"{math.floor(lat/10)}_{math.floor(lon/10)}",str(z.analysis_regime))
    out={}
    for v in VARIANTS:
        frame=trait_table(t,"high_quality_only" if v=="native_high_quality" else "all_recorded")
        idx={s:i for i,s in enumerate(frame.index)}
        _,_,_,_,sc,_=encode(frame,CORE,red_only=v=="native_strict_red")
        H=sc[:,:4].mean(1)-sc[:,4];V=np.column_stack([H,sc])
        gen=pd.factorize(frame.genus,sort=True)[0]
        gg=pd.factorize(pd.MultiIndex.from_frame(frame[["genus","grade_pattern"]]),sort=True)[0]
        all_idx={i:{idx[s] for s in rr if s in idx} for i,rr in full.items()}
        tar_idx=all_idx if v=="all_recorded" else {i:{idx[s] for s in rr if s in idx} for i,rr in native.items()}
        @functools.lru_cache(maxsize=512)
        def pool(ids):return frozenset(idx[s] for e in ids for s in region[e] if s in idx)
        for mode,nr in MODES.items():
            rows={"genus_grade":[],"genus":[]};census={"complete_trait_species":len(frame),"target_names_ge5":0,"source_candidates_available":0,"both_classes_ge5":0}
            for i,(block,reg) in geoc.items():
                tar=tar_idx.get(i,set())
                if len(tar)<5:continue
                census["target_names_ge5"]+=1
                ids=keys.get((i,mode),())
                if not ids:continue
                assert len(ids)==nr
                if any(e not in region for e in ids):continue
                census["source_candidates_available"]+=1
                source=pool(ids)
                present=set(source&tar)
                absent=set(source-all_idx[i])
                # Native primary excludes ALL listed nonnative/unresolved
                # species from both classes, not a false unlisted control.
                assert not present&absent
                if len(present)<5 or len(absent)<5:continue
                census["both_classes_ge5"]+=1
                for label,groups in (("genus_grade",gg),("genus",gen)):
                    res=matched_difference(present,absent,groups,gen,V)
                    if res is None:continue
                    rows[label].append(res|{"block":block,"regime":reg,"first_source":ids[0]})
            for label,rr in rows.items():
                key=f"{v}_{mode}_{label}"
                p=summarize(rr);p["support_census"]=census
                out[key]=p
                print("PLANT_Q4_PANEL",key,json.dumps(p,allow_nan=False),flush=True)
        pool.cache_clear()
    assert len(out)==16
    return out

def self_test():
    g=np.repeat([0,1],6);gen=g.copy()
    vals=np.arange(12,dtype=float)[:,None]
    # In each stratum the listed class has H lower by 3.
    p={0,1,2,6,7,8};a={3,4,5,9,10,11}
    q=matched_difference(p,a,g,gen,vals)
    assert q["matched_genera"]==2 and q["matched_present"]==6 and np.allclose(q["delta"],[-3])
    assert matched_difference(p,a,g,np.zeros(12,int),vals) is None
    assert matched_difference({0,1},a,g,gen,vals) is None
    print("PLANT_Q4_MATCHED_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--archive",type=pathlib.Path);p.add_argument("--geography",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    if not all((a.archive,a.geography,a.out)):p.error("archive, geography, out required")
    try:
        f,t,s,m,g,md=read_archive(a.archive,a.geography)
        panels=analyze(f,t,s,m,g)
        out={"original_botanical_parent":4505,"original_geography_parent":8264,"mainland_candidate_regions":748,
             "archive_sha256":ARCHIVE_SHA,"mainland_sha256":md,"source_file_pins":PIN|{"geography":GEO_SHA},
             "primary_q4":PRIMARY,"panels":panels,"warning":"Post-pilot exploratory Q4. Contemporary native candidate pools, listed vs entirely unlisted species. Nonlisted not absence or failed establishment. Equal overlapping genus-grade strata, then equal islands within operational groups. Shared species/regions across islands remain dependent. Not within-population evolution."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    except Exception as e:
        print("PLANT_Q4_ORIGINAL_FAILED",type(e).__name__,flush=True);sys.exit(1)
if __name__=="__main__":main()
