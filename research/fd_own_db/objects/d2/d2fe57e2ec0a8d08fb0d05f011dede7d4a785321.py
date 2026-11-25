#!/usr/bin/env python3
"""Post-pilot pairwise-complete floral syndrome covariance atlas.

Keep every two-axis comparison on one identical carrier roster for all
applicable templates. Never impute missing states or compare template
differences across different pair-specific carrier sets.
"""
import argparse,functools,json,math,pathlib,sys
from itertools import combinations
import numpy as np
import pandas as pd
from plant_syndrome_core import CORE,SIZE,VALID,TEMPLATES,components
from plant_syndrome_q2_original import read_archive,bootstrap,genus,MODES,ARCHIVE_SHA,PIN,GEO_SHA,SEED

AXES=CORE+(SIZE,)
PAIRS=tuple(combinations(AXES,2))
SCOPES=("all_recorded","original_native")
METRICS=("total","within_genus_grade","between_genus_grade","within_genus","between_genus")
FOCUS={tuple(sorted((SIZE,CORE[0]))),tuple(sorted((SIZE,CORE[1])))}

def pair_table(t,pair):
    a=t[t.trait_name.isin(pair)&t.evidence_scope.eq("species_direct")&
        t.evidence_quality.isin(["high","medium"])&t.normalized_value.notna()].copy()
    a["normalized_value"]=a.normalized_value.str.strip()
    assert a.accepted_species.notna().all() and not a.duplicated(["accepted_species","trait_name"]).any()
    keep=np.zeros(len(a),bool)
    for k in pair:keep|=(a.trait_name.eq(k)&a.normalized_value.isin(VALID[k])).to_numpy()
    a=a[keep]
    v=a.pivot(index="accepted_species",columns="trait_name",values="normalized_value").reindex(columns=pair).dropna()
    q=a.pivot(index="accepted_species",columns="trait_name",values="evidence_quality").reindex(index=v.index,columns=pair)
    assert q.notna().all().all()
    v["genus"]=[genus(s) for s in v.index]
    v["grade_pattern"]=q.astype(str).agg("|".join,axis=1)
    return v[v.genus.notna()].sort_index()

def pair_stat(frame,pair):
    templates=[k for k,v in TEMPLATES.items() if all(a in v for a in pair)]
    assert templates
    X=np.column_stack([frame[a].isin(set(TEMPLATES[k][a].split())).to_numpy(float) for k in templates for a in pair])
    pp=np.arange(X.shape[1]).reshape(-1,2)
    gg=pd.factorize(pd.MultiIndex.from_frame(frame[["genus","grade_pattern"]]),sort=True)[0]
    g=pd.factorize(frame.genus,sort=True)[0]
    def calc(ids):
        ix=np.array(sorted(ids),int);assert len(ix)>=10
        a=components(X[ix],pp,gg[ix]);b=components(X[ix],pp,g[ix])
        v=np.stack([a["total"],a["within"],a["between"],b["within"],b["between"]],axis=1)
        assert np.allclose(v[:,0],v[:,1]+v[:,2],atol=1e-12)
        return v,a["phi"],int(len(np.unique(g[ix])))
    return templates,calc

def summary(rows,templates,focus):
    out={"eligible_islands":len(rows),"templates":templates,"metric_order":list(METRICS)}
    if not rows:return out|{"status":"not_estimable"}
    d=np.stack([v["delta"] for v in rows])
    assert np.allclose(d[:,:,0],d[:,:,1]+d[:,:,2],atol=1e-12)
    assert np.allclose(d[:,:,0],d[:,:,3]+d[:,:,4],atol=1e-12)
    x=d.reshape(len(d),-1)
    bg=[r["block"] for r in rows];sg=[r["first_source"] for r in rows]
    s=bootstrap(x,bg);f=bootstrap(x,sg)
    gate=all(v["groups"]>=30 and v["islands"]>=200 and v["valid_draws"]>=1800 for v in (s,f))
    out|={"status":"exploratory_estimable" if gate else "insufficient_support",
        "spatial_block":s,"first_source_group":f,
        "median_target_scored":float(np.median([r["target_n"] for r in rows])),
        "median_source_scored":float(np.median([r["source_n"] for r in rows])),
        "median_target_genera":float(np.median([r["target_genera"] for r in rows])),
        "distinct_target_taxa":len(set().union(*(r["target_idx"] for r in rows))),
        "distinct_source_taxa":len(set().union(*(r["source_idx"] for r in rows))),
        "phi_by_template":{}}
    for j,k in enumerate(templates):
        rr=[r for r in rows if np.isfinite(r["phi"][j])]
        out["phi_by_template"][k]={"islands":len(rr),
          "spatial_block":bootstrap(np.array([r["phi"][j] for r in rr]),[r["block"] for r in rr]),
          "first_source_group":bootstrap(np.array([r["phi"][j] for r in rr]),[r["first_source"] for r in rr])}
    if focus:
        assert "bird" in templates and len(templates)>1
        # Both focal contrasts were defined after seeing the five-axis
        # pilot; this is not an independent replication or primary Q2.
        j=templates.index("bird")
        contrast=d[:,j,:]-np.delete(d,j,axis=1).mean(1)
        a=bootstrap(contrast,bg);b=bootstrap(contrast,sg)
        phi_rows=[r for r in rows if np.isfinite(r["phi"]).all()]
        pc=np.array([r["phi"][j]-np.delete(r["phi"],j).mean() for r in phi_rows])
        out["bird_minus_other_templates"]={
            "spatial_block":a,"first_source_group":b,
            "both_intervals_positive":bool(gate and a["interval"][0][0]>0 and b["interval"][0][0]>0),
            "common_nonconstant_phi":{
                "spatial_block":bootstrap(pc,[r["block"] for r in phi_rows]),
                "first_source_group":bootstrap(pc,[r["first_source"] for r in phi_rows])}}
    # Regional diagnostics are descriptive and not independent primary tests.
    out["regional_descriptive"]={}
    for reg in sorted(set(r["regime"] for r in rows)):
        rr=[r for r in rows if r["regime"]==reg]
        z=np.stack([r["delta"][:,0] for r in rr])
        out["regional_descriptive"][reg]={"islands":len(rr),"spatial_block":bootstrap(z,[r["block"] for r in rr])}
    return out

def analyze(f,t,a,m,g):
    full=f.groupby("island_id").accepted_species.agg(set).to_dict()
    nat=f[f.origin_status.eq("native")].groupby("island_id").accepted_species.agg(set).to_dict()
    region=m.groupby("entity_ID").work_species.agg(set).to_dict()
    keys=a.sort_values("source_rank").groupby(["island_id","source_mode"]).entity_ID.agg(tuple).to_dict()
    geo=g.set_index("island_id")
    geoc={}
    for i in sorted(full.keys()&set(geo.index)):
        z=geo.loc[i];lat=pd.to_numeric(z.island_latitude,errors="coerce");lon=pd.to_numeric(z.island_longitude,errors="coerce")
        if not np.isfinite([lat,lon]).all() or pd.isna(z.analysis_regime):continue
        geoc[i]=(f"{math.floor(lat/10)}_{math.floor(lon/10)}",str(z.analysis_regime))
    results={}
    for pair in PAIRS:
        frame=pair_table(t,pair);lookup={s:i for i,s in enumerate(frame.index)}
        templates,calc=pair_stat(frame,pair)
        focus=tuple(sorted(pair)) in FOCUS
        full_idx={i:{lookup[s] for s in v if s in lookup} for i,v in full.items()}
        native_idx={i:{lookup[s] for s in v if s in lookup} for i,v in nat.items()}
        @functools.lru_cache(maxsize=512)
        def pool(ids):return frozenset(lookup[s] for e in ids for s in region[e] if s in lookup)
        @functools.lru_cache(maxsize=512)
        def pooled(ids):return calc(pool(ids))
        for scope in SCOPES:
            target=full_idx if scope=="all_recorded" else native_idx
            for mode,nr in MODES.items():
                rows={"inclusive":[],"disjoint":[]};census={"complete_pair_species":len(frame),"target_geographic_eligible":0,
                    "missing_donor":0,"paired_donor_eligible":0}
                for i,(block,reg) in geoc.items():
                    tr=target.get(i,set())
                    if len(tr)<10:continue
                    census["target_geographic_eligible"]+=1
                    ids=keys.get((i,mode),())
                    if not ids: census["missing_donor"]+=1;continue
                    assert len(ids)==nr
                    if any(e not in region for e in ids):census["missing_donor"]+=1;continue
                    inc=pool(ids);exc=inc-full_idx[i]
                    if len(exc)<10:continue
                    census["paired_donor_eligible"]+=1
                    assert not tr&exc
                    u,uphi,ng=calc(tr)
                    for tag,src,res in (("inclusive",inc,pooled(ids)),("disjoint",exc,calc(exc))):
                        v,vphi,_=res
                        rows[tag].append({"delta":u-v,"phi":uphi-vphi,"target_n":len(tr),"source_n":len(src),
                            "target_genera":ng,"target_idx":tr,"source_idx":src,
                            "block":block,"regime":reg,"first_source":ids[0]})
                assert len(rows["inclusive"])==len(rows["disjoint"])
                for tag,rr in rows.items():
                    key="__".join((*pair,scope,mode,tag))
                    p=summary(rr,templates,focus);p["support_census"]=census
                    results[key]=p
                    print("PLANT_PAIR_ATLAS_PANEL",key,json.dumps(p,allow_nan=False),flush=True)
        pool.cache_clear();pooled.cache_clear()
    assert len(results)==80
    return results

def self_test():
    # Same recorded species across all templates, exactly paired axes.
    n=24;pair=(CORE[1],SIZE)
    frame=pd.DataFrame({pair[0]:["tubular","open_radial","bell_campanulate","bilabiate"]*6,
                        pair[1]:["large","small","large","medium"]*6,
                        "genus":np.repeat(["A","B","C"],8),"grade_pattern":"high|high"})
    names,stat=pair_stat(frame,pair)
    assert names==["specialist_bee","bird","bat","open_generalist"]
    d,phi,g=stat(set(range(n)))
    assert d.shape==(4,5) and np.allclose(d[:,0],d[:,1]+d[:,2])
    assert np.isclose(d[1,0],n/(n-1)*(.25-.25*.5))
    print("PLANT_PAIR_ATLAS_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--archive",type=pathlib.Path);p.add_argument("--geography",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    try:
        f,t,s,m,g,md=read_archive(a.archive,a.geography)
        panels=analyze(f,t,s,m,g)
        out={"original_botanical_parent":4505,"original_geography_parent":8264,"mainland_candidate_regions":748,
             "archive_sha256":ARCHIVE_SHA,"mainland_sha256":md,"source_file_pins":PIN|{"geography":GEO_SHA},
             "panels":panels,"warning":"Post-outcome all-pair atlas, not independent replication of complete-core Q2. Every comparison uses genuine same-species two-axis Direct data and a common carrier roster for applicable templates. Different pairs use different carriers and cannot be pooled as independent. Bird-minus-others contrasts are exploratory, and groups share taxa/regions. No actual guild or selective evolutionary claim."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    except Exception as e:
        print("PLANT_PAIR_ATLAS_FAILED",type(e).__name__,flush=True);sys.exit(1)
if __name__=="__main__":main()
