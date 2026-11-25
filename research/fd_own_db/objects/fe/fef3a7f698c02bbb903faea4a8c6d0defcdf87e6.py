#!/usr/bin/env python3
"""FD-owned Q2, complete same-species GIFT traits, anonymous summaries only."""
import argparse,functools,hashlib,io,json,math,pathlib,re,sys,zipfile
import numpy as np
import pandas as pd
from plant_syndrome_core import CORE,SIZE,VALID,TEMPLATES,encode,components,common_phi,self_test as core_test

SEED=20261125
ARCHIVE_SHA="b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999"
PIN={
 "fixed/canonical/input/chapter1_status_flora.csv.gz":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
 "fixed/source/current/evidence/direct_only_trait_ledger.csv.gz":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
 "fixed/source/result/island_source_assignments.csv.gz":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44"}
GEO_SHA="08218617b5880153986012e79817de65a74dd4e44b5a65e8a537dab8cb2cd0a7"
MAINLAND="fixed/source/source/gift_native_mainland_flora.csv.gz"
MODES={"geo_k5":5,"geo50_climate10":10}
VARIANTS=("all_recorded","original_native","high_quality_only","strict_red","size_core","size_expanded")
PRIMARY="all_recorded_geo_k5_inclusive"
METRICS=("total","within","between","within_genus","between_genus")

def sha(b):return hashlib.sha256(b).hexdigest()
def file_sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()

def read_archive(path,geo_path):
    assert file_sha(path)==ARCHIVE_SHA
    def load(z,suffix,cols):
        names=[n for n in z.namelist() if n.endswith("/"+suffix) or n==suffix]
        assert len(names)==1
        raw=z.read(names[0]);digest=sha(raw)
        if suffix in PIN:assert digest==PIN[suffix]
        return pd.read_csv(io.BytesIO(raw),compression="gzip",dtype=str,usecols=cols),digest
    with zipfile.ZipFile(path) as z:
        f,_=load(z,next(iter(PIN)),["island_id","accepted_species","origin_status"])
        t,_=load(z,list(PIN)[1],["accepted_species","trait_name","normalized_value","evidence_scope","evidence_quality"])
        a,_=load(z,list(PIN)[2],["island_id","source_mode","source_rank","entity_ID"])
        m,md=load(z,MAINLAND,["entity_ID","work_species"])
    raw=geo_path.read_bytes();assert sha(raw)==GEO_SHA
    g=pd.read_csv(io.BytesIO(raw),dtype={"island_id":str})
    assert f.island_id.nunique()==4505 and not f.duplicated(["island_id","accepted_species"]).any()
    assert g.island_id.nunique()==8264 and g.island_id.is_unique
    assert m.entity_ID.nunique()==748 and not m.duplicated(["entity_ID","work_species"]).any()
    assert not f[["island_id","accepted_species"]].isna().any().any()
    assert not m[["entity_ID","work_species"]].isna().any().any()
    a=a[a.source_mode.isin(MODES)].copy()
    a.source_rank=pd.to_numeric(a.source_rank,errors="raise").astype(int)
    assert not a.duplicated(["island_id","source_mode","source_rank"]).any()
    assert a.source_mode.value_counts().to_dict()=={"geo50_climate10":82190,"geo_k5":41325}
    return f,t,a,m,g,md

def genus(s):
    a=str(s).strip().split()
    return a[0] if len(a)>=2 and a[0][0].isupper() and a[0].replace("-","").isalpha() and a[1].lower() not in {"sp.","spp.","cf.","aff.","?"} else None

def trait_table(t,variant):
    axes=CORE+(SIZE,) if variant in ("size_core","size_expanded") else CORE
    a=t[t.trait_name.isin(axes)&t.evidence_scope.eq("species_direct")&t.evidence_quality.isin(["high","medium"])&t.normalized_value.notna()].copy()
    a["normalized_value"]=a.normalized_value.str.strip()
    assert a.accepted_species.notna().all() and not a.duplicated(["accepted_species","trait_name"]).any()
    if variant=="high_quality_only":a=a[a.evidence_quality.eq("high")]
    valid=np.zeros(len(a),bool)
    for axis in axes:valid|=(a.trait_name.eq(axis)&a.normalized_value.isin(VALID[axis])).to_numpy()
    a=a[valid]
    v=a.pivot(index="accepted_species",columns="trait_name",values="normalized_value").reindex(columns=axes).dropna()
    q=a.pivot(index="accepted_species",columns="trait_name",values="evidence_quality").reindex(index=v.index,columns=axes)
    assert q.notna().all().all()
    v["genus"]=[genus(s) for s in v.index]
    v["grade_pattern"]=q.astype(str).agg("|".join,axis=1)
    return v[v.genus.notna()].sort_index()

def bootstrap(vals,labels,draws=1999):
    x=np.asarray(vals,float)
    if x.ndim==1:x=x[:,None]
    assert len(x)==len(labels)
    if len(x)==0:return {"islands":0,"groups":0,"mean":None,"interval":None,"tail_negative":None,"valid_draws":0,"nondegenerate":False}
    _,idx=np.unique(np.asarray(labels,str),return_inverse=True)
    k=int(idx.max())+1
    n=np.bincount(idx)
    a=np.column_stack([np.bincount(idx,weights=x[:,j])/n for j in range(x.shape[1])])
    mu=a.mean(0)
    if k<2:return {"islands":len(x),"groups":k,"mean":mu.tolist(),"interval":None,"tail_negative":None,"valid_draws":0,"nondegenerate":False}
    rng=np.random.default_rng(SEED)
    boot=np.stack([a[rng.integers(k,size=k)].mean(0) for _ in range(draws)])
    ci=np.quantile(boot,[.025,.975],axis=0).T
    tail=(1+(boot>=0).sum(0))/(draws+1)
    return {"islands":len(x),"groups":k,"mean":mu.tolist(),"interval":ci.tolist(),"tail_negative":tail.tolist(),
            "valid_draws":len(boot),"nondegenerate":bool(ci[0,1]-ci[0,0]>1e-12)}

def summarize(rows):
    ans={"eligible_islands":len(rows),"metric_order":list(METRICS)}
    if not rows:return ans|{"status":"not_estimable"}
    z=pd.DataFrame(rows)
    x=np.vstack(z.delta)
    assert np.allclose(x[:,0],x[:,1]+x[:,2],atol=1e-12)
    assert np.allclose(x[:,0],x[:,3]+x[:,4],atol=1e-12)
    g=bootstrap(x,z.block);s=bootstrap(x,z.first_source)
    support=all(v["groups"]>=30 and v["islands"]>=200 and v["valid_draws"]>=1800 and v["nondegenerate"] for v in (g,s))
    p=max(g["tail_negative"][0],s["tail_negative"][0]) if support else None
    bound=min(1.,4*p) if p is not None else None
    ok=bool(support and bound<.05 and g["interval"][0][1]<0 and s["interval"][0][1]<0)
    phi=z[z.phi.notna()]
    # Q1 was frozen before Q2 inspection. Formal H inference uses exactly
    # the original paired island cohort and BOTH grouping definitions.
    comp=np.column_stack([z.h_delta.to_numpy(float),np.vstack(z.score_delta),np.vstack(z.strict_delta),z.unclassified_delta])
    cg=bootstrap(comp,z.block);cs=bootstrap(comp,z.first_source)
    cg_support=all(v["groups"]>=30 and v["islands"]>=200 and v["valid_draws"]>=1800 and v["nondegenerate"] for v in (cg,cs))
    q1p=max(cg["tail_negative"][0],cs["tail_negative"][0]) if cg_support else None
    q1={"status":"exploratory_estimable" if cg_support else "insufficient_support",
        "endpoint_order":["specialist_minus_generalist",*["soft_"+k for k in TEMPLATES],*["strict_"+k for k in TEMPLATES],"unclassified"],
        "spatial_block":cg,"first_source_group":cs,
        "worst_directional_bootstrap_tail":q1p,
        "interim_four_question_p_upper_bound":min(1.,4*q1p) if q1p is not None else None,
        "interim_directional_candidate":bool(cg_support and 4*q1p<.05 and cg["interval"][0][1]<0 and cs["interval"][0][1]<0)}
    return ans|{"status":"exploratory_estimable" if support else "insufficient_support",
        "spatial_block":g,"first_source_group":s,
        "interim_four_question_p_upper_bound":bound,"interim_directional_candidate":ok,
        "median_target_scored":float(z.target_n.median()),"median_source_scored":float(z.source_n.median()),
        "median_shared_scored_names":float(z.shared_n.median()),
        "target_coherence_island_equal":float(z.target_coherence.mean()),
        "source_coherence_island_equal":float(z.source_coherence.mean()),
        "composition_q1_inference":q1,
        "composition_q1_descriptive_only":{"template_order":list(TEMPLATES),
            "target_minus_source_soft_scores":np.vstack(z.score_delta).mean(0).tolist(),
            "target_minus_source_strict_matches":np.vstack(z.strict_delta).mean(0).tolist(),
            "target_minus_source_h":float(z.h_delta.mean()),
            "target_minus_source_unclassified":float(z.unclassified_delta.mean())},
        "phi_common_pair_sensitivity":{"islands":len(phi),"spatial_block":bootstrap(phi.phi.to_numpy(float),phi.block),
                "first_source_group":bootstrap(phi.phi.to_numpy(float),phi.first_source)}}

def variant_stats(frame,variant):
    axes=CORE+(SIZE,) if variant=="size_expanded" else CORE
    X,pairs,weights,labels,score,strict=encode(frame,axes,red_only=variant=="strict_red")
    # Both genus and genus x exact per-axis grade patterns on the same carriers.
    gen=pd.factorize(frame.genus,sort=True)[0]
    gg=pd.factorize(pd.MultiIndex.from_frame(frame[["genus","grade_pattern"]]),sort=True)[0]
    h=score[:,:4].mean(1)-score[:,4]
    if variant=="size_expanded":
        from itertools import combinations
        is_size=[];core_w=[]
        for k,spec in TEMPLATES.items():
            a=[v for v in axes if v in spec]
            cp=list(combinations([v for v in a if v!=SIZE],2))
            for j,l in combinations(a,2):
                z=SIZE in (j,l);is_size.append(z)
                core_w.append(0. if z else 1/(len(TEMPLATES)*len(cp)))
        is_size=np.asarray(is_size,bool);core_w=np.asarray(core_w,float)
        assert len(is_size)==len(weights) and np.isclose(core_w.sum(),1)
        assert np.isclose(weights[is_size].sum()+weights[~is_size].sum(),1)
    def stat(ids):
        ix=np.asarray(sorted(ids),int)
        n=len(ix);assert n>=10
        a=components(X[ix],pairs,gg[ix]);b=components(X[ix],pairs,gen[ix])
        vals=np.array([weights@a["total"],weights@a["within"],weights@a["between"],weights@b["within"],weights@b["between"]])
        out={"metrics":vals,"phi":a["phi"],"score":score[ix].mean(0),"strict":strict[ix].mean(0),
                "h":float(h[ix].mean()),"unclassified":float((strict[ix].sum(1)==0).mean()),"n":n}
        if variant=="size_expanded":
            mat=np.stack([a["total"],a["within"],a["between"],b["within"],b["between"]])
            c0=mat@core_w
            c1=mat@(weights*(~is_size))
            c2=mat@(weights*is_size)
            assert np.allclose(c1+c2,vals,atol=1e-12)
            assert np.allclose(vals-c0,(c1-c0)+c2,atol=1e-12)
            out["size_parts"]=np.stack([c0,c1,c2,c1-c0,vals-c0])
            # For each template, isolate the new size-pair contribution
            # under the expanded weights. All five templates are returned.
            out["size_by_template"]=np.array([mat[0]@(weights*is_size*(labels==k)) for k in TEMPLATES])
        return out
    return stat,labels

def analyze(f,t,a,m,g):
    full=f.groupby("island_id").accepted_species.agg(set).to_dict()
    region=m.groupby("entity_ID").work_species.agg(set).to_dict()
    keys=a.sort_values("source_rank").groupby(["island_id","source_mode"]).entity_ID.agg(tuple).to_dict()
    geo=g.set_index("island_id")
    # Only geographic groups used, no individual coordinates in reports.
    geoc={}
    for i in sorted(full.keys()&set(geo.index)):
        z=geo.loc[i];lat=pd.to_numeric(z.island_latitude,errors="coerce");lon=pd.to_numeric(z.island_longitude,errors="coerce")
        if not np.isfinite([lat,lon]).all() or pd.isna(z.analysis_regime):continue
        geoc[i]=(f"{math.floor(lat/10)}_{math.floor(lon/10)}",str(z.analysis_regime))
    panels={}
    for variant in VARIANTS:
        frame=trait_table(t,variant);lookup={s:i for i,s in enumerate(frame.index)}
        stat,labels=variant_stats(frame,variant)
        f0=f[f.origin_status.eq("native")] if variant=="original_native" else f
        tar={i:{lookup[s] for s in v if s in lookup} for i,v in f0.groupby("island_id").accepted_species.agg(set).items()}
        all_idx={i:{lookup[s] for s in v if s in lookup} for i,v in full.items()}
        @functools.lru_cache(maxsize=512)
        def pool(ids):return frozenset(lookup[s] for e in ids for s in region[e] if s in lookup)
        @functools.lru_cache(maxsize=512)
        def inclstat(ids):return stat(pool(ids))
        for mode,nr in MODES.items():
            rows={"inclusive":[],"disjoint":[]};missing_assignment=missing_donors=0
            target_geographic_eligible=0;donor_inclusive_eligible=0;common_eligible=0
            for i in sorted(geoc):
                tr=tar.get(i,set())
                if len(tr)<10:continue
                target_geographic_eligible+=1
                ids=keys.get((i,mode),())
                if not ids:missing_assignment+=1;continue
                assert len(ids)==nr
                if any(e not in region for e in ids):missing_donors+=1;continue
                inc=pool(ids);exc=inc-all_idx[i]
                if len(inc)<10:continue
                donor_inclusive_eligible+=1
                if len(exc)<10:continue
                common_eligible+=1
                assert not tr&exc
                u=stat(tr)
                for tag,roster,v in (("inclusive",inc,inclstat(ids)),("disjoint",exc,stat(exc))):
                    block,reg=geoc[i]
                    rows[tag].append({"block":block,"regime":reg,"first_source":ids[0],
                        "delta":u["metrics"]-v["metrics"],"target_n":u["n"],"source_n":v["n"],
                        "shared_n":len(tr&roster),"phi":common_phi(u["phi"],v["phi"],labels),
                        "target_coherence":u["metrics"][0],"source_coherence":v["metrics"][0],
                        "score_delta":u["score"]-v["score"],"strict_delta":u["strict"]-v["strict"],
                        "h_delta":u["h"]-v["h"],"unclassified_delta":u["unclassified"]-v["unclassified"]})
                    if variant=="size_expanded":
                        rows[tag][-1]["size_parts"]=u["size_parts"]-v["size_parts"]
                        rows[tag][-1]["size_by_template"]=u["size_by_template"]-v["size_by_template"]
            assert len(rows["inclusive"])==len(rows["disjoint"])==common_eligible
            for tag,rr in rows.items():
                key=f"{variant}_{mode}_{tag}"
                p=summarize(rr)
                if variant=="size_expanded" and rr:
                    parts=np.stack([x["size_parts"] for x in rr])
                    bt=np.stack([x["size_by_template"] for x in rr])
                    assert np.allclose(parts[:,4],parts[:,3]+parts[:,2],atol=1e-12)
                    assert np.allclose(bt.sum(1),parts[:,2,0],atol=1e-12)
                    X=np.column_stack([parts.reshape(len(rr),-1),bt])
                    p["size_decomposition"]={"part_order":["core_original_weights","core_expanded_weights","size_pairs_expanded_weights","core_reweighting","expanded_minus_core"],
                        "component_order":list(METRICS),"size_template_order":list(TEMPLATES),
                        "spatial_block":bootstrap(X,[x["block"] for x in rr]),
                        "first_source_group":bootstrap(X,[x["first_source"] for x in rr]),
                        "warning":"Post-outcome paired algebraic diagnostic; size contribution and reweighting are not causal effects."}
                p["support_census"]={"complete_trait_species":len(frame),"target_geographic_eligible":target_geographic_eligible,
                    "missing_assignment":missing_assignment,"missing_candidate_flora":missing_donors,
                    "inclusive_donor_eligible":donor_inclusive_eligible,"paired_donor_eligible":common_eligible}
                p["regional_descriptive"]={}
                for reg in sorted({v[1] for v in geoc.values()}):
                    v=[x for x in rr if x["regime"]==reg]
                    if not v:p["regional_descriptive"][reg]={"eligible_islands":0}
                    else:
                        p["regional_descriptive"][reg]={"eligible_islands":len(v),
                            "total_delta":bootstrap(np.array([x["delta"][0] for x in v]),[x["block"] for x in v]),
                            "island_equal_components":np.vstack([x["delta"] for x in v]).mean(0).tolist()}
                panels[key]=p
                # Only scalar aggregate support and effect; no roster or group identifiers.
                print("PLANT_Q2_PANEL",key,json.dumps(p,allow_nan=False),flush=True)
        pool.cache_clear();inclstat.cache_clear()
    assert len(panels)==24
    for mode in MODES:
        for tag in ("inclusive","disjoint"):
            p0=panels[f"size_core_{mode}_{tag}"];p1=panels[f"size_expanded_{mode}_{tag}"]
            assert p0["eligible_islands"]==p1["eligible_islands"]
            if p0["eligible_islands"]:
                for gr in ("spatial_block","first_source_group"):
                    assert np.allclose(p1["size_decomposition"][gr]["mean"][:5],p0[gr]["mean"],atol=1e-12)
                    assert np.allclose(np.array(p1["size_decomposition"][gr]["mean"][20:25]),np.array(p1[gr]["mean"])-np.array(p0[gr]["mean"]),atol=1e-12)
    return panels

def self_test():
    core_test()
    rng=np.random.default_rng(SEED)
    tab=pd.DataFrame({axis:rng.choice(sorted(VALID[axis]),120) for axis in CORE})
    tab["genus"]=np.repeat(np.arange(30),4).astype(str);tab["grade_pattern"]="high|high|high|high"
    stat,labels=variant_stats(tab,"all_recorded")
    a=stat(set(range(60)));b=stat(set(range(60,120)))
    assert np.isclose(a["metrics"][0],sum(a["metrics"][1:3]))
    assert common_phi(a["phi"],b["phi"],labels) is not None
    x=bootstrap(np.full(240,-.1)+rng.normal(0,.01,240),np.repeat(np.arange(40),6),draws=1999)
    assert x["groups"]==40 and x["interval"][0][1]<0 and x["tail_negative"][0]<.01
    print("PLANT_Q2_ORIGINAL_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--archive",type=pathlib.Path);p.add_argument("--geography",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:self_test();return
    if not all((a.archive,a.geography,a.out)):p.error("archive, geography and out are required")
    try:
        f,t,s,m,g,md=read_archive(a.archive,a.geography)
        panels=analyze(f,t,s,m,g)
        out={"original_botanical_parent":4505,"original_geography_parent":8264,"mainland_candidate_regions":748,
            "archive_sha256":ARCHIVE_SHA,"mainland_sha256":md,"source_file_pins":PIN|{"geography":GEO_SHA},
            "primary_q2":PRIMARY,"panels":panels,"holm_four_primary_questions":"pending_other_questions",
            "warning":"Plant morphology proxies only; candidate mainland floras are not ancestors. Q2 signed excess covariance depends on marginal variance; phi sensitivity uses common nonconstant pairs. Bootstrap tails exploratory; source and geographic blocks are not independent evolutionary lineages. No individual botanical/geographic information exported."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    except Exception as e:
        # Never expose source row values or paths from parser exceptions.
        print("PLANT_Q2_ORIGINAL_FAILED",type(e).__name__,flush=True)
        sys.exit(1)
if __name__=="__main__":main()
