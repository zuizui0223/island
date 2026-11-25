#!/usr/bin/env python3
"""Post-atlas diagnostics: same-island origin-status mixture and colour definition.

Primary atlas is frozen. This is a follow-up on colour/size and form/size,
not independent confirmation. Compare native and full rosters on exactly
the same islands, candidate sources, and scoring definition.
"""
import argparse,functools,json,math,pathlib,sys
import numpy as np
import pandas as pd
from plant_syndrome_core import CORE,SIZE,TEMPLATES
from plant_syndrome_pair_atlas import pair_table
from plant_syndrome_q2_original import read_archive,bootstrap,MODES,ARCHIVE_SHA,PIN,GEO_SHA

PAIRS=((CORE[0],SIZE),(CORE[1],SIZE))
STATUS=("native","introduced","unresolved")
PARTS=("all_minus_source","native_minus_source","all_minus_native",
       "introduced_within_status_mixture","unresolved_within_status_mixture",
       "between_status_covariance","finite_n_correction")
METRICS=("specialist_bee","bird","bat","open_generalist","bird_minus_others")

def vectors(frame,pair,variant):
    names=[k for k,v in TEMPLATES.items() if all(j in v for j in pair)]
    assert names==list(METRICS[:4])
    X=[];Y=[]
    for k in names:
        spec=TEMPLATES[k]
        states=set(spec[pair[0]].split())
        if k=="bird" and pair[0]==CORE[0]:
            if variant=="red_pink_only":states={"red_pink"}
            if variant=="yellow_orange_only":states={"yellow_orange"}
        X.append(frame[pair[0]].isin(states).to_numpy(float))
        Y.append(frame[SIZE].isin(set(spec[SIZE].split())).to_numpy(float))
    return np.column_stack(X),np.column_stack(Y)

def augment(a):
    # Equal weight for the three comparison templates.
    return np.append(a,a[1]-np.mean(a[[0,2,3]]))

def population(x,y,ix,weights=None):
    z=x[ix];w=y[ix]
    if weights is None:p=np.repeat(1/len(ix),len(ix))
    else:
        p=np.asarray(weights,float);assert (p>0).all()
        p=p/p.sum()
    a=p@z;b=p@w
    return p@(z*w)-a*b,a,b

def mixture(x,y,sets):
    ids=[np.array(sorted(s),int) for s in sets]
    assert len(ids[0])>=10
    assert not(sets[0]&sets[1] or sets[0]&sets[2] or sets[1]&sets[2])
    full=np.concatenate(ids);n=len(full);nn=len(ids[0])
    w=np.array([len(i)/n for i in ids])
    C,A,B=population(x,y,full)
    cn,an,bn=population(x,y,ids[0])
    f=n/(n-1);fn=nn/(nn-1)
    within=np.zeros((3,x.shape[1]))
    between=np.zeros(x.shape[1])
    for j,ix in enumerate(ids):
        if not len(ix):continue
        cj,aj,bj=population(x,y,ix)
        within[j]=f*w[j]*(cj-cn)
        between+=f*w[j]*(aj-A)*(bj-B)
    finite=(f-fn)*cn
    all_c=f*C;native_c=fn*cn
    assert np.allclose(all_c-native_c,within.sum(0)+between+finite,atol=1e-12)
    return all_c,native_c,within[1],within[2],between,finite,w

def genus_equal_cov(x,y,ids,genus):
    ix=np.asarray(sorted(ids),int)
    _,g=np.unique(genus[ix],return_inverse=True)
    n=np.bincount(g)
    # One vote per genus, equal name weights within genus.
    return population(x,y,ix,weights=1/n[g])[0]

def analyze(f,t,a,m,g):
    assert f.origin_status.isin(STATUS).all()
    full=f.groupby("island_id").accepted_species.agg(set).to_dict()
    origin={k:f[f.origin_status.eq(k)].groupby("island_id").accepted_species.agg(set).to_dict() for k in STATUS}
    source=m.groupby("entity_ID").work_species.agg(set).to_dict()
    keys=a.sort_values("source_rank").groupby(["island_id","source_mode"]).entity_ID.agg(tuple).to_dict()
    geoc={}
    geo=g.set_index("island_id")
    for i in sorted(full.keys()&set(geo.index)):
        z=geo.loc[i];lat=pd.to_numeric(z.island_latitude,errors="coerce");lon=pd.to_numeric(z.island_longitude,errors="coerce")
        if np.isfinite([lat,lon]).all() and pd.notna(z.analysis_regime):
            geoc[i]=(f"{math.floor(lat/10)}_{math.floor(lon/10)}",str(z.analysis_regime))
    panels={}
    for pair in PAIRS:
        frame=pair_table(t,pair);idx={s:i for i,s in enumerate(frame.index)}
        gen=pd.factorize(frame.genus,sort=True)[0]
        all_idx={i:{idx[s] for s in v if s in idx} for i,v in full.items()}
        groups={k:{i:{idx[s] for s in v if s in idx} for i,v in origin[k].items()} for k in STATUS}
        variants=("original","red_pink_only","yellow_orange_only") if pair[0]==CORE[0] else ("original",)
        @functools.lru_cache(maxsize=512)
        def pool(ids):return frozenset(idx[s] for e in ids for s in source[e] if s in idx)
        for mode,nr in MODES.items():
            eligible=[]
            census={"complete_pair_species":len(frame),"native_ge10":0,"common_donor_ge10":0}
            for i,(block,reg) in geoc.items():
                parts=[groups[k].get(i,set()) for k in STATUS]
                if len(parts[0])<10:continue
                census["native_ge10"]+=1
                assert all_idx[i]==set().union(*parts)
                ids=keys.get((i,mode),())
                if not ids:continue
                assert len(ids)==nr
                if any(k not in source for k in ids):continue
                inc=pool(ids);exc=inc-all_idx[i]
                if len(exc)<10:continue
                census["common_donor_ge10"]+=1
                eligible.append((i,parts,block,reg,ids[0],inc,exc))
            for variant in variants:
                x,y=vectors(frame,pair,variant)
                rr={"inclusive":[],"disjoint":[]}
                for i,parts,block,reg,first,inc,exc in eligible:
                    ac,nc,intro,unknown,between,finite,w=mixture(x,y,parts)
                    full_i=all_idx[i]
                    ap=population(x,y,np.array(sorted(full_i)))[0]
                    npop=population(x,y,np.array(sorted(parts[0])))[0]
                    ag=genus_equal_cov(x,y,full_i,gen);ng=genus_equal_cov(x,y,parts[0],gen)
                    for tag,roster in (("inclusive",inc),("disjoint",exc)):
                        sp=population(x,y,np.array(sorted(roster)))[0]
                        sc=sp*len(roster)/(len(roster)-1)
                        sg=genus_equal_cov(x,y,roster,gen)
                        val=np.stack(list(map(augment,[ac-sc,nc-sc,ac-nc,intro,unknown,between,finite])))
                        assert np.allclose(val[0]-val[1],val[2],atol=1e-12)
                        assert np.allclose(val[2],val[3:].sum(0),atol=1e-12)
                        sensitivity=np.stack(list(map(augment,[ap-sp,npop-sp,ag-sg,ng-sg])))
                        rr[tag].append({"block":block,"first_source":first,"regime":reg,"val":val,
                           "sensitivity":sensitivity,"fractions":w,"native_n":len(parts[0]),"all_n":len(full_i),
                           "native_idx":parts[0],"all_idx":full_i,"source_idx":roster})
                for tag,rows in rr.items():
                    key="__".join((*pair,variant,mode,tag))
                    out={"eligible_islands":len(rows),"support_census":census,
                         "part_order":list(PARTS),"metric_order":list(METRICS),
                         "sensitivity_order":["all_species_equal_population","native_species_equal_population",
                                              "all_genus_equal_population","native_genus_equal_population"]}
                    if not rows:
                        panels[key]=out|{"status":"not_estimable"}
                        print("PLANT_ORIGIN_SIZE_PANEL",key,json.dumps(panels[key]),flush=True);continue
                    bg=[r["block"] for r in rows];sg=[r["first_source"] for r in rows]
                    v=np.stack([r["val"] for r in rows]);s=np.stack([r["sensitivity"] for r in rows])
                    bg0=bootstrap(v.reshape(len(v),-1),bg);sg0=bootstrap(v.reshape(len(v),-1),sg)
                    gate=all(z["groups"]>=30 and z["islands"]>=200 and z["valid_draws"]>=1800 for z in (bg0,sg0))
                    out|={"status":"exploratory_estimable" if gate else "insufficient_support",
                           "spatial_block":bg0,"first_source_group":sg0,
                           "population_and_genus_equal":{"spatial_block":bootstrap(s.reshape(len(s),-1),bg),
                                                       "first_source_group":bootstrap(s.reshape(len(s),-1),sg)},
                           "median_native_scored":float(np.median([r["native_n"] for r in rows])),
                           "median_all_scored":float(np.median([r["all_n"] for r in rows])),
                           "distinct_native_taxa":len(set().union(*(r["native_idx"] for r in rows))),
                           "distinct_all_taxa":len(set().union(*(r["all_idx"] for r in rows))),
                           "distinct_source_taxa":len(set().union(*(r["source_idx"] for r in rows))),
                           "origin_status_order":list(STATUS),
                           "status_fraction_island_equal":np.stack([r["fractions"] for r in rows]).mean(0).tolist(),
                           "regional_descriptive":{}}
                    for reg in sorted(set(r["regime"] for r in rows)):
                        z=[r for r in rows if r["regime"]==reg]
                        out["regional_descriptive"][reg]=bootstrap(np.stack([r["val"][:3,-1] for r in z]),[r["block"] for r in z])
                    panels[key]=out
                    print("PLANT_ORIGIN_SIZE_PANEL",key,json.dumps(out,allow_nan=False),flush=True)
            # All definitions use the same target and donor cohort.
            for tag in ("inclusive","disjoint"):
                assert len({panels["__".join((*pair,v,mode,tag))]["eligible_islands"] for v in variants})==1
        pool.cache_clear()
    assert len(panels)==16
    return panels

def self_test():
    rng=np.random.default_rng(20261125)
    x=rng.integers(0,2,(40,4));y=rng.integers(0,2,(40,4))
    sets=[set(range(20)),set(range(20,30)),set(range(30,40))]
    a,n,i,u,b,c,w=mixture(x,y,sets)
    assert np.allclose(a-n,i+u+b+c)
    assert np.allclose(mixture(x,y,[sets[0],set(),set()])[0],n)
    assert np.isfinite(genus_equal_cov(x,y,sets[0],np.repeat(np.arange(20),2))).all()
    print("PLANT_ORIGIN_SIZE_MIXTURE_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--archive",type=pathlib.Path);p.add_argument("--geography",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    try:
        f,t,s,m,g,md=read_archive(a.archive,a.geography)
        panels=analyze(f,t,s,m,g)
        out={"original_botanical_parent":4505,"original_geography_parent":8264,"mainland_candidate_regions":748,
            "archive_sha256":ARCHIVE_SHA,"mainland_sha256":md,"source_file_pins":PIN|{"geography":GEO_SHA},
            "panels":panels,"warning":"Post-atlas diagnostic, same-island full/native source comparison. Exact algebraic mixture components are not causal contributions of introduction. Unresolved origin is not nonnative. Red/pink and yellow/orange are coarse inseparable categories, no actual pollinator labels. Genus-equal sensitivity uses population covariance and same species; shared taxa across islands remain dependent. No protected identifiers."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    except Exception as e:
        print("PLANT_ORIGIN_SIZE_FAILED",type(e).__name__,flush=True);sys.exit(1)
if __name__=="__main__":main()
