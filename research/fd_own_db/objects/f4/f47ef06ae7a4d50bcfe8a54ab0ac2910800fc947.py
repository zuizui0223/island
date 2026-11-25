#!/usr/bin/env python3
"""Focused colour-size diagnostic: paired colour specificity and coherent genus weights.

Post-outcome exploration. Same native-eligible islands and identical candidate
mainland pools as the origin-size diagnostic. Population covariance for all
weighted estimands; no species or genus identities in output. Exponential
genus weights are a sensitivity to shared taxonomy, not a calibrated sampling
distribution or phylogenetic model.
"""
import argparse,functools,json,math,pathlib,sys
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from plant_syndrome_core import CORE,SIZE,components
from plant_syndrome_pair_atlas import pair_table
from plant_syndrome_origin_size_diagnostic import vectors
from plant_syndrome_q2_original import read_archive,MODES,ARCHIVE_SHA,PIN,GEO_SHA,SEED

PAIR=(CORE[0],SIZE)
RAW=("bee","bird_broad","bat","generalist","red_pink_large","yellow_orange_large")
ENDS=RAW+("red_minus_other_templates","yellow_minus_other_templates","red_minus_yellow")
GROUPINGS=("spatial_block","first_source_group")
COHORTS=("all_on_native_eligible_islands","literal_native")
DRAWS=999

def endpoints(c):
    other=c[..., [0,2,3]].mean(axis=-1)
    return np.concatenate([c, (c[...,4]-other)[...,None],
        (c[...,5]-other)[...,None],(c[...,4]-c[...,5])[...,None]],axis=-1)

def membership(rosters,p):
    rr=[];cc=[]
    for j,s in enumerate(rosters):
        rr.extend([j]*len(s));cc.extend(sorted(s))
    return csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(rosters),p))

def cov(A,w,stat):
    mass=A@w
    assert (mass>0).all()
    q=(A@(w[:,None]*stat))/mass[:,None]
    p=q.shape[1]//3
    c=q[:,2*p:]-q[:,:p]*q[:,p:2*p]
    assert np.allclose(c[:,1],c[:,4]+c[:,5],atol=1e-12)
    return endpoints(c)

def grouping(labels):
    _,g=np.unique(np.asarray(labels,str),return_inverse=True)
    n=np.bincount(g)
    # group-equal average of island contrasts, not species records
    W=csr_matrix((1/n[g],(g,np.arange(len(g)))),shape=(len(n),len(g)))
    return W,len(n)

def interval(draw,point):
    a=np.asarray(draw,float)
    return {"mean":point.tolist(),"interval":np.quantile(a,[.025,.975],axis=0).T.tolist(),
            "draws":len(a)}

def test_cohort(tar,src,stat,x,y,gen,grade,labels,draws=DRAWS):
    A=membership(tar,len(gen));M=membership(src,len(gen))
    w=np.ones(len(gen));delta=cov(A,w,stat)-cov(M,w,stat)
    # Exact within/between-genus-grade population covariance partition.
    pairs=np.column_stack([np.arange(x.shape[1]),np.arange(x.shape[1])+x.shape[1]])
    z=np.column_stack([x,y])
    def partition(roster):
        out=[]
        for s in roster:
            ix=np.asarray(sorted(s),int)
            q=components(z[ix],pairs,grade[ix])
            scale=(len(ix)-1)/len(ix)
            out.append(np.stack([endpoints(q[k]*scale) for k in ("total","within","between")]))
        return np.stack(out)
    parts=partition(tar)-partition(src)
    assert np.allclose(parts[:,0,:],delta,atol=1e-12)
    assert np.allclose(parts[:,1,:]+parts[:,2,:],delta,atol=1e-12)
    gr=[grouping(v) for v in labels]
    result={}
    all_taxa=set().union(*tar,*src)
    used=np.asarray(sorted(all_taxa),int)
    taxon_groups=np.unique(gen[used])
    # One positive random weight per genus, common to all accepted names,
    # all repeated occurrences, and both source and target.
    rng=np.random.default_rng(SEED)
    tax=[[] for _ in gr];joint=[[] for _ in gr];spatial=[[] for _ in gr]
    for j in range(draws):
        gw=rng.exponential(size=int(gen.max())+1)
        v=cov(A,gw[gen],stat)-cov(M,gw[gen],stat)
        for k,(W,G) in enumerate(gr):
            q=W@v;z0=W@delta
            sampled=rng.integers(G,size=G)
            tax[k].append(q.mean(0))
            joint[k].append(q[sampled].mean(0))
            spatial[k].append(z0[sampled].mean(0))
    for k,(W,G) in enumerate(gr):
        p=(W@delta).mean(0)
        flat=parts.reshape(len(parts),-1)
        decomp=(W@flat).mean(0).reshape(3,-1)
        assert np.allclose(decomp[0],p,atol=1e-12)
        result[GROUPINGS[k]]={"groups":G,"spatial_only":interval(spatial[k],p),
            "genus_weight_only":interval(tax[k],p),
            "joint_group_and_genus_weight":interval(joint[k],p),
            "decomposition_order":["total","within_genus_grade","between_genus_grade"],
            "decomposition_means":decomp.tolist()}
    return {"islands":len(tar),"distinct_taxa":len(used),"distinct_genera":len(taxon_groups),
            "groupings":result,
            "support_gate":bool(len(tar)>=200 and len(taxon_groups)>=30 and all(G>=30 for _,G in gr))}

def analyze(f,t,a,m,g):
    frame=pair_table(t,PAIR);lookup={s:i for i,s in enumerate(frame.index)}
    assert len(frame)==3676
    bx,by=vectors(frame,PAIR,"original")
    rx,ry=vectors(frame,PAIR,"red_pink_only")
    yx,yy=vectors(frame,PAIR,"yellow_orange_only")
    x=np.column_stack([bx,rx[:,1],yx[:,1]]);y=np.column_stack([by,ry[:,1],yy[:,1]])
    assert np.array_equal(x[:,1],x[:,4]+x[:,5]) and np.array_equal(y[:,1],y[:,4]) and np.array_equal(y[:,1],y[:,5])
    stat=np.column_stack([x,y,x*y])
    gen=pd.factorize(frame.genus,sort=True)[0]
    grade=pd.factorize(pd.MultiIndex.from_frame(frame[["genus","grade_pattern"]]),sort=True)[0]
    full=f.groupby("island_id").accepted_species.agg(set).to_dict()
    native=f[f.origin_status.eq("native")].groupby("island_id").accepted_species.agg(set).to_dict()
    full_idx={i:{lookup[s] for s in v if s in lookup} for i,v in full.items()}
    nat_idx={i:{lookup[s] for s in v if s in lookup} for i,v in native.items()}
    regions=m.groupby("entity_ID").work_species.agg(set).to_dict()
    keys=a.sort_values("source_rank").groupby(["island_id","source_mode"]).entity_ID.agg(tuple).to_dict()
    geo=g.set_index("island_id")
    geoc={}
    for i in sorted(full.keys()&set(geo.index)):
        z=geo.loc[i];lat=pd.to_numeric(z.island_latitude,errors="coerce");lon=pd.to_numeric(z.island_longitude,errors="coerce")
        if np.isfinite([lat,lon]).all() and pd.notna(z.analysis_regime):
            geoc[i]=f"{math.floor(lat/10)}_{math.floor(lon/10)}"
    @functools.lru_cache(maxsize=512)
    def pool(ids):return frozenset(lookup[s] for j in ids for s in regions[j] if s in lookup)
    panels={}
    for mode,nr in MODES.items():
        rows=[]
        for i,block in geoc.items():
            nat=nat_idx.get(i,set())
            if len(nat)<10:continue
            ids=keys.get((i,mode),())
            if not ids:continue
            assert len(ids)==nr
            if any(e not in regions for e in ids):continue
            inc=pool(ids);exc=inc-full_idx[i]
            if len(exc)<10:continue
            assert nat<=full_idx[i] and not nat&exc
            rows.append((full_idx[i],nat,inc,exc,block,ids[0]))
        assert len(rows)==(226 if mode=="geo_k5" else 227)
        labels=[[v[4] for v in rows],[v[5] for v in rows]]
        for donor,j in (("inclusive",2),("disjoint",3)):
            src=[v[j] for v in rows]
            for cohort,k in zip(COHORTS,(0,1)):
                p=test_cohort([v[k] for v in rows],src,stat,x,y,gen,grade,labels)
                p["endpoint_order"]=list(ENDS)
                key="__".join((mode,donor,cohort))
                panels[key]=p
                print("PLANT_COLOR_TAXON_PANEL",key,json.dumps(p,allow_nan=False),flush=True)
    assert len(panels)==8
    return panels

def self_test():
    rng=np.random.default_rng(SEED);n=60
    red=rng.integers(0,2,n);yellow=(1-red)*rng.integers(0,2,n)
    x=rng.integers(0,2,(n,6)).astype(float);x[:,4]=red;x[:,5]=yellow;x[:,1]=red+yellow
    y=rng.integers(0,2,(n,6)).astype(float);y[:,4]=y[:,1];y[:,5]=y[:,1]
    stat=np.column_stack([x,y,x*y]);g=np.repeat(np.arange(30),2)
    tar=[set(rng.choice(n,20,replace=False)) for _ in range(30)]
    src=[set(rng.choice(n,20,replace=False)) for _ in range(30)]
    q=test_cohort(tar,src,stat,x,y,g,g,[np.arange(30)%6,np.arange(30)%7],draws=19)
    assert q["islands"]==30 and q["groupings"]["spatial_block"]["joint_group_and_genus_weight"]["draws"]==19
    q=test_cohort(tar,tar,stat,x,y,g,g,[np.arange(30)%6,np.arange(30)%7],draws=19)
    assert np.allclose(q["groupings"]["spatial_block"]["joint_group_and_genus_weight"]["interval"],0)
    print("PLANT_COLOR_COHERENT_GENUS_WEIGHT_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--archive",type=pathlib.Path);p.add_argument("--geography",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    try:
        f,t,s,m,g,md=read_archive(a.archive,a.geography)
        panels=analyze(f,t,s,m,g)
        out={"original_botanical_parent":4505,"original_geography_parent":8264,"mainland_candidate_regions":748,
             "archive_sha256":ARCHIVE_SHA,"mainland_sha256":md,"source_file_pins":PIN|{"geography":GEO_SHA},
             "panels":panels,"warning":"Post-outcome exploratory colour-size specificity sensitivity. Red/pink and yellow/orange are observed categories, not pollinator guilds. Coherent positive genus weights and paired geographic/source group resampling are sensitivity distributions, not calibrated confidence intervals. Same native-eligible islands, shared species and genus labels; no independent phylogenetic contrast or causal mechanism claim."}
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    except Exception as e:
        print("PLANT_COLOR_TAXON_FAILED",type(e).__name__,flush=True);sys.exit(1)
if __name__=="__main__":main()
