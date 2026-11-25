#!/usr/bin/env python3
"""Post-result FD diagnostic: factor joint deep-flower composition and
compare thermal regions on one common eligible geography cohort.
No change to prior 64-panel test or protected data custody.
"""
import argparse, collections, functools, hashlib, io, json, pathlib
import numpy as np
import pandas as pd
from scipy import sparse
from fd_depth_conditional_protected import Stats,table,membership,group_matrix,read_pin,PIN,MODES,SEED,MIN_EACH
from fd_depth_conditional_architecture import SYMMETRY

DRAWS=1999
REGIONS=("north","tropical")
SCOPES=("all_recorded","original_native","high_quality_only")
VARIANTS=[(m,v) for m in MODES for v in ("inclusive","disjoint")]
TERMS=("joint_delta","depth_prevalence_component","conditional_symmetry_component","depth_prevalence_delta","conditional_symmetry_delta","nondeep_conditional_delta","deep_minus_nondeep")
NAMES=[f"{k}__{s}" for k in TERMS for s in SYMMETRY]
OUT_NAMES=("north","tropical","tropical_minus_north")
GROUPINGS=("geographic_block","first_source_group")

def factors(T,S):
    T=np.asarray(T,float).reshape(-1,3,4)
    S=np.asarray(S,float).reshape(-1,3,4)
    T=T/T.sum((1,2))[:,None,None];S=S/S.sum((1,2))[:,None,None]
    d=T[:,:,3].sum(1);e=S[:,:,3].sum(1)
    if np.any((d<=0)|(d>=1)|(e<=0)|(e>=1)):raise ValueError("undefined conditional denominator")
    q=T[:,:,3]/d[:,None];v=S[:,:,3]/e[:,None]
    u=T[:,:,:3].sum(2)/(1-d[:,None]);z=S[:,:,:3].sum(2)/(1-e[:,None])
    dj=T[:,:,3]-S[:,:,3]
    prevalence=(d-e)[:,None]*(q+v)/2
    composition=(d+e)[:,None]*(q-v)/2
    assert np.allclose(dj,prevalence+composition,atol=1e-12)
    assert np.allclose(composition.sum(1),0,atol=1e-12)
    return np.column_stack((dj,prevalence,composition,np.repeat((d-e)[:,None],3,axis=1),q-v,u-z,(q-v)-(u-z)))

def paired_group_design(rows,group):
    keys=sorted({r[group] for r in rows});lookup={k:i for i,k in enumerate(keys)}
    rids=[];mat=[]
    for reg in REGIONS:
        ix=np.array([j for j,r in enumerate(rows) if r["region"]==reg],int)
        if not len(ix):return None
        G=group_matrix([rows[j][group] for j in ix])
        labs=sorted({rows[j][group] for j in ix})
        rids.append(np.array([lookup[k] for k in labs],int))
        mat.append((ix,G))
    return keys,rids,mat

def get_means(Y,design,w):
    keys,rids,mat=design
    out=[]
    for ids,(ix,H) in zip(rids,mat):
        a=w[ids]
        if a.sum()==0:return None
        out.append(np.average(H@Y[ix],axis=0,weights=a))
    return np.stack([*out,out[1]-out[0]])

def quant(v):
    return np.quantile(v,[.025,.975],axis=0).transpose(1,2,0).tolist() if len(v)>=.9*DRAWS else None

def analyze(rows,st):
    if not rows:return {"status":"not_estimable","n":0}
    A=membership([r["t"] for r in rows],len(st.tab))
    B=membership([r["s"] for r in rows],len(st.tab))
    def calc(w):
        V=w[st.g,None]*st.V
        return factors(A@V,B@V)
    obs=calc(np.ones(st.G))
    groups={}
    for tag in GROUPINGS:
        d=paired_group_design(rows,tag)
        if d is not None:groups[tag]=d
    if not groups:return {"status":"not_estimable","n":len(rows)}
    rng=np.random.default_rng(SEED)
    points={k:get_means(obs,d,np.ones(len(d[0]))) for k,d in groups.items()}
    raw={k:[] for k in groups};tax={k:[] for k in groups};joint={k:[] for k in groups}
    for j in range(DRAWS):
        Y=calc(np.maximum(rng.exponential(size=st.G),1e-12))
        for k,d in groups.items():
            n=len(d[0]);w=np.bincount(rng.integers(n,size=n),minlength=n)
            for coll,v in ((raw,get_means(obs,d,w)),(tax,get_means(Y,d,np.ones(n))),(joint,get_means(Y,d,w))):
                if v is not None:coll[k].append(v)
    n_by={r:sum(x["region"]==r for x in rows) for r in REGIONS}
    out={}
    for k,d in groups.items():
        support={reg:{"islands":n_by[reg],"groups":len(ids)} for reg,ids in zip(REGIONS,d[1])}
        out[k]={"support":support,"groups_union":len(d[0]),"groups_shared":len(set(d[1][0])&set(d[1][1])),
                "point":points[k].tolist(),"interval":quant(raw[k]),"genus_weight_interval":quant(tax[k]),
                "joint_genus_group_interval":quant(joint[k]),"draws":[len(raw[k]),len(tax[k]),len(joint[k])]}
    gate=len(groups)==2 and all(v["islands"]>=200 and v["groups"]>=30 for a in out.values() for v in a["support"].values()) and all(min(a["draws"])>=1800 for a in out.values())
    return {"status":"exploratory_estimable" if gate else "insufficient_support","n":len(rows),"regions":n_by,"groupings":out}

def self_test():
    T=np.zeros((1,3,4));S=np.zeros_like(T)
    T[0,0,3]=30;T[0,1,3]=10;T[0,0,0]=60
    S[0,0,3]=10;S[0,1,3]=10;S[0,0,0]=80
    z=factors(T,S)[0]
    assert np.isclose(z[0],.2)
    assert np.isclose(z[3],.125)
    assert np.isclose(z[6],.075)
    assert np.isclose(z[12],.25)
    assert np.isclose(z[18],.25)
    assert np.allclose(factors(S,T),-factors(T,S))
    rows=[{"region":"north","geographic_block":"a"},{"region":"north","geographic_block":"b"},{"region":"tropical","geographic_block":"b"}]
    a=paired_group_design(rows,"geographic_block");assert np.allclose(get_means(np.array([[1.],[3.],[5.]]),a,np.ones(2)).ravel(),[2,5,3])
    print("FD_DEPTH_FACTOR_SYNTHETIC_PASS",flush=True)

def inputs(root,geofile):
    f=pd.read_csv(io.BytesIO(read_pin(root/"fixed/canonical/input/chapter1_status_flora.csv.gz","flora")),compression="gzip",dtype=str,usecols=["island_id","accepted_species","origin_status"])
    t=pd.read_csv(io.BytesIO(read_pin(root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz","traits")),compression="gzip",dtype=str)
    ass=pd.read_csv(io.BytesIO(read_pin(root/"fixed/source/result/island_source_assignments.csv.gz","assign")),compression="gzip",dtype=str)
    mb=(root/"fixed/source/source/gift_native_mainland_flora.csv.gz").read_bytes()
    m=pd.read_csv(io.BytesIO(mb),compression="gzip",dtype=str,usecols=["entity_ID","work_species"])
    geo=pd.read_csv(io.BytesIO(read_pin(geofile,"geography")),dtype={"island_id":str})
    assert f.island_id.nunique()==4505 and not f.duplicated(["island_id","accepted_species"]).any()
    assert m.entity_ID.nunique()==748 and geo.island_id.nunique()==8264 and geo.island_id.is_unique
    st=Stats(table(t));assert len(st.tab)==2180
    full=f.groupby("island_id").accepted_species.agg(set).to_dict()
    tar={i:st.ids(v) for i,v in full.items()}
    nat={i:st.ids(v) for i,v in f[f.origin_status.eq("native")].groupby("island_id").accepted_species.agg(set).items()}
    hi={i:frozenset(j for j in v if st.grade_high[j]) for i,v in tar.items()}
    regions={i:st.ids(v) for i,v in m.groupby("entity_ID").work_species.agg(set).items()}
    ass=ass[ass.source_mode.isin(MODES)].copy();ass.source_rank=pd.to_numeric(ass.source_rank,errors="raise")
    assert not ass.duplicated(["island_id","source_mode","source_rank"]).any()
    keys=ass.sort_values("source_rank").groupby(["island_id","source_mode"]).entity_ID.agg(tuple).to_dict()
    @functools.lru_cache(maxsize=512)
    def pool(k):return frozenset().union(*(regions[i] for i in k))
    plans={s:{v:{} for v in VARIANTS} for s in SCOPES};wide={s:{v:{} for v in VARIANTS} for s in SCOPES}
    for z in geo.itertuples():
        i=z.island_id
        if i not in full or not np.isfinite([z.island_latitude,z.island_longitude]).all():continue
        reg="tropical" if z.analysis_regime=="tropical" else "north" if z.island_latitude>23.5 else None
        if reg is None:continue
        block=f"{int(np.floor(z.island_latitude/10))}_{int(np.floor(z.island_longitude/10))}"
        for mode,nr in MODES.items():
            ids=keys.get((i,mode),())
            if len(ids)!=nr or any(j not in regions for j in ids):continue
            inc=pool(ids);dis=inc-tar[i]
            base={"region":reg,"geographic_block":block,"first_source_group":ids[0]}
            for scope,T,S,E in (("all_recorded",tar[i],inc,dis),("original_native",nat.get(i,frozenset()),inc,dis),
                 ("high_quality_only",hi[i],frozenset(j for j in inc if st.grade_high[j]),frozenset(j for j in dis if st.grade_high[j]))):
                if min(len(T),len(S),len(E))>=5:
                    for v,R in (("inclusive",S),("disjoint",E)):
                        wide[scope][mode,v][i]={"t":T,"s":R,**base}
                if not(st.eligible(T) and st.eligible(S) and st.eligible(E)):continue
                for v,R in (("inclusive",S),("disjoint",E)):
                    plans[scope][mode,v][i]={"t":T,"s":R,**base}
    return st,plans,wide,hashlib.sha256(mb).hexdigest()

def wide_joint(rows,st):
    if not rows:return {}
    A=membership([r["t"] for r in rows],len(st.tab));B=membership([r["s"] for r in rows],len(st.tab))
    T=A@st.V;S=B@st.V;T=T/T.sum(1)[:,None];S=S/S.sum(1)[:,None]
    Y=(T-S)[:,[3,7,11]]
    out={}
    for tag in GROUPINGS:
        d=paired_group_design(rows,tag)
        if d is not None:out[tag]=get_means(Y,d,np.ones(len(d[0]))).tolist()
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument("--original-root",type=pathlib.Path);p.add_argument("--geography",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    st,plans,wide,md=inputs(a.original_root,a.geography)
    panels={}
    for scope in SCOPES:
        common=set.intersection(*(set(plans[scope][v]) for v in VARIANTS))
        for mode,variant in VARIANTS:
            rr=plans[scope][mode,variant]
            z=[rr[i] for i in sorted(common)]
            orig=[rr[i] for i in sorted(rr)]
            y=wide[scope][mode,variant]
            key="__".join((scope,mode,variant))
            panel=analyze(z,st)
            panel["mode_specific_conditional_eligible"]=collections.Counter(r["region"] for r in orig)
            panel["wide_n_ge5"]=collections.Counter(r["region"] for r in y.values())
            panel["wide_raw_joint_delta"]=wide_joint(list(y.values()),st)
            panel["mode_specific_conditional_joint_delta"]=wide_joint(orig,st)
            panel["common_cohort_raw_joint_delta"]=wide_joint(z,st)
            panels[key]=panel
            print("FD_DEPTH_FACTOR_PANEL",key,json.dumps(panel,allow_nan=False,separators=(",",":")),flush=True)
    out={"botanical_parent":4505,"geographic_parent":8264,"mainland_candidate_regions":748,"complete_pair_species":2180,
         "mainland_sha256":md,"source_pins":PIN,"metric_order":NAMES,"region_order":OUT_NAMES,
         "panels":panels,"warning":"Post-result diagnostic, not independent replication or formal confirmation. All four donor variants use the intersection of the conditional eligible target islands within scope. Symmetric product decomposition is exact per-island and aggregated linearly. Shared group resampling across north and tropical; coherent genus weight intervals are sensitivity distributions, not calibrated sampling CIs. Joint differences on wider cohorts are descriptive selection diagnostics, not effects of eligibility. No actual animal guild, causal or evolutionary inference."}
    assert len(panels)==12
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    print("FD_DEPTH_FACTOR_COMPLETE",json.dumps({"panels":len(panels),"support_qualified":sum(v["status"]=="exploratory_estimable" for v in panels.values()),"metric_order":NAMES,"region_order":OUT_NAMES}),flush=True)
if __name__=="__main__":main()
