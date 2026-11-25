#!/usr/bin/env python3
"""Protected execution of FD PR78's depth-conditional architecture contract.
No animal inputs or protected IDs exported. All scientific results exploratory.
"""
import argparse, collections, functools, hashlib, io, json, pathlib
import numpy as np
import pandas as pd
from scipy import sparse
from fd_original_tropical_pollinator_bridge import PIN, read_pin, accepted_genus
from fd_depth_conditional_architecture import SYMMETRY, DEPTH, CELLS, VECTORS, pair as exact_pair

SEED=20261125
MODES={"geo_k5":5,"geo50_climate10":10}
SCOPES=("all_recorded","original_native","all_on_native_eligible","high_quality_only")
REGIONS=("all","north","tropical","south")
MIN_EACH=10
DRAWS=1999
CONTRACT="https://github.com/zuizui0223/fd/blob/b3d6587b5996b76aedd341b6865699b4a363d0ff/shared_db_research/FD_SYNDROME_ARCHITECTURE_CONTRACT_20261125.md"
GROUP_NAMES=("geographic_block","first_source_group")

def sha(b):return hashlib.sha256(b).hexdigest()

def table(t):
    axes=("floral_symmetry","tube_depth_class")
    a=t[t.trait_name.isin(axes)&t.evidence_scope.eq("species_direct")&
        t.evidence_quality.isin(["high","medium"])&t.normalized_value.notna()].copy()
    assert a.accepted_species.notna().all() and not a.duplicated(["accepted_species","trait_name"]).any()
    a["normalized_value"]=a.normalized_value.str.strip()
    a=a[(a.trait_name.eq(axes[0])&a.normalized_value.isin(SYMMETRY))|
        (a.trait_name.eq(axes[1])&a.normalized_value.isin(DEPTH))]
    v=a.pivot(index="accepted_species",columns="trait_name",values="normalized_value").reindex(columns=axes).dropna()
    q=a.pivot(index="accepted_species",columns="trait_name",values="evidence_quality").loc[v.index,list(axes)]
    assert q.notna().all().all()
    v["grade_pattern"]=q.astype(str).agg("|".join,axis=1)
    v["genus"]=[accepted_genus(s) for s in v.index]
    return v[v.genus.notna()].sort_index()

class Stats:
    def __init__(self,tab):
        self.tab=tab
        self.names={s:i for i,s in enumerate(tab.index)}
        self.s=np.array([SYMMETRY.index(x) for x in tab.floral_symmetry])
        self.d=np.array([DEPTH.index(x) for x in tab.tube_depth_class])
        self.g=pd.factorize(tab.genus,sort=True)[0]
        self.sg=pd.factorize(pd.MultiIndex.from_frame(tab[["genus","grade_pattern"]]),sort=True)[0]
        self.K=int(self.sg.max())+1
        self.G=int(self.g.max())+1
        self.c=self.s*4+self.d
        self.V=np.eye(12)[self.c]
        self.deep=(self.d==3)
        self.grade_high=tab.grade_pattern.eq("high|high").to_numpy()
    def ids(self,names):return frozenset(self.names[s] for s in names if s in self.names)
    def eligible(self,ids):
        k=sum(self.deep[i] for i in ids);return k>=MIN_EACH and len(ids)-k>=MIN_EACH
    def stat(self,ids):
        ix=np.asarray(sorted(ids),int);n=len(ix)
        assert n>0
        count=np.bincount(self.sg[ix]*12+self.c[ix],minlength=self.K*12).reshape(self.K,3,4)
        total=count.sum(0)
        p=total/n
        h=count.sum((1,2))
        use=h>0
        expected=np.einsum("gi,gj,g->ij",count[use].sum(2),count[use].sum(1),1/h[use])/n
        marg=p.sum(1)[:,None]*p.sum(0)[None,:]
        return {"count":count,"joint":p,"independent_joint":marg,"between":expected-marg,"within":p-expected,
                "n":n,"genera":len(set(self.g[ix])),"high_share":float(self.grade_high[ix].mean())}

def cond_parts(t,s,deep):
    A=t["count"][:,:,3] if deep else t["count"][:,:,:3].sum(2)
    B=s["count"][:,:,3] if deep else s["count"][:,:,:3].sum(2)
    n=A.sum(1);m=B.sum(1)
    assert min(n.sum(),m.sum())>=MIN_EACH
    w=n/n.sum();v=m/m.sum()
    qa=np.divide(A,n[:,None],out=np.zeros_like(A,dtype=float),where=n[:,None]>0)
    qb=np.divide(B,m[:,None],out=np.zeros_like(B,dtype=float),where=m[:,None]>0)
    both=(n>0)&(m>0);ot=(n>0)&(m==0);os=(m>0)&(n==0)
    tq=w@qa;sq=v@qb;diff=tq-sq
    part=np.stack([diff,((w[both]+v[both])/2)@(qa[both]-qb[both]),
        ((w[both]-v[both])/2)@(qa[both]+qb[both]),w[ot]@qa[ot],-v[os]@qb[os]])
    assert np.allclose(part[0],part[1:].sum(0),atol=1e-12)
    assert np.isclose(tq.sum(),1) and np.isclose(sq.sum(),1)
    diag={"common_strata":int(both.sum()),"target_only_strata":int(ot.sum()),"source_only_strata":int(os.sum()),
          "common_target_mass":float(w[both].sum()),"common_source_mass":float(v[both].sum()),
          "common_variable_strata":int((((A+B)>0).sum(1)>1)[both].sum()),
          "common_changed_distribution_strata":int((np.max(np.abs(qa-qb),axis=1)>1e-12)[both].sum()),
          "target_carriers":int(n.sum()),"source_carriers":int(m.sum())}
    return tq,sq,part,diag

def compare(t,s):
    names=[];values=[]
    def add(prefix,a,labels):
        names.extend(prefix+"__"+x for x in labels)
        values.extend(np.asarray(a,float).ravel())
    cellnames=[a+"__"+b for a,b in CELLS]
    add("target_joint",t["joint"],cellnames);add("source_joint",s["joint"],cellnames)
    for k in ("joint","independent_joint","between","within"):
        add("delta_"+k,t[k]-s[k],cellnames)
    assert np.allclose(t["joint"]-s["joint"],(t["independent_joint"]-s["independent_joint"])+(t["between"]-s["between"])+(t["within"]-s["within"]),atol=1e-12)
    for k,a in (("symmetry",0),("depth",1)):
        axis=1 if a==0 else 0;labels=SYMMETRY if a==0 else DEPTH
        add("target_"+k,t["joint"].sum(axis),labels);add("source_"+k,s["joint"].sum(axis),labels)
        add("delta_"+k,(t["joint"]-s["joint"]).sum(axis),labels)
    diagnostics={};parts=[]
    for tag,deep in (("deep",True),("nondeep",False)):
        tq,sq,p,diag=cond_parts(t,s,deep)
        add(tag+"__target_conditional",tq,SYMMETRY);add(tag+"__source_conditional",sq,SYMMETRY)
        for j,k in enumerate(VECTORS):add(tag+"__"+k,p[j],SYMMETRY)
        parts.append(p);diagnostics[tag]=diag
    for j,k in enumerate(VECTORS):add("deep_minus_nondeep__"+k,parts[0][j]-parts[1][j],SYMMETRY)
    add("target_grade",t["high_share"],["high_high_share"])
    add("source_grade",s["high_share"],["high_high_share"])
    return names,np.asarray(values,float),diagnostics

def group_matrix(labels):
    _,g=np.unique(np.asarray(labels,str),return_inverse=True);n=np.bincount(g)
    return sparse.csr_matrix((1/n[g],(g,np.arange(len(g)))),shape=(len(n),len(g)))

def summarize(rows,labels):
    if not rows:return {"islands":0,"groups":0,"status":"not_estimable"}
    x=np.stack([r["v"] for r in rows]);G=group_matrix(labels);z=G@x;n=len(z)
    mu=z.mean(0);rng=np.random.default_rng(SEED)
    boot=np.stack([z[rng.integers(n,size=n)].mean(0) for j in range(DRAWS)]) if n>=2 else None
    lohi=np.quantile(boot,[.025,.975],axis=0).T.tolist() if boot is not None else None
    return {"islands":len(rows),"groups":n,"status":"estimated" if n>=2 else "descriptive_one_group",
        "mean":mu.tolist(),"interval":lohi,"draws":DRAWS if boot is not None else 0}

def membership(rosters,p):
    r=[];c=[]
    for j,ids in enumerate(rosters):
        r.extend([j]*len(ids));c.extend(sorted(ids))
    return sparse.csr_matrix((np.ones(len(r)),(r,c)),shape=(len(rosters),p))

def conditional6(C):
    D=C.reshape(-1,3,4)
    x=D[:,:,3];y=D[:,:,:3].sum(2)
    if np.any(x.sum(1)<=0) or np.any(y.sum(1)<=0):raise ValueError("zero weighted conditional denominator")
    x=x/x.sum(1)[:,None];y=y/y.sum(1)[:,None]
    return np.column_stack([x,y,x-y])

def genus_draws(rows,st):
    # Shared genus weights over all focal and mainland names; no raw IDs
    # or individual-level values leave the custodian runner.
    if not rows:return {}
    A=membership([r["t"] for r in rows],len(st.tab));B=membership([r["s"] for r in rows],len(st.tab))
    def effect(w):
        return conditional6(A@(w[st.g,None]*st.V))-conditional6(B@(w[st.g,None]*st.V))
    base=effect(np.ones(st.G));assert np.allclose(base.sum(1),0,atol=1e-12)
    sub={}
    for reg in REGIONS:
        ix=np.array([i for i,r in enumerate(rows) if reg=="all" or r["region"]==reg])
        if not len(ix):continue
        for label in GROUP_NAMES:
            H=group_matrix([rows[i][label] for i in ix]);sub[(reg,label)]=(ix,H)
    rng=np.random.default_rng(SEED);draw={k:[] for k in sub};tax={k:[] for k in sub}
    for j in range(DRAWS):
        w=np.maximum(rng.exponential(size=st.G),1e-12);d=effect(w)
        for k,(ix,H) in sub.items():
            v=H@d[ix];n=H.shape[0]
            tax[k].append(v.mean(0));draw[k].append(v[rng.integers(n,size=n)].mean(0))
    out={}
    # Exact first-order derivative of conditional probabilities at unit
    # genus weights, used ONLY to rank genera for deletion sensitivity.
    # Then recompute ratios with actual deleted genera; report attrition.
    for reg in REGIONS:
        z={}
        for label in GROUP_NAMES:
            k=(reg,label)
            if k not in sub:continue
            ix,H=sub[k];b=(H@base[ix]).mean(0)
            # Analytic derivative of mean deep radial contrast for each
            # genus; equal geographic/source-group and within-group weights.
            w=np.asarray(H.sum(0)).ravel()/H.shape[0]
            infl=np.zeros(st.G)
            for j in (0,1):
                R=(A if j==0 else B)[ix].tocsr()
                deepn=np.asarray(R@st.deep.astype(float)).ravel()
                radial=(st.deep&(st.s==0)).astype(float)
                p=np.asarray(R@radial).ravel()/deepn
                v0=np.asarray(R.T@(w/deepn)).ravel()*radial
                v1=np.asarray(R.T@(w*p/deepn)).ravel()*st.deep
                infl+=(1 if j==0 else -1)*np.bincount(st.g,weights=v0-v1,minlength=st.G)
            assert abs(infl.sum())<1e-10
            active=np.flatnonzero(np.abs(infl)>1e-14)
            ranked=active[np.argsort(-np.abs(infl[active]),kind="stable")]
            total=float(np.abs(infl).sum())
            deletions={}
            for m in (1,3,5):
                drop=ranked[:m];ww=np.ones(st.G);ww[drop]=0
                T=A[ix]@(ww[st.g,None]*st.V);S=B[ix]@(ww[st.g,None]*st.V)
                def ok(C):
                    D=C.reshape(-1,3,4)
                    return (D[:,:,3].sum(1)>=MIN_EACH)&(D[:,:,:3].sum((1,2))>=MIN_EACH)
                keep=ok(T)&ok(S)
                if keep.any():
                    GG=group_matrix([rows[i][label] for i in ix[keep]])
                    de=conditional6(T[keep])-conditional6(S[keep])
                    val=(GG@de).mean(0).tolist()
                    orig=(GG@base[ix[keep]]).mean(0).tolist()
                else:val=orig=None
                deletions[str(m)]={"removed_genera":len(drop),"retained_islands":int(keep.sum()),"recomputed_mean":val,"original_mean_on_retained":orig}
            z[label]={"islands":len(ix),"groups":H.shape[0],"draws":DRAWS,"mean":b.tolist(),
                "genus_weight_interval":np.quantile(tax[k],[.025,.975],axis=0).T.tolist(),
                "joint_genus_group_interval":np.quantile(draw[k],[.025,.975],axis=0).T.tolist(),
                "radial_deep_influence_nonzero_genera":len(active),
                "top1_abs_derivative_share":float(abs(infl[ranked[0]])/total) if total else None,
                "top5_abs_derivative_share":float(np.abs(infl[ranked[:5]]).sum()/total) if total else None,
                "radial_deep_ranked_genus_deletion":deletions}
        if z:out[reg]=z
    return out

def self_test():
    rng=np.random.default_rng(SEED)
    rows=[]
    for j in range(180):
        rows.append({"floral_symmetry":SYMMETRY[rng.integers(3)],"tube_depth_class":DEPTH[rng.integers(4)],
                     "genus":"G"+str(rng.integers(7)),"grade_pattern":"high|high" if j%2 else "medium|high"})
    df=pd.DataFrame(rows,index=["P "+str(i) for i in range(len(rows))]);st=Stats(df)
    t=frozenset(range(110));s=frozenset(range(70,180))
    assert st.eligible(t) and st.eligible(s)
    a,b=st.stat(t),st.stat(s);nm,v,d=compare(a,b)
    def rr(ids):return [tuple(df.iloc[i][["floral_symmetry","tube_depth_class","genus","grade_pattern"]]) for i in sorted(ids)]
    e=exact_pair(rr(t),rr(s));assert e["eligible"]
    z=dict(zip(nm,v))
    for depth in ("deep","nondeep"):
        for k in VECTORS:
            for j,sym in enumerate(SYMMETRY):assert np.isclose(z[depth+"__"+k+"__"+sym],float(e["conditional"][depth][k][j]))
    for k in ("joint","independent_joint","between","within"):
        assert np.allclose([z["delta_"+k+"__"+a+"__"+b] for a,b in CELLS],np.array(e["delta"][k],float))
    global DRAWS
    old=DRAWS;DRAWS=19
    try:
        q=genus_draws([{"t":t,"s":t,"region":"tropical","geographic_block":"a","first_source_group":"b"}],st)
        assert np.allclose(q["tropical"]["geographic_block"]["joint_genus_group_interval"],0)
    finally:DRAWS=old
    print("FD_DEPTH_PROTECTED_SYNTHETIC_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--original-root",type=pathlib.Path);p.add_argument("--geography",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    r=a.original_root
    f=pd.read_csv(io.BytesIO(read_pin(r/"fixed/canonical/input/chapter1_status_flora.csv.gz","flora")),compression="gzip",dtype=str,usecols=["island_id","accepted_species","origin_status"])
    t=pd.read_csv(io.BytesIO(read_pin(r/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz","traits")),compression="gzip",dtype=str)
    ass=pd.read_csv(io.BytesIO(read_pin(r/"fixed/source/result/island_source_assignments.csv.gz","assign")),compression="gzip",dtype=str)
    mb=(r/"fixed/source/source/gift_native_mainland_flora.csv.gz").read_bytes()
    m=pd.read_csv(io.BytesIO(mb),compression="gzip",dtype=str,usecols=["entity_ID","work_species"])
    geo=pd.read_csv(io.BytesIO(read_pin(a.geography,"geography")),dtype={"island_id":str})
    assert f.island_id.nunique()==4505 and not f.duplicated(["island_id","accepted_species"]).any()
    assert m.entity_ID.nunique()==748 and geo.island_id.nunique()==8264 and geo.island_id.is_unique
    assert f.origin_status.isin(["native","introduced","unresolved"]).all()
    st=Stats(table(t))
    full=f.groupby("island_id").accepted_species.agg(set).to_dict()
    tar={i:st.ids(v) for i,v in full.items()}
    nat={i:st.ids(v) for i,v in f[f.origin_status.eq("native")].groupby("island_id").accepted_species.agg(set).items()}
    high={i:frozenset(j for j in v if st.grade_high[j]) for i,v in tar.items()}
    regions={i:st.ids(v) for i,v in m.groupby("entity_ID").work_species.agg(set).items()}
    aa=ass[ass.source_mode.isin(MODES)].copy();aa.source_rank=pd.to_numeric(aa.source_rank,errors="raise")
    assert not aa.duplicated(["island_id","source_mode","source_rank"]).any()
    keys=aa.sort_values("source_rank").groupby(["island_id","source_mode"]).entity_ID.agg(tuple).to_dict()
    @functools.lru_cache(maxsize=512)
    def pool(k):return frozenset().union(*(regions[i] for i in k))
    geoc={}
    for z in geo.itertuples():
        i=z.island_id
        if i not in full or not np.isfinite([z.island_latitude,z.island_longitude]).all():continue
        reg="tropical" if z.analysis_regime=="tropical" else "north" if z.island_latitude>23.5 else "south" if z.island_latitude< -23.5 else None
        if reg is not None:geoc[i]=(reg,f"{int(np.floor(z.island_latitude/10))}_{int(np.floor(z.island_longitude/10))}")
    panels={};metric_order=None
    for mode,nr in MODES.items():
        plans={k:[] for k in SCOPES};census={k:collections.Counter() for k in SCOPES}
        for i,(reg,block) in geoc.items():
            ids=keys.get((i,mode),())
            if len(ids)!=nr or any(k not in regions for k in ids):continue
            inc=pool(ids);dis=inc-tar[i] # all focal island names, restricted to known pair carriers
            assert not dis&tar[i]
            base={"region":reg,"geographic_block":block,"first_source_group":ids[0]}
            for scope,T,S,E in (("all_recorded",tar[i],inc,dis),("original_native",nat.get(i,frozenset()),inc,dis),
                 ("high_quality_only",high[i],frozenset(j for j in inc if st.grade_high[j]),frozenset(j for j in dis if st.grade_high[j]))):
                c=census[scope];c["geography_and_source_available"]+=1
                if len(T)>=MIN_EACH and min(len(S),len(E))>=MIN_EACH:c["joint_ge10"]+=1
                if not st.eligible(T):continue
                c["target_deep_and_nondeep_ge10"]+=1
                if not(st.eligible(S) and st.eligible(E)):continue
                c["both_source_variants_conditional_eligible"]+=1
                plans[scope].append((i,T,S,E,base))
                if scope=="original_native":
                    assert st.eligible(tar[i]);plans["all_on_native_eligible"].append((i,tar[i],S,E,base))
        census["all_on_native_eligible"]=census["original_native"].copy()
        for scope in SCOPES:
            rows={"inclusive":[],"disjoint":[]}
            # Source moments can be reused across islands.
            @functools.lru_cache(maxsize=128)
            def calc(ids):return st.stat(ids)
            for i,T,S,E,meta in plans[scope]:
                u=calc(T)
                for donor,R in (("inclusive",S),("disjoint",E)):
                    names,v,diag=compare(u,calc(R))
                    if metric_order is None:metric_order=names
                    assert metric_order==names
                    rows[donor].append({"t":T,"s":R,"v":v,"diag":diag,**meta})
            calc.cache_clear()
            assert len(rows["inclusive"])==len(rows["disjoint"])
            for donor,rr in rows.items():
                G=genus_draws(rr,st)
                for reg in REGIONS:
                    z=[v for v in rr if reg=="all" or v["region"]==reg]
                    group={label:summarize(z,[r[label] for r in z]) for label in GROUP_NAMES}
                    support=len(z)>=200 and all(v["groups"]>=30 and v.get("draws",0)>=1800 for v in group.values())
                    key="__".join((scope,mode,donor,reg))
                    q={"status":"exploratory_estimable" if support else "insufficient_support" if z else "not_estimable",
                       "islands":len(z),"groupings":group,"scope_support_census":dict(census[scope]),
                       "genus_sensitivity":G.get(reg,{}),"distinct_target_taxa":len(set().union(*(v["t"] for v in z))),
                       "distinct_source_taxa":len(set().union(*(v["s"] for v in z))),
                       "conditional_stratum_diagnostics":{d:{k:float(np.median([v["diag"][d][k] for v in z])) for k in z[0]["diag"][d]} for d in ("deep","nondeep")} if z else {}}
                    panels[key]=q
                    print("FD_DEPTH_CONDITIONAL_PANEL",key,json.dumps(q,allow_nan=False,separators=(",",":")),flush=True)
        pool.cache_clear()
    assert len(panels)==64
    out={"contract":CONTRACT,"botanical_parent":4505,"geographic_parent":8264,"mainland_candidate_regions":748,
         "source_pins":PIN,"mainland_sha256":sha(mb),"complete_pair_species":len(st.tab),
         "metric_order":metric_order,"genus_metric_order":[d+"__"+s for d in ("deep","nondeep","deep_minus_nondeep") for s in SYMMETRY],
         "panels":panels,"warning":"Post-pilot source-relative descriptive architecture, exact paired Direct states. Fixed 10 deep and 10 nondeep target and both donor variants; condition on native/high availability. Genus weights are sensitivity distributions not calibrated CIs; derivative-ranked deletion changes eligible cohort, original on retained also reported. No effective pollinator, adaptation, evolutionary or causal claims."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    print("FD_DEPTH_CONDITIONAL_COMPLETE",json.dumps({"panels":len(panels),"support_qualified":sum(v["status"]=="exploratory_estimable" for v in panels.values()),"complete_pair_species":len(st.tab),"metric_order":metric_order,"genus_metric_order":out["genus_metric_order"]}),flush=True)

if __name__=="__main__":main()
