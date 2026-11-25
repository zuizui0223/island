#!/usr/bin/env python3
"""FD-owned post-result functional-diversity diagnostic, protected custodian.
Exact ordered-distinct-species pair probabilities, not KDE volume.
Reuse the frozen morphology table and the previous 465-island common cohort.
No protected identifiers or raw records in outputs.
"""
import argparse, collections, json, pathlib
import numpy as np
from fd_depth_conditional_protected import membership,group_matrix,PIN,SEED
from fd_depth_factor_region_comparison import inputs,paired_group_design,get_means,VARIANTS,SCOPES,REGIONS

DRAWS=1999
NAMES=["joint_Q","depth_Q","symmetry_Q","symmetry_difference_at_same_depth",
       *["same_depth_symmetry_difference__"+d for d in ("absent_or_open","shallow","intermediate","deep")],
       *["same_depth_pair_probability__"+d for d in ("absent_or_open","shallow","intermediate","deep")],
       *["common_depth_exposure_component__"+d for d in ("absent_or_open","shallow","intermediate","deep")],
       *["common_within_depth_diversity_component__"+d for d in ("absent_or_open","shallow","intermediate","deep")],
       *["target_only_pair_support__"+d for d in ("absent_or_open","shallow","intermediate","deep")],
       *["source_only_pair_support__"+d for d in ("absent_or_open","shallow","intermediate","deep")]]
GROUPINGS=("geographic_block","first_source_group")
REGION_ORDER=("north","tropical","tropical_minus_north")

def components(C,SS):
    # For general positive taxon weights, subtract self-pair mass.
    # Q is the probability of a difference when two distinct species
    # are sampled proportional to product weights.
    C=np.asarray(C,float).reshape(-1,3,4)
    SS=np.asarray(SS,float)
    total=C.sum((1,2));den=total**2-SS.sum(1)
    assert np.all(den>0)
    D=C.sum(1);S=C.sum(2)
    # ordered unlike-symmetry pairs within each depth class
    B=2*(C[:,0,:]*C[:,1,:]+C[:,0,:]*C[:,2,:]+C[:,1,:]*C[:,2,:])/den[:,None]
    d=2*sum(D[:,i]*D[:,j] for i in range(4) for j in range(i+1,4))/den
    sym=2*(S[:,0]*S[:,1]+S[:,0]*S[:,2]+S[:,1]*S[:,2])/den
    exposure=(D**2-SS)/den[:,None]
    assert np.all(exposure>=-1e-10)
    exposure=np.maximum(exposure,0)
    exposure[exposure<1e-14]=0
    assert np.allclose(d+exposure.sum(1),1,atol=1e-10)
    assert np.all(B<=exposure+1e-10)
    q=d+B.sum(1)
    return np.column_stack((q,d,sym,B.sum(1),B,exposure))

def difference(T,S):
    delta=T-S
    A=T[:,8:12];B=S[:,8:12]
    F=T[:,4:8];G=S[:,4:8]
    # Conditional Q at a depth is defined only if two distinct species
    # carry that depth. Retain one-sided support as its own component;
    # never substitute zero for an undefined conditional Q.
    ct=np.divide(F,A,out=np.zeros_like(F),where=A>0)
    cs=np.divide(G,B,out=np.zeros_like(G),where=B>0)
    common=(A>0)&(B>0)
    exposure=np.where(common,(A-B)*(ct+cs)/2,0)
    within=np.where(common,(A+B)*(ct-cs)/2,0)
    target=np.where((A>0)&(B==0),F,0)
    source=np.where((B>0)&(A==0),-G,0)
    assert np.allclose(delta[:,4:8],exposure+within+target+source,atol=1e-10)
    assert np.allclose(delta[:,0],delta[:,1]+delta[:,3],atol=1e-10)
    return np.column_stack((delta,exposure,within,target,source))

def group_design_all(rows,tag):
    return paired_group_design(rows,tag)

def quant(v):
    return np.quantile(v,[.025,.975],axis=0).transpose(1,2,0).tolist() if len(v)>=.9*DRAWS else None

def analyze(rows,st,weighted=True):
    if not rows:return {"status":"not_estimable","n":0}
    A=membership([r["t"] for r in rows],len(st.tab));B=membership([r["s"] for r in rows],len(st.tab))
    V=st.V
    depthV=np.eye(4)[st.d]
    def calc(w):
        sw=w[st.g]
        # Same-depth self-pair mass, not count squared.
        X=sw[:,None]*V;Y=(sw**2)[:,None]*depthV
        T=components(A@X,A@Y);S=components(B@X,B@Y)
        return difference(T,S),T,S
    obs,T,S=calc(np.ones(st.G))
    groups={}
    for tag in GROUPINGS:
        d=group_design_all(rows,tag)
        if d is not None:groups[tag]=d
    if not groups:return {"status":"not_estimable","n":len(rows)}
    rng=np.random.default_rng(SEED)
    points={k:get_means(obs,d,np.ones(len(d[0]))) for k,d in groups.items()}
    raw={k:[] for k in groups};tax={k:[] for k in groups};joint={k:[] for k in groups}
    for j in range(DRAWS):
        if weighted:
            z,_,_=calc(np.maximum(rng.exponential(size=st.G),1e-8))
        for k,d in groups.items():
            n=len(d[0]);w=np.bincount(rng.integers(n,size=n),minlength=n)
            v=get_means(obs,d,w)
            if v is not None:raw[k].append(v)
            if weighted:
                tax[k].append(get_means(z,d,np.ones(n)))
                v=get_means(z,d,w)
                if v is not None:joint[k].append(v)
    n_by={r:sum(x["region"]==r for x in rows) for r in REGIONS}
    out={}
    for k,d in groups.items():
        support={reg:{"islands":n_by[reg],"groups":len(ids)} for reg,ids in zip(REGIONS,d[1])}
        out[k]={"support":support,"point":points[k].tolist(),"interval":quant(raw[k]),
           "target_mean":get_means(T,d,np.ones(len(d[0]))).tolist(),
           "source_mean":get_means(S,d,np.ones(len(d[0]))).tolist(),
           "genus_weight_interval":quant(tax[k]) if weighted else None,
           "joint_genus_group_interval":quant(joint[k]) if weighted else None,
           "draws":[len(raw[k]),len(tax[k]),len(joint[k])]}
    gate=len(groups)==2 and all(v["islands"]>=200 and v["groups"]>=30 for a in out.values() for v in a["support"].values()) and all(a["draws"][0]>=1800 for a in out.values())
    return {"status":"exploratory_estimable" if gate else "insufficient_support","n":len(rows),"regions":n_by,"groupings":out,
            "target_depth_pair_support":[int(x) for x in (T[:,8:12]>0).sum(0)],"source_depth_pair_support":[int(x) for x in (S[:,8:12]>0).sum(0)]}

def self_test():
    from itertools import permutations
    from collections import Counter
    rng=np.random.default_rng(SEED)
    for _ in range(50):
        a=rng.integers(0,12,size=20);w=rng.uniform(.1,2,size=len(a))
        C=np.bincount(a,weights=w,minlength=12)
        SS=np.bincount(a%4,weights=w*w,minlength=4)
        x=components(C[None,:],SS[None,:])[0]
        pairs=list(permutations(range(len(a)),2))
        weights=np.array([w[i]*w[j] for i,j in pairs])
        den=weights.sum()
        vals=np.array([[a[i]!=a[j],a[i]%4!=a[j]%4,a[i]//4!=a[j]//4,a[i]%4==a[j]%4 and a[i]//4!=a[j]//4] for i,j in pairs])
        assert np.allclose(x[:4],weights@vals/den)
        b=rng.integers(0,12,size=15)
        Y=components(np.bincount(b,minlength=12)[None,:],np.bincount(b%4,minlength=4)[None,:])
        forward=difference(x[None,:],Y)
        reverse=difference(Y,x[None,:])
        # Joint/marginal contrasts, common exposure and common within
        # terms are individually antisymmetric. One-sided support terms
        # exchange roles under T <-> S, in addition to reversing sign.
        assert np.allclose(forward[:,:20],-reverse[:,:20])
        assert np.allclose(forward[:,20:24],-reverse[:,24:28])
        assert np.allclose(forward[:,24:28],-reverse[:,20:24])
        assert np.allclose(forward[:,20:24]+forward[:,24:28],-(reverse[:,20:24]+reverse[:,24:28]))
    # One-sided depth support is not assigned a fabricated Q.
    T=components(np.array([[1,0,0,1,0,0,0,1,0,0,0,0]]),np.array([[1,0,0,2]]))
    S=components(np.array([[2,0,0,0,0,0,0,1,0,0,0,0]]),np.array([[2,0,0,1]]))
    d=difference(T,S)[0];rev=difference(S,T)[0]
    assert d[23]>0 and np.isclose(d[19],0)
    assert np.isclose(d[23],-rev[27]) and np.isclose(d[27],0)
    assert np.allclose(d[:20],-rev[:20])
    assert np.allclose(d[20:24],-rev[24:28])
    assert np.allclose(d[24:28],-rev[20:24])
    print("FD_DEPTH_RAO_SYNTHETIC_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--original-root",type=pathlib.Path);p.add_argument("--geography",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    st,plans,wide,md=inputs(a.original_root,a.geography)
    panels={}
    for scope in SCOPES:
        common=set.intersection(*(set(plans[scope][v]) for v in VARIANTS))
        broad=set.intersection(*(set(wide[scope][v]) for v in VARIANTS))
        assert common<=broad
        for mode,variant in VARIANTS:
            for cohort,pop,ids in (("conditional_common",plans,common),("wide_common_n5",wide,broad)):
                rows=[pop[scope][mode,variant][i] for i in sorted(ids)]
                k="__".join((scope,mode,variant,cohort))
                q=analyze(rows,st)
                if scope=="all_recorded" and cohort=="conditional_common":
                    assert q["regions"]=={"north":264,"tropical":201}
                panels[k]=q
                print("FD_DEPTH_RAO_PANEL",k,json.dumps(q,allow_nan=False,separators=(",",":")),flush=True)
    out={"botanical_parent":4505,"geographic_parent":8264,"mainland_candidate_regions":748,"complete_pair_species":2180,
         "mainland_sha256":md,"source_pins":PIN,"metric_order":NAMES,"region_order":REGION_ORDER,"panels":panels,
         "warning":"Post-result two-axis categorical diversity diagnostic. Ordered distinct-species pair probability (weighted sensitivity proportional to product weights). Exact joint Q = depth Q + symmetry differences at equal depth. Conditional same-depth differences decomposed only where both rosters have pair support; one-sided terms separate. Geographic/source group bootstrap conditional on fixed rosters; genus-weight intervals not calibrated CIs. Nested cohort variants not independent. No causal pollinator, realized function, evolutionary or total trait-space volume claim."}
    assert len(panels)==24
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    print("FD_DEPTH_RAO_COMPLETE",json.dumps({"panels":len(panels),"support_qualified":sum(q["status"]=="exploratory_estimable" for q in panels.values()),"metric_order":NAMES,"region_order":REGION_ORDER}),flush=True)

if __name__=="__main__":main()
