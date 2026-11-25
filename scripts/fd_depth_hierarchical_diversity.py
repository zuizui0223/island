#!/usr/bin/env python3
"""FD-owned hierarchical categorical diversity, fixed cohorts.

Alpha / within-group beta / between-group beta / gamma are probability
diversities of occurrence-profile distributions, using with-replacement
sampling. They are not the distinct-species Rao Q of the prior analysis.
Mainland comparisons are distributions of paired candidate pools, not
independent mainland assemblages or reconstructed colonist sources.
"""
import argparse, json, pathlib, itertools
import numpy as np
from fd_depth_factor_region_comparison import inputs, VARIANTS, SCOPES, REGIONS, SEED
from fd_depth_conditional_protected import membership, group_matrix, PIN

DRAWS=999
AXES=("joint","depth","symmetry")
PARTS=("alpha","beta_within_group","beta_between_group","gamma")
GROUPINGS=("geographic_block","first_source_group")
REGION_ORDER=("north","tropical","tropical_minus_north")
SLICES=(slice(0,12),slice(12,16),slice(16,19))
METRICS=[a+"__"+p for a in AXES for p in PARTS]

def expand(C):
    J=C.reshape(-1,3,4)
    return np.column_stack((C,J.sum(1),J.sum(2)))

def q(p):
    return np.stack([1-np.sum(p[...,sl]**2,axis=-1) for sl in SLICES],axis=-1)

def design(rows,tag):
    labels=sorted({r[tag] for r in rows})
    index={x:i for i,x in enumerate(labels)}
    regs=[]
    for reg in REGIONS:
        ix=np.array([i for i,r in enumerate(rows) if r["region"]==reg],int)
        if not len(ix):return None
        local=[rows[i][tag] for i in ix]
        keys=sorted(set(local))
        H=group_matrix(local)
        # group_matrix uses np.unique strings in sorted order.
        inv=np.array([keys.index(x) for x in local])
        n=np.bincount(inv)
        assert np.allclose(np.asarray(H.sum(1)).ravel(),1)
        regs.append((ix,H,np.array([index[x] for x in keys]),n))
    return labels,regs

def stats(P,var,design):
    # var = E squared norm of sample profile error for matched-size
    # independent SRSWOR at each island/source. Zero for observed P.
    out=[]
    for ix,H,ids,n in design[1]:
        means=H@P[ix]
        alpha=H@(q(P[ix])-var[ix])
        group_error=(H@var[ix])/n[:,None]
        group_gamma=q(means)-group_error
        within=group_gamma-alpha
        assert np.min(within)>-1e-9
        out.append((means,alpha,within,group_gamma,group_error,ids))
    return out

def average(st,w):
    vals=[]
    for means,alpha,within,gg,err,ids in st:
        z=w[...,ids]
        den=z.sum(-1,keepdims=True)
        if np.any(den<=0):raise ValueError("zero regional group weights")
        z=z/den
        a=z@alpha;b=z@within
        # A repeated bootstrap group is the same observed group. For
        # rarefaction expectation its coefficient is squared as a whole.
        gamma=q(z@means)-(z*z)@err
        c=gamma-z@gg
        assert np.min(c)>-1e-8
        assert np.allclose(a+b+c,gamma,atol=1e-9)
        vals.append(np.stack((a,b,c,gamma),axis=-1).reshape(*a.shape[:-1],12))
    return np.stack((vals[0],vals[1],vals[1]-vals[0]),axis=-2)

def interval(x):
    return np.moveaxis(np.quantile(x,[.025,.975],axis=0),0,-1).tolist() if len(x)>=.9*DRAWS else None

def rare_var(P,n,m):
    assert np.all(n>=m) and np.all(m>=2)
    return q(P)*((n-m)/(m*(n-1)))[:,None]

def analyze(rows,st):
    if not rows:return {"status":"not_estimable","n":0}
    A=membership([r["t"] for r in rows],len(st.tab))
    B=membership([r["s"] for r in rows],len(st.tab))
    V=expand(st.V)
    nt=np.asarray(A.sum(1)).ravel();ns=np.asarray(B.sum(1)).ravel()
    m=np.minimum(nt,ns);assert np.min(m)>=5
    def probs(w):
        sw=w[st.g];X=sw[:,None]*V
        return (A@X)/(A@sw)[:,None],(B@X)/(B@sw)[:,None]
    T,S=probs(np.ones(st.G))
    vT=rare_var(T,nt,m);vS=rare_var(S,ns,m)
    zero=np.zeros_like(vT)
    ds={k:d for k in GROUPINGS if (d:=design(rows,k)) is not None}
    if not ds:return {"status":"not_estimable","n":len(rows)}
    rng=np.random.default_rng(SEED)
    boot={};baseline={};rar={};raw_draw={};rare_draw={}
    for tag,d in ds.items():
        n=len(d[0]);w=np.ones(n)
        z=np.stack([np.bincount(rng.integers(n,size=n),minlength=n) for _ in range(DRAWS)])
        # Pair geographic/source bootstrap draws across both regions,
        # dropping a draw only if it loses an entire regional group set.
        keep=np.ones(len(z),bool)
        for ix,H,ids,nn in d[1]:keep&=z[:,ids].sum(1)>0
        z=z[keep];boot[tag]=z
        a,b=stats(T,zero,d),stats(S,zero,d)
        c,e=stats(T,vT,d),stats(S,vS,d)
        baseline[tag]=(average(a,w),average(b,w))
        rar[tag]=(average(c,w),average(e,w))
        raw_draw[tag]=average(a,z)-average(b,z)
        rare_draw[tag]=average(c,z)-average(e,z)
    genus={k:[] for k in ds};joint={k:[] for k in ds}
    for j in range(DRAWS):
        U,W=probs(np.maximum(rng.exponential(size=st.G),1e-8))
        for tag,d in ds.items():
            a,b=stats(U,zero,d),stats(W,zero,d)
            n=len(d[0]);ones=np.ones(n)
            genus[tag].append(average(a,ones)-average(b,ones))
            z=boot[tag]
            if len(z):joint[tag].append(average(a,z[j%len(z)])-average(b,z[j%len(z)]))
    out={}
    for tag,d in ds.items():
        support={reg:{"islands":len(ix),"groups":len(ids)} for reg,(ix,H,ids,n) in zip(REGIONS,d[1])}
        a,b=baseline[tag];c,e=rar[tag]
        out[tag]={"support":support,
            "observed":{"target":a.tolist(),"source":b.tolist(),"delta":(a-b).tolist(),"interval":interval(raw_draw[tag])},
            "matched_size_expected":{"target":c.tolist(),"source":e.tolist(),"delta":(c-e).tolist(),"interval":interval(rare_draw[tag])},
            "observed_genus_weight_interval":interval(genus[tag]),
            "observed_joint_genus_group_interval":interval(joint[tag]),
            "valid_draws":[len(boot[tag]),len(genus[tag]),len(joint[tag])]}
    gate=len(ds)==2 and all(v["islands"]>=200 and v["groups"]>=30 for a in out.values() for v in a["support"].values()) and all(min(a["valid_draws"])>=900 for a in out.values())
    return {"status":"exploratory_estimable" if gate else "insufficient_support","n":len(rows),"regions":{r:sum(x["region"]==r for x in rows) for r in REGIONS},
            "median_target_carriers":float(np.median(nt)),"median_source_carriers":float(np.median(ns)),
            "median_matched_sample_size":float(np.median(m)),"groupings":out}

def self_test():
    # Compare analytic rarefaction expectation with exhaustive finite
    # samples. Use two northern islands in one group, one in another,
    # and one tropical island, including a repeated group bootstrap.
    codes=[[0,3,7,8],[0,1,4],[3,7,11],[0,1,2]]
    rows=[{"region":"north","geographic_block":"a"},{"region":"north","geographic_block":"a"},
          {"region":"north","geographic_block":"b"},{"region":"tropical","geographic_block":"c"}]
    n=np.array([len(x) for x in codes],float);m=np.full(len(codes),2.)
    P=expand(np.array([np.bincount(x,minlength=12)/len(x) for x in codes]))
    d=design(rows,"geographic_block")
    for w in (np.ones(3),np.array([2.,1.,1.])):
        expected=average(stats(P,rare_var(P,n,m),d),w)
        choices=[list(itertools.combinations(x,2)) for x in codes]
        draws=[]
        for z in itertools.product(*choices):
            R=expand(np.array([np.bincount(x,minlength=12)/2 for x in z]))
            draws.append(average(stats(R,np.zeros((4,3)),d),w))
        assert np.allclose(expected,np.mean(draws,axis=0),atol=1e-12)
        assert np.allclose(expected.reshape(3,3,4)[...,:3].sum(-1),expected.reshape(3,3,4)[...,3])
    print("FD_HIERARCHICAL_Q_SYNTHETIC_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--original-root",type=pathlib.Path);p.add_argument("--geography",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    st,plans,wide,md=inputs(a.original_root,a.geography)
    panels={}
    for scope in SCOPES:
        for label,pop in (("wide_common_n5",wide),("conditional_common",plans)):
            ids=sorted(set.intersection(*(set(pop[scope][v]) for v in VARIANTS)))
            if scope=="all_recorded":assert len(ids)==(2244 if label=="wide_common_n5" else 465)
            for mode,variant in VARIANTS:
                rows=[pop[scope][mode,variant][i] for i in ids]
                key="__".join((scope,mode,variant,label))
                z=analyze(rows,st);panels[key]=z
                print("FD_HIERARCHICAL_Q_PANEL",key,json.dumps(z,allow_nan=False,separators=(",",":")),flush=True)
    assert len(panels)==24
    out={"botanical_parent":4505,"geographic_parent":8264,"mainland_candidate_regions":748,"complete_pair_species":2180,
         "mainland_sha256":md,"source_pins":PIN,"metric_order":METRICS,"region_order":REGION_ORDER,
         "panels":panels,"warning":"Post-result hierarchical diversity of categorical occurrence profiles, with replacement; distinct-species Rao Q is a different estimand. Equal groups, equal islands within groups. Paired source pools are assigned candidate pools, not independent mainland communities. Match sample size m=min(target,source) per island; exact expectation under independent within-roster sampling without replacement, not a missing-data correction. Rarefaction variance changes beta by sampling variation; neither equal-size nor genus sensitivity establishes ecological causality. Genus weights are sensitivity distributions, not calibrated CIs. All panels retained; no significance-selected endpoints."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    print("FD_HIERARCHICAL_Q_COMPLETE",json.dumps({"panels":len(panels),"support_qualified":sum(x["status"]=="exploratory_estimable" for x in panels.values()),"metric_order":METRICS,"region_order":REGION_ORDER}),flush=True)
if __name__=="__main__":main()
