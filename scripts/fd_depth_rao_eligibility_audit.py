#!/usr/bin/env python3
"""FD post-result eligibility audit. Exact nested-cohort bookkeeping.

Uses existing protected morphology rosters and no new endpoints, thresholds,
source or geography. Separate group support, carrier richness and depth
balance eligibility, not causal recording-bias correction.
"""
import argparse,collections,json,pathlib
import numpy as np
from fd_depth_factor_region_comparison import inputs,VARIANTS,SCOPES,REGIONS,SEED
from fd_depth_conditional_protected import membership,group_matrix,PIN,MIN_EACH
from fd_depth_rao_diversity import components,difference,NAMES as RAO_NAMES

DRAWS=1999
IX=[0,1,2,3,7,15,19,23,27]
METRICS=[RAO_NAMES[i] for i in IX]
STAGES=("wide_all_groups","wide_on_conditional_groups","richness_on_conditional_groups","conditional",
        "group_support_change","richness_selection_change","depth_balance_selection_change","total_change")
GROUPINGS=("geographic_block","first_source_group")
REGION_ORDER=("north","tropical","tropical_minus_north")

def design(rows,rich,cond,tag):
    keys=sorted({r[tag] for r in rows})
    lookup={k:i for i,k in enumerate(keys)}
    out=[]
    for reg in REGIONS:
        w=np.array([i for i,r in enumerate(rows) if r["region"]==reg],int)
        C=w[cond[w]]
        if not len(C):return None
        g={rows[i][tag] for i in C}
        wc=np.array([i for i in w if rows[i][tag] in g],int)
        R=wc[rich[wc]]
        assert set(C)<=set(R)<=set(wc)<=set(w)
        sizes=[]
        for ix in (w,wc,R,C):
            labels=[rows[i][tag] for i in ix]
            labs=sorted(set(labels))
            H=group_matrix(labels)
            sizes.append({"ix":ix,"H":H,"keys":np.array([lookup[k] for k in labs],int),
                          "islands":len(ix),"groups":len(labs)})
        assert all(np.array_equal(sizes[1]["keys"],z["keys"]) for z in sizes[2:])
        out.append(sizes)
    return keys,out

def means(Y,d,w):
    z=[]
    for reg in d[1]:
        a=[]
        for stage in reg:
            weights=w[stage["keys"]]
            if weights.sum()==0:return None
            a.append(np.average(stage["H"]@Y[stage["ix"]],axis=0,weights=weights))
        p=np.stack(a)
        dx=np.diff(p,axis=0)
        total=p[3]-p[0]
        assert np.allclose(dx.sum(0),total,atol=1e-12)
        z.append(np.concatenate((p,dx,total[None,:])))
    return np.stack([z[0],z[1],z[1]-z[0]])

def quant(draw):
    return np.moveaxis(np.quantile(draw,[.025,.975],axis=0),0,-1).tolist() if len(draw)>=.9*DRAWS else None

def analyze(rows,rich,cond,st):
    if not rows:return {"status":"not_estimable","n":0}
    assert len(rows)==len(rich)==len(cond) and np.all(~cond|rich)
    A=membership([r["t"] for r in rows],len(st.tab));B=membership([r["s"] for r in rows],len(st.tab))
    V=st.V;Dv=np.eye(4)[st.d]
    def calc(w):
        x=w[st.g]
        X=x[:,None]*V;Y=(x*x)[:,None]*Dv
        T=components(A@X,A@Y);S=components(B@X,B@Y)
        v=difference(T,S)[:,IX]
        # Separate observed target and donor for core Q to demonstrate
        # that the selection contrast need not be a change in target Q.
        return np.column_stack([v,T[:,0],S[:,0]])
    obs=calc(np.ones(st.G))
    ds={k:d for k in GROUPINGS if (d:=design(rows,rich,cond,k)) is not None}
    if not ds:return {"status":"not_estimable","n":len(rows)}
    points={k:means(obs,d,np.ones(len(d[0]))) for k,d in ds.items()}
    raw={k:[] for k in ds};tax={k:[] for k in ds};joint={k:[] for k in ds}
    rng=np.random.default_rng(SEED)
    for j in range(DRAWS):
        Y=calc(np.maximum(rng.exponential(size=st.G),1e-8))
        for k,d in ds.items():
            n=len(d[0]);w=np.bincount(rng.integers(n,size=n),minlength=n)
            for a,z in ((raw,means(obs,d,w)),(tax,means(Y,d,np.ones(n))),(joint,means(Y,d,w))):
                if z is not None:a[k].append(z)
    out={}
    for k,d in ds.items():
        support={reg:{name:{"islands":s["islands"],"groups":s["groups"]}
                       for name,s in zip(STAGES[:4],stages)} for reg,stages in zip(REGIONS,d[1])}
        assert np.allclose(points[k][...,0],points[k][...,-2]-points[k][...,-1],atol=1e-12)
        out[k]={"support":support,"point":points[k].tolist(),"interval":quant(raw[k]),
                "genus_weight_interval":quant(tax[k]),"joint_genus_group_interval":quant(joint[k]),
                "valid_draws":[len(raw[k]),len(tax[k]),len(joint[k])]}
    gate=len(ds)==2 and all(v["conditional"]["islands"]>=200 and v["conditional"]["groups"]>=30 for a in out.values() for v in a["support"].values()) and all(min(a["valid_draws"])>=1800 for a in out.values())
    return {"status":"exploratory_estimable" if gate else "insufficient_support","n":len(rows),"groupings":out}

def self_test():
    rows=[{"region":"north","geographic_block":"a"},{"region":"north","geographic_block":"b"},{"region":"north","geographic_block":"b"},
          {"region":"north","geographic_block":"b"},{"region":"tropical","geographic_block":"c"},{"region":"tropical","geographic_block":"c"},
          {"region":"tropical","geographic_block":"c"}]
    R=np.array([False,False,True,True,False,True,True])
    C=np.array([False,False,False,True,False,False,True])
    y=np.array([100,0,10,20,0,10,20],float)[:,None]
    d=design(rows,R,C,"geographic_block")
    a=means(y,d,np.ones(len(d[0])))
    assert np.allclose(a[0,:,0],[55,10,15,20,-45,5,5,-35])
    assert np.allclose(a[1,:,0],[10,10,15,20,0,5,5,10])
    assert np.allclose(a[2],a[1]-a[0])
    print("FD_RAO_ELIGIBILITY_SYNTHETIC_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--original-root",type=pathlib.Path);p.add_argument("--geography",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    st,plans,wide,md=inputs(a.original_root,a.geography)
    panels={}
    for scope in SCOPES:
        B=set.intersection(*(set(wide[scope][v]) for v in VARIANTS))
        C=set.intersection(*(set(plans[scope][v]) for v in VARIANTS))
        # Richness-only necessary condition for >=10 deep and >=10 nondeep
        # in both rosters and all four donor variants.
        R={i for i in B if all(min(len(wide[scope][v][i]["t"]),len(wide[scope][v][i]["s"]))>=2*MIN_EACH for v in VARIANTS)}
        assert C<=R<=B
        ids=sorted(B);rich=np.array([i in R for i in ids]);cond=np.array([i in C for i in ids])
        for mode,variant in VARIANTS:
            rows=[wide[scope][mode,variant][i] for i in ids]
            for i in C:
                assert wide[scope][mode,variant][i]==plans[scope][mode,variant][i]
            q=analyze(rows,rich,cond,st)
            q["global_nested_census"]={"wide":len(B),"richness":len(R),"conditional":len(C)}
            if scope=="all_recorded":
                assert len(B)==2244 and len(C)==465
                assert q["groupings"]["geographic_block"]["support"]["north"]["conditional"]["islands"]==264
                assert q["groupings"]["geographic_block"]["support"]["tropical"]["conditional"]["islands"]==201
            k="__".join((scope,mode,variant));panels[k]=q
            print("FD_RAO_ELIGIBILITY_PANEL",k,json.dumps(q,allow_nan=False,separators=(",",":")),flush=True)
    assert len(panels)==12
    out={"botanical_parent":4505,"geographic_parent":8264,"mainland_candidate_regions":748,"complete_pair_species":2180,
         "mainland_sha256":md,"source_pins":PIN,"metric_order":METRICS+["target_joint_Q","source_joint_Q"],
         "stage_order":STAGES,"region_order":REGION_ORDER,"panels":panels,
         "warning":"Post-result nested cohort selection accounting, not new ecological outcomes or a missing-data correction. Within each scope all four variants share wide, richness and conditional island sets. Group support then total-carrier richness then depth balance is a chosen order, not unique causal components. Joint Q is same observed roster for each island in all stages; coherent positive genus weights do not redefine inclusion. Geographic and first-source groups separately equal weighted; shared resampling across nested stages and regions. Genus-weight intervals are sensitivity distributions, not calibrated CIs. No pollinator mechanism."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    print("FD_RAO_ELIGIBILITY_COMPLETE",json.dumps({"panels":len(panels),"support_qualified":sum(v["status"]=="exploratory_estimable" for v in panels.values()),"metric_order":out["metric_order"],"stage_order":STAGES,"region_order":REGION_ORDER}),flush=True)
if __name__=="__main__":main()
