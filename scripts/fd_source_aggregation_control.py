#!/usr/bin/env python3
"""FD post-hierarchy diagnostic: first-ranked vs pooled mainland reference.
Use exactly the existing wide-common paired carriers, nested support, same
target and same matched sample size for every compared reference.
No raw records or identifiers exported. Not a new independent ecological test.
"""
import argparse, collections, hashlib, io, json, pathlib
import numpy as np
import pandas as pd
from fd_depth_factor_region_comparison import inputs, VARIANTS, SCOPES, REGIONS, SEED
from fd_depth_conditional_protected import membership, read_pin, PIN, MODES
from fd_depth_hierarchical_diversity import expand, q, design, stats, average, rare_var, METRICS, GROUPINGS, REGION_ORDER

DRAWS=999
CONTRASTS=("pooled_full","pooled_restricted","rank1_restricted","eligibility_change","reference_change","total_change","pooled_matched_size","rank1_matched_size","matched_reference_change")
MIN_N=5

def quant(x):
    return np.moveaxis(np.quantile(x,[.025,.975],axis=0),0,-1).tolist() if len(x)>=.9*DRAWS else None

def point(a,b,c,d):
    # a,b,c are target-source contrasts for wide union, restricted union,
    # and restricted rank1; d gives matched-size pooled/rank1 contrasts.
    return np.stack((a,b,c,b-a,c-b,c-a,d[0],d[1],d[1]-d[0]),axis=-2)

def make_probs(rows,st):
    V=expand(st.V)
    mats={k:membership([r[k] for r in rows],len(st.tab)) for k in ("t","s")}
    n={k:np.asarray(M.sum(1)).ravel() for k,M in mats.items()}
    assert min(n["t"])>=MIN_N and min(n["s"])>=MIN_N
    def prob(w):
        sw=w[st.g];X=sw[:,None]*V
        return tuple((M@X)/(M@sw)[:,None] for M in mats.values())
    return prob,n

def analyze(wide,restricted,rank1,ms,st):
    assert len(restricted)==len(rank1)==len(ms)
    if not restricted:return {"status":"not_estimable","n_wide":len(wide),"n_restricted":0}
    assert all(r["t"]==s["t"] and r["region"]==s["region"] and
        r["first_source_group"]==s["first_source_group"] and s["s"]<=r["s"]
        for r,s in zip(restricted,rank1))
    W,nw=make_probs(wide,st);U,nu=make_probs(restricted,st);R,nr=make_probs(rank1,st)
    ms=np.asarray(ms,float)
    assert np.min(ms)>=MIN_N and all(np.all(n>=ms) for n in (*nu.values(),*nr.values()))
    one=np.ones(st.G)
    w,u,r=W(one),U(one),R(one)
    assert np.allclose(u[0],r[0])
    zeros=lambda x:np.zeros((len(x),3))
    vu=[rare_var(p,n,ms) for p,n in zip(u,nu.values())]
    vr=[rare_var(p,n,ms) for p,n in zip(r,nr.values())]
    assert np.allclose(vu[0],vr[0])
    ds={}
    for tag in GROUPINGS:
        a=design(wide,tag);b=design(restricted,tag)
        if a is None or b is None:continue
        labels=a[0];ix={x:i for i,x in enumerate(labels)}
        # Use the same bootstrap weight for a group whenever it appears
        # in either cohort. First-source labels are unchanged by pool choice.
        ind=np.array([ix[x] for x in b[0]],int)
        ds[tag]=(a,b,ind)
    if not ds:return {"status":"not_estimable","n_wide":len(wide),"n_restricted":len(restricted)}
    def prepared(a,b,pw,pu,pr):
        return (tuple(stats(x,zeros(x),a) for x in pw),
                tuple(stats(x,zeros(x),b) for x in pu),
                tuple(stats(x,zeros(x),b) for x in pr))
    def deltas(sts,v):
        return average(sts[0],v)-average(sts[1],v)
    def obs_values(p,weights,subweights):
        return [deltas(x,z) for x,z in zip(p,(weights,subweights,subweights))]
    rng=np.random.default_rng(SEED)
    configs={};out={}
    raw={};matched={};genus={};joint={}
    for tag,(a,b,ix) in ds.items():
        n=len(a[0]);ow=np.ones(n);sub=ow[ix]
        boot=np.stack([np.bincount(rng.integers(n,size=n),minlength=n) for _ in range(DRAWS)])
        keep=np.ones(DRAWS,bool)
        for dd,z in ((a,boot),(b,boot[:,ix])):
            for rr in dd[1]:keep&=z[:,rr[2]].sum(1)>0
        boot=boot[keep];bb=boot[:,ix]
        p=prepared(a,b,w,u,r)
        rm=(tuple(stats(x,v,b) for x,v in zip(u,vu)),tuple(stats(x,v,b) for x,v in zip(r,vr)))
        base=obs_values(p,ow,sub)
        match=[deltas(z,sub) for z in rm]
        values=point(*base,match)
        assert np.allclose(values[...,5,:],values[...,3,:]+values[...,4,:])
        # Regional and cohort differences computed on identical draws.
        raw[tag]=point(*obs_values(p,boot,bb),[deltas(z,bb) for z in rm])
        configs[tag]=(a,b,ix,boot)
        genus[tag]=[];joint[tag]=[]
        support={reg:{"wide_islands":len(aa[0]),"wide_groups":len(aa[2]),
                      "restricted_islands":len(cc[0]),"restricted_groups":len(cc[2])}
                 for reg,aa,cc in zip(REGIONS,a[1],b[1])}
        out[tag]={"support":support,"point":values.tolist(),"geographic_group_interval":quant(raw[tag]),"valid_bootstrap":len(boot)}
    # Genus sensitivity uses actual weighted roster probabilities only.
    # Do not combine the analytic unweighted SRSWOR correction with weights.
    for j in range(DRAWS):
        gw=np.maximum(rng.exponential(size=st.G),1e-8);pw,pu,pr=W(gw),U(gw),R(gw)
        for tag,(a,b,ix,boot) in configs.items():
            pp=prepared(a,b,pw,pu,pr)
            def six(z):
                x,y,v=obs_values(pp,z,z[ix])
                return np.stack((x,y,v,y-x,v-y,v-x),axis=-2)
            genus[tag].append(six(np.ones(len(a[0]))))
            if len(boot):joint[tag].append(six(boot[j%len(boot)]))
    for tag,x in out.items():
        x["observed_genus_weight_interval"]=quant(genus[tag])
        x["observed_joint_genus_group_interval"]=quant(joint[tag])
        x["valid_genus_draws"]=len(genus[tag]);x["valid_joint_draws"]=len(joint[tag])
    gate=len(out)==2 and all(v["restricted_islands"]>=200 and v["restricted_groups"]>=30
        for a in out.values() for v in a["support"].values()) and all(min(a["valid_bootstrap"],a["valid_genus_draws"],a["valid_joint_draws"])>=900 for a in out.values())
    return {"status":"exploratory_estimable" if gate else "insufficient_support",
            "n_wide":len(wide),"n_restricted":len(restricted),"median_common_sample_size":float(np.median(ms)),
            "median_pooled_carriers":float(np.median(nu["s"])),"median_rank1_carriers":float(np.median(nr["s"])),
            "groupings":out}

def self_test():
    a=np.arange(36).reshape(3,12);b=a+1;c=a+3
    z=point(a,b,c,[b+10,c+11])
    assert z.shape==(3,9,12)
    assert np.all(z[:,3,:]==1) and np.all(z[:,4,:]==2) and np.all(z[:,5,:]==3)
    assert np.all(z[:,8,:]==3)
    # Identity when pooled reference equals rank1 reference.
    z=point(a,b,b,[b,b]);assert not z[:,4,:].any() and not z[:,8,:].any()
    print("FD_SOURCE_AGGREGATION_CONTROL_SYNTHETIC_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--original-root",type=pathlib.Path);p.add_argument("--geography",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    st,plans,wide,md=inputs(a.original_root,a.geography)
    ass=pd.read_csv(io.BytesIO(read_pin(a.original_root/"fixed/source/result/island_source_assignments.csv.gz","assign")),compression="gzip",dtype=str)
    raw=(a.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz").read_bytes();assert hashlib.sha256(raw).hexdigest()==md
    m=pd.read_csv(io.BytesIO(raw),compression="gzip",dtype=str,usecols=["entity_ID","work_species"])
    pools={k:st.ids(v) for k,v in m.groupby("entity_ID").work_species.agg(set).items()}
    ass=ass[ass.source_mode.isin(MODES)].copy();ass.source_rank=pd.to_numeric(ass.source_rank,errors="raise")
    keys=ass.sort_values("source_rank").groupby(["island_id","source_mode"]).entity_ID.first().to_dict()
    panels={}
    for scope in SCOPES:
        B=set.intersection(*(set(wide[scope][v]) for v in VARIANTS))
        rank={v:{} for v in VARIANTS}
        counts={}
        for i in sorted(B):
            sizes=[]
            for mode,variant in VARIANTS:
                v=(mode,variant);r=wide[scope][v][i]
                first=keys[i,mode];assert first==r["first_source_group"]
                names=pools[first]
                # Native scope restricts only target. Disjoint means removal
                # of all focal pair-recorded species, not just native/high.
                if variant=="disjoint":names=names-wide["all_recorded"][v][i]["t"]
                if scope=="high_quality_only":names=frozenset(j for j in names if st.grade_high[j])
                assert names<=r["s"]
                rank[v][i]={**r,"s":names}
                sizes.extend((len(r["t"]),len(r["s"]),len(names)))
            counts[i]=min(sizes)
        C=sorted(i for i in B if counts[i]>=MIN_N)
        if scope=="all_recorded":assert len(B)==2244
        for mode,variant in VARIANTS:
            v=(mode,variant)
            z=analyze([wide[scope][v][i] for i in sorted(B)],[wide[scope][v][i] for i in C],
                      [rank[v][i] for i in C],[counts[i] for i in C],st)
            z["wide_common_count"]=len(B);z["all_variant_rank1_common_count"]=len(C)
            key="__".join((scope,mode,variant));panels[key]=z
            print("FD_SOURCE_AGG_CONTROL_PANEL",key,json.dumps(z,allow_nan=False,separators=(",",":")),flush=True)
    assert len(panels)==12
    out={"botanical_parent":4505,"geographic_parent":8264,"mainland_candidate_regions":748,"complete_pair_species":2180,
         "mainland_sha256":md,"source_pins":PIN,"metric_order":METRICS,"region_order":REGION_ORDER,"contrast_order":CONTRASTS,
         "genus_contrast_order":CONTRASTS[:6],"panels":panels,
         "warning":"Post-result reference aggregation diagnostic. Rank1 is first-ranked geographic or climate-matched candidate, not actual source population. Same retained islands and sample size across all modes, overlap and both reference types; equal-size expected diversity under independent SRSWOR per roster. Observed genus sensitivity only. Reference vs eligibility change is algebraic, not causal. No proof of comparable geographic units, measured pollinators, or evolutionary change."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    print("FD_SOURCE_AGG_CONTROL_COMPLETE",json.dumps({"panels":len(panels),"support_qualified":sum(z["status"]=="exploratory_estimable" for z in panels.values()),"metric_order":METRICS,"region_order":REGION_ORDER,"contrast_order":CONTRASTS,"genus_contrast_order":CONTRASTS[:6]}),flush=True)
if __name__=="__main__":main()
