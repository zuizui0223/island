#!/usr/bin/env python3
"""FD diagnostic: architecture-dependent color documentation, not a biological mechanism."""
import argparse,functools,hashlib,io,json,pathlib
import numpy as np
import pandas as pd
from fd_four_axis_space import PIN,read_pin,states,pair_states,REGIONS
from fd_complete_case_reversal import tier_sets
from fd_global_source_architecture import region_design
from fd_original_source_relative_pollinator import membership
SEED=20261125
PRIMARY="geo_k5__disjoint"

def audit(d,tar,src,pair,color,perms=999):
    if not len(d):return {"status":"empty"}
    ids=d.island_id.tolist();rr=[tar[i] for i in ids];ss=[src[i] for i in ids]
    names=sorted(set().union(*rr,*ss));ix={s:i for i,s in enumerate(names)}
    A,n=membership(rr,ix);M,m=membership(ss,ix)
    t=pair.loc[names];v=np.array([s in color for s in names],float)
    rmask=np.array([a=="actinomorphic" and b=="deep" for a,b in t.pair])
    bmask=np.array([a=="zygomorphic" and b=="deep" for a,b in t.pair])
    G,reg,bid,blocks=region_design(d)
    groups=t.groupby(["genus","grade"],sort=True).indices
    groups=[np.array(z,int) for z in groups.values()]
    movable=np.zeros(len(t),bool);E=v.copy()
    for g in groups:
        E[g]=v[g].mean()
        if np.ptp(v[g])>0:movable[g]=True
    weights={};details={}
    for j,name in enumerate(REGIONS):
        mask=reg==j
        if not mask.any():continue
        W=np.asarray(G[mask].mean(0)).ravel()
        details[name]={"islands":int(d.region.eq(name).sum()),"blocks":int(mask.sum()),"roles":{}}
        for role,mat in [("target",A),("donor",M)]:
            w=np.asarray(mat.T@W).ravel()
            parts={};co=[]
            for label,k in [("radial_deep",rmask),("bilateral_deep",bmask)]:
                h=w*k;den=h.sum();present=h>0
                row={"pair_share":float(den),"distinct_species":int(present.sum()),"color_documented_species":int((present&(v>0)).sum()),"movable_species":int((present&movable).sum()),"movable_genera":int(t.loc[present&movable,"genus"].nunique())}
                if den>0:
                    h=h/den;co.append(h)
                    row|={"weighted_color_coverage":float(h@v),"genus_grade_expected_coverage":float(h@E)}
                parts[label]=row
            details[name]["roles"][role]=parts
            if len(co)==2:weights[name+"__"+role]=co[0]-co[1]
        a=name+"__target";b=name+"__donor"
        if a in weights and b in weights:weights[name+"__target_minus_donor"]=weights[a]-weights[b]
    if "north__target_minus_donor" in weights and "tropical__target_minus_donor" in weights:
        weights["tropical_minus_north__target_minus_donor"]=weights["tropical__target_minus_donor"]-weights["north__target_minus_donor"]
    if not weights:return {"status":"insufficient","coverage":details}
    keys=list(weights);H=np.vstack([weights[k] for k in keys])
    obs=H@v;ex=H@E;res=obs-ex
    rng=np.random.default_rng(SEED);draw=np.empty((perms,len(keys)))
    mv=[g for g in groups if np.ptp(v[g])>0]
    for k in range(perms):
        u=v.copy()
        for g in mv:u[g]=rng.permutation(v[g])
        draw[k]=H@u
    out={}
    for i,k in enumerate(keys):
        p=min(1.,2*min((1+int(np.sum(draw[:,i]<=obs[i]+1e-12)))/(perms+1),(1+int(np.sum(draw[:,i]>=obs[i]-1e-12)))/(perms+1)))
        relevant=["north","tropical"] if k.startswith("tropical_minus_north") else [k.split("__")[0]]
        gate=all(details[r]["islands"]>=200 and details[r]["blocks"]>=30 for r in relevant)
        active=(np.abs(H[i])>1e-14)&movable
        ng=int(t.loc[active,"genus"].nunique())
        gate=gate and active.sum()>=30 and ng>=10
        out[k]={"observed_coverage_contrast":float(obs[i]),"genus_grade_expected":float(ex[i]),"residual":float(res[i]),"p_coherent":float(p),"null_interval":np.quantile(draw[:,i],[.025,.975]).tolist(),"movable_weighted_species":int(active.sum()),"movable_weighted_genera":ng,"gate":bool(gate)}
    return {"status":"estimated","islands":len(d),"distinct_pair_species":len(names),"coverage":details,"tests":out}

def self_test():
    names=["Plant "+str(i) for i in range(8)]
    tab=pd.DataFrame({"pair":[("actinomorphic","deep")]*4+[("zygomorphic","deep")]*4,"genus":"Plant","grade":"high"},index=names)
    d=pd.DataFrame({"island_id":[1,2],"region":["north","tropical"],"block":["a","b"]})
    tar={i:set(names) for i in [1,2]};src=tar
    z=audit(d,tar,src,tab,set(names[:4]),perms=19)
    assert z["tests"]["north__target"]["observed_coverage_contrast"]==1
    assert z["tests"]["north__target_minus_donor"]["observed_coverage_contrast"]==0
    print("FD_COLOR_DOCUMENTATION_SELF_TEST_PASS",flush=True)
def main():
    p=argparse.ArgumentParser();p.add_argument("--self-test",action="store_true");p.add_argument("--original-root",type=pathlib.Path);p.add_argument("--geography",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);a=p.parse_args()
    if a.self_test:self_test();return
    r=a.original_root
    f=pd.read_csv(io.BytesIO(read_pin(r/"fixed/canonical/input/chapter1_status_flora.csv.gz","flora")),compression="gzip",dtype=str,usecols=["island_id","accepted_species","origin_status"])
    t=pd.read_csv(io.BytesIO(read_pin(r/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz","traits")),compression="gzip",dtype=str)
    ass=pd.read_csv(io.BytesIO(read_pin(r/"fixed/source/result/island_source_assignments.csv.gz","assign")),compression="gzip",dtype=str)
    raw=(r/"fixed/source/source/gift_native_mainland_flora.csv.gz").read_bytes()
    mainland=pd.read_csv(io.BytesIO(raw),compression="gzip",dtype=str,usecols=["entity_ID","work_species"])
    geo=pd.read_csv(io.BytesIO(read_pin(a.geography,"geography")),dtype={"island_id":str})
    assert f.island_id.nunique()==4505 and not f.duplicated(["island_id","accepted_species"]).any()
    assert mainland.entity_ID.nunique()==748 and geo.island_id.is_unique and len(geo)==8264
    full=f.groupby("island_id").accepted_species.agg(set).to_dict();pools=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    ass=ass[ass.source_mode.isin(["geo_k5","geo50_climate10"])].copy();ass.source_rank=pd.to_numeric(ass.source_rank)
    assert not ass.duplicated(["island_id","source_mode","source_rank"]).any()
    keys=ass.sort_values("source_rank").groupby(["island_id","source_mode"]).entity_ID.agg(tuple).to_dict()
    @functools.lru_cache(maxsize=128)
    def pool(ids):return frozenset().union(*(pools[i] for i in ids))
    g=geo[geo.island_id.isin(full)].dropna(subset=["island_latitude","island_longitude"]).copy()
    g["region"]=np.where(g.analysis_regime.eq("tropical"),"tropical",np.where(g.island_latitude.gt(23.5),"north",np.where(g.island_latitude.lt(-23.5),"south","unclassified")))
    g=g[g.region.isin(REGIONS)].copy();g["block"]=np.floor(g.island_latitude/10).astype(int).astype(str)+"_"+np.floor(g.island_longitude/10).astype(int).astype(str)
    pair=pair_states(t);scored=set(pair.index);color=tier_sets(t,pair,"all_recorded")["color"]
    tar=f[f.accepted_species.isin(scored)].groupby("island_id").accepted_species.agg(set).to_dict()
    panels={}
    for mode,nr in [("geo_k5",5),("geo50_climate10",10)]:
        incl={};disj={}
        for i in g.island_id:
            if len(tar.get(i,set()))<5:continue
            k=keys.get((i,mode),())
            if len(k)!=nr or any(x not in pools for x in k):continue
            s=pool(k)&scored;v=s-full[i]
            if len(v)<5:continue
            incl[i]=s;disj[i]=v
        d=g[g.island_id.isin(disj)].copy()
        for label,src in [("inclusive",incl),("disjoint",disj)]:
            key=mode+"__"+label
            panels[key]=audit(d,tar,src,pair,color)
    names=[r+"__"+k for r in ["north","tropical"] for k in ["target","donor","target_minus_donor"]]+["tropical_minus_north__target_minus_donor"]
    rows=[]
    for k in names:
        row=panels[PRIMARY].get("tests",{}).get(k,{})
        rows.append((row.get("p_coherent",1.) if row.get("gate") else 1.,row))
    last=0.
    for j,(pv,row) in enumerate(sorted(rows,key=lambda q:q[0])):
        last=max(last,min(1.,pv*(len(rows)-j)));row["primary_holm_p"]=last
    out={"botanical_parent":4505,"geography_parent":8264,"source_pins":PIN,"mainland_sha256":hashlib.sha256(raw).hexdigest(),"primary":PRIMARY,"primary_family_size":7,"panels":panels,"warning":"Post-result documentation audit conditional on already pair-documented species. Coverage is ratio of block-equal incidence weights, not mean of island-specific ratios. Genus and pair evidence grade conditional species-coherent color-availability permutations; does not identify recording bias versus biology, missing trait values, or actual pollination. No raw IDs exported."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False,indent=2))
    for k,v in panels.items():print("FD_COLOR_DOCUMENTATION_PANEL",k,json.dumps(v,allow_nan=False),flush=True)
    print("FD_COLOR_DOCUMENTATION_COMPLETE",json.dumps({"panels":len(panels),"primary_family_size":7}),flush=True)
if __name__=="__main__":main()
