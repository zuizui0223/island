#!/usr/bin/env python3
"""Protected GIFT4505 source-relative genuine floral pair vs disjoint animal guilds.
All original plant/source/geographic identities remain inside the protected runner."""
import argparse,collections,functools,hashlib,io,json,pathlib
import numpy as np
import pandas as pd
from scipy import sparse
from fd_original_tropical_pollinator_bridge import (
    PIN,read_pin,public_animals,pair_states,fit_block,null_p)
ENDPOINTS=("radial_deep","bilateral_deep","radial_deep_covariance","joint_Q")
def membership(rosters,index):
    rr=[];cc=[];n=np.array([len(s) for s in rosters])
    for i,s in enumerate(rosters):
        rr.extend([i]*len(s));cc.extend(index[v] for v in sorted(s))
    a=sparse.csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(rosters),len(index)))
    return sparse.diags(1/n)@a,n
def source_fit(d,pair,rosters,source,features,controls,nperm,seed):
    cols=features+controls
    d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=cols+["group"]).reset_index(drop=True)
    result={"n":len(d),"groups":int(d.group.nunique()),"features":features,"controls":controls}
    if len(d)<8:return result|{"status":"insufficient"}
    design=fit_block(d,cols)
    if design is None:return result|{"status":"insufficient_rank"}
    B,groups,dof,ninfo=design
    B=B[:len(features)]
    labels=list(features)
    contrast=None
    if "log_specialist" in features and "log_partial" in features:
        contrast=np.zeros(len(features));contrast[features.index("log_specialist")]=1;contrast[features.index("log_partial")]=-1
        B=np.vstack([B,contrast@B]);labels.append("specialist_minus_partial")
    tar=[rosters[i] for i in d.island_id];src=[source[i] for i in d.island_id]
    assert all(len(s)>=5 for s in src) and all(not a&b for a,b in zip(tar,src))
    species=sorted(set().union(*tar,*src));index={s:i for i,s in enumerate(species)}
    A,n=membership(tar,index);M,m=membership(src,index)
    tab=pair.loc[species].reset_index(drop=True)
    cats=sorted(set(tab.pair));ci={c:i for i,c in enumerate(cats)}
    code=np.array([ci[c] for c in tab.pair]);eye=np.eye(len(cats))
    T=np.array([[c==("actinomorphic","deep"),c==("zygomorphic","deep"),c[0]=="actinomorphic",c[1]=="deep"] for c in cats],float)
    def outcomes(a,size,c):
        p=a@eye[c];v=p@T
        return np.column_stack([v[:,:2],v[:,0]-v[:,2]*v[:,3],(1-np.square(p).sum(1))*size/(size-1)])
    y=outcomes(A,n,code);z=outcomes(M,m,code);delta=y-z
    obs=B@delta;bt=B@y;bs=B@z
    assert np.allclose(obs,bt-bs)
    rng=np.random.default_rng(seed)
    blocks=[np.array(v,int) for v in tab.groupby(["genus","grade"],sort=True).indices.values() if len(v)>1 and len(np.unique(code[v]))>1]
    if not blocks:return result|{"status":"no_movable_taxa"}
    draws=np.empty((nperm,len(labels),len(ENDPOINTS)))
    for k in range(nperm):
        c=code.copy()
        for ids in blocks:c[ids]=rng.permutation(c[ids])
        draws[k]=B@(outcomes(A,n,c)-outcomes(M,m,c))
    drops=[]
    for g in groups:
        keep=d.group.ne(g).to_numpy()
        other=fit_block(d.loc[keep].reset_index(drop=True),cols)
        if other is None:continue
        b=other[0][:len(features)]
        if contrast is not None:b=np.vstack([b,contrast@b])
        drops.append(b@delta[keep])
    rows={}
    for j,name in enumerate(labels):
        rows[name]={}
        for i,end in enumerate(ENDPOINTS):
            v=draws[:,j,i];actual=float(obs[j,i]);lo=[float(q[j,i]) for q in drops]
            rows[name][end]={"target_coefficient":float(bt[j,i]),"source_coefficient":float(bs[j,i]),"delta_coefficient":actual,"null_mean":float(v.mean()),"delta_minus_null":float(actual-v.mean()),"null_q025_q975":np.quantile(v,[.025,.975]).tolist(),"p":null_p(v,actual),"leave_one_group_range":[min(lo),max(lo)] if lo else None,"n_leave_one_group":len(lo)}
    return result|{"status":"exploratory","informative_groups":len(groups),"informative_islands":ninfo,"dof":dof,"passes_support_gate":len(groups)>=10 and ninfo>=30,"distinct_target_source_taxa":len(species),"movable_taxa":sum(map(len,blocks)),"median_target_scored":float(np.median(n)),"median_source_scored":float(np.median(m)),"mean_observed_delta":dict(zip(ENDPOINTS,delta.mean(0).tolist())),"outcomes":rows}
def main():
    p=argparse.ArgumentParser();p.add_argument("--original-root",required=True,type=pathlib.Path);p.add_argument("--geography",required=True,type=pathlib.Path);p.add_argument("--bridge",required=True,type=pathlib.Path);p.add_argument("--out",required=True,type=pathlib.Path);p.add_argument("--permutations",type=int,default=199);a=p.parse_args()
    r=a.original_root
    flora=pd.read_csv(io.BytesIO(read_pin(r/"fixed/canonical/input/chapter1_status_flora.csv.gz","flora")),compression="gzip",dtype=str,usecols=["island_id","accepted_species","origin_status"])
    t=pd.read_csv(io.BytesIO(read_pin(r/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz","traits")),compression="gzip",dtype=str,usecols=["accepted_species","trait_name","normalized_value","evidence_scope","evidence_quality"])
    assigned=pd.read_csv(io.BytesIO(read_pin(r/"fixed/source/result/island_source_assignments.csv.gz","assign")),compression="gzip",dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    raw=(r/"fixed/source/source/gift_native_mainland_flora.csv.gz").read_bytes();mainland_sha=hashlib.sha256(raw).hexdigest()
    mainland=pd.read_csv(io.BytesIO(raw),compression="gzip",dtype=str,usecols=["entity_ID","work_species"])
    geo=pd.read_csv(io.BytesIO(read_pin(a.geography,"geography")),dtype={"island_id":str})
    assert flora.island_id.nunique()==4505 and not flora.duplicated(["island_id","accepted_species"]).any() and geo.island_id.nunique()==8264 and geo.island_id.is_unique
    assert mainland.entity_ID.nunique()==748 and not mainland.duplicated(["entity_ID","work_species"]).any()
    groups=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    full=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    ass=assigned[assigned.source_mode.isin(["geo_k5","geo50_climate10"])].copy()
    ass.source_rank=pd.to_numeric(ass.source_rank,errors="raise").astype(int)
    assert not ass.duplicated(["island_id","source_mode","source_rank"]).any()
    assert {k:int(ass.source_mode.eq(k).sum()) for k in ("geo_k5","geo50_climate10")}=={"geo_k5":41325,"geo50_climate10":82190}
    keys=ass.sort_values("source_rank").groupby(["island_id","source_mode"]).entity_ID.agg(tuple).to_dict()
    @functools.lru_cache(maxsize=64)
    def pool(ids):return frozenset().union(*(groups.get(i,set()) for i in ids))
    tropical=set(geo.loc[geo.analysis_regime.eq("tropical"),"island_id"])
    assert len(set(full)&tropical)==1558
    d=public_animals(a.bridge);d=d[d.island_id.isin(tropical)].copy()
    for col in ("nectar_gt0","nectar_ge50","all_birds"):d[col]=pd.to_numeric(d[col],errors="coerce")
    ok=d.nectar_gt0.notna()
    d["specialist"]=d.nectar_ge50;d["partial"]=d.nectar_gt0-d.nectar_ge50;d["other"]=d.all_birds-d.nectar_gt0
    assert (d.loc[ok,["specialist","partial","other"]]>=0).all().all()
    for c in ("specialist","partial","other"):d["log_"+c]=np.log1p(d[c])
    pa=pair_states(t);panels={};support={}
    for scope,minimum in (("all_recorded",5),("all_recorded",30),("original_native",5),("high_quality_only",5)):
        f=flora[flora.island_id.isin(tropical)]
        if scope=="original_native":f=f[f.origin_status.eq("native")]
        pair=pa if scope!="high_quality_only" else pa[pa.grade.eq("high")]
        scored=set(pair.index)
        tar=f[f.island_id.isin(d.island_id)&f.accepted_species.isin(scored)].groupby("island_id").accepted_species.agg(set).to_dict()
        for mode,nr in (("geo_k5",5),("geo50_climate10",10)):
            src={}
            for i in tar:
                ix=keys.get((i,mode),())
                if ix:
                    assert len(ix)==nr
                    # identical published source construction: union mainland candidates, remove entire island flora.
                    src[i]=(set(pool(ix))-full[i])&scored
            selected=[i for i,s in tar.items() if len(s)>=minimum and len(src.get(i,set()))>=5]
            q=d[d.island_id.isin(selected)].copy()
            q["log_n"]=np.log2(q.island_id.map(lambda i:len(tar[i])))
            q["log_source_n"]=np.log2(q.island_id.map(lambda i:len(src[i])))
            cv=["log_area","log_distance","log_n","log_source_n"]
            name=f"{scope}_n{minimum}_{mode}"
            support[name]={"matched_pair_islands":len(q),"bird_candidates":int(q.specialist.notna().sum())}
            cases=[("bee",q,["log_bee"],cv),("bird_disjoint",q,["log_bee","log_specialist","log_partial","log_other"],cv)]
            for mask in ("all","precision_le1km","no_gbifid_with_reference"):
                b=q[q[mask+"_records"].ge(10)&q[mask+"_recorded_species"].gt(0)].copy()
                b["log_butterfly"]=np.log1p(b[mask+"_recorded_species"]);b["log_effort"]=np.log1p(b[mask+"_records"])
                cases.append(("butterfly_"+mask,b,["log_bee","log_butterfly"],cv+["log_effort"]))
            for label,b,features,controls in cases:
                k=name+"_"+label
                panels[k]=source_fit(b,pair,tar,src,features,controls,a.permutations,20261124)
    tests=[v for p in panels.values() for o in p.get("outcomes",{}).values() for v in o.values()]
    order=sorted(tests,key=lambda x:x["p"]);last=1.
    for j in range(len(order)-1,-1,-1):
        last=min(last,order[j]["p"]*len(order)/(j+1));order[j]["BH_q_exploratory"]=float(last)
    out={"status":"ORIGINAL_GIFT4505_SOURCE_RELATIVE_ANIMAL_GUILD_EXPLORATORY","original_botanical_parent":4505,"original_tropical_parent":1558,"original_geographic_parent":8264,"mainland_regions":748,"mainland_input_sha256":mainland_sha,"source_pins":PIN,"support":support,"panel_count":len(panels),"test_count":len(tests),"panels":panels,"warning":"Potential mainland source roster minus whole island flora, not historical ancestors. Joint same-species Direct categories, not local measured phenotypes or pollination. Bird crosswalk provisional, butterfly positive-only record richness. Shared distinct-species genus x quality permutation across island and source; controls include scored island and source richness. All tests exploratory. No original taxon/source/geographic identities exported."}
    assert len(panels)==40
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    print("SOURCE_RELATIVE_SUPPORT",json.dumps(support),flush=True)
    for k,v in panels.items():print("SOURCE_RELATIVE_PANEL",k,json.dumps(v),flush=True)
if __name__=="__main__":main()
