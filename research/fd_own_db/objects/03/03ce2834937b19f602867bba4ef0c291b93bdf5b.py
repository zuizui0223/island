#!/usr/bin/env python3
"""Bird-species diet-label null, conditional on fixed original checklist incidence.
Tests a specific PUBLIC tubular flower association; not realized pollination."""
import argparse,json,math,pathlib,urllib.parse
import numpy as np
import pandas as pd
from tropical_pollinator_source_intake import MATRICES,matrix,get_diet,norm_name
from tropical_pollinator_source_preflight import BASE,get,blobsha
from tropical_bird_hwi_occupancy import load_hwi
from tropical_public_flora_guild_association import read_artifact,get_public_flora
from tropical_three_guild_candidate_audit import fetch

PANELS=(("direct",5),("direct",50),("direct_single_state",5),("direct_single_state",50))
ENDPOINTS=("tubular_only","extended_access_shape","rao_q_without_replacement")
def prepare_incidence(b):
    diet,_=get_diet();hwi,_=load_hwi();present={}
    for group,v in b.groupby("group"):
        path,expected=MATRICES[group]
        raw=get(BASE+urllib.parse.quote(path));assert blobsha(raw)==expected
        islands,sp,area,footer=matrix(raw)
        ii={name:j for j,name in enumerate(islands)}
        for r in v.itertuples():
            j=ii[r.bird_name]
            ss={norm_name(s) for s,p in sp.items() if p[j]}
            assert len(ss)==r.all_birds and ss<=diet.keys() and ss<=hwi.keys()
            assert sum(diet[s]>=50 for s in ss)==r.nectar_ge50
            assert sum(diet[s]>0 for s in ss)==r.nectar_gt0
            present[r.bee_id]=ss
    species=sorted(set().union(*present.values()))
    labels=np.array([2 if diet[s]>=50 else 1 if diet[s]>0 else 0 for s in species])
    fam=np.array([hwi[s][1] for s in species]);gen=np.array([s.split()[0] for s in species])
    return present,species,labels,fam,gen

def model_design(d,counts,mode):
    base=np.column_stack([np.log1p(d.native_bees),np.log2(d.bee_area_km2),np.log2(1+d.dist_km),np.log2(d.n_names)])
    other=counts[:,0];partial=counts[:,1];special=counts[:,2]
    if mode=="total":return np.column_stack([base,np.log1p(other+partial+special)])
    if mode=="broad":return np.column_stack([base,np.log1p(partial+special),np.log1p(other)])
    return np.column_stack([base,np.log1p(special),np.log1p(partial),np.log1p(other)])

def fit(X,Y,g):
    xs=[];ys=[]
    for idx in g:
        if len(idx)<2:continue
        x=X[idx];y=Y[idx];w=1/math.sqrt(len(idx))
        xs.append((x-x.mean(0))*w);ys.append((y-y.mean(0))*w)
    if not xs:return None
    x=np.vstack(xs);y=np.vstack(ys)
    if len(x)-len(xs)-x.shape[1]<4 or np.linalg.matrix_rank(x)<x.shape[1]:return None
    return np.linalg.lstsq(x,y,rcond=None)[0]
def holdout(d,Y,counts,groups):
    allidx=np.arange(len(d));loss=[]
    X=[model_design(d,counts,z) for z in ("total","broad","split")]
    for idx in groups:
        if len(idx)<2:continue
        tr=np.setdiff1d(allidx,idx)
        train_group=[np.flatnonzero(d.group.iloc[tr].to_numpy()==k) for k in sorted(d.group.iloc[tr].unique())]
        preds=[]
        for x in X:
            b=fit(x[tr],Y[tr],train_group)
            if b is None:break
            p=x[idx];preds.append((p-p.mean(0))@b)
        if len(preds)!=3:continue
        y=Y[idx]-Y[idx].mean(0);loss.append(np.stack([np.square(y-v).mean(0) for v in preds]))
    if not loss:return None
    a=np.stack(loss);m=a.mean(0)
    return {"groups":len(a),"endpoints":{k:{"total":float(m[0,j]),"broad":float(m[1,j]),"split":float(m[2,j]),"broad_minus_split":float(m[1,j]-m[2,j]),"total_minus_split":float(m[0,j]-m[2,j]),"split_beats_broad_groups":int((a[:,1,j]>a[:,2,j]).sum())} for j,k in enumerate(ENDPOINTS)}}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=pathlib.Path);p.add_argument("--draws",type=int,default=1999);a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    bird=fetch("bee_bird");fl,sha=get_public_flora(read_artifact())
    d0=bird[bird.group!="Ryukyu sensitivity"].copy()
    present,sp,labels,fam,gen=prepare_incidence(d0);ix={s:i for i,s in enumerate(sp)}
    out={};rng=np.random.default_rng(20261124);eye=np.eye(3,dtype=float)
    for scope,nmin in PANELS:
        f=fl[scope];d=f[(f.analysis_regime=="tropical")&(f.n_names>=nmin)].merge(d0,on="bee_id",validate="one_to_one").reset_index(drop=True)
        groups=[v.index.to_numpy() for _,v in d.groupby("group")]
        A=np.zeros((len(d),len(sp)))
        for j,id in enumerate(d.bee_id):A[j,[ix[s] for s in present[id]]]=1.
        counts=A@eye[labels]
        assert np.array_equal(counts[:,2],d.nectar_ge50) and np.array_equal(counts[:,1]+counts[:,2],d.nectar_gt0)
        Y=d[list(ENDPOINTS)].to_numpy(float)
        b=fit(model_design(d,counts,"split"),Y,groups)
        result={"n":len(d),"groups":sum(len(g)>=2 for g in groups),"source_distinct_birds":len(sp)}
        if b is None:out[f"{scope}_n{nmin}"]=result|{"status":"insufficient"};continue
        obs=b[4]-b[5];null={}
        for tax,key in (("family",fam),("genus",gen)):
            blocks=[np.flatnonzero(key==k) for k in sorted(set(key))]
            blocks=[v for v in blocks if len(v)>1 and len(np.unique(labels[v]))>1]
            draws=[]
            for i in range(a.draws):
                lab=labels.copy()
                for v in blocks:lab[v]=rng.permutation(lab[v])
                beta=fit(model_design(d,A@eye[lab],"split"),Y,groups)
                if beta is not None:draws.append(beta[4]-beta[5])
            z=np.asarray(draws)
            null[tax]={"movable_species":sum(map(len,blocks)),"blocks":len(blocks),"valid_draws":len(z),"endpoints":{}}
            if len(z):
                for j,k in enumerate(ENDPOINTS):
                    v=z[:,j];pval=(1+int((v>=obs[j]).sum()))/(len(v)+1)
                    null[tax]["endpoints"][k]={"null_q025_q50_q975":np.quantile(v,[.025,.5,.975]).tolist(),"one_sided_positive_p":pval}
        result|={"status":"exploratory","specialist_minus_partial":dict(zip(ENDPOINTS,obs.tolist())),"null":null,"holdout_within_archipelago":holdout(d,Y,counts,groups)}
        out[f"{scope}_n{nmin}"]=result
    result={"profile_sha256":sha,"panels":out,"warning":"Original fixed bird checklist species incidence, global family/genus diet-label shuffles preserve source-pool class counts, not local class counts or actual pollination. Valid rank-conditional draws only; sparse overlapping cohorts, one-sided descriptive p-values unadjusted. Heldout compares within-group centered outcomes, not absolute floras."}
    (a.out/"diet_label_null.json").write_text(json.dumps(result,indent=2,allow_nan=False))
    print("DIET_LABEL_NULL",json.dumps(result),flush=True)
if __name__=="__main__":main()
