"""Region-specific, globally coherent within-genus whole-profile null.

Matched to the corrected isolation triple endpoint. No original identities leave
the runner. The two focal tests were selected after the previous FD result:
this is a retrospective falsification, not a preregistered discovery.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import json
from math import log
from pathlib import Path
import re

import numpy as np
import pandas as pd
import yaml

from scripts.fd_self_db_isolation import (
    C, F, S, DIST, COV, LOCK, METRICS, REGIONS, SCOPES, VARIANTS, REPEATS,
    direct_states, digest, cell_expectation, seed_for, fit_cluster, adjust_fixed,
)
from island_v2.chapter1_h2_complete_selfing_sensitivity import build_species_score

N = 30
B = 399
SEED = 20261126
ENDPOINT = "triple_cfs_n30"
FOCAL = (("northern_midlatitude", "excess_mi"), ("tropical", "marginal_q"))


def genus_key(name):
    """Use accepted binomial's genus token; unresolved/hybrid tokens singleton."""
    parts = str(name).split()
    if len(parts) >= 2 and re.fullmatch(r"[A-Z][a-z]+", parts[0]):
        return parts[0]
    return "unresolved:" + str(name)


def random_assignments(codes, genera, draws, seed):
    """Row 0 observed; a SINGLE species->whole-profile map per regional draw."""
    codes = np.asarray(codes, dtype=np.int32)
    groups = defaultdict(list)
    for i, g in enumerate(genera):
        groups[g].append(i)
    result = np.tile(codes, (draws + 1, 1))
    variable = np.zeros(len(codes), dtype=bool)
    rng = np.random.default_rng(seed)
    for ids in groups.values():
        if len(ids) < 2 or len(set(codes[ids])) < 2:
            continue
        variable[ids] = True
        # Independent uniform permutations of species positions, no replacement.
        orders = np.argsort(rng.random((draws, len(ids))), axis=1)
        result[1:, ids] = codes[np.asarray(ids)[orders]]
    return result, variable


def dense_counts(codes, k):
    x = np.asarray(codes, dtype=np.int64)
    r = x.shape[0]
    return np.bincount((x + np.arange(r)[:,None]*k).ravel(),
                       minlength=r*k).reshape(r,k)


def coding(profiles):
    # The labels are indices into actual observed whole profiles, not hypothetical
    # Cartesian combinations. Marginal indices refer to those same profiles.
    axes = [np.unique(profiles[:,j], return_inverse=True)[1] for j in range(3)]
    arch = np.unique(profiles[:,1:],axis=0,return_inverse=True)[1]
    return axes, arch


def all_metrics(p, profiles, sample_indices, n=N):
    """Exact expected rarefied Q and sampled MI for all coherent assignments.

    p rows are observed + regional permutations; columns are the identical
    sorted target species. Reuse the exact same 30 parent rarefaction subsets
    in every assignment, hence preserve their sampled genus counts as well.
    """
    nrows, total = p.shape
    if total < n:
        raise ValueError("insufficient same-species carriers")
    axes, arch = coding(profiles)
    factor=(n-1)/n*total/(total-1)
    qs=[]
    for axis in axes:
        cnt=dense_counts(axis[p], int(axis.max())+1)
        qs.append(factor*(1-np.square(cnt/total).sum(axis=1)))
    joint_counts=dense_counts(p,len(profiles))
    joint=factor*(1-np.square(joint_counts/total).sum(axis=1))
    # Chunk by permutations to bound memory use.
    table=np.zeros((n+1,n+1))
    for a in range(1,n+1):
        for b in range(1,n+1):
            table[a,b]=cell_expectation(n,a,b)
    clog=np.zeros(n+1)
    c=np.arange(1,n+1)
    clog[1:]=c*np.log(c)
    mi=np.empty(nrows)
    for start in range(0,nrows,40):
        stop=min(start+40,nrows)
        samples=p[start:stop, sample_indices].reshape(-1,n)
        col=dense_counts(axes[0][samples],int(axes[0].max())+1)
        ar=dense_counts(arch[samples],int(arch.max())+1)
        joint_s=dense_counts(samples,len(profiles))
        info=(clog[joint_s].sum(1)-clog[col].sum(1)-clog[ar].sum(1)+n*log(n))/n
        expect=table[col[:,:,None],ar[:,None,:]].sum(axis=(1,2))
        mi[start:stop]=(info-expect).reshape(stop-start,-1).mean(1)
    return np.column_stack((np.mean(qs,axis=0),joint,mi))


def parent_draws(total, island, scope, n=N):
    if total == n:
        return np.arange(n)[None,:]
    rng=np.random.default_rng(seed_for(island,ENDPOINT,scope))
    return np.array([rng.choice(total,size=n,replace=False) for _ in range(REPEATS)])


def linear_weights(work,predictors):
    x=work[predictors].to_numpy(float)
    x=(x-x.mean(0))/x.std(0)
    design=np.column_stack((np.ones(len(x)),x))
    return np.linalg.pinv(design)[1,:]


def perm_p(observed, null):
    # Inclusive two-tail randomization test; +1 correction, no zero p.
    tol=1e-12
    return float(min(1,2*min(1+int(np.sum(null>=observed-tol)),
                             1+int(np.sum(null<=observed+tol)))/(len(null)+1)))


def null_qvalues(rows):
    # Same fixed 12-slot family for each scope x variant, independent of support.
    for key in {(r["scope"],r["variant"]) for r in rows}:
        group=[r for r in rows if (r["scope"],r["variant"])==key]
        assert len(group)==12
        order=sorted(group,key=lambda r:r.get("perm_p",1.))
        prev=1.
        harmonic=sum(1/i for i in range(1,13))
        for rank,r in reversed(list(enumerate(order,1))):
            prev=min(prev,12*r.get("perm_p",1.)/rank)
            if "perm_p" in r:
                r["perm_q_bh"]=float(prev)
                r["perm_q_by"]=float(min(1,prev*harmonic))
    targets=[r for r in rows if r["scope"]=="all_recorded"
             and r["variant"]=="coverage_adjusted" and (r["region"],r["metric"]) in FOCAL]
    assert len(targets)==2
    prev=0.
    for i,r in enumerate(sorted(targets,key=lambda x:x.get("perm_p",1.))):
        prev=max(prev,min(1,(2-i)*r.get("perm_p",1.)))
        if "perm_p" in r:
            r["focal_holm2"]=float(prev)


def run(flora,traits,geo,scores,parent,draws=B):
    if flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("duplicate island/species")
    complete=traits.dropna(subset=[C,F,S])[[C,F,S]]
    if len(complete)!=2234 or flora.island_id.nunique()!=4505:
        raise ValueError("locked morphology universe changed")
    if geo.island_id.duplicated().any() or len(geo)!=8264:
        raise ValueError("corrected geography changed")
    geo=geo[["island_id","analysis_regime","spatial_block",DIST,*COV]].copy()
    for c in (DIST,*COV):
        geo[c]=pd.to_numeric(geo[c],errors="coerce")
    score=scores.set_index("accepted_species").score.dropna()
    rows=[]; support=[]; validations=[]
    parent_lookup={(r["scope"],r["region"],r["variant"],r["metric"]):r
                   for r in parent["models"] if r["endpoint"]==ENDPOINT}
    for scope in SCOPES:
        f=flora if scope=="all_recorded" else flora.loc[flora.origin_status.eq("native")]
        totals=f.groupby("island_id").size()
        rep=f[["island_id","accepted_species"]].copy()
        rep["score"]=rep.accepted_species.map(score)
        rep=rep.groupby("island_id").score.agg(["mean","count"]).rename(columns={"mean":"selfing_core","count":"selfing_n"})
        joined=f[["island_id","accepted_species"]].merge(complete,left_on="accepted_species",right_index=True,validate="many_to_one")
        counts=joined.groupby("island_id").size()
        if int((counts>=N).sum())!=(928 if scope=="all_recorded" else 87):
            raise ValueError("original triple support changed")
        # Pool contains every Direct-complete species recorded in ANY island in
        # this region/scope, including islands below the target n=30 gate.
        joined=joined.merge(geo[["island_id","analysis_regime"]],on="island_id",validate="many_to_one")
        for region in REGIONS:
            reg=joined.loc[joined.analysis_regime.eq(region)]
            eligible=[(i,p.sort_values("accepted_species")) for i,p in reg.groupby("island_id",sort=True) if len(p)>=N]
            names=sorted(reg.accepted_species.unique())
            profiles, codes=np.unique(complete.loc[names].to_numpy(str),axis=0,return_inverse=True)
            genera=[genus_key(s) for s in names]
            region_seed=int.from_bytes(sha256(f"{SEED}:{scope}:{region}".encode()).digest()[:8],"big")
            assignments,variable=random_assignments(codes,genera,draws,region_seed)
            index={s:i for i,s in enumerate(names)}
            rec=[]; outcomes=[]
            for island,p in eligible:
                ids=np.array([index[s] for s in p.accepted_species],dtype=int)
                own=all_metrics(assignments[:,ids],profiles,parent_draws(len(p),island,scope))
                shared=p.accepted_species.map(score).dropna()
                rec.append({"island_id":island,"log_carriers":log(len(p)),
                    "coverage":len(p)/totals[island],"same_carrier_selfing":float(shared.mean()) if len(shared) else np.nan,
                    "same_carrier_n":len(shared),"variable_fraction":float(variable[ids].mean()),
                    **{m:float(own[0,j]) for j,m in enumerate(METRICS)}})
                outcomes.append(own)
            support.append({"scope":scope,"region":region,"pool_species":len(names),
                "pool_genera":len(set(genera)),"variable_pool_species":int(variable.sum()),
                "unresolved_genus_singletons":sum(x.startswith("unresolved:") for x in genera),
                "eligible_islands":len(rec),
                "median_exchangeable_carrier_fraction":float(np.median([r["variable_fraction"] for r in rec])) if rec else None})
            if rec:
                private=pd.DataFrame(rec).merge(geo,on="island_id",how="left",validate="one_to_one").merge(rep,left_on="island_id",right_index=True,how="left",validate="one_to_one")
                cube=np.stack(outcomes)
            else:
                private=pd.DataFrame(columns=["island_id","log_carriers","coverage","same_carrier_selfing","same_carrier_n","selfing_core","selfing_n","spatial_block",DIST,*COV,*METRICS])
                cube=np.empty((0,draws+1,3))
            base=private.dropna(subset=[DIST,*COV,"spatial_block",*METRICS])
            coverage=base.dropna(subset=["log_carriers","coverage"])
            h2=coverage.loc[coverage.selfing_n.ge(10)].dropna(subset=["selfing_core"])
            same=coverage.loc[coverage.same_carrier_n.ge(10)].dropna(subset=["same_carrier_selfing"])
            baseline=[DIST,*COV]; control=baseline+["log_carriers","coverage"]
            designs={"baseline":(base,baseline),"coverage_adjusted":(coverage,control),
                     "h2_matched_base":(h2,control),"h2_matched_adjusted":(h2,control+["selfing_core"]),
                     "h2_same_carrier_base":(same,control),"h2_same_carrier_adjusted":(same,control+["same_carrier_selfing"])}
            for variant in VARIANTS:
                sample,preds=designs[variant]
                for j,m in enumerate(METRICS):
                    original=fit_cluster(sample,m,preds)
                    prev=parent_lookup[(scope,region,variant,m)]
                    if original["n"]!=prev["n"] or original["blocks"]!=prev["blocks"] or original["status"]!=prev["status"]:
                        raise ValueError("parent model support/status mismatch: "+str((scope,region,variant,m)))
                    row={"scope":scope,"region":region,"variant":variant,"metric":m,"family":"triple",
                         "endpoint":ENDPOINT,"n":original["n"],"blocks":original["blocks"],"status":original["status"]}
                    if original["status"]=="fit":
                        err=abs(original["beta"]-prev["beta"])
                        if err>1e-10:
                            raise ValueError("observed parent slope not reproduced")
                        work=sample[[m,*preds,"spatial_block"]].replace([np.inf,-np.inf],np.nan).dropna()
                        work=work.loc[work.spatial_block.astype(str).str.strip().ne("")]
                        a=cube[work.index.to_numpy(),:,j]
                        w=linear_weights(work,preds)
                        slopes=w@a
                        if abs(slopes[0]-original["beta"])>1e-10:
                            raise ValueError("permutation regression alignment failure")
                        res=work.copy(); res["residual"]=a[:,0]-a[:,1:].mean(1)
                        fit=fit_cluster(res,"residual",preds)
                        row.update(observed_beta=float(slopes[0]),expected_beta=float(slopes[1:].mean()),
                            null_lo=float(np.quantile(slopes[1:],.025)),null_hi=float(np.quantile(slopes[1:],.975)),
                            null_sd=float(slopes[1:].std()),perm_p=perm_p(slopes[0],slopes[1:]),
                            residual_beta=float(slopes[0]-slopes[1:].mean()),
                            residual_fit_status=fit["status"])
                        if fit["status"]=="fit":
                            if abs(fit["beta"]-row["residual_beta"])>1e-10:
                                raise ValueError("additive slope decomposition failed")
                            row.update(beta=fit["beta"],se=fit["se"],lo=fit["lo"],hi=fit["hi"],p=fit["p"])
                        else:
                            row["status"]=fit["status"]
                        validations.append(err)
                    rows.append(row)
    adjust_fixed(rows)
    null_qvalues(rows)
    return {"status":"COHERENT_REGIONAL_GENUS_NULL_EXPLORATORY","draws":draws,
            "source_sha256":LOCK,"original_slots":len(rows),"support":support,"models":rows,
            "n_reproduced_parent_fits":len(validations),"max_parent_beta_difference":max(validations,default=0.),
            "not_causal":True,"no_row_level_export":True}


def main():
    p=argparse.ArgumentParser()
    for key in ("flora","ledger","geography","state-audit","syndrome-config","parent","output"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    for k,v in LOCK.items():
        if digest(getattr(a,k))!=v:
            raise ValueError("source SHA mismatch: "+k)
    if digest(a.parent)!="d06af1e649bd7070f251927dfafe1000ac08f951355c4bf514cf87950a692871":
        raise ValueError("parent full aggregate changed")
    parent=json.loads(a.parent.read_text())
    for k,path in (("state_audit",a.state_audit),("syndrome_config",a.syndrome_config)):
        if digest(path)!=parent["source_sha256"][k]:
            raise ValueError("parent selfing source changed: "+k)
    flora=pd.read_csv(a.flora,dtype=str,usecols=["island_id","accepted_species","origin_status"])
    traits=direct_states(pd.read_csv(a.ledger,dtype=str))
    geo=pd.read_csv(a.geography,dtype={"island_id":str,"spatial_block":str})
    audit=pd.read_csv(a.state_audit,dtype=str)
    cfg=yaml.safe_load(a.syndrome_config.read_text())
    scores=build_species_score(audit,cfg["syndromes"]["selfing_core"],minimum_informative_traits=int(cfg["score_definition"]["minimum_informative_traits"]),require_all=False)
    result=run(flora,traits,geo,scores,parent)
    assert len(result["models"])==144
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("GENUS_SUPPORT "+json.dumps(result["support"]))
    for row in result["models"]:
        print("GENUS_MODEL "+json.dumps(row,allow_nan=False))
    print("GENUS_COMPLETE "+json.dumps({k:result[k] for k in ("original_slots","n_reproduced_parent_fits","max_parent_beta_difference")}))


if __name__=="__main__":
    main()
