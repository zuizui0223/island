"""Species-name crossfit of whole floral profiles and corrected FD slopes.

All prediction distributions exclude the entire target name fold. Known genus
and observed incidence are oracle inputs, not predicted colonization. Source
data and every row-level prediction remain private. Predictive simulation
intervals are not permutation p values or independent biological uncertainty.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import json
from math import log
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from scripts.fd_self_db_isolation import (
    C,F,S,DIST,COV,LOCK,METRICS,REGIONS,SCOPES,VARIANTS,
    direct_states,digest,fit_cluster,adjust_fixed,
)
from scripts.fd_self_db_genus_null import genus_key,all_metrics,parent_draws,linear_weights
from island_v2.chapter1_h2_complete_selfing_sensitivity import build_species_score

N=30
DRAWS=199
SEEDS=(20261125,20261126)
TAUS=(1.,8.)
MODELS=("blind","genus_tau1","genus_tau8")
ENDPOINT="triple_cfs_n30"


def fold(name,seed):
    return int.from_bytes(sha256(f"{seed}:{name}".encode()).digest()[:8],"big")%5


def probabilities(names,codes,genera,k,seed):
    """Fixed global name folds, no target-fold labels even in the global prior.

    The outcome alphabet is the archived ontology/profile space. Its existence
    does not provide frequencies: tiny symmetric smoothing is identical for
    both references. No same-genus donors means the identical blind fallback.
    """
    folds=np.array([fold(s,seed) for s in names])
    codes=np.asarray(codes)
    out=np.zeros((3,len(names),k))
    support=np.zeros(len(names),dtype=int)
    genus=np.asarray(genera)
    for f in range(5):
        train=folds!=f
        if not train.any():
            raise ValueError("empty training fold")
        prior=np.bincount(codes[train],minlength=k).astype(float)+.02
        prior/=prior.sum()
        targets=np.flatnonzero(~train)
        by_genus={}
        for g in np.unique(genus[train]):
            counts=np.bincount(codes[train & (genus==g)],minlength=k)
            by_genus[g]=counts
        for i in targets:
            counts=by_genus.get(genus[i],np.zeros(k))
            n=int(counts.sum())
            support[i]=n
            out[0,i]=prior
            for j,tau in enumerate(TAUS,1):
                out[j,i]=(counts+tau*prior)/(n+tau)
    if not np.allclose(out.sum(2),1):
        raise ValueError("invalid categorical distributions")
    return out,support,folds


def simulate(p,uniform):
    # Shared random numbers across references and one draw per species globally
    # within the regional pool. No independent draws per island occurrence.
    cdf=np.cumsum(p,axis=1)
    return np.array([np.searchsorted(cdf[i],uniform[:,i],side="right")
                     for i in range(len(p))],dtype=np.int32).T.clip(0,p.shape[1]-1)


def proper_scores(p,truth):
    pr=p[np.arange(len(truth)),truth]
    return np.column_stack((1-2*pr+np.square(p).sum(1),-np.log(np.maximum(pr,1e-300))))


def grouped_gain_interval(gain,groups,seed):
    """Equal-group conditional bootstrap, not a refit or phylogeny bootstrap."""
    a=pd.DataFrame({"gain":gain,"g":groups}).groupby("g").gain.mean().to_numpy(float)
    if len(a)<10:
        return None
    rng=np.random.default_rng(seed)
    d=a[rng.integers(len(a),size=(1999,len(a)))].mean(1)
    return [float(a.mean()),float(np.quantile(d,.025)),float(np.quantile(d,.975))]


def run(flora,traits,geo,scores,parent):
    complete=traits.dropna(subset=[C,F,S])[[C,F,S]]
    if len(complete)!=2234 or flora.island_id.nunique()!=4505 or len(geo)!=8264:
        raise ValueError("source universe changed")
    if flora.duplicated(["island_id","accepted_species"]).any() or geo.island_id.duplicated().any():
        raise ValueError("duplicate source keys")
    geo=geo[["island_id","analysis_regime","spatial_block",DIST,*COV]].copy()
    for c in (DIST,*COV):
        geo[c]=pd.to_numeric(geo[c],errors="coerce")
    score=scores.set_index("accepted_species").score.dropna()
    parent_lookup={(r["scope"],r["region"],r["variant"],r["metric"]):r for r in parent["models"] if r["endpoint"]==ENDPOINT}
    rows=[]; skills=[]; support=[]; errors=[]
    # Scope-specific reference; names may occur across regions, but a name has
    # the same fold everywhere and cannot be used in its own regional training.
    profiles=np.unique(complete.to_numpy(str),axis=0)
    codebook={tuple(p):i for i,p in enumerate(profiles)}
    for scope in SCOPES:
        f=flora if scope=="all_recorded" else flora.loc[flora.origin_status.eq("native")]
        totals=f.groupby("island_id").size()
        rep=f[["island_id","accepted_species"]].copy()
        rep["score"]=rep.accepted_species.map(score)
        rep=rep.groupby("island_id").score.agg(["mean","count"]).rename(columns={"mean":"selfing_core","count":"selfing_n"})
        joined=f[["island_id","accepted_species"]].merge(complete,left_on="accepted_species",right_index=True,validate="many_to_one")
        if int((joined.groupby("island_id").size()>=N).sum())!=(928 if scope=="all_recorded" else 87):
            raise ValueError("triple support changed")
        joined=joined.merge(geo[["island_id","analysis_regime"]],on="island_id",validate="many_to_one")
        for region in REGIONS:
            reg=joined.loc[joined.analysis_regime.eq(region)]
            names=sorted(reg.accepted_species.unique())
            if len(names)<5:
                # No eligible target islands in this tiny scope; fail if that changes.
                if (reg.groupby("island_id").size()>=N).any():
                    raise ValueError("unexpected tiny reference with eligible islands")
                continue
            codes=np.array([codebook[tuple(p)] for p in complete.loc[names].to_numpy(str)])
            genera=[genus_key(s) for s in names]
            idx={s:i for i,s in enumerate(names)}
            eligible=[(i,p.sort_values("accepted_species")) for i,p in reg.groupby("island_id",sort=True) if len(p)>=N]
            target_ids=np.array(sorted({idx[s] for _,p in eligible for s in p.accepted_species}),dtype=int)
            for seed in SEEDS:
                pred,ndonor,folds=probabilities(names,codes,genera,len(profiles),seed)
                support.append({"scope":scope,"region":region,"seed":seed,"pool_species":len(names),"eligible_islands":len(eligible),
                    "unique_target_species":len(target_ids),"target_species_with_congener":int(np.sum(ndonor[target_ids]>0)),
                    "target_species_with_2plus_congeners":int(np.sum(ndonor[target_ids]>=2))})
                if not eligible:
                    continue
                scores_by_model=[proper_scores(p[target_ids],codes[target_ids]) for p in pred]
                for j,model in enumerate(MODELS):
                    for k,metric in enumerate(("brier","logloss")):
                        obs=scores_by_model[j][:,k]
                        gain=scores_by_model[0][:,k]-obs
                        valid=ndonor[target_ids]>0
                        skills.append({"scope":scope,"region":region,"seed":seed,"model":model,"metric":metric,
                            "n_target_species":len(target_ids),"loss":float(obs.mean()),
                            "gain_vs_blind":float(gain.mean()),
                            "gain_supported_species":float(gain[valid].mean()) if valid.any() else None,
                            "genus_equal_gain_interval":grouped_gain_interval(gain,[genera[i] for i in target_ids],seed+100+k)})
                rng=np.random.default_rng(int.from_bytes(sha256(f"{seed}:{scope}:{region}".encode()).digest()[:8],"big"))
                u=rng.random((DRAWS,len(names)))
                assignments=np.vstack([codes[None,:],*[simulate(p,u) for p in pred]])
                rec=[]; cubes=[]
                for island,p in eligible:
                    ids=np.array([idx[s] for s in p.accepted_species],dtype=int)
                    z=all_metrics(assignments[:,ids],profiles,parent_draws(len(p),island,scope))
                    shared=p.accepted_species.map(score).dropna()
                    rec.append({"island_id":island,"log_carriers":log(len(p)),"coverage":len(p)/totals[island],
                        "same_carrier_selfing":float(shared.mean()) if len(shared) else np.nan,"same_carrier_n":len(shared),
                        "donor_coverage":float((ndonor[ids]>0).mean()),**{m:float(z[0,j]) for j,m in enumerate(METRICS)}})
                    cubes.append(z)
                priv=pd.DataFrame(rec).merge(geo,on="island_id",how="left",validate="one_to_one").merge(rep,left_on="island_id",right_index=True,how="left",validate="one_to_one")
                cube=np.stack(cubes)
                base=priv.dropna(subset=[DIST,*COV,"spatial_block",*METRICS])
                coverage=base.dropna(subset=["log_carriers","coverage"])
                h2=coverage.loc[coverage.selfing_n.ge(10)].dropna(subset=["selfing_core"])
                same=coverage.loc[coverage.same_carrier_n.ge(10)].dropna(subset=["same_carrier_selfing"])
                baseline=[DIST,*COV]; control=baseline+["log_carriers","coverage"]
                designs={"baseline":(base,baseline),"coverage_adjusted":(coverage,control),
                    "h2_matched_base":(h2,control),"h2_matched_adjusted":(h2,control+["selfing_core"]),
                    "h2_same_carrier_base":(same,control),"h2_same_carrier_adjusted":(same,control+["same_carrier_selfing"])}
                for variant,(part,predictors) in designs.items():
                    for j,metric in enumerate(METRICS):
                        orig=fit_cluster(part,metric,predictors)
                        prev=parent_lookup[(scope,region,variant,metric)]
                        if any(orig[k]!=prev[k] for k in ("status","n","blocks")):
                            raise ValueError("parent support/status mismatch")
                        if orig["status"]!="fit":
                            continue
                        delta=abs(orig["beta"]-prev["beta"])
                        if delta>1e-10:
                            raise ValueError("parent observed slope mismatch")
                        errors.append(delta)
                        work=part[[metric,*predictors,"spatial_block"]].replace([np.inf,-np.inf],np.nan).dropna()
                        work=work.loc[work.spatial_block.astype(str).str.strip().ne("")]
                        a=cube[work.index.to_numpy(),:,j]
                        w=linear_weights(work,predictors)
                        slopes=w@a
                        if abs(slopes[0]-orig["beta"])>1e-10:
                            raise ValueError("observed slope misalignment")
                        y=a[:,0]
                        blind=a[:,1:1+DRAWS].mean(1)
                        for k,model in enumerate(MODELS):
                            start=1+k*DRAWS
                            draws=slopes[start:start+DRAWS]
                            means=a[:,start:start+DRAWS].mean(1)
                            gain=np.square(y-blind)-np.square(y-means)
                            re=work.copy()
                            re["residual"]=y-means
                            rf=fit_cluster(re,"residual",predictors)
                            row={"scope":scope,"region":region,"seed":seed,"variant":variant,"metric":metric,"model":model,
                                "n":orig["n"],"blocks":orig["blocks"],"observed_beta":orig["beta"],
                                "predicted_beta":float(draws.mean()),"predicted_beta_draw_lo":float(np.quantile(draws,.025)),
                                "predicted_beta_draw_hi":float(np.quantile(draws,.975)),
                                "residual_beta":float(orig["beta"]-draws.mean()),
                                "island_mse":float(np.mean(np.square(y-means))),
                                "island_mse_gain_vs_blind":float(gain.mean()),
                                "block_equal_gain_interval":grouped_gain_interval(gain,work.spatial_block,seed+200+j),
                                "mean_donor_coverage":float(priv.loc[work.index,"donor_coverage"].mean()),
                                "residual_status":rf["status"]}
                            if rf["status"]=="fit":
                                if abs(rf["beta"]-row["residual_beta"])>1e-10:
                                    raise ValueError("residual slope identity failed")
                                row.update(residual_lo=rf["lo"],residual_hi=rf["hi"],residual_p=rf["p"])
                            rows.append(row)
    if len(errors)!=108:
        raise ValueError("54 parent fits x 2 seeds not reproduced")
    return {"status":"SPECIES_NAME_CROSSFIT_ORACLE_GENUS_NOT_CAUSAL","draws_per_reference":DRAWS,
        "seeds":list(SEEDS),"tau":list(TAUS),"source_sha256":LOCK,"support":support,"proper_scores":skills,
        "models":rows,"n_reproduced_parent_fits":len(errors),"max_parent_beta_difference":max(errors),
        "no_row_level_export":True,"not_geographic_holdout":True}


def main():
    p=argparse.ArgumentParser()
    for key in ("flora","ledger","geography","state-audit","syndrome-config","parent","output"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    for k,v in LOCK.items():
        if digest(getattr(a,k))!=v:
            raise ValueError("source SHA mismatch: "+k)
    if digest(a.parent)!="d06af1e649bd7070f251927dfafe1000ac08f951355c4bf514cf87950a692871":
        raise ValueError("parent result changed")
    parent=json.loads(a.parent.read_text())
    for k,path in (("state_audit",a.state_audit),("syndrome_config",a.syndrome_config)):
        if digest(path)!=parent["source_sha256"][k]:
            raise ValueError("parent selfing source changed")
    flora=pd.read_csv(a.flora,dtype=str,usecols=["island_id","accepted_species","origin_status"])
    traits=direct_states(pd.read_csv(a.ledger,dtype=str))
    geo=pd.read_csv(a.geography,dtype={"island_id":str,"spatial_block":str})
    audit=pd.read_csv(a.state_audit,dtype=str)
    cfg=yaml.safe_load(a.syndrome_config.read_text())
    scores=build_species_score(audit,cfg["syndromes"]["selfing_core"],minimum_informative_traits=int(cfg["score_definition"]["minimum_informative_traits"]),require_all=False)
    d=run(flora,traits,geo,scores,parent)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
    print("CROSS_SUPPORT "+json.dumps(d["support"]))
    for r in d["proper_scores"]:
        print("CROSS_SKILL "+json.dumps(r))
    for r in d["models"]:
        print("CROSS_MODEL "+json.dumps(r))
    print("CROSS_COMPLETE "+json.dumps({k:d[k] for k in ("n_reproduced_parent_fits","max_parent_beta_difference")}))


if __name__=="__main__":
    main()
