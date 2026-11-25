"""Geological-origin-stratified binary profile-margin tests, including origin contrasts.

Directly tests oceanic-versus-continental differences in observed-minus-null
functional metrics WITHIN biogeographic region, rather than juxtaposing
significant/nonsignificant results from separate strata. Exploratory and
source-unmatched: profile binary margins alone cannot identify selection.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from .core import ProfileSpace, csv_rows, load_profiles, load_carriers
from .g1_fixed_margins import (
    REGIMES, METRICS, _curveball_once, _margins,
    _local_q, _pair_nn, _mixing, _p_value, _bh,
)
from .pipeline import SCOPES, REF_INPUT, _write_csv, file_sha256

TYPES=("continental", "oceanic")
STATUSES=("ok", "no_complete_islands", "insufficient_islands", "insufficient_pairs", "mixing_failed")
GROUP_COLUMNS=(
    "scope", "regime", "origin", "metric", "status", "n_islands", "n_pairs_all",
    "n_local_scored", "n_pairs_scored", "observed", "null_median",
    "delta", "rhat", "max_abs_lag1", "n_distinct", "p_approx",
)
CONTRAST_COLUMNS=(
    "scope", "regime", "metric", "status", "n_continental", "n_oceanic",
    "n_pairs_continental", "n_pairs_oceanic", "observed_difference_oceanic_minus_continental",
    "null_median_difference", "excess_difference", "null_q025", "null_q975",
    "p_approx", "q_bh_32", "q_by_32",
)


def _seed(seed: int, *terms: str):
    return int.from_bytes(hashlib.sha256((":".join([str(seed),*terms])).encode()).digest()[:8], "little")


def _subsample(elements, max_n, seed):
    ordered=sorted(elements)
    if len(ordered)<=max_n:
        return ordered
    rng=np.random.default_rng(seed)
    indexes=sorted(rng.choice(len(ordered),size=max_n,replace=False))
    return [ordered[int(i)] for i in indexes]


def run_origin_null(inputs: Path, g0: Path, crosswalk: Path, output: Path, *,
                    draws_per_chain: int = 200, max_islands: int = 36,
                    max_pairs: int = 36, min_islands: int = 20,
                    min_pairs: int = 20, seed: int = 20261108) -> dict:
    if draws_per_chain < 20 or min_islands < 2 or min_pairs < 1:
        raise ValueError("MCMC support below permitted screening threshold")
    inputs, g0, crosswalk, output=map(Path,(inputs,g0,crosswalk,output))
    output.mkdir(parents=True,exist_ok=True)
    classifications = {r["island_id"]:r for r in csv_rows(crosswalk)}
    space=ProfileSpace(load_profiles(inputs/REF_INPUT))
    groups=[]
    traces={}
    scored_counts={}
    for scope in SCOPES:
        carrier = load_carriers(
            inputs/"floral_space_20261106"/scope/"island_profile_carriers.csv.gz")
        metrics = {r["island_id"]:r for r in csv_rows(g0/f"{scope}_islands.csv")}
        pairs=list(csv_rows(g0/f"{scope}_pairs.csv"))
        for regime in REGIMES:
            for origin in TYPES:
                ids=sorted(i for i in metrics
                           if i in classifications
                           and metrics[i]["analysis_regime"] == regime
                           and classifications[i]["geological_origin"] == origin)
                pair_ids=sorted(set(tuple(sorted((r["island_a"],r["island_b"])))
                          for r in pairs if r["analysis_regime"] == regime
                          and r["island_a"] in ids and r["island_b"] in ids))
                if len(ids) < min_islands:
                    for metric in METRICS:
                        groups.append({"scope":scope,"regime":regime,"origin":origin,"metric":metric,
                                       "status":"insufficient_islands","n_islands":len(ids),
                                       "n_pairs_all":len(pair_ids), "n_local_scored":0,
                                       "n_pairs_scored":0})
                    continue
                local_ids=_subsample(ids,max_islands,_seed(seed,scope,regime,origin,"local"))
                selected_pairs=_subsample(pair_ids,max_pairs,
                                          _seed(seed,scope,regime,origin,"pairs"))
                index={identifier:i for i,identifier in enumerate(ids)}
                rows=[set(space.lookup[p] for p in carrier[identifier]) for identifier in ids]
                row_margins=_margins(rows)
                island_ix=[index[i] for i in local_ids]
                pair_ix=[(index[a],index[b]) for a,b in selected_pairs]
                observed={METRICS[0]:_local_q(rows,island_ix,space.matrix)}
                if len(pair_ix) >= min_pairs:
                    observed[METRICS[1]]=_pair_nn(rows,pair_ix,space.matrix)
                samples={metric:[] for metric in observed}
                for chain in range(2):
                    rng=np.random.default_rng(_seed(seed,scope,regime,origin,"chain",str(chain)))
                    current=[r.copy() for r in rows]
                    for _ in range(12*len(rows)):
                        _curveball_once(current,rng)
                    ts={metric:np.empty(draws_per_chain,dtype=float) for metric in observed}
                    for t in range(draws_per_chain):
                        for _ in range(3*len(rows)):
                            _curveball_once(current,rng)
                        ts[METRICS[0]][t]=_local_q(current,island_ix,space.matrix)
                        if METRICS[1] in observed:
                            ts[METRICS[1]][t]=_pair_nn(current,pair_ix,space.matrix)
                    if _margins(current)!=row_margins:
                        raise AssertionError("Curveball profile margins changed")
                    for metric in observed:
                        samples[metric].append(ts[metric])
                for metric in METRICS:
                    base={"scope":scope,"regime":regime,"origin":origin,
                          "metric":metric,"n_islands":len(ids),
                          "n_pairs_all":len(pair_ids),"n_local_scored":len(local_ids),
                          "n_pairs_scored":len(pair_ix)}
                    if metric not in observed:
                        groups.append({**base,"status":"insufficient_pairs"})
                        continue
                    values=np.concatenate(samples[metric])
                    rhat,acf,distinct=_mixing(samples[metric])
                    mixing=bool(np.isfinite(rhat) and rhat<=1.05 and acf<=0.30 and distinct>=8)
                    status="ok" if mixing else "mixing_failed"
                    groups.append({**base,"status":status,"observed":observed[metric],
                                   "null_median":float(np.median(values)),
                                   "delta":float(observed[metric]-np.median(values)),
                                   "rhat":rhat,"max_abs_lag1":acf,"n_distinct":distinct,
                                   "p_approx":_p_value(observed[metric],values) if mixing else None})
                    traces[(scope,regime,origin,metric)]=(observed[metric],values,mixing)
    contrasts=[]
    for scope in SCOPES:
        for regime in REGIMES:
            for metric in METRICS:
                a=traces.get((scope,regime,"continental",metric))
                b=traces.get((scope,regime,"oceanic",metric))
                row_a=next(r for r in groups if r["scope"]==scope and r["regime"]==regime
                            and r["metric"]==metric and r["origin"]=="continental")
                row_b=next(r for r in groups if r["scope"]==scope and r["regime"]==regime
                            and r["metric"]==metric and r["origin"]=="oceanic")
                base={"scope":scope,"regime":regime,"metric":metric,
                      "n_continental":row_a["n_islands"],"n_oceanic":row_b["n_islands"],
                      "n_pairs_continental":row_a["n_pairs_all"],
                      "n_pairs_oceanic":row_b["n_pairs_all"]}
                if not a or not b or not a[2] or not b[2]:
                    contrasts.append({**base,"status":"insufficient_or_unmixed_component"})
                    continue
                obs=float(b[0]-a[0])
                null=b[1]-a[1] # matched-index independent chain draws
                ref=float(np.median(null))
                contrasts.append({**base,"status":"exploratory_binary_margin",
                                  "observed_difference_oceanic_minus_continental":obs,
                                  "null_median_difference":ref,
                                  "excess_difference":obs-ref,
                                  "null_q025":float(np.quantile(null,0.025)),
                                  "null_q975":float(np.quantile(null,0.975)),
                                  "p_approx":_p_value(obs,null)})
    if len(contrasts)!=32 or len(groups)!=64:
        raise AssertionError("Unexpected number of origin-null hypothesis slots")
    p=np.asarray([r.get("p_approx",1.0) if r.get("p_approx") is not None else 1.0
                  for r in contrasts])
    bh=_bh(p)
    H=float(sum(1/i for i in range(1,33)))
    for i,row in enumerate(contrasts):
        if row["status"]=="exploratory_binary_margin":
            row["q_bh_32"]=float(bh[i])
            row["q_by_32"]=float(min(1,bh[i]*H))
    _write_csv(output/"origin_g1_groups.csv",groups,GROUP_COLUMNS)
    _write_csv(output/"origin_g1_contrasts.csv",contrasts,CONTRAST_COLUMNS)
    audit={
        "status":"posthoc_exploratory_direct_origin_contrasts",
        "null":"fixed within-origin, within-regime island profile count and prevalence",
        "not_conditioned_on":["source pool","genus","shared taxa","origin-label selection",
                              "observation effort","island history"],
        "metric_weighting":"equal recorded profile (not species-carrier weighted)",
        "n_draws_per_group":2*draws_per_chain,"n_chains":2,
        "comparisons_family":32,"tested_contrasts":sum(r["status"]=="exploratory_binary_margin" for r in contrasts),
        "bh_supported":sum(r.get("q_bh_32",1)>=0 and r.get("q_bh_32",1)<0.05 for r in contrasts),
        "by_supported":sum(r.get("q_by_32",1)>=0 and r.get("q_by_32",1)<0.05 for r in contrasts),
        "notes":"MCMC approximate; overlapped evidence scopes and island-pair reuse; non-causal",
        "inputs_sha256":{"crosswalk":file_sha256(crosswalk),
                         "profile_reference":file_sha256(inputs/REF_INPUT)},
    }
    (output/"manifest.json").write_text(json.dumps(audit,indent=2)+"\n")
    print(json.dumps(audit))
    return audit


def main():
    ap=argparse.ArgumentParser()
    for name in ("input","g0","crosswalk","output"):
        ap.add_argument("--"+name,type=Path,required=True)
    ap.add_argument("--draws",type=int,default=200)
    args=ap.parse_args()
    run_origin_null(args.input,args.g0,args.crosswalk,args.output,
                    draws_per_chain=args.draws)


if __name__=="__main__":
    main()
