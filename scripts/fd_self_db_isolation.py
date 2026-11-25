"""Exploratory FD extension of Chapter 1, not a submission selector.

Own GBIF occurrence/status roster, species-Direct states, corrected geography.
No row-level outputs. Q uses the exact expected with-replacement statistic in
a fixed-size draw WITHOUT replacement; excess MI uses exact conditional
fixed-margin expectation within each seeded rarefaction draw. Method adapted
from fd PR14/15, with the architecture vector kept intact.
"""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, log
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t
import yaml

from island_v2.chapter1_h2_complete_selfing_sensitivity import build_species_score

C, F, S, D = "flower_primary_color", "floral_form", "floral_symmetry", "tube_depth_class"
ENDPOINTS = {
    "triple_cfs_n30": ((C, F, S), 30),
    "colour_form_n5": ((C, F), 5),
    "colour_symmetry_n5": ((C, S), 5),
    "form_symmetry_n5": ((F, S), 5),
    "symmetry_depth_n5": ((S, D), 5),
    "form_depth_n5": ((F, D), 5),
}
REGIONS = ("northern_midlatitude", "northern_high_latitude", "tropical", "southern_extratropical")
SCOPES = ("all_recorded", "strict_native")
METRICS = ("marginal_q", "joint_q", "excess_mi")
DIST = "log_distance_to_continent_km"
COV = ("log_island_area_km2", "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4")
VARIANTS = ("baseline", "coverage_adjusted", "h2_matched_base", "h2_matched_adjusted",
            "h2_same_carrier_base", "h2_same_carrier_adjusted")
SEED = 20261125
REPEATS = 30
LOCK = {
    "flora": "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
    "ledger": "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
    "geography": "08218617b5880153986012e79817de65a74dd4e44b5a65e8a537dab8cb2cd0a7",
}


def digest(path):
    h = sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def direct_states(ledger):
    x = ledger.loc[
        ledger.trait_name.isin({C, F, S, D})
        & ledger.evidence_scope.eq("species_direct")
        & ledger.evidence_quality.isin(("high", "medium"))
        & ledger.normalized_value.notna()
        & ~ledger.normalized_value.astype(str).str.contains("|", regex=False)
    ].copy()
    if x.duplicated(["accepted_species", "trait_name"]).any():
        raise ValueError("duplicate Direct species-trait evidence")
    return x.pivot(index="accepted_species", columns="trait_name", values="normalized_value")


def expected_q(counts, n):
    """Exact E[1-sum (sample category proportions)^2], equal-size rarefaction."""
    counts = np.asarray(list(counts), dtype=float)
    N = counts.sum()
    if not 2 <= n <= N:
        raise ValueError("rarefaction support insufficient")
    return float((n - 1) / n * N / (N - 1) * (1 - np.square(counts / N).sum()))


@lru_cache(maxsize=100000)
def cell_expectation(n, a, b):
    den = comb(n, b)
    return sum(comb(a, k) * comb(n-a, b-k) / den * (k/n) * log(k*n/(a*b))
               for k in range(max(1, a+b-n), min(a,b)+1))


def excess_mi(sample):
    """First trait vs intact remaining trait tuple; exact fixed margins."""
    n = len(sample)
    a = Counter(r[0] for r in sample)
    b = Counter(tuple(r[1:]) for r in sample)
    cells = Counter((r[0], tuple(r[1:])) for r in sample)
    obs = sum((v/n) * log(v*n/(a[c]*b[z])) for (c,z),v in cells.items())
    exp = sum(cell_expectation(n, u, v) for u in a.values() for v in b.values())
    return float(obs-exp)


def morphology(arr, n, seed, repeats=REPEATS):
    N = len(arr)
    if N < n:
        raise ValueError("insufficient carriers")
    marginal = np.mean([expected_q(Counter(arr[:, j]).values(), n) for j in range(arr.shape[1])])
    joint = expected_q(Counter(map(tuple, arr)).values(), n)
    rng = np.random.default_rng(seed)
    if N == n:
        mi = excess_mi(arr)
    else:
        mi = np.mean([excess_mi(arr[rng.choice(N, size=n, replace=False)]) for _ in range(repeats)])
    return {"marginal_q": float(marginal), "joint_q": float(joint), "excess_mi": float(mi)}


def seed_for(island, endpoint, scope):
    return int.from_bytes(sha256(f"{SEED}:{island}:{endpoint}:{scope}".encode()).digest()[:8], "big")


def fit_cluster(part, response, predictors):
    cols = [response, *predictors, "spatial_block"]
    work = part[cols].replace([np.inf, -np.inf], np.nan).dropna()
    work = work.loc[work.spatial_block.astype(str).str.strip().ne("")]
    n = len(work)
    g = work.spatial_block.nunique()
    out = {"n": int(n), "blocks": int(g), "status": "not_testable"}
    if n < 50 or g < 10:
        return out
    X0 = work[list(predictors)].astype(float).to_numpy()
    sd = X0.std(axis=0)
    if np.any(sd <= 1e-12):
        return {**out, "status": "constant_predictor"}
    X = np.column_stack([np.ones(n), (X0-X0.mean(axis=0))/sd])
    k = X.shape[1]
    if np.linalg.matrix_rank(X) != k or n <= k or g <= k:
        return {**out, "status": "rank_or_df_failure"}
    y = work[response].astype(float).to_numpy()
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    residual = y-X@b
    inv = np.linalg.inv(X.T@X)
    labels = work.spatial_block.to_numpy()
    meat = np.zeros((k,k))
    for label in pd.unique(labels):
        v = X[labels == label].T @ residual[labels == label]
        meat += np.outer(v, v)
    vc = g/(g-1)*(n-1)/(n-k)*inv@meat@inv
    j = list(predictors).index(DIST)+1
    se = float(np.sqrt(max(0, vc[j,j])))
    if not np.isfinite(se) or se <= 1e-14:
        return {**out, "status": "degenerate_se"}
    beta = float(b[j])
    crit = float(t.ppf(.975,g-1))
    return {**out, "status":"fit", "beta":beta, "se":se, "lo":beta-crit*se,
            "hi":beta+crit*se, "p":float(2*t.sf(abs(beta/se),g-1)),
            "distance_sd_raw":float(sd[j-1])}


def adjust_fixed(rows):
    """Within scope, variant, and endpoint family; retain untestable slots as p=1."""
    for key in {(r["scope"],r["variant"],r["family"]) for r in rows}:
        ix = [i for i,r in enumerate(rows) if (r["scope"],r["variant"],r["family"]) == key]
        ordered = sorted(ix,key=lambda i: rows[i].get("p",1.0))
        m = len(ix)
        harmonic = sum(1/j for j in range(1,m+1))
        prev = 1.
        for rank,i in reversed(list(enumerate(ordered,1))):
            prev = min(prev, rows[i].get("p",1.)*m/rank)
            if rows[i]["status"] == "fit":
                rows[i]["q_bh"] = float(prev)
                rows[i]["q_by"] = float(min(1,prev*harmonic))
            rows[i]["family_slots"] = m


def models(flora, traits, geo, scores):
    if flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("duplicate island/species")
    if flora[["island_id","accepted_species","origin_status"]].isna().any().any():
        raise ValueError("missing flora identity/status")
    if not flora.origin_status.isin(["native","introduced","unresolved","endemic"]).all():
        raise ValueError("unknown origin status")
    if geo.island_id.duplicated().any():
        raise ValueError("duplicate geography")
    geo = geo[["island_id","analysis_regime","spatial_block",DIST,*COV]].copy()
    for c in (DIST,*COV):
        geo[c] = pd.to_numeric(geo[c],errors="coerce")
    score_map = scores.set_index("accepted_species").score.dropna()
    out = []
    support = []
    for scope in SCOPES:
        f = flora if scope == "all_recorded" else flora.loc[flora.origin_status.eq("native")]
        totals = f.groupby("island_id").size()
        rep = f[["island_id","accepted_species"]].copy()
        rep["selfing"] = rep.accepted_species.map(score_map)
        rep = rep.groupby("island_id").selfing.agg(["mean","count"]).rename(columns={"mean":"selfing_core","count":"selfing_n"})
        for endpoint,(columns,nmin) in ENDPOINTS.items():
            complete = traits.dropna(subset=list(columns))[list(columns)]
            joined = f[["island_id","accepted_species"]].merge(complete,left_on="accepted_species",right_index=True,validate="many_to_one")
            records=[]
            for island,part in joined.groupby("island_id",sort=True):
                if len(part) < nmin:
                    continue
                part = part.sort_values("accepted_species")
                row=morphology(part[list(columns)].to_numpy(str),nmin,seed_for(island,endpoint,scope))
                shared = part.accepted_species.map(score_map).dropna()
                records.append({"island_id":island,**row,"carrier_n":len(part),
                                "log_carriers":log(len(part)),"coverage":len(part)/totals[island],
                                "same_carrier_selfing":float(shared.mean()) if len(shared) else np.nan,
                                "same_carrier_n":len(shared)})
            private=pd.DataFrame(records,columns=["island_id",*METRICS,"carrier_n","log_carriers","coverage","same_carrier_selfing","same_carrier_n"])
            if scope=="all_recorded" and endpoint=="triple_cfs_n30" and (len(complete)!=2234 or len(private)!=928):
                raise ValueError(f"Original FD triple source drift: {len(complete)} species/{len(private)} islands")
            pre=len(private)
            private=private.merge(geo,on="island_id",how="inner",validate="one_to_one").merge(rep,left_on="island_id",right_index=True,how="left",validate="one_to_one")
            support.append({"scope":scope,"endpoint":endpoint,"direct_species":len(complete),"eligible_pre_geography":pre,"matched_geography":len(private)})
            for region in REGIONS:
                region_rows=private.loc[private.analysis_regime.eq(region)].copy()
                base=region_rows.dropna(subset=[DIST,*COV,"spatial_block",*METRICS])
                coverage=base.dropna(subset=["log_carriers","coverage"])
                h2=coverage.loc[coverage.selfing_n.ge(10)].dropna(subset=["selfing_core"])
                same=coverage.loc[coverage.same_carrier_n.ge(10)].dropna(subset=["same_carrier_selfing"])
                baseline=[DIST,*COV]
                control=baseline+["log_carriers","coverage"]
                designs = {
                    "baseline":(base,baseline),
                    "coverage_adjusted":(coverage,control),
                    "h2_matched_base":(h2,control),
                    "h2_matched_adjusted":(h2,control+["selfing_core"]),
                    "h2_same_carrier_base":(same,control),
                    "h2_same_carrier_adjusted":(same,control+["same_carrier_selfing"]),
                }
                for variant in VARIANTS:
                    # H2 is restricted in advance to the triple cohort.
                    if endpoint!="triple_cfs_n30" and variant.startswith("h2_"):
                        continue
                    sample,preds=designs[variant]
                    for metric in METRICS:
                        out.append({"scope":scope,"endpoint":endpoint,"region":region,
                                    "family":"triple" if nmin==30 else "pair",
                                    "variant":variant,"metric":metric,**fit_cluster(sample,metric,preds)})
    adjust_fixed(out)
    return {"status":"EXPLORATORY_DIRECT_FD_ISOLATION_NOT_CAUSAL",
            "n_original_flora_rows":len(flora),"n_original_flora_islands":flora.island_id.nunique(),
            "n_geography":len(geo),"support":support,"models":out,
            "seed":SEED,"mi_repeats":REPEATS,"n_selfing_score_species":len(score_map),
            "uncertainty":"CR1 spatial cluster t(G-1); repeated taxa and shared histories not fully modeled",
            "not_independent_replication":True,"no_row_level_export":True}


def main():
    p=argparse.ArgumentParser()
    for key in ("flora","ledger","geography","state-audit","syndrome-config","output"):
        p.add_argument("--"+key,required=True,type=Path)
    a=p.parse_args()
    for k,expected in LOCK.items():
        if digest(getattr(a,k))!=expected:
            raise ValueError("locked source SHA mismatch: "+k)
    flora=pd.read_csv(a.flora,dtype=str,usecols=["island_id","accepted_species","origin_status"])
    if flora.island_id.nunique()!=4505:
        raise ValueError("canonical full observed island denominator changed")
    traits=direct_states(pd.read_csv(a.ledger,dtype=str))
    geo=pd.read_csv(a.geography,dtype={"island_id":str,"spatial_block":str})
    if len(geo)!=8264:
        raise ValueError("corrected geography universe changed")
    audit=pd.read_csv(a.state_audit,dtype=str)
    cfg=yaml.safe_load(a.syndrome_config.read_text())
    scores=build_species_score(audit,cfg["syndromes"]["selfing_core"],
        minimum_informative_traits=int(cfg["score_definition"]["minimum_informative_traits"]),require_all=False)
    result=models(flora,traits,geo,scores)
    result["source_sha256"]={**LOCK,"state_audit":digest(a.state_audit),"syndrome_config":digest(a.syndrome_config)}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("FD_SUPPORT "+json.dumps(result["support"]))
    for row in result["models"]:
        if row["endpoint"]=="triple_cfs_n30" or (row["scope"]=="all_recorded" and row["variant"]=="baseline"):
            print("FD_MODEL "+json.dumps(row,allow_nan=False))
    print("FD_COMPLETE",len(result["models"]),"anonymous aggregate model slots")


if __name__=="__main__":
    main()
