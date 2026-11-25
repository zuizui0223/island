"""One frozen tropical island floral-sector x independent guild-map test.

Requires PRIVATE island/source joins. Never exports site/taxon/source identifiers.
This is not a direct pollinator visitation, fitness or adaptation test.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import f as f_distribution

POPULATION = 4505
PRIMARY_REGION = "tropical"
PRIMARY_SCOPE = "all_recorded"
PRIMARY_SOURCE_MODE = "geo_k5"
PRIMARY_SECTOR = "actinomorphic_deep_tube"
MIN_ISLANDS = 60
MIN_CLUSTERS = 15
COVARIATES = (
    "log1p_distance_km", "log1p_area_km2", "absolute_latitude",
    "temperature_pc1", "precipitation_pc1",
    "genus_reference_expected_delta",
    "island_pair_coverage", "source_pair_coverage",
)
REQUIRED = (
    "island_id", "region", "scope", "source_mode", "sector",
    "spatial_block", "first_source", "island_total", "source_total",
    "island_pair_n", "source_pair_n", "island_pair_positive",
    "source_pair_positive", "bird_island_richness", "bird_source_richness",
    "sphingid_island_richness", "sphingid_source_richness",
    "bird_range_qualified", "sphingid_model_qualified",
    "genus_reference_qualified", "genus_reference_expected_delta",
    "log1p_distance_km", "log1p_area_km2", "absolute_latitude",
    "temperature_pc1", "precipitation_pc1",
)
NUMERIC = (
    "island_total", "source_total", "island_pair_n", "source_pair_n",
    "island_pair_positive", "source_pair_positive",
    "bird_island_richness", "bird_source_richness",
    "sphingid_island_richness", "sphingid_source_richness",
    "log1p_distance_km", "log1p_area_km2", "absolute_latitude",
    "temperature_pc1", "precipitation_pc1",
    "genus_reference_expected_delta",
)


def _select(df: pd.DataFrame, require_population: bool) -> tuple[pd.DataFrame, dict]:
    missing = sorted(set(REQUIRED) - set(df.columns))
    if missing:
        raise ValueError(f"missing original source/mapping columns: {missing}")
    if df.island_id.isna().any() or df.island_id.duplicated().any():
        raise ValueError("private input must have one non-null original island_id per row")
    if require_population and df.island_id.nunique() != POPULATION:
        raise ValueError("source population not equal to 4505 original GIFT floras")
    frame = df.loc[
        df.region.eq(PRIMARY_REGION) & df.scope.eq(PRIMARY_SCOPE)
        & df.sector.eq(PRIMARY_SECTOR) & df.source_mode.eq(PRIMARY_SOURCE_MODE)
    ].copy()
    for key in NUMERIC:
        frame[key] = pd.to_numeric(frame[key], errors="coerce")
    for pos, denom in (("island_pair_positive", "island_pair_n"),
                       ("source_pair_positive", "source_pair_n"),
                       ("island_pair_n", "island_total"),
                       ("source_pair_n", "source_total")):
        bad = (frame[pos].notna() & frame[denom].notna()
               & ((frame[pos] < 0) | (frame[denom] < frame[pos])))
        if bad.any():
            raise ValueError(f"invalid nested original recorded counts: {pos} / {denom}")
    for key in ("bird_island_richness", "bird_source_richness",
                "sphingid_island_richness", "sphingid_source_richness"):
        if frame[key].lt(0).any():
            raise ValueError(f"negative modeled guild richness: {key}")
    for key in ("bird_range_qualified", "sphingid_model_qualified",
                "genus_reference_qualified"):
        if not frame[key].dropna().isin([True, False, 0, 1]).all():
            raise ValueError(f"invalid qualified evidence flag: {key}")
    if frame.genus_reference_expected_delta.abs().gt(1 + 1e-9).any():
        raise ValueError("genus expected island-minus-source fraction outside [-1,1]")
    audit = {
        "n_source_parent_floras": int(df.island_id.nunique()),
        "parent_complete_4505": bool(df.island_id.nunique() == POPULATION),
        "n_tropical_candidate": int(len(frame)),
        "n_independent_name_disjoint_genus_reference": int(
            frame.genus_reference_qualified.eq(True).sum()),
    }
    frame = frame.loc[
        frame.island_pair_n.ge(5) & frame.source_pair_n.ge(5)
        & frame.island_total.gt(0) & frame.source_total.gt(0)
        & frame.bird_range_qualified.eq(True) & frame.sphingid_model_qualified.eq(True)
        & frame.genus_reference_qualified.eq(True)
    ].dropna(subset=[*NUMERIC, "spatial_block", "first_source"])
    audit["n_map_and_direct_qualified"] = int(len(frame))
    audit["n_10deg_blocks"] = int(frame.spatial_block.nunique())
    audit["n_first_source_groups"] = int(frame.first_source.nunique())
    if not frame.empty:
        frame["island_pair_coverage"] = frame.island_pair_n / frame.island_total
        frame["source_pair_coverage"] = frame.source_pair_n / frame.source_total
        frame["response"] = (frame.island_pair_positive / frame.island_pair_n
                              - frame.source_pair_positive / frame.source_pair_n)
        # Both guild components are source-relative, equal scale, and derived
        # exclusively from independently constructed range/SDM maps.
        for name, label in (("bird", "bird"), ("sphingid", "sphingid")):
            delta = (np.log1p(frame[f"{label}_island_richness"])
                     - np.log1p(frame[f"{label}_source_richness"]))
            sd = float(delta.std(ddof=0))
            if sd < 1e-10 or not np.isfinite(sd):
                audit[f"{name}_guild_contrast_has_variation"] = False
            else:
                frame[f"{name}_contrast_z"] = (delta - float(delta.mean())) / sd
        # Conservative binary-sector bounds: all unknown Direct assignments
        # could be positive or negative. These bounds do not assume MAR.
        frame["delta_lower"] = (frame.island_pair_positive / frame.island_total
            - (frame.source_pair_positive + frame.source_total
               - frame.source_pair_n) / frame.source_total)
        frame["delta_upper"] = ((frame.island_pair_positive + frame.island_total
              - frame.island_pair_n) / frame.island_total
            - frame.source_pair_positive / frame.source_total)
        audit["median_island_pair_coverage"] = float(frame.island_pair_coverage.median())
        audit["median_mainland_pair_coverage"] = float(frame.source_pair_coverage.median())
    return frame, audit


def _cluster_fit(x: np.ndarray, y: np.ndarray, labels: np.ndarray,
                 draws: int, seed: int) -> dict:
    """Block-equal regression, CR1 joint Wald F, group bootstrap intervals."""
    group_ids = sorted(pd.unique(labels).tolist())
    n, p = x.shape
    g = len(group_ids)
    moments = []
    rhs = []
    for label in group_ids:
        select = np.flatnonzero(labels == label)
        z = x[select]
        q = y[select]
        moments.append(z.T @ z / len(select))
        rhs.append(z.T @ q / len(select))
    A = np.stack(moments)
    b = np.stack(rhs)
    gram = A.sum(axis=0)
    if np.linalg.matrix_rank(gram, tol=1e-9) < p:
        return {"status": "RANK_DEFICIENT"}
    beta = np.linalg.solve(gram, b.sum(axis=0))
    score = b - np.einsum("gij,j->gi", A, beta)
    inv = np.linalg.inv(gram)
    cov = inv @ (score.T @ score) @ inv * (g / (g - 1)) * ((n - 1) / (n - p))
    v = cov[1:3, 1:3]
    if np.linalg.matrix_rank(v) < 2 or not np.isfinite(v).all():
        return {"status": "GUILD_COVARIANCE_UNIDENTIFIABLE"}
    f_value = float(beta[1:3] @ np.linalg.solve(v, beta[1:3]) / 2)
    p_joint = float(f_distribution.sf(f_value, 2, g - 1))
    rng = np.random.default_rng(seed)
    draws_ok = []
    for _ in range(draws):
        ix = rng.integers(0, g, g)
        a = A[ix].sum(axis=0)
        if np.linalg.matrix_rank(a, tol=1e-9) < p:
            continue
        draws_ok.append(np.linalg.solve(a, b[ix].sum(axis=0))[1:3])
    if len(draws_ok) < max(50, int(.8 * draws)):
        return {"status": "BOOTSTRAP_RANK_UNSTABLE", "n_valid_draws":len(draws_ok)}
    ci = np.quantile(np.asarray(draws_ok), [.025, .975], axis=0)
    return {"status":"ESTIMABLE", "n_groups":g,
            "bird_beta_z":float(beta[1]), "sphingid_beta_z":float(beta[2]),
            "bird_ci95":[float(ci[0,0]),float(ci[1,0])],
            "sphingid_ci95":[float(ci[0,1]),float(ci[1,1])],
            "joint_wald_F":f_value, "joint_p_exploratory":p_joint,
            "n_valid_draws":len(draws_ok)}


def run(df: pd.DataFrame, *, require_population: bool=True,
        bootstrap: int=499, seed: int=20261119) -> dict:
    if bootstrap < 100:
        raise ValueError("bootstrap >= 100 draws required")
    frame, audit = _select(df, require_population)
    out = {"contract":"fd_tropical_guild_flower_alignment_v1",
           "role":"exploratory_frozen_after_prior_descriptive_results_not_causal",
           "region":PRIMARY_REGION, "source_mode":PRIMARY_SOURCE_MODE,
           "scope":PRIMARY_SCOPE, "same_species_pair_sector":PRIMARY_SECTOR,
           "parent_count_matches_4505":bool(audit["parent_complete_4505"]),
           "support":audit,
           "conditions_on_external_name_disjoint_genus_expected_delta":True,
           "scientific_limit":"correlation of mapped potential guild availability with recorded flora; not visitation, adaptation, colonization or fitness"}
    if (len(frame) < MIN_ISLANDS or
        frame.spatial_block.nunique() < MIN_CLUSTERS or
        frame.first_source.nunique() < MIN_CLUSTERS or
        not {"bird_contrast_z", "sphingid_contrast_z"}.issubset(frame.columns)):
        out["status"] = "NOT_IDENTIFIABLE_GUILD_MAP_OR_CLUSTER_SUPPORT"
        return out
    fields = ["bird_contrast_z", "sphingid_contrast_z", *COVARIATES]
    X = np.column_stack([np.ones(len(frame)), frame[fields].to_numpy(float)])
    y = frame.response.to_numpy(float)
    if not np.isfinite(X).all() or not np.isfinite(y).all():
        raise ValueError("nonfinite original predictors or outcomes")
    # Report uncertainty twice, from physical blocks and independent
    # candidate-source groups. Do not misstate either as experimental replicates.
    fits = {}
    for ix, key in enumerate(("spatial_block", "first_source")):
        fits[key] = _cluster_fit(X, y, frame[key].to_numpy(), bootstrap, seed+ix)
    out["fits"] = fits
    out["status"] = ("EXPLORATORY_GUILD_ASSOCIATION_ESTIMABLE"
                     if all(v["status"]=="ESTIMABLE" for v in fits.values())
                     else "NOT_IDENTIFIABLE_CLUSTER_MODEL")
    if out["status"].startswith("EXPLORATORY"):
        # One preselected 2-df joint test, no post-hoc choice of guild.
        directional = all(
            v["bird_beta_z"] >= 0 and v["sphingid_beta_z"] >= 0
            and v["joint_p_exploratory"] < .05
            for v in fits.values())
        out["predeclared_joint_positive_gate"] = bool(directional and any(
            all(fits[g][f"{t}_ci95"][0]>0 for g in fits)
            for t in ("bird", "sphingid")))
        out["missingness_bounds_on_mean_delta"] = {
            "lower_unweighted":float(frame.delta_lower.mean()),
            "upper_unweighted":float(frame.delta_upper.mean()),
            "sign_robust_to_arbitrary_missing_binary_states":bool(
                frame.delta_lower.mean() > 0 or frame.delta_upper.mean() < 0),
        }
        out["missingness_warning"] = ("Coverage covariates are not an MNAR correction; "
            "missingness bounds must accompany, not be suppressed by, a positive joint score")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--private-island-join",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--bootstrap",type=int,default=499)
    args=parser.parse_args()
    result=run(pd.read_csv(args.private_island_join),bootstrap=args.bootstrap)
    body=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    forbidden=('"island_id"','"spatial_block"','"first_source"',
               '"species"','"latitude"','"longitude"','"source_region_id"')
    if any(key in body for key in forbidden):
        raise ValueError("nonredistributable source identifier escaped aggregate output")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],
        "n_tropical_candidates":result["support"]["n_tropical_candidate"],
        "n_qualified":result["support"]["n_map_and_direct_qualified"]},sort_keys=True))


if __name__ == "__main__":
    main()
