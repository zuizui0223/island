#!/usr/bin/env python3
"""Repeated island flower-colour resampling diagnostic.

This is an exploratory bridge from the RGFCA repeated-atlas idea into the
island Chapter 1 data model. It does not search for a common geographic
boundary. Instead it asks whether colour-composition slopes along the frozen
distance-to-continent gradient recur after equal-depth random resampling of
directly colour-resolved native species within islands.

The diagnostic is deliberately separate from the canonical Chapter 1
when/where inference.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

PRIMARY_CONTEXTS = ("northern_midlatitude", "tropical")
OUTCOMES = ("white_present", "muted_only", "conspicuous_present")
CONTROL_COLUMNS = (
    "log_island_area_km2",
    "climate_pc1",
    "climate_pc2",
    "climate_pc3",
    "climate_pc4",
)
REQUIRED_COVARIATES = (
    "island_id",
    "analysis_regime",
    "spatial_block",
    "log_distance_to_continent_km",
    *CONTROL_COLUMNS,
)

PLAIN = {"white", "green_brown_inconspicuous"}
CONSPICUOUS = {"yellow_orange", "red_pink", "blue_purple"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--status-flora-csv", type=Path, required=True)
    p.add_argument("--trait-audit-csv", type=Path, required=True)
    p.add_argument("--covariates-csv", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--draws", type=int, default=1000)
    p.add_argument("--species-per-island", type=int, default=20)
    p.add_argument("--depth-draws", type=int, default=500)
    p.add_argument("--depth-common-min", type=int, default=50)
    p.add_argument("--depth-grid", type=int, nargs="+", default=[10, 20, 50])
    p.add_argument("--seed", type=int, default=20260908)
    return p.parse_args()


def _tokens(value: object) -> set[str]:
    if pd.isna(value):
        return set()
    return {v.strip() for v in str(value).split("|") if v.strip()}


def build_colour_species(trait_audit: pd.DataFrame) -> pd.DataFrame:
    required = {
        "accepted_species",
        "trait_name",
        "resolved_for_primary",
        "canonical_signature",
    }
    missing = required - set(trait_audit.columns)
    if missing:
        raise ValueError(f"trait audit missing columns: {sorted(missing)}")

    frame = trait_audit.loc[
        trait_audit["trait_name"].eq("flower_primary_color")
        & trait_audit["resolved_for_primary"].eq(True),
        ["accepted_species", "canonical_signature"],
    ].copy()
    if frame["accepted_species"].duplicated().any():
        dup = frame.loc[frame["accepted_species"].duplicated(), "accepted_species"].head().tolist()
        raise ValueError(f"resolved colour audit has duplicated species: {dup}")

    frame["tokens"] = frame["canonical_signature"].map(_tokens)
    frame["white_present"] = frame["tokens"].map(lambda x: "white" in x)
    frame["muted_only"] = frame["tokens"].map(lambda x: bool(x) and x.issubset(PLAIN))
    frame["conspicuous_present"] = frame["tokens"].map(lambda x: bool(x & CONSPICUOUS))
    return frame[["accepted_species", "white_present", "muted_only", "conspicuous_present"]]


def stratum_rows(status_flora: pd.DataFrame, stratum: str) -> pd.DataFrame:
    required = {"island_id", "accepted_species", "origin_status", "floristic_status"}
    missing = required - set(status_flora.columns)
    if missing:
        raise ValueError(f"status flora missing columns: {sorted(missing)}")
    if stratum == "all_native":
        mask = status_flora["origin_status"].eq("native")
    elif stratum == "native_nonendemic":
        mask = status_flora["floristic_status"].eq("native_nonendemic")
    else:
        raise ValueError(f"unsupported stratum: {stratum}")
    return status_flora.loc[mask, ["island_id", "accepted_species"]].drop_duplicates()


def build_joined(
    status_flora: pd.DataFrame,
    colour_species: pd.DataFrame,
    covariates: pd.DataFrame,
    stratum: str,
) -> pd.DataFrame:
    missing = set(REQUIRED_COVARIATES) - set(covariates.columns)
    if missing:
        raise ValueError(f"covariates missing columns: {sorted(missing)}")
    cov = covariates[list(REQUIRED_COVARIATES)].drop_duplicates("island_id").copy()
    if cov["island_id"].duplicated().any():
        raise ValueError("covariate table contains duplicate island_id")
    return (
        stratum_rows(status_flora, stratum)
        .merge(colour_species, on="accepted_species", how="inner", validate="many_to_one")
        .merge(cov, on="island_id", how="left", validate="many_to_one")
    )


def _design(cov: pd.DataFrame) -> tuple[np.ndarray, float]:
    x = cov["log_distance_to_continent_km"].to_numpy(float)
    xz = (x - x.mean()) / x.std(ddof=0)
    controls = cov[list(CONTROL_COLUMNS)].to_numpy(float)
    controls = (controls - controls.mean(axis=0)) / controls.std(axis=0, ddof=0)
    z = np.column_stack([np.ones(len(cov)), controls])
    residual_x = xz - z @ np.linalg.lstsq(z, xz, rcond=None)[0]
    return residual_x, float(residual_x @ residual_x)


def _sample_realisations(
    data: pd.DataFrame,
    cov: pd.DataFrame,
    *,
    n_draws: int,
    k: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    groups = {
        island_id: group[list(OUTCOMES)].astype(float).to_numpy()
        for island_id, group in data.groupby("island_id", sort=False)
    }
    ids = cov["island_id"].tolist()
    out = np.empty((n_draws, len(ids), len(OUTCOMES)), dtype=float)
    for r in range(n_draws):
        for i, island_id in enumerate(ids):
            arr = groups[island_id]
            if arr.shape[0] < k:
                raise RuntimeError("eligibility drift: island has fewer rows than fixed k")
            idx = rng.choice(arr.shape[0], size=k, replace=False)
            out[r, i, :] = arr[idx].mean(axis=0)
    return out


def _block_null(
    realisations: np.ndarray,
    cov: pd.DataFrame,
    residual_x: np.ndarray,
    denominator: float,
    *,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    blocks = cov["spatial_block"].astype(str).to_numpy()
    unique_blocks = np.unique(blocks)
    null = np.empty((realisations.shape[0], realisations.shape[2]), dtype=float)
    for r in range(realisations.shape[0]):
        perm = np.arange(len(cov))
        for block in unique_blocks:
            idx = np.flatnonzero(blocks == block)
            if len(idx) > 1:
                perm[idx] = rng.permutation(idx)
        null[r, :] = residual_x @ realisations[r, perm, :] / denominator
    return null


def repeated_summary(
    joined: pd.DataFrame,
    *,
    stratum: str,
    context: str,
    k: int,
    n_draws: int,
    seed: int,
) -> list[dict[str, object]]:
    counts = (
        joined.groupby(["island_id", "analysis_regime"], as_index=False)
        .agg(n_colour_resolved_species=("accepted_species", "nunique"))
    )
    eligible = counts.loc[
        counts["analysis_regime"].eq(context)
        & counts["n_colour_resolved_species"].ge(k),
        "island_id",
    ]
    cov = (
        joined.loc[joined["island_id"].isin(eligible), list(REQUIRED_COVARIATES)]
        .drop_duplicates("island_id")
        .dropna(subset=["log_distance_to_continent_km", *CONTROL_COLUMNS])
        .sort_values("island_id")
        .reset_index(drop=True)
    )
    data = joined.loc[joined["island_id"].isin(cov["island_id"])].copy()
    residual_x, denominator = _design(cov)
    realisations = _sample_realisations(data, cov, n_draws=n_draws, k=k, seed=seed)
    beta = np.einsum("i,rij->rj", residual_x, realisations) / denominator
    null = _block_null(
        realisations,
        cov,
        residual_x,
        denominator,
        seed=seed + 1000003,
    )

    rows: list[dict[str, object]] = []
    for j, outcome in enumerate(OUTCOMES):
        b = beta[:, j]
        n = null[:, j]
        rows.append(
            {
                "stratum": stratum,
                "context": context,
                "outcome": outcome,
                "n_islands": int(len(cov)),
                "species_per_island": int(k),
                "draws": int(n_draws),
                "adjusted_beta_mean": float(b.mean()),
                "adjusted_beta_median": float(np.median(b)),
                "adjusted_beta_q025": float(np.quantile(b, 0.025)),
                "adjusted_beta_q975": float(np.quantile(b, 0.975)),
                "realisations_beta_positive": float(np.mean(b > 0)),
                "block_null_mean": float(n.mean()),
                "block_null_q025": float(np.quantile(n, 0.025)),
                "block_null_q975": float(np.quantile(n, 0.975)),
                "paired_observed_gt_block_null": float(np.mean(b > n)),
                "paired_observed_lt_block_null": float(np.mean(b < n)),
            }
        )
    return rows


def depth_sensitivity(
    joined: pd.DataFrame,
    *,
    stratum: str,
    context: str,
    common_min: int,
    depth_grid: list[int],
    n_draws: int,
    seed: int,
) -> list[dict[str, object]]:
    counts = (
        joined.groupby(["island_id", "analysis_regime"], as_index=False)
        .agg(n_colour_resolved_species=("accepted_species", "nunique"))
    )
    eligible = counts.loc[
        counts["analysis_regime"].eq(context)
        & counts["n_colour_resolved_species"].ge(common_min),
        "island_id",
    ]
    cov = (
        joined.loc[joined["island_id"].isin(eligible), list(REQUIRED_COVARIATES)]
        .drop_duplicates("island_id")
        .dropna(subset=["log_distance_to_continent_km", *CONTROL_COLUMNS])
        .sort_values("island_id")
        .reset_index(drop=True)
    )
    data = joined.loc[joined["island_id"].isin(cov["island_id"])].copy()
    residual_x, denominator = _design(cov)
    rows: list[dict[str, object]] = []
    for k in depth_grid:
        if k > common_min:
            raise ValueError("depth grid cannot exceed common eligibility minimum")
        realisations = _sample_realisations(
            data,
            cov,
            n_draws=n_draws,
            k=k,
            seed=seed + k * 101,
        )
        beta = np.einsum("i,rij->rj", residual_x, realisations) / denominator
        for j, outcome in enumerate(OUTCOMES):
            b = beta[:, j]
            rows.append(
                {
                    "stratum": stratum,
                    "context": context,
                    "outcome": outcome,
                    "common_eligibility_min_species": int(common_min),
                    "n_common_eligible_islands": int(len(cov)),
                    "species_per_island": int(k),
                    "draws": int(n_draws),
                    "adjusted_beta_mean": float(b.mean()),
                    "adjusted_beta_q025": float(np.quantile(b, 0.025)),
                    "adjusted_beta_q975": float(np.quantile(b, 0.975)),
                    "realisations_beta_positive": float(np.mean(b > 0)),
                }
            )
    return rows


def main() -> None:
    args = parse_args()
    if args.draws < 1 or args.depth_draws < 1:
        raise ValueError("draw counts must be positive")
    if args.species_per_island < 1:
        raise ValueError("species-per-island must be positive")

    status = pd.read_csv(args.status_flora_csv)
    audit = pd.read_csv(args.trait_audit_csv)
    covariates = pd.read_csv(args.covariates_csv)
    colours = build_colour_species(audit)

    summary_rows: list[dict[str, object]] = []
    depth_rows: list[dict[str, object]] = []
    stratum_counts: dict[str, object] = {}
    for s_idx, stratum in enumerate(("all_native", "native_nonendemic")):
        joined = build_joined(status, colours, covariates, stratum)
        stratum_counts[stratum] = {
            "joined_rows": int(len(joined)),
            "distinct_species": int(joined["accepted_species"].nunique()),
            "distinct_islands": int(joined["island_id"].nunique()),
        }
        for c_idx, context in enumerate(PRIMARY_CONTEXTS):
            summary_rows.extend(
                repeated_summary(
                    joined,
                    stratum=stratum,
                    context=context,
                    k=args.species_per_island,
                    n_draws=args.draws,
                    seed=args.seed + s_idx * 10000 + c_idx * 1000,
                )
            )
            depth_rows.extend(
                depth_sensitivity(
                    joined,
                    stratum=stratum,
                    context=context,
                    common_min=args.depth_common_min,
                    depth_grid=list(args.depth_grid),
                    n_draws=args.depth_draws,
                    seed=args.seed + 500000 + s_idx * 10000 + c_idx * 1000,
                )
            )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(summary_rows)
    depth = pd.DataFrame(depth_rows)
    summary.to_csv(args.output_dir / "rgfca_island_colour_repeated_summary.csv", index=False)
    depth.to_csv(args.output_dir / "rgfca_island_colour_depth_sensitivity.csv", index=False)

    manifest = {
        "contract": "rgfca_island_colour_repeated_exploratory_v1",
        "status": "exploratory_not_canonical_chapter1_inference",
        "question": "Does island flower-colour composition become whiter or more muted along the frozen distance-to-continent gradient under repeated equal-depth species resampling?",
        "contexts": list(PRIMARY_CONTEXTS),
        "strata": ["all_native", "native_nonendemic"],
        "outcomes": {
            "white_present": "resolved species colour signature contains white; may also contain conspicuous states",
            "muted_only": "all resolved colour states are in {white, green_brown_inconspicuous}",
            "conspicuous_present": "resolved signature contains at least one of {yellow_orange, red_pink, blue_purple}",
        },
        "resampling": {
            "draws": int(args.draws),
            "species_per_island": int(args.species_per_island),
            "sampling": "without replacement within each island; island-equal depth and island-equal regression weight",
            "seed": int(args.seed),
        },
        "model": {
            "response": "random-realisation island species share",
            "distance_predictor": "within-context standardized log_distance_to_continent_km",
            "controls": list(CONTROL_COLUMNS),
            "coefficient": "OLS proportion-scale beta via Frisch-Waugh-Lovell residualization",
            "spatial_null": "one paired permutation per realization of the island response vector within spatial_block",
        },
        "interpretation_guards": [
            "resampling quantiles describe stability to equal-depth species resampling; they are not confidence intervals for all plants",
            "realisations_beta_positive is a reproducibility fraction, not a p-value",
            "paired block-null dominance fractions are descriptive, not confirmatory permutation p-values",
            "native_nonendemic persistence supports assemblage/source-pool filtering but does not establish within-lineage evolutionary whitening",
            "white_present and conspicuous_present can overlap for multistate species; muted_only is the non-overlapping strict muted endpoint",
            "this lane does not replace the frozen Chapter 1 when/where analysis",
        ],
        "input_counts": stratum_counts,
    }
    (args.output_dir / "rgfca_island_colour_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
