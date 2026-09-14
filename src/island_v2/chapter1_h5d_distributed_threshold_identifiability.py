"""Outcome-closed identifiability audit for distributed lineage thresholds.

Only the realized Palearctic island distance distribution and genus-support counts are
used. Observed floral/reproductive scores and observed genus-entry outcomes are never
opened. The audit asks whether a distributed-step generator can be distinguished from
a heterogeneous smooth-cline generator before any observed threshold distribution is fit.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml
from scipy.optimize import minimize
from scipy.special import expit, ndtr

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _clip(p: np.ndarray, eps: float) -> np.ndarray:
    return np.clip(np.asarray(p, dtype=float), float(eps), 1.0 - float(eps))


def _binomial_loglik(y: np.ndarray, n: np.ndarray, p: np.ndarray, eps: float) -> float:
    q = _clip(p, eps)
    return float(np.sum(y * np.log(q) + (n - y) * np.log1p(-q)))


def _fit_probit_threshold(x: np.ndarray, y: np.ndarray, n: np.ndarray, eps: float) -> tuple[float, float]:
    def objective(theta: np.ndarray) -> float:
        mu = float(theta[0])
        sigma = math.exp(float(theta[1]))
        p = ndtr((x - mu) / sigma)
        return -_binomial_loglik(y, n, p, eps)

    result = minimize(objective, np.array([0.0, math.log(0.5)]), method="L-BFGS-B", bounds=[(-4, 4), (math.log(0.05), math.log(4.0))])
    if not result.success:
        raise RuntimeError("probit threshold fit failed")
    return float(result.x[0]), float(math.exp(result.x[1]))


def _fit_logistic_cline(x: np.ndarray, y: np.ndarray, n: np.ndarray, eps: float) -> tuple[float, float]:
    def objective(theta: np.ndarray) -> float:
        p = expit(float(theta[0]) + float(theta[1]) * x)
        return -_binomial_loglik(y, n, p, eps)

    result = minimize(objective, np.array([0.0, 1.0]), method="L-BFGS-B", bounds=[(-8, 8), (-8, 8)])
    if not result.success:
        raise RuntimeError("logistic cline fit failed")
    return float(result.x[0]), float(result.x[1])


def _predict_probit(x: np.ndarray, mu: float, sigma: float, eps: float) -> np.ndarray:
    return _clip(ndtr((x - float(mu)) / float(sigma)), eps)


def _predict_logistic(x: np.ndarray, a: float, b: float, eps: float) -> np.ndarray:
    return _clip(expit(float(a) + float(b) * x), eps)


def _standardize(x: pd.Series) -> np.ndarray:
    values = pd.to_numeric(x, errors="coerce").to_numpy(float)
    mean = float(np.nanmean(values))
    sd = float(np.nanstd(values))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant distance design")
    return (values - mean) / sd


def build_design(artifact_root: Path, *, stratum: str, source_mode: str, context: str) -> pd.DataFrame:
    decomposition = pd.read_csv(artifact_root / "taxonomic-depth/all_analysis_eligible/decomposition.csv")
    covariates = pd.read_csv(artifact_root / "fixed/isolation/results/purpose_shortest_island_data.csv")
    realm = pd.read_csv(artifact_root / "fixed/realm/realm/island_biogeographic_realm_assignment.csv")
    work = decomposition.loc[
        decomposition["syndrome"].astype(str).eq("generalized_accessible")
        & decomposition["stratum"].astype(str).eq(str(stratum))
        & decomposition["source_mode"].astype(str).eq(str(source_mode))
    ][["island_id", "n_genera"]].copy()
    cov = covariates[["island_id", "log_distance_to_continent_km", "spatial_block"]].drop_duplicates("island_id")
    work = work.merge(cov, on="island_id", how="left", validate="one_to_one")
    work = work.merge(
        realm[["island_id", "biogeographic_realm"]].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="one_to_one",
    )
    work["n_genera"] = pd.to_numeric(work["n_genera"], errors="coerce")
    work["log_distance_to_continent_km"] = pd.to_numeric(work["log_distance_to_continent_km"], errors="coerce")
    work["spatial_block"] = work["spatial_block"].fillna("").astype(str)
    work = work.loc[work["biogeographic_realm"].astype(str).eq(str(context))]
    work = work.dropna(subset=["n_genera", "log_distance_to_continent_km"])
    work = work.loc[work["n_genera"].ge(2) & work["spatial_block"].ne("")].copy()
    work["z_distance"] = _standardize(work["log_distance_to_continent_km"])
    return work.reset_index(drop=True)


def _distributed_probability(x: np.ndarray, rng: np.random.Generator, *, sigma: float, latent_genera: int, eps: float) -> np.ndarray:
    thresholds = rng.normal(0.0, float(sigma), int(latent_genera))
    p = np.mean(x[:, None] > thresholds[None, :], axis=1)
    return _clip(p, eps)


def _smooth_probability(
    x: np.ndarray,
    rng: np.random.Generator,
    *,
    midpoint_sigma: float,
    log_slope_mean: float,
    log_slope_sd: float,
    latent_genera: int,
    eps: float,
) -> np.ndarray:
    midpoint = rng.normal(0.0, float(midpoint_sigma), int(latent_genera))
    slope = np.exp(rng.normal(float(log_slope_mean), float(log_slope_sd), int(latent_genera)))
    p = np.mean(expit(slope[None, :] * (x[:, None] - midpoint[None, :])), axis=1)
    return _clip(p, eps)


def _split_indices(n_rows: int, train_fraction: float, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    order = rng.permutation(n_rows)
    n_train = max(10, min(n_rows - 5, int(round(n_rows * float(train_fraction)))))
    return order[:n_train], order[n_train:]


def _classify_replicate(
    *,
    x: np.ndarray,
    n: np.ndarray,
    y: np.ndarray,
    train_fraction: float,
    rng: np.random.Generator,
    eps: float,
) -> tuple[str, float]:
    train, test = _split_indices(len(x), train_fraction, rng)
    mu, sigma = _fit_probit_threshold(x[train], y[train], n[train], eps)
    a, b = _fit_logistic_cline(x[train], y[train], n[train], eps)
    ll_distributed = _binomial_loglik(y[test], n[test], _predict_probit(x[test], mu, sigma, eps), eps)
    ll_smooth = _binomial_loglik(y[test], n[test], _predict_logistic(x[test], a, b, eps), eps)
    label = "distributed_threshold" if ll_distributed > ll_smooth else "smooth_cline"
    return label, sigma


def simulate_design(design: pd.DataFrame, config: dict[str, Any], *, seed_offset: int = 0) -> pd.DataFrame:
    sim = config["simulation"]
    rng = np.random.default_rng(int(sim["seed"]) + int(seed_offset))
    x = design["z_distance"].to_numpy(float)
    n = design["n_genera"].to_numpy(int)
    reps = int(sim["replicates"])
    eps = float(sim["probability_clip"])
    rows: list[dict[str, Any]] = []
    for sigma in [float(v) for v in sim["distributed_threshold_sigma"]]:
        for rep in range(reps):
            p = _distributed_probability(x, rng, sigma=sigma, latent_genera=int(sim["latent_genera"]), eps=eps)
            y = rng.binomial(n, p)
            predicted, fitted_sigma = _classify_replicate(
                x=x,
                n=n,
                y=y,
                train_fraction=float(sim["train_fraction"]),
                rng=rng,
                eps=eps,
            )
            rows.append(
                {
                    "generator": "distributed_threshold",
                    "true_sigma": sigma,
                    "replicate": rep,
                    "predicted": predicted,
                    "fitted_sigma": fitted_sigma,
                    "sigma_relative_error": abs(fitted_sigma - sigma) / sigma,
                }
            )
    for rep in range(reps):
        p = _smooth_probability(
            x,
            rng,
            midpoint_sigma=float(sim["smooth_cline_midpoint_sigma"]),
            log_slope_mean=float(sim["smooth_cline_log_slope_mean"]),
            log_slope_sd=float(sim["smooth_cline_log_slope_sd"]),
            latent_genera=int(sim["latent_genera"]),
            eps=eps,
        )
        y = rng.binomial(n, p)
        predicted, fitted_sigma = _classify_replicate(
            x=x,
            n=n,
            y=y,
            train_fraction=float(sim["train_fraction"]),
            rng=rng,
            eps=eps,
        )
        rows.append(
            {
                "generator": "smooth_cline",
                "true_sigma": float("nan"),
                "replicate": rep,
                "predicted": predicted,
                "fitted_sigma": fitted_sigma,
                "sigma_relative_error": float("nan"),
            }
        )
    return pd.DataFrame(rows)


def summarize_simulation(results: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    distributed = results.loc[results["generator"].eq("distributed_threshold")].copy()
    smooth = results.loc[results["generator"].eq("smooth_cline")].copy()
    correct_distributed = float((distributed["predicted"] == "distributed_threshold").mean())
    correct_smooth = float((smooth["predicted"] == "smooth_cline").mean())
    accuracy = float(
        (
            (results["generator"].eq("distributed_threshold") & results["predicted"].eq("distributed_threshold"))
            | (results["generator"].eq("smooth_cline") & results["predicted"].eq("smooth_cline"))
        ).mean()
    )
    false_distributed = float((smooth["predicted"] == "distributed_threshold").mean())
    sigma_recovery = float(
        pd.to_numeric(distributed["sigma_relative_error"], errors="coerce")
        .le(float(config["qualification"]["maximum_relative_sigma_error"]))
        .mean()
    )
    qualified = (
        accuracy >= float(config["qualification"]["minimum_classification_accuracy"])
        and false_distributed <= float(config["qualification"]["maximum_false_distributed_selection_under_smooth_clines"])
        and sigma_recovery >= float(config["qualification"]["minimum_sigma_recovery_fraction"])
    )
    return {
        "classification_accuracy": accuracy,
        "distributed_correct_fraction": correct_distributed,
        "smooth_correct_fraction": correct_smooth,
        "false_distributed_selection_under_smooth": false_distributed,
        "sigma_recovery_fraction": sigma_recovery,
        "qualified": bool(qualified),
    }


def run_audit(*, artifact_root: Path, config_path: Path, output_dir: Path) -> dict[str, Any]:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if config.get("contract") != "chapter1_h5d_distributed_threshold_identifiability_v1":
        raise ValueError("unexpected H5d contract")
    all_results: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []
    cell = 0
    for stratum in config["scope"]["strata"]:
        for source_mode in config["scope"]["source_modes"]:
            design = build_design(
                artifact_root,
                stratum=str(stratum),
                source_mode=str(source_mode),
                context=str(config["scope"]["context"]),
            )
            if len(design) < 30 or design["spatial_block"].nunique() < 10:
                summary_rows.append(
                    {
                        "stratum": str(stratum),
                        "source_mode": str(source_mode),
                        "n_islands": int(len(design)),
                        "n_blocks": int(design["spatial_block"].nunique()),
                        "status": "insufficient_design",
                        "qualified": False,
                    }
                )
                cell += 1
                continue
            results = simulate_design(design, config, seed_offset=cell * 10000)
            results.insert(0, "source_mode", str(source_mode))
            results.insert(0, "stratum", str(stratum))
            all_results.append(results)
            summary = summarize_simulation(results, config)
            summary_rows.append(
                {
                    "stratum": str(stratum),
                    "source_mode": str(source_mode),
                    "n_islands": int(len(design)),
                    "n_blocks": int(design["spatial_block"].nunique()),
                    "status": "simulated",
                    **summary,
                }
            )
            cell += 1
    result_df = pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()
    summary_df = pd.DataFrame(summary_rows)
    all_modes_required = bool(config["qualification"]["all_source_modes_required"])
    overall_qualified = bool(summary_df["qualified"].all()) if all_modes_required and len(summary_df) else bool(summary_df["qualified"].any())
    output_dir.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(output_dir / "h5d_identifiability_replicates.csv.gz", index=False, compression="gzip")
    summary_df.to_csv(output_dir / "h5d_identifiability_summary.csv", index=False)
    manifest = {
        "contract": config["contract"],
        "observed_threshold_distribution_opened": False,
        "n_design_cells": int(len(summary_df)),
        "n_qualified_cells": int(summary_df["qualified"].fillna(False).sum()),
        "overall_qualified": overall_qualified,
        "failure_action": config["failure_action"],
        "claim_ceiling": config["claim_ceiling"],
    }
    (output_dir / "h5d_identifiability_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("run")
def run_command(
    artifact_root: Path = typer.Option(..., exists=True, file_okay=False),
    config_path: Path = typer.Option(
        Path("config/chapter1_h5d_distributed_threshold_identifiability.yml"), exists=True
    ),
    output_dir: Path = typer.Option(...),
) -> None:
    typer.echo(json.dumps(run_audit(artifact_root=artifact_root, config_path=config_path, output_dir=output_dir), indent=2))


if __name__ == "__main__":
    app()
