"""Prospective N1 pollination-channel retention model for the Chapter 1 NEE challenge.

The model is deliberately frozen before N1 biological outcomes are inspected.  It
uses only source-available, effort-evaluable island x channel rows and tests the
single global null that isolation-retention slopes are exchangeable among qualified
channels.  Pairwise channel slopes are descriptive after that global test.
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
from scipy.stats import chi2

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_PROJECTED = {
    "island_id",
    "source_region_id",
    "channel_id",
    "source_state",
    "channel_state",
    "background_record_count",
    "background_spatial_units",
    "background_temporal_units",
    "distinct_dataset_count",
    "latest_background_year",
}
REQUIRED_QUALIFICATION = {"channel_id", "support_tier", "N1_gate_eligible"}
REQUIRED_COVARIATES = {
    "island_id",
    "log1p_distance_to_continent_km",
    "log_area",
    "climate_pc1",
    "climate_pc2",
    "climate_pc3",
    "climate_pc4",
    "spatial_block",
}
CONTINUOUS = [
    "isolation",
    "log_area",
    "climate_pc1",
    "climate_pc2",
    "climate_pc3",
    "climate_pc4",
    "effort_log_records",
    "effort_log_spatial",
    "effort_log_temporal",
    "effort_log_datasets",
    "effort_recency",
]


def load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("N1 model config must be a mapping")
    if payload.get("contract") != "chapter1_nee_n1_model_v1":
        raise typer.BadParameter("unexpected N1 model contract")
    return payload


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def _require(table: pd.DataFrame, required: set[str], label: str) -> None:
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"{label} missing columns: {sorted(missing)}")


def _standardize(series: pd.Series) -> tuple[pd.Series, dict[str, float]]:
    values = pd.to_numeric(series, errors="coerce")
    mean = float(values.mean())
    sd = float(values.std(ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError(f"constant or invalid continuous predictor: {series.name}")
    return (values - mean) / sd, {"mean": mean, "sd": sd}


def qualified_channels(qualification: pd.DataFrame, config: dict[str, Any]) -> list[str]:
    _require(qualification, REQUIRED_QUALIFICATION, "channel qualification")
    allowed = set(config["channel_gate"]["allowed_channels"])
    eligible = qualification.loc[
        qualification["channel_id"].astype(str).isin(allowed)
        & qualification["support_tier"].astype(str).eq("confirmatory")
        & qualification["N1_gate_eligible"].map(_as_bool),
        "channel_id",
    ].astype(str).drop_duplicates().tolist()
    return sorted(eligible)


def build_primary_frame(
    projected: pd.DataFrame,
    qualification: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Build the complete primary N1 support without inspecting plant outcomes."""
    _require(projected, REQUIRED_PROJECTED, "projected channel states")
    _require(covariates, REQUIRED_COVARIATES, "N1 covariates")
    if covariates["island_id"].astype(str).duplicated().any():
        raise ValueError("N1 covariates must contain one row per island_id")

    channels = qualified_channels(qualification, config)
    minimum = int(config["channel_gate"]["min_confirmatory_channels_for_primary_N1"])
    if len(channels) < minimum:
        return projected.head(0).copy(), {
            "status": "N1_not_evaluable_for_NEE_gate",
            "confirmatory_channels": channels,
            "n_confirmatory_channels": len(channels),
            "required_confirmatory_channels": minimum,
        }

    work = projected.loc[
        projected["channel_id"].astype(str).isin(channels)
        & projected["source_state"].astype(str).eq("available")
        & projected["channel_state"].astype(str).isin({"retained", "disrupted"})
    ].copy()
    if work[["island_id", "channel_id"]].duplicated().any():
        raise ValueError("primary N1 support contains duplicate island x channel rows")

    work["retained"] = work["channel_state"].astype(str).eq("retained").astype(int)
    for column in [
        "background_record_count",
        "background_spatial_units",
        "background_temporal_units",
        "distinct_dataset_count",
        "latest_background_year",
    ]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    effort_columns = list(config["observation_effort"]["required_columns"])
    if work[effort_columns].isna().any().any():
        raise ValueError("evaluable N1 row has missing observation-effort fields")
    if (work[effort_columns[:-1]] < 0).any().any():
        raise ValueError("N1 observation-effort counts must be non-negative")

    reference_year = int(config["observation_effort"]["reference_year"])
    work["effort_log_records"] = np.log1p(work["background_record_count"])
    work["effort_log_spatial"] = np.log1p(work["background_spatial_units"])
    work["effort_log_temporal"] = np.log1p(work["background_temporal_units"])
    work["effort_log_datasets"] = np.log1p(work["distinct_dataset_count"])
    work["effort_recency"] = reference_year - work["latest_background_year"]
    if (work["effort_recency"] < 0).any():
        raise ValueError("latest_background_year exceeds the frozen reference year")

    cov = covariates.copy()
    cov["island_id"] = cov["island_id"].astype(str)
    work["island_id"] = work["island_id"].astype(str)
    work = work.merge(cov, on="island_id", how="left", validate="many_to_one")
    work = work.rename(columns={"log1p_distance_to_continent_km": "isolation"})

    complete_columns = [*CONTINUOUS, "spatial_block", "source_region_id"]
    missing_mask = work[complete_columns].isna().any(axis=1)
    blank_mask = (
        work[["spatial_block", "source_region_id"]]
        .fillna("")
        .astype(str)
        .eq("")
        .any(axis=1)
    )
    work = work.loc[~missing_mask & ~blank_mask].copy()
    if work.empty:
        raise ValueError("N1 complete-case support is empty")

    # Pool sparse source regions without consulting the response.
    min_rows = int(
        config["primary_model"]["source_region_fixed_effects"]
        ["levels_with_fewer_than_5_primary_rows"]
    )
    counts = work["source_region_id"].astype(str).value_counts()
    sparse = set(counts.loc[counts < min_rows].index)
    work["source_region_model"] = work["source_region_id"].astype(str).where(
        ~work["source_region_id"].astype(str).isin(sparse), "other_source_region"
    )

    scaling: dict[str, dict[str, float]] = {}
    for column in CONTINUOUS:
        work[f"z_{column}"], scaling[column] = _standardize(work[column])

    reference = str(config["channel_gate"]["reference_channel"])
    if reference not in channels:
        reference = sorted(channels)[0]
    metadata = {
        "status": "qualified",
        "confirmatory_channels": channels,
        "reference_channel": reference,
        "n_rows": int(len(work)),
        "n_islands": int(work["island_id"].nunique()),
        "n_spatial_blocks": int(work["spatial_block"].astype(str).nunique()),
        "n_source_regions_raw": int(work["source_region_id"].astype(str).nunique()),
        "n_source_regions_model": int(work["source_region_model"].astype(str).nunique()),
        "continuous_scaling": scaling,
    }
    return work.reset_index(drop=True), metadata


def design_matrix(
    frame: pd.DataFrame,
    channels: list[str],
    reference_channel: str,
    heterogeneous: bool,
) -> tuple[np.ndarray, list[str], list[str]]:
    """Create the frozen numeric design; returns X, names, interaction names."""
    if reference_channel not in channels:
        raise ValueError("reference channel absent from model channels")
    names = ["intercept", *[f"z_{c}" for c in CONTINUOUS]]
    columns = [
        np.ones(len(frame)),
        *[frame[f"z_{c}"].to_numpy(float) for c in CONTINUOUS],
    ]

    channel_dummies: list[str] = []
    for channel in sorted(channels):
        if channel == reference_channel:
            continue
        name = f"channel[{channel}]"
        values = frame["channel_id"].astype(str).eq(channel).astype(float).to_numpy()
        columns.append(values)
        names.append(name)
        channel_dummies.append(channel)

    source_levels = sorted(frame["source_region_model"].astype(str).unique())
    if source_levels:
        for region in source_levels[1:]:
            columns.append(
                frame["source_region_model"].astype(str).eq(region).astype(float).to_numpy()
            )
            names.append(f"source_region[{region}]")
    interaction_names: list[str] = []
    if heterogeneous:
        isolation = frame["z_isolation"].to_numpy(float)
        for channel in channel_dummies:
            dummy = frame["channel_id"].astype(str).eq(channel).astype(float).to_numpy()
            name = f"z_isolation:channel[{channel}]"
            columns.append(isolation * dummy)
            names.append(name)
            interaction_names.append(name)
    return np.column_stack(columns), names, interaction_names


def _expit(x: np.ndarray) -> np.ndarray:
    x = np.clip(x, -35.0, 35.0)
    return 1.0 / (1.0 + np.exp(-x))


def fit_clustered_logit(
    frame: pd.DataFrame,
    X: np.ndarray,
    names: list[str],
    cluster_column: str = "spatial_block",
    max_iter: int = 100,
    tolerance: float = 1e-9,
) -> dict[str, Any]:
    """IRLS logistic regression with finite-sample cluster sandwich covariance."""
    y = frame["retained"].to_numpy(float)
    if len(y) <= X.shape[1] + 2:
        raise ValueError("insufficient N1 rows for design dimension")
    beta = np.zeros(X.shape[1], dtype=float)
    converged = False
    for iteration in range(1, max_iter + 1):
        eta = X @ beta
        p = np.clip(_expit(eta), 1e-8, 1 - 1e-8)
        w = p * (1 - p)
        z = eta + (y - p) / w
        xtwx = X.T @ (w[:, None] * X)
        beta_new = np.linalg.pinv(xtwx) @ (X.T @ (w * z))
        if float(np.max(np.abs(beta_new - beta))) < tolerance:
            beta = beta_new
            converged = True
            break
        beta = beta_new

    p = np.clip(_expit(X @ beta), 1e-10, 1 - 1e-10)
    w = p * (1 - p)
    bread = np.linalg.pinv(X.T @ (w[:, None] * X))
    residual = y - p
    clusters = frame[cluster_column].fillna("missing").astype(str).to_numpy()
    unique = np.unique(clusters)
    meat = np.zeros_like(bread)
    for cluster in unique:
        mask = clusters == cluster
        score = X[mask].T @ residual[mask]
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    n, k = X.shape
    g = len(unique)
    if g > 1 and n > k:
        covariance *= (g / (g - 1.0)) * ((n - 1.0) / (n - k))
    se = np.sqrt(np.clip(np.diag(covariance), 0, None))
    ll = float(np.sum(y * np.log(p) + (1 - y) * np.log(1 - p)))
    return {
        "beta": beta,
        "covariance": covariance,
        "standard_errors": se,
        "names": names,
        "converged": converged,
        "iterations": iteration,
        "log_likelihood": ll,
        "n_rows": int(n),
        "n_clusters": int(g),
    }


def global_heterogeneity_test(
    fit: dict[str, Any], interaction_names: list[str]
) -> dict[str, Any]:
    """Joint Wald test of H0: every isolation x channel contrast is zero."""
    if not interaction_names:
        raise ValueError("global heterogeneity test requires interaction terms")
    names = list(fit["names"])
    idx = [names.index(name) for name in interaction_names]
    b = np.asarray(fit["beta"], dtype=float)[idx]
    cov = np.asarray(fit["covariance"], dtype=float)[np.ix_(idx, idx)]
    statistic = float(b.T @ np.linalg.pinv(cov) @ b)
    df = len(idx)
    p_value = float(chi2.sf(statistic, df))
    return {
        "Wald_statistic": statistic,
        "df": df,
        "p_value": p_value,
        "interaction_terms": interaction_names,
    }


def channel_slopes(
    fit: dict[str, Any], channels: list[str], reference_channel: str
) -> pd.DataFrame:
    """Return standardized isolation slopes and cluster-robust SEs for each channel."""
    names = list(fit["names"])
    beta = np.asarray(fit["beta"], dtype=float)
    cov = np.asarray(fit["covariance"], dtype=float)
    isolation_idx = names.index("z_isolation")
    rows: list[dict[str, Any]] = []
    contrasts: dict[str, np.ndarray] = {}
    for channel in sorted(channels):
        c = np.zeros(len(names), dtype=float)
        c[isolation_idx] = 1.0
        if channel != reference_channel:
            interaction = f"z_isolation:channel[{channel}]"
            c[names.index(interaction)] = 1.0
        contrasts[channel] = c
        estimate = float(c @ beta)
        variance = float(c @ cov @ c)
        rows.append(
            {
                "channel_id": channel,
                "is_reference": channel == reference_channel,
                "isolation_slope_log_odds_per_sd": estimate,
                "cluster_robust_se": math.sqrt(max(variance, 0.0)),
            }
        )
    result = pd.DataFrame(rows)
    # Pairwise differences are intentionally descriptive after the global gate.
    differences: list[str] = []
    for i, left in enumerate(sorted(channels)):
        for right in sorted(channels)[i + 1 :]:
            c = contrasts[left] - contrasts[right]
            estimate = float(c @ beta)
            se = math.sqrt(max(float(c @ cov @ c), 0.0))
            differences.append(f"{left}-{right}:{estimate:.8g}|se={se:.8g}")
    result.attrs["pairwise_slope_differences"] = differences
    return result


def fit_n1_primary(
    projected: pd.DataFrame,
    qualification: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any], dict[str, Any]]:
    frame, metadata = build_primary_frame(projected, qualification, covariates, config)
    if metadata["status"] != "qualified":
        return frame, pd.DataFrame(), {}, metadata
    channels = list(metadata["confirmatory_channels"])
    reference = str(metadata["reference_channel"])
    X0, names0, _ = design_matrix(frame, channels, reference, heterogeneous=False)
    X1, names1, interactions = design_matrix(
        frame, channels, reference, heterogeneous=True
    )
    common = fit_clustered_logit(frame, X0, names0)
    heterogeneous = fit_clustered_logit(frame, X1, names1)
    global_test = global_heterogeneity_test(heterogeneous, interactions)
    slopes = channel_slopes(heterogeneous, channels, reference)
    coefficients = pd.DataFrame(
        {
            "predictor": heterogeneous["names"],
            "estimate_log_odds": heterogeneous["beta"],
            "cluster_robust_se": heterogeneous["standard_errors"],
        }
    )
    comparison = {
        "common_slope_log_likelihood": float(common["log_likelihood"]),
        "heterogeneous_slope_log_likelihood": float(
            heterogeneous["log_likelihood"]
        ),
        "heterogeneous_minus_common_log_likelihood": float(
            heterogeneous["log_likelihood"] - common["log_likelihood"]
        ),
        "same_primary_rows": True,
        "n_rows": int(len(frame)),
        "reference_channel": reference,
    }
    result_meta = {
        **metadata,
        "global_test": global_test,
        "model_comparison": comparison,
    }
    return frame, coefficients, global_test, {"slopes": slopes, **result_meta}


def deletion_robustness(
    frame: pd.DataFrame,
    channels: list[str],
    reference_channel: str,
    delete_column: str,
) -> pd.DataFrame:
    """Refit the frozen heterogeneous model after deleting each spatial/source block."""
    rows: list[dict[str, Any]] = []
    levels = sorted(frame[delete_column].astype(str).unique())
    for level in levels:
        subset = frame.loc[~frame[delete_column].astype(str).eq(level)].copy()
        if len(subset) <= 20 or subset["channel_id"].nunique() < 3:
            rows.append(
                {"deleted": level, "status": "not_evaluable", "p_value": np.nan}
            )
            continue
        try:
            X, names, interactions = design_matrix(
                subset, channels, reference_channel, heterogeneous=True
            )
            fit = fit_clustered_logit(subset, X, names)
            test = global_heterogeneity_test(fit, interactions)
            rows.append(
                {"deleted": level, "status": "fitted", "p_value": test["p_value"]}
            )
        except Exception as exc:  # noqa: BLE001
            rows.append(
                {
                    "deleted": level,
                    "status": "fit_failed",
                    "p_value": np.nan,
                    "error": str(exc),
                }
            )
    return pd.DataFrame(rows)


def gate_receipt(
    global_test: dict[str, Any],
    comparison: dict[str, Any],
    block_deletions: pd.DataFrame,
    source_deletions: pd.DataFrame,
    n_channels: int,
    config: dict[str, Any],
) -> dict[str, Any]:
    alpha = float(config["primary_global_test"]["alpha"])

    def reversal_fraction(table: pd.DataFrame) -> float:
        evaluable = table.loc[
            table["status"].eq("fitted") & table["p_value"].notna()
        ]
        if evaluable.empty:
            return 1.0
        return float((evaluable["p_value"] >= alpha).mean())

    block_reverse = reversal_fraction(block_deletions)
    source_reverse = reversal_fraction(source_deletions)
    pass_gate = (
        n_channels
        >= int(config["channel_gate"]["min_confirmatory_channels_for_primary_N1"])
        and float(global_test["p_value"]) < alpha
        and float(comparison["heterogeneous_minus_common_log_likelihood"]) > 0
        and block_reverse <= 0.20
        and source_reverse <= 0.20
    )
    return {
        "contract": config["contract"],
        "N1_pass": bool(pass_gate),
        "n_confirmatory_channels": int(n_channels),
        "global_Wald_p_value": float(global_test["p_value"]),
        "heterogeneous_minus_common_log_likelihood": float(
            comparison["heterogeneous_minus_common_log_likelihood"]
        ),
        "spatial_block_reversal_fraction": block_reverse,
        "source_region_reversal_fraction": source_reverse,
        "failure_action": (
            "stop_before_N2_and_keep_frozen_Chapter1" if not pass_gate else "N2_may_open"
        ),
        "claim_ceiling": (
            config["claim_ceiling"]["pass"] if pass_gate else "N1_not_promoted"
        ),
    }


@app.command("validate-config")
def validate_config_command(
    config_path: Path = typer.Option(
        Path("config/chapter1_nee_n1_model.yml"), exists=True
    ),
) -> None:
    config = load_config(config_path)
    required = int(config["channel_gate"]["min_confirmatory_channels_for_primary_N1"])
    if required < 3:
        raise typer.BadParameter(
            "primary N1 must require at least three confirmatory channels"
        )
    if (
        config["primary_global_test"]["method"]
        != "cluster_robust_joint_Wald_test_of_isolation_x_channel_terms"
    ):
        raise typer.BadParameter("N1 global test must remain the frozen joint Wald test")
    if config["N1_pass_rule"]["no_posthoc_channel_relabel_or_substitution"] is not True:
        raise typer.BadParameter("posthoc channel rescue must remain prohibited")
    typer.echo(json.dumps({"contract": config["contract"], "status": "valid"}))


if __name__ == "__main__":
    app()
