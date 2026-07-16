"""Repeated spatial-block cross-validation for incremental M0-M4 model evidence.

The core question is predictive and nested: does adding Bombus environmental
compatibility improve held-out geographic-block prediction beyond the frozen
baseline, and does the M4 interaction improve prediction beyond its additive
components? Whole 10-degree spatial blocks are assigned to folds, so no block is
split between training and testing.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)


@dataclass
class LogisticModel:
    beta: np.ndarray
    predictors: list[str]
    means: pd.Series
    sds: pd.Series


def _load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("full-analysis config must be a mapping")
    return payload


def _expit(values: np.ndarray) -> np.ndarray:
    values = np.clip(values, -35.0, 35.0)
    return 1.0 / (1.0 + np.exp(-values))


def _complete_rows(
    data: pd.DataFrame,
    successes: str,
    trials: str,
    predictors: list[str],
) -> pd.DataFrame:
    columns = ["island_id", "spatial_block", successes, trials, *predictors]
    missing = set(columns).difference(data.columns)
    if missing:
        raise typer.BadParameter(f"spatial CV input missing columns: {sorted(missing)}")
    work = data[columns].copy()
    for column in [successes, trials, *predictors]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.dropna()
    work = work.loc[
        (work[trials] > 0)
        & (work[successes] >= 0)
        & (work[successes] <= work[trials])
        & work["spatial_block"].astype(str).ne("missing")
    ].copy()
    return work


def fit_grouped_binomial(
    train: pd.DataFrame,
    successes: str,
    trials: str,
    predictors: list[str],
    max_iter: int = 100,
    tolerance: float = 1e-9,
) -> LogisticModel:
    y_success = train[successes].to_numpy(dtype=float)
    n_trials = train[trials].to_numpy(dtype=float)
    y = y_success / n_trials
    X_raw = train[predictors].astype(float)
    means = X_raw.mean(axis=0)
    sds = X_raw.std(axis=0, ddof=0)
    bad = sds.index[(~np.isfinite(sds)) | (sds <= 0)].tolist()
    if bad:
        raise ValueError(f"constant or invalid training predictors: {bad}")
    X_std = (X_raw - means) / sds
    X = np.column_stack([np.ones(len(train)), X_std.to_numpy(dtype=float)])
    beta = np.zeros(X.shape[1], dtype=float)

    for _ in range(max_iter):
        eta = X @ beta
        probability = np.clip(_expit(eta), 1e-8, 1.0 - 1e-8)
        variance = probability * (1.0 - probability)
        weights = n_trials * variance
        working_response = eta + (y - probability) / variance
        beta_new = np.linalg.pinv(X.T @ (weights[:, None] * X)) @ (
            X.T @ (weights * working_response)
        )
        if float(np.max(np.abs(beta_new - beta))) < tolerance:
            beta = beta_new
            break
        beta = beta_new
    return LogisticModel(beta=beta, predictors=predictors, means=means, sds=sds)


def predict_probability(model: LogisticModel, test: pd.DataFrame) -> np.ndarray:
    X_raw = test[model.predictors].astype(float)
    X_std = (X_raw - model.means) / model.sds
    X = np.column_stack([np.ones(len(test)), X_std.to_numpy(dtype=float)])
    return np.clip(_expit(X @ model.beta), 1e-10, 1.0 - 1e-10)


def heldout_metrics(
    test: pd.DataFrame,
    successes: str,
    trials: str,
    probability: np.ndarray,
) -> dict[str, float]:
    success = test[successes].to_numpy(dtype=float)
    total = test[trials].to_numpy(dtype=float)
    observed = success / total
    log_likelihood = float(
        np.sum(success * np.log(probability) + (total - success) * np.log(1.0 - probability))
    )
    brier = float(np.sum(total * (observed - probability) ** 2) / np.sum(total))
    return {
        "log_likelihood": log_likelihood,
        "weighted_brier": brier,
        "n_islands": int(len(test)),
        "n_trials": float(total.sum()),
    }


def assign_balanced_block_folds(
    data: pd.DataFrame,
    n_folds: int,
    seed: int,
) -> dict[str, int]:
    block_counts = data.groupby("spatial_block", sort=True).size().rename("n_islands")
    if len(block_counts) < n_folds:
        raise typer.BadParameter(
            f"fewer spatial blocks ({len(block_counts)}) than folds ({n_folds})"
        )
    rng = np.random.default_rng(seed)
    order = pd.DataFrame(
        {
            "spatial_block": block_counts.index.astype(str),
            "n_islands": block_counts.to_numpy(dtype=int),
            "tie_break": rng.random(len(block_counts)),
        }
    ).sort_values(["n_islands", "tie_break"], ascending=[False, True])
    fold_loads = [0 for _ in range(n_folds)]
    mapping: dict[str, int] = {}
    for row in order.itertuples(index=False):
        fold = min(range(n_folds), key=lambda index: (fold_loads[index], index))
        mapping[str(row.spatial_block)] = fold
        fold_loads[fold] += int(row.n_islands)
    return mapping


def _comparison_specs(baseline: list[str]) -> list[dict[str, Any]]:
    no_latitude = [column for column in baseline if column != "abs_latitude"]
    return [
        {
            "comparison": "M0_vs_M1_full_baseline",
            "successes": "floral_phenotype_successes",
            "trials": "floral_phenotype_trials",
            "reduced": baseline,
            "full": [*baseline, "bombus_channel_state"],
        },
        {
            "comparison": "M0_vs_M1_no_abs_latitude",
            "successes": "floral_phenotype_successes",
            "trials": "floral_phenotype_trials",
            "reduced": no_latitude,
            "full": [*no_latitude, "bombus_channel_state"],
        },
        {
            "comparison": "M2_vs_M3_full_baseline",
            "successes": "floral_phenotype_successes",
            "trials": "floral_phenotype_trials",
            "reduced": [*baseline, "mixed_mating_share", "self_compatibility_share"],
            "full": [
                *baseline,
                "mixed_mating_share",
                "self_compatibility_share",
                "bombus_channel_state",
            ],
        },
        {
            "comparison": "M2_vs_M3_no_abs_latitude",
            "successes": "floral_phenotype_successes",
            "trials": "floral_phenotype_trials",
            "reduced": [*no_latitude, "mixed_mating_share", "self_compatibility_share"],
            "full": [
                *no_latitude,
                "mixed_mating_share",
                "self_compatibility_share",
                "bombus_channel_state",
            ],
        },
        {
            "comparison": "M4_additive_vs_interaction",
            "successes": "floral_phenotype_successes",
            "trials": "floral_phenotype_trials",
            "reduced": [*baseline, "bombus_channel_state", "alternative_pollen_vector"],
            "full": [
                *baseline,
                "bombus_channel_state",
                "alternative_pollen_vector",
                "bombus_x_alternative",
            ],
        },
    ]


def run_spatial_cv(
    data: pd.DataFrame,
    baseline: list[str],
    n_folds: int,
    repeats: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if n_folds < 5:
        raise typer.BadParameter("spatial CV requires at least 5 folds")
    if repeats < 1:
        raise typer.BadParameter("spatial CV repeats must be positive")

    fold_rows: list[dict[str, Any]] = []
    for spec in _comparison_specs(baseline):
        union_predictors = list(dict.fromkeys([*spec["reduced"], *spec["full"]]))
        complete = _complete_rows(
            data,
            spec["successes"],
            spec["trials"],
            union_predictors,
        )
        for repeat in range(repeats):
            repeat_seed = seed + repeat * 1009
            mapping = assign_balanced_block_folds(complete, n_folds, repeat_seed)
            fold_id = complete["spatial_block"].astype(str).map(mapping)
            for fold in range(n_folds):
                test = complete.loc[fold_id.eq(fold)].copy()
                train = complete.loc[~fold_id.eq(fold)].copy()
                if test.empty or train.empty:
                    continue
                reduced_model = fit_grouped_binomial(
                    train,
                    spec["successes"],
                    spec["trials"],
                    spec["reduced"],
                )
                full_model = fit_grouped_binomial(
                    train,
                    spec["successes"],
                    spec["trials"],
                    spec["full"],
                )
                reduced_metrics = heldout_metrics(
                    test,
                    spec["successes"],
                    spec["trials"],
                    predict_probability(reduced_model, test),
                )
                full_metrics = heldout_metrics(
                    test,
                    spec["successes"],
                    spec["trials"],
                    predict_probability(full_model, test),
                )
                fold_rows.append(
                    {
                        "comparison": spec["comparison"],
                        "repeat": repeat + 1,
                        "fold": fold + 1,
                        "seed": repeat_seed,
                        "n_train_islands": int(len(train)),
                        "n_test_islands": int(len(test)),
                        "n_test_spatial_blocks": int(test["spatial_block"].nunique()),
                        "n_test_trials": full_metrics["n_trials"],
                        "reduced_log_likelihood": reduced_metrics["log_likelihood"],
                        "full_log_likelihood": full_metrics["log_likelihood"],
                        "delta_log_likelihood_full_minus_reduced": (
                            full_metrics["log_likelihood"] - reduced_metrics["log_likelihood"]
                        ),
                        "reduced_weighted_brier": reduced_metrics["weighted_brier"],
                        "full_weighted_brier": full_metrics["weighted_brier"],
                        "delta_brier_full_minus_reduced": (
                            full_metrics["weighted_brier"] - reduced_metrics["weighted_brier"]
                        ),
                    }
                )

    folds = pd.DataFrame(fold_rows)
    repeat_rows: list[dict[str, Any]] = []
    for (comparison, repeat), group in folds.groupby(["comparison", "repeat"], sort=True):
        repeat_rows.append(
            {
                "comparison": comparison,
                "repeat": int(repeat),
                "n_folds": int(len(group)),
                "total_test_trials": float(group["n_test_trials"].sum()),
                "delta_log_likelihood_full_minus_reduced": float(
                    group["delta_log_likelihood_full_minus_reduced"].sum()
                ),
                "delta_log_likelihood_per_1000_trials": float(
                    1000.0
                    * group["delta_log_likelihood_full_minus_reduced"].sum()
                    / group["n_test_trials"].sum()
                ),
                "mean_delta_brier_full_minus_reduced": float(
                    np.average(
                        group["delta_brier_full_minus_reduced"],
                        weights=group["n_test_trials"],
                    )
                ),
            }
        )
    repeats_table = pd.DataFrame(repeat_rows)

    summary_rows: list[dict[str, Any]] = []
    for comparison, group in repeats_table.groupby("comparison", sort=True):
        ll = group["delta_log_likelihood_per_1000_trials"].to_numpy(dtype=float)
        brier = group["mean_delta_brier_full_minus_reduced"].to_numpy(dtype=float)
        summary_rows.append(
            {
                "comparison": comparison,
                "n_repeats": int(len(group)),
                "median_delta_log_likelihood_per_1000_trials": float(np.median(ll)),
                "min_delta_log_likelihood_per_1000_trials": float(np.min(ll)),
                "max_delta_log_likelihood_per_1000_trials": float(np.max(ll)),
                "proportion_repeats_full_model_better_loglik": float(np.mean(ll > 0)),
                "median_delta_brier_full_minus_reduced": float(np.median(brier)),
                "proportion_repeats_full_model_better_brier": float(np.mean(brier < 0)),
                "predictive_gain_consistent": bool(np.all(ll > 0) and np.all(brier < 0)),
            }
        )
    return folds, repeats_table, pd.DataFrame(summary_rows)


@app.command("run")
def run(
    island_data_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    analysis_tier: str = typer.Option(...),
    config_path: Path = typer.Option(Path("config/m0_m4_full_analysis.yml"), exists=True),
    n_folds: int = typer.Option(10, min=5),
    repeats: int = typer.Option(5, min=1),
    seed: int = typer.Option(20260713),
) -> None:
    config = _load_config(config_path)
    data = pd.read_csv(island_data_csv)
    baseline = [str(value) for value in config["baseline_covariates"]]
    folds, repeat_table, summary = run_spatial_cv(
        data,
        baseline,
        n_folds=n_folds,
        repeats=repeats,
        seed=seed,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    folds.to_csv(output_dir / "spatial_cv_fold_results.csv", index=False)
    repeat_table.to_csv(output_dir / "spatial_cv_repeat_results.csv", index=False)
    summary.to_csv(output_dir / "spatial_cv_summary.csv", index=False)
    manifest = {
        "contract": "m0_m4_repeated_spatial_block_cv_v1",
        "analysis_tier": analysis_tier,
        "n_folds": n_folds,
        "repeats": repeats,
        "seed": seed,
        "spatial_unit": "frozen_10_degree_spatial_block",
        "n_comparisons": int(summary["comparison"].nunique()),
        "comparisons": summary.to_dict("records"),
        "interpretation": (
            "Positive held-out delta log-likelihood and negative delta Brier indicate that the full model "
            "adds out-of-region predictive information beyond the nested reduced model."
        ),
    }
    (output_dir / "spatial_cv_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
