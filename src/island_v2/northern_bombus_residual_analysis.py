"""Test whether Bombus hypervolume compatibility adds information beyond lineage.

The response is the observed-minus-genus-fixed-null trait composition produced by
``genus_fixed_trait_null``. Only the northern mid-latitude regime is confirmatory
for this diagnostic. Environmental compatibility is a climatic niche proxy only;
it is not occurrence, abundance, visitation, pollination service, or historical
Bombus loss.

Three nested models are compared:

M0: distance + area + climate PCs
M1: M0 + Bombus environmental compatibility
M2: M1 + distance x compatibility

M2 is a sensitivity model for the specific "remote but environmentally compatible"
arrival-limitation prediction. Incremental value is evaluated with deterministic
spatial-block cross-validation in addition to full-data cluster-robust coefficients.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer

from island_v2.status_stratified_lineage_analysis import fit_weighted_linear_clustered

app = typer.Typer(add_completion=False, no_args_is_help=True)

BASELINE = [
    "log_distance_to_continent_km",
    "log_island_area_km2",
    "climate_pc1",
    "climate_pc2",
    "climate_pc3",
    "climate_pc4",
]
COMPATIBILITY = "environmental_compatibility_max"
INTERACTION = "distance_x_bombus_compatibility"
MODELS = {
    "M0_lineage_geography_climate": BASELINE,
    "M1_add_bombus_compatibility": [*BASELINE, COMPATIBILITY],
    "M2_add_distance_x_compatibility": [*BASELINE, COMPATIBILITY, INTERACTION],
}
DEFAULT_STRATA = ("all_native", "native_nonendemic")


def _add_interaction(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    distance = pd.to_numeric(out["log_distance_to_continent_km"], errors="coerce")
    compatibility = pd.to_numeric(out[COMPATIBILITY], errors="coerce")
    d_sd = float(distance.std(ddof=0))
    c_sd = float(compatibility.std(ddof=0))
    if not math.isfinite(d_sd) or d_sd <= 0 or not math.isfinite(c_sd) or c_sd <= 0:
        out[INTERACTION] = np.nan
        return out
    z_d = (distance - float(distance.mean())) / d_sd
    z_c = (compatibility - float(compatibility.mean())) / c_sd
    out[INTERACTION] = z_d * z_c
    return out


def _fold_for_block(block: object, n_folds: int) -> int:
    digest = hashlib.sha256(str(block).encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % n_folds


def _fit_predict_wls(
    train: pd.DataFrame,
    test: pd.DataFrame,
    predictors: list[str],
    response: str,
    weight: str,
) -> np.ndarray:
    train_columns: list[np.ndarray] = []
    test_columns: list[np.ndarray] = []
    for predictor in predictors:
        tr = pd.to_numeric(train[predictor], errors="coerce").to_numpy(dtype=float)
        te = pd.to_numeric(test[predictor], errors="coerce").to_numpy(dtype=float)
        mean = float(np.mean(tr))
        sd = float(np.std(tr, ddof=0))
        if not math.isfinite(sd) or sd <= 0:
            raise ValueError(f"constant predictor in CV training fold: {predictor}")
        train_columns.append((tr - mean) / sd)
        test_columns.append((te - mean) / sd)
    X_train = np.column_stack([np.ones(len(train)), *train_columns])
    X_test = np.column_stack([np.ones(len(test)), *test_columns])
    y = pd.to_numeric(train[response], errors="coerce").to_numpy(dtype=float)
    w = pd.to_numeric(train[weight], errors="coerce").to_numpy(dtype=float)
    beta = np.linalg.pinv(X_train.T @ (w[:, None] * X_train)) @ (X_train.T @ (w * y))
    return X_test @ beta


def spatial_block_cv(
    data: pd.DataFrame,
    *,
    predictors: list[str],
    response: str = "deviation_observed_minus_null",
    weight: str = "observed_n_species",
    cluster: str = "spatial_block",
    n_folds: int = 5,
) -> dict[str, Any]:
    required = {response, weight, cluster, *predictors}
    work = data[list(required)].copy()
    for column in [response, weight, *predictors]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.dropna()
    work = work.loc[work[weight] > 0].copy()
    if len(work) < 30 or work[cluster].nunique() < n_folds:
        return {"status": "insufficient_support", "n_islands": int(len(work))}
    work["cv_fold"] = work[cluster].map(lambda value: _fold_for_block(value, n_folds))
    records: list[pd.DataFrame] = []
    for fold in range(n_folds):
        test = work.loc[work["cv_fold"].eq(fold)].copy()
        train = work.loc[~work["cv_fold"].eq(fold)].copy()
        if test.empty or train.empty:
            continue
        try:
            prediction = _fit_predict_wls(train, test, predictors, response, weight)
        except ValueError:
            return {"status": "constant_predictor_in_fold", "n_islands": int(len(work))}
        part = pd.DataFrame(
            {
                "observed": test[response].to_numpy(dtype=float),
                "predicted": prediction,
                "weight": test[weight].to_numpy(dtype=float),
            }
        )
        records.append(part)
    if not records:
        return {"status": "no_valid_folds", "n_islands": int(len(work))}
    scored = pd.concat(records, ignore_index=True)
    residual = scored["observed"] - scored["predicted"]
    weights = scored["weight"]
    wmse = float(np.average(np.square(residual), weights=weights))
    wmae = float(np.average(np.abs(residual), weights=weights))
    return {
        "status": "scored",
        "n_islands": int(len(scored)),
        "n_spatial_blocks": int(work[cluster].nunique()),
        "weighted_rmse": math.sqrt(max(wmse, 0.0)),
        "weighted_mae": wmae,
    }


def fit_northern_bombus_models(
    genus_null: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    strata: tuple[str, ...] = DEFAULT_STRATA,
    cluster: str = "spatial_block",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    required_null = {
        "island_id",
        "outcome",
        "stratum",
        "observed_n_species",
        "deviation_observed_minus_null",
    }
    missing = required_null - set(genus_null.columns)
    if missing:
        raise typer.BadParameter(f"genus-null table missing columns: {sorted(missing)}")
    required_cov = {"island_id", "analysis_regime", cluster, COMPATIBILITY, *BASELINE}
    missing_cov = required_cov - set(covariates.columns)
    if missing_cov:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing_cov)}")

    data = genus_null.merge(
        covariates[list(required_cov)].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    data = data.loc[data["analysis_regime"].eq("northern_midlatitude")].copy()
    data = _add_interaction(data)

    coefficient_tables: list[pd.DataFrame] = []
    cv_rows: list[dict[str, Any]] = []
    support_rows: list[dict[str, Any]] = []
    for stratum in strata:
        outcomes = sorted(data.loc[data["stratum"].eq(stratum), "outcome"].dropna().unique())
        for outcome in outcomes:
            subset = data.loc[
                data["stratum"].eq(stratum) & data["outcome"].eq(outcome)
            ].copy()
            support_rows.append(
                {
                    "stratum": stratum,
                    "outcome": outcome,
                    "n_islands": int(len(subset)),
                    "n_spatial_blocks": int(subset[cluster].dropna().nunique()),
                    "median_direct_trials": (
                        float(subset["observed_n_species"].median()) if len(subset) else np.nan
                    ),
                }
            )
            for model, predictors in MODELS.items():
                coefficients, fit = fit_weighted_linear_clustered(
                    subset,
                    response_column="deviation_observed_minus_null",
                    weight_column="observed_n_species",
                    predictors=predictors,
                    cluster_column=cluster,
                )
                if not coefficients.empty:
                    coefficients.insert(0, "model", model)
                    coefficients.insert(0, "outcome", outcome)
                    coefficients.insert(0, "stratum", stratum)
                    coefficient_tables.append(coefficients)
                cv = spatial_block_cv(subset, predictors=predictors, cluster=cluster)
                cv_rows.append({"stratum": stratum, "outcome": outcome, "model": model, **cv})

    coefficients = (
        pd.concat(coefficient_tables, ignore_index=True)
        if coefficient_tables
        else pd.DataFrame()
    )
    cv = pd.DataFrame(cv_rows)
    if not cv.empty:
        baseline_rmse = (
            cv.loc[cv["model"].eq("M0_lineage_geography_climate"),
                   ["stratum", "outcome", "weighted_rmse"]]
            .rename(columns={"weighted_rmse": "m0_weighted_rmse"})
        )
        cv = cv.merge(baseline_rmse, on=["stratum", "outcome"], how="left")
        cv["rmse_improvement_vs_m0"] = cv["m0_weighted_rmse"] - cv["weighted_rmse"]
    return coefficients, cv, pd.DataFrame(support_rows)


@app.command("run")
def run(
    genus_null_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    coefficients, cv, support = fit_northern_bombus_models(
        pd.read_csv(genus_null_csv), pd.read_csv(covariates_csv)
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    coefficients.to_csv(output_dir / "northern_bombus_residual_coefficients.csv", index=False)
    cv.to_csv(output_dir / "northern_bombus_residual_spatial_cv.csv", index=False)
    support.to_csv(output_dir / "northern_bombus_residual_support.csv", index=False)
    focal = coefficients.loc[
        coefficients["predictor"].isin([COMPATIBILITY, INTERACTION])
    ].copy() if not coefficients.empty else coefficients.copy()
    focal.to_csv(output_dir / "northern_bombus_focal_coefficients.csv", index=False)
    manifest = {
        "contract": "northern_bombus_lineage_residual_v1",
        "regime": "northern_midlatitude",
        "response": "observed-minus-genus-fixed-null direct trait composition",
        "compatibility": COMPATIBILITY,
        "models": MODELS,
        "guardrails": [
            "environmental compatibility is not observed Bombus occurrence or loss",
            "compatibility is climate-derived, so incremental value must exceed flexible climate controls",
            "the distance x compatibility model is an arrival-limitation sensitivity, not causal proof",
            "spatial-block CV is used to judge incremental predictive value",
        ],
    }
    (output_dir / "northern_bombus_residual_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
