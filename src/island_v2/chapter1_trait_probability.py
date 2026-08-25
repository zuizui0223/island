"""Trait-centric probability surfaces for Chapter 1.

The estimand is deliberately conditional on the directly trait-resolved part of the
observed realised island flora. For each floristic-status stratum, biogeographic
context and atomic trait state, fit

    P(trait state | directly trait-resolved observed flora,
      mainland-distance/source-pool-accessibility gradient, area, climate).

This is an assembly/composition probability, not the probability that a trait evolved
on an island and not a species-occurrence probability relative to a complete source
pool. Islands/species without direct trait evidence are unresolved here; the full
island observation process is modelled separately by chapter1_observation_bias.py.
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

from island_v2.chapter1_context_analysis import _fit_grouped_binomial_design

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _expit(value: np.ndarray | float) -> np.ndarray:
    x = np.asarray(value, dtype=float)
    return 1.0 / (1.0 + np.exp(-np.clip(x, -35.0, 35.0)))


def _z(values: pd.Series) -> tuple[np.ndarray, float, float]:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid continuous predictor")
    return (x - mean) / sd, mean, sd


def _fit_curve(
    rows: pd.DataFrame,
    *,
    baseline: list[str],
    distance: str,
    cluster: str,
    quantiles: list[float],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    work = rows.copy()
    numeric = ["successes", "trials", distance, *baseline]
    for column in numeric:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work[cluster] = work[cluster].fillna("").astype(str)
    work = work.dropna(subset=numeric)
    work = work.loc[
        (work["trials"] > 0)
        & (work["successes"] >= 0)
        & (work["successes"] <= work["trials"])
        & work[cluster].ne("")
    ].copy()
    if work.empty:
        return pd.DataFrame(), pd.DataFrame(), {"status": "empty"}

    names = ["intercept"]
    columns = [np.ones(len(work), dtype=float)]
    scaling: dict[str, tuple[float, float]] = {}
    for predictor in [*baseline, distance]:
        z, mean, sd = _z(work[predictor])
        scaling[predictor] = (mean, sd)
        names.append(f"z_{predictor}")
        columns.append(z)
    X = np.column_stack(columns)
    coef, fit, covariance = _fit_grouped_binomial_design(
        work["successes"].to_numpy(float),
        work["trials"].to_numpy(float),
        X,
        names,
        work[cluster].to_numpy(str),
    )
    coef = coef.copy()
    coef["odds_ratio"] = np.exp(coef["estimate_log_odds"])

    beta = coef.set_index("predictor")["estimate_log_odds"]
    intercept = float(beta["intercept"])
    distance_name = f"z_{distance}"
    distance_beta = float(beta[distance_name])
    intercept_idx = names.index("intercept")
    distance_idx = names.index(distance_name)
    distance_values = pd.to_numeric(work[distance], errors="coerce")

    prediction_rows: list[dict[str, Any]] = []
    for q in quantiles:
        raw = float(distance_values.quantile(q))
        mean, sd = scaling[distance]
        z_distance = (raw - mean) / sd
        # Baseline covariates are held at their within-fit means (z=0).
        linear = intercept + distance_beta * z_distance
        design = np.zeros(len(names), dtype=float)
        design[intercept_idx] = 1.0
        design[distance_idx] = z_distance
        variance = float(design @ covariance @ design)
        se_linear = math.sqrt(max(variance, 0.0))
        prediction_rows.append(
            {
                "distance_quantile": float(q),
                "distance_gradient_value": raw,
                "distance_gradient_z": float(z_distance),
                "predicted_probability": float(_expit(linear)),
                "predicted_probability_lcl": float(_expit(linear - 1.96 * se_linear)),
                "predicted_probability_ucl": float(_expit(linear + 1.96 * se_linear)),
            }
        )

    metadata = {
        **fit,
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_trait_resolved_species_equivalents": int(work["trials"].sum()),
        "distance_mean": scaling[distance][0],
        "distance_sd": scaling[distance][1],
        "distance_slope_log_odds_per_sd": distance_beta,
        "distance_slope_or_per_sd": float(math.exp(distance_beta)),
    }
    return coef, pd.DataFrame(prediction_rows), metadata


def run_trait_probability(
    composition: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    required_comp = {
        "island_id",
        "stratum",
        "trait_name",
        "trait_state",
        "successes",
        "trials",
    }
    missing = required_comp - set(composition.columns)
    if missing:
        raise typer.BadParameter(f"composition missing columns: {sorted(missing)}")

    baseline = [str(x) for x in config["baseline_covariates"]]
    distance = str(config["distance_gradient"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    contexts = [str(x) for x in config["contexts"]]
    strata = [str(x) for x in config["strata"]]
    min_islands = int(config["min_islands_per_curve"])
    confirmatory = int(config["confirmatory_islands_per_curve"])
    quantiles = [float(x) for x in config["distance_grid_quantiles"]]

    needed_cov = {"island_id", context, cluster, distance, *baseline}
    missing_cov = needed_cov - set(covariates.columns)
    if missing_cov:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing_cov)}")

    joined = composition.merge(
        covariates[["island_id", context, cluster, distance, *baseline]],
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    joined[context] = joined[context].fillna("").astype(str)

    coefficient_parts: list[pd.DataFrame] = []
    probability_parts: list[pd.DataFrame] = []
    fit_rows: list[dict[str, Any]] = []

    category_keys = (
        joined[["stratum", "trait_name", "trait_state"]]
        .drop_duplicates()
        .sort_values(["stratum", "trait_name", "trait_state"])
    )
    for key in category_keys.itertuples(index=False):
        stratum = str(key.stratum)
        trait_name = str(key.trait_name)
        trait_state = str(key.trait_state)
        if stratum not in strata:
            continue
        category = joined.loc[
            joined["stratum"].astype(str).eq(stratum)
            & joined["trait_name"].astype(str).eq(trait_name)
            & joined["trait_state"].astype(str).eq(trait_state)
        ].copy()
        for focal_context in contexts:
            rows = category.loc[category[context].eq(focal_context)].copy()
            n_islands = int(rows["island_id"].nunique())
            if n_islands < min_islands:
                fit_rows.append(
                    {
                        "stratum": stratum,
                        "trait_name": trait_name,
                        "trait_state": trait_state,
                        "context": focal_context,
                        "status": "below_min_islands",
                        "support_tier": "not_testable",
                        "n_unique_islands": n_islands,
                    }
                )
                continue
            try:
                coef, probabilities, metadata = _fit_curve(
                    rows,
                    baseline=baseline,
                    distance=distance,
                    cluster=cluster,
                    quantiles=quantiles,
                )
            except (ValueError, np.linalg.LinAlgError) as exc:
                fit_rows.append(
                    {
                        "stratum": stratum,
                        "trait_name": trait_name,
                        "trait_state": trait_state,
                        "context": focal_context,
                        "status": f"fit_error:{type(exc).__name__}",
                        "support_tier": "not_testable",
                        "n_unique_islands": n_islands,
                    }
                )
                continue
            tier = "confirmatory" if n_islands >= confirmatory else "pilot"
            if not coef.empty:
                coef.insert(0, "context", focal_context)
                coef.insert(0, "trait_state", trait_state)
                coef.insert(0, "trait_name", trait_name)
                coef.insert(0, "stratum", stratum)
                coefficient_parts.append(coef)
            if not probabilities.empty:
                probabilities.insert(0, "context", focal_context)
                probabilities.insert(0, "trait_state", trait_state)
                probabilities.insert(0, "trait_name", trait_name)
                probabilities.insert(0, "stratum", stratum)
                probabilities["support_tier"] = tier
                probability_parts.append(probabilities)
            fit_rows.append(
                {
                    "stratum": stratum,
                    "trait_name": trait_name,
                    "trait_state": trait_state,
                    "context": focal_context,
                    "support_tier": tier,
                    **metadata,
                }
            )

    coefficients = (
        pd.concat(coefficient_parts, ignore_index=True)
        if coefficient_parts
        else pd.DataFrame()
    )
    probabilities = (
        pd.concat(probability_parts, ignore_index=True)
        if probability_parts
        else pd.DataFrame()
    )
    fits = pd.DataFrame(fit_rows)
    return coefficients, probabilities, fits


@app.command("run")
def run(
    composition_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(
        Path("config/chapter1_trait_probability.yml"), exists=True
    ),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    coefficients, probabilities, fits = run_trait_probability(
        pd.read_csv(composition_csv), pd.read_csv(covariates_csv), config
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    coefficients.to_csv(output_dir / "trait_probability_coefficients.csv", index=False)
    probabilities.to_csv(output_dir / "trait_probability_curves.csv", index=False)
    fits.to_csv(output_dir / "trait_probability_fit_support.csv", index=False)

    confirmatory = fits.loc[fits.get("support_tier", pd.Series(dtype=str)).eq("confirmatory")]
    manifest = {
        "contract": "chapter1_trait_probability_v1",
        "estimand": "P(trait_state | directly_trait_resolved_observed_realised_flora, geography)",
        "n_fitted_curves": int(fits["status"].eq("fit").sum()) if "status" in fits else 0,
        "n_confirmatory_curves": int(len(confirmatory)),
        "missing_trait_is_zero": False,
        "missing_flora_is_zero": False,
        "distance_gradient_interpretation": (
            "dispersal limitation plus mainland/source-pool accessibility"
        ),
        "source_pool_occurrence_probability_claimed": False,
        "pollinator_predictors": False,
    }
    (output_dir / "chapter1_trait_probability_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
