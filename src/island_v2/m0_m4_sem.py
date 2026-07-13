"""Fit the preregistered M0-M4 island-level piecewise SEM ladder.

This module operates on island-level trait-composition proxies and Bombus channel
metrics. It does not infer realised pollinators from floral traits. Baseline
covariates are used only when supplied; missing baseline blocks are reported in
metadata rather than fabricated.
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

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or "proxies" not in payload or "models" not in payload:
        raise typer.BadParameter("SEM config must contain proxies and models")
    return payload


def _share(table: pd.DataFrame, trait_name: str, positive_values: list[str]) -> pd.DataFrame:
    work = table.loc[table["trait_name"].eq(trait_name)].copy()
    if work.empty:
        return pd.DataFrame(columns=["island_id", "value"])
    work["species_share_within_trait"] = pd.to_numeric(
        work["species_share_within_trait"], errors="coerce"
    )
    mask = work["filled_value"].astype(str).isin([str(v) for v in positive_values])
    result = (
        work.loc[mask]
        .groupby("island_id", sort=True)["species_share_within_trait"]
        .sum(min_count=1)
        .rename("value")
        .reset_index()
    )
    islands = pd.DataFrame({"island_id": sorted(work["island_id"].astype(str).unique())})
    return islands.merge(result, on="island_id", how="left").fillna({"value": 0.0})


def build_sem_island_data(
    composition: pd.DataFrame,
    island_metrics: pd.DataFrame,
    config: dict[str, Any],
    covariates: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    required_comp = {"island_id", "trait_name", "filled_value", "species_share_within_trait"}
    missing = required_comp.difference(composition.columns)
    if missing:
        raise typer.BadParameter(f"trait composition missing columns: {sorted(missing)}")
    if "island_id" not in island_metrics.columns:
        raise typer.BadParameter("island metrics missing island_id")

    proxies = config["proxies"]
    floral = proxies["floral_phenotype_proxy"]
    floral_share = _share(composition, floral["trait_name"], floral["positive_values"]).rename(
        columns={"value": "floral_phenotype_proxy"}
    )

    mediator_parts: list[pd.DataFrame] = []
    for idx, component in enumerate(proxies["reproductive_assurance_proxy"]["components"]):
        mediator_parts.append(
            _share(composition, component["trait_name"], component["positive_values"]).rename(
                columns={"value": f"reproductive_component_{idx}"}
            )
        )
    mediator = mediator_parts[0]
    for part in mediator_parts[1:]:
        mediator = mediator.merge(part, on="island_id", how="outer", validate="one_to_one")
    mediator_cols = [c for c in mediator.columns if c.startswith("reproductive_component_")]
    mediator["reproductive_assurance_proxy"] = mediator[mediator_cols].mean(axis=1, skipna=True)
    mediator = mediator[["island_id", "reproductive_assurance_proxy"]]

    alt = proxies["alternative_functional_replacement_proxy"]
    alternative = _share(composition, alt["trait_name"], alt["positive_values"]).rename(
        columns={"value": "alternative_functional_replacement_proxy"}
    )
    alternative["alternative_functional_replacement_proxy"] = alternative[
        "alternative_functional_replacement_proxy"
    ].clip(0.0, 1.0)

    base = island_metrics.drop_duplicates("island_id").copy()
    result = base.merge(floral_share, on="island_id", how="left", validate="one_to_one")
    result = result.merge(mediator, on="island_id", how="left", validate="one_to_one")
    result = result.merge(alternative, on="island_id", how="left", validate="one_to_one")

    bombus_column = str(config["bombus_channel"]["column"])
    if bombus_column not in result.columns:
        raise typer.BadParameter(f"island metrics missing Bombus column {bombus_column}")
    result["bombus_channel_state"] = pd.to_numeric(result[bombus_column], errors="coerce")

    if covariates is not None:
        if "island_id" not in covariates.columns:
            raise typer.BadParameter("covariate table missing island_id")
        result = result.merge(
            covariates.drop_duplicates("island_id"), on="island_id", how="left", validate="one_to_one"
        )

    candidate_baseline = [str(v) for v in config.get("baseline_candidates", [])]
    available_baseline = [column for column in candidate_baseline if column in result.columns]
    missing_baseline = [column for column in candidate_baseline if column not in result.columns]
    metadata = {
        "baseline_candidates": candidate_baseline,
        "available_baseline_covariates": available_baseline,
        "missing_baseline_covariates": missing_baseline,
        "baseline_complete": not missing_baseline,
    }
    return result, metadata


def _standardize(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    sd = float(values.std(ddof=0))
    if not math.isfinite(sd) or sd == 0:
        return pd.Series(np.zeros(len(values)), index=values.index, dtype=float)
    return (values - float(values.mean())) / sd


def _fit_ols(data: pd.DataFrame, outcome: str, predictors: list[str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    columns = [outcome, *predictors]
    work = data[columns].copy()
    for column in columns:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.dropna()
    n = int(len(work))
    p = int(len(predictors) + 1)
    if n <= p or n < 3:
        return pd.DataFrame(), {"n": n, "r2": np.nan, "aic": np.nan, "status": "insufficient_complete_cases"}

    y = _standardize(work[outcome]).to_numpy(dtype=float)
    xcols = [_standardize(work[pred]).to_numpy(dtype=float) for pred in predictors]
    X = np.column_stack([np.ones(n), *xcols])
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    fitted = X @ beta
    resid = y - fitted
    sse = float(np.sum(resid**2))
    sst = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - sse / sst if sst > 0 else np.nan
    sigma2 = max(sse / n, np.finfo(float).tiny)
    aic = float(n * np.log(sigma2) + 2 * p)

    xtx_inv = np.linalg.pinv(X.T @ X)
    dof = max(n - p, 1)
    mse = sse / dof
    se = np.sqrt(np.diag(xtx_inv) * mse)
    rows = []
    names = ["intercept", *predictors]
    for name, estimate, stderr in zip(names, beta, se, strict=True):
        rows.append(
            {
                "outcome": outcome,
                "predictor": name,
                "standardized_estimate": float(estimate),
                "standard_error": float(stderr),
                "n": n,
            }
        )
    return pd.DataFrame(rows), {"n": n, "r2": float(r2), "aic": aic, "status": "fitted"}


def fit_model_ladder(data: pd.DataFrame, baseline: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    equations: dict[str, list[tuple[str, list[str]]]] = {
        "M0": [("floral_phenotype_proxy", baseline)],
        "M1": [("floral_phenotype_proxy", [*baseline, "bombus_channel_state"])],
        "M2": [
            ("reproductive_assurance_proxy", [*baseline, "bombus_channel_state"]),
            ("floral_phenotype_proxy", [*baseline, "reproductive_assurance_proxy"]),
        ],
        "M3": [
            ("reproductive_assurance_proxy", [*baseline, "bombus_channel_state"]),
            (
                "floral_phenotype_proxy",
                [*baseline, "bombus_channel_state", "reproductive_assurance_proxy"],
            ),
        ],
        "M4": [
            (
                "floral_phenotype_proxy",
                [
                    *baseline,
                    "bombus_channel_state",
                    "alternative_functional_replacement_proxy",
                    "bombus_x_alternative",
                ],
            )
        ],
    }
    work = data.copy()
    work["bombus_x_alternative"] = pd.to_numeric(
        work["bombus_channel_state"], errors="coerce"
    ) * pd.to_numeric(work["alternative_functional_replacement_proxy"], errors="coerce")

    path_parts: list[pd.DataFrame] = []
    fit_rows: list[dict[str, Any]] = []
    coefficient_lookup: dict[tuple[str, str, str], float] = {}
    for model, model_equations in equations.items():
        for equation_index, (outcome, predictors) in enumerate(model_equations, start=1):
            coeffs, fit = _fit_ols(work, outcome, predictors)
            if not coeffs.empty:
                coeffs.insert(0, "equation", equation_index)
                coeffs.insert(0, "model", model)
                path_parts.append(coeffs)
                for row in coeffs.itertuples(index=False):
                    coefficient_lookup[(model, outcome, row.predictor)] = float(row.standardized_estimate)
            fit_rows.append(
                {
                    "model": model,
                    "equation": equation_index,
                    "outcome": outcome,
                    "predictors": " + ".join(predictors) if predictors else "intercept_only",
                    **fit,
                }
            )

    indirect_rows: list[dict[str, Any]] = []
    for model in ("M2", "M3"):
        a = coefficient_lookup.get((model, "reproductive_assurance_proxy", "bombus_channel_state"))
        b = coefficient_lookup.get((model, "floral_phenotype_proxy", "reproductive_assurance_proxy"))
        indirect_rows.append(
            {
                "model": model,
                "path": "bombus_channel_state -> reproductive_assurance_proxy -> floral_phenotype_proxy",
                "indirect_effect": float(a * b) if a is not None and b is not None else np.nan,
            }
        )

    paths = pd.concat(path_parts, ignore_index=True) if path_parts else pd.DataFrame()
    return paths, pd.DataFrame(fit_rows), pd.DataFrame(indirect_rows)


@app.command("run")
def run(
    trait_composition_csv: Path = typer.Option(..., exists=True),
    island_metrics_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/m0_m4_sem.yml"), exists=True),
    covariates_csv: Path | None = typer.Option(None),
    analysis_tier: str = typer.Option("primary"),
) -> None:
    config = _load_config(config_path)
    composition = pd.read_csv(trait_composition_csv)
    metrics = pd.read_csv(island_metrics_csv)
    covariates = pd.read_csv(covariates_csv) if covariates_csv else None
    island_data, metadata = build_sem_island_data(composition, metrics, config, covariates)

    gate_column = str(config["bombus_channel"]["applicability_gate_column"])
    if gate_column in island_data.columns:
        gate = island_data[gate_column].fillna(False).astype(bool)
        primary_data = island_data.loc[gate].copy()
    else:
        primary_data = island_data.copy()

    quality_column = str(config["bombus_channel"].get("extrapolation_quality_column", ""))
    allowed_quality = set(config.get("reporting", {}).get("primary_extrapolation_quality", []))
    if quality_column and quality_column in primary_data.columns and allowed_quality:
        primary_data = primary_data.loc[
            primary_data[quality_column].fillna("unknown").astype(str).isin(allowed_quality)
        ].copy()

    baseline = metadata["available_baseline_covariates"]
    paths, fit, indirect = fit_model_ladder(primary_data, baseline)

    output_dir.mkdir(parents=True, exist_ok=True)
    island_data.to_csv(output_dir / "sem_island_data.csv", index=False)
    paths.to_csv(output_dir / "sem_path_coefficients.csv", index=False)
    fit.to_csv(output_dir / "sem_model_fit.csv", index=False)
    indirect.to_csv(output_dir / "sem_indirect_effects.csv", index=False)
    manifest = {
        "contract": "m0_m4_piecewise_sem_v1",
        "analysis_tier": analysis_tier,
        "n_islands_input": int(island_data["island_id"].nunique()),
        "n_islands_primary_sem": int(primary_data["island_id"].nunique()),
        **metadata,
        "m0_m4_definition": "config/analysis_models.yml",
        "proxy_contract": str(config_path),
        "warning": (
            "SEM paths use island-level trait-composition proxies. They do not identify realised pollinators. "
            "Any missing baseline blocks are reported explicitly and must not be interpreted as adjusted effects."
        ),
    }
    (output_dir / "sem_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
