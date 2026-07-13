"""Run the preregistered M0-M4 SEM ladder with model-specific applicability scopes.

M0 is descriptive over resolved applicability strata; M1-M3 are restricted to the
predeclared primary Bombus channel; M4 is the compensation model over resolved
weak/absent/not-applicable cases. Floral trait proxies remain proxies and are never
interpreted as realised pollinator observations.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer

from island_v2.m0_m4_sem import _fit_ols, _load_config, build_sem_island_data

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.fillna("").astype(str).str.lower().isin({"true", "1", "yes", "y"})


def select_model_scopes(data: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Return explicit, non-overwritten analysis scopes for M0-M4."""
    work = data.copy()
    applicability = work.get("applicability", pd.Series("unresolved", index=work.index)).fillna("unresolved").astype(str)
    resolved = applicability.isin({"applicable", "structurally_not_applicable"})
    primary = _as_bool(work.get("primary_bombus_channel_analysis", pd.Series(False, index=work.index)))
    bombus = pd.to_numeric(work.get("bombus_channel_state", pd.Series(np.nan, index=work.index)), errors="coerce")

    # M4 follows the preregistered scope: weak/absent Bombus channel or structural non-applicability.
    m4_scope = resolved & (
        applicability.eq("structurally_not_applicable") | bombus.lt(0.5)
    )
    return {
        "M0": work.loc[resolved].copy(),
        "M1": work.loc[primary].copy(),
        "M2": work.loc[primary].copy(),
        "M3": work.loc[primary].copy(),
        "M4": work.loc[m4_scope].copy(),
    }


def fit_scoped_model_ladder(
    data: pd.DataFrame,
    baseline: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
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

    scopes = select_model_scopes(data)
    path_parts: list[pd.DataFrame] = []
    fit_rows: list[dict[str, Any]] = []
    scope_rows: list[dict[str, Any]] = []
    coefficient_lookup: dict[tuple[str, str, str], float] = {}

    for model, model_equations in equations.items():
        model_data = scopes[model].copy()
        model_data["bombus_x_alternative"] = pd.to_numeric(
            model_data["bombus_channel_state"], errors="coerce"
        ) * pd.to_numeric(
            model_data["alternative_functional_replacement_proxy"], errors="coerce"
        )
        scope_rows.append(
            {
                "model": model,
                "n_islands_scope": int(model_data["island_id"].nunique()),
                "scope_rule": {
                    "M0": "resolved applicability strata",
                    "M1": "primary Bombus channel only",
                    "M2": "primary Bombus channel only",
                    "M3": "primary Bombus channel only",
                    "M4": "resolved and (structurally not applicable or Bombus compatibility < 0.5)",
                }[model],
            }
        )
        for equation_index, (outcome, predictors) in enumerate(model_equations, start=1):
            coeffs, fit = _fit_ols(model_data, outcome, predictors)
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
    return paths, pd.DataFrame(fit_rows), pd.DataFrame(indirect_rows), pd.DataFrame(scope_rows)


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

    quality_column = str(config["bombus_channel"].get("extrapolation_quality_column", ""))
    allowed_quality = set(config.get("reporting", {}).get("primary_extrapolation_quality", []))
    if quality_column and quality_column in island_data.columns and allowed_quality:
        quality = island_data[quality_column].fillna("unknown").astype(str)
        applicable = island_data.get("applicability", pd.Series("unresolved", index=island_data.index)).astype(str).eq("applicable")
        # Keep structural comparison islands for M4; quality-filter only the applicable Bombus channel.
        island_data = island_data.loc[(~applicable) | quality.isin(allowed_quality)].copy()

    baseline = metadata["available_baseline_covariates"]
    paths, fit, indirect, scopes = fit_scoped_model_ladder(island_data, baseline)

    output_dir.mkdir(parents=True, exist_ok=True)
    island_data.to_csv(output_dir / "sem_island_data.csv", index=False)
    paths.to_csv(output_dir / "sem_path_coefficients.csv", index=False)
    fit.to_csv(output_dir / "sem_model_fit.csv", index=False)
    indirect.to_csv(output_dir / "sem_indirect_effects.csv", index=False)
    scopes.to_csv(output_dir / "sem_model_scopes.csv", index=False)
    manifest = {
        "contract": "m0_m4_piecewise_sem_scoped_v1",
        "analysis_tier": analysis_tier,
        "n_islands_input": int(island_data["island_id"].nunique()),
        "model_scope_counts": {
            row.model: int(row.n_islands_scope) for row in scopes.itertuples(index=False)
        },
        **metadata,
        "m0_m4_definition": "config/analysis_models.yml",
        "proxy_contract": str(config_path),
        "warning": (
            "SEM paths use island-level trait-composition proxies, not realised pollinator identities. "
            "M0-M4 use model-specific applicability scopes; missing baseline blocks remain explicit."
        ),
    }
    (output_dir / "sem_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
