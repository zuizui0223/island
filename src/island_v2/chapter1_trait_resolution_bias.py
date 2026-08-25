"""Trait-resolution bias audit and coverage-adjusted WHEN/WHERE sensitivity.

This layer separates two questions:

1. Among observed species in a floristic-status stratum, what fraction has direct
   evidence for each Chapter 1 trait domain, and does that resolution probability
   vary with geography?
2. Do the headline response-vector results persist after adding a trait-specific
   resolution-coverage term to each atomic response model?

The observed flora denominator is not assumed complete. This layer addresses the
second-stage evidence-selection process conditional on the opportunistically observed
flora; no missing species or traits are recoded as biological zeros.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_context_analysis import _fit_grouped_binomial_design
from island_v2.chapter1_when_where_omnibus import run_when_where_omnibus

app = typer.Typer(add_completion=False, no_args_is_help=True)

TRAITS = ("flower_primary_color", "floral_form", "self_incompatibility")
STRATA = ("all_native", "native_nonendemic", "endemic")


def _stratum_mask(frame: pd.DataFrame, stratum: str) -> pd.Series:
    if stratum == "all_native":
        return frame["origin_status"].astype(str).eq("native")
    if stratum == "native_nonendemic":
        return frame["floristic_status"].astype(str).eq("native_nonendemic")
    if stratum == "endemic":
        return frame["floristic_status"].astype(str).eq("endemic")
    raise typer.BadParameter(f"unknown stratum: {stratum}")


def build_trait_resolution_coverage(
    status_flora: pd.DataFrame,
    trait_audit: pd.DataFrame,
) -> pd.DataFrame:
    required_flora = {
        "island_id",
        "accepted_species",
        "origin_status",
        "floristic_status",
    }
    required_audit = {"accepted_species", "trait_name", "resolved_for_primary"}
    miss_flora = required_flora - set(status_flora.columns)
    miss_audit = required_audit - set(trait_audit.columns)
    if miss_flora:
        raise typer.BadParameter(f"status flora missing columns: {sorted(miss_flora)}")
    if miss_audit:
        raise typer.BadParameter(f"trait audit missing columns: {sorted(miss_audit)}")

    resolved = trait_audit.copy()
    resolved["resolved_for_primary"] = resolved["resolved_for_primary"].astype(str).str.lower().isin(
        {"true", "1", "yes"}
    )
    resolved = resolved.loc[
        resolved["resolved_for_primary"]
        & resolved["trait_name"].astype(str).isin(TRAITS)
    ]
    direct_species = {
        trait: set(
            resolved.loc[resolved["trait_name"].astype(str).eq(trait), "accepted_species"]
            .dropna()
            .astype(str)
        )
        for trait in TRAITS
    }

    rows: list[dict[str, Any]] = []
    for stratum in STRATA:
        flora = status_flora.loc[_stratum_mask(status_flora, stratum)].copy()
        for island_id, group in flora.groupby("island_id", sort=False):
            species = set(group["accepted_species"].dropna().astype(str))
            total = len(species)
            if total <= 0:
                continue
            for trait in TRAITS:
                covered = len(species & direct_species[trait])
                fraction = covered / total
                rows.append(
                    {
                        "island_id": island_id,
                        "stratum": stratum,
                        "trait_name": trait,
                        "n_observed_stratum_species": total,
                        "n_direct_trait_species": covered,
                        "direct_trait_fraction": fraction,
                        "coverage_logit_smoothed": float(
                            np.log((covered + 0.5) / (total - covered + 0.5))
                        ),
                    }
                )
    return pd.DataFrame(rows)


def fit_resolution_selection(
    coverage: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    baseline = [str(x) for x in config["baseline_covariates"]]
    gradient = str(config["distance_gradient"])
    context_col = str(config["context_column"])
    cluster_col = str(config["cluster_column"])
    reference = str(config["reference_context"])

    joined = coverage.merge(covariates, on="island_id", how="left", validate="many_to_one")
    coef_parts: list[pd.DataFrame] = []
    fit_rows: list[dict[str, Any]] = []
    support_rows: list[dict[str, Any]] = []

    for (stratum, trait), work in joined.groupby(["stratum", "trait_name"], sort=True):
        numeric = [
            "n_direct_trait_species",
            "n_observed_stratum_species",
            gradient,
            *baseline,
        ]
        work = work.copy()
        for col in numeric:
            work[col] = pd.to_numeric(work[col], errors="coerce")
        work[context_col] = work[context_col].fillna("").astype(str)
        work[cluster_col] = work[cluster_col].fillna("").astype(str)
        work = work.dropna(subset=numeric)
        work = work.loc[
            work["n_observed_stratum_species"].gt(0)
            & work[cluster_col].ne("")
            & work[context_col].ne("")
        ].copy()
        if len(work) < 30:
            continue

        contexts = sorted(work[context_col].unique())
        if reference not in contexts:
            reference_here = contexts[0]
        else:
            reference_here = reference
        contrasts = [x for x in contexts if x != reference_here]

        names = ["intercept"]
        cols = [np.ones(len(work))]
        for predictor in [*baseline, gradient, "log1p_observed_stratum_species"]:
            if predictor == "log1p_observed_stratum_species":
                values = np.log1p(work["n_observed_stratum_species"].to_numpy(float))
            else:
                values = work[predictor].to_numpy(float)
            mean = float(values.mean())
            sd = float(values.std(ddof=0))
            if not np.isfinite(sd) or sd <= 0:
                continue
            names.append(f"z_{predictor}")
            cols.append((values - mean) / sd)
        distance_index = names.index(f"z_{gradient}")
        z_distance = cols[distance_index]
        for ctx in contrasts:
            dummy = work[context_col].eq(ctx).astype(float).to_numpy()
            names.append(f"context[{ctx}]")
            cols.append(dummy)
            names.append(f"z_{gradient}:context[{ctx}]")
            cols.append(z_distance * dummy)

        X = np.column_stack(cols)
        coef, fit, _ = _fit_grouped_binomial_design(
            work["n_direct_trait_species"].to_numpy(float),
            work["n_observed_stratum_species"].to_numpy(float),
            X,
            names,
            work[cluster_col].to_numpy(str),
        )
        coef.insert(0, "trait_name", trait)
        coef.insert(0, "stratum", stratum)
        coef["odds_ratio"] = np.exp(coef["estimate_log_odds"])
        coef_parts.append(coef)
        fit_rows.append({"stratum": stratum, "trait_name": trait, **fit})

        for ctx, group in work.groupby(context_col):
            support_rows.append(
                {
                    "stratum": stratum,
                    "trait_name": trait,
                    "context": ctx,
                    "n_islands": int(group["island_id"].nunique()),
                    "median_direct_trait_fraction": float(
                        group["direct_trait_fraction"].median()
                    ),
                    "mean_direct_trait_fraction": float(
                        group["direct_trait_fraction"].mean()
                    ),
                    "median_observed_stratum_species": float(
                        group["n_observed_stratum_species"].median()
                    ),
                }
            )

    return (
        pd.concat(coef_parts, ignore_index=True) if coef_parts else pd.DataFrame(),
        pd.DataFrame(fit_rows),
        pd.DataFrame(support_rows),
    )


def run_coverage_adjusted_omnibus(
    composition: pd.DataFrame,
    coverage: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    adjusted = composition.merge(
        coverage[
            [
                "island_id",
                "stratum",
                "trait_name",
                "coverage_logit_smoothed",
            ]
        ],
        on=["island_id", "stratum", "trait_name"],
        how="left",
        validate="many_to_one",
    )
    if adjusted["coverage_logit_smoothed"].isna().any():
        raise RuntimeError("coverage-adjusted omnibus has missing trait-resolution coverage")

    cov = covariates.merge(
        adjusted[
            ["island_id", "stratum", "trait_name", "coverage_logit_smoothed"]
        ].drop_duplicates(),
        on="island_id",
        how="right",
    )
    # run_when_where_omnibus joins covariates on island only; trait-specific coverage
    # therefore cannot be passed as an island-level covariate directly. Encode it in
    # composition via a unique row key and use a dedicated copy of the omnibus logic
    # would be unnecessarily fragile. Instead construct a composite island id that
    # preserves each island × stratum × trait-specific coverage while retaining the
    # original spatial block for cluster-robust covariance.
    adjusted = adjusted.copy()
    adjusted["original_island_id"] = adjusted["island_id"]
    adjusted["island_id"] = (
        adjusted["island_id"].astype(str)
        + "||"
        + adjusted["stratum"].astype(str)
        + "||"
        + adjusted["trait_name"].astype(str)
    )
    cov_trait = adjusted[
        ["island_id", "original_island_id", "coverage_logit_smoothed"]
    ].drop_duplicates("island_id")
    cov_trait = cov_trait.merge(
        covariates,
        left_on="original_island_id",
        right_on="island_id",
        how="left",
        suffixes=("", "_original"),
        validate="many_to_one",
    )
    cov_trait = cov_trait.drop(columns=["island_id_original"]).rename(
        columns={"island_id": "composite_island_id"}
    )
    cov_trait = cov_trait.rename(columns={"composite_island_id": "island_id"})

    omnibus_config = dict(config["omnibus"])
    omnibus_config["baseline_covariates"] = [
        *[str(x) for x in omnibus_config["baseline_covariates"]],
        "coverage_logit_smoothed",
    ]
    within, between, _ = run_when_where_omnibus(
        adjusted, cov_trait, omnibus_config
    )
    return within, between


@app.command("run")
def run(
    status_flora_csv: Path = typer.Option(..., exists=True),
    trait_audit_csv: Path = typer.Option(..., exists=True),
    composition_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(
        Path("config/chapter1_trait_resolution_bias.yml"), exists=True
    ),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    status_flora = pd.read_csv(status_flora_csv)
    trait_audit = pd.read_csv(trait_audit_csv)
    composition = pd.read_csv(composition_csv)
    covariates = pd.read_csv(covariates_csv)

    coverage = build_trait_resolution_coverage(status_flora, trait_audit)
    coefficients, fits, support = fit_resolution_selection(coverage, covariates, config)
    within, between = run_coverage_adjusted_omnibus(
        composition, coverage, covariates, config
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    coverage.to_csv(output_dir / "trait_resolution_coverage.csv.gz", index=False)
    coefficients.to_csv(output_dir / "trait_resolution_selection_coefficients.csv", index=False)
    fits.to_csv(output_dir / "trait_resolution_selection_fit.csv", index=False)
    support.to_csv(output_dir / "trait_resolution_support_by_context.csv", index=False)
    within.to_csv(output_dir / "coverage_adjusted_when_where_within.csv", index=False)
    between.to_csv(output_dir / "coverage_adjusted_when_where_between.csv", index=False)

    manifest = {
        "contract": "chapter1_trait_resolution_bias_v1",
        "observed_flora_assumed_complete": False,
        "missing_trait_is_zero": False,
        "n_coverage_rows": int(len(coverage)),
        "coverage_adjustment": "trait_specific_smoothed_logit_as_response_specific_baseline_covariate",
    }
    (output_dir / "chapter1_trait_resolution_bias_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
