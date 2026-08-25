"""Trait-resolution bias audit and coverage-adjusted WHEN/WHERE sensitivity.

This layer separates two questions:

1. Among opportunistically observed species in a floristic-status stratum, what
   fraction has direct evidence for each Chapter 1 trait domain, and does that
   resolution probability vary with geography?
2. Do the northern-midlatitude WHERE, tropical WHERE, and north-versus-tropical
   BETWEEN-WHERE results persist when each atomic response includes its own
   trait-resolution coverage covariate?

The observed flora is not assumed complete. Missing traits/species are never recoded
as biological trait zeros.
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

from island_v2.chapter1_context_analysis import (
    _chi_square_sf_integer_df,
    _fit_grouped_binomial_design,
)

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


def _bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.fillna("").astype(str).str.lower().isin({"true", "1", "yes"})


def _bh(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=values.index, dtype=float)
    ok = p.notna()
    if not ok.any():
        return out
    x = p.loc[ok].to_numpy(float)
    order = np.argsort(x)
    ranked = x[order]
    n = len(ranked)
    adjusted = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    restored = np.empty(n, dtype=float)
    restored[order] = np.clip(adjusted, 0, 1)
    out.loc[ok] = restored
    return out


def _joint_wald(beta: np.ndarray, covariance: np.ndarray) -> tuple[float, int, float]:
    rank = int(np.linalg.matrix_rank(covariance))
    if rank <= 0:
        return float("nan"), 0, float("nan")
    stat = float(beta @ np.linalg.pinv(covariance) @ beta)
    return stat, rank, _chi_square_sf_integer_df(stat, rank)


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
    resolved["resolved_for_primary"] = _bool(resolved["resolved_for_primary"])
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
                rows.append(
                    {
                        "island_id": island_id,
                        "stratum": stratum,
                        "trait_name": trait,
                        "n_observed_stratum_species": total,
                        "n_direct_trait_species": covered,
                        "direct_trait_fraction": covered / total,
                        "coverage_logit_smoothed": float(
                            math.log((covered + 0.5) / (total - covered + 0.5))
                        ),
                    }
                )
    return pd.DataFrame(rows)


def _standardize(values: pd.Series) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant predictor")
    return (x - mean) / sd


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
        reference_here = reference if reference in contexts else contexts[0]
        contrasts = [x for x in contexts if x != reference_here]

        names = ["intercept"]
        cols = [np.ones(len(work), dtype=float)]
        for predictor in baseline:
            names.append(f"z_{predictor}")
            cols.append(_standardize(work[predictor]))
        names.append(f"z_{gradient}")
        z_distance = _standardize(work[gradient])
        cols.append(z_distance)
        richness = np.log1p(work["n_observed_stratum_species"].to_numpy(float))
        names.append("z_log1p_observed_stratum_species")
        mean = float(richness.mean())
        sd = float(richness.std(ddof=0))
        cols.append((richness - mean) / sd)

        for ctx in contrasts:
            dummy = work[context_col].eq(ctx).astype(float).to_numpy()
            names.append(f"context[{ctx}]")
            cols.append(dummy)
            names.append(f"z_{gradient}:context[{ctx}]")
            cols.append(z_distance * dummy)

        coef, fit, _ = _fit_grouped_binomial_design(
            work["n_direct_trait_species"].to_numpy(float),
            work["n_observed_stratum_species"].to_numpy(float),
            np.column_stack(cols),
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
                    "median_direct_trait_fraction": float(group["direct_trait_fraction"].median()),
                    "mean_direct_trait_fraction": float(group["direct_trait_fraction"].mean()),
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


def _prepare_adjusted(
    composition: pd.DataFrame,
    coverage: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    stratum: str,
    config: dict[str, Any],
) -> pd.DataFrame:
    baseline = [str(x) for x in config["baseline_covariates"]]
    gradient = str(config["distance_gradient"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    contexts = set(str(x) for x in config["contexts"])

    comp = composition.loc[composition["stratum"].astype(str).eq(stratum)].copy()
    comp["response_id"] = comp["trait_name"].astype(str) + "::" + comp["trait_state"].astype(str)
    comp = comp.merge(
        coverage[
            ["island_id", "stratum", "trait_name", "coverage_logit_smoothed"]
        ],
        on=["island_id", "stratum", "trait_name"],
        how="left",
        validate="many_to_one",
    )
    comp = comp.merge(
        covariates[["island_id", gradient, context, cluster, *baseline]],
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    numeric = ["successes", "trials", "coverage_logit_smoothed", gradient, *baseline]
    for col in numeric:
        comp[col] = pd.to_numeric(comp[col], errors="coerce")
    comp[context] = comp[context].fillna("").astype(str)
    comp[cluster] = comp[cluster].fillna("").astype(str)
    comp = comp.dropna(subset=numeric)
    return comp.loc[
        comp[context].isin(contexts)
        & comp[cluster].ne("")
        & comp["trials"].gt(0)
        & comp["successes"].ge(0)
        & comp["successes"].le(comp["trials"])
    ].copy()


def _response_specific_z(
    work: pd.DataFrame,
    response_mask: np.ndarray,
    predictor: str,
) -> np.ndarray:
    values = work.loc[response_mask, predictor]
    z = _standardize(values)
    full = np.zeros(len(work), dtype=float)
    full[response_mask] = z
    return full


def _within_adjusted(
    rows: pd.DataFrame,
    context_name: str,
    *,
    config: dict[str, Any],
    threshold: int,
) -> dict[str, Any] | None:
    baseline = [str(x) for x in config["baseline_covariates"]]
    gradient = str(config["distance_gradient"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])

    work = rows.loc[rows[context].eq(context_name)].copy()
    counts = work.groupby("response_id")["island_id"].nunique()
    responses = sorted(counts.loc[counts >= threshold].index)
    work = work.loc[work["response_id"].isin(responses)].copy()
    if len(responses) < 2:
        return None

    names: list[str] = []
    cols: list[np.ndarray] = []
    iso_names: list[str] = []
    for response in responses:
        response_mask = work["response_id"].eq(response).to_numpy()
        mask = response_mask.astype(float)
        names.append(f"response[{response}]")
        cols.append(mask)
        for predictor in [*baseline, "coverage_logit_smoothed"]:
            names.append(f"response[{response}]:z_{predictor}")
            cols.append(_response_specific_z(work, response_mask, predictor))
        iso_name = f"response[{response}]:z_{gradient}"
        names.append(iso_name)
        cols.append(_response_specific_z(work, response_mask, gradient))
        iso_names.append(iso_name)

    coef, fit, covariance = _fit_grouped_binomial_design(
        work["successes"].to_numpy(float),
        work["trials"].to_numpy(float),
        np.column_stack(cols),
        names,
        work[cluster].to_numpy(str),
    )
    beta = coef.set_index("predictor")["estimate_log_odds"]
    indices = [names.index(x) for x in iso_names]
    vector = np.array([float(beta[x]) for x in iso_names])
    stat, df, p = _joint_wald(vector, covariance[np.ix_(indices, indices)])
    return {
        "context": context_name,
        "n_responses": len(responses),
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "joint_wald_chisq": stat,
        "joint_df": df,
        "p_value": p,
    }


def _between_adjusted(
    rows: pd.DataFrame,
    context_a: str,
    context_b: str,
    *,
    config: dict[str, Any],
    threshold: int,
) -> dict[str, Any] | None:
    baseline = [str(x) for x in config["baseline_covariates"]]
    gradient = str(config["distance_gradient"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])

    work = rows.loc[rows[context].isin([context_a, context_b])].copy()
    counts = work.groupby(["response_id", context])["island_id"].nunique().unstack(fill_value=0)
    if context_a not in counts.columns or context_b not in counts.columns:
        return None
    responses = sorted(
        counts.index[counts[context_a].ge(threshold) & counts[context_b].ge(threshold)]
    )
    work = work.loc[work["response_id"].isin(responses)].copy()
    if len(responses) < 2:
        return None

    names: list[str] = []
    cols: list[np.ndarray] = []
    interaction_names: list[str] = []
    ctx_indicator = work[context].eq(context_b).to_numpy(float)
    for response in responses:
        response_mask = work["response_id"].eq(response).to_numpy()
        mask = response_mask.astype(float)
        names.append(f"response[{response}]")
        cols.append(mask)
        for predictor in [*baseline, "coverage_logit_smoothed"]:
            names.append(f"response[{response}]:z_{predictor}")
            cols.append(_response_specific_z(work, response_mask, predictor))
        names.append(f"response[{response}]:context[{context_b}]")
        cols.append(mask * ctx_indicator)
        full_iso = _response_specific_z(work, response_mask, gradient)
        names.append(f"response[{response}]:z_{gradient}")
        cols.append(full_iso)
        interaction = f"response[{response}]:z_{gradient}:context[{context_b}]"
        names.append(interaction)
        cols.append(full_iso * ctx_indicator)
        interaction_names.append(interaction)

    coef, fit, covariance = _fit_grouped_binomial_design(
        work["successes"].to_numpy(float),
        work["trials"].to_numpy(float),
        np.column_stack(cols),
        names,
        work[cluster].to_numpy(str),
    )
    beta = coef.set_index("predictor")["estimate_log_odds"]
    indices = [names.index(x) for x in interaction_names]
    vector = np.array([float(beta[x]) for x in interaction_names])
    stat, df, p = _joint_wald(vector, covariance[np.ix_(indices, indices)])
    return {
        "context_a": context_a,
        "context_b": context_b,
        "n_responses": len(responses),
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "joint_wald_chisq": stat,
        "joint_df": df,
        "p_value": p,
    }


def run_coverage_adjusted_omnibus(
    composition: pd.DataFrame,
    coverage: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    threshold = int(config["confirmatory_islands_per_response"])
    contexts = [str(x) for x in config["headline_contexts"]]
    strata = [str(x) for x in config["headline_strata"]]
    within_rows: list[dict[str, Any]] = []
    between_rows: list[dict[str, Any]] = []

    for stratum in strata:
        rows = _prepare_adjusted(
            composition, coverage, covariates, stratum=stratum, config=config
        )
        for ctx in contexts:
            result = _within_adjusted(rows, ctx, config=config, threshold=threshold)
            if result is not None:
                within_rows.append({"stratum": stratum, **result})
        if len(contexts) >= 2:
            result = _between_adjusted(
                rows, contexts[0], contexts[1], config=config, threshold=threshold
            )
            if result is not None:
                between_rows.append({"stratum": stratum, **result})

    within = pd.DataFrame(within_rows)
    between = pd.DataFrame(between_rows)
    if not within.empty:
        within["q_value"] = within.groupby("stratum", group_keys=False)["p_value"].apply(_bh)
        within["supported"] = within["q_value"].le(0.05)
    if not between.empty:
        between["q_value"] = _bh(between["p_value"])
        between["supported"] = between["q_value"].le(0.05)

    summary_rows = []
    for stratum in strata:
        w = within.loc[within["stratum"].eq(stratum)]
        b = between.loc[between["stratum"].eq(stratum)]
        north = w.loc[w["context"].eq(contexts[0])]
        tropical = w.loc[w["context"].eq(contexts[1])] if len(contexts) > 1 else pd.DataFrame()
        summary_rows.append(
            {
                "stratum": stratum,
                "north_supported": bool(not north.empty and north.iloc[0]["supported"]),
                "tropical_supported": bool(not tropical.empty and tropical.iloc[0]["supported"]),
                "between_supported": bool(not b.empty and b.iloc[0]["supported"]),
                "headline_replication": bool(
                    not north.empty
                    and not tropical.empty
                    and not b.empty
                    and north.iloc[0]["supported"]
                    and tropical.iloc[0]["supported"]
                    and b.iloc[0]["supported"]
                ),
            }
        )
    return within, between, pd.DataFrame(summary_rows)


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
    within, between, summary = run_coverage_adjusted_omnibus(
        composition, coverage, covariates, config
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    coverage.to_csv(output_dir / "trait_resolution_coverage.csv.gz", index=False)
    coefficients.to_csv(output_dir / "trait_resolution_selection_coefficients.csv", index=False)
    fits.to_csv(output_dir / "trait_resolution_selection_fit.csv", index=False)
    support.to_csv(output_dir / "trait_resolution_support_by_context.csv", index=False)
    within.to_csv(output_dir / "coverage_adjusted_when_where_within.csv", index=False)
    between.to_csv(output_dir / "coverage_adjusted_when_where_between.csv", index=False)
    summary.to_csv(output_dir / "coverage_adjusted_headline_summary.csv", index=False)

    manifest = {
        "contract": "chapter1_trait_resolution_bias_v2",
        "observed_flora_assumed_complete": False,
        "missing_trait_is_zero": False,
        "n_coverage_rows": int(len(coverage)),
        "coverage_adjustment": "trait-specific smoothed-logit coverage with response-specific coefficient",
        "headline_replications": int(summary["headline_replication"].fillna(False).sum()),
        "headline_opportunities": int(len(summary)),
    }
    (output_dir / "chapter1_trait_resolution_bias_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
