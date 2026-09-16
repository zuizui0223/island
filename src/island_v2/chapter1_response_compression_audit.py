"""Audit whether the two-family summary preserves the six-axis response difference.

The all-data probability analysis estimates a six-component North--Tropical distance
interaction vector. This module decomposes that frozen six-dimensional estimand into:

1. the full six-axis vector;
2. two equal-weight family means (accessibility/generalisation and reproductive assurance);
3. four within-family contrasts that are discarded by the two-family compression.

This is a linear-contrast audit of the same fitted model and same island support. It
introduces no new biological outcome definition after seeing the signs.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_all_data_probability import (
    _assemble_cluster_covariance,
    _chi_square_sf_integer_df,
    _fit_single_beta_binomial,
    _prepare,
    _standardize,
    build_broad_counts,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _wald(vector: np.ndarray, covariance: np.ndarray) -> tuple[float, int, float]:
    rank = int(np.linalg.matrix_rank(covariance))
    if rank <= 0:
        return float("nan"), 0, float("nan")
    stat = float(vector @ np.linalg.pinv(covariance) @ vector)
    return stat, rank, _chi_square_sf_integer_df(stat, rank)


def _contrast_matrices(outcomes: list[str], config: dict[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    index = {name: i for i, name in enumerate(outcomes)}
    family_rows: list[np.ndarray] = []
    within_rows: list[np.ndarray] = []
    for family in ("accessibility_generalization", "reproductive_assurance"):
        members = [str(x) for x in config["response_families"][family] if str(x) in index]
        if len(members) < 2:
            continue
        row = np.zeros(len(outcomes), dtype=float)
        for member in members:
            row[index[member]] = 1.0 / len(members)
        family_rows.append(row)
        anchor = members[0]
        for member in members[1:]:
            contrast = np.zeros(len(outcomes), dtype=float)
            contrast[index[anchor]] = 1.0
            contrast[index[member]] = -1.0
            within_rows.append(contrast)
    return np.vstack(family_rows), np.vstack(within_rows)


def fit_compression_audit(
    counts: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
    *,
    stratum: str = "all_observed",
    context_a: str = "northern_midlatitude",
    context_b: str = "tropical",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    data = _prepare(counts, covariates, config)
    context = str(config["context_column"])
    geography = str(config["geography_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    threshold = int(config["minimum_islands_per_outcome"])
    work = data.loc[
        data["stratum"].eq(stratum) & data[context].isin([context_a, context_b])
    ].copy()
    support = work.groupby(["outcome", context])["island_id"].nunique().unstack(fill_value=0)
    for value in (context_a, context_b):
        if value not in support.columns:
            support[value] = 0
    outcomes = [
        str(x) for x in config["model_outcomes"]
        if str(x) in support.index
        and int(support.loc[str(x), context_a]) >= threshold
        and int(support.loc[str(x), context_b]) >= threshold
    ]
    if len(outcomes) != len(config["model_outcomes"]):
        raise typer.BadParameter(
            f"compression audit requires all declared outcomes; retained={outcomes}"
        )

    fits: list[dict[str, Any]] = []
    cluster_parts: list[np.ndarray] = []
    interaction_indices: list[int] = []
    component_rows: list[dict[str, Any]] = []
    offset = 0
    for outcome in outcomes:
        part = work.loc[work["outcome"].eq(outcome)].copy()
        b = part[context].eq(context_b).to_numpy(float)
        columns = [np.ones(len(part), dtype=float), b]
        names = [f"{outcome}:intercept", f"{outcome}:context[{context_b}]"]
        for predictor in baseline:
            z = _standardize(part[predictor])
            columns.extend([z, z * b])
            names.extend([
                f"{outcome}:z_{predictor}",
                f"{outcome}:z_{predictor}:context[{context_b}]",
            ])
        z_geo = _standardize(part[geography])
        interaction = f"{outcome}:z_{geography}:context[{context_b}]"
        columns.extend([z_geo, z_geo * b])
        names.extend([f"{outcome}:z_{geography}", interaction])
        fit = _fit_single_beta_binomial(
            part["successes"].to_numpy(float),
            part["trials"].to_numpy(float),
            np.column_stack(columns),
            names,
            max_iter=int(config.get("max_iter", 1000)),
        )
        fits.append(fit)
        cluster_parts.append(part[cluster].to_numpy(str))
        interaction_indices.append(offset + names.index(interaction))
        component_rows.append({
            "outcome": outcome,
            "n_islands_context_a": int(support.loc[outcome, context_a]),
            "n_islands_context_b": int(support.loc[outcome, context_b]),
            "optimizer_success": bool(fit["success"]),
        })
        offset += len(fit["names"])

    covariance, _, theta = _assemble_cluster_covariance(fits, cluster_parts)
    vector = theta[interaction_indices]
    vcov = covariance[np.ix_(interaction_indices, interaction_indices)]
    se = np.sqrt(np.clip(np.diag(vcov), 0.0, None))
    for row, estimate, stderr in zip(component_rows, vector, se, strict=True):
        row["slope_difference_tropical_minus_north"] = float(estimate)
        row["cluster_robust_se"] = float(stderr)

    family, within = _contrast_matrices(outcomes, config)
    full_stat, full_df, full_p = _wald(vector, vcov)
    family_vector = family @ vector
    family_cov = family @ vcov @ family.T
    fam_stat, fam_df, fam_p = _wald(family_vector, family_cov)
    within_vector = within @ vector
    within_cov = within @ vcov @ within.T
    within_stat, within_df, within_p = _wald(within_vector, within_cov)

    projection = family.T @ np.linalg.pinv(family @ family.T) @ family
    total_ss = float(vector @ vector)
    family_ss = float((projection @ vector) @ (projection @ vector))
    family_fraction = family_ss / total_ss if total_ss > 0 else float("nan")

    summary = pd.DataFrame([
        {
            "contrast": "full_six_axis",
            "df": full_df,
            "wald_chisq": full_stat,
            "p_value": full_p,
            "euclidean_fraction_retained_by_two_family_subspace": family_fraction,
        },
        {
            "contrast": "two_family_means",
            "df": fam_df,
            "wald_chisq": fam_stat,
            "p_value": fam_p,
            "euclidean_fraction_retained_by_two_family_subspace": family_fraction,
        },
        {
            "contrast": "within_family_structure_discarded_by_compression",
            "df": within_df,
            "wald_chisq": within_stat,
            "p_value": within_p,
            "euclidean_fraction_retained_by_two_family_subspace": family_fraction,
        },
    ])
    return pd.DataFrame(component_rows), summary


@app.command("run")
def run(
    status_flora_csv: Path = typer.Option(..., exists=True),
    state_audit_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    counts = build_broad_counts(
        pd.read_csv(status_flora_csv), pd.read_csv(state_audit_csv), config
    )
    components, summary = fit_compression_audit(counts, pd.read_csv(covariates_csv), config)
    components.insert(0, "evidence_scope", evidence_scope)
    summary.insert(0, "evidence_scope", evidence_scope)
    output_dir.mkdir(parents=True, exist_ok=True)
    components.to_csv(output_dir / "compression_component_estimates.csv", index=False)
    summary.to_csv(output_dir / "compression_audit_summary.csv", index=False)
    typer.echo(summary.to_csv(index=False))


if __name__ == "__main__":
    app()
