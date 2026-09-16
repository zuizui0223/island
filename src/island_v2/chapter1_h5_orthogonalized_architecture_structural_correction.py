"""Post-hoc structural correction for the Chapter 1 H5 architecture decomposition.

The frozen v1 analysis exposed two structural mismatches: removing one shared factor
from three template dimensions leaves a residual subspace of intrinsic dimension two,
and the reused multivariate helper intentionally cannot test a one-response shared
factor.  This module corrects only those two mathematical/interface defects.  It does
not change the source-trained factor, evidence hierarchy, all-observed flora scope,
North/Tropical contexts, per-island >=50-species gate, controls, or claim ceiling.

This is explicitly post hoc and cannot promote realized pollinator identity or causal
pollination-service mechanisms.
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

from island_v2.chapter1_h5_orthogonalized_architecture_specificity import (
    IDENTITY_AXIS_SET,
    run_scope,
)
from island_v2.chapter1_pr138_syndrome_analysis import _fit_weighted_clustered_design

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT = "chapter1_h5_orthogonalized_architecture_structural_correction_v2"


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected orthogonalized structural-correction contract")
    return config


def _normal_two_sided_p(z_value: float) -> float:
    return math.erfc(abs(float(z_value)) / math.sqrt(2.0)) if math.isfinite(z_value) else float("nan")


def _z(values: pd.Series, *, name: str) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise typer.BadParameter(f"constant or invalid scalar H2 predictor: {name}")
    return (x - mean) / sd


def summarize_rank2_residual_h2(
    between: pd.DataFrame,
    *,
    evidence_scope: str,
    config: dict[str, Any],
) -> dict[str, Any]:
    correction = config["residual_family_correction"]
    contexts = [str(x) for x in config["unchanged_from_v1"]["contexts"]]
    required_labels = sorted(str(x) for x in correction["labels"])
    part = between.loc[
        between["context_layer"].astype(str).eq("analysis_regime")
        & between["axis_set"].astype(str).eq(IDENTITY_AXIS_SET)
        & between["stratum"].astype(str).eq("all_observed")
        & between["support_tier"].astype(str).eq("confirmatory")
        & between["context_a"].astype(str).eq(contexts[0])
        & between["context_b"].astype(str).eq(contexts[1])
    ]
    if len(part) != 1:
        raise ValueError(f"expected one residual H2 row, found {len(part)}")
    row = part.iloc[0]
    observed_labels = sorted(
        x for x in str(row.get("retained_syndromes", "")).split("|") if x
    )
    observed_df = (
        int(row["joint_context_difference_df"])
        if pd.notna(row.get("joint_context_difference_df"))
        else 0
    )
    p_value = float(row["p_value"]) if pd.notna(row.get("p_value")) else None
    evaluable = bool(
        str(row.get("status", "")) == "fit"
        and observed_labels == required_labels
        and int(row.get("n_retained_syndromes", 0) or 0)
        == int(correction["required_labels_present"])
        and observed_df == int(correction["required_joint_wald_df"])
        and p_value is not None
        and math.isfinite(p_value)
    )
    if not evaluable:
        p_value = None
    return {
        "evidence_scope": evidence_scope,
        "evaluable": evaluable,
        "joint_chisq": float(row["joint_context_difference_chisq"])
        if evaluable and pd.notna(row.get("joint_context_difference_chisq"))
        else None,
        "joint_df": observed_df if evaluable else None,
        "joint_p": p_value,
        "retained_labels": observed_labels,
        "n_unique_islands": int(row["n_unique_islands"])
        if pd.notna(row.get("n_unique_islands"))
        else None,
        "n_clusters": int(row["n_clusters"])
        if pd.notna(row.get("n_clusters"))
        else None,
    }


def _unique_covariates(covariates: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    required = {"island_id", *columns}
    if missing := required - set(covariates.columns):
        raise typer.BadParameter(f"scalar H2 covariates missing columns: {sorted(missing)}")
    work = covariates[["island_id", *columns]].copy()
    multiplicity = work.groupby("island_id", dropna=False)[columns].nunique(dropna=False)
    conflicts = multiplicity.gt(1).any(axis=1)
    if bool(conflicts.any()):
        examples = [str(x) for x in conflicts.index[conflicts][:5]]
        raise typer.BadParameter(f"conflicting scalar H2 covariate rows for island_id: {examples}")
    return work.drop_duplicates("island_id")


def fit_scalar_shared_h2(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> dict[str, Any]:
    shared = config["shared_factor_correction"]
    unchanged = config["unchanged_from_v1"]
    contexts = [str(x) for x in unchanged["contexts"]]
    if contexts != ["northern_midlatitude", "tropical"]:
        raise typer.BadParameter("scalar shared H2 requires frozen northern_midlatitude/tropical contexts")
    response_name = str(shared["response"])
    context_column = "analysis_regime"
    cluster_column = str(unchanged["cluster_column"])
    geography = str(unchanged["geography_predictor"])
    controls = [str(x) for x in unchanged["controls"]]
    required_score = {"island_id", "stratum", "syndrome", "syndrome_score", "n_species"}
    if missing := required_score - set(island_scores.columns):
        raise typer.BadParameter(f"island component table missing columns: {sorted(missing)}")

    part = island_scores.loc[
        island_scores["stratum"].astype(str).eq("all_observed")
        & island_scores["syndrome"].astype(str).eq(response_name)
    ].copy()
    threshold_species = int(unchanged["minimum_scored_species_per_island_component"])
    part["n_species"] = pd.to_numeric(part["n_species"], errors="coerce")
    part = part.loc[part["n_species"].ge(threshold_species)].copy()

    cov_columns = [context_column, cluster_column, geography, *controls]
    cov = _unique_covariates(covariates, cov_columns)
    work = part.merge(cov, on="island_id", how="left", validate="one_to_one")
    numeric = ["syndrome_score", geography, *controls]
    for column in numeric:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.loc[work[context_column].astype(str).isin(contexts)].copy()
    work = work.dropna(subset=numeric)
    work = work.loc[work[cluster_column].fillna("").astype(str).ne("")].copy()

    counts = work.groupby(context_column)["island_id"].nunique().to_dict()
    minimum = int(shared["minimum_islands_per_context"])
    support_passed = all(int(counts.get(context, 0)) >= minimum for context in contexts)
    base = {
        "support_passed": support_passed,
        "northern_midlatitude_n_islands": int(counts.get(contexts[0], 0)),
        "tropical_n_islands": int(counts.get(contexts[1], 0)),
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(work[cluster_column].nunique()),
    }
    if not support_passed:
        return {
            **base,
            "evaluable": False,
            "distance_by_tropical_estimate": None,
            "distance_by_tropical_se": None,
            "distance_by_tropical_p": None,
            "northern_distance_slope": None,
            "northern_distance_slope_p": None,
            "tropical_distance_slope": None,
            "tropical_distance_slope_p": None,
        }

    names = ["intercept"]
    columns: list[np.ndarray] = [np.ones(len(work), dtype=float)]
    for control in controls:
        names.append(f"z_{control}")
        columns.append(_z(work[control], name=control))
    tropical = work[context_column].astype(str).eq(contexts[1]).to_numpy(float)
    names.append("context_tropical")
    columns.append(tropical)
    z_geo = _z(work[geography], name=geography)
    names.append(f"z_{geography}")
    columns.append(z_geo)
    interaction_name = f"z_{geography}:context_tropical"
    names.append(interaction_name)
    columns.append(z_geo * tropical)

    coefficients, covariance, fit = _fit_weighted_clustered_design(
        work["syndrome_score"].to_numpy(float),
        np.ones(len(work), dtype=float),
        np.column_stack(columns),
        names,
        work[cluster_column].astype(str).to_numpy(),
    )
    if coefficients.empty or str(fit.get("status")) != "fit":
        return {
            **base,
            "evaluable": False,
            "fit_status": str(fit.get("status", "fit_failed")),
            "distance_by_tropical_estimate": None,
            "distance_by_tropical_se": None,
            "distance_by_tropical_p": None,
            "northern_distance_slope": None,
            "northern_distance_slope_p": None,
            "tropical_distance_slope": None,
            "tropical_distance_slope_p": None,
        }

    coef = coefficients.set_index("predictor")
    geo_name = f"z_{geography}"
    geo_index = names.index(geo_name)
    interaction_index = names.index(interaction_name)
    north = float(coef.loc[geo_name, "estimate"])
    interaction = float(coef.loc[interaction_name, "estimate"])
    north_se = float(coef.loc[geo_name, "cluster_robust_se"])
    north_z = north / north_se if north_se > 0 else float("nan")
    tropical_slope = north + interaction
    tropical_var = float(
        covariance[geo_index, geo_index]
        + covariance[interaction_index, interaction_index]
        + 2.0 * covariance[geo_index, interaction_index]
    )
    tropical_se = math.sqrt(max(tropical_var, 0.0))
    tropical_z = tropical_slope / tropical_se if tropical_se > 0 else float("nan")
    return {
        **base,
        "evaluable": True,
        "fit_status": "fit",
        "distance_by_tropical_estimate": interaction,
        "distance_by_tropical_se": float(coef.loc[interaction_name, "cluster_robust_se"]),
        "distance_by_tropical_p": float(coef.loc[interaction_name, "p_value"]),
        "northern_distance_slope": north,
        "northern_distance_slope_se": north_se,
        "northern_distance_slope_p": _normal_two_sided_p(north_z),
        "tropical_distance_slope": tropical_slope,
        "tropical_distance_slope_se": tropical_se,
        "tropical_distance_slope_p": _normal_two_sided_p(tropical_z),
    }


def integrate_structural_correction(
    all_result: dict[str, Any],
    direct_result: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    residual_alpha = float(config["residual_family_correction"]["alpha"])
    shared_alpha = float(config["shared_factor_correction"]["alpha"])

    def supported(scope: dict[str, Any], key: str, p_key: str, alpha: float) -> bool:
        part = scope[key]
        p_value = part.get(p_key)
        return bool(part.get("evaluable") and p_value is not None and float(p_value) <= alpha)

    all_residual = supported(all_result, "residual", "joint_p", residual_alpha)
    direct_residual = supported(direct_result, "residual", "joint_p", residual_alpha)
    all_shared = supported(all_result, "shared", "distance_by_tropical_p", shared_alpha)
    direct_shared = supported(direct_result, "shared", "distance_by_tropical_p", shared_alpha)
    robust_residual = all_residual and direct_residual
    robust_shared = all_shared and direct_shared

    if robust_residual and robust_shared:
        classification = "shared_and_identity_specific_architecture_H2_supported_posthoc"
    elif robust_residual:
        classification = "identity_specific_residual_H2_without_shared_factor_posthoc"
    elif robust_shared:
        classification = "shared_architecture_H2_without_residual_specificity_posthoc"
    elif all_residual or direct_residual:
        classification = (
            "all_analysis_only_residual_signal_scope_sensitive"
            if all_residual
            else "direct_only_residual_signal_scope_sensitive"
        )
    elif all_shared or direct_shared:
        classification = (
            "all_analysis_only_shared_signal_scope_sensitive"
            if all_shared
            else "direct_only_shared_signal_scope_sensitive"
        )
    else:
        classification = "no_robust_named_architecture_H2_after_structural_correction"

    return {
        "contract": CONTRACT,
        "status": "completed_posthoc_structural_correction",
        "all_analysis": all_result,
        "direct_only": direct_result,
        "robust_identity_specific_residual_H2": robust_residual,
        "robust_shared_architecture_H2": robust_shared,
        "classification": classification,
        "confirmatory_specificity_test": False,
        "pollinator_identity_promoted": False,
        "causal_pollinator_mechanism_promoted": False,
    }


def run_analysis(
    all_species_scores: pd.DataFrame,
    direct_species_scores: pd.DataFrame,
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    branching_config: dict[str, Any],
    v1_config: dict[str, Any],
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    scope_rows = []
    score_parts = []
    between_parts = []
    for evidence_scope, scores in (
        ("all_analysis_eligible", all_species_scores),
        ("direct_only", direct_species_scores),
    ):
        _, island_scores, _, between, v1_decision, _ = run_scope(
            scores,
            status_flora,
            gift_flora,
            covariates,
            realm_assignment,
            pattern_config,
            branching_config,
            v1_config,
            evidence_scope=evidence_scope,
        )
        residual = summarize_rank2_residual_h2(
            between,
            evidence_scope=evidence_scope,
            config=config,
        )
        shared = fit_scalar_shared_h2(island_scores, covariates, config)
        scope_rows.append(
            {
                "evidence_scope": evidence_scope,
                "factor_variance_fraction": float(v1_decision["factor_variance_fraction"]),
                "residual": residual,
                "shared": shared,
            }
        )
        score_parts.append(island_scores)
        between_parts.append(between)

    integrated = integrate_structural_correction(scope_rows[0], scope_rows[1], config)
    return (
        pd.concat(score_parts, ignore_index=True),
        pd.concat(between_parts, ignore_index=True),
        integrated,
    )


@app.command("run")
def run(
    all_species_scores_csv: Path = typer.Option(..., exists=True),
    direct_species_scores_csv: Path = typer.Option(..., exists=True),
    status_flora_csv: Path = typer.Option(..., exists=True),
    gift_flora_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    realm_assignment_csv: Path = typer.Option(..., exists=True),
    pattern_config_path: Path = typer.Option(..., exists=True),
    branching_config_path: Path = typer.Option(..., exists=True),
    v1_config_path: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    v1_config = yaml.safe_load(v1_config_path.read_text(encoding="utf-8"))
    pattern_config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    branching_config = yaml.safe_load(branching_config_path.read_text(encoding="utf-8"))
    island_scores, between, decision = run_analysis(
        pd.read_csv(all_species_scores_csv),
        pd.read_csv(direct_species_scores_csv),
        pd.read_csv(status_flora_csv),
        pd.read_csv(gift_flora_csv),
        pd.read_csv(covariates_csv),
        pd.read_csv(realm_assignment_csv),
        pattern_config,
        branching_config,
        v1_config,
        config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    island_scores.to_csv(
        output_dir / "strict_support_architecture_components.csv.gz",
        index=False,
        compression="gzip",
    )
    between.to_csv(output_dir / "rank_audited_residual_between_context.csv", index=False)
    (output_dir / "structural_correction_decision.json").write_text(
        json.dumps(decision, indent=2) + "\n",
        encoding="utf-8",
    )
    all_result = decision["all_analysis"]
    direct_result = decision["direct_only"]
    lines = [
        "# H5 orthogonalized architecture — post-hoc structural correction",
        "",
        f"- all-analysis PCA1 shared variance: {all_result['factor_variance_fraction']:.6g}",
        f"- Direct PCA1 shared variance: {direct_result['factor_variance_fraction']:.6g}",
        f"- all-analysis rank-2 residual H2 p: {all_result['residual']['joint_p']}",
        f"- Direct rank-2 residual H2 p: {direct_result['residual']['joint_p']}",
        f"- all-analysis shared-factor H2 interaction p: {all_result['shared']['distance_by_tropical_p']}",
        f"- Direct shared-factor H2 interaction p: {direct_result['shared']['distance_by_tropical_p']}",
        "",
        f"**Classification:** {decision['classification']}",
        "",
        "This is a post-hoc structural correction of the frozen v1 rank/interface defects, not a confirmatory pollinator-identity test.",
    ]
    (output_dir / "RESULT_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.echo("\n".join(lines))


if __name__ == "__main__":
    app()
