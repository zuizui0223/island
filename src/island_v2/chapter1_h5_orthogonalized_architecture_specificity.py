"""Focused H5 test of shared versus template-specific floral architecture.

This module reuses the prospectively defined V4 source-trained PCA decomposition but
applies it to the current all-observed Chapter 1 route.  It tests the predeclared
northern-midlatitude versus tropical distance-response contrast separately for the
shared architecture factor and for the three template-specific residuals.

Residual labels are plant-trait contrasts.  They do not identify realized pollinators,
visitation, effective service, historical loss, or causal selection.
"""
from __future__ import annotations

import json
import math
from copy import deepcopy
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.chapter1_global_branching import run_global_branching
from island_v2.chapter1_pollination_architecture_factor import fit_and_project_source_factor
from island_v2.chapter1_pr138_syndrome_analysis import build_island_syndrome_scores

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT = "chapter1_h5_orthogonalized_architecture_specificity_v1"
IDENTITY_AXIS_SET = "identity_specific_residuals"
SHARED_AXIS_SET = "shared_architecture"


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected orthogonalized architecture specificity contract")
    return config


def filter_island_component_support(
    island_scores: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Apply the frozen species-per-island support gate before H2 modelling."""

    if "n_species" not in island_scores.columns:
        raise typer.BadParameter("island component table missing n_species")
    threshold = int(config["island_aggregation"]["minimum_scored_species_per_component"])
    out = island_scores.copy()
    out["n_species"] = pd.to_numeric(out["n_species"], errors="coerce")
    return out.loc[out["n_species"].ge(threshold)].copy()


def focused_configs(
    pattern_config: dict[str, Any],
    branching_config: dict[str, Any],
    config: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Restrict the existing branching machinery to the frozen H5 specificity test."""

    pattern = deepcopy(pattern_config)
    branching = deepcopy(branching_config)

    contexts = [str(x) for x in config["model"]["contexts"]]
    stratum = str(config["flora_scope"]["primary"])
    threshold = int(config["island_aggregation"]["minimum_scored_species_per_component"])
    residuals = [str(x) for x in config["projected_components"][IDENTITY_AXIS_SET]]
    shared = [str(x) for x in config["projected_components"]["shared"]]

    pattern["contexts"] = contexts
    pattern["strata"] = [stratum]
    pattern["support_tiers"] = {"confirmatory": threshold}
    # The shared architecture family is intentionally one-dimensional.  The primary
    # residual family is separately required to retain all three axes in summarize_scope.
    pattern["minimum_outcomes_per_vector"] = 1

    branching["contract"] = CONTRACT
    branching["hypothesis_name"] = "H5 orthogonalized floral-architecture specificity"
    components = [*shared, *residuals]
    branching["branch_axes"] = {
        component: {"components": {component: 1.0}} for component in components
    }
    branching["axis_sets"] = {
        SHARED_AXIS_SET: {
            "axes": shared,
            "role": "decomposition_context_not_identity_specificity_gate",
            "classify": False,
        },
        IDENTITY_AXIS_SET: {
            "axes": residuals,
            "role": "secondary_H5_identity_specificity_falsification",
            "classify": False,
        },
    }
    context_column = str(config["model"]["context_column"])
    branching["context_layers"] = {
        "analysis_regime": {
            "column": context_column,
            "contexts": contexts,
            "role": "frozen_H2_north_tropical_layer",
        }
    }
    branching["alpha"] = float(config["primary_tests"]["identity_specific_H2"]["alpha"])
    return pattern, branching


def _select_between_row(
    between: pd.DataFrame,
    *,
    axis_set: str,
    config: dict[str, Any],
) -> pd.Series:
    contexts = [str(x) for x in config["model"]["contexts"]]
    if len(contexts) != 2:
        raise ValueError("orthogonalized specificity contract requires exactly two contexts")
    stratum = str(config["flora_scope"]["primary"])
    part = between.loc[
        between["context_layer"].astype(str).eq("analysis_regime")
        & between["axis_set"].astype(str).eq(axis_set)
        & between["stratum"].astype(str).eq(stratum)
        & between["support_tier"].astype(str).eq("confirmatory")
        & between["context_a"].astype(str).eq(contexts[0])
        & between["context_b"].astype(str).eq(contexts[1])
    ]
    if len(part) != 1:
        raise ValueError(f"expected one frozen between-context row for {axis_set}, found {len(part)}")
    return part.iloc[0]


def _context_slopes(
    slopes: pd.DataFrame,
    *,
    axis_set: str,
    config: dict[str, Any],
) -> dict[str, dict[str, float | None]]:
    if slopes.empty:
        return {}
    contexts = [str(x) for x in config["model"]["contexts"]]
    stratum = str(config["flora_scope"]["primary"])
    part = slopes.loc[
        slopes["context_layer"].astype(str).eq("analysis_regime")
        & slopes["axis_set"].astype(str).eq(axis_set)
        & slopes["stratum"].astype(str).eq(stratum)
        & slopes["support_tier"].astype(str).eq("confirmatory")
        & slopes["context"].astype(str).isin(contexts)
    ].copy()
    out: dict[str, dict[str, float | None]] = {}
    for row in part.itertuples(index=False):
        p_value = float(row.p_value) if pd.notna(row.p_value) else None
        out.setdefault(str(row.context), {})[str(row.syndrome)] = float(row.distance_slope)
        out[str(row.context)][f"{row.syndrome}__p"] = p_value
    return out


def summarize_scope(
    between: pd.DataFrame,
    slopes: pd.DataFrame,
    *,
    evidence_scope: str,
    factor_audit: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    identity = _select_between_row(between, axis_set=IDENTITY_AXIS_SET, config=config)
    shared = _select_between_row(between, axis_set=SHARED_AXIS_SET, config=config)

    required_identity = int(
        config["primary_tests"]["identity_specific_H2"]["required_retained_components"]
    )
    required_df = int(config["primary_tests"]["identity_specific_H2"]["joint_df_target"])
    identity_status = str(identity.get("status", ""))
    identity_n = int(identity.get("n_retained_syndromes", 0) or 0)
    identity_df = (
        int(identity["joint_context_difference_df"])
        if pd.notna(identity.get("joint_context_difference_df"))
        else 0
    )
    identity_p = (
        float(identity["p_value"])
        if pd.notna(identity.get("p_value"))
        else None
    )
    identity_evaluable = bool(
        identity_status == "fit"
        and identity_n == required_identity
        and identity_df == required_df
        and identity_p is not None
        and math.isfinite(identity_p)
    )
    if not identity_evaluable:
        identity_p = None

    shared_status = str(shared.get("status", ""))
    shared_n = int(shared.get("n_retained_syndromes", 0) or 0)
    shared_df = (
        int(shared["joint_context_difference_df"])
        if pd.notna(shared.get("joint_context_difference_df"))
        else 0
    )
    shared_p = float(shared["p_value"]) if pd.notna(shared.get("p_value")) else None
    shared_evaluable = bool(
        shared_status == "fit"
        and shared_n == 1
        and shared_df == 1
        and shared_p is not None
        and math.isfinite(shared_p)
    )
    if not shared_evaluable:
        shared_p = None

    return {
        "evidence_scope": evidence_scope,
        "identity_specific_evaluable": identity_evaluable,
        "identity_specific_joint_p": identity_p,
        "identity_specific_joint_df": identity_df if identity_evaluable else None,
        "identity_specific_n_retained_components": identity_n,
        "identity_specific_retained_components": str(identity.get("retained_syndromes", "")),
        "identity_specific_n_unique_islands": int(identity["n_unique_islands"])
        if pd.notna(identity.get("n_unique_islands"))
        else None,
        "identity_specific_n_clusters": int(identity["n_clusters"])
        if pd.notna(identity.get("n_clusters"))
        else None,
        "shared_architecture_evaluable": shared_evaluable,
        "shared_architecture_joint_p": shared_p,
        "shared_architecture_joint_df": shared_df if shared_evaluable else None,
        "shared_architecture_n_unique_islands": int(shared["n_unique_islands"])
        if pd.notna(shared.get("n_unique_islands"))
        else None,
        "shared_architecture_n_clusters": int(shared["n_clusters"])
        if pd.notna(shared.get("n_clusters"))
        else None,
        "factor_variance_fraction": float(factor_audit["variance_fraction"]),
        "n_complete_scored_species": int(factor_audit["n_complete_scored_species"]),
        "n_complete_gift_source_species": int(factor_audit["n_complete_gift_source_species"]),
        "source_factor_residual_covariances": factor_audit[
            "source_factor_residual_covariances"
        ],
        "identity_specific_context_slopes": _context_slopes(
            slopes, axis_set=IDENTITY_AXIS_SET, config=config
        ),
        "shared_architecture_context_slopes": _context_slopes(
            slopes, axis_set=SHARED_AXIS_SET, config=config
        ),
    }


def integrate_scope_decisions(
    primary: dict[str, Any],
    direct: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    alpha = float(config["primary_tests"]["identity_specific_H2"]["alpha"])
    primary_ok = bool(primary.get("identity_specific_evaluable"))
    direct_ok = bool(direct.get("identity_specific_evaluable"))
    p_primary = primary.get("identity_specific_joint_p")
    p_direct = direct.get("identity_specific_joint_p")
    primary_supported = bool(primary_ok and p_primary is not None and float(p_primary) <= alpha)
    direct_supported = bool(direct_ok and p_direct is not None and float(p_direct) <= alpha)

    if not primary_ok or not direct_ok:
        classification = "identity_specific_architecture_not_evaluable_in_both_scopes"
        robust = False
    elif primary_supported and direct_supported:
        classification = "robust_identity_specific_architecture_structure"
        robust = True
    elif primary_supported:
        classification = "all_analysis_only_identity_specific_signal_not_promoted"
        robust = False
    elif direct_supported:
        classification = "direct_only_identity_specific_signal_not_promoted"
        robust = False
    else:
        classification = "identity_specific_residual_H2_not_supported"
        robust = False

    return {
        "contract": CONTRACT,
        "status": "completed",
        "primary_scope": "all_analysis_eligible",
        "sensitivity_scope": "direct_only",
        "all_analysis": primary,
        "direct_only": direct,
        "robust_identity_specific_structure": robust,
        "classification": classification,
        "pollinator_identity_promoted": False,
        "causal_pollinator_mechanism_promoted": False,
        "maximum_claim": (
            "The all-observed North-Tropical H2 contrast contains plant-trait geometry beyond the dominant shared floral-architecture factor."
            if robust
            else "The current test does not promote robust identity-specific plant-architecture geometry across both evidence scopes."
        ),
    }


def run_scope(
    species_scores: pd.DataFrame,
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    branching_config: dict[str, Any],
    config: dict[str, Any],
    *,
    evidence_scope: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any], pd.DataFrame]:
    templates = tuple(str(x) for x in config["sampled_templates"])
    projected, factor_model, factor_audit = fit_and_project_source_factor(
        species_scores,
        gift_flora,
        templates=templates,
        minimum_complete_source_species=int(
            config["factor_training"]["minimum_complete_source_species"]
        ),
    )
    island_scores = build_island_syndrome_scores(
        status_flora,
        projected,
        [str(config["flora_scope"]["primary"])],
    )
    island_scores = filter_island_component_support(island_scores, config)
    focused_pattern, focused_branching = focused_configs(
        pattern_config, branching_config, config
    )
    _, slopes, within, between, _ = run_global_branching(
        island_scores,
        covariates,
        realm_assignment,
        focused_pattern,
        focused_branching,
    )
    decision = summarize_scope(
        between,
        slopes,
        evidence_scope=evidence_scope,
        factor_audit=factor_audit,
        config=config,
    )
    factor_model = factor_model.copy()
    factor_model.insert(0, "evidence_scope", evidence_scope)
    island_scores = island_scores.copy()
    island_scores.insert(0, "evidence_scope", evidence_scope)
    slopes = slopes.copy()
    slopes.insert(0, "evidence_scope", evidence_scope)
    within = within.copy()
    within.insert(0, "evidence_scope", evidence_scope)
    between = between.copy()
    between.insert(0, "evidence_scope", evidence_scope)
    return factor_model, island_scores, slopes, between, decision, within


def run_analysis(
    all_species_scores: pd.DataFrame,
    direct_species_scores: pd.DataFrame,
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    branching_config: dict[str, Any],
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    outputs = []
    decisions = []
    for scope, scores in (
        ("all_analysis_eligible", all_species_scores),
        ("direct_only", direct_species_scores),
    ):
        outputs.append(
            run_scope(
                scores,
                status_flora,
                gift_flora,
                covariates,
                realm_assignment,
                pattern_config,
                branching_config,
                config,
                evidence_scope=scope,
            )
        )
        decisions.append(outputs[-1][4])

    factor_models = pd.concat([x[0] for x in outputs], ignore_index=True)
    island_scores = pd.concat([x[1] for x in outputs], ignore_index=True)
    slopes = pd.concat([x[2] for x in outputs], ignore_index=True)
    between = pd.concat([x[3] for x in outputs], ignore_index=True)
    within = pd.concat([x[5] for x in outputs], ignore_index=True)
    integrated = integrate_scope_decisions(decisions[0], decisions[1], config)
    return factor_models, island_scores, slopes, within, between, integrated


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
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    pattern_config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    branching_config = yaml.safe_load(branching_config_path.read_text(encoding="utf-8"))
    factor_models, island_scores, slopes, within, between, decision = run_analysis(
        pd.read_csv(all_species_scores_csv),
        pd.read_csv(direct_species_scores_csv),
        pd.read_csv(status_flora_csv),
        pd.read_csv(gift_flora_csv),
        pd.read_csv(covariates_csv),
        pd.read_csv(realm_assignment_csv),
        pattern_config,
        branching_config,
        config,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    factor_models.to_csv(output_dir / "factor_models.csv", index=False)
    island_scores.to_csv(
        output_dir / "island_architecture_components.csv.gz",
        index=False,
        compression="gzip",
    )
    slopes.to_csv(output_dir / "architecture_context_slopes.csv", index=False)
    within.to_csv(output_dir / "architecture_within_context_omnibus.csv", index=False)
    between.to_csv(output_dir / "architecture_between_context_omnibus.csv", index=False)
    (output_dir / "orthogonalized_architecture_decision.json").write_text(
        json.dumps(decision, indent=2) + "\n", encoding="utf-8"
    )

    all_result = decision["all_analysis"]
    direct_result = decision["direct_only"]
    lines = [
        "# H5 orthogonalized architecture specificity",
        "",
        f"- all-analysis source PCA1 variance fraction: {all_result['factor_variance_fraction']:.6g}",
        f"- Direct source PCA1 variance fraction: {direct_result['factor_variance_fraction']:.6g}",
        f"- all-analysis identity-specific residual H2 p: {all_result['identity_specific_joint_p']}",
        f"- Direct identity-specific residual H2 p: {direct_result['identity_specific_joint_p']}",
        f"- all-analysis shared-architecture H2 p: {all_result['shared_architecture_joint_p']}",
        f"- Direct shared-architecture H2 p: {direct_result['shared_architecture_joint_p']}",
        "",
        f"**Classification:** {decision['classification']}",
        "",
        "Residual labels are plant-trait contrasts, not realized pollinator identities or effective-service measurements.",
    ]
    (output_dir / "RESULT_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.echo("\n".join(lines))


if __name__ == "__main__":
    app()
