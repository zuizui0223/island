"""Common-support adjudication of the H5 named-architecture null.

The orthogonalized named-architecture analysis uses a narrow complete-case support:
species must score all three sampled named templates and islands must contain at least
50 such species.  Before interpreting a null architecture H2 result, this module asks
whether the primary six-atomic North--Tropical H2 itself remains on exactly that
species/island support.

This is post-hoc adjudication.  It does not replace the full-data H2 result and it does
not observe pollinator identity or effective pollination service.
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

from island_v2.chapter1_pr138_biogeographic_pattern import (
    build_observed_broad_counts,
    run_observed_pattern,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT = "chapter1_h5_architecture_common_support_h2_v1"
SHARED_COMPONENT = "shared_architecture_factor"


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected H5 architecture common-support contract")
    return config


def complete_template_species(
    species_scores: pd.DataFrame,
    config: dict[str, Any],
) -> set[str]:
    required = {"accepted_species", "syndrome", "syndrome_concordance"}
    if missing := required - set(species_scores.columns):
        raise typer.BadParameter(f"species syndrome scores missing columns: {sorted(missing)}")
    templates = [str(x) for x in config["named_template_complete_case_basis"]]
    work = species_scores.loc[
        species_scores["syndrome"].astype(str).isin(templates),
        ["accepted_species", "syndrome", "syndrome_concordance"],
    ].copy()
    work["accepted_species"] = work["accepted_species"].fillna("").astype(str)
    work["syndrome_concordance"] = pd.to_numeric(
        work["syndrome_concordance"], errors="coerce"
    )
    work = work.loc[
        work["accepted_species"].ne("") & work["syndrome_concordance"].notna()
    ].drop_duplicates(["accepted_species", "syndrome"])
    counts = work.groupby("accepted_species")["syndrome"].nunique()
    return set(counts.index[counts.eq(len(templates))].astype(str))


def _strict_support_islands(
    strict_support: pd.DataFrame,
    config: dict[str, Any],
    *,
    evidence_scope: str,
) -> set[str]:
    required = {"evidence_scope", "island_id", "syndrome", "n_species", "stratum"}
    if missing := required - set(strict_support.columns):
        raise typer.BadParameter(f"strict architecture support missing columns: {sorted(missing)}")
    minimum = int(config["minimum_complete_template_species_per_island"])
    work = strict_support.copy()
    work["island_id"] = work["island_id"].astype(str)
    work["n_species"] = pd.to_numeric(work["n_species"], errors="coerce")
    work = work.loc[
        work["evidence_scope"].astype(str).eq(evidence_scope)
        & work["stratum"].astype(str).eq(str(config["flora_scope"]))
        & work["syndrome"].astype(str).eq(SHARED_COMPONENT)
        & work["n_species"].ge(minimum)
    ]
    return set(work["island_id"].astype(str))


def restrict_common_support(
    status_flora: pd.DataFrame,
    state_audit: pd.DataFrame,
    covariates: pd.DataFrame,
    strict_support: pd.DataFrame,
    species_scores: pd.DataFrame,
    *,
    config: dict[str, Any],
    evidence_scope: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    required_flora = {"island_id", "accepted_species"}
    required_audit = {"accepted_species", "trait_name", "resolved_for_primary", "canonical_signature"}
    if missing := required_flora - set(status_flora.columns):
        raise typer.BadParameter(f"status flora missing columns: {sorted(missing)}")
    if missing := required_audit - set(state_audit.columns):
        raise typer.BadParameter(f"state audit missing columns: {sorted(missing)}")
    if "island_id" not in covariates.columns:
        raise typer.BadParameter("covariates missing island_id")

    species = complete_template_species(species_scores, config)
    islands = _strict_support_islands(
        strict_support,
        config,
        evidence_scope=evidence_scope,
    )
    contexts = {str(x) for x in config["contexts"]}
    context_column = str(config["model"]["context_column"])
    if context_column not in covariates.columns:
        raise typer.BadParameter(f"covariates missing context column: {context_column}")

    cov = covariates.copy()
    cov["island_id"] = cov["island_id"].astype(str)
    cov = cov.loc[
        cov["island_id"].isin(islands)
        & cov[context_column].fillna("").astype(str).isin(contexts)
    ].copy()
    allowed_islands = set(cov["island_id"].astype(str))

    flora = status_flora.copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    flora = flora.loc[
        flora["island_id"].isin(allowed_islands)
        & flora["accepted_species"].isin(species)
    ].copy()

    audit = state_audit.copy()
    audit["accepted_species"] = audit["accepted_species"].astype(str)
    audit = audit.loc[audit["accepted_species"].isin(species)].copy()
    return flora, audit, cov


def focused_atomic_config(
    pattern_config: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    focused = deepcopy(pattern_config)
    outcomes = [str(x) for x in config["atomic_outcomes"]]
    missing = [x for x in outcomes if x not in pattern_config["broad_outcomes"]]
    if missing:
        raise typer.BadParameter(f"atomic outcomes absent from pattern contract: {missing}")
    focused["contexts"] = [str(x) for x in config["contexts"]]
    focused["strata"] = [str(config["flora_scope"])]
    focused["support_tiers"] = {
        "confirmatory": int(config["minimum_islands_per_context_per_outcome"])
    }
    focused["minimum_outcomes_per_vector"] = int(config["required_atomic_outcomes"])
    focused["broad_outcomes"] = {
        outcome: deepcopy(pattern_config["broad_outcomes"][outcome])
        for outcome in outcomes
    }
    return focused


def _summarize_between(
    between: pd.DataFrame,
    *,
    evidence_scope: str,
    config: dict[str, Any],
) -> dict[str, Any]:
    contexts = [str(x) for x in config["contexts"]]
    part = between.loc[
        between["stratum"].astype(str).eq(str(config["flora_scope"]))
        & between["support_tier"].astype(str).eq("confirmatory")
        & between["context_a"].astype(str).eq(contexts[0])
        & between["context_b"].astype(str).eq(contexts[1])
    ]
    if len(part) != 1:
        raise ValueError(f"expected one common-support North-Tropical H2 row, found {len(part)}")
    row = part.iloc[0]
    retained = [x for x in str(row.get("retained_outcomes", "")).split("|") if x]
    expected = [str(x) for x in config["atomic_outcomes"]]
    p_value = float(row["p_value"]) if pd.notna(row.get("p_value")) else None
    evaluable = bool(
        str(row.get("status", "")) == "fit"
        and int(row.get("n_retained_outcomes", 0) or 0) == int(config["required_atomic_outcomes"])
        and set(retained) == set(expected)
        and p_value is not None
        and math.isfinite(p_value)
    )
    if not evaluable:
        p_value = None
    return {
        "evidence_scope": evidence_scope,
        "evaluable": evaluable,
        "n_retained_outcomes": int(row.get("n_retained_outcomes", 0) or 0),
        "retained_outcomes": retained,
        "n_unique_islands": int(row["n_unique_islands"])
        if pd.notna(row.get("n_unique_islands"))
        else None,
        "n_clusters": int(row["n_clusters"])
        if pd.notna(row.get("n_clusters"))
        else None,
        "joint_chisq": float(row["joint_wald_chisq"])
        if evaluable and pd.notna(row.get("joint_wald_chisq"))
        else None,
        "joint_df": int(row["joint_df"])
        if evaluable and pd.notna(row.get("joint_df"))
        else None,
        "joint_p": p_value,
    }


def integrate_decisions(
    all_result: dict[str, Any],
    direct_result: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    alpha = float(config["alpha"])
    required = int(config["required_atomic_outcomes"])

    def supported(result: dict[str, Any]) -> bool:
        p_value = result.get("joint_p")
        return bool(
            result.get("evaluable")
            and int(result.get("n_retained_outcomes", 0)) == required
            and p_value is not None
            and float(p_value) <= alpha
        )

    all_ok = supported(all_result)
    direct_ok = supported(direct_result)
    if all_ok and direct_ok:
        classification = "common_support_H2_retained_architecture_null_informative"
        informative = True
    elif all_ok or direct_ok:
        classification = "common_support_H2_not_robust_across_evidence_scopes"
        informative = False
    elif not all_result.get("evaluable") or not direct_result.get("evaluable"):
        classification = "common_support_H2_not_evaluable_in_both_scopes"
        informative = False
    else:
        classification = "common_support_H2_not_supported_in_either_scope"
        informative = False
    return {
        "contract": CONTRACT,
        "status": "completed_posthoc_common_support_adjudication",
        "all_analysis": all_result,
        "direct_only": direct_result,
        "architecture_null_informative_about_H2_representation": informative,
        "classification": classification,
        "full_data_H2_replaced": False,
        "pollinator_mechanism_identified": False,
    }


def run_scope(
    species_scores: pd.DataFrame,
    state_audit: pd.DataFrame,
    status_flora: pd.DataFrame,
    strict_support: pd.DataFrame,
    covariates: pd.DataFrame,
    pattern_config: dict[str, Any],
    config: dict[str, Any],
    *,
    evidence_scope: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any], dict[str, Any]]:
    flora, audit, cov = restrict_common_support(
        status_flora,
        state_audit,
        covariates,
        strict_support,
        species_scores,
        config=config,
        evidence_scope=evidence_scope,
    )
    focused = focused_atomic_config(pattern_config, config)
    counts = build_observed_broad_counts(flora, audit, focused)
    within_slopes, between_slopes, within, between = run_observed_pattern(
        counts,
        cov,
        focused,
    )
    decision = _summarize_between(
        between,
        evidence_scope=evidence_scope,
        config=config,
    )
    context_column = str(config["model"]["context_column"])
    support_audit = {
        "evidence_scope": evidence_scope,
        "n_complete_template_species": len(complete_template_species(species_scores, config)),
        "n_support_islands": int(cov["island_id"].astype(str).nunique()),
        "support_islands_by_context": {
            str(key): int(value)
            for key, value in cov.groupby(context_column)["island_id"].nunique().items()
        },
        "n_restricted_flora_rows": len(flora),
        "n_restricted_species": int(flora["accepted_species"].nunique()),
    }
    counts = counts.copy()
    counts.insert(0, "evidence_scope", evidence_scope)
    between_slopes = between_slopes.copy()
    between_slopes.insert(0, "evidence_scope", evidence_scope)
    within = within.copy()
    within.insert(0, "evidence_scope", evidence_scope)
    between = between.copy()
    between.insert(0, "evidence_scope", evidence_scope)
    return counts, between_slopes, within, between, decision, support_audit


@app.command("run")
def run(
    all_species_scores_csv: Path = typer.Option(..., exists=True),
    direct_species_scores_csv: Path = typer.Option(..., exists=True),
    all_state_audit_csv: Path = typer.Option(..., exists=True),
    direct_state_audit_csv: Path = typer.Option(..., exists=True),
    status_flora_csv: Path = typer.Option(..., exists=True),
    strict_support_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    pattern_config_path: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    pattern_config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    status_flora = pd.read_csv(status_flora_csv)
    strict_support = pd.read_csv(strict_support_csv)
    covariates = pd.read_csv(covariates_csv)

    parts = []
    decisions = []
    audits = []
    for evidence_scope, species_path, audit_path in (
        ("all_analysis_eligible", all_species_scores_csv, all_state_audit_csv),
        ("direct_only", direct_species_scores_csv, direct_state_audit_csv),
    ):
        result = run_scope(
            pd.read_csv(species_path),
            pd.read_csv(audit_path),
            status_flora,
            strict_support,
            covariates,
            pattern_config,
            config,
            evidence_scope=evidence_scope,
        )
        parts.append(result[:4])
        decisions.append(result[4])
        audits.append(result[5])

    integrated = integrate_decisions(decisions[0], decisions[1], config)
    integrated["support_audit"] = audits
    output_dir.mkdir(parents=True, exist_ok=True)
    pd.concat([x[0] for x in parts], ignore_index=True).to_csv(
        output_dir / "common_support_atomic_counts.csv.gz",
        index=False,
        compression="gzip",
    )
    pd.concat([x[1] for x in parts], ignore_index=True).to_csv(
        output_dir / "common_support_between_outcome_slope_differences.csv",
        index=False,
    )
    pd.concat([x[2] for x in parts], ignore_index=True).to_csv(
        output_dir / "common_support_within_region_omnibus.csv",
        index=False,
    )
    pd.concat([x[3] for x in parts], ignore_index=True).to_csv(
        output_dir / "common_support_between_region_omnibus.csv",
        index=False,
    )
    (output_dir / "common_support_decision.json").write_text(
        json.dumps(integrated, indent=2) + "\n",
        encoding="utf-8",
    )
    lines = [
        "# H5 architecture common-support six-atomic H2 gate",
        "",
        f"- all-analysis H2 p: {decisions[0]['joint_p']}",
        f"- Direct H2 p: {decisions[1]['joint_p']}",
        f"- all-analysis retained atomic outcomes: {decisions[0]['n_retained_outcomes']}",
        f"- Direct retained atomic outcomes: {decisions[1]['n_retained_outcomes']}",
        f"- all-analysis H2 islands/blocks: {decisions[0]['n_unique_islands']} / {decisions[0]['n_clusters']}",
        f"- Direct H2 islands/blocks: {decisions[1]['n_unique_islands']} / {decisions[1]['n_clusters']}",
        "",
        f"**Classification:** {integrated['classification']}",
        "",
        "This post-hoc gate adjudicates whether the named-architecture null is informative on the broad H2 support; it does not replace full-data H2 or identify pollinator mechanism.",
    ]
    (output_dir / "RESULT_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.echo("\n".join(lines))


if __name__ == "__main__":
    app()
