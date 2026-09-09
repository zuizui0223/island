"""Prospective MNAR tipping-point sensitivity for Chapter 1 trait resolution.

The observed syndrome score on an island is the mean among species with enough
trait information to score that syndrome.  This module asks how the completed
island mean would change if focal-state species had different odds of being
trait-resolved than non-focal species.  It never codes unresolved species as
zeros and it retains the observed number of resolved species as the regression
information weight.

The odds-ratio grid and context contrasts are read from the prospectively frozen
``chapter1_explanation_gap_validation_v1`` contract.  Biological estimates are
allowed to disappear or reverse; those outcomes are recorded rather than used
to retune the grid.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Annotated, Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_global_branching import run_global_branching
from island_v2.flora_status_support import stratum_mask

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT = "chapter1_trait_resolution_mnar_tipping_v1"
VALIDATION_KEY = "V5_trait_resolution_MNAR_tipping_point"


def completed_prevalence_from_resolution_or(
    *,
    total_species: float,
    resolved_species: float,
    resolved_focal_count: float,
    resolution_odds_ratio: float,
) -> float:
    """Identify completed prevalence under a binary selection-odds model.

    ``resolution_odds_ratio`` is

    ``odds(R=1 | focal) / odds(R=1 | nonfocal)``.

    The formula solves the two observed-data equations exactly.  It supports
    fractional focal counts because island syndrome means are soft membership
    scores.  OR=1 exactly reproduces the observed resolved prevalence.
    """

    total = float(total_species)
    resolved = float(resolved_species)
    focal = float(resolved_focal_count)
    odds_ratio = float(resolution_odds_ratio)
    if not math.isfinite(total) or total <= 0:
        raise ValueError("total_species must be finite and positive")
    if not math.isfinite(resolved) or resolved <= 0 or resolved > total:
        raise ValueError("resolved_species must be in (0, total_species]")
    if not math.isfinite(focal) or focal < 0 or focal > resolved:
        raise ValueError("resolved_focal_count must be in [0, resolved_species]")
    if not math.isfinite(odds_ratio) or odds_ratio <= 0:
        raise ValueError("resolution_odds_ratio must be finite and positive")

    observed_focal = focal / total
    observed_nonfocal = (resolved - focal) / total
    if focal <= 0:
        return 0.0
    if resolved - focal <= 0:
        return 1.0

    nonfocal_resolution = (
        observed_focal / odds_ratio + observed_nonfocal
    ) / (1.0 - observed_focal + observed_focal / odds_ratio)
    focal_resolution = (
        odds_ratio * nonfocal_resolution
        / (1.0 + (odds_ratio - 1.0) * nonfocal_resolution)
    )
    prevalence = observed_focal / focal_resolution
    return float(np.clip(prevalence, 0.0, 1.0))


def build_total_species_counts(
    status_flora: pd.DataFrame,
    strata: list[str],
) -> pd.DataFrame:
    required = {
        "island_id",
        "accepted_species",
        "origin_status",
        "endemic_status",
        "floristic_status",
    }
    missing = required - set(status_flora.columns)
    if missing:
        raise typer.BadParameter(f"status flora missing columns: {sorted(missing)}")
    flora = status_flora.copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    rows: list[pd.DataFrame] = []
    for stratum in strata:
        subset = flora.loc[stratum_mask(flora, stratum)].copy()
        subset = subset.loc[subset["accepted_species"].ne("")]
        counts = (
            subset.drop_duplicates(["island_id", "accepted_species"])
            .groupby("island_id", as_index=False)["accepted_species"]
            .nunique()
            .rename(columns={"accepted_species": "n_total_stratum_species"})
        )
        counts["stratum"] = str(stratum)
        rows.append(counts)
    if not rows:
        return pd.DataFrame(
            columns=["island_id", "n_total_stratum_species", "stratum"]
        )
    return pd.concat(rows, ignore_index=True)


def _context_assignments(
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
) -> pd.DataFrame:
    required_cov = {"island_id", "analysis_regime"}
    required_realm = {"island_id", "biogeographic_realm"}
    missing_cov = required_cov - set(covariates.columns)
    missing_realm = required_realm - set(realm_assignment.columns)
    if missing_cov:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing_cov)}")
    if missing_realm:
        raise typer.BadParameter(
            f"realm assignment missing columns: {sorted(missing_realm)}"
        )
    cov = covariates[["island_id", "analysis_regime"]].copy()
    realm = realm_assignment[["island_id", "biogeographic_realm"]].copy()
    for frame in (cov, realm):
        frame["island_id"] = frame["island_id"].astype(str)
    cov = cov.drop_duplicates("island_id")
    realm = realm.drop_duplicates("island_id")
    return cov.merge(realm, on="island_id", how="left", validate="one_to_one")


def _scenario_rows(v5_config: dict[str, Any]) -> list[dict[str, Any]]:
    grid = [float(x) for x in v5_config["freeze_before_execution"]["resolution_odds_ratio_grid"]]
    if 1.0 not in grid:
        raise ValueError("MNAR odds-ratio grid must contain the OR=1 baseline")
    families = v5_config["freeze_before_execution"]["scenario_families"]
    rows: list[dict[str, Any]] = []
    for family, spec in families.items():
        contexts = [str(x) for x in spec.get("contrast_contexts", [])]
        if contexts and len(contexts) != 2:
            raise ValueError(f"scenario family {family} must declare zero or two contexts")
        for odds_ratio in grid:
            rows.append(
                {
                    "scenario_id": f"{family}__or_{odds_ratio:.10g}".replace(".", "p"),
                    "scenario_family": str(family),
                    "scenario_type": "selection_grid",
                    "context_layer": str(spec["context_layer"]),
                    "context_a": contexts[0] if contexts else "",
                    "context_b": contexts[1] if contexts else "",
                    "resolution_odds_ratio": odds_ratio,
                    "log_resolution_odds_ratio": math.log(odds_ratio),
                    "bound_mode": "",
                }
            )
        if not contexts:
            bound_modes = ["all_unresolved_nonfocal", "all_unresolved_focal"]
        else:
            bound_modes = [
                "first_context_nonfocal_second_context_focal",
                "first_context_focal_second_context_nonfocal",
            ]
        for bound_mode in bound_modes:
            rows.append(
                {
                    "scenario_id": f"{family}__{bound_mode}",
                    "scenario_family": str(family),
                    "scenario_type": "partial_identification_bound",
                    "context_layer": str(spec["context_layer"]),
                    "context_a": contexts[0] if contexts else "",
                    "context_b": contexts[1] if contexts else "",
                    "resolution_odds_ratio": float("nan"),
                    "log_resolution_odds_ratio": float("nan"),
                    "bound_mode": bound_mode,
                }
            )
    return rows


def _local_resolution_or(
    assignments: pd.DataFrame,
    scenario: dict[str, Any],
) -> np.ndarray:
    odds_ratio = float(scenario["resolution_odds_ratio"])
    layer = str(scenario["context_layer"])
    if layer == "all_islands":
        return np.full(len(assignments), odds_ratio, dtype=float)
    first = str(scenario["context_a"])
    second = str(scenario["context_b"])
    log_or = math.log(odds_ratio)
    local = np.ones(len(assignments), dtype=float)
    local[assignments[layer].astype(str).eq(first).to_numpy()] = math.exp(log_or / 2.0)
    local[assignments[layer].astype(str).eq(second).to_numpy()] = math.exp(-log_or / 2.0)
    return local


def _bound_prevalence(
    *,
    observed_focal: np.ndarray,
    unresolved: np.ndarray,
    total: np.ndarray,
    assignments: pd.DataFrame,
    scenario: dict[str, Any],
) -> np.ndarray:
    lower = observed_focal / total
    upper = (observed_focal + unresolved) / total
    mode = str(scenario["bound_mode"])
    if mode == "all_unresolved_nonfocal":
        return lower
    if mode == "all_unresolved_focal":
        return upper
    layer = str(scenario["context_layer"])
    first = assignments[layer].astype(str).eq(str(scenario["context_a"])).to_numpy()
    second = assignments[layer].astype(str).eq(str(scenario["context_b"])).to_numpy()
    resolved = total - unresolved
    prevalence = np.divide(
        observed_focal,
        resolved,
        out=np.zeros_like(observed_focal, dtype=float),
        where=resolved > 0,
    )
    if mode == "first_context_nonfocal_second_context_focal":
        prevalence[first] = lower[first]
        prevalence[second] = upper[second]
        return prevalence
    if mode == "first_context_focal_second_context_nonfocal":
        prevalence[first] = upper[first]
        prevalence[second] = lower[second]
        return prevalence
    raise ValueError(f"unknown partial-identification bound: {mode}")


def adjust_island_scores(
    island_scores: pd.DataFrame,
    total_species: pd.DataFrame,
    assignments: pd.DataFrame,
    scenario: dict[str, Any],
    affected_axes: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {"island_id", "stratum", "syndrome", "syndrome_score", "n_species"}
    missing = required - set(island_scores.columns)
    if missing:
        raise typer.BadParameter(f"island scores missing columns: {sorted(missing)}")
    scores = island_scores.copy()
    scores["island_id"] = scores["island_id"].astype(str)
    scores["stratum"] = scores["stratum"].astype(str)
    scores["syndrome"] = scores["syndrome"].astype(str)
    scores["syndrome_score"] = pd.to_numeric(scores["syndrome_score"], errors="coerce")
    scores["n_species"] = pd.to_numeric(scores["n_species"], errors="coerce")
    scores = scores.merge(
        total_species,
        on=["island_id", "stratum"],
        how="left",
        validate="many_to_one",
    ).merge(assignments, on="island_id", how="left", validate="many_to_one")
    affected = scores["syndrome"].isin(affected_axes)
    if scores.loc[affected, "n_total_stratum_species"].isna().any():
        raise ValueError("affected island scores lack total stratum species counts")
    if (scores.loc[affected, "n_species"] > scores.loc[affected, "n_total_stratum_species"]).any():
        raise ValueError("resolved syndrome species exceed total stratum species")

    diagnostics = scores.loc[affected].copy()
    if diagnostics.empty:
        raise ValueError(f"no island scores found for affected axes: {sorted(affected_axes)}")
    resolved = diagnostics["n_species"].to_numpy(float)
    total = diagnostics["n_total_stratum_species"].to_numpy(float)
    observed_membership = np.clip(
        (diagnostics["syndrome_score"].to_numpy(float) + 1.0) / 2.0,
        0.0,
        1.0,
    )
    observed_focal = observed_membership * resolved
    unresolved = total - resolved

    if str(scenario["scenario_type"]) == "selection_grid":
        local_or = _local_resolution_or(diagnostics, scenario)
        completed = np.array(
            [
                completed_prevalence_from_resolution_or(
                    total_species=n_total,
                    resolved_species=n_resolved,
                    resolved_focal_count=n_focal,
                    resolution_odds_ratio=or_value,
                )
                for n_total, n_resolved, n_focal, or_value in zip(
                    total, resolved, observed_focal, local_or, strict=True
                )
            ],
            dtype=float,
        )
    else:
        local_or = np.full(len(diagnostics), np.nan)
        completed = _bound_prevalence(
            observed_focal=observed_focal,
            unresolved=unresolved,
            total=total,
            assignments=diagnostics,
            scenario=scenario,
        )

    adjusted_score = 2.0 * completed - 1.0
    scores.loc[affected, "syndrome_score"] = adjusted_score
    diagnostics["observed_soft_membership"] = observed_membership
    diagnostics["completed_soft_membership"] = completed
    diagnostics["observed_syndrome_score"] = 2.0 * observed_membership - 1.0
    diagnostics["completed_syndrome_score"] = adjusted_score
    diagnostics["n_unresolved_species"] = unresolved
    diagnostics["local_resolution_odds_ratio"] = local_or
    for key in (
        "scenario_id",
        "scenario_family",
        "scenario_type",
        "context_layer",
        "context_a",
        "context_b",
        "resolution_odds_ratio",
        "log_resolution_odds_ratio",
        "bound_mode",
    ):
        diagnostics[key] = scenario[key]
    return scores[list(island_scores.columns)], diagnostics


def _add_scenario_columns(
    frame: pd.DataFrame,
    scenario: dict[str, Any],
    evidence_scope: str,
) -> pd.DataFrame:
    out = frame.copy()
    metadata = {
        "evidence_scope": evidence_scope,
        **{
            key: scenario[key]
            for key in (
                "scenario_id",
                "scenario_family",
                "scenario_type",
                "resolution_odds_ratio",
                "log_resolution_odds_ratio",
                "bound_mode",
            )
        },
    }
    for key, value in reversed(list(metadata.items())):
        out.insert(0, key, value)
    return out


def _as_bool(value: object) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes"}


def _tipping_table(
    frame: pd.DataFrame,
    *,
    result_type: str,
    key_columns: list[str],
    support_column: str,
    q_column: str,
    estimate_column: str | None,
) -> pd.DataFrame:
    if frame.empty or "scenario_type" not in frame.columns:
        return pd.DataFrame()
    grid = frame.loc[frame["scenario_type"].eq("selection_grid")].copy()
    if grid.empty:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    group_columns = ["evidence_scope", "scenario_family", *key_columns]
    for keys, group in grid.groupby(group_columns, dropna=False, sort=True):
        key_values = keys if isinstance(keys, tuple) else (keys,)
        identity = dict(zip(group_columns, key_values, strict=True))
        baseline = group.loc[np.isclose(group["resolution_odds_ratio"], 1.0)]
        if len(baseline) != 1:
            raise ValueError(f"expected one OR=1 baseline for {identity}, found {len(baseline)}")
        base = baseline.iloc[0]
        baseline_supported = _as_bool(base.get(support_column, False))
        baseline_estimate = (
            float(base[estimate_column]) if estimate_column is not None else float("nan")
        )
        for direction, candidates in (
            ("below_one", group.loc[group["resolution_odds_ratio"].lt(1.0)]),
            ("above_one", group.loc[group["resolution_odds_ratio"].gt(1.0)]),
        ):
            candidates = candidates.assign(
                absolute_log_or=candidates["log_resolution_odds_ratio"].abs()
            ).sort_values("absolute_log_or")
            break_row: pd.Series | None = None
            break_events: list[str] = []
            for _, candidate in candidates.iterrows():
                events: list[str] = []
                candidate_supported = _as_bool(candidate.get(support_column, False))
                if baseline_supported and not candidate_supported:
                    events.append("support_lost")
                if estimate_column is not None:
                    estimate = float(candidate[estimate_column])
                    if (
                        math.isfinite(baseline_estimate)
                        and math.isfinite(estimate)
                        and baseline_estimate != 0
                        and estimate * baseline_estimate <= 0
                    ):
                        events.append("sign_changed")
                if events:
                    break_row = candidate
                    break_events = events
                    break
            if not baseline_supported:
                verdict = "baseline_not_supported"
            elif break_row is None:
                verdict = "no_break_in_prespecified_grid"
            else:
                verdict = "break_detected"
            rows.append(
                {
                    "result_type": result_type,
                    **identity,
                    "odds_ratio_direction": direction,
                    "baseline_supported": baseline_supported,
                    "baseline_q_value": float(base.get(q_column, float("nan"))),
                    "baseline_estimate": baseline_estimate,
                    "tipping_resolution_odds_ratio": (
                        float(break_row["resolution_odds_ratio"])
                        if break_row is not None
                        else float("nan")
                    ),
                    "tipping_absolute_log_or": (
                        float(abs(break_row["log_resolution_odds_ratio"]))
                        if break_row is not None
                        else float("nan")
                    ),
                    "tipping_event": "|".join(break_events),
                    "verdict": verdict,
                }
            )
    return pd.DataFrame(rows)


def build_tipping_summary(
    slopes: pd.DataFrame,
    within: pd.DataFrame,
    between: pd.DataFrame,
) -> pd.DataFrame:
    parts = [
        _tipping_table(
            slopes,
            result_type="within_context_slope",
            key_columns=[
                "context_layer",
                "axis_set",
                "stratum",
                "support_tier",
                "context",
                "syndrome",
            ],
            support_column="axis_supported",
            q_column="q_axis_family",
            estimate_column="distance_slope",
        ),
        _tipping_table(
            within,
            result_type="within_context_vector",
            key_columns=[
                "context_layer",
                "axis_set",
                "stratum",
                "support_tier",
                "context",
            ],
            support_column="context_vector_supported",
            q_column="q_context_vector_family",
            estimate_column=None,
        ),
        _tipping_table(
            between,
            result_type="between_context_vector",
            key_columns=[
                "context_layer",
                "axis_set",
                "stratum",
                "support_tier",
                "context_a",
                "context_b",
            ],
            support_column="context_vector_difference_supported",
            q_column="q_between_context_family",
            estimate_column=None,
        ),
    ]
    usable = [part for part in parts if not part.empty]
    return pd.concat(usable, ignore_index=True) if usable else pd.DataFrame()


def run_mnar_sensitivity(
    *,
    island_scores: pd.DataFrame,
    status_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    branching_config: dict[str, Any],
    explanation_config: dict[str, Any],
    evidence_scope: str,
) -> dict[str, pd.DataFrame | dict[str, Any]]:
    if explanation_config.get("contract") != "chapter1_explanation_gap_validation_v1":
        raise ValueError("unexpected explanation-gap contract")
    v5_config = explanation_config["validations"][VALIDATION_KEY]
    fixed = v5_config["freeze_before_execution"]
    affected_axes = {str(x) for x in fixed["affected_syndrome_axes"]}
    strata = [str(x) for x in pattern_config["strata"]]
    totals = build_total_species_counts(status_flora, strata)
    assignments = _context_assignments(covariates, realm_assignment)
    scenarios = _scenario_rows(v5_config)

    diagnostics_parts: list[pd.DataFrame] = []
    slope_parts: list[pd.DataFrame] = []
    within_parts: list[pd.DataFrame] = []
    between_parts: list[pd.DataFrame] = []
    state_parts: list[pd.DataFrame] = []
    baseline_score_differences: list[float] = []

    baseline_scores = island_scores[
        ["island_id", "stratum", "syndrome", "syndrome_score"]
    ].copy()
    baseline_scores["island_id"] = baseline_scores["island_id"].astype(str)
    baseline_scores["stratum"] = baseline_scores["stratum"].astype(str)
    baseline_scores["syndrome"] = baseline_scores["syndrome"].astype(str)
    baseline_scores = baseline_scores.set_index(["island_id", "stratum", "syndrome"])

    for scenario in scenarios:
        adjusted, diagnostics = adjust_island_scores(
            island_scores,
            totals,
            assignments,
            scenario,
            affected_axes,
        )
        diagnostics.insert(0, "evidence_scope", evidence_scope)
        diagnostics_parts.append(diagnostics)
        if (
            scenario["scenario_type"] == "selection_grid"
            and math.isclose(float(scenario["resolution_odds_ratio"]), 1.0)
        ):
            check = adjusted[
                ["island_id", "stratum", "syndrome", "syndrome_score"]
            ].copy()
            for column in ("island_id", "stratum", "syndrome"):
                check[column] = check[column].astype(str)
            check = check.set_index(["island_id", "stratum", "syndrome"])
            difference = pd.to_numeric(
                check["syndrome_score"] - baseline_scores["syndrome_score"],
                errors="coerce",
            ).abs()
            baseline_score_differences.append(float(difference.max()))

        _, slopes, within, between, states = run_global_branching(
            adjusted,
            covariates,
            realm_assignment,
            pattern_config,
            branching_config,
        )
        slope_parts.append(_add_scenario_columns(slopes, scenario, evidence_scope))
        within_parts.append(_add_scenario_columns(within, scenario, evidence_scope))
        between_parts.append(_add_scenario_columns(between, scenario, evidence_scope))
        state_parts.append(_add_scenario_columns(states, scenario, evidence_scope))

    slopes_all = pd.concat(slope_parts, ignore_index=True)
    within_all = pd.concat(within_parts, ignore_index=True)
    between_all = pd.concat(between_parts, ignore_index=True)
    states_all = pd.concat(state_parts, ignore_index=True)
    tipping = build_tipping_summary(slopes_all, within_all, between_all)
    maximum_baseline_difference = max(baseline_score_differences, default=float("nan"))
    if not math.isfinite(maximum_baseline_difference) or maximum_baseline_difference > 1e-12:
        raise ValueError(
            f"OR=1 failed to reproduce baseline island scores: {maximum_baseline_difference}"
        )
    manifest = {
        "contract": CONTRACT,
        "status": "post_baseline_prospective_falsification",
        "evidence_scope": evidence_scope,
        "affected_syndrome_axes": sorted(affected_axes),
        "resolution_odds_ratio_grid": [
            float(x) for x in fixed["resolution_odds_ratio_grid"]
        ],
        "scenario_families": list(fixed["scenario_families"]),
        "n_scenarios": len(scenarios),
        "n_tipping_rows": len(tipping),
        "maximum_OR1_score_difference": maximum_baseline_difference,
        "unresolved_species_are_biological_zeros": False,
        "regression_information_is_imputed": False,
        "grid_selected_from_outcomes": False,
        "favorable_results_promote_H5": False,
    }
    return {
        "scenario_table": pd.DataFrame(scenarios),
        "score_diagnostics": pd.concat(diagnostics_parts, ignore_index=True),
        "slopes": slopes_all,
        "within": within_all,
        "between": between_all,
        "states": states_all,
        "tipping": tipping,
        "manifest": manifest,
    }


@app.command("run")
def run(
    island_scores_csv: Annotated[Path, typer.Option(exists=True)],
    status_flora_csv: Annotated[Path, typer.Option(exists=True)],
    covariates_csv: Annotated[Path, typer.Option(exists=True)],
    realm_assignment_csv: Annotated[Path, typer.Option(exists=True)],
    pattern_config_path: Annotated[Path, typer.Option(exists=True)],
    branching_config_path: Annotated[Path, typer.Option(exists=True)],
    explanation_config_path: Annotated[Path, typer.Option(exists=True)],
    output_dir: Annotated[Path, typer.Option()],
    evidence_scope: Annotated[str, typer.Option()] = "all_analysis_eligible",
) -> None:
    result = run_mnar_sensitivity(
        island_scores=pd.read_csv(island_scores_csv),
        status_flora=pd.read_csv(status_flora_csv),
        covariates=pd.read_csv(covariates_csv),
        realm_assignment=pd.read_csv(realm_assignment_csv),
        pattern_config=yaml.safe_load(pattern_config_path.read_text(encoding="utf-8")),
        branching_config=yaml.safe_load(
            branching_config_path.read_text(encoding="utf-8")
        ),
        explanation_config=yaml.safe_load(
            explanation_config_path.read_text(encoding="utf-8")
        ),
        evidence_scope=evidence_scope,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    frame_outputs = {
        "scenario_table": "mnar_scenarios.csv",
        "score_diagnostics": "mnar_score_diagnostics.csv.gz",
        "slopes": "mnar_branch_slopes.csv.gz",
        "within": "mnar_branch_within.csv.gz",
        "between": "mnar_branch_between.csv.gz",
        "states": "mnar_branch_states.csv.gz",
        "tipping": "mnar_tipping_points.csv",
    }
    for key, filename in frame_outputs.items():
        frame = result[key]
        assert isinstance(frame, pd.DataFrame)
        compression = "gzip" if filename.endswith(".gz") else None
        frame.to_csv(output_dir / filename, index=False, compression=compression)
    manifest = result["manifest"]
    assert isinstance(manifest, dict)
    (output_dir / "chapter1_trait_resolution_mnar_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
