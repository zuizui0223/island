#!/usr/bin/env python3
"""Validate the v2 layered, conflict-safe, joint-composition analysis contract."""

from __future__ import annotations

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "config" / "v2_workflow_registry.yml"
CONTRACT = ROOT / "config" / "v2_analysis_contract.yml"
LAYERED_SCRIPT = ROOT / "analysis" / "v2" / "run_inla_category_preserving_north.R"
JOINT_SCRIPT = ROOT / "analysis" / "v2" / "run_inla_joint_color_form_composition.R"
STATE_BUILDER = ROOT / "src" / "island_v2" / "build_species_trait_states.py"
AGGREGATOR = ROOT / "src" / "island_v2" / "aggregate_island_trait_evidence.py"


def fail(message: str) -> None:
    raise SystemExit(f"v2 workflow registry validation failed: {message}")


def read_text(path: Path) -> str:
    if not path.is_file():
        fail(f"required file does not exist: {path.relative_to(ROOT)}")
    return path.read_text(encoding="utf-8")


def require_tokens(text: str, tokens: dict[str, str], context: str) -> None:
    missing = [name for name, token in tokens.items() if token not in text]
    if missing:
        fail(f"{context} is missing: {', '.join(sorted(missing))}")


def main() -> None:
    registry = yaml.safe_load(REGISTRY.read_text(encoding="utf-8"))
    contract = yaml.safe_load(CONTRACT.read_text(encoding="utf-8"))
    current = registry["canonical"]["current"]
    workflow_path = ROOT / current["workflow"]
    workflow_text = read_text(workflow_path)
    layered_text = read_text(LAYERED_SCRIPT)
    joint_text = read_text(JOINT_SCRIPT)
    state_text = read_text(STATE_BUILDER)
    aggregate_text = read_text(AGGREGATOR)
    tier = str(current.get("evidence_tier", ""))

    if current.get("engine") != "INLA":
        fail("canonical workflow must use INLA")
    if f"bombus_join/{tier}/" not in workflow_text:
        fail(f"workflow does not read evidence tier {tier!r}")
    if f"--analysis-tier {tier}" not in workflow_text:
        fail(f"workflow does not pass evidence tier {tier!r}")
    require_tokens(
        workflow_text,
        {
            "layered analysis": "run_inla_category_preserving_north.R",
            "joint analysis": "run_inla_joint_color_form_composition.R",
            "positive-distance support": "distance_to_continent_km'] > 0",
            "joint trial audit": "full[joint_cols].fillna(0).sum(axis=1)",
            "conflict audit": "n_species_with_any_trait_conflict",
            "four regimes": "expected_regimes.issubset(regimes)",
        },
        "canonical workflow",
    )

    canonical = contract.get("canonical_analysis", {})
    if canonical.get("workflow") != current["workflow"]:
        fail("analysis contract and registry disagree on workflow")
    if canonical.get("engine") != "INLA":
        fail("analysis contract must use INLA")
    if canonical.get("current_evidence_tier") != tier:
        fail("analysis contract and registry disagree on evidence tier")
    if canonical.get("planned_confirmatory_evidence_tier") != "primary":
        fail("planned confirmatory evidence tier must remain primary")
    if canonical.get("layered_all_data") is not True:
        fail("canonical analysis must remain layered_all_data")
    if canonical.get("joint_colour_form_composition") is not True:
        fail("canonical analysis must include joint colour-form composition")

    expected_colours = {"plain", "yellow_orange", "red_pink", "blue_purple"}
    expected_forms = {
        "open_generalized",
        "tubular_trumpet",
        "zygomorphic_specialized",
        "composite_brush",
    }
    retained = contract.get("retained_categories", {})
    if set(retained.get("flower_color", [])) != expected_colours:
        fail("retained flower-colour categories are incomplete")
    if set(retained.get("floral_form", [])) != expected_forms:
        fail("retained floral-form categories are incomplete")
    if retained.get("joint_colour_form_cells") != 16:
        fail("joint colour-form composition must contain 16 cells")

    principles = contract.get("principles", {})
    required_true = {
        "conflicting_species_states_are_unresolved",
        "separate_binary_models_are_not_joint_composition",
        "joint_colour_form_counts_must_exhaust_joint_trials",
        "exclude_zero_distance_from_fitted_models",
        "use_equation_specific_maximum_support",
        "compare_models_only_on_matched_support",
        "preserve_bombus_channel_components",
        "do_not_force_collinear_bombus_components_into_one_equation",
        "alternative_guilds_are_global_descriptive_outcomes",
    }
    for name in required_true:
        if principles.get(name) is not True:
            fail(f"analysis contract must require {name}")
    for name in (
        "category_omission_allowed",
        "outcome_layer_omission_allowed",
        "alternative_pollinator_covariates_in_primary_models",
    ):
        if principles.get(name) is not False:
            fail(f"analysis contract must prohibit {name}")

    require_tokens(
        state_text,
        {
            "conflict-safe resolution": "count(DISTINCT filled_value)",
            "conflict flag": "conflict_flag",
            "mixed pollen-vector state": "pollen_vector_mode = 'mixed'",
        },
        "species-state builder",
    )
    if "max(CASE WHEN trait_name='flower_primary_color' THEN filled_value END)" in state_text:
        fail("species-state builder still selects raw conflicting values with max()")

    require_tokens(
        aggregate_text,
        {
            "joint trial count": "joint_color_form_trials",
            "joint cell counts": "joint_{color}__{form}",
            "conflict totals": "n_species_with_any_trait_conflict",
            "mixed species": "n_mixed_species",
        },
        "island trait aggregator",
    )

    require_tokens(
        layered_text,
        {
            "all four colours": "color_blue_purple",
            "all four forms": "form_composite_brush",
            "wind-adjusted path": "N4_wind_adjusted",
            "four regimes": '"northern_high_latitude"',
        },
        "layered INLA analysis",
    )
    require_tokens(
        joint_text,
        {
            "16-cell input": "joint_color_form_trials",
            "shared island effect": "island_comp_id",
            "Poisson log-linear family": 'family = "poisson"',
            "joint Bombus interaction": "joint_category:z_bombus_deficit",
            "joint SC interaction": "joint_category:z_sc_share",
            "joint wind interaction": "joint_category:z_wind_mixed_share",
            "joint total check": "joint colour-form cells do not sum to joint_color_form_trials",
        },
        "joint colour-form INLA analysis",
    )

    joint_models = {
        "J0_joint_geo",
        "J1_joint_bombus",
        "J2_joint_sc_wind",
        "J3_joint_full",
    }
    if not joint_models.issubset(set(contract.get("models", {}))):
        fail("analysis contract is missing joint composition model ladder")

    expected_regimes = {
        "northern_midlatitude",
        "northern_high_latitude",
        "tropical",
        "southern_extratropical",
    }
    if set(canonical.get("population", {}).get("global_regimes", [])) != expected_regimes:
        fail("global regime set is incomplete")

    print(
        {
            "registry": str(REGISTRY.relative_to(ROOT)),
            "contract": str(CONTRACT.relative_to(ROOT)),
            "workflow": str(workflow_path.relative_to(ROOT)),
            "engine": "INLA",
            "evidence_tier": tier,
            "joint_colour_form_cells": 16,
            "conflict_safe_species_states": True,
            "status": "ok",
        }
    )


if __name__ == "__main__":
    main()
