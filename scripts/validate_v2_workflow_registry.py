#!/usr/bin/env python3
"""Validate the v2 layered all-data analysis contract and workflow."""

from __future__ import annotations

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "config" / "v2_workflow_registry.yml"
CONTRACT = ROOT / "config" / "v2_analysis_contract.yml"
INLA_SCRIPT = ROOT / "analysis" / "v2" / "run_inla_category_preserving_north.R"


def fail(message: str) -> None:
    raise SystemExit(f"v2 workflow registry validation failed: {message}")


def read_text(relative: str) -> str:
    path = ROOT / relative
    if not path.is_file():
        fail(f"registered workflow does not exist: {relative}")
    return path.read_text(encoding="utf-8")


def require_tokens(text: str, tokens: dict[str, str], context: str) -> None:
    missing = [name for name, token in tokens.items() if token not in text]
    if missing:
        fail(f"{context} is missing: {', '.join(sorted(missing))}")


def main() -> None:
    registry = yaml.safe_load(REGISTRY.read_text(encoding="utf-8"))
    contract = yaml.safe_load(CONTRACT.read_text(encoding="utf-8"))
    current = registry["canonical"]["current"]
    current_path = current["workflow"]
    workflow_text = read_text(current_path)
    inla_text = INLA_SCRIPT.read_text(encoding="utf-8")
    tier = str(current.get("evidence_tier", ""))

    if f"bombus_join/{tier}/" not in workflow_text:
        fail(f"current workflow does not read registry tier {tier!r}")
    if f"--analysis-tier {tier}" not in workflow_text:
        fail(f"current workflow does not pass analysis tier {tier!r}")

    if registry["guardrails"].get("exclude_zero_distance_from_fitted_models"):
        positive_distance_markers = (
            "distance_to_continent_km'].gt(0)",
            "distance_to_continent_km'] > 0",
            'distance_to_continent_km"] > 0',
        )
        if not any(marker in workflow_text for marker in positive_distance_markers):
            fail("workflow does not construct positive-distance model support")

    if registry["guardrails"].get("category_engine") != "INLA":
        fail("category_engine must be INLA")
    if current.get("engine") != "INLA":
        fail("canonical workflow must use INLA")
    if "run_inla_category_preserving_north.R" not in workflow_text:
        fail("canonical workflow does not run the INLA analysis")

    brms_role = registry.get("engine_roles", {}).get("bayesian_brms", {})
    if registry["guardrails"].get("prohibit_brms_category_models"):
        required_exclusions = {
            "category_preserving_flower_colour",
            "category_preserving_floral_form",
            "multinomial_category_models",
        }
        if not required_exclusions.issubset(set(brms_role.get("excludes", []))):
            fail("brms role must exclude category-preserving models")

    canonical = contract.get("canonical_analysis", {})
    if canonical.get("workflow") != current_path:
        fail("analysis contract and registry disagree on canonical workflow")
    if canonical.get("engine") != current.get("engine"):
        fail("analysis contract and registry disagree on engine")
    if canonical.get("current_evidence_tier") != tier:
        fail("analysis contract and registry disagree on evidence tier")
    if canonical.get("planned_confirmatory_evidence_tier") != "primary":
        fail("planned confirmatory evidence tier must remain primary")
    if canonical.get("layered_all_data") is not True:
        fail("canonical analysis must be declared layered_all_data")

    principles = contract.get("principles", {})
    required_true = {
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
    required_false = {
        "category_omission_allowed",
        "outcome_layer_omission_allowed",
        "alternative_pollinator_covariates_in_primary_models",
    }
    for name in required_false:
        if principles.get(name) is not False:
            fail(f"analysis contract must prohibit {name}")

    expected_colours = {"plain", "yellow_orange", "red_pink", "blue_purple"}
    expected_forms = {
        "open_generalized",
        "tubular_trumpet",
        "zygomorphic_specialized",
        "composite_brush",
    }
    retained = contract.get("retained_categories", {})
    contract_colours = set(retained.get("flower_color", []))
    contract_forms = set(retained.get("floral_form", []))
    if contract_colours != expected_colours:
        fail(f"retained flower colours are incomplete: {sorted(contract_colours)}")
    if contract_forms != expected_forms:
        fail(f"retained floral forms are incomplete: {sorted(contract_forms)}")

    require_tokens(
        inla_text,
        {
            "plain": 'plain = c("color_plain", "color_trials")',
            "yellow_orange": 'yellow_orange = c("color_yellow_orange", "color_trials")',
            "red_pink": 'red_pink = c("color_red_pink", "color_trials")',
            "blue_purple": 'blue_purple = c("color_blue_purple", "color_trials")',
            "open_generalized": 'open_generalized = c("form_open_generalized", "form_trials")',
            "tubular_trumpet": 'tubular_trumpet = c("form_tubular_trumpet", "form_trials")',
            "zygomorphic_specialized": 'zygomorphic_specialized = c("form_zygomorphic_specialized", "form_trials")',
            "composite_brush": 'composite_brush = c("form_composite_brush", "form_trials")',
        },
        "canonical INLA category implementation",
    )

    require_tokens(
        inla_text,
        {
            "colour total check": "colour categories do not sum to color_trials",
            "form total check": "form categories do not sum to form_trials",
            "animal/wind total check": "animal + wind species do not sum to observed pollen-vector trials",
            "wind trial identity": "wind_mixed_trials and animal_status_observed disagree",
            "wind success identity": "wind_mixed_successes and n_wind_species disagree",
        },
        "canonical INLA input checks",
    )
    require_tokens(
        workflow_text,
        {
            "colour workflow audit": "full[colours].fillna(0).sum(axis=1)",
            "form workflow audit": "full[forms].fillna(0).sum(axis=1)",
            "animal/wind workflow audit": "full['n_animal_species'] + full['n_wind_species']",
            "wind trial workflow audit": "full['wind_mixed_trials'].fillna(0) == full['animal_status_observed'].fillna(0)",
            "four-regime workflow audit": "expected_regimes.issubset(regimes)",
        },
        "workflow pre-fit audit",
    )

    expected_regimes = {
        "northern_midlatitude",
        "northern_high_latitude",
        "tropical",
        "southern_extratropical",
    }
    contract_regimes = set(canonical.get("population", {}).get("global_regimes", []))
    if contract_regimes != expected_regimes:
        fail(f"global regimes are incomplete: {sorted(contract_regimes)}")
    for regime in expected_regimes:
        if f'"{regime}"' not in inla_text:
            fail(f"INLA script omits regime {regime}")

    expected_full_flora = {"animal_pollinated", "wind_mixed", "self_compatibility"}
    if set(contract.get("full_flora_outcomes", [])) != expected_full_flora:
        fail("full-flora outcomes are incomplete")
    require_tokens(
        inla_text,
        {
            "animal-pollinated outcome": 'animal_pollinated = c("animal_pollinated_successes", "animal_pollinated_trials")',
            "wind/mixed outcome": 'wind_mixed = c("wind_mixed_successes", "wind_mixed_trials")',
            "SC outcome": 'self_compatibility = c("sc_successes", "sc_trials")',
            "wind-adjusted path": "N4_wind_adjusted",
        },
        "layered full-flora and path analysis",
    )

    expected_all_categories = expected_colours | expected_forms
    n3_categories = set(contract["models"]["N3_direct_indirect"].get("applies_to", []))
    if n3_categories != expected_all_categories:
        fail("N3 paths must apply to every retained flower category")
    if "N4_wind_adjusted" not in contract.get("models", {}):
        fail("contract must define wind-adjusted path model")

    expected_global = (
        expected_all_categories
        | expected_full_flora
        | {"showy_alternative_guild", "other_bee_guild", "generalist_insect_guild"}
    )
    global_responses = set(contract["models"]["global_falsification"].get("responses", []))
    if not expected_global.issubset(global_responses):
        fail("global falsification does not include all available outcome layers")

    forbidden_predictors = {
        "alternative_pollinator_guild",
        "showy_alternative_guild",
        "other_bee_guild",
        "generalist_insect_guild",
    }
    global_predictors = set(contract["models"]["global_falsification"].get("predictors", []))
    if global_predictors & forbidden_predictors:
        fail("alternative guilds must not be causal predictors")

    registered: set[str] = {current_path}
    for group in registry.get("robustness", {}).values():
        registered.update(group.get("workflows", []))
    for group in registry.get("replication", {}).values():
        workflow = group.get("workflow")
        if workflow:
            registered.add(workflow)
    registered.update(registry.get("legacy", {}).get("workflows", []))
    missing_files = [path for path in sorted(registered) if not (ROOT / path).is_file()]
    if missing_files:
        fail("registered workflow files are missing: " + ", ".join(missing_files))

    print(
        {
            "registry": str(REGISTRY.relative_to(ROOT)),
            "contract": str(CONTRACT.relative_to(ROOT)),
            "current_workflow": current_path,
            "current_engine": current.get("engine"),
            "current_evidence_tier": tier,
            "retained_colour_categories": sorted(contract_colours),
            "retained_form_categories": sorted(contract_forms),
            "full_flora_outcomes": sorted(expected_full_flora),
            "global_regimes": sorted(contract_regimes),
            "status": "ok",
        }
    )


if __name__ == "__main__":
    main()
