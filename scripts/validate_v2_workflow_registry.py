#!/usr/bin/env python3
"""Validate the scientific role registry and analysis contract for v2 workflows."""

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


def main() -> None:
    registry = yaml.safe_load(REGISTRY.read_text(encoding="utf-8"))
    contract = yaml.safe_load(CONTRACT.read_text(encoding="utf-8"))
    current = registry["canonical"]["current"]
    current_path = current["workflow"]
    current_text = read_text(current_path)
    inla_text = INLA_SCRIPT.read_text(encoding="utf-8")
    tier = str(current.get("evidence_tier", ""))

    expected_path = f"bombus_join/{tier}/"
    expected_arg = f"--analysis-tier {tier}"
    if expected_path not in current_text:
        fail(f"current workflow does not explicitly read registry tier {tier!r}: {current_path}")
    if expected_arg not in current_text:
        fail(f"current workflow does not pass {expected_arg!r}: {current_path}")

    if registry["guardrails"].get("exclude_zero_distance_from_fitted_models"):
        markers = (
            "distance_to_continent_km'].gt(0)",
            "distance_to_continent_km'] > 0",
            'distance_to_continent_km"] > 0',
        )
        if not any(marker in current_text for marker in markers):
            fail("current workflow does not explicitly construct positive-distance model support")

    if registry["guardrails"].get("category_engine") != "INLA":
        fail("category_engine must be INLA")
    if current.get("engine") != "INLA":
        fail("canonical category-preserving workflow must use INLA")
    if "run_inla_category_preserving_north.R" not in current_text:
        fail("canonical workflow does not run the category-preserving INLA analysis")

    brms_role = registry.get("engine_roles", {}).get("bayesian_brms", {})
    if registry["guardrails"].get("prohibit_brms_category_models"):
        required_exclusions = {
            "category_preserving_flower_colour",
            "category_preserving_floral_form",
            "multinomial_category_models",
        }
        if not required_exclusions.issubset(set(brms_role.get("excludes", []))):
            fail("brms role must explicitly exclude category-preserving/multinomial models")

    for relative in (
        registry["guardrails"]["require_analysis_contract"],
        registry["guardrails"]["require_input_validation"],
    ):
        if not (ROOT / relative).is_file():
            fail(f"required contract file is missing: {relative}")

    canonical_contract = contract.get("canonical_analysis", {})
    if canonical_contract.get("workflow") != current_path:
        fail("analysis contract and registry disagree on canonical workflow")
    if canonical_contract.get("engine") != current.get("engine"):
        fail("analysis contract and registry disagree on canonical engine")
    if canonical_contract.get("current_evidence_tier") != tier:
        fail("analysis contract and registry disagree on current evidence tier")
    if canonical_contract.get("planned_confirmatory_evidence_tier") != "primary":
        fail("planned confirmatory evidence tier must remain primary")

    guardrails = registry["guardrails"]
    if guardrails.get("planned_confirmatory_evidence_tier") != "primary":
        fail("registry planned confirmatory evidence tier must remain primary")
    if guardrails.get("use_equation_specific_maximum_support") is not True:
        fail("registry must require equation-specific maximum support")
    if guardrails.get("prohibit_alternative_pollinator_primary_covariates") is not True:
        fail("registry must prohibit alternative-pollinator primary covariates")

    contract_models = set(contract.get("models", {}))
    if not set(current.get("models", [])).issubset(contract_models):
        fail("registry canonical models are not all defined in the analysis contract")
    if "global_falsification" not in contract_models:
        fail("analysis contract must define the global falsification model")

    principles = contract.get("principles", {})
    required_true = {
        "exclude_zero_distance_from_fitted_models",
        "use_equation_specific_maximum_support",
        "compare_models_only_on_matched_support",
        "verify_category_counts_against_trials",
    }
    for name in required_true:
        if principles.get(name) is not True:
            fail(f"analysis contract must require {name}")
    if principles.get("category_omission_allowed") is not False:
        fail("analysis contract must prohibit category omission")
    if principles.get("alternative_pollinator_covariates_in_primary_models") is not False:
        fail("alternative pollinator covariates must remain outside primary models")

    forbidden = {
        "alternative_pollinator_guild",
        "showy_alt_guild",
        "other_bee_guild",
        "generalist_insect_guild",
    }
    global_predictors = set(contract["models"]["global_falsification"].get("predictors", []))
    if global_predictors & forbidden:
        fail("global falsification model must not include alternative-pollinator guild covariates")

    retained = contract.get("retained_categories", {})
    expected_colours = {"plain", "yellow_orange", "red_pink", "blue_purple"}
    expected_forms = {
        "open_generalized",
        "tubular_trumpet",
        "zygomorphic_specialized",
        "composite_brush",
    }
    contract_colours = set(retained.get("flower_color", []))
    contract_forms = set(retained.get("floral_form", []))
    if contract_colours != expected_colours:
        fail(f"retained flower-colour categories changed or are incomplete: {sorted(contract_colours)}")
    if contract_forms != expected_forms:
        fail(f"retained floral-form categories changed or are incomplete: {sorted(contract_forms)}")

    required_script_tokens = {
        "plain": 'plain = c("color_plain", "color_trials")',
        "yellow_orange": 'yellow_orange = c("color_yellow_orange", "color_trials")',
        "red_pink": 'red_pink = c("color_red_pink", "color_trials")',
        "blue_purple": 'blue_purple = c("color_blue_purple", "color_trials")',
        "open_generalized": 'open_generalized = c("form_open_generalized", "form_trials")',
        "tubular_trumpet": 'tubular_trumpet = c("form_tubular_trumpet", "form_trials")',
        "zygomorphic_specialized": 'zygomorphic_specialized = c("form_zygomorphic_specialized", "form_trials")',
        "composite_brush": 'composite_brush = c("form_composite_brush", "form_trials")',
    }
    missing_script_categories = [
        category for category, token in required_script_tokens.items() if token not in inla_text
    ]
    if missing_script_categories:
        fail(
            "canonical INLA script omits retained categories: "
            + ", ".join(sorted(missing_script_categories))
        )

    if "retained colour categories do not sum to color_trials" not in inla_text:
        fail("canonical INLA script does not verify colour category totals")
    if "retained floral-form categories do not sum to form_trials" not in inla_text:
        fail("canonical INLA script does not verify floral-form category totals")
    if "colour categories do not exhaust color_trials" not in current_text:
        fail("workflow does not verify colour category totals before fitting")
    if "form categories do not exhaust form_trials" not in current_text:
        fail("workflow does not verify floral-form category totals before fitting")

    expected_all = expected_colours | expected_forms
    n3_categories = set(contract["models"]["N3_direct_indirect"].get("applies_to", []))
    if n3_categories != expected_all:
        fail("N3 direct/indirect paths must apply to every retained colour/form category")

    global_responses = set(contract["models"]["global_falsification"].get("responses", []))
    if not ({"self_compatibility"} | expected_all).issubset(global_responses):
        fail("global falsification must include SC and every retained colour/form category")

    registered: set[str] = {current_path}
    for group in registry.get("robustness", {}).values():
        registered.update(group.get("workflows", []))
    for group in registry.get("replication", {}).values():
        workflow = group.get("workflow")
        if workflow:
            registered.add(workflow)
    registered.update(registry.get("legacy", {}).get("workflows", []))

    missing = [relative for relative in sorted(registered) if not (ROOT / relative).is_file()]
    if missing:
        fail("registered workflow files are missing: " + ", ".join(missing))

    print(
        {
            "registry": str(REGISTRY.relative_to(ROOT)),
            "contract": str(CONTRACT.relative_to(ROOT)),
            "current_workflow": current_path,
            "current_engine": current.get("engine"),
            "current_evidence_tier": tier,
            "retained_colour_categories": sorted(contract_colours),
            "retained_form_categories": sorted(contract_forms),
            "category_totals_verified": True,
            "registered_workflows": len(registered),
            "status": "ok",
        }
    )


if __name__ == "__main__":
    main()
