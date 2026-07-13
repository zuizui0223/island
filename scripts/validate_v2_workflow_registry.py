#!/usr/bin/env python3
"""Validate the scientific role registry for v2 analysis workflows."""

from __future__ import annotations

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "config" / "v2_workflow_registry.yml"


def fail(message: str) -> None:
    raise SystemExit(f"v2 workflow registry validation failed: {message}")


def read_text(relative: str) -> str:
    path = ROOT / relative
    if not path.is_file():
        fail(f"registered workflow does not exist: {relative}")
    return path.read_text(encoding="utf-8")


def main() -> None:
    registry = yaml.safe_load(REGISTRY.read_text(encoding="utf-8"))
    current = registry["canonical"]["current"]
    current_path = current["workflow"]
    current_text = read_text(current_path)
    tier = str(current.get("evidence_tier", ""))

    expected_path = f"bombus_join/{tier}/"
    expected_arg = f"--analysis-tier {tier}"
    if expected_path not in current_text:
        fail(f"current workflow does not explicitly read registry tier {tier!r}: {current_path}")
    if expected_arg not in current_text:
        fail(f"current workflow does not pass {expected_arg!r}: {current_path}")

    if registry["guardrails"].get("exclude_zero_distance_from_fitted_models"):
        positive_distance_markers = (
            "distance_to_continent_km'].gt(0)",
            "distance_to_continent_km'] > 0",
            'distance_to_continent_km\"] > 0',
        )
        if not any(marker in current_text for marker in positive_distance_markers):
            fail("current workflow does not explicitly construct positive-distance model support")

    if registry["guardrails"].get("category_engine") != "INLA":
        fail("category_engine must be INLA")
    if current.get("engine") != "INLA":
        fail("canonical category-preserving workflow must use INLA")
    if "run_inla_category_preserving_north.R" not in current_text:
        fail("canonical workflow does not run the category-preserving INLA analysis")

    brms_role = registry.get("engine_roles", {}).get("bayesian_brms", {})
    if registry["guardrails"].get("prohibit_brms_category_models"):
        excluded = set(brms_role.get("excludes", []))
        required_exclusions = {
            "category_preserving_flower_colour",
            "category_preserving_floral_form",
            "multinomial_category_models",
        }
        if not required_exclusions.issubset(excluded):
            fail("brms role must explicitly exclude category-preserving/multinomial models")

    required_contract = registry["guardrails"]["require_analysis_contract"]
    required_validator = registry["guardrails"]["require_input_validation"]
    for relative in (required_contract, required_validator):
        if not (ROOT / relative).is_file():
            fail(f"required contract file is missing: {relative}")

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
            "current_workflow": current_path,
            "current_engine": current.get("engine"),
            "current_evidence_tier": tier,
            "category_engine": registry["guardrails"].get("category_engine"),
            "brms_category_models_allowed": not registry["guardrails"].get("prohibit_brms_category_models", False),
            "registered_workflows": len(registered),
            "status": "ok",
        }
    )


if __name__ == "__main__":
    main()
