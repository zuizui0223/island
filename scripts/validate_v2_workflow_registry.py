#!/usr/bin/env python3
"""Validate the scientific role registry for v2 analysis workflows."""

from __future__ import annotations

import sys
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
        if "distance_to_continent_km'].gt(0)" not in current_text:
            fail("current workflow does not explicitly construct positive-distance model support")

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

    replication = registry.get("replication", {}).get("bayesian_engine", {})
    replication_path = replication.get("workflow")
    if replication_path:
        text = read_text(replication_path)
        if "sensitivity_all" in text and replication.get("required_evidence_tier") == "primary":
            print(
                "WARNING: brms replication still references sensitivity_all; "
                "it is classified as replication, not the current INLA route.",
                file=sys.stderr,
            )

    print(
        {
            "registry": str(REGISTRY.relative_to(ROOT)),
            "current_workflow": current_path,
            "current_evidence_tier": tier,
            "registered_workflows": len(registered),
            "status": "ok",
        }
    )


if __name__ == "__main__":
    main()
