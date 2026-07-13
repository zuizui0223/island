#!/usr/bin/env python3
"""Validate the scientific role registry for v2 analysis workflows.

This is intentionally lightweight: it prevents workflow-role drift without trying to
parse every shell command as a full programming language.
"""

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
    canonical = registry["canonical"]["primary"]
    primary_path = canonical["workflow"]
    primary_text = read_text(primary_path)

    if canonical.get("evidence_tier") != "primary":
        fail("canonical primary evidence_tier must be 'primary'")
    if "bombus_join/primary/" not in primary_text:
        fail(f"canonical workflow does not explicitly read primary tier: {primary_path}")
    if "bombus_join/sensitivity_all/" in primary_text:
        fail(f"canonical workflow reads sensitivity_all: {primary_path}")
    if "--analysis-tier primary" not in primary_text:
        fail(f"canonical workflow does not pass --analysis-tier primary: {primary_path}")

    required_contract = registry["guardrails"]["require_analysis_contract"]
    required_validator = registry["guardrails"]["require_input_validation"]
    for relative in (required_contract, required_validator):
        if not (ROOT / relative).is_file():
            fail(f"required contract file is missing: {relative}")

    registered: set[str] = {primary_path}
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
                "it is classified as replication, not canonical primary inference.",
                file=sys.stderr,
            )

    print(
        {
            "registry": str(REGISTRY.relative_to(ROOT)),
            "canonical_primary": primary_path,
            "registered_workflows": len(registered),
            "status": "ok",
        }
    )


if __name__ == "__main__":
    main()
