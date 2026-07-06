#!/usr/bin/env python3
"""Validate basic GitHub Actions workflow contracts without external tooling.

This deliberately catches the inexpensive, high-frequency failures that Python
unit tests cannot see: malformed YAML, missing workflow names/triggers, and
invalid schedule expressions. GitHub Actions cron uses exactly five fields.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import yaml

WORKFLOW_DIR = Path(".github/workflows")


def _trigger_section(document: dict[Any, Any]) -> Any:
    """Return the Actions trigger block despite PyYAML's YAML 1.1 `on` coercion."""
    return document.get("on", document.get(True))


def _validate_cron(value: Any, path: Path, errors: list[str]) -> None:
    if not isinstance(value, str) or not value.strip():
        errors.append(f"{path}: schedule cron must be a non-empty string")
        return
    fields = value.split()
    if len(fields) != 5:
        errors.append(
            f"{path}: GitHub Actions cron must have 5 fields, found {len(fields)}: {value!r}"
        )


def validate_workflow(path: Path) -> list[str]:
    errors: list[str] = []
    try:
        document = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        return [f"{path}: invalid YAML: {exc}"]

    if not isinstance(document, dict):
        return [f"{path}: top-level workflow document must be a mapping"]
    if not isinstance(document.get("name"), str) or not document["name"].strip():
        errors.append(f"{path}: workflow needs a non-empty name")

    triggers = _trigger_section(document)
    if triggers in (None, "", False):
        errors.append(f"{path}: workflow needs an `on` trigger")
        return errors
    if isinstance(triggers, str):
        return errors
    if not isinstance(triggers, dict):
        errors.append(f"{path}: `on` must be a trigger string or mapping")
        return errors

    schedule = triggers.get("schedule", [])
    if schedule and not isinstance(schedule, list):
        errors.append(f"{path}: schedule must be a list")
        return errors
    for item in schedule:
        if not isinstance(item, dict):
            errors.append(f"{path}: each schedule entry must be a mapping")
            continue
        _validate_cron(item.get("cron"), path, errors)
    return errors


def main() -> int:
    if not WORKFLOW_DIR.is_dir():
        print(f"Workflow directory does not exist: {WORKFLOW_DIR}", file=sys.stderr)
        return 2

    paths = sorted((*WORKFLOW_DIR.glob("*.yml"), *WORKFLOW_DIR.glob("*.yaml")))
    if not paths:
        print("No workflow files found.", file=sys.stderr)
        return 2

    errors = [error for path in paths for error in validate_workflow(path)]
    if errors:
        print("GitHub Actions workflow contract failures:", file=sys.stderr)
        print("\n".join(f"- {error}" for error in errors), file=sys.stderr)
        return 1

    print(f"Validated {len(paths)} GitHub Actions workflow file(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
