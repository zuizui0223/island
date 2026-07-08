"""Derive size / tube-depth classes from raw measurements by a versioned rule.

The raw measurement is stored first; the class is derived here so it is
reproducible and reviewer-independent. The raw fields
(``raw_measurement, raw_unit, measurement_structure, measurement_source_text``)
are never overwritten; this only appends ``derived_class`` and
``classification_rule_version``. A measurement that cannot be parsed or has an
unknown unit yields a blank class (unresolved), never a guessed one.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

# Measurements are non-negative; a hyphen denotes a range, not a negative sign.
_NUMBER = re.compile(r"\d*\.?\d+")


@app.callback()
def main() -> None:
    """Derive size/tube-depth classes from raw measurements (raw kept)."""


def load_rules(path: Path) -> dict[str, Any]:
    rules = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(rules, dict) or "traits" not in rules:
        raise typer.BadParameter("measurement classification config must contain a traits mapping")
    return rules


def _to_mm(raw_measurement: str, raw_unit: str, unit_to_mm: dict[str, float]) -> float | None:
    """Parse a raw measurement to millimetres; a range uses its midpoint."""
    numbers = [float(n) for n in _NUMBER.findall(str(raw_measurement or ""))]
    if not numbers:
        return None
    value = sum(numbers) / len(numbers)  # midpoint of a range, or the single value
    factor = unit_to_mm.get(str(raw_unit or "").strip().lower())
    if factor is None:
        return None
    return value * factor


def derive_class(
    trait: str,
    raw_measurement: str,
    raw_unit: str,
    rules: dict[str, Any],
) -> str:
    """Return the derived class for a size/tube trait, or '' if not derivable."""
    trait_rule = rules.get("traits", {}).get(str(trait).strip())
    if not trait_rule:
        return ""
    mm = _to_mm(raw_measurement, raw_unit, rules.get("unit_to_mm", {}))
    if mm is None:
        return ""
    for band in trait_rule.get("bands", []):
        max_mm = band.get("max_mm")
        if max_mm is None or mm <= float(max_mm):
            return str(band.get("derived_class", ""))
    return ""


def annotate_measurements(table: pd.DataFrame, rules: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Append derived_class + classification_rule_version; never overwrite raw."""
    trait_col = next((c for c in ("trait_name", "trait") if c in table.columns), None)
    if trait_col is None:
        raise typer.BadParameter("measurement table must have a 'trait_name' or 'trait' column")
    for required in ("raw_measurement", "raw_unit"):
        if required not in table.columns:
            raise typer.BadParameter(f"measurement table missing required column: {required}")

    result = table.copy()
    version = str(rules.get("classification_rule_version", ""))
    derived = [
        derive_class(str(t), str(m), str(u), rules)
        for t, m, u in zip(result[trait_col], result["raw_measurement"], result["raw_unit"], strict=False)
    ]
    result["derived_class"] = derived
    result["classification_rule_version"] = version
    n_resolved = int(sum(1 for d in derived if d))
    summary = {
        "note": "Class derived from raw measurement by a versioned rule; raw fields are unchanged.",
        "classification_rule_version": version,
        "n_rows": int(len(result)),
        "n_class_derived": n_resolved,
        "n_unresolved": int(len(result)) - n_resolved,
    }
    return result, summary


@app.command("annotate")
def annotate(
    measurements_csv: Path = typer.Option(..., exists=True, help="Raw measurement CSV (raw_measurement, raw_unit, trait)."),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/measurement_classification.yml")),
) -> None:
    """Write measurement rows with a derived class and the rule version."""
    rules = load_rules(config_path)
    table = pd.read_csv(measurements_csv, dtype=str).fillna("")
    result, summary = annotate_measurements(table, rules)
    output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_dir / "trait_measurements_classified.csv", index=False)
    (output_dir / "measurement_classification_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"Classified {summary['n_rows']} measurement(s) with rule "
        f"{summary['classification_rule_version']}: {summary['n_class_derived']} derived, "
        f"{summary['n_unresolved']} unresolved."
    )


if __name__ == "__main__":
    app()
