"""Validate structured flower-colour rows and flag species-level admissibility.

A single ``flower_primary_color`` value is insufficient: cultivar colours,
anthesis colour change, and geographic colour variation must not be merged into a
species-level primary colour. This module keeps the raw description and the
structured fields, validates the controlled vocabulary, and appends a
``species_level_primary_colour_admissible`` flag. Non-admissible rows are RETAINED
and flagged, never dropped and never silently promoted to a species colour.

This decides no colour value: it only validates and flags. A blank primary_colour
means unresolved, not "no colour".
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.callback()
def main() -> None:
    """Validate structured colour rows and flag species-level admissibility."""


def load_schema(path: Path) -> dict[str, Any]:
    schema = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(schema, dict) or "required_columns" not in schema:
        raise typer.BadParameter("colour schema must contain required_columns")
    return schema


def validate_colour(table: pd.DataFrame, schema: dict[str, Any]) -> list[str]:
    """Return a list of schema/vocabulary problems (empty means valid)."""
    issues: list[str] = []
    required = list(schema.get("required_columns", []))
    missing = [c for c in required if c not in table.columns]
    if missing:
        issues.append(f"missing required columns: {missing}")

    def _bad_values(column: str, allowed_key: str) -> None:
        if column not in table.columns:
            return
        allowed = {str(v) for v in schema.get(allowed_key, [])}
        seen = {str(v).strip() for v in table[column] if str(v).strip()}
        invalid = sorted(seen.difference(allowed))
        if invalid:
            issues.append(f"invalid {column} values: {invalid[:8]}")

    _bad_values("primary_colour", "primary_colour_terms")
    _bad_values("secondary_colour", "primary_colour_terms")
    _bad_values("within_species_variation", "within_species_variation")
    _bad_values("geographic_scope", "geographic_scope")
    _bad_values("wild_or_cultivated", "wild_or_cultivated")
    _bad_values("colour_evidence_type", "colour_evidence_type")
    return issues


def _row_admissible(row: pd.Series, excludes: dict[str, Any]) -> bool:
    """A row is species-level admissible unless it is a cultivar/variation artefact."""
    for column, bad_values in excludes.items():
        value = str(row.get(column, "")).strip()
        if value in {str(v) for v in bad_values}:
            return False
    return True


def annotate_admissibility(
    table: pd.DataFrame, schema: dict[str, Any]
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Append species_level_primary_colour_admissible; retain (never drop) rows."""
    result = table.copy()
    excludes = schema.get("species_level_primary_colour_excludes", {}) or {}
    admissible = [
        _row_admissible(row, excludes) for _, row in result.iterrows()
    ]
    result["species_level_primary_colour_admissible"] = admissible
    n_admissible = int(sum(1 for a in admissible if a))
    summary = {
        "note": (
            "Non-admissible rows (cultivated, anthesis colour change, geographic "
            "variation) are retained and flagged, not merged into a species colour."
        ),
        "n_rows": int(len(result)),
        "n_species_level_admissible": n_admissible,
        "n_flagged_non_admissible": int(len(result)) - n_admissible,
    }
    return result, summary


@app.command("validate")
def validate(
    colour_csv: Path = typer.Option(..., exists=True, help="Structured colour CSV."),
    schema_path: Path = typer.Option(Path("config/colour_schema.yml")),
) -> None:
    """Validate the structured colour table against the schema; exit non-zero on problems."""
    schema = load_schema(schema_path)
    table = pd.read_csv(colour_csv, dtype=str).fillna("")
    issues = validate_colour(table, schema)
    if issues:
        for issue in issues:
            typer.echo(f"- {issue}", err=True)
        raise typer.Exit(code=1)
    typer.echo(f"OK: {len(table)} structured colour row(s) conform to the schema.")


@app.command("annotate")
def annotate(
    colour_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    schema_path: Path = typer.Option(Path("config/colour_schema.yml")),
) -> None:
    """Flag species-level admissibility for each structured colour row (rows kept)."""
    schema = load_schema(schema_path)
    table = pd.read_csv(colour_csv, dtype=str).fillna("")
    issues = validate_colour(table, schema)
    if issues:
        raise typer.BadParameter("; ".join(issues))
    result, summary = annotate_admissibility(table, schema)
    output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_dir / "trait_colour_flagged.csv", index=False)
    (output_dir / "colour_admissibility_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"Flagged {summary['n_rows']} colour row(s): "
        f"{summary['n_species_level_admissible']} species-level admissible, "
        f"{summary['n_flagged_non_admissible']} non-admissible."
    )


if __name__ == "__main__":
    app()
