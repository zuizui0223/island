"""Validation commands used locally and by GitHub Actions."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)


def load_ontology(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or "traits" not in data:
        raise ValueError("Ontology must contain a top-level traits mapping.")
    return data


def validate_ontology(ontology: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    all_values: set[str] = set()
    for trait_name, spec in ontology["traits"].items():
        values = spec.get("allowed_values", [])
        if not values:
            errors.append(f"{trait_name}: allowed_values is empty")
        if len(values) != len(set(values)):
            errors.append(f"{trait_name}: allowed_values contains duplicates")
        all_values.update(values)
    for name, mapping in ontology.get("secondary_v1_summaries", {}).items():
        for group, values in mapping.items():
            unknown = set(values).difference(all_values)
            if unknown:
                errors.append(f"{name}.{group}: unknown ontology values {sorted(unknown)}")
    return errors


def validate_candidate_table(table: pd.DataFrame, ontology: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    required = {
        "accepted_species",
        "trait_name",
        "standardized_value",
        "candidate_kind",
        "evidence_scope",
        "source_type",
        "source_reliability",
        "needs_human_review",
        "review_status",
    }
    missing = required.difference(table.columns)
    if missing:
        return [f"candidate CSV missing columns: {sorted(missing)}"]

    allowed = {name: set(spec["allowed_values"]) for name, spec in ontology["traits"].items()}
    for index, row in table.iterrows():
        label = f"row {index + 2} ({row['accepted_species']})"
        trait_name = row["trait_name"]
        if trait_name not in allowed:
            errors.append(f"{label}: unknown trait_name {trait_name}")
            continue
        if row["standardized_value"] not in allowed[trait_name]:
            errors.append(f"{label}: invalid standardized_value {row['standardized_value']}")
        if str(row["needs_human_review"]).strip().lower() not in {"true", "1"}:
            errors.append(f"{label}: LLM candidate must require human review")
        if row["candidate_kind"] == "hierarchical_inference":
            for field in ("supporting_taxa", "inference_rule"):
                if field not in table.columns or not str(row.get(field, "")).strip():
                    errors.append(f"{label}: hierarchical inference lacks {field}")
        if row["candidate_kind"] == "source_backed" and row["source_type"] == "none_found":
            errors.append(f"{label}: source_backed candidate cannot have source_type none_found")
    return errors


def check_or_raise(errors: list[str]) -> None:
    if errors:
        raise typer.Exit(code=1)


@app.command("ontology")
def ontology_command(
    ontology_path: Path = typer.Option(Path("config/trait_ontology.yml")),
) -> None:
    """Validate the versioned trait ontology."""
    errors = validate_ontology(load_ontology(ontology_path))
    for error in errors:
        typer.echo(f"ERROR: {error}")
    if errors:
        raise typer.Exit(code=1)
    typer.echo("Ontology validation passed.")


@app.command("candidates")
def candidates_command(
    candidate_csv: Path = typer.Option(..., exists=True),
    ontology_path: Path = typer.Option(Path("config/trait_ontology.yml")),
) -> None:
    """Validate a web-search LLM candidate export before review."""
    table = pd.read_csv(candidate_csv, dtype=str).fillna("")
    errors = validate_candidate_table(table, load_ontology(ontology_path))
    for error in errors:
        typer.echo(f"ERROR: {error}")
    if errors:
        raise typer.Exit(code=1)
    typer.echo(f"Candidate validation passed ({len(table)} rows).")


@app.command("print-schema")
def print_schema() -> None:
    """Print a compact JSON description of the accepted curation columns."""
    print(json.dumps({"required_review_status": ["pending", "accepted", "rejected"]}, indent=2))


if __name__ == "__main__":
    app()
