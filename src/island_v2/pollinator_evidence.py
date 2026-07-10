"""Validate the long-format pollinator functional-evidence table and derive an index.

The raw evidence is kept **long** - one row per (species, pollinator_guild,
observation) - and is never collapsed to a single ``mixed_or_generalist`` bucket
at the raw stage. A functional-replacement index may be DERIVED from accepted
rows, but the derivation is separate from and does not mutate the raw table.

This module validates the schema and vocabulary and builds the derived index. It
does not decide any trait value or analysis inclusion; only ``accepted`` rows
feed the derived index, and everything else is reported, not silently dropped.
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
    """Validate and index long-format pollinator functional evidence."""


def load_schema(path: Path) -> dict[str, Any]:
    schema = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(schema, dict) or "required_columns" not in schema:
        raise typer.BadParameter("pollinator evidence schema must contain required_columns")
    return schema


def validate_long_format(table: pd.DataFrame, schema: dict[str, Any]) -> list[str]:
    """Return a list of schema/vocabulary problems (empty means valid)."""
    issues: list[str] = []
    required = list(schema.get("required_columns", []))
    missing = [c for c in required if c not in table.columns]
    if missing:
        issues.append(f"missing required columns: {missing}")
    forbidden = [c for c in schema.get("forbidden_columns", []) if c in table.columns]
    if forbidden:
        issues.append(f"raw table must not collapse guilds; forbidden columns present: {forbidden}")

    def _bad_values(column: str, allowed_key: str) -> None:
        if column not in table.columns:
            return
        allowed = {str(v) for v in schema.get(allowed_key, [])}
        seen = {str(v).strip() for v in table[column] if str(v).strip()}
        invalid = sorted(seen.difference(allowed))
        if invalid:
            issues.append(f"invalid {column} values: {invalid[:8]}")

    _bad_values("evidence_type", "evidence_type")
    _bad_values("effectiveness_class", "effectiveness_class")
    _bad_values("wild_or_cultivated", "wild_or_cultivated")
    _bad_values("review_status", "review_status")
    return issues


def functional_replacement_index(table: pd.DataFrame, schema: dict[str, Any]) -> pd.DataFrame:
    """Derive a per-species guild summary from ACCEPTED rows (raw kept intact).

    Never collapses the raw evidence: this is a separate derived table recording,
    per species, how many distinct guilds have accepted evidence and how many are
    confirmed effective pollinators - the signal for functional replacement.
    """
    columns = [
        "accepted_species",
        "n_guilds_with_accepted_evidence",
        "guilds_with_accepted_evidence",
        "n_confirmed_effective_guilds",
        "confirmed_effective_guilds",
        "strongest_evidence_type",
        "functional_replacement_candidate",
    ]
    if table.empty:
        return pd.DataFrame(columns=columns)

    strength = {str(v): i for i, v in enumerate(schema.get("evidence_type", []))}
    accepted = table.loc[table.get("review_status", "").astype(str).str.strip() == "accepted"].copy()
    if accepted.empty:
        return pd.DataFrame(columns=columns)

    accepted["_species"] = accepted["accepted_species"].astype(str).str.strip()
    accepted["_guild"] = accepted["pollinator_guild"].astype(str).str.strip()
    accepted["_evidence"] = accepted["evidence_type"].astype(str).str.strip()
    accepted["_effective"] = accepted["effectiveness_class"].astype(str).str.strip().eq(
        "confirmed_effective_pollinator"
    )
    accepted["_strength"] = accepted["_evidence"].map(strength).fillna(-1)

    rows: list[dict[str, Any]] = []
    for species, group in accepted.groupby("_species", sort=True):
        guilds = sorted({g for g in group["_guild"] if g})
        confirmed = sorted({g for g, ok in zip(group["_guild"], group["_effective"], strict=True) if ok and g})
        best_idx = int(group["_strength"].max())
        strongest = next((k for k, v in strength.items() if v == best_idx), "")
        rows.append(
            {
                "accepted_species": species,
                "n_guilds_with_accepted_evidence": len(guilds),
                "guilds_with_accepted_evidence": "|".join(guilds),
                "n_confirmed_effective_guilds": len(confirmed),
                "confirmed_effective_guilds": "|".join(confirmed),
                "strongest_evidence_type": strongest,
                # A functional-replacement candidate has 2+ guilds each confirmed as
                # an effective pollinator -- not a single generalist label.
                "functional_replacement_candidate": len(confirmed) >= 2,
            }
        )
    return pd.DataFrame(rows, columns=columns)


@app.command("validate")
def validate(
    evidence_csv: Path = typer.Option(..., exists=True, help="Long-format pollinator evidence CSV."),
    schema_path: Path = typer.Option(Path("config/pollinator_evidence_schema.yml")),
) -> None:
    """Validate the long-format table against the schema; exit non-zero on problems."""
    schema = load_schema(schema_path)
    table = pd.read_csv(evidence_csv, dtype=str).fillna("")
    issues = validate_long_format(table, schema)
    if issues:
        for issue in issues:
            typer.echo(f"- {issue}", err=True)
        raise typer.Exit(code=1)
    typer.echo(f"OK: {len(table)} long-format pollinator-evidence rows conform to the schema.")


@app.command("index")
def index(
    evidence_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    schema_path: Path = typer.Option(Path("config/pollinator_evidence_schema.yml")),
) -> None:
    """Derive the functional-replacement index from accepted rows (raw kept intact)."""
    schema = load_schema(schema_path)
    table = pd.read_csv(evidence_csv, dtype=str).fillna("")
    issues = validate_long_format(table, schema)
    if issues:
        raise typer.BadParameter("; ".join(issues))
    derived = functional_replacement_index(table, schema)
    output_dir.mkdir(parents=True, exist_ok=True)
    derived.to_csv(output_dir / "pollinator_functional_replacement_index.csv", index=False)
    summary = {
        "note": "Derived from accepted rows only; the raw long-format table is not modified.",
        "n_species": int(len(derived)),
        "n_functional_replacement_candidates": int(
            derived["functional_replacement_candidate"].sum()
        ) if len(derived) else 0,
    }
    (output_dir / "pollinator_index_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    typer.echo(
        f"Derived guild index for {summary['n_species']} species; "
        f"{summary['n_functional_replacement_candidates']} functional-replacement candidate(s)."
    )


if __name__ == "__main__":
    app()
