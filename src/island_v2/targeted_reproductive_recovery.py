"""Focused SC/SI and mating-system recovery using exact phrase searches."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

import island_v2.targeted_trait_recovery as targeted

app = typer.Typer(add_completion=False, no_args_is_help=True)

REPRODUCTIVE_QUERY_GROUPS = {
    "self_incompatibility": (
        '"{name}" "self-compatible"',
        '"{name}" "self-incompatible"',
        '"{name}" "self compatibility"',
        '"{name}" "self incompatibility"',
        '"{name}" "self-fertile" "self-sterile"',
    ),
    "mating_system": (
        '"{name}" autogamy autogamous selfing',
        '"{name}" outcrossing outcrossed',
        '"{name}" "mixed mating"',
        '"{name}" "mating system"',
        '"{name}" "breeding system"',
    ),
}

# Broaden biologically interpretable wording while retaining the compact user schema.
targeted.MATING_PATTERNS = {
    "obligate_outcrossing": (
        r"\bobligate(?:ly)? outcross",
        r"\bstrict(?:ly)? outcross",
        r"\bdioecious\b",
    ),
    "mixed_mating": (
        r"\bmixed[- ]mating\b",
        r"\bboth selfing and outcrossing\b",
        r"\bfacultative selfing\b",
        r"\bpartial(?:ly)? selfing\b",
        r"\bpredominantly outcrossing\b",
        r"\bmainly outcrossing\b",
    ),
    "mainly_selfing": (
        r"\bpredominantly selfing\b",
        r"\bmainly selfing\b",
        r"\bprimarily selfing\b",
        r"\bhigh selfing rate\b",
        r"\bautogamous\b",
        r"\bself[- ]pollinating\b",
        r"\bcleistogamous\b",
    ),
    "obligate_selfing": (
        r"\bobligate(?:ly)? selfing\b",
        r"\bstrict(?:ly)? selfing\b",
    ),
}

targeted.SI_PATTERNS = {
    "SI": (
        r"\bself[- ]incompatible\b",
        r"\bself[- ]sterile\b",
        r"\bself incompatibility\b",
        r"\bself-incompatibility\b",
    ),
    "SC": (
        r"\bself[- ]compatible\b",
        r"\bself[- ]fertile\b",
        r"\bself compatibility\b",
        r"\bself-compatibility\b",
        r"\bcapable of self[- ]fertili[sz]ation\b",
    ),
}


@app.command()
def collect(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    max_taxa: int = typer.Option(25, min=1, max=500),
    max_results_per_query: int = typer.Option(8, min=1, max=20),
    workers: int = typer.Option(16, min=1, max=32),
) -> None:
    targeted.QUERY_GROUPS = REPRODUCTIVE_QUERY_GROUPS
    species = pd.read_csv(species_csv, dtype=str).fillna("")
    frame, report = targeted.collect_targeted(
        species,
        max_taxa=max_taxa,
        max_synonyms=0,
        max_results_per_query=max_results_per_query,
        workers=workers,
    )
    report["query_mode"] = "reproductive_exact_phrase_species_plus_genus"
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "targeted_reproductive_recovery.csv", index=False)
    (output_dir / "targeted_reproductive_recovery_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
