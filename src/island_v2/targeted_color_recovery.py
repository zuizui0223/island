"""Focused flower-colour recovery with species/genus fallback and flower-context filtering."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

import island_v2.targeted_trait_recovery as targeted

app = typer.Typer(add_completion=False, no_args_is_help=True)

COLOR_QUERY_GROUPS = {
    "flower_color": (
        '"{name}" "flower color"',
        '"{name}" "flower colour"',
        '"{name}" flowers corolla petals',
        '"{name}" bloom flower',
    )
}

FLOWER_CONTEXT_TERMS = ("flower", "flowers", "floral", "corolla", "petal", "petals", "bloom", "blooms")


@app.command()
def collect(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    max_taxa: int = typer.Option(25, min=1, max=500),
    max_synonyms: int = typer.Option(1, min=0, max=10),
    max_results_per_query: int = typer.Option(8, min=1, max=20),
    workers: int = typer.Option(16, min=1, max=32),
) -> None:
    targeted.QUERY_GROUPS = COLOR_QUERY_GROUPS
    species = pd.read_csv(species_csv, dtype=str).fillna("")
    frame, report = targeted.collect_targeted(
        species,
        max_taxa=max_taxa,
        max_synonyms=max_synonyms,
        max_results_per_query=max_results_per_query,
        workers=workers,
    )
    if not frame.empty:
        context = (
            frame["result_title"].astype(str).str.cat(frame["result_snippet"].astype(str), sep=" ").str.casefold()
        )
        keep = context.map(lambda text: any(term in text for term in FLOWER_CONTEXT_TERMS))
        frame = frame.loc[keep].copy()
    coverage = int(frame["accepted_species"].nunique()) if len(frame) else 0
    report["coverage"] = {"flower_color": coverage}
    report["n_rows"] = int(len(frame))
    report["query_mode"] = "exact_color_species_synonym_genus_with_flower_context"
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "targeted_color_recovery.csv", index=False)
    (output_dir / "targeted_color_recovery_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
