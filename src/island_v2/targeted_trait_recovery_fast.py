"""Fast wrapper for targeted trait recovery using four compressed query families."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

import island_v2.targeted_trait_recovery as targeted

app = typer.Typer(add_completion=False, no_args_is_help=True)

FAST_QUERY_GROUPS = {
    "flower_color": ('"{name}" flower color colour corolla petals blooms',),
    "pollination": ('"{name}" pollinated pollinator bee fly butterfly moth bird wind',),
    "self_incompatibility": ('"{name}" self-compatible self-incompatible self-fertile self-sterile',),
    "mating_system": ('"{name}" autogamy selfing outcrossing mixed mating breeding system',),
}


@app.command()
def collect(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    max_taxa: int = typer.Option(25, min=1, max=500),
    max_synonyms: int = typer.Option(1, min=0, max=10),
    max_results_per_query: int = typer.Option(8, min=1, max=20),
    workers: int = typer.Option(8, min=1, max=32),
) -> None:
    targeted.QUERY_GROUPS = FAST_QUERY_GROUPS
    species = pd.read_csv(species_csv, dtype=str).fillna("")
    frame, report = targeted.collect_targeted(
        species,
        max_taxa=max_taxa,
        max_synonyms=max_synonyms,
        max_results_per_query=max_results_per_query,
        workers=workers,
    )
    report["query_mode"] = "four_compressed_trait_queries"
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "targeted_trait_recovery.csv", index=False)
    (output_dir / "targeted_trait_recovery_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
