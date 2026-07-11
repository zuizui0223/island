"""Parallel execution wrapper for query-driven trait source discovery."""

from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import httpx
import pandas as pd
import typer

from island_v2.trait_dynamic_discovery import (
    SEARCH_COLUMNS,
    SOURCE_COLUMNS,
    USER_AGENT,
    _clean,
    discover_species,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _discover_one(
    species: str,
    max_synonyms: int,
    max_results_per_query: int,
    max_pages_per_species: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    headers = {"User-Agent": USER_AGENT, "Accept-Language": "en;q=0.9,*;q=0.5"}
    with httpx.Client(headers=headers) as client:
        return discover_species(
            client,
            species,
            max_synonyms=max_synonyms,
            max_results_per_query=max_results_per_query,
            max_pages_per_species=max_pages_per_species,
            delay_seconds=0.0,
        )


@app.command()
def discover(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    species_column: str = typer.Option("accepted_species"),
    max_taxa: int = typer.Option(25, min=1, max=500),
    max_synonyms: int = typer.Option(1, min=0, max=20),
    max_results_per_query: int = typer.Option(4, min=1, max=20),
    max_pages_per_species: int = typer.Option(8, min=1, max=50),
    workers: int = typer.Option(8, min=1, max=32),
) -> None:
    table = pd.read_csv(species_csv, dtype=str).fillna("")
    if species_column not in table.columns:
        raise typer.BadParameter(f"missing species column: {species_column}")
    species: list[str] = []
    seen: set[str] = set()
    for value in table[species_column]:
        name = _clean(value)
        if name and name not in seen:
            seen.add(name)
            species.append(name)
        if len(species) >= max_taxa:
            break

    search_rows: list[dict[str, str]] = []
    source_rows: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(
                _discover_one,
                name,
                max_synonyms,
                max_results_per_query,
                max_pages_per_species,
            ): name
            for name in species
        }
        for future in as_completed(futures):
            name = futures[future]
            try:
                searches, sources = future.result()
                search_rows.extend(searches)
                source_rows.extend(sources)
            except Exception as exc:
                errors.append({"accepted_species": name, "error": str(exc)})

    search_frame = pd.DataFrame(search_rows, columns=SEARCH_COLUMNS).drop_duplicates()
    source_frame = pd.DataFrame(source_rows, columns=SOURCE_COLUMNS).drop_duplicates(
        ["accepted_species", "source_url", "content_sha256"]
    )
    error_frame = pd.DataFrame(errors, columns=["accepted_species", "error"])
    output_dir.mkdir(parents=True, exist_ok=True)
    search_frame.to_csv(output_dir / "search_results.csv", index=False)
    source_frame.to_csv(output_dir / "source_texts.csv", index=False)
    error_frame.to_csv(output_dir / "discovery_errors.csv", index=False)
    summary = {
        "n_species_requested": len(species),
        "n_species_with_search_results": int(search_frame["accepted_species"].nunique()) if len(search_frame) else 0,
        "n_species_with_fetched_sources": int(source_frame["accepted_species"].nunique()) if len(source_frame) else 0,
        "n_search_results": len(search_frame),
        "n_fetched_sources": len(source_frame),
        "n_species_errors": len(error_frame),
        "workers": workers,
        "max_synonyms": max_synonyms,
        "max_results_per_query": max_results_per_query,
        "max_pages_per_species": max_pages_per_species,
        "source_type_counts": source_frame["source_type"].value_counts().to_dict() if len(source_frame) else {},
    }
    (output_dir / "dynamic_discovery_manifest.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    typer.echo(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    app()
