"""Build a free, reproducible review packet for staged trait evidence work.

This module deliberately performs no model call and no external API request. It
turns a gated pilot taxon list into a taxon-by-trait review table that can be
filled from transparent web and literature searches carried out in the research
workflow. The packet is not a trait matrix and does not contain biological
claims.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_TAXON_COLUMNS = {"accepted_species", "genus", "family"}
REQUIRED_TRAIT_KEYS = {"domain", "allowed_values"}

OUTPUT_COLUMNS = [
    "accepted_species",
    "genus",
    "family",
    "island_id",
    "canonical_island_name",
    "country_iso3",
    "source_region_id",
    "n_records",
    "trait_candidate_status",
    "release_gate",
    "trait_name",
    "trait_domain",
    "allowed_values_json",
    "recommended_search_query",
    "evidence_status",
    "provisional_standardized_value",
    "source_citation",
    "source_url",
    "evidence_excerpt",
    "evidence_scope",
    "source_type",
    "source_reliability",
    "candidate_kind",
    "review_status",
    "review_note",
]


@app.callback()
def main() -> None:
    """Create an API-free trait evidence review packet from staged taxa."""


def _text(values: pd.Series) -> pd.Series:
    return values.fillna("").astype(str).str.strip()


def load_ontology(path: Path) -> dict[str, dict[str, Any]]:
    """Read and validate the versioned trait ontology."""
    if not path.exists():
        raise typer.BadParameter(f"Ontology file does not exist: {path}")
    document = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    traits = document.get("traits")
    if not isinstance(traits, dict) or not traits:
        raise typer.BadParameter("Trait ontology must contain a non-empty traits mapping")
    for trait_name, definition in traits.items():
        if not isinstance(definition, dict):
            raise typer.BadParameter(f"Trait {trait_name} must be a mapping")
        missing = REQUIRED_TRAIT_KEYS.difference(definition)
        if missing:
            raise typer.BadParameter(f"Trait {trait_name} missing keys: {', '.join(sorted(missing))}")
        if not isinstance(definition["allowed_values"], list):
            raise typer.BadParameter(f"Trait {trait_name} allowed_values must be a list")
    return traits


def build_review_packet(taxa: pd.DataFrame, traits: dict[str, dict[str, Any]]) -> pd.DataFrame:
    """Expand staged taxa into an empty, evidence-first taxon-by-trait packet."""
    missing = REQUIRED_TAXON_COLUMNS.difference(taxa.columns)
    if missing:
        raise typer.BadParameter(f"Taxon input missing columns: {', '.join(sorted(missing))}")
    work = taxa.copy()
    for column in work.columns:
        work[column] = _text(work[column])
    work = work.loc[work["accepted_species"].ne("")].drop_duplicates("accepted_species")
    if work.empty:
        raise typer.BadParameter("Taxon input has no nonblank accepted_species rows")

    rows: list[dict[str, str]] = []
    for taxon in work.sort_values("accepted_species").to_dict("records"):
        species = str(taxon["accepted_species"])
        for trait_name, definition in traits.items():
            query = f'"{species}" {trait_name.replace("_", " ")} flora pollination'
            rows.append(
                {
                    "accepted_species": species,
                    "genus": str(taxon.get("genus", "")),
                    "family": str(taxon.get("family", "")),
                    "island_id": str(taxon.get("island_id", "")),
                    "canonical_island_name": str(taxon.get("canonical_island_name", "")),
                    "country_iso3": str(taxon.get("country_iso3", "")),
                    "source_region_id": str(taxon.get("source_region_id", "")),
                    "n_records": str(taxon.get("n_records", "")),
                    "trait_candidate_status": str(taxon.get("trait_candidate_status", "staged_not_curated")),
                    "release_gate": str(taxon.get("release_gate", "")),
                    "trait_name": trait_name,
                    "trait_domain": str(definition["domain"]),
                    "allowed_values_json": json.dumps(definition["allowed_values"], ensure_ascii=False),
                    "recommended_search_query": query,
                    "evidence_status": "not_searched",
                    "provisional_standardized_value": "",
                    "source_citation": "",
                    "source_url": "",
                    "evidence_excerpt": "",
                    "evidence_scope": "",
                    "source_type": "",
                    "source_reliability": "",
                    "candidate_kind": "",
                    "review_status": "pending",
                    "review_note": "No biological value asserted. Fill only from a cited source; unresolved is valid.",
                }
            )
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


def packet_summary(packet: pd.DataFrame) -> dict[str, Any]:
    """Describe the work queue without treating blank values as biological absence."""
    return {
        "n_staged_taxa": int(packet["accepted_species"].nunique()),
        "n_trait_rows": int(len(packet)),
        "n_traits": int(packet["trait_name"].nunique()),
        "n_islands": int(packet["island_id"].replace("", pd.NA).dropna().nunique()),
        "evidence_status": "not_searched",
        "curation_rule": "Rows are review templates only. Blank values and unresolved values are not biological absences.",
    }


@app.command("build")
def build(
    taxa_csv: Path = typer.Option(..., exists=True, help="Staged pilot taxa CSV."),
    output_dir: Path = typer.Option(..., help="Directory for API-free review packet artifacts."),
    ontology_path: Path = typer.Option(
        Path("config/trait_ontology.yml"), help="Versioned trait ontology."
    ),
) -> None:
    """Write an API-free trait review packet and summary."""
    taxa = pd.read_csv(taxa_csv, dtype=str).fillna("")
    packet = build_review_packet(taxa, load_ontology(ontology_path))
    output_dir.mkdir(parents=True, exist_ok=True)
    packet.to_csv(output_dir / "trait_evidence_review_packet.csv", index=False)
    summary = packet_summary(packet)
    (output_dir / "trait_evidence_review_packet_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"Built {summary['n_trait_rows']} review rows for {summary['n_staged_taxa']} staged taxa."
    )


if __name__ == "__main__":
    app()
