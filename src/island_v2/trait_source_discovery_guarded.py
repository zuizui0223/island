"""Strict CLI wrapper for free trait source discovery.

The legacy discovery functions are retained for offline unit testing, while this
entrypoint enforces the operational contract: a staged-taxa CSV must be for one
nominated Core-pilot island and must retain its curation release gate. It writes
unreviewed source leads only; run ``island-v2-trait-pdf-locator`` separately to
locate candidate pages in declared public PDFs.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

from island_v2 import trait_source_discovery as discovery

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_STAGED_COLUMNS = {
    "accepted_species",
    "island_id",
    "trait_candidate_status",
    "release_gate",
}


@app.callback()
def main() -> None:
    """Discover source leads only for one nominated Core-pilot island."""


def validate_staged_taxa_context(staged_taxa: pd.DataFrame, island_id: str) -> pd.DataFrame:
    """Require a single matching island and preserve all non-biological gates."""
    missing = REQUIRED_STAGED_COLUMNS.difference(staged_taxa.columns)
    if missing:
        raise typer.BadParameter(
            f"Staged Core-pilot taxa missing columns: {', '.join(sorted(missing))}"
        )
    table = staged_taxa.copy().fillna("")
    table["island_id"] = table["island_id"].astype(str).str.strip()
    table["accepted_species"] = table["accepted_species"].astype(str).str.strip()
    table = table.loc[table["accepted_species"].ne("")].copy()
    observed = sorted(set(table["island_id"]))
    if observed != [island_id]:
        raise typer.BadParameter(
            f"staged taxa must contain exactly island_id={island_id!r}; observed {observed}"
        )
    prohibited = discovery.PROHIBITED_OUTPUT_COLUMNS.intersection(table.columns)
    if prohibited:
        raise typer.BadParameter(
            "staged taxa contains prohibited finalized-value columns: "
            f"{', '.join(sorted(prohibited))}"
        )
    return table


def attach_release_gate(frame: pd.DataFrame, staged: pd.DataFrame) -> pd.DataFrame:
    """Attach the pre-existing curation gate without interpreting it."""
    lookup = (
        staged[["accepted_species", "release_gate"]]
        .drop_duplicates("accepted_species")
        .set_index("accepted_species")["release_gate"]
        .to_dict()
    )
    result = frame.copy()
    result["release_gate"] = result["query_taxon"].map(lookup).fillna("")
    return result


@app.command("discover")
def discover(
    staged_taxa_csv: Path = typer.Option(..., exists=True, help="One nominated island's staged taxa CSV."),
    nomination_report: Path = typer.Option(..., exists=True, help="core_pilot_nomination_report.json."),
    island_id: str = typer.Option(..., help="Island ID that must be nominated and match every staged row."),
    contact_email: str = typer.Option(..., help="Contact email for OpenAlex and Unpaywall polite pools."),
    output_dir: Path = typer.Option(..., help="Directory for unreviewed literature/access leads."),
    max_taxa: int = typer.Option(5, min=1, max=25, help="Bounded staged taxa to query."),
    max_seeds_per_taxon: int = typer.Option(5, min=1, max=20, help="Bounded seeds per taxon."),
) -> None:
    """Write M0 literature/access leads; no biological value is assigned."""
    discovery._require_nominated_island(nomination_report, island_id)
    staged = pd.read_csv(staged_taxa_csv, dtype=str).fillna("")
    staged = validate_staged_taxa_context(staged, island_id)
    frame, report = discovery.discover_sources(
        staged,
        discovery._httpx_getter(),
        contact_email,
        max_taxa,
        max_seeds_per_taxon,
    )
    frame = attach_release_gate(frame, staged)
    report.update(
        {
            "island_id": island_id,
            "release_boundary": (
                "Rows are unreviewed M0 literature/access leads. Neither source discovery nor PDF page "
                "location can accept a trait value, native/establishment status, Bombus applicability, or analysis inclusion."
            ),
        }
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "trait_literature_seeds.csv", index=False)
    receipt_columns = [
        column
        for column in [
            "query_taxon", "island_id", "doi", "source_api", "is_open_access", "oa_status",
            "open_access_url", "candidate_pdf_url", "lead_status", "provenance_note", "release_gate",
        ]
        if column in frame.columns
    ]
    frame.loc[:, receipt_columns].to_csv(output_dir / "trait_public_source_receipts.csv", index=False)
    # Backward-compatible name for existing downstream inspection scripts.
    frame.to_csv(output_dir / "trait_source_leads.csv", index=False)
    (output_dir / "trait_source_discovery_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"Wrote {report['n_literature_seeds']} unreviewed literature/access leads for "
        f"{report['n_taxa_queried']} staged taxa on {island_id}."
    )


if __name__ == "__main__":
    app()
