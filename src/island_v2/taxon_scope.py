"""Construct the reviewed angiosperm analysis universe for v2.

The global GBIF campaign intentionally retrieves Tracheophyta as a broad raw
candidate universe. That is appropriate for occurrence recovery, taxonomy
audit, and source-pool documentation, but it cannot be used unchanged as the
denominator of a pollen-vector or floral-trait proportion: ferns and lycophytes
do not have flowers or pollen-vector states comparable to angiosperms.

This module turns a reviewed species-level scope table into an explicit
island-by-species analysis universe. It does not infer that unreviewed taxa are
angiosperms; those rows remain unresolved and excluded from confirmatory
angiosperm denominators until reviewed.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_ISLAND_SPECIES_COLUMNS = {"island_id", "species"}
REQUIRED_SCOPE_COLUMNS = {
    "accepted_species",
    "taxonomic_group",
    "decision_basis",
    "source_citation",
    "source_url",
    "review_status",
    "reviewer",
    "decision_date",
}
VALID_TAXONOMIC_GROUPS = {
    "angiosperm",
    "gymnosperm",
    "non_seed_vascular",
    "unresolved",
}
OUTPUT_COLUMNS = [
    "island_id",
    "species",
    "taxonomic_group",
    "taxon_scope_review_status",
    "taxon_scope_decision_basis",
    "taxon_scope_source_citation",
    "taxon_scope_source_url",
    "taxon_scope_reviewer",
    "taxon_scope_decision_date",
    "angiosperm_analysis_eligible",
]


@app.callback()
def main() -> None:
    """Build an explicit angiosperm species universe before composition models."""


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _require_columns(table: pd.DataFrame, required: set[str], label: str) -> None:
    missing = required.difference(table.columns)
    if missing:
        raise typer.BadParameter(f"{label} missing columns: {', '.join(sorted(missing))}")


def _unique_nonblank(table: pd.DataFrame, column: str, label: str) -> None:
    values = _text(table[column])
    if values.eq("").any():
        raise typer.BadParameter(f"{label} has blank {column} values")
    duplicated = values.loc[values.duplicated()].unique().tolist()
    if duplicated:
        raise typer.BadParameter(f"{label} has duplicate {column} values: {duplicated[:5]}")


def build_taxon_scope(
    island_species: pd.DataFrame,
    scope_decisions: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return full scoped island taxa and the accepted-angiosperm subset.

    A taxon enters the angiosperm analysis denominator only when its decision is
    both `accepted` and `angiosperm`. Non-seed vascular taxa, gymnosperms, and
    all unresolved or not-yet-reviewed taxa are kept in the audit output but
    excluded from the confirmatory angiosperm table.
    """
    _require_columns(island_species, REQUIRED_ISLAND_SPECIES_COLUMNS, "island species table")
    _require_columns(scope_decisions, REQUIRED_SCOPE_COLUMNS, "taxon scope decisions")
    _unique_nonblank(scope_decisions, "accepted_species", "taxon scope decisions")

    raw = island_species.copy()
    raw["island_id"] = _text(raw["island_id"])
    raw["species"] = _text(raw["species"])
    raw = raw.loc[raw["island_id"].ne("") & raw["species"].ne("")].drop_duplicates(
        ["island_id", "species"]
    )

    decisions = scope_decisions.copy()
    for column in REQUIRED_SCOPE_COLUMNS:
        decisions[column] = _text(decisions[column])
    invalid_groups = set(decisions["taxonomic_group"]).difference(VALID_TAXONOMIC_GROUPS)
    if invalid_groups:
        raise typer.BadParameter(f"Invalid taxonomic_group values: {sorted(invalid_groups)}")
    decisions = decisions.rename(
        columns={
            "accepted_species": "species",
            "review_status": "taxon_scope_review_status",
            "decision_basis": "taxon_scope_decision_basis",
            "source_citation": "taxon_scope_source_citation",
            "source_url": "taxon_scope_source_url",
            "reviewer": "taxon_scope_reviewer",
            "decision_date": "taxon_scope_decision_date",
        }
    )
    scoped = raw.merge(
        decisions[
            [
                "species",
                "taxonomic_group",
                "taxon_scope_review_status",
                "taxon_scope_decision_basis",
                "taxon_scope_source_citation",
                "taxon_scope_source_url",
                "taxon_scope_reviewer",
                "taxon_scope_decision_date",
            ]
        ],
        on="species",
        how="left",
        validate="many_to_one",
    )
    scoped["taxonomic_group"] = _text(scoped["taxonomic_group"]).replace("", "unresolved")
    scoped["taxon_scope_review_status"] = _text(scoped["taxon_scope_review_status"]).replace(
        "", "unreviewed"
    )
    for column in [
        "taxon_scope_decision_basis",
        "taxon_scope_source_citation",
        "taxon_scope_source_url",
        "taxon_scope_reviewer",
        "taxon_scope_decision_date",
    ]:
        scoped[column] = _text(scoped[column])
    scoped["angiosperm_analysis_eligible"] = (
        scoped["taxonomic_group"].eq("angiosperm")
        & scoped["taxon_scope_review_status"].str.lower().eq("accepted")
    )
    scoped = scoped[OUTPUT_COLUMNS].sort_values(["island_id", "species"]).reset_index(drop=True)
    return scoped, scoped.loc[scoped["angiosperm_analysis_eligible"]].copy()


def taxon_scope_summary(scoped: pd.DataFrame) -> dict[str, Any]:
    """Summarise how much of each island's raw species universe is review-eligible."""
    by_group = scoped["taxonomic_group"].value_counts().sort_index().to_dict()
    reviewed = scoped["taxon_scope_review_status"].str.lower().eq("accepted")
    eligible = scoped["angiosperm_analysis_eligible"]
    by_island = (
        scoped.groupby("island_id", sort=True)
        .agg(
            n_raw_species=("species", "nunique"),
            n_accepted_angiosperms=("angiosperm_analysis_eligible", "sum"),
        )
        .reset_index()
    )
    by_island["accepted_angiosperm_coverage_fraction"] = (
        by_island["n_accepted_angiosperms"] / by_island["n_raw_species"].where(by_island["n_raw_species"] > 0)
    ).fillna(0.0)
    return {
        "n_island_species_pairs": int(len(scoped)),
        "n_unique_raw_species": int(scoped["species"].nunique()),
        "n_accepted_scope_decisions": int(reviewed.sum()),
        "n_accepted_angiosperm_pairs": int(eligible.sum()),
        "n_islands_with_accepted_angiosperms": int(scoped.loc[eligible, "island_id"].nunique()),
        "taxonomic_group_counts": {str(key): int(value) for key, value in by_group.items()},
        "per_island": by_island.to_dict("records"),
    }


@app.command("build")
def build(
    island_species_csv: Path = typer.Option(..., exists=True, help="Raw island_species_occurrences.csv."),
    taxon_scope_decisions_csv: Path = typer.Option(
        ..., exists=True, help="Reviewed accepted_species taxonomic-scope table."
    ),
    output_dir: Path = typer.Option(..., help="Directory for scoped species products."),
) -> None:
    """Write full taxonomic-scope audit and accepted angiosperm species table."""
    island_species = pd.read_csv(island_species_csv, dtype=str).fillna("")
    decisions = pd.read_csv(taxon_scope_decisions_csv, dtype=str).fillna("")
    scoped, angiosperms = build_taxon_scope(island_species, decisions)
    summary = taxon_scope_summary(scoped)

    output_dir.mkdir(parents=True, exist_ok=True)
    scoped.to_csv(output_dir / "island_species_taxonomic_scope.csv", index=False)
    angiosperms.to_csv(output_dir / "island_angiosperm_species.csv", index=False)
    candidate_taxa = scoped[
        ["species", "taxonomic_group", "taxon_scope_review_status"]
    ].drop_duplicates("species").sort_values("species")
    candidate_taxa.to_csv(output_dir / "taxon_scope_review_queue.csv", index=False)
    (output_dir / "taxon_scope_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"Scoped {summary['n_unique_raw_species']} raw species; "
        f"{summary['n_accepted_angiosperm_pairs']} island-species pairs are accepted angiosperms."
    )


if __name__ == "__main__":
    app()
