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
COVERAGE_COLUMNS = [
    "island_id",
    "n_raw_species",
    "n_raw_records",
    "n_accepted_angiosperm_species",
    "n_accepted_angiosperm_records",
    "accepted_angiosperm_species_fraction",
    "accepted_angiosperm_record_fraction",
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


def _record_count(table: pd.DataFrame) -> pd.Series:
    """Use species-level record counts when available, otherwise one row = one record."""
    if "n_records" not in table.columns:
        return pd.Series(1, index=table.index, dtype="int64")
    values = pd.to_numeric(table["n_records"], errors="coerce").fillna(0)
    return values.clip(lower=0).astype(int)


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
    raw["_n_records"] = _record_count(raw)
    raw = raw.loc[raw["island_id"].ne("") & raw["species"].ne("")].copy()
    raw = raw.sort_values(["island_id", "species"]).drop_duplicates(["island_id", "species"], keep="first")

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
    scoped = scoped.sort_values(["island_id", "species"]).reset_index(drop=True)
    return scoped[OUTPUT_COLUMNS], scoped.loc[scoped["angiosperm_analysis_eligible"], OUTPUT_COLUMNS].copy()


def island_angiosperm_coverage(scoped: pd.DataFrame, island_species: pd.DataFrame) -> pd.DataFrame:
    """Summarise the reviewed angiosperm denominator for each island."""
    records = island_species.copy()
    records["island_id"] = _text(records["island_id"])
    records["species"] = _text(records["species"])
    records["_n_records"] = _record_count(records)
    records = records.loc[records["island_id"].ne("") & records["species"].ne("")].copy()
    records = records.sort_values(["island_id", "species"]).drop_duplicates(["island_id", "species"], keep="first")
    eligibility = scoped[["island_id", "species", "angiosperm_analysis_eligible"]].copy()
    joined = records.merge(eligibility, on=["island_id", "species"], how="left", validate="one_to_one")
    joined["angiosperm_analysis_eligible"] = joined["angiosperm_analysis_eligible"].fillna(False)

    rows: list[dict[str, Any]] = []
    for island_id, group in joined.groupby("island_id", sort=True):
        raw_species = int(group["species"].nunique())
        raw_records = int(group["_n_records"].sum())
        accepted = group.loc[group["angiosperm_analysis_eligible"]]
        accepted_species = int(accepted["species"].nunique())
        accepted_records = int(accepted["_n_records"].sum())
        rows.append(
            {
                "island_id": island_id,
                "n_raw_species": raw_species,
                "n_raw_records": raw_records,
                "n_accepted_angiosperm_species": accepted_species,
                "n_accepted_angiosperm_records": accepted_records,
                "accepted_angiosperm_species_fraction": accepted_species / raw_species if raw_species else 0.0,
                "accepted_angiosperm_record_fraction": accepted_records / raw_records if raw_records else 0.0,
            }
        )
    return pd.DataFrame(rows, columns=COVERAGE_COLUMNS)


def taxon_scope_summary(scoped: pd.DataFrame, coverage: pd.DataFrame) -> dict[str, Any]:
    """Summarise review status and angiosperm coverage without treating unknown taxa as absence."""
    by_group = scoped["taxonomic_group"].value_counts().sort_index().to_dict()
    reviewed = scoped["taxon_scope_review_status"].str.lower().eq("accepted")
    eligible = scoped["angiosperm_analysis_eligible"]
    return {
        "n_island_species_pairs": int(len(scoped)),
        "n_unique_raw_species": int(scoped["species"].nunique()),
        "n_accepted_scope_decisions": int(reviewed.sum()),
        "n_accepted_angiosperm_pairs": int(eligible.sum()),
        "n_islands_with_accepted_angiosperms": int(scoped.loc[eligible, "island_id"].nunique()),
        "taxonomic_group_counts": {str(key): int(value) for key, value in by_group.items()},
        "n_islands_in_coverage_table": int(len(coverage)),
    }


@app.command("build")
def build(
    island_species_csv: Path = typer.Option(..., exists=True, help="Raw island_species_occurrences.csv."),
    taxon_scope_decisions_csv: Path = typer.Option(
        ..., exists=True, help="Reviewed accepted_species taxonomic-scope table."
    ),
    output_dir: Path = typer.Option(..., help="Directory for scoped species products."),
) -> None:
    """Write taxonomic-scope audit, angiosperm table, and coverage-gate input."""
    island_species = pd.read_csv(island_species_csv, dtype=str).fillna("")
    decisions = pd.read_csv(taxon_scope_decisions_csv, dtype=str).fillna("")
    scoped, angiosperms = build_taxon_scope(island_species, decisions)
    coverage = island_angiosperm_coverage(scoped, island_species)
    summary = taxon_scope_summary(scoped, coverage)

    output_dir.mkdir(parents=True, exist_ok=True)
    scoped.to_csv(output_dir / "island_species_taxonomic_scope.csv", index=False)
    angiosperms.to_csv(output_dir / "island_angiosperm_species.csv", index=False)
    coverage.to_csv(output_dir / "island_angiosperm_coverage.csv", index=False)
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
