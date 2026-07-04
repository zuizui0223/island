"""Nominate Core-pilot islands using only outcome-blind eligibility gates.

A Core-pilot island is one on which the Bombus confirmatory workflow may begin.
This nomination is deliberately narrow and outcome-blind: it uses ONLY

1. frozen Bombus applicability (from island-v2-bombus-applicability), which is
   itself decided only from source-region native Bombus status, source pool, and
   reviewed assignments; and
2. whether the current GBIF collection produced enough **raw exact-island
   species labels** for the island (a data-availability gate).

The second gate is deliberately a *raw data-availability* check, NOT a measure
of reviewed flora richness. The count is the number of raw exact-island GBIF
species labels collected so far; it has not passed accepted taxonomic scope,
native-establishment review, or an angiosperm-flora anchor. It only decides
whether there is enough exact-island data to *start* Core-pilot work.

This nomination never uses floral traits, current island Bombus records or
non-detection, reproductive traits, or model results. Passing nomination does
NOT accept any taxon, native status, or trait value, and puts nothing into an
analysis table; a nominated island still requires an island flora anchor and
taxon/establishment review first. When no island is eligible, the gap is
reported explicitly rather than inventing a Core island.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_APPLICABILITY_COLUMNS = {
    "island_id",
    "source_region_id",
    "applicability",
    "regional_analysis_only",
    "native_bombus_status",
}
REQUIRED_EFFORT_COLUMNS = {"island_id", "n_species"}
# Nomination must never be influenced by these; presence is an error.
PROHIBITED_COLUMNS = {
    "flower_primary_color",
    "floral_traits",
    "reproductive_traits",
    "island_Bombus_occurrence_or_non_detection",
    "bombus_occurrence_evidence",
    "model_results",
    "analysis_included",
}


@app.callback()
def main() -> None:
    """Nominate Core-pilot islands from applicability and flora-data availability."""


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _bool(series: pd.Series) -> pd.Series:
    return _text(series).str.lower().isin({"true", "1", "yes", "y"})


def _guard_inputs(table: pd.DataFrame, label: str) -> None:
    prohibited = PROHIBITED_COLUMNS.intersection(table.columns)
    if prohibited:
        raise typer.BadParameter(
            f"{label} contains prohibited outcome/predictor columns: {', '.join(sorted(prohibited))}"
        )


def nominate_core_pilots(
    applicability: pd.DataFrame,
    collected_effort: pd.DataFrame,
    min_raw_exact_island_species_labels: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Return per-island nomination rows and an eligibility/gap report."""
    missing_a = REQUIRED_APPLICABILITY_COLUMNS.difference(applicability.columns)
    if missing_a:
        raise typer.BadParameter(f"applicability table missing columns: {', '.join(sorted(missing_a))}")
    missing_e = REQUIRED_EFFORT_COLUMNS.difference(collected_effort.columns)
    if missing_e:
        raise typer.BadParameter(f"collected effort table missing columns: {', '.join(sorted(missing_e))}")
    _guard_inputs(applicability, "applicability table")
    _guard_inputs(collected_effort, "collected effort table")

    table = applicability.copy()
    table["island_id"] = _text(table["island_id"])
    table["applicability"] = _text(table["applicability"])
    table["is_core_applicable"] = table["applicability"].eq("applicable") & ~_bool(
        table["regional_analysis_only"]
    )

    effort = collected_effort.copy()
    effort["island_id"] = _text(effort["island_id"])
    effort["n_collected_exact_island_raw_species_labels"] = (
        pd.to_numeric(effort["n_species"], errors="coerce").fillna(0).astype(int)
    )
    raw_labels = effort.groupby("island_id")["n_collected_exact_island_raw_species_labels"].max()

    min_labels = int(min_raw_exact_island_species_labels)
    table["n_collected_exact_island_raw_species_labels"] = (
        table["island_id"].map(raw_labels).fillna(0).astype(int)
    )
    table["has_raw_exact_island_taxon_availability"] = (
        table["n_collected_exact_island_raw_species_labels"] >= min_labels
    )
    table["core_pilot_eligible"] = (
        table["is_core_applicable"] & table["has_raw_exact_island_taxon_availability"]
    )
    table["nomination_note"] = [
        (
            "Eligible to START Core-pilot: Bombus-applicable Core island with enough raw "
            "exact-island species labels. This is a data-availability gate, not analysis "
            "inclusion; still requires an island flora anchor and taxon/establishment review."
            if row.core_pilot_eligible
            else "Applicable Core island but not enough raw exact-island species labels yet "
            "(collect its GBIF block)."
            if row.is_core_applicable
            else "Not a Bombus-applicable Core island."
        )
        for row in table.itertuples(index=False)
    ]

    eligible = table.loc[table["core_pilot_eligible"]]
    applicable_no_data = table.loc[
        table["is_core_applicable"] & ~table["has_raw_exact_island_taxon_availability"]
    ]
    report = {
        "gate_meaning": (
            "Raw exact-island GBIF species-label availability to START Core-pilot work; "
            "not reviewed flora richness, native establishment, or analysis inclusion."
        ),
        "min_raw_exact_island_species_labels": min_labels,
        "n_islands_considered": int(len(table)),
        "n_core_applicable_islands": int(table["is_core_applicable"].sum()),
        "n_core_applicable_with_raw_species_labels": int(
            (table["is_core_applicable"] & table["has_raw_exact_island_taxon_availability"]).sum()
        ),
        "n_core_pilot_eligible": int(len(eligible)),
        "eligible_island_ids": eligible["island_id"].tolist(),
        "applicable_but_no_raw_species_labels": applicable_no_data[
            ["island_id", "source_region_id", "n_collected_exact_island_raw_species_labels"]
        ].to_dict("records"),
    }
    return table.reset_index(drop=True), report


@app.command("nominate")
def nominate(
    applicability_csv: Path = typer.Option(..., exists=True, help="bombus_applicability_by_island.csv"),
    collected_effort_csv: Path = typer.Option(..., exists=True, help="island_observation_effort.csv"),
    output_dir: Path = typer.Option(...),
    min_raw_exact_island_species_labels: int = typer.Option(
        5, min=1, help="Min raw exact-island GBIF species labels available to START Core-pilot work."
    ),
) -> None:
    """Write Core-pilot nominations and an explicit eligibility/gap report."""
    applicability = pd.read_csv(applicability_csv, dtype=str).fillna("")
    effort = pd.read_csv(collected_effort_csv, dtype=str).fillna("")
    table, report = nominate_core_pilots(applicability, effort, min_raw_exact_island_species_labels)
    output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_dir / "core_pilot_nomination_table.csv", index=False)
    table.loc[table["core_pilot_eligible"]].to_csv(output_dir / "core_pilot_candidates.csv", index=False)
    (output_dir / "core_pilot_nomination_report.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    if report["n_core_pilot_eligible"]:
        typer.echo(
            f"{report['n_core_pilot_eligible']} Core-pilot candidate(s): {report['eligible_island_ids']}"
        )
    else:
        typer.echo(
            "No eligible Core-pilot island in the collected blocks. "
            f"{report['n_core_applicable_islands']} Core-applicable island(s) present; "
            f"{len(report['applicable_but_no_raw_species_labels'])} await exact-island collection."
        )


if __name__ == "__main__":
    app()
