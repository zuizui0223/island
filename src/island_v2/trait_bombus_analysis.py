"""Aggregate island trait composition and join Bombus island metrics."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_TRAIT_COLUMNS = {
    "island_id",
    "accepted_species",
    "trait_name",
    "filled_value",
}
REQUIRED_BOMBUS_COLUMNS = {
    "island_id",
    "n_candidate_bombus_species",
    "n_scored_bombus_species",
    "n_unscored_bombus_species",
    "n_environmentally_compatible_species",
    "environmentally_compatible_species_share",
    "environmental_compatibility_max",
    "environmental_compatibility_mean",
    "environmental_extrapolation_species_share",
    "bombus_regime",
    "primary_bombus_channel_analysis",
    "global_comparison_only",
}


def aggregate_trait_composition(table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    missing = REQUIRED_TRAIT_COLUMNS.difference(table.columns)
    if missing:
        raise typer.BadParameter(f"trait table missing columns: {sorted(missing)}")

    work = table.copy()
    for column in REQUIRED_TRAIT_COLUMNS:
        work[column] = work[column].fillna("").astype(str).str.strip()
    work = work.loc[
        work["island_id"].ne("")
        & work["accepted_species"].ne("")
        & work["trait_name"].ne("")
        & work["filled_value"].ne("")
    ].copy()
    work = work.drop_duplicates(["island_id", "accepted_species", "trait_name"])

    counts = (
        work.groupby(["island_id", "trait_name", "filled_value"], sort=True)
        .agg(n_species=("accepted_species", "nunique"))
        .reset_index()
    )
    totals = counts.groupby(["island_id", "trait_name"])["n_species"].transform("sum")
    counts["species_share_within_trait"] = counts["n_species"] / totals

    richness = (
        work.groupby("island_id", sort=True)["accepted_species"]
        .nunique()
        .rename("n_trait_species")
        .reset_index()
    )
    return counts, richness


def build_analysis_input(
    trait_table: pd.DataFrame,
    bombus_table: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    missing = REQUIRED_BOMBUS_COLUMNS.difference(bombus_table.columns)
    if missing:
        raise typer.BadParameter(f"Bombus table missing columns: {sorted(missing)}")

    composition, richness = aggregate_trait_composition(trait_table)
    bombus = bombus_table.drop_duplicates("island_id").copy()
    island_input = bombus.merge(richness, on="island_id", how="left", validate="one_to_one")
    island_input["n_trait_species"] = island_input["n_trait_species"].fillna(0).astype(int)

    joined_long = composition.merge(
        island_input,
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    return island_input, joined_long


@app.command("run")
def run(
    island_species_traits_csv: Path = typer.Option(..., exists=True),
    bombus_island_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    traits = pd.read_csv(island_species_traits_csv, dtype=str).fillna("")
    bombus = pd.read_csv(bombus_island_csv)
    island_input, joined_long = build_analysis_input(traits, bombus)

    output_dir.mkdir(parents=True, exist_ok=True)
    island_input.to_csv(output_dir / "island_bombus_trait_analysis_input.csv", index=False)
    joined_long.to_csv(
        output_dir / "island_trait_composition_with_bombus.csv.gz",
        index=False,
        compression="gzip",
    )
    summary = {
        "contract": "island_bombus_trait_analysis_v1",
        "n_islands": int(island_input["island_id"].nunique()),
        "n_trait_composition_rows": int(len(joined_long)),
        "traits": sorted(joined_long["trait_name"].dropna().astype(str).unique().tolist()),
        "primary_bombus_channel_islands": int(
            island_input["primary_bombus_channel_analysis"].fillna(False).astype(bool).sum()
        ),
        "interpretation": (
            "Trait composition is kept separate by trait and value. Bombus metrics are island-level predictors; "
            "structural absence and unresolved applicability remain explicit."
        ),
    }
    (output_dir / "island_bombus_trait_analysis_manifest.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary))


if __name__ == "__main__":
    app()
