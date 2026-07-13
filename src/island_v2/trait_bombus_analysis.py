"""Aggregate island trait composition and join canonical Bombus island metrics."""

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

BOMBUS_COLUMN_ALIASES = {
    "n_candidate_bombus_species": "n_bombus_species_total",
    "n_scored_bombus_species": "n_bombus_species_scored",
    "n_unscored_bombus_species": "n_bombus_species_unscored",
    "environmentally_compatible_species_share": "share_environmentally_compatible_species",
    "environmental_extrapolation_species_share": "environmental_extrapolation_share",
    "bombus_regime": "bombus_regime_set",
}

REQUIRED_BOMBUS_COLUMNS = {
    "island_id",
    "n_bombus_species_total",
    "n_bombus_species_scored",
    "n_bombus_species_unscored",
    "n_environmentally_compatible_species",
    "share_environmentally_compatible_species",
    "environmental_compatibility_max",
    "environmental_compatibility_mean",
    "environmental_extrapolation_share",
    "bombus_regime_set",
    "primary_bombus_channel_analysis",
    "global_comparison_only",
}


def normalize_bombus_columns(table: pd.DataFrame) -> pd.DataFrame:
    """Accept the legacy PR #113 names and the canonical PR #112 names."""
    result = table.copy()
    for legacy, canonical in BOMBUS_COLUMN_ALIASES.items():
        if canonical not in result.columns and legacy in result.columns:
            result = result.rename(columns={legacy: canonical})
    missing = REQUIRED_BOMBUS_COLUMNS.difference(result.columns)
    if missing:
        raise typer.BadParameter(f"Bombus table missing columns: {sorted(missing)}")
    if result["island_id"].duplicated().any():
        raise typer.BadParameter("Bombus island table must contain one row per island_id")
    return result


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
    duplicate_keys = ["island_id", "accepted_species", "trait_name"]
    if "analysis_tier" in work.columns:
        duplicate_keys.append("analysis_tier")
    work = work.drop_duplicates(duplicate_keys)

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
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    bombus = normalize_bombus_columns(bombus_table)
    composition, richness = aggregate_trait_composition(trait_table)
    island_input = bombus.merge(richness, on="island_id", how="left", validate="one_to_one")
    island_input["n_trait_species"] = island_input["n_trait_species"].fillna(0).astype(int)

    joined_long = composition.merge(
        island_input,
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    species_master = trait_table.merge(
        bombus,
        on="island_id",
        how="left",
        validate="many_to_one",
        indicator="bombus_join_status",
    )
    return island_input, joined_long, species_master


@app.command("run")
def run(
    island_species_traits_csv: Path = typer.Option(..., exists=True),
    bombus_island_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    traits = pd.read_csv(island_species_traits_csv, dtype=str).fillna("")
    bombus = pd.read_csv(bombus_island_csv)
    island_input, joined_long, species_master = build_analysis_input(traits, bombus)

    output_dir.mkdir(parents=True, exist_ok=True)
    island_input.to_csv(output_dir / "island_bombus_trait_analysis_input.csv", index=False)
    joined_long.to_csv(
        output_dir / "island_trait_composition_with_bombus.csv.gz",
        index=False,
        compression="gzip",
    )
    species_master.to_csv(
        output_dir / "island_species_trait_bombus_master.csv.gz",
        index=False,
        compression="gzip",
    )

    trait_islands = set(traits["island_id"].astype(str))
    bombus_islands = set(bombus["island_id"].astype(str))
    unmatched_trait_islands = sorted(value for value in trait_islands - bombus_islands if value)
    unmatched_bombus_islands = sorted(value for value in bombus_islands - trait_islands if value)
    pd.DataFrame({"island_id": unmatched_trait_islands}).to_csv(
        output_dir / "trait_islands_missing_bombus_metrics.csv", index=False
    )
    pd.DataFrame({"island_id": unmatched_bombus_islands}).to_csv(
        output_dir / "bombus_islands_missing_trait_occurrences.csv", index=False
    )

    join_counts = species_master["bombus_join_status"].value_counts(dropna=False).to_dict()
    summary = {
        "contract": "island_bombus_trait_analysis_v2",
        "n_islands": int(island_input["island_id"].nunique()),
        "n_species_trait_master_rows": int(len(species_master)),
        "n_trait_composition_rows": int(len(joined_long)),
        "traits": sorted(joined_long["trait_name"].dropna().astype(str).unique().tolist()),
        "primary_bombus_channel_islands": int(
            island_input["primary_bombus_channel_analysis"].fillna(False).astype(bool).sum()
        ),
        "n_trait_islands_missing_bombus_metrics": int(len(unmatched_trait_islands)),
        "n_bombus_islands_missing_trait_occurrences": int(len(unmatched_bombus_islands)),
        "bombus_join_status_counts": {str(key): int(value) for key, value in join_counts.items()},
        "interpretation": (
            "The species-level master keeps trait fill tier and Bombus island predictors together. "
            "Applicability, structural absence, unresolved status, and extrapolation quality remain explicit."
        ),
    }
    (output_dir / "island_bombus_trait_analysis_manifest.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary))


if __name__ == "__main__":
    app()
