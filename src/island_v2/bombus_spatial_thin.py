"""Deterministically spatially thin Bombus occurrence records before niche fitting."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)


def spatial_thin(
    frame: pd.DataFrame,
    *,
    cell_size_degrees: float,
    species_column: str = "bombus_species",
    latitude_column: str = "decimal_latitude",
    longitude_column: str = "decimal_longitude",
) -> pd.DataFrame:
    """Keep one deterministic record per species and geographic grid cell."""
    if cell_size_degrees <= 0:
        raise typer.BadParameter("cell_size_degrees must be positive")
    required = {species_column, latitude_column, longitude_column}
    missing = required.difference(frame.columns)
    if missing:
        raise typer.BadParameter(f"occurrence table missing columns: {sorted(missing)}")

    result = frame.copy()
    result[latitude_column] = pd.to_numeric(result[latitude_column], errors="coerce")
    result[longitude_column] = pd.to_numeric(result[longitude_column], errors="coerce")
    result = result.dropna(subset=[species_column, latitude_column, longitude_column])
    result = result.loc[
        result[latitude_column].between(-90, 90)
        & result[longitude_column].between(-180, 180)
    ].copy()

    result["spatial_thin_lat_cell"] = (
        ((result[latitude_column] + 90.0) / cell_size_degrees).astype(int)
    )
    result["spatial_thin_lon_cell"] = (
        ((result[longitude_column] + 180.0) / cell_size_degrees).astype(int)
    )

    sort_columns = [species_column, "spatial_thin_lat_cell", "spatial_thin_lon_cell"]
    for optional in ("year", "gbif_id"):
        if optional in result.columns:
            sort_columns.append(optional)
    result = result.sort_values(sort_columns, kind="mergesort")
    result = result.drop_duplicates(
        [species_column, "spatial_thin_lat_cell", "spatial_thin_lon_cell"],
        keep="first",
    )
    return result.reset_index(drop=True)


@app.command("run")
def run(
    input_csv: Path = typer.Option(..., exists=True),
    output_csv: Path = typer.Option(...),
    manifest_json: Path = typer.Option(...),
    cell_size_degrees: float = typer.Option(0.25, min=0.01, max=5.0),
) -> None:
    """Thin a Bombus occurrence-environment CSV and record retention diagnostics."""
    source = pd.read_csv(input_csv)
    thinned = spatial_thin(source, cell_size_degrees=cell_size_degrees)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    thinned.to_csv(output_csv, index=False)

    before_by_species = source.groupby("bombus_species").size()
    after_by_species = thinned.groupby("bombus_species").size()
    diagnostics = (
        pd.DataFrame({"n_before": before_by_species, "n_after": after_by_species})
        .fillna(0)
        .astype(int)
    )
    diagnostics["retention_fraction"] = diagnostics["n_after"] / diagnostics["n_before"]
    diagnostics.reset_index().to_csv(
        output_csv.with_name("bombus_spatial_thinning_by_species.csv"), index=False
    )

    manifest = {
        "contract": "bombus_spatial_thinning_v1",
        "cell_size_degrees": cell_size_degrees,
        "n_rows_before": int(len(source)),
        "n_rows_after": int(len(thinned)),
        "retention_fraction": float(len(thinned) / len(source)) if len(source) else None,
        "n_species_before": int(source["bombus_species"].nunique()),
        "n_species_after": int(thinned["bombus_species"].nunique()),
        "selection_rule": "one deterministic record per species x latitude-cell x longitude-cell",
    }
    manifest_json.parent.mkdir(parents=True, exist_ok=True)
    manifest_json.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
