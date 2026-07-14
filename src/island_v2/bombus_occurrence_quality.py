"""Lightweight Bombus occurrence quality filtering and spatial thinning.

This module deliberately depends only on pandas so validation and unit tests do not
need rasterio/geopandas. Raster acquisition remains in ``bombus_niche_inputs``.
"""

from __future__ import annotations

import math

import pandas as pd

EXCLUDED_BASIS = {"FOSSIL_SPECIMEN", "MACHINE_OBSERVATION"}


def quality_filter_and_thin_occurrences(
    frame: pd.DataFrame,
    *,
    grid_degrees: float = 0.25,
    max_coordinate_uncertainty_m: float = 20_000,
    min_year: int = 1950,
) -> pd.DataFrame:
    """Filter low-quality records and retain one deterministic record per species×cell."""
    if grid_degrees <= 0:
        raise ValueError("grid_degrees must be positive")

    required = {"bombus_species", "decimal_latitude", "decimal_longitude"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"occurrence table missing required columns: {sorted(missing)}")

    d = frame.copy()
    d["decimal_latitude"] = pd.to_numeric(d["decimal_latitude"], errors="coerce")
    d["decimal_longitude"] = pd.to_numeric(d["decimal_longitude"], errors="coerce")
    d["year_numeric"] = pd.to_numeric(d.get("year"), errors="coerce")
    d["uncertainty_numeric"] = pd.to_numeric(
        d.get("coordinate_uncertainty_m"), errors="coerce"
    )

    basis = (
        d["basis_of_record"].fillna("")
        if "basis_of_record" in d.columns
        else pd.Series("", index=d.index)
    )
    d = d[
        d["decimal_latitude"].between(-90, 90)
        & d["decimal_longitude"].between(-180, 180)
    ]
    d = d[~basis.loc[d.index].isin(EXCLUDED_BASIS)]
    d = d[d["year_numeric"].isna() | (d["year_numeric"] >= min_year)]
    d = d[
        d["uncertainty_numeric"].isna()
        | (d["uncertainty_numeric"] <= max_coordinate_uncertainty_m)
    ]
    d = d.drop_duplicates(
        ["bombus_species", "decimal_latitude", "decimal_longitude"]
    )

    d["grid_x"] = (d["decimal_longitude"] / grid_degrees).map(math.floor)
    d["grid_y"] = (d["decimal_latitude"] / grid_degrees).map(math.floor)
    d = d.sort_values(
        ["bombus_species", "grid_x", "grid_y", "uncertainty_numeric", "year_numeric"],
        na_position="last",
    )
    return d.drop_duplicates(
        ["bombus_species", "grid_x", "grid_y"], keep="first"
    ).drop(columns=["grid_x", "grid_y", "year_numeric", "uncertainty_numeric"])
