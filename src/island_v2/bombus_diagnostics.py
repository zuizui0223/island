"""Effort-aware island-level Bombus occurrence diagnostics.

This module implements the observation-process half of the v2 Bombus channel.
It deliberately separates a detected Bombus record from interpretable
non-detection. A missing GBIF record is never coded as biological absence.

Input records must already have been assigned to exact island polygons (or to a
separately predeclared island-system spatial rule). The module is intentionally
agnostic about download method so the same diagnostics can be applied to GBIF
exports, curated occurrence compilations, or a future dedicated pollinator
collector.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_RECORD_COLUMNS = {
    "island_id",
    "family",
    "genus",
    "dataset_key",
    "year",
    "decimal_latitude",
    "decimal_longitude",
}

OUTPUT_COLUMNS = [
    "island_id",
    "target_group_definition",
    "bombus_occurrence_evidence",
    "bombus_channel_evaluable",
    "bombus_record_count",
    "target_group_record_count",
    "target_group_spatial_units",
    "target_group_temporal_units",
    "distinct_dataset_count",
    "proportion_preserved_specimens",
    "latest_record_year",
    "n_target_records_before_quality_filter",
    "n_target_records_excluded_unresolved_taxonomy",
    "n_target_records_excluded_invalid_coordinates",
    "n_target_records_excluded_establishment",
    "observation_diagnostic_flags",
    "occupancy_sensitivity_eligible",
]


@app.callback()
def main() -> None:
    """Classify island-level Bombus occurrence evidence without treating non-records as absence."""


def load_config(path: Path) -> dict[str, Any]:
    """Load and minimally validate the versioned Bombus observation policy."""
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise typer.BadParameter("Bombus diagnostics config must be a mapping")
    if "primary_target_group" not in config or "primary_thresholds" not in config:
        raise typer.BadParameter(
            "Bombus diagnostics config must contain primary_target_group and primary_thresholds"
        )
    thresholds = config["primary_thresholds"]
    required = {
        "min_background_records",
        "min_background_spatial_units",
        "min_background_temporal_units",
        "min_distinct_datasets",
        "required_bombus_records_for_detection",
    }
    missing = required.difference(thresholds)
    if missing:
        raise typer.BadParameter(f"Bombus diagnostics config missing thresholds: {sorted(missing)}")
    return config


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _bool_contains(series: pd.Series, token: str) -> pd.Series:
    return _text(series).str.upper().str.contains(token, regex=False)


def _normalise_records(records: pd.DataFrame, grid_degrees: float) -> pd.DataFrame:
    """Normalise a pre-assigned occurrence table and derive audit fields."""
    missing = REQUIRED_RECORD_COLUMNS.difference(records.columns)
    if missing:
        raise typer.BadParameter(
            f"Bombus occurrence table missing columns: {', '.join(sorted(missing))}"
        )
    if grid_degrees <= 0:
        raise typer.BadParameter("spatial_grid_degrees must be greater than zero")

    out = records.copy()
    out["_island_id"] = _text(out["island_id"])
    out["_family"] = _text(out["family"]).str.upper()
    out["_genus"] = _text(out["genus"]).str.upper()
    out["_dataset"] = _text(out["dataset_key"])
    out["_year"] = pd.to_numeric(out["year"], errors="coerce")
    out["_lat"] = pd.to_numeric(out["decimal_latitude"], errors="coerce")
    out["_lon"] = pd.to_numeric(out["decimal_longitude"], errors="coerce")
    out["_valid_coordinates"] = (
        out["_lat"].between(-90, 90, inclusive="both")
        & out["_lon"].between(-180, 180, inclusive="both")
    )
    out["_valid_year"] = out["_year"].between(1500, 2100, inclusive="both")
    out["_has_resolved_genus"] = out["_genus"].ne("") & ~out["_genus"].isin(
        {"UNRESOLVED", "UNKNOWN", "INCERTAE SEDIS"}
    )

    establishment = out.get("establishment_means", pd.Series("", index=out.index))
    out["_excluded_establishment"] = _bool_contains(establishment, "CULTIVATED") | _bool_contains(
        establishment, "CAPTIVE"
    )

    # GBIF downloads may overlap. When a stable GBIF identifier exists, retain
    # one deterministic copy so background effort is not inflated by duplicate
    # request blocks.
    if "gbif_id" in out.columns:
        gbif_id = _text(out["gbif_id"])
        has_id = gbif_id.ne("")
        out["_gbif_id"] = gbif_id
        out["_dedupe_keep"] = ~has_id | ~out.loc[has_id].duplicated("_gbif_id", keep="first")
    else:
        out["_dedupe_keep"] = True

    out["_spatial_unit"] = [
        f"{math.floor(lon / grid_degrees)}:{math.floor(lat / grid_degrees)}"
        if pd.notna(lon) and pd.notna(lat)
        else ""
        for lon, lat in zip(out["_lon"], out["_lat"], strict=True)
    ]
    return out


def _quality_mask(records: pd.DataFrame, target_group: str) -> pd.Series:
    """Records eligible to quantify target-group background observation effort."""
    return (
        records["_island_id"].ne("")
        & records["_family"].eq(target_group.upper())
        & records["_has_resolved_genus"]
        & records["_valid_coordinates"]
        & records["_valid_year"]
        & ~records["_excluded_establishment"]
        & records["_dedupe_keep"]
    )


def compute_bombus_diagnostics(records: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Return per-island Bombus detection/non-detection evidence and audit metrics.

    The `adequate_non_detection` category requires zero quality-filtered Bombus
    records *and* all background, spatial, temporal, and source-dispersion
    thresholds. `detected` requires a quality-filtered Bombus record but is not
    itself a quantitative abundance estimate.
    """
    target_group = str(config["primary_target_group"]["name"]).strip().upper()
    thresholds = config["primary_thresholds"]
    grid_degrees = float(config.get("spatial_grid_degrees", 0.1))
    work = _normalise_records(records, grid_degrees)
    if work.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)

    quality = _quality_mask(work, target_group)
    raw_island_ids = sorted(value for value in work["_island_id"].unique() if value)
    rows: list[dict[str, Any]] = []
    for island_id in raw_island_ids:
        group = work.loc[work["_island_id"].eq(island_id)].copy()
        raw_target = group["_family"].eq(target_group)
        eligible_target = group.loc[quality.loc[group.index]]
        bombus = eligible_target.loc[eligible_target["_genus"].eq("BOMBUS")]

        background_count = int(len(eligible_target))
        bombus_count = int(len(bombus))
        spatial_units = int(eligible_target["_spatial_unit"].replace("", pd.NA).dropna().nunique())
        temporal_units = int(eligible_target["_year"].dropna().astype(int).nunique())
        dataset_count = int(eligible_target["_dataset"].loc[eligible_target["_dataset"].ne("")].nunique())
        latest_year = eligible_target["_year"].dropna()
        latest_record_year: int | None = int(latest_year.max()) if not latest_year.empty else None

        basis = _text(eligible_target.get("basis_of_record", pd.Series("", index=eligible_target.index))).str.upper()
        specimen_fraction = float(basis.eq("PRESERVED_SPECIMEN").mean()) if background_count else 0.0

        unresolved_taxonomy = int((raw_target & ~group["_has_resolved_genus"]).sum())
        invalid_coordinates = int((raw_target & ~group["_valid_coordinates"]).sum())
        excluded_establishment = int((raw_target & group["_excluded_establishment"]).sum())
        raw_target_count = int(raw_target.sum())

        has_detection = bombus_count >= int(thresholds["required_bombus_records_for_detection"])
        adequate_background = background_count >= int(thresholds["min_background_records"])
        adequate_spatial = spatial_units >= int(thresholds["min_background_spatial_units"])
        adequate_temporal = temporal_units >= int(thresholds["min_background_temporal_units"])
        adequate_datasets = dataset_count >= int(thresholds["min_distinct_datasets"])
        adequate_effort = adequate_background and adequate_spatial and adequate_temporal and adequate_datasets

        # Reserve unresolved for a genuinely unusable target-group record set;
        # ordinary sparse coverage remains insufficient effort.
        if has_detection:
            evidence = "detected"
        elif raw_target_count > 0 and background_count == 0 and unresolved_taxonomy == raw_target_count:
            evidence = "unresolved"
        elif adequate_effort:
            evidence = "adequate_non_detection"
        else:
            evidence = "insufficient_effort"

        flags: list[str] = []
        if raw_target_count == 0:
            flags.append("no_target_group_records")
        if background_count == 0:
            flags.append("no_quality_filtered_target_records")
        if unresolved_taxonomy:
            flags.append("unresolved_target_taxonomy")
        if invalid_coordinates:
            flags.append("invalid_target_coordinates")
        if excluded_establishment:
            flags.append("excluded_target_establishment")
        if not adequate_background:
            flags.append("below_background_records")
        if not adequate_spatial:
            flags.append("low_spatial_dispersion")
        if not adequate_temporal:
            flags.append("low_temporal_dispersion")
        if not adequate_datasets:
            flags.append("low_dataset_dispersion")

        rows.append(
            {
                "island_id": island_id,
                "target_group_definition": target_group,
                "bombus_occurrence_evidence": evidence,
                "bombus_channel_evaluable": evidence in {"detected", "adequate_non_detection"},
                "bombus_record_count": bombus_count,
                "target_group_record_count": background_count,
                "target_group_spatial_units": spatial_units,
                "target_group_temporal_units": temporal_units,
                "distinct_dataset_count": dataset_count,
                "proportion_preserved_specimens": specimen_fraction,
                "latest_record_year": latest_record_year,
                "n_target_records_before_quality_filter": raw_target_count,
                "n_target_records_excluded_unresolved_taxonomy": unresolved_taxonomy,
                "n_target_records_excluded_invalid_coordinates": invalid_coordinates,
                "n_target_records_excluded_establishment": excluded_establishment,
                "observation_diagnostic_flags": "|".join(flags),
                "occupancy_sensitivity_eligible": adequate_effort,
            }
        )
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


def diagnostic_summary(diagnostics: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    """Return a compact JSON-safe audit of the evidence categories and flags."""
    evidence_counts = (
        diagnostics["bombus_occurrence_evidence"].value_counts().sort_index().to_dict()
        if not diagnostics.empty
        else {}
    )
    flag_counts: dict[str, int] = {}
    for value in diagnostics.get("observation_diagnostic_flags", pd.Series(dtype=str)):
        for flag in str(value).split("|"):
            if flag:
                flag_counts[flag] = flag_counts.get(flag, 0) + 1
    return {
        "config_version": config.get("version"),
        "target_group": config["primary_target_group"]["name"],
        "n_islands": int(len(diagnostics)),
        "n_bombus_channel_evaluable": int(diagnostics.get("bombus_channel_evaluable", pd.Series(dtype=bool)).sum()),
        "occurrence_evidence_counts": evidence_counts,
        "observation_diagnostic_flag_counts": flag_counts,
    }


@app.command("run")
def run(
    occurrences_csv: Path = typer.Option(..., exists=True, help="Pre-assigned island pollinator occurrence table."),
    output_dir: Path = typer.Option(..., help="Directory for Bombus diagnostics outputs."),
    config_path: Path = typer.Option(Path("config/bombus_observation_diagnostics.yml")),
) -> None:
    """Write effort-aware island-level Bombus observation diagnostics."""
    config = load_config(config_path)
    occurrences = pd.read_csv(occurrences_csv, dtype=str).fillna("")
    diagnostics = compute_bombus_diagnostics(occurrences, config)
    summary = diagnostic_summary(diagnostics, config)

    output_dir.mkdir(parents=True, exist_ok=True)
    diagnostics.to_csv(output_dir / "island_bombus_observation_diagnostics.csv", index=False)
    (output_dir / "bombus_observation_diagnostics_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    typer.echo(
        "Diagnosed "
        f"{summary['n_islands']} islands; {summary['n_bombus_channel_evaluable']} have detected "
        "Bombus or adequate non-detection evidence."
    )


if __name__ == "__main__":
    app()
