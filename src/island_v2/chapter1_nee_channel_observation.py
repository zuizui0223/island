"""Effort-aware island observation states for frozen N1 pollination channels.

Input occurrence rows must already come from a predeclared broad background-group
campaign for one channel and be assigned to exact island polygons.  The target set
is the independently frozen functional channel taxon catalog.  A missing target
record is never an absence unless the broad-group effort gate passes.
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

REQUIRED_RECORDS = {
    "island_id",
    "species",
    "dataset_key",
    "year",
    "decimal_latitude",
    "decimal_longitude",
}
REQUIRED_CATALOG = {"channel_id", "pollinator_species", "catalog_tier"}
OUTPUT_COLUMNS = [
    "island_id",
    "channel_id",
    "observation_state",
    "channel_record_count",
    "background_record_count",
    "background_spatial_units",
    "background_temporal_units",
    "distinct_dataset_count",
    "latest_background_year",
    "evidence_source",
    "quality_flags",
]


def load_policy(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("channel observation policy must be a mapping")
    if payload.get("contract") != "chapter1_nee_channel_observation_policy_v1":
        raise typer.BadParameter("unexpected channel observation policy contract")
    return payload


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _introduced(series: pd.Series) -> pd.Series:
    text = _text(series).str.upper()
    tokens = ("INTRODUCED", "INVASIVE", "NATURALISED", "NATURALIZED", "ALIEN")
    result = pd.Series(False, index=text.index)
    for token in tokens:
        result |= text.str.contains(token, regex=False)
    return result


def _captive(series: pd.Series) -> pd.Series:
    text = _text(series).str.upper()
    return text.str.contains("CAPTIVE", regex=False) | text.str.contains("CULTIVATED", regex=False)


def _prepare_records(records: pd.DataFrame) -> pd.DataFrame:
    missing = REQUIRED_RECORDS.difference(records.columns)
    if missing:
        raise ValueError(f"channel occurrence records missing columns: {sorted(missing)}")
    out = records.copy()
    out["_island"] = _text(out["island_id"])
    out["_species"] = _text(out["species"])
    out["_dataset"] = _text(out["dataset_key"])
    out["_year"] = pd.to_numeric(out["year"], errors="coerce")
    out["_lat"] = pd.to_numeric(out["decimal_latitude"], errors="coerce")
    out["_lon"] = pd.to_numeric(out["decimal_longitude"], errors="coerce")
    out["_valid_coords"] = out["_lat"].between(-90, 90) & out["_lon"].between(-180, 180)
    out["_valid_year"] = out["_year"].between(1500, 2100)
    establishment = out.get("establishment_means", pd.Series("", index=out.index))
    out["_introduced"] = _introduced(establishment)
    out["_captive"] = _captive(establishment)
    if "gbif_id" in out.columns:
        ids = _text(out["gbif_id"])
        has_id = ids.ne("")
        out["_dedupe_keep"] = ~has_id | ~ids.loc[has_id].duplicated(keep="first")
    else:
        out["_dedupe_keep"] = True
    out["_background_eligible"] = (
        out["_island"].ne("")
        & out["_species"].ne("")
        & out["_valid_coords"]
        & out["_valid_year"]
        & ~out["_captive"]
        & out["_dedupe_keep"]
    )
    out["_spatial_unit"] = [
        f"{math.floor(lon / 0.1)}:{math.floor(lat / 0.1)}"
        if pd.notna(lon) and pd.notna(lat)
        else ""
        for lon, lat in zip(out["_lon"], out["_lat"], strict=True)
    ]
    return out


def target_species(catalog: pd.DataFrame, channel_id: str, confirmatory_only: bool = True) -> set[str]:
    missing = REQUIRED_CATALOG.difference(catalog.columns)
    if missing:
        raise ValueError(f"channel taxon catalog missing columns: {sorted(missing)}")
    work = catalog.loc[catalog["channel_id"].astype(str).eq(str(channel_id))].copy()
    if confirmatory_only:
        work = work.loc[work["catalog_tier"].astype(str).eq("confirmatory")]
    return {value for value in _text(work["pollinator_species"]) if value}


def classify_channel_observations(
    records: pd.DataFrame,
    catalog: pd.DataFrame,
    channel_id: str,
    policy: dict[str, Any],
    *,
    effort_tier: str = "primary",
    confirmatory_catalog_only: bool = True,
    evidence_source: str = "GBIF_exact_island_background_campaign",
) -> pd.DataFrame:
    """Classify every represented island for one frozen background campaign."""
    if channel_id not in policy["channels"]:
        raise ValueError(f"unregistered N1 channel: {channel_id}")
    if effort_tier == "primary":
        thresholds = policy["primary_island_effort_gate"]
    elif effort_tier in {"liberal", "strict"}:
        thresholds = policy["sensitivity_effort_gates"][effort_tier]
    else:
        raise ValueError("effort_tier must be primary, liberal or strict")

    targets = target_species(catalog, channel_id, confirmatory_catalog_only)
    if not targets:
        raise ValueError(f"no target catalog taxa for channel {channel_id}")
    work = _prepare_records(records)
    if work.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)

    rows: list[dict[str, Any]] = []
    reference_year = int(policy["primary_island_effort_gate"]["reference_year"])
    max_age = int(thresholds["max_years_since_latest_background"])
    for island_id, group in work.groupby("_island", sort=True):
        background = group.loc[group["_background_eligible"]].copy()
        # Primary target detections exclude known introduced records. Unknown establishment
        # is retained but flagged, matching the prospective policy.
        target_mask = background["_species"].isin(targets) & ~background["_introduced"]
        target = background.loc[target_mask]
        n_target = int(len(target))
        n_background = int(len(background))
        spatial = int(background["_spatial_unit"].replace("", pd.NA).dropna().nunique())
        temporal = int(background["_year"].dropna().astype(int).nunique())
        datasets = int(background["_dataset"].loc[background["_dataset"].ne("")].nunique())
        latest = background["_year"].dropna()
        latest_year = int(latest.max()) if not latest.empty else None
        recent = latest_year is not None and 0 <= reference_year - latest_year <= max_age
        adequate = (
            n_background >= int(thresholds["min_background_records"])
            and spatial >= int(thresholds["min_background_spatial_units"])
            and temporal >= int(thresholds["min_background_temporal_units"])
            and datasets >= int(thresholds["min_distinct_datasets"])
            and recent
        )

        if n_target > 0:
            state = "detected"
        elif adequate:
            state = "adequate_non_detection"
        else:
            state = "insufficient_effort"

        flags: list[str] = []
        known_introduced_targets = int(
            (background["_species"].isin(targets) & background["_introduced"]).sum()
        )
        unknown_establishment_targets = int(
            (
                background["_species"].isin(targets)
                & _text(background.get("establishment_means", pd.Series("", index=background.index))).eq("")
            ).sum()
        )
        if known_introduced_targets:
            flags.append(f"introduced_target_rows_excluded={known_introduced_targets}")
        if unknown_establishment_targets:
            flags.append(f"target_rows_unknown_establishment={unknown_establishment_targets}")
        if n_background < int(thresholds["min_background_records"]):
            flags.append("below_background_records")
        if spatial < int(thresholds["min_background_spatial_units"]):
            flags.append("low_spatial_dispersion")
        if temporal < int(thresholds["min_background_temporal_units"]):
            flags.append("low_temporal_dispersion")
        if datasets < int(thresholds["min_distinct_datasets"]):
            flags.append("low_dataset_dispersion")
        if not recent:
            flags.append("stale_or_missing_background_recency")

        rows.append(
            {
                "island_id": island_id,
                "channel_id": channel_id,
                "observation_state": state,
                "channel_record_count": n_target,
                "background_record_count": n_background,
                "background_spatial_units": spatial,
                "background_temporal_units": temporal,
                "distinct_dataset_count": datasets,
                "latest_background_year": latest_year,
                "evidence_source": evidence_source,
                "quality_flags": "|".join(flags),
            }
        )
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


def observation_receipt(table: pd.DataFrame, channel_id: str, effort_tier: str) -> dict[str, Any]:
    counts = table["observation_state"].value_counts().to_dict() if not table.empty else {}
    return {
        "channel_id": channel_id,
        "effort_tier": effort_tier,
        "n_islands": int(len(table)),
        "n_detected": int(counts.get("detected", 0)),
        "n_adequate_non_detection": int(counts.get("adequate_non_detection", 0)),
        "n_insufficient_effort": int(counts.get("insufficient_effort", 0)),
        "uses_focal_plant_traits": False,
        "detection_claim": "potential_partner_channel_present_not_realized_service",
    }


@app.command("classify")
def classify_command(
    records_csv: Path = typer.Option(..., exists=True),
    catalog_csv: Path = typer.Option(..., exists=True),
    channel_id: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
    effort_tier: str = typer.Option("primary"),
    policy_path: Path = typer.Option(Path("config/chapter1_nee_channel_observation_policy.yml")),
) -> None:
    policy = load_policy(policy_path)
    records = pd.read_csv(records_csv, dtype=str).fillna("")
    catalog = pd.read_csv(catalog_csv, dtype=str).fillna("")
    result = classify_channel_observations(
        records, catalog, channel_id, policy, effort_tier=effort_tier
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_dir / f"{channel_id}_island_observation.csv", index=False)
    receipt = observation_receipt(result, channel_id, effort_tier)
    (output_dir / f"{channel_id}_observation_receipt.json").write_text(
        json.dumps(receipt, indent=2), encoding="utf-8"
    )
    typer.echo(json.dumps(receipt))


if __name__ == "__main__":
    app()
