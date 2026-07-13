"""Validate scientific and structural invariants of the v2 analysis input."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = {
    "island_id",
    "analysis_tier",
    "analysis_regime",
    "area_km2",
    "distance_to_continent_km",
    "distance_zero_flag",
    "primary_distance_support",
    "bombus_deficit",
    "color_trials",
    "color_plain",
    "form_trials",
    "form_open_generalized",
    "sc_trials",
    "sc_successes",
}


def _read_reference_ids(path: Path) -> pd.Index:
    reference = pd.read_csv(path, usecols=["island_id"])
    if reference["island_id"].duplicated().any():
        raise ValueError("reference island universe contains duplicate island_id values")
    return pd.Index(reference["island_id"].astype(str), name="island_id")


def validate_analysis_inputs(
    csv_path: Path,
    *,
    expected_tier: str | None = None,
    reference_islands_csv: Path | None = None,
    expected_islands: int | None = None,
) -> None:
    data = pd.read_csv(csv_path)
    missing = REQUIRED_COLUMNS.difference(data.columns)
    if missing:
        raise ValueError(f"missing required columns: {sorted(missing)}")
    if not data["island_id"].is_unique:
        raise ValueError("island_id must be unique")

    if expected_tier is not None:
        tiers = set(data["analysis_tier"].dropna().astype(str).unique())
        if tiers != {expected_tier}:
            raise ValueError(f"expected analysis_tier={expected_tier!r}, found {sorted(tiers)}")

    if reference_islands_csv is not None:
        reference_ids = _read_reference_ids(reference_islands_csv)
        observed_ids = pd.Index(data["island_id"].astype(str), name="island_id")
        missing_ids = reference_ids.difference(observed_ids)
        extra_ids = observed_ids.difference(reference_ids)
        if len(missing_ids) or len(extra_ids):
            raise ValueError(
                "analysis island universe differs from reference: "
                f"missing={len(missing_ids)}, extra={len(extra_ids)}"
            )
    elif expected_islands is not None and len(data) != expected_islands:
        # Backward-compatible fallback for older workflows. New canonical workflows
        # should validate against the frozen/reference island table instead.
        raise ValueError(f"expected {expected_islands} islands, found {len(data)}")

    if data["distance_to_continent_km"].isna().any():
        raise ValueError("distance_to_continent_km contains missing values")
    if (data["distance_to_continent_km"] < 0).any():
        raise ValueError("distance_to_continent_km must be non-negative")

    expected_zero = data["distance_to_continent_km"].eq(0)
    observed_zero = data["distance_zero_flag"].astype(bool)
    if not expected_zero.equals(observed_zero):
        raise ValueError("distance_zero_flag is inconsistent with distance_to_continent_km")

    expected_primary = data["distance_to_continent_km"].gt(0)
    observed_primary = data["primary_distance_support"].astype(bool)
    if not expected_primary.equals(observed_primary):
        raise ValueError("primary_distance_support is inconsistent with distance_to_continent_km")

    checks = [
        ("color_plain", "color_trials"),
        ("form_open_generalized", "form_trials"),
        ("sc_successes", "sc_trials"),
    ]
    for successes, trials in checks:
        bad = data[successes].fillna(0) > data[trials].fillna(0)
        if bad.any():
            raise ValueError(f"{successes} exceeds {trials} for {int(bad.sum())} islands")

    observed = data["bombus_deficit"].dropna()
    if ((observed < 0) | (observed > 1)).any():
        raise ValueError("bombus_deficit must lie in [0, 1]")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-csv", type=Path, required=True)
    parser.add_argument("--expected-tier")
    parser.add_argument("--reference-islands-csv", type=Path)
    parser.add_argument("--expected-islands", type=int)
    args = parser.parse_args()
    validate_analysis_inputs(
        args.input_csv,
        expected_tier=args.expected_tier,
        reference_islands_csv=args.reference_islands_csv,
        expected_islands=args.expected_islands,
    )


if __name__ == "__main__":
    main()
