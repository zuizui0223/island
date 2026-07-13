"""Validate the scientific and structural invariants of the v2 main-analysis input."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = {
    "island_id",
    "analysis_regime",
    "area_km2",
    "distance_to_continent_km",
    "bombus_deficit",
    "color_trials",
    "color_plain",
    "form_trials",
    "form_open_generalized",
    "sc_trials",
    "sc_successes",
}


def validate_analysis_inputs(csv_path: Path, expected_islands: int | None = None) -> None:
    data = pd.read_csv(csv_path)
    missing = REQUIRED_COLUMNS.difference(data.columns)
    if missing:
        raise ValueError(f"missing required columns: {sorted(missing)}")
    if not data["island_id"].is_unique:
        raise ValueError("island_id must be unique")
    if expected_islands is not None and len(data) != expected_islands:
        raise ValueError(f"expected {expected_islands} islands, found {len(data)}")

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
    parser.add_argument("--expected-islands", type=int)
    args = parser.parse_args()
    validate_analysis_inputs(args.input_csv, args.expected_islands)


if __name__ == "__main__":
    main()
