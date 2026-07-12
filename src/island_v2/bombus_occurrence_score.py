"""Continuous Bombus occurrence availability/deficit scoring.

This module converts island-level Bombus observation diagnostics into continuous
0-1 signals while preserving the component scores. It does not collapse climate
compatibility or functional-service evidence into the occurrence component.
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

REQUIRED_COLUMNS = {
    "island_id",
    "bombus_record_count",
    "target_group_record_count",
    "target_group_spatial_units",
    "target_group_temporal_units",
    "distinct_dataset_count",
    "latest_record_year",
}

SCORE_COLUMNS = [
    "island_id",
    "detection_strength",
    "effort_strength",
    "recency_strength",
    "measurement_certainty",
    "bombus_occurrence_availability_score",
    "bombus_occurrence_deficit_score",
]


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or "thresholds" not in config:
        raise typer.BadParameter("Bombus occurrence score config must contain thresholds")
    return config


def _num(value: object, default: float = 0.0) -> float:
    parsed = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    return default if pd.isna(parsed) else float(parsed)


def _cap_ratio(value: object, threshold: float) -> float:
    if threshold <= 0:
        raise ValueError("threshold must be positive")
    return min(1.0, max(0.0, _num(value) / threshold))


def _geometric_mean(values: list[float]) -> float:
    if not values:
        return 0.0
    if any(value <= 0 for value in values):
        return 0.0
    return math.prod(values) ** (1.0 / len(values))


def score_bombus_occurrence(
    diagnostics: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    missing = REQUIRED_COLUMNS.difference(diagnostics.columns)
    if missing:
        raise typer.BadParameter(
            f"Bombus diagnostics missing columns: {', '.join(sorted(missing))}"
        )

    thresholds = config["thresholds"]
    reference_year = int(config.get("reference_year", 2026))
    saturation = float(thresholds["bombus_records_saturation"])
    max_age = float(thresholds["max_years_since_latest_background"])
    if saturation <= 0 or max_age <= 0:
        raise typer.BadParameter("score saturation and recency window must be positive")

    rows: list[dict[str, float | str]] = []
    for _, row in diagnostics.iterrows():
        bombus_count = max(0.0, _num(row["bombus_record_count"]))
        detection_strength = 1.0 - math.exp(-bombus_count / saturation)

        effort_parts = [
            _cap_ratio(row["target_group_record_count"], float(thresholds["target_group_records"])),
            _cap_ratio(row["target_group_spatial_units"], float(thresholds["spatial_units"])),
            _cap_ratio(row["target_group_temporal_units"], float(thresholds["temporal_units"])),
            _cap_ratio(row["distinct_dataset_count"], float(thresholds["distinct_datasets"])),
        ]
        effort_strength = _geometric_mean(effort_parts)

        latest_year = pd.to_numeric(pd.Series([row["latest_record_year"]]), errors="coerce").iloc[0]
        if pd.isna(latest_year):
            recency_strength = 0.0
        else:
            age = max(0.0, reference_year - float(latest_year))
            recency_strength = max(0.0, 1.0 - age / max_age)

        non_detection_information = effort_strength * recency_strength
        measurement_certainty = max(detection_strength, non_detection_information)

        # With no detection and no usable effort, stay at neutral 0.5 rather than
        # pretending either presence or absence. As observation evidence becomes
        # complete, zero detections push availability toward 0. Direct detections
        # monotonically push availability upward.
        availability = detection_strength + (
            (1.0 - detection_strength) * (1.0 - non_detection_information) * 0.5
        )
        availability = min(1.0, max(0.0, availability))

        rows.append(
            {
                "island_id": str(row["island_id"]),
                "detection_strength": detection_strength,
                "effort_strength": effort_strength,
                "recency_strength": recency_strength,
                "measurement_certainty": measurement_certainty,
                "bombus_occurrence_availability_score": availability,
                "bombus_occurrence_deficit_score": 1.0 - availability,
            }
        )

    return pd.DataFrame(rows, columns=SCORE_COLUMNS)


@app.command("run")
def run(
    diagnostics_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/bombus_occurrence_score.yml")),
) -> None:
    diagnostics = pd.read_csv(diagnostics_csv, dtype=str).fillna("")
    config = load_config(config_path)
    scores = score_bombus_occurrence(diagnostics, config)

    output_dir.mkdir(parents=True, exist_ok=True)
    scores.to_csv(output_dir / "island_bombus_occurrence_scores.csv", index=False)
    summary = {
        "n_islands": int(len(scores)),
        "mean_availability": float(scores["bombus_occurrence_availability_score"].mean())
        if len(scores)
        else None,
        "mean_deficit": float(scores["bombus_occurrence_deficit_score"].mean())
        if len(scores)
        else None,
        "mean_measurement_certainty": float(scores["measurement_certainty"].mean())
        if len(scores)
        else None,
        "config_version": config.get("version"),
    }
    (output_dir / "bombus_occurrence_score_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(summary))


if __name__ == "__main__":
    app()
