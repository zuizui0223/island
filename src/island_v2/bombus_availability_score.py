"""Continuous Bombus occurrence-availability evidence from island diagnostics.

This score is an observation-evidence axis, not a claim of biological occupancy.
It preserves the v2 distinction between detection and observation effort while
providing a v1-like continuous predictor for downstream models.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

SCORE_COLUMNS = [
    "bombus_detection_support",
    "bombus_observation_coverage",
    "bombus_recency_support",
    "bombus_occurrence_availability_score",
    "bombus_occurrence_deficit_score",
    "bombus_occurrence_score_confidence",
]


def _num(frame: pd.DataFrame, column: str) -> pd.Series:
    return pd.to_numeric(frame.get(column, pd.Series(index=frame.index)), errors="coerce")


def _cap_ratio(values: pd.Series, threshold: float) -> pd.Series:
    return (values.fillna(0.0).clip(lower=0.0) / float(threshold)).clip(upper=1.0)


def score_occurrence_evidence(
    diagnostics: pd.DataFrame,
    config: dict,
    detection_scale: float = 3.0,
) -> pd.DataFrame:
    """Add continuous 0-1 Bombus availability and deficit evidence scores.

    Components are kept separate. A detected island is mapped to 0.5-1.0 as
    record support increases. A zero-detection island is mapped from 0.5 toward
    0.0 only as observation coverage becomes strong. Thus no records plus no
    effort remains uncertain rather than being treated as absence.
    """
    out = diagnostics.copy()
    thresholds = config["primary_thresholds"]
    recency = config["observation_recency"]

    bombus_count = _num(out, "bombus_record_count").fillna(0.0).clip(lower=0.0)
    background = _num(out, "target_group_record_count")
    spatial = _num(out, "target_group_spatial_units")
    temporal = _num(out, "target_group_temporal_units")
    datasets = _num(out, "distinct_dataset_count")
    age = _num(out, "latest_record_age_years")

    count_support = _cap_ratio(background, thresholds["min_background_records"])
    spatial_support = _cap_ratio(spatial, thresholds["min_background_spatial_units"])
    temporal_support = _cap_ratio(temporal, thresholds["min_background_temporal_units"])
    dataset_support = _cap_ratio(datasets, thresholds["min_distinct_datasets"])

    # Geometric mean: one severely weak observation dimension cannot be hidden by
    # another very large dimension.
    coverage = pd.Series(
        [
            math.prod(max(float(v), 1e-12) for v in values) ** 0.25
            if any(float(v) > 0 for v in values)
            else 0.0
            for values in zip(
                count_support,
                spatial_support,
                temporal_support,
                dataset_support,
                strict=True,
            )
        ],
        index=out.index,
        dtype=float,
    )

    max_age = float(recency["max_years_since_latest_background"])
    recency_support = (1.0 - age / max_age).clip(lower=0.0, upper=1.0).fillna(0.0)
    nondetection_coverage = coverage * recency_support

    # Saturating evidence from repeated detections.
    detection_support = 1.0 - (-bombus_count / float(detection_scale)).map(math.exp)
    detected = bombus_count.gt(0)

    availability = pd.Series(index=out.index, dtype=float)
    availability.loc[detected] = 0.5 + 0.5 * detection_support.loc[detected]
    availability.loc[~detected] = 0.5 * (1.0 - nondetection_coverage.loc[~detected])
    availability = availability.clip(lower=0.0, upper=1.0)

    confidence = pd.Series(index=out.index, dtype=float)
    confidence.loc[detected] = (0.5 + 0.5 * detection_support.loc[detected]).clip(upper=1.0)
    confidence.loc[~detected] = nondetection_coverage.loc[~detected]

    out["bombus_detection_support"] = detection_support.round(6)
    out["bombus_observation_coverage"] = coverage.round(6)
    out["bombus_recency_support"] = recency_support.round(6)
    out["bombus_occurrence_availability_score"] = availability.round(6)
    out["bombus_occurrence_deficit_score"] = (1.0 - availability).round(6)
    out["bombus_occurrence_score_confidence"] = confidence.round(6)
    return out


def score_summary(scored: pd.DataFrame) -> dict:
    return {
        "n_islands": int(len(scored)),
        "availability_mean": float(scored["bombus_occurrence_availability_score"].mean()),
        "availability_median": float(scored["bombus_occurrence_availability_score"].median()),
        "deficit_mean": float(scored["bombus_occurrence_deficit_score"].mean()),
        "confidence_mean": float(scored["bombus_occurrence_score_confidence"].mean()),
        "interpretation": (
            "0=strong occurrence-deficit evidence; 0.5=uninformative/weak evidence; "
            "1=strong occurrence-availability evidence. This is not biological occupancy probability."
        ),
    }


@app.command("run")
def run(
    diagnostics_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/bombus_observation_diagnostics.yml")),
) -> None:
    diagnostics = pd.read_csv(diagnostics_csv).fillna("")
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    scored = score_occurrence_evidence(diagnostics, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    scored.to_csv(output_dir / "island_bombus_occurrence_scores.csv", index=False)
    (output_dir / "bombus_occurrence_score_summary.json").write_text(
        json.dumps(score_summary(scored), indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(f"Scored {len(scored)} islands on continuous Bombus occurrence evidence.")


if __name__ == "__main__":
    app()
