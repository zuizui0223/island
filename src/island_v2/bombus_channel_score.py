"""Continuous Bombus channel metrics from occurrence effort and climate compatibility.

The score is intentionally component-first. Observation evidence and environmental
compatibility remain separate columns; a combined availability/deficit score is a
secondary summary for v1-style continuous analyses.
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


@app.callback()
def main() -> None:
    """Build continuous, component-first Bombus channel metrics."""


def _clip01(value: float) -> float:
    return max(0.0, min(1.0, float(value)))


def _ratio(value: object, threshold: float) -> float:
    number = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(number) or threshold <= 0:
        return 0.0
    return _clip01(float(number) / threshold)


def aggregate_environmental_compatibility(species_scores: pd.DataFrame) -> pd.DataFrame:
    """Aggregate island x source-pool Bombus compatibility without hiding species scores.

    Input columns: island_id, bombus_species, environmental_compatibility (0..1).
    Maximum compatibility is the primary summary because species-level niche
    similarities are not calibrated independent probabilities. Mean and noisy-OR
    are retained as explicitly named sensitivity summaries. Missing species models
    remain visible in coverage counts instead of silently removing the island.
    """
    required = {"island_id", "bombus_species", "environmental_compatibility"}
    missing = required.difference(species_scores.columns)
    if missing:
        raise typer.BadParameter(
            f"environmental score table missing columns: {sorted(missing)}"
        )

    work = species_scores.copy()
    for column in ("island_id", "bombus_species"):
        work[column] = work[column].fillna("").astype(str).str.strip()
    if work[["island_id", "bombus_species"]].eq("").any(axis=None):
        raise typer.BadParameter(
            "environmental score table contains blank island_id or bombus_species"
        )
    if work.duplicated(["island_id", "bombus_species"]).any():
        raise typer.BadParameter(
            "environmental score table must have one row per island_id x bombus_species"
        )

    raw_scores = work["environmental_compatibility"].fillna("").astype(str).str.strip()
    numeric_scores = pd.to_numeric(raw_scores, errors="coerce")
    invalid_numeric = raw_scores.ne("") & numeric_scores.isna()
    if invalid_numeric.any():
        raise typer.BadParameter(
            "environmental_compatibility contains non-numeric non-missing values"
        )
    outside_unit_interval = numeric_scores.notna() & ~numeric_scores.between(0, 1)
    if outside_unit_interval.any():
        raise typer.BadParameter("environmental_compatibility must be between 0 and 1")
    work["environmental_compatibility"] = numeric_scores

    rows: list[dict[str, Any]] = []
    for island_id, group in work.groupby("island_id", sort=True):
        scored = group["environmental_compatibility"].dropna().astype(float)
        values = scored.tolist()
        n_total = int(group["bombus_species"].nunique())
        n_scored = int(scored.size)
        if values:
            maximum = max(values)
            mean = sum(values) / len(values)
            noisy_or = 1.0 - math.prod(1.0 - value for value in values)
            status = "scored_all" if n_scored == n_total else "scored_partial"
        else:
            maximum = pd.NA
            mean = pd.NA
            noisy_or = pd.NA
            status = "no_scored_species"
        rows.append(
            {
                "island_id": island_id,
                "n_source_pool_bombus_species": n_total,
                "n_scored_source_pool_bombus_species": n_scored,
                "n_unscored_source_pool_bombus_species": n_total - n_scored,
                "environmental_compatibility": maximum,
                "environmental_compatibility_max": maximum,
                "environmental_compatibility_mean": mean,
                "environmental_compatibility_noisy_or": noisy_or,
                "environmental_compatibility_status": status,
            }
        )
    return pd.DataFrame(rows)


def score_observation_evidence(
    diagnostics: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Convert Bombus records plus effort/recency into a neutral-at-low-effort score.

    Low observation effort keeps a zero-detection score near 0.5 rather than
    treating it as absence. Any quality-filtered detection remains above 0.5
    regardless of background effort, and repeated detections approach 1.
    """
    thresholds = config["primary_thresholds"]
    recency = config["observation_recency"]
    count_scale = float(
        config.get("continuous_score", {}).get("bombus_record_scale", 3.0)
    )
    max_age = float(recency["max_years_since_latest_background"])

    rows: list[dict[str, Any]] = []
    for _, row in diagnostics.iterrows():
        bombus_value = pd.to_numeric(
            pd.Series([row.get("bombus_record_count", 0)]),
            errors="coerce",
        ).fillna(0)
        bombus_count = max(0.0, float(bombus_value.iloc[0]))
        occurrence_support = 1.0 - math.exp(
            -bombus_count / max(count_scale, 1e-9)
        )

        effort_parts = [
            _ratio(
                row.get("target_group_record_count", 0),
                float(thresholds["min_background_records"]),
            ),
            _ratio(
                row.get("target_group_spatial_units", 0),
                float(thresholds["min_background_spatial_units"]),
            ),
            _ratio(
                row.get("target_group_temporal_units", 0),
                float(thresholds["min_background_temporal_units"]),
            ),
            _ratio(
                row.get("distinct_dataset_count", 0),
                float(thresholds["min_distinct_datasets"]),
            ),
        ]
        structural_effort = math.prod(effort_parts) ** (1.0 / len(effort_parts))

        age = pd.to_numeric(
            pd.Series([row.get("latest_record_age_years", pd.NA)]),
            errors="coerce",
        ).iloc[0]
        if pd.isna(age) or age < 0:
            recency_support = 0.0
        else:
            recency_support = _clip01(1.0 - float(age) / max(max_age, 1.0))

        effort_quality = math.sqrt(structural_effort * recency_support)
        if bombus_count > 0:
            # A valid detection is positive occurrence evidence even when the
            # target-group background is sparse or old. Effort/recency governs
            # how strongly a zero can support non-detection, not whether an
            # observed Bombus record counts as a detection.
            observation_availability = 0.5 + 0.5 * occurrence_support
            score_confidence = observation_availability
        else:
            observation_availability = 0.5 * (1.0 - effort_quality)
            score_confidence = effort_quality

        rows.append(
            {
                "island_id": row["island_id"],
                "bombus_record_count": int(bombus_count),
                "occurrence_support": _clip01(occurrence_support),
                "observation_effort_quality": _clip01(effort_quality),
                "recency_support": _clip01(recency_support),
                "observation_availability": _clip01(observation_availability),
                "observation_deficit": 1.0 - _clip01(observation_availability),
                "observation_score_confidence": _clip01(score_confidence),
                "bombus_occurrence_evidence": row.get(
                    "bombus_occurrence_evidence", ""
                ),
            }
        )
    return pd.DataFrame(rows)


def combine_channel_components(
    observation: pd.DataFrame,
    environmental: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Join components and create a secondary combined availability/deficit score."""
    out = observation.copy()
    if environmental is None or environmental.empty:
        out["environmental_compatibility"] = pd.NA
        out["bombus_channel_availability"] = pd.NA
        out["bombus_channel_deficit"] = pd.NA
        return out

    out = out.merge(environmental, on="island_id", how="left", validate="one_to_one")
    env = pd.to_numeric(out["environmental_compatibility"], errors="coerce")
    obs = pd.to_numeric(out["observation_availability"], errors="coerce")
    combined = (env * obs).pow(0.5)
    out["bombus_channel_availability"] = combined
    out["bombus_channel_deficit"] = 1.0 - combined
    return out


@app.command("run")
def run(
    diagnostics_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/bombus_observation_diagnostics.yml")),
    environmental_scores_csv: Path | None = typer.Option(
        None,
        help="Optional island x source-pool Bombus environmental compatibility table.",
    ),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    diagnostics = pd.read_csv(diagnostics_csv, dtype=str).fillna("")
    observation = score_observation_evidence(diagnostics, config)
    environmental = None
    if environmental_scores_csv is not None:
        species_scores = pd.read_csv(
            environmental_scores_csv,
            dtype=str,
        ).fillna("")
        environmental = aggregate_environmental_compatibility(species_scores)

    result = combine_channel_components(observation, environmental)
    output_dir.mkdir(parents=True, exist_ok=True)
    observation.to_csv(
        output_dir / "bombus_observation_continuous_scores.csv",
        index=False,
    )
    if environmental is not None:
        environmental.to_csv(
            output_dir / "bombus_environmental_compatibility.csv",
            index=False,
        )
    result.to_csv(
        output_dir / "bombus_channel_continuous_scores.csv",
        index=False,
    )
    summary = {
        "n_islands": int(len(result)),
        "component_policy": (
            "keep environmental compatibility and observation evidence separate; "
            "combined score is secondary"
        ),
        "environmental_primary_aggregation": "maximum species compatibility",
        "environmental_sensitivity_aggregations": ["mean", "noisy_or"],
        "combined_formula": (
            "sqrt(environmental_compatibility * observation_availability)"
        ),
        "observation_formula": {
            "detected": "0.5 + 0.5 * occurrence_support",
            "zero_detection": "0.5 * (1 - effort_quality)",
        },
    }
    (output_dir / "bombus_channel_score_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary))


if __name__ == "__main__":
    app()
