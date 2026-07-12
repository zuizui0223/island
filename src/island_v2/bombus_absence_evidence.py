"""Recency-aware final Bombus occurrence evidence for island analyses.

This module wraps the base effort diagnostics and adds the preregistered recency
requirement for interpreting zero Bombus records. Historical detections remain
valid detections; only adequate non-detection requires recent target-group
observation effort.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2.bombus_diagnostics import (
    compute_bombus_diagnostics,
    diagnostic_summary,
    load_config,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

RECENCY_COLUMNS = ["latest_record_age_years", "background_recent_enough"]


@app.callback()
def main() -> None:
    """Build recency-aware final Bombus occurrence evidence."""


def apply_recency_gate(diagnostics: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Require recent target-group effort before calling adequate non-detection.

    Detection is not downgraded by record age. For zero-detection islands, the
    base diagnostic must first satisfy count/spatial/temporal/dataset thresholds;
    a stale or missing latest target-group record then downgrades the state to
    insufficient_effort.
    """
    out = diagnostics.copy()
    recency = config.get("observation_recency", {})
    reference_year = int(recency["reference_year"])
    max_age = int(recency["max_years_since_latest_background"])

    latest = pd.to_numeric(
        out.get("latest_record_year", pd.Series(index=out.index)), errors="coerce"
    )
    age = reference_year - latest
    recent = latest.notna() & age.ge(0) & age.le(max_age)

    out["latest_record_age_years"] = age.astype("Int64")
    out["background_recent_enough"] = recent

    stale_non_detection = out["bombus_occurrence_evidence"].eq("adequate_non_detection") & ~recent
    out.loc[stale_non_detection, "bombus_occurrence_evidence"] = "insufficient_effort"
    out.loc[stale_non_detection, "bombus_channel_evaluable"] = False
    out.loc[stale_non_detection, "occupancy_sensitivity_eligible"] = False

    flags = (
        out.get("observation_diagnostic_flags", pd.Series("", index=out.index))
        .fillna("")
        .astype(str)
    )
    addition = pd.Series("", index=out.index, dtype=str)
    addition.loc[latest.isna()] = "no_latest_target_group_year"
    addition.loc[latest.notna() & ~recent] = "stale_background_effort"
    out["observation_diagnostic_flags"] = [
        "|".join(part for part in (base, extra) if part)
        for base, extra in zip(flags, addition, strict=True)
    ]
    return out


def final_summary(diagnostics: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    summary = diagnostic_summary(diagnostics, config)
    recency = config["observation_recency"]
    summary["observation_recency"] = {
        "reference_year": int(recency["reference_year"]),
        "max_years_since_latest_background": int(
            recency["max_years_since_latest_background"]
        ),
        "n_background_recent_enough": int(
            diagnostics.get("background_recent_enough", pd.Series(dtype=bool))
            .fillna(False)
            .sum()
        ),
    }
    return summary


@app.command("run")
def run(
    occurrences_csv: Path = typer.Option(
        ..., exists=True, help="Pre-assigned island pollinator occurrence table."
    ),
    output_dir: Path = typer.Option(
        ..., help="Directory for final Bombus occurrence evidence outputs."
    ),
    config_path: Path = typer.Option(Path("config/bombus_observation_diagnostics.yml")),
) -> None:
    """Write final effort- and recency-aware Bombus occurrence evidence."""
    config = load_config(config_path)
    if "observation_recency" not in config:
        raise typer.BadParameter("Bombus diagnostics config must contain observation_recency")
    occurrences = pd.read_csv(occurrences_csv, dtype=str).fillna("")
    base = compute_bombus_diagnostics(occurrences, config)
    diagnostics = apply_recency_gate(base, config)
    summary = final_summary(diagnostics, config)

    output_dir.mkdir(parents=True, exist_ok=True)
    diagnostics.to_csv(output_dir / "island_bombus_occurrence_evidence.csv", index=False)
    (output_dir / "bombus_occurrence_evidence_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"Classified {summary['n_islands']} islands; "
        f"{summary['n_bombus_channel_evaluable']} have evaluable Bombus occurrence evidence."
    )


if __name__ == "__main__":
    app()
