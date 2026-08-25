"""Materialize a review-ready, outcome-blind PR136 source-exposure curation sheet.

The scaffold converts the support-only source-exposure queue into a human-reviewable
sheet with geography needed to identify an island/source region and empty evidence
fields needed by the existing Bombus applicability curation contract. It never fills
or guesses source exposure.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

PROHIBITED_PRIORITY_COLUMNS = {
    "deviation_observed_minus_null",
    "trait_state",
    "trait_state_value",
    "effect_direction",
    "p_value",
    "q_value",
    "pollination_syndrome_match",
}

GEO_COLUMNS = [
    "island_latitude",
    "island_longitude",
    "area_km2",
    "distance_to_continent_km",
    "log_distance_to_continent_km",
    "spatial_block",
    "analysis_regime",
]

REVIEW_COLUMNS = [
    "source_region_id",
    "biogeographic_region",
    "floristic_source_pool",
    "plausible_dispersal_connection_or_gateway",
    "documented_native_bombus_status",
    "proposed_applicability",
    "source_region_evidence_id",
    "assignment_evidence_id",
    "source_citation",
    "source_url",
    "source_locator",
    "evidence_excerpt",
    "review_status",
    "reviewer",
    "decision_date",
]


def build_review_scaffold(
    queue: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    wave: int | None = 1,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    missing = {"island_id", "curation_wave"} - set(queue.columns)
    if missing:
        raise typer.BadParameter(f"queue missing columns: {sorted(missing)}")
    prohibited = PROHIBITED_PRIORITY_COLUMNS.intersection(queue.columns)
    if prohibited:
        raise typer.BadParameter(
            f"queue contains prohibited outcome columns: {sorted(prohibited)}"
        )
    if "island_id" not in covariates.columns:
        raise typer.BadParameter("covariates missing island_id")

    work = queue.copy()
    if wave is not None:
        work = work.loc[pd.to_numeric(work["curation_wave"], errors="coerce").eq(wave)].copy()

    available_geo = [column for column in GEO_COLUMNS if column in covariates.columns]
    geo = covariates[["island_id", *available_geo]].drop_duplicates("island_id")
    if geo["island_id"].duplicated().any():
        raise typer.BadParameter("covariates must contain one row per island_id")
    # Avoid duplicate columns already retained by the queue; the canonical covariates
    # are authoritative for the review scaffold.
    for column in available_geo:
        if column in work.columns:
            work = work.drop(columns=column)
    work = work.merge(geo, on="island_id", how="left", validate="one_to_one")

    for column in REVIEW_COLUMNS:
        work[column] = "pending_review" if column == "review_status" else ""
    work["source_exposure_decision_is_outcome_blind"] = True
    work["trait_values_available_to_reviewer"] = False
    work["review_instruction"] = (
        "Assign source region and native Bombus status from biogeographic/floristic "
        "evidence only; do not inspect floral/reproductive outcomes."
    )

    manifest = {
        "contract": "chapter1_pr136_source_exposure_review_scaffold_v1",
        "wave": wave,
        "n_rows": int(len(work)),
        "geography_columns_attached": available_geo,
        "review_fields": REVIEW_COLUMNS,
        "trait_values_exposed_to_reviewer": False,
        "residual_effects_exposed_to_reviewer": False,
        "applicability_auto_assigned": False,
        "interpretation": (
            "Review scaffold only. Every source-region/applicability decision still "
            "requires independent evidence and review status."
        ),
    }
    return work, manifest


@app.command("run")
def run(
    queue_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    wave: int = typer.Option(1, min=1),
) -> None:
    scaffold, manifest = build_review_scaffold(
        pd.read_csv(queue_csv), pd.read_csv(covariates_csv), wave=wave
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    scaffold.to_csv(output_dir / "source_exposure_review_scaffold.csv", index=False)
    (output_dir / "source_exposure_review_scaffold_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
