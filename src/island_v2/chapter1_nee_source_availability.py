"""Build outcome-blind source-region pollination-channel availability for N1.

Positive evidence can establish that a frozen functional channel is available in a
source region. Structural absence is deliberately harder: it requires explicit,
accepted biogeographic evidence. Zero occurrence records, island non-detection,
climate compatibility and floral traits never create source absence.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

CATALOG_REQUIRED = {"channel_id", "pollinator_species", "catalog_tier"}
OUTPUT_COLUMNS = [
    "island_id",
    "source_region_id",
    "channel_id",
    "source_state",
    "evidence_id",
    "evidence_type",
    "source_citation",
    "source_url",
    "review_status",
]


def load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("source availability config must be a mapping")
    if payload.get("contract") != "chapter1_nee_source_availability_v1":
        raise typer.BadParameter("unexpected source availability contract")
    return payload


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _require(table: pd.DataFrame, columns: list[str] | set[str], label: str) -> None:
    missing = set(columns).difference(table.columns)
    if missing:
        raise ValueError(f"{label} missing columns: {sorted(missing)}")


def _normalize(table: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    result = table.copy()
    for column in columns:
        result[column] = _text(result[column])
    return result


def _join_unique(values: pd.Series) -> str:
    return "|".join(sorted({value for value in _text(values) if value}))


def confirmatory_species_by_channel(catalog: pd.DataFrame, config: dict[str, Any]) -> dict[str, set[str]]:
    _require(catalog, CATALOG_REQUIRED, "channel catalog")
    allowed = set(config["channels"])
    work = catalog.loc[
        catalog["catalog_tier"].astype(str).eq("confirmatory")
        & catalog["channel_id"].astype(str).isin(allowed)
    ].copy()
    result: dict[str, set[str]] = {channel: set() for channel in sorted(allowed)}
    for channel, group in work.groupby("channel_id", sort=True):
        result[str(channel)] = {value for value in _text(group["pollinator_species"]) if value}
    return result


def validate_assignments(assignments: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    required = list(config["assignment"]["required_columns"])
    _require(assignments, required, "source assignments")
    result = _normalize(assignments[required], required)
    if result["island_id"].eq("").any():
        raise ValueError("source assignments contain blank island_id")
    if result["island_id"].duplicated().any():
        raise ValueError("source assignments must contain one row per island_id")
    valid_review = set(config["review_status"].values())
    invalid = sorted(set(result["review_status"]).difference(valid_review))
    if invalid:
        raise ValueError(f"source assignments contain invalid review status: {invalid}")
    return result


def validate_positive_evidence(table: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    required = list(config["positive_evidence"]["required_columns"])
    _require(table, required, "positive source evidence")
    result = _normalize(table[required], required)
    accepted_types = set(config["positive_evidence"]["accepted_types"])
    accepted_review = str(config["review_status"]["accepted"])
    result = result.loc[
        result["review_status"].eq(accepted_review)
        & result["evidence_type"].isin(accepted_types)
    ].copy()
    if result[["source_region_id", "pollinator_species", "evidence_id"]].eq("").any().any():
        raise ValueError("accepted positive source evidence requires region, species and evidence_id")
    if result["source_citation"].eq("").any():
        raise ValueError("accepted positive source evidence requires source_citation")
    return result


def validate_structural_absence(table: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    required = list(config["structural_absence_evidence"]["required_columns"])
    _require(table, required, "structural source-absence evidence")
    result = _normalize(table[required], required)
    accepted_types = set(config["structural_absence_evidence"]["accepted_types"])
    accepted_review = str(config["review_status"]["accepted"])
    # Prohibited types are not merely ignored: an accepted row using them is an invalid
    # attempt to turn missing/indirect evidence into structural absence.
    prohibited = set(config["structural_absence_evidence"]["prohibited_types"])
    invalid_accepted = result.loc[
        result["review_status"].eq(accepted_review) & result["evidence_type"].isin(prohibited)
    ]
    if not invalid_accepted.empty:
        raise ValueError("prohibited evidence type cannot establish structural source absence")
    result = result.loc[
        result["review_status"].eq(accepted_review)
        & result["evidence_type"].isin(accepted_types)
    ].copy()
    if result[["source_region_id", "channel_id", "evidence_id"]].eq("").any().any():
        raise ValueError("accepted structural absence requires region, channel and evidence_id")
    if result["source_citation"].eq("").any():
        raise ValueError("accepted structural absence requires source_citation")
    unknown_channels = sorted(set(result["channel_id"]).difference(set(config["channels"])))
    if unknown_channels:
        raise ValueError(f"structural absence contains unregistered channels: {unknown_channels}")
    return result


def _source_channel_evidence(
    catalog: pd.DataFrame,
    positive: pd.DataFrame,
    structural: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    targets = confirmatory_species_by_channel(catalog, config)
    species_to_channels: dict[str, set[str]] = {}
    for channel, species_set in targets.items():
        for species in species_set:
            species_to_channels.setdefault(species, set()).add(channel)

    positive_rows: list[dict[str, str]] = []
    for row in positive.to_dict("records"):
        species = str(row["pollinator_species"])
        for channel in sorted(species_to_channels.get(species, set())):
            positive_rows.append(
                {
                    "source_region_id": str(row["source_region_id"]),
                    "channel_id": channel,
                    "evidence_id": str(row["evidence_id"]),
                    "evidence_type": str(row["evidence_type"]),
                    "source_citation": str(row["source_citation"]),
                    "source_url": str(row["source_url"]),
                    "kind": "positive",
                }
            )
    positive_channel = pd.DataFrame(positive_rows)
    structural_rows = structural.assign(kind="structural")[
        [
            "source_region_id",
            "channel_id",
            "evidence_id",
            "evidence_type",
            "source_citation",
            "source_url",
            "kind",
        ]
    ] if not structural.empty else pd.DataFrame(
        columns=[
            "source_region_id",
            "channel_id",
            "evidence_id",
            "evidence_type",
            "source_citation",
            "source_url",
            "kind",
        ]
    )
    combined = pd.concat([positive_channel, structural_rows], ignore_index=True)
    if combined.empty:
        return pd.DataFrame(
            columns=[
                "source_region_id",
                "channel_id",
                "source_state",
                "evidence_id",
                "evidence_type",
                "source_citation",
                "source_url",
                "conflict",
            ]
        )

    rows: list[dict[str, Any]] = []
    for (region, channel), group in combined.groupby(["source_region_id", "channel_id"], sort=True):
        has_positive = group["kind"].eq("positive").any()
        has_structural = group["kind"].eq("structural").any()
        conflict = bool(has_positive and has_structural)
        if conflict:
            state = "unresolved"
        elif has_positive:
            state = "available"
        elif has_structural:
            state = "structurally_absent"
        else:
            state = "unresolved"
        rows.append(
            {
                "source_region_id": region,
                "channel_id": channel,
                "source_state": state,
                "evidence_id": _join_unique(group["evidence_id"]),
                "evidence_type": _join_unique(group["evidence_type"]),
                "source_citation": _join_unique(group["source_citation"]),
                "source_url": _join_unique(group["source_url"]),
                "conflict": conflict,
            }
        )
    return pd.DataFrame(rows)


def build_island_source_availability(
    assignments: pd.DataFrame,
    catalog: pd.DataFrame,
    positive_evidence: pd.DataFrame,
    structural_absence: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return one source-state row per assigned island x frozen N0 channel."""
    assignment = validate_assignments(assignments, config)
    positive = validate_positive_evidence(positive_evidence, config)
    structural = validate_structural_absence(structural_absence, config)
    evidence = _source_channel_evidence(catalog, positive, structural, config)
    accepted_review = str(config["review_status"]["accepted"])
    pending_review = str(config["output"]["unresolved_review_status"])

    lookup = {
        (str(row.source_region_id), str(row.channel_id)): row
        for row in evidence.itertuples(index=False)
    }
    output: list[dict[str, str]] = []
    conflicts: list[dict[str, str]] = []
    for assignment_row in assignment.itertuples(index=False):
        island_id = str(assignment_row.island_id)
        region = str(assignment_row.source_region_id)
        assignment_accepted = str(assignment_row.review_status) == accepted_review and bool(region)
        for channel in config["channels"]:
            record = lookup.get((region, channel)) if assignment_accepted else None
            if not assignment_accepted or record is None:
                state = "unresolved"
                evidence_id = evidence_type = citation = url = ""
                review = pending_review
            else:
                state = str(record.source_state)
                evidence_id = str(record.evidence_id)
                evidence_type = str(record.evidence_type)
                citation = str(record.source_citation)
                url = str(record.source_url)
                review = accepted_review if state in {"available", "structurally_absent"} else pending_review
                if bool(record.conflict):
                    conflicts.append(
                        {
                            "island_id": island_id,
                            "source_region_id": region,
                            "channel_id": channel,
                            "evidence_id": evidence_id,
                        }
                    )
            output.append(
                {
                    "island_id": island_id,
                    "source_region_id": region,
                    "channel_id": channel,
                    "source_state": state,
                    "evidence_id": evidence_id,
                    "evidence_type": evidence_type,
                    "source_citation": citation,
                    "source_url": url,
                    "review_status": review,
                }
            )
    result = pd.DataFrame(output, columns=OUTPUT_COLUMNS)
    conflict_table = pd.DataFrame(
        conflicts, columns=["island_id", "source_region_id", "channel_id", "evidence_id"]
    )
    return result, conflict_table


def source_availability_receipt(table: pd.DataFrame, conflicts: pd.DataFrame) -> dict[str, Any]:
    counts = (
        table.groupby(["channel_id", "source_state"]).size().rename("n").reset_index().to_dict("records")
        if not table.empty
        else []
    )
    return {
        "n_islands": int(table["island_id"].nunique()) if not table.empty else 0,
        "n_rows": int(len(table)),
        "n_conflicts": int(len(conflicts)),
        "state_counts": counts,
        "uses_focal_plant_traits": False,
        "uses_island_channel_observation": False,
        "zero_occurrence_can_create_structural_absence": False,
    }


@app.command("build")
def build_command(
    assignments_csv: Path = typer.Option(..., exists=True),
    catalog_csv: Path = typer.Option(..., exists=True),
    positive_evidence_csv: Path = typer.Option(..., exists=True),
    structural_absence_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/chapter1_nee_source_availability.yml")),
) -> None:
    config = load_config(config_path)
    assignments = pd.read_csv(assignments_csv, dtype=str).fillna("")
    catalog = pd.read_csv(catalog_csv, dtype=str).fillna("")
    positive = pd.read_csv(positive_evidence_csv, dtype=str).fillna("")
    structural = pd.read_csv(structural_absence_csv, dtype=str).fillna("")
    result, conflicts = build_island_source_availability(
        assignments, catalog, positive, structural, config
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_dir / "island_channel_source_availability.csv", index=False)
    conflicts.to_csv(output_dir / "source_channel_evidence_conflicts.csv", index=False)
    receipt = source_availability_receipt(result, conflicts)
    (output_dir / "source_availability_receipt.json").write_text(
        json.dumps(receipt, indent=2), encoding="utf-8"
    )
    typer.echo(json.dumps(receipt))


if __name__ == "__main__":
    app()
