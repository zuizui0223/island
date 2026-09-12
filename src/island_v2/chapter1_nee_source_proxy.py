"""Aggregate frozen PR138 mainland source assignments into N1 channel source states.

The source assignment modes predate pollinator outcomes.  Primary N1 uses geo_k5:
a channel is available if at least one of the five selected mainland source entities
has accepted positive channel evidence. Structural absence requires explicit absence
for every selected source entity. Partial or incomplete information remains unresolved.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

ASSIGNMENT_REQUIRED = {"island_id", "source_mode", "source_rank", "entity_ID"}
ENTITY_STATE_REQUIRED = {
    "entity_ID",
    "channel_id",
    "source_state",
    "evidence_id",
    "evidence_type",
    "source_citation",
    "source_url",
    "review_status",
}
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
AUDIT_COLUMNS = [
    "island_id",
    "source_mode",
    "channel_id",
    "source_context_entity_ID",
    "selected_entity_IDs",
    "n_selected_entities",
    "n_available_entities",
    "n_structurally_absent_entities",
    "n_unresolved_entities",
    "source_state",
]


def load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("source proxy config must be a mapping")
    if payload.get("contract") != "chapter1_nee_source_proxy_v1":
        raise typer.BadParameter("unexpected source proxy contract")
    return payload


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _require(table: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = columns.difference(table.columns)
    if missing:
        raise ValueError(f"{label} missing columns: {sorted(missing)}")


def _join(values: pd.Series) -> str:
    return "|".join(sorted({value for value in _text(values) if value}))


def source_mode_spec(config: dict[str, Any], mode: str) -> dict[str, Any]:
    primary = config["primary_source_proxy"]
    if mode == str(primary["source_mode"]):
        return {"expected_k": int(primary["expected_k"]), "role": "primary"}
    sensitivity = config["sensitivity_source_proxies"]
    if mode not in sensitivity or not isinstance(sensitivity[mode], dict):
        raise ValueError(f"source mode is not frozen in the N1 source proxy: {mode}")
    return {"expected_k": int(sensitivity[mode]["expected_k"]), "role": "sensitivity"}


def validate_assignments(assignments: pd.DataFrame, config: dict[str, Any], mode: str) -> pd.DataFrame:
    _require(assignments, ASSIGNMENT_REQUIRED, "source assignments")
    spec = source_mode_spec(config, mode)
    work = assignments.loc[assignments["source_mode"].astype(str).eq(mode)].copy()
    for column in ["island_id", "source_mode", "entity_ID"]:
        work[column] = _text(work[column])
    work["source_rank"] = pd.to_numeric(work["source_rank"], errors="coerce")
    if work.empty:
        raise ValueError(f"no frozen source assignments for mode {mode}")
    if work[["island_id", "source_mode", "source_rank"]].duplicated().any():
        raise ValueError("duplicate island x source_mode x source_rank assignments")
    if work[["island_id", "entity_ID"]].duplicated().any():
        raise ValueError("same source entity appears more than once within an island source set")
    expected = set(range(1, int(spec["expected_k"]) + 1))
    for island_id, group in work.groupby("island_id", sort=False):
        observed = {int(value) for value in group["source_rank"].dropna()}
        if observed != expected:
            raise ValueError(
                f"island {island_id} has incomplete/unexpected {mode} ranks: "
                f"expected {sorted(expected)}, found {sorted(observed)}"
            )
        if group["entity_ID"].eq("").any():
            raise ValueError(f"island {island_id} has blank source entity ID")
    return work.sort_values(["island_id", "source_rank"]).reset_index(drop=True)


def validate_entity_states(entity_states: pd.DataFrame, channels: list[str]) -> pd.DataFrame:
    _require(entity_states, ENTITY_STATE_REQUIRED, "source-entity channel states")
    work = entity_states[list(ENTITY_STATE_REQUIRED)].copy()
    for column in ENTITY_STATE_REQUIRED:
        work[column] = _text(work[column])
    if work[["entity_ID", "channel_id"]].duplicated().any():
        raise ValueError("source-entity states must be unique by entity_ID x channel_id")
    invalid_channels = sorted(set(work["channel_id"]).difference(set(channels)))
    if invalid_channels:
        raise ValueError(f"unregistered source-entity channels: {invalid_channels}")
    allowed_states = {"available", "structurally_absent", "unresolved"}
    invalid_states = sorted(set(work["source_state"]).difference(allowed_states))
    if invalid_states:
        raise ValueError(f"invalid source-entity state: {invalid_states}")
    decisive = work["source_state"].isin({"available", "structurally_absent"})
    if (decisive & ~work["review_status"].eq("accepted")).any():
        raise ValueError("decisive source-entity states require accepted review")
    return work


def aggregate_source_mode(
    assignments: pd.DataFrame,
    entity_states: pd.DataFrame,
    config: dict[str, Any],
    mode: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate one frozen source mode to island x channel source availability."""
    selected = validate_assignments(assignments, config, mode)
    channels = list(config["channels"] if "channels" in config else [
        "bombus", "non_bombus_bees", "lepidoptera", "flower_visiting_birds", "diptera"
    ])
    states = validate_entity_states(entity_states, channels)
    lookup = states.set_index(["entity_ID", "channel_id"], drop=False)
    pending_review = "pending"

    output_rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, Any]] = []
    for island_id, group in selected.groupby("island_id", sort=True):
        group = group.sort_values("source_rank")
        entities = group["entity_ID"].astype(str).tolist()
        context = str(group.iloc[0]["entity_ID"])
        for channel in channels:
            records: list[dict[str, str]] = []
            for entity in entities:
                key = (entity, channel)
                if key in lookup.index:
                    row = lookup.loc[key]
                    if isinstance(row, pd.DataFrame):
                        raise ValueError("duplicate source-entity state after validation")
                    records.append(row.to_dict())
                else:
                    records.append(
                        {
                            "entity_ID": entity,
                            "channel_id": channel,
                            "source_state": "unresolved",
                            "evidence_id": "",
                            "evidence_type": "",
                            "source_citation": "",
                            "source_url": "",
                            "review_status": pending_review,
                        }
                    )
            evidence = pd.DataFrame(records)
            n_available = int(evidence["source_state"].eq("available").sum())
            n_absent = int(evidence["source_state"].eq("structurally_absent").sum())
            n_unresolved = int(evidence["source_state"].eq("unresolved").sum())
            if n_available >= 1:
                state = "available"
                review = "accepted"
            elif n_absent == len(entities):
                state = "structurally_absent"
                review = "accepted"
            else:
                state = "unresolved"
                review = pending_review
            relevant = evidence.loc[
                evidence["source_state"].eq("available")
                if state == "available"
                else evidence["source_state"].eq("structurally_absent")
                if state == "structurally_absent"
                else evidence["source_state"].ne("")
            ]
            output_rows.append(
                {
                    "island_id": island_id,
                    "source_region_id": context,
                    "channel_id": channel,
                    "source_state": state,
                    "evidence_id": _join(relevant["evidence_id"]),
                    "evidence_type": _join(relevant["evidence_type"]),
                    "source_citation": _join(relevant["source_citation"]),
                    "source_url": _join(relevant["source_url"]),
                    "review_status": review,
                }
            )
            audit_rows.append(
                {
                    "island_id": island_id,
                    "source_mode": mode,
                    "channel_id": channel,
                    "source_context_entity_ID": context,
                    "selected_entity_IDs": "|".join(entities),
                    "n_selected_entities": len(entities),
                    "n_available_entities": n_available,
                    "n_structurally_absent_entities": n_absent,
                    "n_unresolved_entities": n_unresolved,
                    "source_state": state,
                }
            )
    return (
        pd.DataFrame(output_rows, columns=OUTPUT_COLUMNS),
        pd.DataFrame(audit_rows, columns=AUDIT_COLUMNS),
    )


def receipt(output: pd.DataFrame, audit: pd.DataFrame, mode: str, config: dict[str, Any]) -> dict[str, Any]:
    spec = source_mode_spec(config, mode)
    counts = (
        output.groupby(["channel_id", "source_state"]).size().rename("n").reset_index().to_dict("records")
        if not output.empty
        else []
    )
    return {
        "contract": config["contract"],
        "source_mode": mode,
        "role": spec["role"],
        "expected_k": spec["expected_k"],
        "n_islands": int(output["island_id"].nunique()) if not output.empty else 0,
        "n_rows": int(len(output)),
        "state_counts": counts,
        "source_context_definition": "rank1_entity_ID",
        "uses_pollinator_retention_outcomes": False,
        "uses_focal_plant_traits": False,
        "all_selected_source_entities_preserved_in_audit": bool(
            not audit.empty and audit["n_selected_entities"].eq(int(spec["expected_k"])).all()
        ),
    }


@app.command("aggregate")
def aggregate_command(
    assignments_csv: Path = typer.Option(..., exists=True),
    entity_states_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    source_mode: str = typer.Option("geo_k5"),
    config_path: Path = typer.Option(Path("config/chapter1_nee_source_proxy.yml")),
) -> None:
    config = load_config(config_path)
    assignments = pd.read_csv(assignments_csv, dtype=str).fillna("")
    entity_states = pd.read_csv(entity_states_csv, dtype=str).fillna("")
    output, audit = aggregate_source_mode(assignments, entity_states, config, source_mode)
    output_dir.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_dir / f"{source_mode}_island_channel_source_availability.csv", index=False)
    audit.to_csv(output_dir / f"{source_mode}_source_set_audit.csv", index=False)
    report = receipt(output, audit, source_mode, config)
    (output_dir / f"{source_mode}_source_proxy_receipt.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    typer.echo(json.dumps(report))


if __name__ == "__main__":
    app()
