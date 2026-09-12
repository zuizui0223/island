"""Standardize outcome-blind pollinator-channel evidence for the Chapter 1 NEE challenge.

This module keeps three layers separate:

1. source availability: was a channel available in the island's source region?
2. island observation: was the channel detected, adequately not detected, or not evaluable?
3. projected channel state: retained / disrupted / structurally absent / unresolved.

The separation is deliberate. Island non-detection can never create structural
absence, and focal floral traits never enter any state projection.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

SOURCE_COLUMNS = [
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
OBSERVATION_COLUMNS = [
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
PROJECTED_COLUMNS = [
    *SOURCE_COLUMNS,
    "observation_state",
    "channel_record_count",
    "background_record_count",
    "background_spatial_units",
    "background_temporal_units",
    "distinct_dataset_count",
    "latest_background_year",
    "evidence_source",
    "quality_flags",
    "channel_state",
    "projection_flags",
]

BOMBUS_SOURCE_REQUIRED = {
    "island_id",
    "source_region_id",
    "applicability",
    "source_region_evidence_id",
    "source_region_review_status",
    "assignment_review_status",
}
BOMBUS_EVIDENCE_REQUIRED = {
    "source_region_evidence_id",
    "evidence_type",
    "source_citation",
    "source_url",
    "review_status",
}
BOMBUS_OBSERVATION_REQUIRED = {
    "island_id",
    "bombus_occurrence_evidence",
    "bombus_record_count",
    "target_group_record_count",
    "target_group_spatial_units",
    "target_group_temporal_units",
    "distinct_dataset_count",
    "latest_record_year",
    "observation_diagnostic_flags",
}


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise typer.BadParameter("NEE channel-input config must be a mapping")
    if config.get("contract") != "chapter1_nee_channel_inputs_v1":
        raise typer.BadParameter("unexpected NEE channel-input contract")
    return config


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _require_columns(table: pd.DataFrame, columns: list[str] | set[str], label: str) -> None:
    missing = set(columns).difference(table.columns)
    if missing:
        raise ValueError(f"{label} missing required columns: {sorted(missing)}")


def _require_unique_pair(table: pd.DataFrame, left: str, right: str, label: str) -> None:
    if table[[left, right]].duplicated().any():
        raise ValueError(f"{label} contains duplicate {left} x {right} rows")


def _numeric_nonnegative(table: pd.DataFrame, columns: list[str], label: str) -> None:
    for column in columns:
        values = pd.to_numeric(table[column], errors="coerce")
        invalid = values.notna() & values.lt(0)
        if invalid.any():
            raise ValueError(f"{label}.{column} contains negative values")


def validate_source_availability(table: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Validate one outcome-blind source-state row per island x channel."""
    required = list(config["required_source_columns"])
    _require_columns(table, required, "source availability")
    result = table[required].copy()
    for column in required:
        result[column] = _text(result[column])

    if result["island_id"].eq("").any():
        raise ValueError("source availability contains blank island_id")
    if result["channel_id"].eq("").any():
        raise ValueError("source availability contains blank channel_id")
    _require_unique_pair(result, "island_id", "channel_id", "source availability")

    channels = set(config["channels"])
    invalid_channels = sorted(set(result["channel_id"]).difference(channels))
    if invalid_channels:
        raise ValueError(f"source availability contains unregistered channels: {invalid_channels}")

    source_states = set(config["source_states"])
    invalid_states = sorted(set(result["source_state"]).difference(source_states))
    if invalid_states:
        raise ValueError(f"source availability contains invalid source_state: {invalid_states}")

    review_values = {
        config["review_status"]["accepted"],
        config["review_status"]["pending"],
        config["review_status"]["rejected"],
    }
    invalid_reviews = sorted(set(result["review_status"]).difference(review_values))
    if invalid_reviews:
        raise ValueError(f"source availability contains invalid review_status: {invalid_reviews}")

    decisive = result["source_state"].isin({"available", "structurally_absent"})
    accepted = result["review_status"].eq(config["review_status"]["accepted"])
    if (decisive & ~accepted).any():
        raise ValueError("available/structurally_absent source states require accepted review")
    if (decisive & result["source_region_id"].eq("")).any():
        raise ValueError("decisive source states require source_region_id")
    if (decisive & result["evidence_id"].eq("")).any():
        raise ValueError("decisive source states require evidence_id")
    if (decisive & result["source_citation"].eq("")).any():
        raise ValueError("decisive source states require source_citation")
    return result


def validate_island_observation(table: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Validate one effort-aware observation row per island x channel."""
    required = list(config["required_observation_columns"])
    _require_columns(table, required, "island observation")
    result = table[required].copy()
    text_columns = [
        "island_id",
        "channel_id",
        "observation_state",
        "evidence_source",
        "quality_flags",
    ]
    for column in text_columns:
        result[column] = _text(result[column])

    if result["island_id"].eq("").any():
        raise ValueError("island observation contains blank island_id")
    if result["channel_id"].eq("").any():
        raise ValueError("island observation contains blank channel_id")
    _require_unique_pair(result, "island_id", "channel_id", "island observation")

    channels = set(config["channels"])
    invalid_channels = sorted(set(result["channel_id"]).difference(channels))
    if invalid_channels:
        raise ValueError(f"island observation contains unregistered channels: {invalid_channels}")

    states = set(config["observation_states"])
    invalid_states = sorted(set(result["observation_state"]).difference(states))
    if invalid_states:
        raise ValueError(f"island observation contains invalid observation_state: {invalid_states}")

    count_columns = [
        "channel_record_count",
        "background_record_count",
        "background_spatial_units",
        "background_temporal_units",
        "distinct_dataset_count",
    ]
    _numeric_nonnegative(result, count_columns, "island observation")
    channel_count = pd.to_numeric(result["channel_record_count"], errors="coerce").fillna(0)
    detected = result["observation_state"].eq("detected")
    if (detected & channel_count.lt(1)).any():
        raise ValueError("detected channel observations require channel_record_count >= 1")
    adequate_zero = result["observation_state"].eq("adequate_non_detection")
    if (adequate_zero & channel_count.ne(0)).any():
        raise ValueError("adequate_non_detection requires channel_record_count == 0")
    background = pd.to_numeric(result["background_record_count"], errors="coerce").fillna(0)
    if (adequate_zero & background.lt(1)).any():
        raise ValueError("adequate_non_detection requires positive background effort")
    return result


def project_channel_states(
    source: pd.DataFrame,
    observation: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Project source + observation evidence into N1 channel states."""
    source_valid = validate_source_availability(source, config)
    observation_valid = validate_island_observation(observation, config)
    merged = source_valid.merge(
        observation_valid,
        on=["island_id", "channel_id"],
        how="left",
        validate="one_to_one",
    )
    for column in OBSERVATION_COLUMNS:
        if column in {"island_id", "channel_id"}:
            continue
        if column not in merged.columns:
            merged[column] = ""
    merged["observation_state"] = _text(merged["observation_state"]).replace("", "unresolved")
    merged["quality_flags"] = _text(merged["quality_flags"])
    merged["evidence_source"] = _text(merged["evidence_source"])

    channel_states: list[str] = []
    projection_flags: list[str] = []
    for row in merged.to_dict("records"):
        source_state = str(row["source_state"])
        observation_state = str(row["observation_state"])
        flags: list[str] = []
        if source_state == "structurally_absent":
            channel_state = "structurally_absent"
            if observation_state == "detected":
                flags.append("detected_despite_structural_source_absence")
        elif source_state == "unresolved":
            channel_state = "unresolved"
            flags.append("source_state_unresolved")
        elif observation_state == "detected":
            channel_state = "retained"
        elif observation_state == "adequate_non_detection":
            channel_state = "disrupted"
        else:
            channel_state = "unresolved"
            flags.append("observation_not_evaluable")
        channel_states.append(channel_state)
        projection_flags.append("|".join(flags))

    merged["channel_state"] = channel_states
    merged["projection_flags"] = projection_flags
    projected_states = set(config["projected_states"])
    invalid = sorted(set(channel_states).difference(projected_states))
    if invalid:
        raise ValueError(f"projection generated invalid channel states: {invalid}")
    return merged[PROJECTED_COLUMNS].copy()


def qualification_receipt(projected: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Return outcome-blind support status for every frozen N0 channel."""
    _require_columns(projected, PROJECTED_COLUMNS, "projected channel state")
    pilot = config["support_tiers"]["pilot"]
    confirmatory = config["support_tiers"]["confirmatory"]
    rows: list[dict[str, Any]] = []
    for channel in config["channels"]:
        group = projected.loc[projected["channel_id"].eq(channel)].copy()
        available = group["source_state"].eq("available")
        retained = group["channel_state"].eq("retained")
        disrupted = group["channel_state"].eq("disrupted")
        evaluable = available & (retained | disrupted)
        source_regions = int(
            group.loc[available, "source_region_id"].replace("", pd.NA).dropna().nunique()
        )
        metrics = {
            "n_rows": int(len(group)),
            "n_source_available": int(available.sum()),
            "n_structurally_absent": int(group["source_state"].eq("structurally_absent").sum()),
            "n_source_unresolved": int(group["source_state"].eq("unresolved").sum()),
            "n_retained": int(retained.sum()),
            "n_disrupted": int(disrupted.sum()),
            "n_observation_unresolved": int((available & group["channel_state"].eq("unresolved")).sum()),
            "n_evaluable_source_available": int(evaluable.sum()),
            "n_source_regions": source_regions,
        }

        def passes(thresholds: dict[str, Any]) -> bool:
            return (
                metrics["n_evaluable_source_available"]
                >= int(thresholds["min_evaluable_source_available_islands"])
                and metrics["n_retained"] >= int(thresholds["min_retained"])
                and metrics["n_disrupted"] >= int(thresholds["min_disrupted"])
                and metrics["n_source_regions"] >= int(thresholds["min_source_regions"])
            )

        if passes(confirmatory):
            tier = "confirmatory"
            reason = ""
        elif passes(pilot):
            tier = "pilot"
            reason = "below_confirmatory_support"
        else:
            tier = "not_qualified"
            failed: list[str] = []
            if metrics["n_evaluable_source_available"] < int(
                pilot["min_evaluable_source_available_islands"]
            ):
                failed.append("insufficient_evaluable_islands")
            if metrics["n_retained"] < int(pilot["min_retained"]):
                failed.append("insufficient_retained")
            if metrics["n_disrupted"] < int(pilot["min_disrupted"]):
                failed.append("insufficient_disrupted")
            if metrics["n_source_regions"] < int(pilot["min_source_regions"]):
                failed.append("insufficient_source_regions")
            reason = "|".join(failed) or "data_quality_support_failure"
        rows.append(
            {
                "channel_id": channel,
                **metrics,
                "support_tier": tier,
                "N1_gate_eligible": tier == config["support_tiers"]["N1_gate_eligible_tier"],
                "exclusion_reason": reason,
            }
        )
    return pd.DataFrame(rows)


def adapt_bombus_source_availability(
    applicability: pd.DataFrame,
    source_evidence: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Adapt the existing outcome-blind Bombus applicability registry to N1 source rows."""
    _require_columns(applicability, BOMBUS_SOURCE_REQUIRED, "Bombus applicability")
    _require_columns(source_evidence, BOMBUS_EVIDENCE_REQUIRED, "Bombus source evidence")
    evidence = source_evidence[list(BOMBUS_EVIDENCE_REQUIRED)].copy()
    for column in BOMBUS_EVIDENCE_REQUIRED:
        evidence[column] = _text(evidence[column])
    if evidence["source_region_evidence_id"].duplicated().any():
        raise ValueError("Bombus source evidence contains duplicate source_region_evidence_id")

    work = applicability.copy()
    for column in BOMBUS_SOURCE_REQUIRED:
        work[column] = _text(work[column])
    work = work.merge(
        evidence,
        on="source_region_evidence_id",
        how="left",
        validate="many_to_one",
        suffixes=("", "_evidence"),
    )
    mapping = config["bombus_adapter"]["source_mapping"]
    accepted_label = config["review_status"]["accepted"]
    rows: list[dict[str, Any]] = []
    for row in work.to_dict("records"):
        reviews_accepted = (
            str(row.get("source_region_review_status", "")).lower() == "accepted"
            and str(row.get("assignment_review_status", "")).lower() == "accepted"
            and str(row.get("review_status", "")).lower() == "accepted"
        )
        mapped = str(mapping.get(str(row.get("applicability", "")), "unresolved"))
        source_state = mapped if reviews_accepted else "unresolved"
        review_status = accepted_label if reviews_accepted else config["review_status"]["pending"]
        rows.append(
            {
                "island_id": str(row.get("island_id", "")),
                "source_region_id": str(row.get("source_region_id", "")),
                "channel_id": "bombus",
                "source_state": source_state,
                "evidence_id": str(row.get("source_region_evidence_id", "")),
                "evidence_type": str(row.get("evidence_type", "")),
                "source_citation": str(row.get("source_citation", "")),
                "source_url": str(row.get("source_url", "")),
                "review_status": review_status,
            }
        )
    return validate_source_availability(pd.DataFrame(rows, columns=SOURCE_COLUMNS), config)


def adapt_bombus_observation(
    diagnostics: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Adapt the existing recency-aware Bombus evidence to generic N1 observations."""
    _require_columns(diagnostics, BOMBUS_OBSERVATION_REQUIRED, "Bombus occurrence evidence")
    work = diagnostics.copy()
    for column in BOMBUS_OBSERVATION_REQUIRED:
        work[column] = _text(work[column])
    mapping = config["bombus_adapter"]["observation_mapping"]
    rows: list[dict[str, Any]] = []
    for row in work.to_dict("records"):
        state = str(mapping.get(str(row.get("bombus_occurrence_evidence", "")), "unresolved"))
        rows.append(
            {
                "island_id": str(row.get("island_id", "")),
                "channel_id": "bombus",
                "observation_state": state,
                "channel_record_count": str(row.get("bombus_record_count", "")),
                "background_record_count": str(row.get("target_group_record_count", "")),
                "background_spatial_units": str(row.get("target_group_spatial_units", "")),
                "background_temporal_units": str(row.get("target_group_temporal_units", "")),
                "distinct_dataset_count": str(row.get("distinct_dataset_count", "")),
                "latest_background_year": str(row.get("latest_record_year", "")),
                "evidence_source": "bombus_occurrence_evidence_v1",
                "quality_flags": str(row.get("observation_diagnostic_flags", "")),
            }
        )
    return validate_island_observation(pd.DataFrame(rows, columns=OBSERVATION_COLUMNS), config)


def _write_outputs(
    output_dir: Path,
    source: pd.DataFrame,
    observation: pd.DataFrame,
    projected: pd.DataFrame,
    receipt: pd.DataFrame,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    source.to_csv(output_dir / "nee_source_availability.csv", index=False)
    observation.to_csv(output_dir / "nee_island_channel_observation.csv", index=False)
    projected.to_csv(output_dir / "nee_projected_channel_state.csv", index=False)
    receipt.to_csv(output_dir / "nee_channel_qualification_receipt.csv", index=False)
    summary = {
        "n_source_rows": int(len(source)),
        "n_observation_rows": int(len(observation)),
        "n_projected_rows": int(len(projected)),
        "channels": receipt.to_dict("records"),
    }
    (output_dir / "nee_channel_qualification_receipt.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


@app.command("validate")
def validate_command(
    source_csv: Path = typer.Option(..., exists=True),
    observation_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/chapter1_nee_channel_inputs.yml")),
) -> None:
    """Validate standardized channel inputs, project states, and write a support receipt."""
    config = load_config(config_path)
    source = validate_source_availability(pd.read_csv(source_csv, dtype=str).fillna(""), config)
    observation = validate_island_observation(
        pd.read_csv(observation_csv, dtype=str).fillna(""), config
    )
    projected = project_channel_states(source, observation, config)
    receipt = qualification_receipt(projected, config)
    _write_outputs(output_dir, source, observation, projected, receipt)
    typer.echo(receipt.to_string(index=False))


@app.command("adapt-bombus")
def adapt_bombus_command(
    applicability_csv: Path = typer.Option(..., exists=True),
    source_evidence_csv: Path = typer.Option(..., exists=True),
    occurrence_evidence_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/chapter1_nee_channel_inputs.yml")),
) -> None:
    """Project existing outcome-blind Bombus assets into the generic N1 evidence layer."""
    config = load_config(config_path)
    source = adapt_bombus_source_availability(
        pd.read_csv(applicability_csv, dtype=str).fillna(""),
        pd.read_csv(source_evidence_csv, dtype=str).fillna(""),
        config,
    )
    observation = adapt_bombus_observation(
        pd.read_csv(occurrence_evidence_csv, dtype=str).fillna(""), config
    )
    projected = project_channel_states(source, observation, config)
    receipt = qualification_receipt(projected, config)
    _write_outputs(output_dir, source, observation, projected, receipt)
    bombus = receipt.loc[receipt["channel_id"].eq("bombus")].iloc[0].to_dict()
    typer.echo(json.dumps(bombus, indent=2))


if __name__ == "__main__":
    app()
