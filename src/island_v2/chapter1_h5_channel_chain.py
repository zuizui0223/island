"""Fail-closed H5 evidence chain for independent pollination-channel measurements.

This module keeps source availability, island retention/disruption, visitation,
single-visit effectiveness, and effective service as separate evidence layers.
It does not infer channel state from floral phenotype or climate compatibility.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

SOURCE_STATES = {"available", "structurally_absent", "unresolved"}
RETENTION_STATES = {"retained", "disrupted", "structurally_absent", "unresolved"}
YES_NO = {"yes", "no"}
KEY = ["system_id", "island_id", "channel_id"]


def _required_columns(config: dict[str, Any]) -> list[str]:
    return [str(value) for value in config["required_columns"]]


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _number(series: pd.Series, *, label: str) -> pd.Series:
    raw = _text(series)
    numeric = pd.to_numeric(raw, errors="coerce")
    invalid = raw.ne("") & numeric.isna()
    if invalid.any():
        raise ValueError(f"{label} contains non-numeric values")
    return numeric


def validate_channel_panel(frame: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Validate one independent channel-evidence row per system x island x channel."""
    required = _required_columns(config)
    missing = set(required).difference(frame.columns)
    if missing:
        raise ValueError(f"H5 channel panel missing columns: {sorted(missing)}")

    work = frame.copy()
    for column in required:
        work[column] = _text(work[column])

    if work[KEY].eq("").any(axis=None):
        raise ValueError("H5 channel panel contains blank unit keys")
    if work.duplicated(KEY).any():
        raise ValueError("H5 channel panel must have one row per system x island x channel")

    source = work["source_channel_state"]
    retention = work["island_retention_state"]
    if not source.isin(SOURCE_STATES).all():
        raise ValueError("invalid source_channel_state")
    if not retention.isin(RETENTION_STATES).all():
        raise ValueError("invalid island_retention_state")

    structural = source.eq("structurally_absent")
    if (structural & ~retention.eq("structurally_absent")).any():
        raise ValueError("structurally absent source channel cannot be labelled as island loss")
    if ((~structural) & retention.eq("structurally_absent")).any():
        raise ValueError("island structurally_absent requires source structurally_absent")

    resolved_source = source.isin({"available", "structurally_absent"})
    if (resolved_source & work["source_evidence_id"].eq("")).any():
        raise ValueError("resolved source state requires source_evidence_id")

    resolved_retention = retention.isin({"retained", "disrupted"})
    if (resolved_retention & work["retention_evidence_type"].eq("")).any():
        raise ValueError("retained/disrupted state requires retention_evidence_type")
    if (resolved_retention & work["retention_evidence_id"].eq("")).any():
        raise ValueError("retained/disrupted state requires retention_evidence_id")

    controls = work["no_visit_control_present"]
    if not controls.isin(YES_NO).all():
        raise ValueError("no_visit_control_present must be yes or no")

    effort = _number(work["observation_effort_flower_hours"], label="observation effort")
    bouts = _number(work["visit_bouts"], label="visit_bouts")
    if (effort.dropna() < 0).any():
        raise ValueError("observation effort must be non-negative")
    if (bouts.dropna() < 0).any():
        raise ValueError("visit_bouts must be non-negative")
    non_integer_bouts = bouts.dropna().map(lambda value: float(value).is_integer()).eq(False)
    if non_integer_bouts.any():
        raise ValueError("visit_bouts must be integer-valued")

    _number(
        work["mean_single_visit_conspecific_pollen"],
        label="mean_single_visit_conspecific_pollen",
    )
    _number(
        work["mean_no_visit_conspecific_pollen"],
        label="mean_no_visit_conspecific_pollen",
    )
    return work


def evaluate_channel_chain(frame: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Evaluate evidence gates and compute realized effective pollen service."""
    work = validate_channel_panel(frame, config)
    effort = _number(work["observation_effort_flower_hours"], label="observation effort")
    bouts = _number(work["visit_bouts"], label="visit_bouts")
    single = _number(
        work["mean_single_visit_conspecific_pollen"],
        label="mean_single_visit_conspecific_pollen",
    )
    control = _number(
        work["mean_no_visit_conspecific_pollen"],
        label="mean_no_visit_conspecific_pollen",
    )

    source_gate = (
        work["source_channel_state"].isin({"available", "structurally_absent"})
        & work["source_evidence_id"].ne("")
    )
    loss_eligible = work["source_channel_state"].eq("available")
    retention_gate = (
        work["island_retention_state"].isin({"retained", "disrupted"})
        & work["retention_evidence_type"].ne("")
        & work["retention_evidence_id"].ne("")
        & loss_eligible
    )
    visitation_gate = effort.gt(0) & bouts.notna()
    visit_rate = bouts / effort.where(effort.gt(0))

    positive_visits = visitation_gate & bouts.gt(0)
    zero_visits = visitation_gate & bouts.eq(0)
    svd_gate = (
        positive_visits
        & work["no_visit_control_present"].eq("yes")
        & single.notna()
        & control.notna()
    )
    adjusted_svd = single - control
    adjusted_svd = adjusted_svd.where(svd_gate)

    service = pd.Series(pd.NA, index=work.index, dtype="Float64")
    service.loc[zero_visits] = 0.0
    service.loc[svd_gate] = (
        visit_rate.loc[svd_gate].astype(float) * adjusted_svd.loc[svd_gate].astype(float)
    )
    service_gate = visitation_gate & (zero_visits | svd_gate)
    full_observed_chain = source_gate & retention_gate & visitation_gate & svd_gate & service_gate

    work["source_gate_pass"] = source_gate
    work["loss_contrast_eligible_unit"] = loss_eligible & retention_gate
    work["retention_gate_pass"] = retention_gate
    work["visitation_gate_pass"] = visitation_gate
    work["visit_bouts_per_flower_hour"] = visit_rate
    work["single_visit_effectiveness_gate_pass"] = svd_gate
    work["mean_background_adjusted_svd"] = adjusted_svd
    work["single_visit_effectiveness_status"] = "missing_or_not_evaluable"
    work.loc[zero_visits, "single_visit_effectiveness_status"] = (
        "not_applicable_zero_visitation"
    )
    work.loc[svd_gate, "single_visit_effectiveness_status"] = "measured_with_control"
    work["effective_service_gate_pass"] = service_gate
    work["effective_pollen_delivery_per_flower_hour"] = service
    work["effective_service_basis"] = "withheld"
    work.loc[zero_visits, "effective_service_basis"] = "adequate_zero_visitation"
    work.loc[svd_gate, "effective_service_basis"] = "visit_rate_x_background_adjusted_svd"
    work["full_positive_visit_chain_observed"] = full_observed_chain
    work["structural_absence_not_loss"] = work["source_channel_state"].eq(
        "structurally_absent"
    )
    return work


def contrast_readiness(audit: pd.DataFrame) -> pd.DataFrame:
    """Report whether a channel has retained/disrupted service contrasts without pooling guilds."""
    rows: list[dict[str, Any]] = []
    for (system_id, channel_id), group in audit.groupby(["system_id", "channel_id"], sort=True):
        usable = group.loc[
            group["loss_contrast_eligible_unit"].astype(bool)
            & group["effective_service_gate_pass"].astype(bool)
        ]
        retained = usable["island_retention_state"].eq("retained").sum()
        disrupted = usable["island_retention_state"].eq("disrupted").sum()
        positive_full = group["full_positive_visit_chain_observed"].astype(bool).sum()
        rows.append(
            {
                "system_id": system_id,
                "channel_id": channel_id,
                "n_service_estimable_units": int(len(usable)),
                "n_retained_service_units": int(retained),
                "n_disrupted_service_units": int(disrupted),
                "n_full_positive_visit_chain_units": int(positive_full),
                "retained_disrupted_service_contrast_ready": bool(
                    retained > 0 and disrupted > 0 and positive_full > 0
                ),
            }
        )
    return pd.DataFrame(rows)


def build_outputs(
    frame: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    audit = evaluate_channel_chain(frame, config)
    readiness = contrast_readiness(audit)
    ready = readiness["retained_disrupted_service_contrast_ready"].astype(bool)
    summary = {
        "contract": config["contract"],
        "n_units": int(len(audit)),
        "n_source_resolved": int(audit["source_gate_pass"].astype(bool).sum()),
        "n_retention_resolved": int(audit["retention_gate_pass"].astype(bool).sum()),
        "n_visitation_estimable": int(audit["visitation_gate_pass"].astype(bool).sum()),
        "n_single_visit_effectiveness_measured": int(
            audit["single_visit_effectiveness_gate_pass"].astype(bool).sum()
        ),
        "n_effective_service_estimable": int(
            audit["effective_service_gate_pass"].astype(bool).sum()
        ),
        "n_full_positive_visit_chain": int(
            audit["full_positive_visit_chain_observed"].astype(bool).sum()
        ),
        "n_channel_contrasts_ready": int(ready.sum()),
        "H5_channel_side_status": "contrast_ready" if ready.any() else "not_evaluable",
        "H5_full_status": "requires_post_H2_H4_plant_residual_join",
        "zero_visit_policy": (
            "adequate zero visitation identifies zero realized service but not zero per-visit "
            "effectiveness"
        ),
    }
    return audit, readiness, summary


@app.command("run")
def run(
    channel_panel_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/chapter1_h5_channel_chain.yml")),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    frame = pd.read_csv(channel_panel_csv, dtype=str).fillna("")
    audit, readiness, summary = build_outputs(frame, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output_dir / "h5_channel_chain_audit.csv", index=False)
    service = audit.loc[audit["effective_service_gate_pass"].astype(bool)].copy()
    service.to_csv(output_dir / "h5_effective_service_panel.csv", index=False)
    readiness.to_csv(output_dir / "h5_channel_contrast_readiness.csv", index=False)
    (output_dir / "h5_channel_chain_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary))


if __name__ == "__main__":
    app()
