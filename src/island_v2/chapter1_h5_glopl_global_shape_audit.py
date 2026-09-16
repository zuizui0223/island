"""Post-hoc shape diagnostic for the supported all-site GloPL distance association.

This module never changes the parent model, p-value, or promotion decision. It reuses the
frozen parent publication x site x measurement cells and their weights, decomposing the
continuous distance association into a mainland-to-average-offshore level shift and a
within-offshore distance slope.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_h5_glopl_global_distance import (
    MEASUREMENT_COLUMNS,
    _clustered_wls,
    _context_dummies,
    _measurement_dummies,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)
CONTRACT = "chapter1_h5_glopl_global_shape_audit_v1"


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected GloPL global shape-audit contract")
    return config


def _two_sided_p(z_value: float) -> float:
    if not math.isfinite(z_value):
        return float("nan")
    return math.erfc(abs(float(z_value)) / math.sqrt(2.0))


def prepare_shape_cells(
    cells: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    required = {
        "study_key",
        "site_key",
        "analysis_regime",
        "log1p_distance_to_major_continent_km",
        "PL_Effect_Size",
        "analysis_weight",
        *MEASUREMENT_COLUMNS,
    }
    if missing := required - set(cells.columns):
        raise typer.BadParameter(f"shape audit parent cells missing columns: {sorted(missing)}")
    work = cells.copy().reset_index(drop=True)
    distance = pd.to_numeric(
        work["log1p_distance_to_major_continent_km"], errors="coerce"
    ).to_numpy(float)
    outcome = pd.to_numeric(work["PL_Effect_Size"], errors="coerce").to_numpy(float)
    weights = pd.to_numeric(work["analysis_weight"], errors="coerce").to_numpy(float)
    if (
        not np.isfinite(distance).all()
        or not np.isfinite(outcome).all()
        or not np.isfinite(weights).all()
        or np.any(weights <= 0)
    ):
        raise typer.BadParameter("nonfinite parent input in GloPL shape audit")
    if np.any(distance < 0):
        raise typer.BadParameter("negative log-distance in GloPL shape audit")

    offshore = distance > 0.0
    if int(offshore.sum()) < 2 or int((~offshore).sum()) < 2:
        raise typer.BadParameter("shape audit requires both mainland and offshore cells")
    positive = distance[offshore]
    mean = float(np.mean(positive))
    sd = float(np.std(positive, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise typer.BadParameter("offshore distance has zero or invalid variance")

    work["offshore"] = offshore.astype(float)
    z = np.zeros(len(work), dtype=float)
    z[offshore] = (positive - mean) / sd
    work["offshore_distance_z"] = z
    work["PL_Effect_Size"] = outcome
    work["analysis_weight"] = weights

    audit = {
        "n_cells": int(len(work)),
        "n_mainland_cells": int((~offshore).sum()),
        "n_offshore_cells": int(offshore.sum()),
        "n_mainland_sites": int(work.loc[~offshore, "site_key"].nunique()),
        "n_offshore_sites": int(work.loc[offshore, "site_key"].nunique()),
        "n_publications": int(work["study_key"].nunique()),
        "n_offshore_publications": int(work.loc[offshore, "study_key"].nunique()),
        "offshore_log_distance_mean": mean,
        "offshore_log_distance_sd": sd,
        "parent_weights_reused_without_renormalization": True,
    }
    return work, audit


def build_shape_design(frame: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    context_cols, context_names, _ = _context_dummies(frame)
    measure_cols, measure_names = _measurement_dummies(frame)
    columns: list[np.ndarray] = [
        np.ones(len(frame), dtype=float),
        *context_cols,
        pd.to_numeric(frame["offshore"], errors="coerce").to_numpy(float),
        pd.to_numeric(frame["offshore_distance_z"], errors="coerce").to_numpy(float),
        *measure_cols,
    ]
    names = [
        "intercept",
        *context_names,
        "offshore_indicator",
        "offshore_distance_z",
        *measure_names,
    ]
    return np.column_stack(columns), names


def build_offshore_only_design(frame: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    context_cols, context_names, _ = _context_dummies(frame)
    measure_cols, measure_names = _measurement_dummies(frame)
    columns: list[np.ndarray] = [
        np.ones(len(frame), dtype=float),
        *context_cols,
        pd.to_numeric(frame["offshore_distance_z"], errors="coerce").to_numpy(float),
        *measure_cols,
    ]
    names = ["intercept", *context_names, "offshore_distance_z", *measure_names]
    return np.column_stack(columns), names


def fit_two_part_shape(frame: pd.DataFrame) -> dict[str, Any]:
    design, names = build_shape_design(frame)
    fit = _clustered_wls(frame, design, names)
    if not fit.get("evaluable"):
        return fit
    index = {name: i for i, name in enumerate(fit["names"])}
    step_i = index["offshore_indicator"]
    slope_i = index["offshore_distance_z"]
    step = float(fit["beta"][step_i])
    step_se = float(fit["se"][step_i])
    slope = float(fit["beta"][slope_i])
    slope_se = float(fit["se"][slope_i])
    step_z = step / step_se if step_se > 0 else float("nan")
    slope_z = slope / slope_se if slope_se > 0 else float("nan")
    return {
        "evaluable": True,
        "offshore_step": step,
        "offshore_step_se": step_se,
        "offshore_step_two_sided_p": _two_sided_p(step_z),
        "within_offshore_slope": slope,
        "within_offshore_slope_se": slope_se,
        "within_offshore_slope_two_sided_p": _two_sided_p(slope_z),
        "n_cells": int(fit["n_cells"]),
        "n_sites": int(fit["n_sites"]),
        "n_publications": int(fit["n_publications"]),
    }


def fit_offshore_only(frame: pd.DataFrame) -> dict[str, Any]:
    offshore = frame.loc[frame["offshore"].eq(1.0)].copy()
    design, names = build_offshore_only_design(offshore)
    fit = _clustered_wls(offshore, design, names)
    if not fit.get("evaluable"):
        return fit
    idx = fit["names"].index("offshore_distance_z")
    slope = float(fit["beta"][idx])
    se = float(fit["se"][idx])
    z = slope / se if se > 0 else float("nan")
    return {
        "evaluable": True,
        "distance_slope": slope,
        "distance_slope_se": se,
        "distance_slope_two_sided_p": _two_sided_p(z),
        "n_cells": int(fit["n_cells"]),
        "n_sites": int(fit["n_sites"]),
        "n_publications": int(fit["n_publications"]),
    }


def classify_shape_audit(
    two_part: dict[str, Any],
    offshore_only: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    alpha = 0.05
    if not two_part.get("evaluable") or not offshore_only.get("evaluable"):
        label = "mixed_or_unresolved"
    else:
        step_supported = float(two_part.get("offshore_step_two_sided_p", 1.0)) <= alpha
        two_part_slope_supported = float(
            two_part.get("within_offshore_slope_two_sided_p", 1.0)
        ) <= alpha
        offshore_only_supported = float(
            offshore_only.get("distance_slope_two_sided_p", 1.0)
        ) <= alpha
        positive_both = (
            float(two_part.get("within_offshore_slope", float("nan"))) > 0
            and float(offshore_only.get("distance_slope", float("nan"))) > 0
        )
        if step_supported and not two_part_slope_supported and not offshore_only_supported:
            label = "step_dominated"
        elif positive_both and two_part_slope_supported and offshore_only_supported:
            label = "within_offshore_gradient_present"
        else:
            label = "mixed_or_unresolved"
    return {
        "shape_classification": label,
        "can_change_parent_promotion": False,
        "role": str(config["inferential_role"]),
    }


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, allow_nan=True) + "\n", encoding="utf-8")


@app.command()
def run(
    cells_csv: Path = typer.Option(..., exists=True),
    parent_result_json: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    parent = json.loads(parent_result_json.read_text(encoding="utf-8"))
    expected = config["parent"]
    if int(expected["workflow_run_id"]) != 35090599662:
        raise typer.BadParameter("unexpected frozen parent run")
    observed_parent_slope = float(parent["global_gradient"]["distance_slope"])
    expected_parent_slope = float(expected["parent_global_slope"])
    if not math.isclose(observed_parent_slope, expected_parent_slope, rel_tol=0.0, abs_tol=1e-12):
        raise typer.BadParameter("parent global slope does not match frozen shape-audit parent")
    if str(parent["decision"]["classification"]) != str(expected["parent_classification"]):
        raise typer.BadParameter("parent classification does not match frozen shape-audit parent")

    cells = pd.read_csv(cells_csv)
    prepared, support = prepare_shape_cells(cells, config)
    two_part = fit_two_part_shape(prepared)
    offshore_only = fit_offshore_only(prepared)
    decision = classify_shape_audit(two_part, offshore_only, config)
    result = {
        "contract": CONTRACT,
        "status": "completed_posthoc_shape_audit",
        "parent": {
            "workflow_run_id": int(expected["workflow_run_id"]),
            "artifact_id": int(expected["artifact_id"]),
            "parent_global_slope": observed_parent_slope,
            "parent_one_sided_p": float(expected["parent_one_sided_p"]),
            "parent_classification": str(parent["decision"]["classification"]),
            "parent_promotion_immutable": True,
        },
        "support": support,
        "two_part_global": two_part,
        "offshore_only": offshore_only,
        "decision": decision,
        "claim_ceiling": config["claim_ceiling"],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    prepared.to_csv(output_dir / "SHAPE_AUDIT_CELLS.csv.gz", index=False, compression="gzip")
    _write_json(output_dir / "RESULT.json", result)
    lines = [
        "# GloPL global distance shape audit",
        "",
        f"- parent global slope: {observed_parent_slope}",
        f"- mainland cells: {support['n_mainland_cells']}",
        f"- offshore cells: {support['n_offshore_cells']}",
        f"- mainland-to-average-offshore step: {two_part.get('offshore_step')}",
        f"- step two-sided p: {two_part.get('offshore_step_two_sided_p')}",
        f"- two-part within-offshore slope: {two_part.get('within_offshore_slope')}",
        f"- within-offshore p: {two_part.get('within_offshore_slope_two_sided_p')}",
        f"- offshore-only slope: {offshore_only.get('distance_slope')}",
        f"- offshore-only p: {offshore_only.get('distance_slope_two_sided_p')}",
        f"- shape classification: {decision['shape_classification']}",
        "- parent promotion changed: false",
    ]
    (output_dir / "RESULT_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.echo(json.dumps(result, indent=2))


if __name__ == "__main__":
    app()
