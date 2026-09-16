"""Post-hoc influence audit for the frozen Chapter 1 GloPL H5 result.

This module cannot change the parent inference. It only asks whether the positive but
unsupported northern distance coefficient, or the near-zero context interaction, is
driven by one island or one publication. It reuses the already standardized and
publication-island aggregated cells from the parent artifact and never searches subsets
for significance.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.chapter1_h5_glopl_pollen_limitation import _fit_clustered_wls

app = typer.Typer(add_completion=False, no_args_is_help=True)
CONTRACT = "chapter1_h5_glopl_influence_audit_v1"


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected GloPL influence-audit contract")
    return config


def recompute_publication_weights(cells: pd.DataFrame) -> pd.DataFrame:
    required = {"study_key", "island_id"}
    if missing := required - set(cells.columns):
        raise typer.BadParameter(f"influence cells missing columns: {sorted(missing)}")
    out = cells.copy()
    if out.empty:
        out["analysis_weight"] = pd.Series(dtype=float)
        return out
    n_cells = out.groupby("study_key")["island_id"].transform("size").astype(float)
    out["analysis_weight"] = 1.0 / n_cells
    return out


def _run_leave_one_out(
    cells: pd.DataFrame,
    *,
    unit_column: str,
    model_config: dict[str, Any],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for unit in sorted(cells[unit_column].astype(str).unique()):
        kept = cells.loc[~cells[unit_column].astype(str).eq(unit)].copy()
        kept = recompute_publication_weights(kept)
        result = _fit_clustered_wls(
            kept,
            kept["analysis_weight"].to_numpy(float),
            model_config,
        )
        rows.append(
            {
                "removed_unit": unit,
                "evaluable": bool(result.get("evaluable")),
                "reason": str(result.get("reason", "")),
                "northern_distance_slope": result.get("northern_distance_slope"),
                "interaction_tropical_minus_northern": result.get(
                    "interaction_tropical_minus_northern"
                ),
                "n_cells": int(len(kept)),
                "n_islands": int(kept["island_id"].nunique()),
                "n_publications": int(kept["study_key"].nunique()),
            }
        )
    return pd.DataFrame(rows)


def summarize_influence(
    refits: pd.DataFrame,
    *,
    full_north: float,
    full_interaction: float,
) -> dict[str, Any]:
    required = {
        "removed_unit",
        "evaluable",
        "northern_distance_slope",
        "interaction_tropical_minus_northern",
    }
    if missing := required - set(refits.columns):
        raise typer.BadParameter(f"influence refits missing columns: {sorted(missing)}")
    n_total = int(len(refits))
    evaluable_mask = refits["evaluable"].fillna(False).astype(bool)
    valid = refits.loc[evaluable_mask].copy()
    n_eval = int(len(valid))
    base: dict[str, Any] = {
        "n_refits": n_total,
        "n_evaluable": n_eval,
        "evaluable_fraction": float(n_eval / n_total) if n_total else 0.0,
    }
    if valid.empty:
        return {
            **base,
            "north_slope_min": None,
            "north_slope_max": None,
            "north_slope_sign_positive_fraction": 0.0,
            "max_absolute_north_slope_delta_from_full": None,
            "interaction_min": None,
            "interaction_max": None,
            "interaction_sign_negative_fraction": 0.0,
            "max_absolute_interaction_delta_from_full": None,
            "most_influential_removed_unit_for_north_slope": None,
            "most_influential_removed_unit_for_interaction": None,
        }
    valid["northern_distance_slope"] = pd.to_numeric(
        valid["northern_distance_slope"], errors="coerce"
    )
    valid["interaction_tropical_minus_northern"] = pd.to_numeric(
        valid["interaction_tropical_minus_northern"], errors="coerce"
    )
    valid = valid.dropna(
        subset=["northern_distance_slope", "interaction_tropical_minus_northern"]
    )
    if valid.empty:
        return {**base, "evaluable_fraction": 0.0}
    north_delta = (valid["northern_distance_slope"] - float(full_north)).abs()
    interaction_delta = (
        valid["interaction_tropical_minus_northern"] - float(full_interaction)
    ).abs()
    north_idx = north_delta.idxmax()
    interaction_idx = interaction_delta.idxmax()
    return {
        **base,
        "north_slope_min": float(valid["northern_distance_slope"].min()),
        "north_slope_max": float(valid["northern_distance_slope"].max()),
        "north_slope_sign_positive_fraction": float(
            valid["northern_distance_slope"].gt(0).mean()
        ),
        "max_absolute_north_slope_delta_from_full": float(north_delta.max()),
        "interaction_min": float(valid["interaction_tropical_minus_northern"].min()),
        "interaction_max": float(valid["interaction_tropical_minus_northern"].max()),
        "interaction_sign_negative_fraction": float(
            valid["interaction_tropical_minus_northern"].lt(0).mean()
        ),
        "max_absolute_interaction_delta_from_full": float(interaction_delta.max()),
        "most_influential_removed_unit_for_north_slope": str(
            valid.loc[north_idx, "removed_unit"]
        ),
        "most_influential_removed_unit_for_interaction": str(
            valid.loc[interaction_idx, "removed_unit"]
        ),
    }


def classify_influence(
    island_summary: dict[str, Any],
    publication_summary: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    required_fraction = float(
        config["leave_one_out"]["minimum_evaluable_fraction"]
    )
    north_diffuse = bool(
        float(island_summary.get("evaluable_fraction", 0.0)) >= required_fraction
        and float(publication_summary.get("evaluable_fraction", 0.0))
        >= required_fraction
        and math.isclose(
            float(island_summary.get("north_slope_sign_positive_fraction", 0.0)),
            1.0,
        )
        and math.isclose(
            float(publication_summary.get("north_slope_sign_positive_fraction", 0.0)),
            1.0,
        )
    )
    interaction_diffuse = bool(
        float(island_summary.get("evaluable_fraction", 0.0)) >= required_fraction
        and float(publication_summary.get("evaluable_fraction", 0.0))
        >= required_fraction
        and math.isclose(
            float(island_summary.get("interaction_sign_negative_fraction", 0.0)),
            1.0,
        )
        and math.isclose(
            float(publication_summary.get("interaction_sign_negative_fraction", 0.0)),
            1.0,
        )
    )
    return {
        "northern_direction": (
            "northern_direction_diffuse"
            if north_diffuse
            else "northern_direction_influence_fragile"
        ),
        "context_interaction": (
            "context_interaction_sign_diffuse"
            if interaction_diffuse
            else "context_interaction_sign_fragile"
        ),
        "parent_classification": "pollen_limitation_gradient_not_supported",
        "can_promote_parent_result": False,
        "audit_role": "posthoc_influence_only",
    }


def _validate_parent_baseline(
    cells: pd.DataFrame,
    parent_result: dict[str, Any],
    model_config: dict[str, Any],
    audit_config: dict[str, Any],
) -> dict[str, Any]:
    weighted = recompute_publication_weights(cells)
    refit = _fit_clustered_wls(
        weighted,
        weighted["analysis_weight"].to_numpy(float),
        model_config,
    )
    if not refit.get("evaluable"):
        raise typer.BadParameter("parent primary cells are no longer evaluable")
    expected_north = float(audit_config["parent_result"]["frozen_primary_north_slope"])
    expected_interaction = float(
        audit_config["parent_result"]["frozen_primary_interaction"]
    )
    parent_primary = parent_result.get("primary", {})
    for observed, expected, label in (
        (refit["northern_distance_slope"], expected_north, "north slope"),
        (
            refit["interaction_tropical_minus_northern"],
            expected_interaction,
            "interaction",
        ),
        (
            parent_primary.get("northern_distance_slope"),
            expected_north,
            "parent JSON north slope",
        ),
        (
            parent_primary.get("interaction_tropical_minus_northern"),
            expected_interaction,
            "parent JSON interaction",
        ),
    ):
        if observed is None or not math.isclose(
            float(observed), float(expected), rel_tol=1e-10, abs_tol=1e-10
        ):
            raise typer.BadParameter(
                f"frozen parent baseline mismatch for {label}: {observed} vs {expected}"
            )
    return refit


@app.command("run")
def run(
    cells_csv: Path = typer.Option(..., exists=True),
    parent_result_json: Path = typer.Option(..., exists=True),
    parent_model_config: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    model_config = yaml.safe_load(parent_model_config.read_text(encoding="utf-8"))
    cells = pd.read_csv(cells_csv, dtype={"study_key": str, "island_id": str})
    parent = json.loads(parent_result_json.read_text(encoding="utf-8"))
    baseline = _validate_parent_baseline(cells, parent, model_config, config)

    island_refits = _run_leave_one_out(
        cells,
        unit_column="island_id",
        model_config=model_config,
    )
    publication_refits = _run_leave_one_out(
        cells,
        unit_column="study_key",
        model_config=model_config,
    )
    full_north = float(baseline["northern_distance_slope"])
    full_interaction = float(baseline["interaction_tropical_minus_northern"])
    island_summary = summarize_influence(
        island_refits,
        full_north=full_north,
        full_interaction=full_interaction,
    )
    publication_summary = summarize_influence(
        publication_refits,
        full_north=full_north,
        full_interaction=full_interaction,
    )
    decision = classify_influence(island_summary, publication_summary, config)
    result = {
        "contract": CONTRACT,
        "status": "completed_posthoc_influence_audit",
        "parent": {
            "workflow_run_id": config["parent_result"]["workflow_run_id"],
            "artifact_id": config["parent_result"]["artifact_id"],
            "full_northern_distance_slope": full_north,
            "full_interaction_tropical_minus_northern": full_interaction,
            "classification": "pollen_limitation_gradient_not_supported",
        },
        "leave_one_island": island_summary,
        "leave_one_publication": publication_summary,
        "decision": decision,
        "claim_ceiling": config["claim_ceiling"],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    island_refits.to_csv(output_dir / "leave_one_island_refits.csv", index=False)
    publication_refits.to_csv(
        output_dir / "leave_one_publication_refits.csv", index=False
    )
    (output_dir / "RESULT.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    lines = [
        "# H5 GloPL post-hoc influence audit",
        "",
        f"- parent North slope: {full_north:.8f}",
        f"- parent interaction: {full_interaction:.8f}",
        f"- leave-one-island North range: {island_summary.get('north_slope_min')} to {island_summary.get('north_slope_max')}",
        f"- leave-one-island positive fraction: {island_summary.get('north_slope_sign_positive_fraction')}",
        f"- leave-one-publication North range: {publication_summary.get('north_slope_min')} to {publication_summary.get('north_slope_max')}",
        f"- leave-one-publication positive fraction: {publication_summary.get('north_slope_sign_positive_fraction')}",
        f"- island interaction negative fraction: {island_summary.get('interaction_sign_negative_fraction')}",
        f"- publication interaction negative fraction: {publication_summary.get('interaction_sign_negative_fraction')}",
        "",
        f"**North direction:** {decision['northern_direction']}",
        f"**Context interaction:** {decision['context_interaction']}",
        "",
        "The parent result remains `pollen_limitation_gradient_not_supported`; this audit cannot promote it.",
    ]
    (output_dir / "RESULT_SUMMARY.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(result, indent=2))


if __name__ == "__main__":
    app()
