"""H3 sensitivities for the mainland-zero concern in corrected GloPL geography.

Tests whether the continuous pollen-limitation gradient persists among positive-distance
(offshore/island) measurement cells, whether a binary mainland/offshore contrast suffices,
and whether the positive-distance result is robust to deleting one publication at a time.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer

from island_v2.chapter1_h5_glopl_global_distance import (
    _clustered_wls,
    _context_dummies,
    _measurement_dummies,
    _p2,
    _standardize_distance,
    build_global_design,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _binary_design(frame: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    context_cols, context_names, _ = _context_dummies(frame)
    measure_cols, measure_names = _measurement_dummies(frame)
    offshore = (
        pd.to_numeric(frame["log1p_distance_to_major_continent_km"], errors="coerce")
        .gt(0)
        .to_numpy(float)
    )
    columns = [np.ones(len(frame)), *context_cols, offshore, *measure_cols]
    names = ["intercept", *context_names, "offshore_indicator", *measure_names]
    return np.column_stack(columns), names


def _fit_distance(frame: pd.DataFrame) -> dict[str, Any]:
    work, _, _ = _standardize_distance(frame)
    X, names = build_global_design(work)
    fit = _clustered_wls(work, X, names)
    if not fit.get("evaluable"):
        return fit
    idx = fit["names"].index("z_distance")
    est = float(fit["beta"][idx])
    se = float(fit["se"][idx])
    z = est / se if se > 0 else float("nan")
    return {
        "evaluable": True,
        "estimate": est,
        "se": se,
        "two_sided_p": _p2(z),
        "n_cells": fit["n_cells"],
        "n_publications": fit["n_publications"],
        "n_sites": fit["n_sites"],
    }


def _fit_binary(frame: pd.DataFrame) -> dict[str, Any]:
    X, names = _binary_design(frame)
    fit = _clustered_wls(frame, X, names)
    if not fit.get("evaluable"):
        return fit
    idx = fit["names"].index("offshore_indicator")
    est = float(fit["beta"][idx])
    se = float(fit["se"][idx])
    z = est / se if se > 0 else float("nan")
    return {
        "evaluable": True,
        "estimate": est,
        "se": se,
        "two_sided_p": _p2(z),
        "n_cells": fit["n_cells"],
        "n_publications": fit["n_publications"],
        "n_sites": fit["n_sites"],
    }


def run_sensitivity(cells: pd.DataFrame) -> tuple[dict[str, Any], pd.DataFrame]:
    required = {
        "study_key",
        "site_key",
        "analysis_regime",
        "log1p_distance_to_major_continent_km",
        "PL_Effect_Size",
        "analysis_weight",
        "PL_Effect_Size_Type1",
        "PL_Effect_Size_Type2",
        "Constant_added",
        "Level_of_Supplementation",
    }
    if missing := required - set(cells.columns):
        raise typer.BadParameter(f"H3 cells missing columns: {sorted(missing)}")
    frame = cells.copy()
    dist = pd.to_numeric(frame["log1p_distance_to_major_continent_km"], errors="coerce")
    frame = frame.loc[dist.notna()].copy()
    offshore = frame.loc[
        pd.to_numeric(frame["log1p_distance_to_major_continent_km"], errors="coerce").gt(0)
    ].copy()

    island_result = _fit_distance(offshore)
    binary_result = _fit_binary(frame)

    jackknife_rows: list[dict[str, Any]] = []
    for study in sorted(offshore["study_key"].astype(str).unique()):
        part = offshore.loc[offshore["study_key"].astype(str).ne(study)].copy()
        result = _fit_distance(part)
        jackknife_rows.append(
            {
                "deleted_publication": study,
                "evaluable": bool(result.get("evaluable", False)),
                "estimate": result.get("estimate"),
                "se": result.get("se"),
                "two_sided_p": result.get("two_sided_p"),
                "n_cells": result.get("n_cells"),
                "n_publications": result.get("n_publications"),
                "n_sites": result.get("n_sites"),
            }
        )
    jackknife = pd.DataFrame(jackknife_rows)
    ok = jackknife.loc[jackknife["evaluable"].astype(bool)].copy()
    summary = {
        "contract": "chapter1_h3_offshore_gradient_sensitivity_v1",
        "offshore_definition": "corrected log1p distance > 0",
        "offshore_continuous_gradient": island_result,
        "mainland_vs_offshore_indicator": binary_result,
        "jackknife": {
            "n_deletions": int(len(ok)),
            "all_estimates_positive": bool(pd.to_numeric(ok["estimate"], errors="coerce").gt(0).all()),
            "worst_two_sided_p": float(pd.to_numeric(ok["two_sided_p"], errors="coerce").max()),
            "min_estimate": float(pd.to_numeric(ok["estimate"], errors="coerce").min()),
            "max_estimate": float(pd.to_numeric(ok["estimate"], errors="coerce").max()),
        },
        "interpretation": (
            "The offshore-only model tests whether the global H3 result is more than a mainland/offshore "
            "step. The binary model asks whether the step alone explains the pattern. These are post-hoc "
            "robustness analyses and do not establish causation."
        ),
    }
    return summary, jackknife


@app.command("run")
def run(
    measurement_cells_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    cells = pd.read_csv(measurement_cells_csv)
    summary, jackknife = run_sensitivity(cells)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "h3_offshore_gradient_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    jackknife.to_csv(output_dir / "h3_offshore_leave_one_publication.csv", index=False)
    typer.echo(json.dumps(summary, indent=2))


if __name__ == "__main__":
    app()
# rerun-marker: 2026-10-03 latest-head persistence retry
# rerun-marker: 2026-10-03 after-h2-persist
