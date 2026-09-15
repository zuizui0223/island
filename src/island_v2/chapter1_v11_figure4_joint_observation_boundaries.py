"""Render Chapter 1 v11 Figure 4 from frozen P3 and mechanism-boundary artifacts.

Panels A/B show the jointly frozen V5×V6 observation-bias surface for the
primary direct-only native-nonendemic Palearctic and Tropical response vectors.
Each heatmap cell is the fraction of the predeclared C0×OR_C subgrid that remains
FDR-supported at the corresponding OR_R×OR_D combination. These fractions are
geometry of an assumption grid, never probabilities.

Panel C shows the partial-identification envelope for accessibility slopes plus
vector-support labels. Panel D preserves the pre-existing response-geometry and
H5 mechanism-identifiability boundaries. No biological model is fitted here.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import typer  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

app = typer.Typer(add_completion=False, no_args_is_help=True)
CONTRACT = "chapter1_v11_figure4_joint_observation_boundaries_v1"


class FigureInputError(ValueError):
    pass


def _require(path: Path) -> Path:
    if not path.is_file():
        raise FigureInputError(f"missing frozen input: {path}")
    return path


def load_inputs(
    joint_root: Path,
    geometry_root: Path,
    h5c_root: Path,
    h5d_root: Path,
) -> dict[str, Any]:
    surface = pd.read_csv(_require(joint_root / "joint_surface_classification.csv.gz"))
    summary = pd.read_csv(_require(joint_root / "joint_robustness_summary.csv"))
    envelope = pd.read_csv(_require(joint_root / "partial_identification_envelope.csv"))
    manifest = json.loads(_require(joint_root / "joint_observation_bias_manifest.json").read_text())
    geometry = pd.read_csv(_require(geometry_root / "observed_geometry_cross_scope.csv"))
    h5c = json.loads(_require(h5c_root / "h5c_observed_result.json").read_text())
    h5d = pd.read_csv(_require(h5d_root / "h5d_identifiability_summary.csv"))

    surface_required = {
        "evidence_scope",
        "target",
        "stratum",
        "trait_resolution_odds_ratio",
        "median_distance_completeness",
        "distance_completeness_odds_ratio",
        "state_recording_odds_ratio",
        "robust_cell",
    }
    summary_required = {
        "evidence_scope",
        "target",
        "stratum",
        "n_fit_cells",
        "n_robust_cells",
        "robust_fraction_of_fit_grid",
    }
    envelope_required = {
        "evidence_scope",
        "target",
        "stratum",
        "estimate_lower",
        "estimate_upper",
        "expected_sign_identified",
        "support_identified_across_envelope",
        "envelope_robust",
    }
    for frame, required, label in (
        (surface, surface_required, "joint surface"),
        (summary, summary_required, "joint summary"),
        (envelope, envelope_required, "partial-identification envelope"),
    ):
        missing = required - set(frame.columns)
        if missing:
            raise FigureInputError(f"{label} missing columns: {sorted(missing)}")

    if manifest.get("contract") != "chapter1_joint_observation_bias_v1":
        raise FigureInputError("joint manifest contract mismatch")
    if manifest.get("grid_fraction_is_probability") is not False:
        raise FigureInputError("joint grid fraction must remain non-probabilistic")
    if len(geometry) != 12 or set(geometry["classification"].astype(str)) != {
        "monotonic_or_unresolved"
    }:
        raise FigureInputError("geometry artifact differs from frozen 12-cell boundary")
    if h5c.get("classification") != "no_pollination_mode_specificity_support":
        raise FigureInputError("H5c classification mismatch")
    if len(h5d) != 8 or pd.Series(h5d["qualified"]).astype(bool).any():
        raise FigureInputError("H5d artifact differs from frozen 8-cell boundary")
    return {
        "surface": surface,
        "summary": summary,
        "envelope": envelope,
        "manifest": manifest,
        "geometry": geometry,
        "h5c": h5c,
        "h5d": h5d,
    }


def _primary_surface_matrix(surface: pd.DataFrame, target: str) -> tuple[np.ndarray, list[float], list[float]]:
    sub = surface.loc[
        surface["evidence_scope"].astype(str).eq("direct_only")
        & surface["stratum"].astype(str).eq("native_nonendemic")
        & surface["target"].astype(str).eq(target)
        & surface["surface_type"].astype(str).eq("joint_selection_grid")
        & surface["status"].astype(str).eq("fit")
    ].copy()
    if sub.empty:
        raise FigureInputError(f"no primary joint cells for {target}")
    sub["trait_resolution_odds_ratio"] = pd.to_numeric(
        sub["trait_resolution_odds_ratio"], errors="raise"
    )
    sub["state_recording_odds_ratio"] = pd.to_numeric(
        sub["state_recording_odds_ratio"], errors="raise"
    )
    sub["robust_cell"] = sub["robust_cell"].astype(bool)
    ors_r = sorted(sub["trait_resolution_odds_ratio"].unique().tolist())
    ors_d = sorted(sub["state_recording_odds_ratio"].unique().tolist())
    grouped = (
        sub.groupby(["trait_resolution_odds_ratio", "state_recording_odds_ratio"], sort=True)[
            "robust_cell"
        ]
        .mean()
        .unstack()
        .reindex(index=ors_r, columns=ors_d)
    )
    if grouped.isna().any().any():
        raise FigureInputError(f"incomplete OR_R×OR_D surface for {target}")
    return grouped.to_numpy(float), ors_r, ors_d


def _heatmap_panel(ax: plt.Axes, surface: pd.DataFrame, target: str, title: str) -> dict[str, Any]:
    matrix, ors_r, ors_d = _primary_surface_matrix(surface, target)
    image = ax.imshow(matrix, vmin=0.0, vmax=1.0, aspect="auto", origin="lower")
    ax.set_xticks(np.arange(len(ors_d)))
    ax.set_xticklabels([f"{x:g}" for x in ors_d], rotation=45, ha="right", fontsize=7.2)
    ax.set_yticks(np.arange(len(ors_r)))
    ax.set_yticklabels([f"{x:g}" for x in ors_r], fontsize=7.2)
    ax.set_xlabel("V6 state-recording OR_D", fontsize=8.2)
    ax.set_ylabel("V5 trait-resolution OR_R", fontsize=8.2)
    ax.set_title(title, loc="left", fontsize=10.3, fontweight="bold")
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            value = matrix[row, col]
            label = "R" if np.isclose(value, 1.0) else f"{value:.2f}"
            ax.text(col, row, label, ha="center", va="center", fontsize=6.2)
    return {
        "target": target,
        "min_subgrid_robust_fraction": float(matrix.min()),
        "max_subgrid_robust_fraction": float(matrix.max()),
        "n_fully_robust_or_cells": int(np.isclose(matrix, 1.0).sum()),
        "n_or_cells": int(matrix.size),
        "image": image,
    }


def _envelope_panel(ax: plt.Axes, envelope: pd.DataFrame) -> dict[str, Any]:
    primary = envelope.loc[
        envelope["evidence_scope"].astype(str).eq("direct_only")
        & envelope["stratum"].astype(str).eq("native_nonendemic")
    ].copy()
    scalar_targets = ["Palearctic_accessibility", "tropical_accessibility"]
    labels = ["Palearctic accessibility", "Tropical accessibility"]
    rows = []
    for target in scalar_targets:
        hit = primary.loc[primary["target"].astype(str).eq(target)]
        if len(hit) != 1:
            raise FigureInputError(f"partial-identification envelope missing {target}")
        rows.append(hit.iloc[0])
    y = np.array([1.0, 0.0])
    ax.axvline(0, linestyle="--", linewidth=0.9)
    for pos, row in zip(y, rows, strict=True):
        lo = float(row["estimate_lower"])
        hi = float(row["estimate_upper"])
        mid = (lo + hi) / 2.0
        ax.errorbar([mid], [pos], xerr=[[mid - lo], [hi - mid]], fmt="o", capsize=3)
        robust = bool(row["envelope_robust"])
        ax.text(hi, pos + 0.13, "ROBUST" if robust else "FRAGILE", fontsize=7.5, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Partial-identification slope envelope", fontsize=8.2)
    ax.set_title("C  Partial identification exposes residual fragility", loc="left", fontsize=10.3, fontweight="bold")

    vector_status: dict[str, bool] = {}
    for target, label in (("Palearctic_vector", "Pal vector"), ("tropical_vector", "Tropical vector"), ("north_tropical_vector_difference", "North–Trop difference")):
        hit = primary.loc[primary["target"].astype(str).eq(target)]
        if len(hit) != 1:
            raise FigureInputError(f"partial-identification envelope missing {target}")
        vector_status[target] = bool(hit.iloc[0]["envelope_robust"])
    text = "  |  ".join(
        f"{label}: {'robust' if vector_status[target] else 'fragile'}"
        for target, label in (("Palearctic_vector", "Pal vector"), ("tropical_vector", "Trop vector"), ("north_tropical_vector_difference", "N–T"))
    )
    ax.text(0.0, -0.24, text, transform=ax.transAxes, fontsize=7.2, va="top")
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="x", labelsize=7.5)
    return {
        "scalar_envelope_robust": {
            target: bool(row["envelope_robust"])
            for target, row in zip(scalar_targets, rows, strict=True)
        },
        "vector_envelope_robust": vector_status,
    }


def _claim_boundary_panel(
    ax: plt.Axes,
    geometry: pd.DataFrame,
    h5c: dict[str, Any],
    h5d: pd.DataFrame,
) -> dict[str, Any]:
    ax.axis("off")
    est = float(h5c["interaction_estimate"])
    lo = float(h5c["interaction_ci_low"])
    hi = float(h5c["interaction_ci_high"])
    lines = [
        "D  Mechanism / geometry gates remain closed",
        "",
        f"Response geometry: 0/{len(geometry)} promoted",
        f"H5c pollination specificity: {est:.3f} [{lo:.3f}, {hi:.3f}]",
        f"H5c p = {float(h5c['interaction_p_value']):.3f}",
        f"H5d threshold identifiability: 0/{len(h5d)} qualified",
        "",
        "Observation robustness cannot promote a pollinator mechanism.",
        "Non-identification is retained as a result, not inverted into absence.",
    ]
    ax.text(0.0, 1.0, "\n".join(lines), transform=ax.transAxes, va="top", fontsize=9.0, linespacing=1.5)
    return {
        "geometry_promoted": 0,
        "geometry_cells": int(len(geometry)),
        "h5c_interaction_estimate": est,
        "h5c_ci": [lo, hi],
        "h5c_p_value": float(h5c["interaction_p_value"]),
        "h5d_qualified": 0,
        "h5d_cells": int(len(h5d)),
    }


def render_figure(
    *,
    joint_root: Path,
    geometry_root: Path,
    h5c_root: Path,
    h5d_root: Path,
    output_dir: Path,
    source_receipts: dict[str, Any],
) -> dict[str, Any]:
    tables = load_inputs(joint_root, geometry_root, h5c_root, h5d_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(13.4, 9.2), constrained_layout=False)
    grid = GridSpec(2, 2, figure=fig, left=0.10, right=0.96, bottom=0.10, top=0.88, hspace=0.42, wspace=0.35)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    pal = _heatmap_panel(ax_a, tables["surface"], "Palearctic_vector", "A  Palearctic joint robustness surface")
    trop = _heatmap_panel(ax_b, tables["surface"], "tropical_vector", "B  Tropical joint robustness surface")
    cbar = fig.colorbar(pal.pop("image"), ax=[ax_a, ax_b], fraction=0.025, pad=0.02)
    cbar.set_label("Fraction robust across frozen C0 × OR_C subgrid\n(assumption-grid geometry, not probability)", fontsize=7.8)
    trop.pop("image")
    envelope = _envelope_panel(ax_c, tables["envelope"])
    boundaries = _claim_boundary_panel(ax_d, tables["geometry"], tables["h5c"], tables["h5d"])

    fig.suptitle(
        "Joint observation bias separates robust pattern from fragile identification",
        x=0.10,
        y=0.975,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.10,
        0.943,
        "V5 trait-resolution MNAR × V6 species-list detection; robust/fragile regions are predeclared assumption domains, not posterior probabilities.",
        fontsize=8.8,
    )
    basename = "chapter1_v11_figure4_joint_observation_boundaries"
    outputs = []
    for suffix in ("png", "svg", "pdf"):
        path = output_dir / f"{basename}.{suffix}"
        if suffix == "png":
            fig.savefig(path, dpi=300, bbox_inches="tight")
        else:
            fig.savefig(path, bbox_inches="tight")
        outputs.append(path.name)
    plt.close(fig)

    manifest = {
        "contract": CONTRACT,
        "status": "rendered_from_frozen_joint_P3_and_existing_identification_boundaries",
        "sources": source_receipts,
        "primary_profile": "direct_only_native_nonendemic",
        "palearctic_surface": pal,
        "tropical_surface": trop,
        "partial_identification": envelope,
        "claim_boundaries": boundaries,
        "grid_fraction_is_probability": False,
        "new_biological_models_fitted": False,
        "new_p_values_generated": False,
        "outputs": outputs,
        "claim_boundary": (
            "Figure 4 shows where the frozen plant-side result survives a joint observation-bias domain, "
            "where partial identification remains fragile, and why observation robustness does not promote "
            "response geometry or pollination-mechanism claims."
        ),
    }
    (output_dir / "chapter1_v11_figure4_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("render")
def render_command(
    joint_root: Path = typer.Option(..., exists=True, file_okay=False),
    geometry_root: Path = typer.Option(..., exists=True, file_okay=False),
    h5c_root: Path = typer.Option(..., exists=True, file_okay=False),
    h5d_root: Path = typer.Option(..., exists=True, file_okay=False),
    output_dir: Path = typer.Option(...),
    source_receipts_json: Path = typer.Option(..., exists=True, dir_okay=False),
) -> None:
    receipts = json.loads(source_receipts_json.read_text(encoding="utf-8"))
    typer.echo(
        json.dumps(
            render_figure(
                joint_root=joint_root,
                geometry_root=geometry_root,
                h5c_root=h5c_root,
                h5d_root=h5d_root,
                output_dir=output_dir,
                source_receipts=receipts,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    app()
