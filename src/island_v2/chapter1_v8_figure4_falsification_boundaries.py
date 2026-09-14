"""Render Chapter 1 v8 Figure 4 from four frozen falsification artifacts.

Presentation only: no biological model is refit and no new inferential quantity is
used beyond transparent summaries of frozen outputs.
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


class FigureInputError(ValueError):
    pass


def _require(path: Path) -> Path:
    if not path.is_file():
        raise FigureInputError(f"missing frozen input: {path}")
    return path


def load_inputs(v6_root: Path, geometry_root: Path, h5c_root: Path, h5d_root: Path) -> dict[str, Any]:
    v6 = pd.read_csv(_require(v6_root / "species_detection_tipping_surface.csv"))
    geom = pd.read_csv(_require(geometry_root / "observed_geometry_cross_scope.csv"))
    h5c = json.loads(_require(h5c_root / "h5c_observed_result.json").read_text())
    h5d = pd.read_csv(_require(h5d_root / "h5d_identifiability_summary.csv"))

    v6_required = {"target", "baseline_supported", "tipping_event"}
    geom_required = {
        "all_observed_D",
        "all_critical_D",
        "direct_observed_D",
        "direct_critical_D",
        "classification",
    }
    h5d_required = {
        "classification_accuracy",
        "false_distributed_selection_under_smooth",
        "qualified",
    }
    for frame, required, label in (
        (v6, v6_required, "V6"),
        (geom, geom_required, "geometry"),
        (h5d, h5d_required, "H5d"),
    ):
        missing = required - set(frame.columns)
        if missing:
            raise FigureInputError(f"{label} missing columns: {sorted(missing)}")

    expected_targets = {
        "Palearctic_accessibility",
        "tropical_accessibility",
        "north_tropical_vector_difference",
    }
    if set(v6["target"].astype(str)) != expected_targets:
        raise FigureInputError("V6 target set differs from frozen contract")
    if len(geom) != 12 or set(geom["classification"].astype(str)) != {
        "monotonic_or_unresolved"
    }:
        raise FigureInputError(
            "geometry artifact must contain 12 monotonic_or_unresolved cross-scope cells"
        )
    if len(h5d) != 8 or pd.Series(h5d["qualified"]).astype(bool).any():
        raise FigureInputError("H5d artifact must contain 8 non-qualified cells")
    if h5c.get("classification") != "no_pollination_mode_specificity_support":
        raise FigureInputError("H5c classification differs from frozen result")
    return {"v6": v6, "geometry": geom, "h5c": h5c, "h5d": h5d}


def _panel_v6(ax: plt.Axes, v6: pd.DataFrame) -> dict[str, tuple[int, int]]:
    labels = [
        "Palearctic accessibility",
        "North–Tropical vector",
        "Tropical accessibility",
    ]
    targets = [
        "Palearctic_accessibility",
        "north_tropical_vector_difference",
        "tropical_accessibility",
    ]
    summary: dict[str, tuple[int, int]] = {}
    survive = []
    for target in targets:
        sub = v6.loc[
            v6["target"].astype(str).eq(target)
            & v6["baseline_supported"].astype(bool)
        ]
        total = len(sub)
        kept = int(
            sub["tipping_event"].astype(str).eq("none_on_frozen_grid").sum()
        )
        summary[target] = (kept, total)
        survive.append(kept / total)
    y = np.arange(3)
    ax.barh(y, survive)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.2)
    ax.set_xlim(0, 1.02)
    ax.invert_yaxis()
    ax.set_xlabel("Fraction surviving frozen species-detection grid", fontsize=8.5)
    ax.set_title(
        "A  Observation-bias defense is asymmetric",
        loc="left",
        fontsize=10.5,
        fontweight="bold",
    )
    for i, (target, frac) in enumerate(zip(targets, survive, strict=True)):
        kept, total = summary[target]
        ax.text(
            min(frac + 0.02, 0.90),
            i,
            f"{kept}/{total}",
            va="center",
            fontsize=8.5,
            fontweight="bold",
        )
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="x", labelsize=8)
    return summary


def _panel_geometry(ax: plt.Axes, geom: pd.DataFrame) -> None:
    x = pd.to_numeric(geom["all_observed_D"]) - pd.to_numeric(geom["all_critical_D"])
    y = pd.to_numeric(geom["direct_observed_D"]) - pd.to_numeric(
        geom["direct_critical_D"]
    )
    ax.axhline(0, color="0.55", linestyle="--", linewidth=0.9)
    ax.axvline(0, color="0.55", linestyle="--", linewidth=0.9)
    ax.scatter(x, y, s=32)
    ax.set_xlabel("All-analysis nonlinear margin (D − critical D)", fontsize=8.2)
    ax.set_ylabel("Direct-only nonlinear margin", fontsize=8.2)
    ax.set_title(
        "B  A nonlinear shape required both scopes",
        loc="left",
        fontsize=10.5,
        fontweight="bold",
    )
    ax.text(
        0.98,
        0.96,
        "promotion quadrant",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7.7,
        color="0.35",
    )
    ax.text(
        0.03,
        0.05,
        "0/12 promoted",
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
    )
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(labelsize=7.8)


def _panel_h5c(ax: plt.Axes, result: dict[str, Any]) -> None:
    est = float(result["interaction_estimate"])
    lo = float(result["interaction_ci_low"])
    hi = float(result["interaction_ci_high"])
    ax.axvline(0, color="0.55", linestyle="--", linewidth=0.9)
    ax.errorbar(
        [est],
        [0],
        xerr=[[est - lo], [hi - est]],
        fmt="o",
        capsize=3,
        markersize=6,
    )
    ax.set_yticks([0])
    ax.set_yticklabels(["distance × biotic"], fontsize=8.5)
    ax.set_ylim(-0.8, 0.8)
    ax.set_xlabel("Interaction estimate", fontsize=8.5)
    ax.set_title(
        "C  Independent pollination-mode specificity",
        loc="left",
        fontsize=10.5,
        fontweight="bold",
    )
    ax.text(
        0.03,
        0.12,
        f"p = {float(result['interaction_p_value']):.3f}\nnot promoted",
        transform=ax.transAxes,
        fontsize=8.5,
        fontweight="bold",
    )
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="x", labelsize=8)


def _panel_h5d(ax: plt.Axes, h5d: pd.DataFrame) -> None:
    acc = pd.to_numeric(h5d["classification_accuracy"])
    fp = pd.to_numeric(h5d["false_distributed_selection_under_smooth"])
    ax.axvline(0.80, color="0.55", linestyle="--", linewidth=0.9)
    ax.axhline(0.10, color="0.55", linestyle="--", linewidth=0.9)
    ax.scatter(acc, fp, s=34)
    ax.set_xlabel("Threshold-vs-cline classification accuracy", fontsize=8.2)
    ax.set_ylabel("False threshold selection under smooth clines", fontsize=8.2)
    ax.set_title(
        "D  Distributed thresholds were not identifiable",
        loc="left",
        fontsize=10.5,
        fontweight="bold",
    )
    ax.text(
        0.03,
        0.94,
        "required: accuracy ≥0.80\nand false selection ≤0.10",
        transform=ax.transAxes,
        va="top",
        fontsize=7.7,
        color="0.35",
    )
    ax.text(
        0.97,
        0.06,
        "0/8 qualified",
        transform=ax.transAxes,
        ha="right",
        fontsize=9,
        fontweight="bold",
    )
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(labelsize=7.8)


def render_figure(
    *,
    v6_root: Path,
    geometry_root: Path,
    h5c_root: Path,
    h5d_root: Path,
    output_dir: Path,
    source_receipts: dict[str, Any],
) -> dict[str, Any]:
    tables = load_inputs(v6_root, geometry_root, h5c_root, h5d_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(13.2, 8.8), constrained_layout=False)
    grid = GridSpec(
        2,
        2,
        figure=fig,
        left=0.09,
        right=0.97,
        bottom=0.10,
        top=0.88,
        hspace=0.40,
        wspace=0.34,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])
    v6_summary = _panel_v6(ax_a, tables["v6"])
    _panel_geometry(ax_b, tables["geometry"])
    _panel_h5c(ax_c, tables["h5c"])
    _panel_h5d(ax_d, tables["h5d"])
    fig.suptitle(
        "Cross-examination defines what the global data can and cannot support",
        x=0.09,
        y=0.975,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.09,
        0.944,
        "Frozen falsification and identifiability audits; adverse evidence is shown rather than rescued.",
        fontsize=9,
        color="0.35",
    )
    basename = "chapter1_v8_figure4_falsification_boundaries"
    png = output_dir / f"{basename}.png"
    svg = output_dir / f"{basename}.svg"
    pdf = output_dir / f"{basename}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    manifest = {
        "contract": "chapter1_v8_figure4_falsification_boundaries_v1",
        "status": "rendered_from_four_frozen_audits",
        "sources": source_receipts,
        "v6_survival": {
            key: {"survived": value[0], "baseline_supported": value[1]}
            for key, value in v6_summary.items()
        },
        "geometry_promoted": 0,
        "geometry_cells": int(len(tables["geometry"])),
        "h5c_interaction_estimate": float(tables["h5c"]["interaction_estimate"]),
        "h5c_p_value": float(tables["h5c"]["interaction_p_value"]),
        "h5d_qualified": int(
            pd.Series(tables["h5d"]["qualified"]).astype(bool).sum()
        ),
        "h5d_cells": int(len(tables["h5d"])),
        "new_biological_models_fitted": False,
        "new_p_values_generated": False,
        "claim_boundary": (
            "Figure 4 visualizes frozen stress tests and failed/non-identified "
            "mechanism gates; it does not promote a mechanism from negative evidence."
        ),
        "outputs": [png.name, svg.name, pdf.name],
    }
    (output_dir / "chapter1_v8_figure4_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest


@app.command("render")
def render_command(
    v6_root: Path = typer.Option(..., exists=True, file_okay=False),
    geometry_root: Path = typer.Option(..., exists=True, file_okay=False),
    h5c_root: Path = typer.Option(..., exists=True, file_okay=False),
    h5d_root: Path = typer.Option(..., exists=True, file_okay=False),
    output_dir: Path = typer.Option(...),
    source_receipts_json: Path = typer.Option(..., exists=True, dir_okay=False),
) -> None:
    receipts = json.loads(source_receipts_json.read_text())
    typer.echo(
        json.dumps(
            render_figure(
                v6_root=v6_root,
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
