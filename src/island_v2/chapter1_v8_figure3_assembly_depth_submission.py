"""Submission layout for Chapter 1 v8 Figure 3.

Uses the validated frozen-input loaders and panel primitives from
``chapter1_v8_figure3_assembly_depth``. Only composition and label spacing differ.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import typer  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

from island_v2.chapter1_v8_figure3_assembly_depth import (  # noqa: E402
    SCOPES,
    SCOPE_LABELS,
    STRATA,
    STRATUM_LABELS,
    _attenuation_panel,
    _axis_stage_panel,
    load_inputs,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _support_panel_compact(ax: plt.Axes, classification) -> None:
    rows = [(scope, stratum) for stratum in STRATA for scope in SCOPES]
    matrix = []
    for scope, stratum in rows:
        row = classification.loc[
            classification["evidence_scope"].astype(str).eq(scope)
            & classification["stratum"].astype(str).eq(stratum)
        ].iloc[0]
        matrix.append(
            [
                int(row["observed_supported_modes"]),
                int(row["after_family_supported_modes"]),
                int(row["after_genus_supported_modes"]),
            ]
        )
    array = np.asarray(matrix, dtype=float)
    ax.imshow(array, vmin=0, vmax=4, cmap="Greys", aspect="auto")
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            value = int(array[i, j])
            ax.text(
                j,
                i,
                f"{value}/4",
                ha="center",
                va="center",
                fontsize=10,
                fontweight="bold",
                color="white" if value >= 3 else "black",
            )
    short_scope = {"all_analysis_eligible": "All", "direct_only": "Direct"}
    short_stratum = {"all_native": "native", "native_nonendemic": "NNE"}
    ax.set_xticks(np.arange(3))
    ax.set_xticklabels(["Observed", "Family", "Genus"], fontsize=8.2)
    ax.set_yticks(np.arange(4))
    ax.set_yticklabels(
        [f"{short_scope[s]} · {short_stratum[t]}" for s, t in rows],
        fontsize=7.7,
    )
    ax.set_title("D  Vector gate", loc="left", fontsize=10.5, fontweight="bold")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(
        0.5,
        -0.19,
        "4/4 → 4/4 → 0/4\nin every cell",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=8.1,
        fontweight="bold",
    )


def render_submission(
    *,
    effect_root: Path,
    pr142_root: Path,
    output_dir: Path,
    effect_run_id: int,
    effect_artifact_id: int,
    effect_artifact_digest: str,
    primary_run_id: int,
    primary_artifact_id: int,
    primary_artifact_digest: str,
) -> dict[str, Any]:
    tables = load_inputs(effect_root, pr142_root)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(14.3, 9.4), constrained_layout=False)
    outer = GridSpec(
        2,
        1,
        figure=fig,
        height_ratios=[1.08, 1.0],
        hspace=0.38,
        top=0.90,
        bottom=0.09,
        left=0.075,
        right=0.975,
    )
    ax_a = fig.add_subplot(outer[0, 0])
    bottom = outer[1, 0].subgridspec(1, 3, width_ratios=[1.55, 1.20, 0.95], wspace=0.36)
    ax_b = fig.add_subplot(bottom[0, 0])
    ax_c = fig.add_subplot(bottom[0, 1])
    ax_d = fig.add_subplot(bottom[0, 2])

    _attenuation_panel(ax_a, tables["attenuation"])
    _axis_stage_panel(ax_b, tables["effects"], "generalized_accessible", "B")
    _axis_stage_panel(ax_c, tables["effects"], "selfing_core", "C")
    _support_panel_compact(ax_d, tables["classification"])

    handles, labels = ax_b.get_legend_handles_labels()
    ax_b.get_legend().remove() if ax_b.get_legend() else None
    fig.legend(
        handles,
        labels,
        frameon=False,
        fontsize=7.4,
        loc="upper center",
        bbox_to_anchor=(0.51, 0.505),
        ncol=4,
        handlelength=2.3,
        columnspacing=1.2,
    )

    fig.suptitle(
        "The strongest Palearctic floral syndrome has a taxonomic assembly depth",
        x=0.075,
        y=0.975,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.075,
        0.945,
        "Frozen source-matched decomposition; attenuation localizes hierarchical expression and is not causal mediation.",
        fontsize=9,
        color="0.35",
    )

    basename = "chapter1_v8_figure3_assembly_depth_submission"
    png_path = output_dir / f"{basename}.png"
    svg_path = output_dir / f"{basename}.svg"
    pdf_path = output_dir / f"{basename}.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    tables["attenuation"].to_csv(output_dir / "figure3_taxonomic_attenuation_profiles.csv", index=False)
    tables["effects"].to_csv(output_dir / "figure3_canonical_stage_effects.csv", index=False)
    tables["classification"].to_csv(output_dir / "figure3_palearctic_vector_gate.csv", index=False)

    manifest = {
        "contract": "chapter1_v8_figure3_assembly_depth_submission_v1",
        "status": "submission_layout_rendered_from_frozen_taxonomic_depth_artifacts",
        "effect_fingerprint": {
            "workflow_run_id": int(effect_run_id),
            "artifact_id": int(effect_artifact_id),
            "artifact_digest": str(effect_artifact_digest),
        },
        "primary_analysis": {
            "workflow_run_id": int(primary_run_id),
            "artifact_id": int(primary_artifact_id),
            "artifact_digest": str(primary_artifact_digest),
        },
        "n_attenuation_profiles": int(len(tables["attenuation"])),
        "n_canonical_stage_effect_rows": int(len(tables["effects"])),
        "n_classification_cells": int(len(tables["classification"])),
        "frozen_support_ladder": "4/4 -> 4/4 -> 0/4",
        "new_biological_models_fitted": False,
        "new_p_values_generated": False,
        "layout_revision": "compact_bottom_row_no_label_overlap_v1",
        "claim_boundary": (
            "Submission Figure 3 visualizes the same frozen taxonomic-depth quantities as "
            "the validated renderer. No biological estimand or test is changed."
        ),
        "outputs": [png_path.name, svg_path.name, pdf_path.name],
    }
    (output_dir / "chapter1_v8_figure3_submission_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("render")
def render_command(
    effect_root: Path = typer.Option(..., exists=True, file_okay=False),
    pr142_root: Path = typer.Option(..., exists=True, file_okay=False),
    output_dir: Path = typer.Option(...),
    effect_run_id: int = typer.Option(...),
    effect_artifact_id: int = typer.Option(...),
    effect_artifact_digest: str = typer.Option(...),
    primary_run_id: int = typer.Option(...),
    primary_artifact_id: int = typer.Option(...),
    primary_artifact_digest: str = typer.Option(...),
) -> None:
    typer.echo(
        json.dumps(
            render_submission(
                effect_root=effect_root,
                pr142_root=pr142_root,
                output_dir=output_dir,
                effect_run_id=effect_run_id,
                effect_artifact_id=effect_artifact_id,
                effect_artifact_digest=effect_artifact_digest,
                primary_run_id=primary_run_id,
                primary_artifact_id=primary_artifact_id,
                primary_artifact_digest=primary_artifact_digest,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    app()
