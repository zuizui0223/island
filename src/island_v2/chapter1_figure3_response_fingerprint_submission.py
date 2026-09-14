"""Submission-layout renderer for Chapter 1 Figure 3.

This module reuses the fail-closed data validation and panel primitives from
``chapter1_figure3_response_fingerprint``. It changes only composition/heading placement
for the manuscript-ready output; no data transformation or inferential quantity differs.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import typer
from matplotlib.gridspec import GridSpec

from island_v2.chapter1_figure3_response_fingerprint import (
    _render_angle_panel,
    _render_attenuation_panel,
    _render_atomic_panel,
    load_inputs,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def render_submission_figure(
    *,
    input_dir: Path,
    output_dir: Path,
    source_run_id: int,
    source_artifact_id: int,
    source_artifact_digest: str,
) -> dict:
    tables = load_inputs(input_dir)
    atomic = tables["atomic"]
    angles = tables["angles"]
    attenuation = tables["attenuation"]
    output_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(15.8, 10.5), constrained_layout=False)
    grid = GridSpec(
        2,
        2,
        figure=fig,
        width_ratios=[2.35, 1.0],
        height_ratios=[1.15, 1.0],
        hspace=0.31,
        wspace=0.25,
    )
    atomic_axes = _render_atomic_panel(fig, grid[:, 0], atomic)
    for text_artist in list(atomic_axes[0].texts):
        if text_artist.get_text() in {"A", "Component response fingerprints"}:
            text_artist.remove()
    _render_angle_panel(fig, grid[0, 1], angles)
    _render_attenuation_panel(fig, grid[1, 1], attenuation)

    fig.suptitle(
        "Isolation-response fingerprints differ in component composition and hierarchical expression",
        x=0.055,
        y=0.985,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.055,
        0.956,
        "Effect-size synthesis of the frozen Chapter 1 dataset; descriptive panels do not replace the preregistered multivariate tests.",
        fontsize=9.2,
        color="0.35",
    )
    fig.text(
        0.055,
        0.918,
        "A  Component response fingerprints",
        fontsize=11.5,
        fontweight="bold",
        ha="left",
        va="top",
    )

    basename = "chapter1_figure3_response_fingerprint_submission"
    png_path = output_dir / f"{basename}.png"
    svg_path = output_dir / f"{basename}.svg"
    pdf_path = output_dir / f"{basename}.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    atomic.to_csv(output_dir / "figure3_panel_a_atomic_effects.csv", index=False)
    angles.to_csv(output_dir / "figure3_panel_b_vector_angles.csv", index=False)
    attenuation.to_csv(output_dir / "figure3_panel_c_taxonomic_attenuation.csv", index=False)

    manifest = {
        "contract": "chapter1_figure3_response_fingerprint_submission_v1",
        "status": "submission_figure_rendered_from_frozen_effect_fingerprint",
        "source_workflow_run_id": int(source_run_id),
        "source_artifact_id": int(source_artifact_id),
        "source_artifact_digest": str(source_artifact_digest),
        "layout_revision": "submission_spacing_v1",
        "n_atomic_rows": int(len(atomic)),
        "n_vector_angle_rows": int(len(angles)),
        "n_taxonomic_profiles": int(len(attenuation)),
        "angle_degree_range": [
            float(pd.to_numeric(angles["vector_angle_degrees"]).min()),
            float(pd.to_numeric(angles["vector_angle_degrees"]).max()),
        ],
        "family_attenuation_fraction_range": [
            float(pd.to_numeric(attenuation["family_attenuation_fraction"]).min()),
            float(pd.to_numeric(attenuation["family_attenuation_fraction"]).max()),
        ],
        "genus_attenuation_fraction_range": [
            float(pd.to_numeric(attenuation["genus_attenuation_fraction"]).min()),
            float(pd.to_numeric(attenuation["genus_attenuation_fraction"]).max()),
        ],
        "new_biological_models_fitted": False,
        "new_p_values_generated": False,
        "claim_boundary": (
            "Submission Figure 3 is a visualization of frozen estimates. Atomic intervals "
            "are descriptive decomposition, vector angles are descriptive geometry, and "
            "attenuation is not causal mediation. Frozen H1-H3 tests remain inferential."
        ),
        "outputs": [png_path.name, svg_path.name, pdf_path.name],
    }
    (output_dir / "chapter1_figure3_response_fingerprint_submission_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("render")
def render_command(
    input_dir: Path = typer.Option(..., exists=True, file_okay=False),
    output_dir: Path = typer.Option(...),
    source_run_id: int = typer.Option(...),
    source_artifact_id: int = typer.Option(...),
    source_artifact_digest: str = typer.Option(...),
) -> None:
    manifest = render_submission_figure(
        input_dir=input_dir,
        output_dir=output_dir,
        source_run_id=source_run_id,
        source_artifact_id=source_artifact_id,
        source_artifact_digest=source_artifact_digest,
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
