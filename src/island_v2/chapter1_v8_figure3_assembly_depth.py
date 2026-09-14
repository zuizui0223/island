"""Render Chapter 1 v8 Figure 3: taxonomic assembly depth.

The renderer consumes two already-frozen artifacts: the PR142 primary-analysis artifact
and the post-freeze effect-fingerprint synthesis artifact. It performs no biological
refitting and generates no new p-values. Taxonomic attenuation is visualized as an
estimand, not causal mediation.
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

SCOPES = ["all_analysis_eligible", "direct_only"]
SCOPE_LABELS = {"all_analysis_eligible": "All-analysis", "direct_only": "Direct-only"}
SCOPE_COLOURS = {"all_analysis_eligible": "#3569a8", "direct_only": "#d0702f"}
STRATA = ["all_native", "native_nonendemic"]
STRATUM_LABELS = {"all_native": "All native", "native_nonendemic": "Native non-endemic"}
LINESTYLES = {"all_native": "-", "native_nonendemic": "--"}
SOURCE_MODES = ["geo_k5", "geo_k10", "geo_k20", "geo50_climate10"]
STAGES = ["observed_score", "after_family_residual", "after_genus_residual"]
STAGE_LABELS = ["Observed", "After family", "After genus"]
AXES = ["generalized_accessible", "selfing_core"]
AXIS_LABELS = {
    "generalized_accessible": "Accessibility / generalization",
    "selfing_core": "Reproductive assurance",
}
CANONICAL_SOURCE_MODE = "geo50_climate10"


class FigureInputError(ValueError):
    """Raised when a frozen figure input does not match the expected contract."""


def _require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    if missing := required - set(frame.columns):
        raise FigureInputError(f"{label} missing columns: {sorted(missing)}")


def load_inputs(effect_root: Path, pr142_root: Path) -> dict[str, pd.DataFrame]:
    attenuation_path = effect_root / "taxonomic_vector_attenuation.csv"
    effects_path = effect_root / "taxonomic_effects_long.csv"
    for path in (attenuation_path, effects_path):
        if not path.is_file():
            raise FigureInputError(f"missing frozen effect-fingerprint table: {path}")

    attenuation = pd.read_csv(attenuation_path)
    effects = pd.read_csv(effects_path)
    _require_columns(
        attenuation,
        {
            "evidence_scope",
            "source_mode",
            "stratum",
            "context",
            "observed_score_vector_norm",
            "after_family_residual_vector_norm",
            "after_genus_residual_vector_norm",
            "family_attenuation_fraction",
            "genus_attenuation_fraction",
            "conditional_genus_attenuation_fraction",
        },
        "taxonomic attenuation",
    )
    expected_attenuation = {
        (scope, stratum, mode)
        for scope in SCOPES
        for stratum in STRATA
        for mode in SOURCE_MODES
    }
    found_attenuation = set(
        attenuation[["evidence_scope", "stratum", "source_mode"]]
        .astype(str)
        .itertuples(index=False, name=None)
    )
    if found_attenuation != expected_attenuation or len(attenuation) != 16:
        raise FigureInputError("attenuation table must contain all 16 frozen profiles")
    if set(attenuation["context"].astype(str)) != {"Palearctic"}:
        raise FigureInputError("attenuation table must contain only the frozen Palearctic profiles")

    _require_columns(
        effects,
        {
            "evidence_scope",
            "source_mode",
            "stratum",
            "context",
            "stage",
            "axis",
            "distance_slope",
            "ci_low",
            "ci_high",
        },
        "taxonomic effects",
    )
    effects = effects.loc[
        effects["source_mode"].astype(str).eq(CANONICAL_SOURCE_MODE)
        & effects["context"].astype(str).eq("Palearctic")
        & effects["evidence_scope"].astype(str).isin(SCOPES)
        & effects["stratum"].astype(str).isin(STRATA)
        & effects["stage"].astype(str).isin(STAGES)
        & effects["axis"].astype(str).isin(AXES)
    ].copy()
    expected_effects = {
        (scope, stratum, stage, axis)
        for scope in SCOPES
        for stratum in STRATA
        for stage in STAGES
        for axis in AXES
    }
    found_effects = set(
        effects[["evidence_scope", "stratum", "stage", "axis"]]
        .astype(str)
        .itertuples(index=False, name=None)
    )
    if found_effects != expected_effects or len(effects) != 24:
        raise FigureInputError("canonical taxonomic effects must contain exactly 24 frozen rows")

    classification_frames: list[pd.DataFrame] = []
    for scope in SCOPES:
        path = pr142_root / f"taxonomic-depth/{scope}/classification.csv"
        if not path.is_file():
            raise FigureInputError(f"missing frozen taxonomic classification: {path}")
        frame = pd.read_csv(path)
        _require_columns(
            frame,
            {
                "context_layer",
                "stratum",
                "context",
                "n_source_modes_expected",
                "observed_supported_modes",
                "after_family_supported_modes",
                "after_genus_supported_modes",
                "V2_taxonomic_depth_classification",
            },
            f"{scope} classification",
        )
        frame = frame.loc[
            frame["context_layer"].astype(str).eq("biogeographic_realm")
            & frame["context"].astype(str).eq("Palearctic")
            & frame["stratum"].astype(str).isin(STRATA)
        ].copy()
        frame.insert(0, "evidence_scope", scope)
        classification_frames.append(frame)
    classification = pd.concat(classification_frames, ignore_index=True)
    expected_classification = {(scope, stratum) for scope in SCOPES for stratum in STRATA}
    found_classification = set(
        classification[["evidence_scope", "stratum"]]
        .astype(str)
        .itertuples(index=False, name=None)
    )
    if found_classification != expected_classification or len(classification) != 4:
        raise FigureInputError("Palearctic taxonomic classification must contain four scope×stratum cells")
    for column in (
        "n_source_modes_expected",
        "observed_supported_modes",
        "after_family_supported_modes",
        "after_genus_supported_modes",
    ):
        classification[column] = pd.to_numeric(classification[column], errors="coerce")
    if not classification["n_source_modes_expected"].eq(4).all():
        raise FigureInputError("taxonomic classification must expect four frozen source modes")
    if not classification["observed_supported_modes"].eq(4).all():
        raise FigureInputError("frozen Palearctic observed stage must be 4/4 in all cells")
    if not classification["after_family_supported_modes"].eq(4).all():
        raise FigureInputError("frozen Palearctic family stage must be 4/4 in all cells")
    if not classification["after_genus_supported_modes"].eq(0).all():
        raise FigureInputError("frozen Palearctic genus stage must be 0/4 in all cells")

    numeric_att = [
        "observed_score_vector_norm",
        "after_family_residual_vector_norm",
        "after_genus_residual_vector_norm",
        "family_attenuation_fraction",
        "genus_attenuation_fraction",
        "conditional_genus_attenuation_fraction",
    ]
    for column in numeric_att:
        attenuation[column] = pd.to_numeric(attenuation[column], errors="coerce")
    for column in ("distance_slope", "ci_low", "ci_high"):
        effects[column] = pd.to_numeric(effects[column], errors="coerce")
    if attenuation[numeric_att].isna().any().any() or effects[["distance_slope", "ci_low", "ci_high"]].isna().any().any():
        raise FigureInputError("frozen Figure 3 inputs contain non-numeric plotted values")

    return {"attenuation": attenuation, "effects": effects, "classification": classification}


def _range_percent(series: pd.Series) -> str:
    values = pd.to_numeric(series, errors="coerce").to_numpy(float)
    return f"{100 * values.min():.1f}–{100 * values.max():.1f}%"


def _attenuation_panel(ax: plt.Axes, attenuation: pd.DataFrame) -> None:
    x = np.arange(3)
    normalized: list[np.ndarray] = []
    for row in attenuation.itertuples(index=False):
        observed = float(row.observed_score_vector_norm)
        values = np.asarray(
            [
                1.0,
                float(row.after_family_residual_vector_norm) / observed,
                float(row.after_genus_residual_vector_norm) / observed,
            ]
        )
        normalized.append(values)
        ax.plot(
            x,
            values,
            color=SCOPE_COLOURS[str(row.evidence_scope)],
            linestyle=LINESTYLES[str(row.stratum)],
            linewidth=1.0,
            alpha=0.35,
            marker="o",
            markersize=2.5,
        )
    matrix = np.vstack(normalized)
    median = np.median(matrix, axis=0)
    q25 = np.quantile(matrix, 0.25, axis=0)
    q75 = np.quantile(matrix, 0.75, axis=0)
    ax.fill_between(x, q25, q75, alpha=0.15, color="0.25")
    ax.plot(x, median, color="0.10", linewidth=2.8, marker="o", markersize=6, zorder=5)
    ax.set_xticks(x)
    ax.set_xticklabels(STAGE_LABELS, fontsize=9)
    ax.set_ylabel("Vector magnitude retained (observed = 1)", fontsize=9)
    ax.set_ylim(-0.02, 1.08)
    ax.set_title("A  Most attenuation occurs from family to genus", loc="left", fontsize=11, fontweight="bold")
    ax.text(
        0.02,
        0.06,
        f"Family attenuation: {_range_percent(attenuation['family_attenuation_fraction'])}\n"
        f"Genus attenuation: {_range_percent(attenuation['genus_attenuation_fraction'])}\n"
        f"Family→genus of remainder: {_range_percent(attenuation['conditional_genus_attenuation_fraction'])}",
        transform=ax.transAxes,
        fontsize=8.2,
        va="bottom",
    )
    ax.grid(axis="y", color="0.9", linewidth=0.7)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def _axis_stage_panel(ax: plt.Axes, effects: pd.DataFrame, axis: str, letter: str) -> None:
    subset = effects.loc[effects["axis"].astype(str).eq(axis)]
    x = np.arange(3)
    for scope in SCOPES:
        for stratum in STRATA:
            cell = (
                subset.loc[
                    subset["evidence_scope"].astype(str).eq(scope)
                    & subset["stratum"].astype(str).eq(stratum)
                ]
                .set_index("stage")
                .reindex(STAGES)
            )
            estimate = cell["distance_slope"].to_numpy(float)
            low = cell["ci_low"].to_numpy(float)
            high = cell["ci_high"].to_numpy(float)
            ax.errorbar(
                x,
                estimate,
                yerr=np.vstack([estimate - low, high - estimate]),
                color=SCOPE_COLOURS[scope],
                linestyle=LINESTYLES[stratum],
                marker="o",
                markersize=4,
                linewidth=1.1,
                capsize=2,
                alpha=0.85,
                label=f"{SCOPE_LABELS[scope]} · {STRATUM_LABELS[stratum]}",
            )
    ax.axhline(0, color="0.65", linestyle="--", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(STAGE_LABELS, fontsize=8.5)
    ax.set_ylabel("Isolation slope", fontsize=8.8)
    ax.set_title(f"{letter}  {AXIS_LABELS[axis]}", loc="left", fontsize=10.5, fontweight="bold")
    ax.tick_params(axis="y", labelsize=8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def _support_panel(ax: plt.Axes, classification: pd.DataFrame) -> None:
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
    ax.set_xticks(np.arange(3))
    ax.set_xticklabels(STAGE_LABELS, fontsize=8.5)
    ax.set_yticks(np.arange(4))
    ax.set_yticklabels(
        [f"{SCOPE_LABELS[s]} · {STRATUM_LABELS[t]}" for s, t in rows],
        fontsize=7.7,
    )
    ax.set_title("D  Frozen vector gate across source modes", loc="left", fontsize=10.5, fontweight="bold")
    ax.set_xlabel("supported source modes / 4", fontsize=8)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(
        0.5,
        -0.23,
        "Observed 4/4 → family 4/4 → genus 0/4 in all four cells",
        transform=ax.transAxes,
        ha="center",
        fontsize=8.2,
        fontweight="bold",
    )


def render_figure(
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

    fig = plt.figure(figsize=(13.4, 9.5), constrained_layout=False)
    grid = GridSpec(2, 2, figure=fig, hspace=0.34, wspace=0.34)
    ax_a = fig.add_subplot(grid[0, :])
    ax_b = fig.add_subplot(grid[1, 0])
    inner = grid[1, 1].subgridspec(1, 2, wspace=0.52)
    ax_c = fig.add_subplot(inner[0, 0])
    ax_d = fig.add_subplot(inner[0, 1])

    _attenuation_panel(ax_a, tables["attenuation"])
    _axis_stage_panel(ax_b, tables["effects"], "generalized_accessible", "B")
    _axis_stage_panel(ax_c, tables["effects"], "selfing_core", "C")
    _support_panel(ax_d, tables["classification"])

    handles, labels = ax_b.get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=7.5, loc="upper right", bbox_to_anchor=(0.96, 0.53))

    fig.suptitle(
        "The strongest Palearctic floral syndrome has a taxonomic assembly depth",
        x=0.055,
        y=0.985,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.055,
        0.955,
        "Frozen source-matched decomposition: taxonomic attenuation is an estimand of hierarchical expression, not causal mediation.",
        fontsize=9,
        color="0.35",
    )

    basename = "chapter1_v8_figure3_assembly_depth"
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
        "contract": "chapter1_v8_figure3_assembly_depth_v1",
        "status": "rendered_from_frozen_taxonomic_depth_artifacts",
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
        "family_attenuation_fraction_range": [
            float(tables["attenuation"]["family_attenuation_fraction"].min()),
            float(tables["attenuation"]["family_attenuation_fraction"].max()),
        ],
        "genus_attenuation_fraction_range": [
            float(tables["attenuation"]["genus_attenuation_fraction"].min()),
            float(tables["attenuation"]["genus_attenuation_fraction"].max()),
        ],
        "conditional_genus_attenuation_fraction_range": [
            float(tables["attenuation"]["conditional_genus_attenuation_fraction"].min()),
            float(tables["attenuation"]["conditional_genus_attenuation_fraction"].max()),
        ],
        "frozen_support_ladder": "4/4 -> 4/4 -> 0/4",
        "new_biological_models_fitted": False,
        "new_p_values_generated": False,
        "claim_boundary": (
            "Figure 3 localizes hierarchical expression of the frozen Palearctic response. "
            "Attenuation after family/genus composition is not proof of dispersal, absence of "
            "within-lineage evolution, or causal mediation."
        ),
        "outputs": [png_path.name, svg_path.name, pdf_path.name],
    }
    (output_dir / "chapter1_v8_figure3_manifest.json").write_text(
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
    manifest = render_figure(
        effect_root=effect_root,
        pr142_root=pr142_root,
        output_dir=output_dir,
        effect_run_id=effect_run_id,
        effect_artifact_id=effect_artifact_id,
        effect_artifact_digest=effect_artifact_digest,
        primary_run_id=primary_run_id,
        primary_artifact_id=primary_artifact_id,
        primary_artifact_digest=primary_artifact_digest,
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
