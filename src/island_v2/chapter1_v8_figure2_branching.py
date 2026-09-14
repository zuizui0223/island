"""Render Chapter 1 v8 Figure 2 from the frozen PR142 artifact.

This renderer is presentation-only. It does not refit any model or generate new
inferential quantities. The primary H2 response is kept distinct from the broader
8-component effect fingerprint used in Extended Data.
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
SCOPE_MARKERS = {"all_analysis_eligible": "o", "direct_only": "s"}
STRATA = ["all_native", "native_nonendemic"]
STRATUM_LABELS = {"all_native": "All native", "native_nonendemic": "Native non-endemic"}
REGIME_CONTEXTS = ["northern_midlatitude", "tropical"]
REGIME_LABELS = {"northern_midlatitude": "Northern mid-latitude", "tropical": "Tropical"}
REGIME_COLOURS = {"northern_midlatitude": "#3569a8", "tropical": "#d0702f"}
RESPONSES = ["accessibility_generalization", "reproductive_assurance"]
RESPONSE_LABELS = {
    "accessibility_generalization": "Accessibility / generalization",
    "reproductive_assurance": "Reproductive assurance",
}
RESPONSE_COLOURS = {
    "accessibility_generalization": "#3569a8",
    "reproductive_assurance": "#8f4c97",
}
SOURCE_MODES = ["geo_k5", "geo_k10", "geo_k20", "geo50_climate10"]
SOURCE_LABELS = {
    "geo_k5": "geo k=5",
    "geo_k10": "geo k=10",
    "geo_k20": "geo k=20",
    "geo50_climate10": "geo50 + climate10",
}


class FigureInputError(ValueError):
    """Raised when a frozen artifact does not match the expected figure contract."""


def _require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    if missing := required - set(frame.columns):
        raise FigureInputError(f"{label} missing columns: {sorted(missing)}")


def load_inputs(artifact_root: Path) -> dict[str, pd.DataFrame]:
    """Load and validate only the frozen rows needed for v8 Figure 2."""
    slope_frames: list[pd.DataFrame] = []
    for short, scope in (("all", "all_analysis_eligible"), ("direct", "direct_only")):
        path = artifact_root / f"branching/{short}/global_branch_distance_slopes.csv"
        if not path.is_file():
            raise FigureInputError(f"missing primary branching table: {path}")
        frame = pd.read_csv(path)
        _require_columns(
            frame,
            {
                "context_layer",
                "axis_set",
                "stratum",
                "support_tier",
                "context",
                "syndrome",
                "distance_slope",
                "cluster_robust_se",
                "q_axis_family",
                "n_islands",
            },
            f"{scope} branching slopes",
        )
        frame = frame.loc[
            frame["axis_set"].astype(str).eq("universal_plant_response")
            & frame["support_tier"].astype(str).eq("confirmatory")
            & frame["stratum"].astype(str).isin(STRATA)
            & frame["syndrome"].astype(str).isin(RESPONSES)
        ].copy()
        frame.insert(0, "evidence_scope", scope)
        slope_frames.append(frame)
    slopes = pd.concat(slope_frames, ignore_index=True)

    regime = slopes.loc[
        slopes["context_layer"].astype(str).eq("analysis_regime")
        & slopes["context"].astype(str).isin(REGIME_CONTEXTS)
    ].copy()
    expected_regime = {
        (scope, stratum, context, response)
        for scope in SCOPES
        for stratum in STRATA
        for context in REGIME_CONTEXTS
        for response in RESPONSES
    }
    found_regime = set(
        regime[["evidence_scope", "stratum", "context", "syndrome"]]
        .astype(str)
        .itertuples(index=False, name=None)
    )
    if found_regime != expected_regime:
        raise FigureInputError("analysis-regime primary response cells are incomplete")

    palearctic = slopes.loc[
        slopes["context_layer"].astype(str).eq("biogeographic_realm")
        & slopes["context"].astype(str).eq("Palearctic")
    ].copy()
    expected_palearctic = {
        (scope, stratum, response)
        for scope in SCOPES
        for stratum in STRATA
        for response in RESPONSES
    }
    found_palearctic = set(
        palearctic[["evidence_scope", "stratum", "syndrome"]]
        .astype(str)
        .itertuples(index=False, name=None)
    )
    if found_palearctic != expected_palearctic:
        raise FigureInputError("Palearctic primary response cells are incomplete")

    pathway_path = artifact_root / "source-adjusted/source_adjusted_pathway_realm.csv"
    if not pathway_path.is_file():
        raise FigureInputError(f"missing frozen pathway table: {pathway_path}")
    pathway = pd.read_csv(pathway_path)
    _require_columns(
        pathway,
        {
            "source_mode",
            "stratum",
            "context",
            "support_tier",
            "model",
            "status",
            "distance_estimate",
            "distance_se",
            "distance_q",
            "n_unique_islands",
        },
        "source-adjusted pathway",
    )
    pathway = pathway.loc[
        pathway["context"].astype(str).eq("Palearctic")
        & pathway["support_tier"].astype(str).eq("confirmatory")
        & pathway["model"].astype(str).eq("attraction_conditional_on_selfing_core")
        & pathway["status"].astype(str).eq("fit")
        & pathway["source_mode"].astype(str).isin(SOURCE_MODES)
        & pathway["stratum"].astype(str).isin(STRATA)
    ].copy()
    expected_pathway = {(mode, stratum) for mode in SOURCE_MODES for stratum in STRATA}
    found_pathway = set(
        pathway[["source_mode", "stratum"]].astype(str).itertuples(index=False, name=None)
    )
    if found_pathway != expected_pathway:
        raise FigureInputError("Palearctic conditional-pathway source-mode cells are incomplete")

    for frame, columns in (
        (regime, ["distance_slope", "cluster_robust_se"]),
        (palearctic, ["distance_slope", "cluster_robust_se"]),
        (pathway, ["distance_estimate", "distance_se"]),
    ):
        for column in columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
        if frame[columns].isna().any().any():
            raise FigureInputError("non-numeric plotted values in frozen Figure 2 input")

    return {"regime": regime, "palearctic": palearctic, "pathway": pathway}


def _response_row(frame: pd.DataFrame, response: str) -> pd.Series:
    rows = frame.loc[frame["syndrome"].astype(str).eq(response)]
    if len(rows) != 1:
        raise FigureInputError(f"expected one row for {response}, found {len(rows)}")
    return rows.iloc[0]


def _phase_panel(ax: plt.Axes, regime: pd.DataFrame, stratum: str, letter: str) -> None:
    subset = regime.loc[regime["stratum"].astype(str).eq(stratum)]
    ax.axhline(0, color="0.75", linewidth=0.8)
    ax.axvline(0, color="0.75", linewidth=0.8)
    for context in REGIME_CONTEXTS:
        for scope in SCOPES:
            cell = subset.loc[
                subset["context"].astype(str).eq(context)
                & subset["evidence_scope"].astype(str).eq(scope)
            ]
            xrow = _response_row(cell, "accessibility_generalization")
            yrow = _response_row(cell, "reproductive_assurance")
            x = float(xrow["distance_slope"])
            y = float(yrow["distance_slope"])
            xse = float(xrow["cluster_robust_se"])
            yse = float(yrow["cluster_robust_se"])
            ax.errorbar(
                x,
                y,
                xerr=1.96 * xse,
                yerr=1.96 * yse,
                fmt=SCOPE_MARKERS[scope],
                markersize=6,
                capsize=2.5,
                linewidth=1.0,
                color=REGIME_COLOURS[context],
                markerfacecolor=REGIME_COLOURS[context] if scope == "direct_only" else "white",
                markeredgewidth=1.4,
                zorder=3,
            )
    ax.set_xlabel("Accessibility/generalization slope", fontsize=9)
    ax.set_ylabel("Reproductive-assurance slope", fontsize=9)
    ax.set_title(f"{letter}  {STRATUM_LABELS[stratum]}", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(labelsize=8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.text(
        0.03,
        0.97,
        "Open = all-analysis\nFilled = direct-only",
        transform=ax.transAxes,
        va="top",
        fontsize=7.8,
        color="0.35",
    )
    for context, ypos in (("northern_midlatitude", 0.17), ("tropical", 0.09)):
        ax.text(
            0.98,
            ypos,
            REGIME_LABELS[context],
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
            color=REGIME_COLOURS[context],
            fontweight="bold",
        )


def _palearctic_panel(ax: plt.Axes, palearctic: pd.DataFrame) -> None:
    rows = [(scope, stratum) for stratum in STRATA for scope in SCOPES]
    ybase = np.arange(len(rows))
    offsets = {"accessibility_generalization": -0.12, "reproductive_assurance": 0.12}
    for response in RESPONSES:
        estimates = []
        lows = []
        highs = []
        for scope, stratum in rows:
            row = _response_row(
                palearctic.loc[
                    palearctic["evidence_scope"].astype(str).eq(scope)
                    & palearctic["stratum"].astype(str).eq(stratum)
                ],
                response,
            )
            est = float(row["distance_slope"])
            se = float(row["cluster_robust_se"])
            estimates.append(est)
            lows.append(est - 1.96 * se)
            highs.append(est + 1.96 * se)
        estimates_arr = np.asarray(estimates)
        ax.errorbar(
            estimates_arr,
            ybase + offsets[response],
            xerr=np.vstack([estimates_arr - np.asarray(lows), np.asarray(highs) - estimates_arr]),
            fmt="o",
            markersize=4.5,
            capsize=2,
            linewidth=1,
            color=RESPONSE_COLOURS[response],
            label=RESPONSE_LABELS[response],
        )
    ax.axvline(0, color="0.65", linestyle="--", linewidth=0.8)
    labels = [f"{SCOPE_LABELS[s]} · {STRATUM_LABELS[t]}" for s, t in rows]
    ax.set_yticks(ybase)
    ax.set_yticklabels(labels, fontsize=7.8)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xlabel("Palearctic isolation slope", fontsize=9)
    ax.set_title("C  Palearctic replication", loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=7.7, loc="lower right")
    ax.tick_params(axis="x", labelsize=8)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)


def _conditional_panel(ax: plt.Axes, pathway: pd.DataFrame) -> None:
    y = np.arange(len(SOURCE_MODES))
    for stratum, offset, marker in (("all_native", -0.10, "o"), ("native_nonendemic", 0.10, "s")):
        subset = pathway.loc[pathway["stratum"].astype(str).eq(stratum)].set_index("source_mode")
        estimates = np.asarray([float(subset.loc[m, "distance_estimate"]) for m in SOURCE_MODES])
        ses = np.asarray([float(subset.loc[m, "distance_se"]) for m in SOURCE_MODES])
        ax.errorbar(
            estimates,
            y + offset,
            xerr=1.96 * ses,
            fmt=marker,
            markersize=4.5,
            capsize=2,
            linewidth=1,
            label=STRATUM_LABELS[stratum],
        )
    ax.axvline(0, color="0.65", linestyle="--", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([SOURCE_LABELS[m] for m in SOURCE_MODES], fontsize=7.8)
    ax.set_ylim(len(SOURCE_MODES) - 0.5, -0.5)
    ax.set_xlabel("Distance effect after conditioning on selfing core", fontsize=8.6)
    ax.set_title("D  Floral architecture remains after selfing-core conditioning", loc="left", fontsize=10.3, fontweight="bold")
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")
    ax.tick_params(axis="x", labelsize=8)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)


def render_figure(
    *,
    artifact_root: Path,
    output_dir: Path,
    source_run_id: int,
    source_artifact_id: int,
    source_artifact_digest: str,
) -> dict[str, Any]:
    tables = load_inputs(artifact_root)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(13.2, 9.4), constrained_layout=False)
    grid = GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.34)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    _phase_panel(ax_a, tables["regime"], "all_native", "A")
    _phase_panel(ax_b, tables["regime"], "native_nonendemic", "B")
    _palearctic_panel(ax_c, tables["palearctic"])
    _conditional_panel(ax_d, tables["pathway"])

    fig.suptitle(
        "Floral and reproductive responses branch among biogeographic contexts",
        x=0.06,
        y=0.985,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.06,
        0.955,
        "Primary H2 responses are pollinator-name-free; intervals are cluster-robust 95% CIs from frozen estimates.",
        fontsize=9,
        color="0.35",
    )

    basename = "chapter1_v8_figure2_biogeographic_branching"
    png_path = output_dir / f"{basename}.png"
    svg_path = output_dir / f"{basename}.svg"
    pdf_path = output_dir / f"{basename}.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    tables["regime"].to_csv(output_dir / "figure2_regime_primary_slopes.csv", index=False)
    tables["palearctic"].to_csv(output_dir / "figure2_palearctic_primary_slopes.csv", index=False)
    tables["pathway"].to_csv(output_dir / "figure2_palearctic_conditional_attraction.csv", index=False)

    manifest = {
        "contract": "chapter1_v8_figure2_biogeographic_branching_v1",
        "status": "rendered_from_frozen_pr142_estimates",
        "source_workflow_run_id": int(source_run_id),
        "source_artifact_id": int(source_artifact_id),
        "source_artifact_digest": str(source_artifact_digest),
        "new_biological_models_fitted": False,
        "new_p_values_generated": False,
        "n_regime_primary_rows": int(len(tables["regime"])),
        "n_palearctic_primary_rows": int(len(tables["palearctic"])),
        "n_conditional_pathway_rows": int(len(tables["pathway"])),
        "claim_boundary": (
            "Figure 2 visualizes frozen H2 and pathway-decomposition estimates only. "
            "The phase panels do not create a new between-context test, and the conditional "
            "pathway panel is decomposition rather than causal mediation."
        ),
        "outputs": [png_path.name, svg_path.name, pdf_path.name],
    }
    (output_dir / "chapter1_v8_figure2_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("render")
def render_command(
    artifact_root: Path = typer.Option(..., exists=True, file_okay=False),
    output_dir: Path = typer.Option(...),
    source_run_id: int = typer.Option(...),
    source_artifact_id: int = typer.Option(...),
    source_artifact_digest: str = typer.Option(...),
) -> None:
    manifest = render_figure(
        artifact_root=artifact_root,
        output_dir=output_dir,
        source_run_id=source_run_id,
        source_artifact_id=source_artifact_id,
        source_artifact_digest=source_artifact_digest,
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
