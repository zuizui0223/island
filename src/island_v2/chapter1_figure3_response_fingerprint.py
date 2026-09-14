"""Render the Chapter 1 response-fingerprint atlas from the frozen synthesis artifact.

The figure is a presentation layer only. It reads the post-freeze effect-fingerprint
artifact and creates three panels:

A. atomic colour/structure/reproduction response fingerprints;
B. descriptive cross-context vector orientation;
C. Palearctic observed -> family -> genus response-vector attenuation.

No biological model is refit and no new p-value is generated.
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
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec  # noqa: E402

app = typer.Typer(add_completion=False, no_args_is_help=True)

OUTCOME_ORDER = [
    "plain_colour",
    "generalized_form",
    "actinomorphic_symmetry",
    "shallow_open_tube",
    "small_flower",
    "self_compatibility",
    "selfing_mating_system",
    "autonomous_selfing",
]
OUTCOME_LABELS = {
    "plain_colour": "Plain colour",
    "generalized_form": "Generalized form",
    "actinomorphic_symmetry": "Actinomorphic",
    "shallow_open_tube": "Shallow/open tube",
    "small_flower": "Small flower",
    "self_compatibility": "Self-compatible",
    "selfing_mating_system": "Selfing mating",
    "autonomous_selfing": "Autonomous selfing",
}
CONTEXTS = ["northern_midlatitude", "tropical"]
CONTEXT_LABELS = {
    "northern_midlatitude": "Northern mid-latitude",
    "tropical": "Tropical",
}
STRATA = ["all_native", "native_nonendemic"]
STRATUM_LABELS = {
    "all_native": "All native",
    "native_nonendemic": "Native non-endemic",
}
SCOPES = ["all_analysis_eligible", "direct_only"]
SCOPE_LABELS = {
    "all_analysis_eligible": "All-analysis",
    "direct_only": "Direct-only",
}
SCOPE_COLOURS = {
    "all_analysis_eligible": "#3569a8",
    "direct_only": "#d0702f",
}
SCOPE_MARKERS = {
    "all_analysis_eligible": "o",
    "direct_only": "s",
}


class FigureInputError(ValueError):
    """Raised when the frozen figure input does not match the expected schema."""


def _require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    if missing := required - set(frame.columns):
        raise FigureInputError(f"{label} missing columns: {sorted(missing)}")


def load_inputs(input_dir: Path) -> dict[str, pd.DataFrame]:
    """Load and fail-close validate the frozen effect-fingerprint tables."""
    paths = {
        "atomic": input_dir / "atomic_response_fingerprint.csv",
        "angles": input_dir / "atomic_cross_context_vector_geometry.csv",
        "attenuation": input_dir / "taxonomic_vector_attenuation.csv",
    }
    for label, path in paths.items():
        if not path.is_file():
            raise FigureInputError(f"missing {label} input: {path}")

    atomic = pd.read_csv(paths["atomic"])
    angles = pd.read_csv(paths["angles"])
    attenuation = pd.read_csv(paths["attenuation"])

    _require_columns(
        atomic,
        {
            "evidence_scope",
            "stratum",
            "context",
            "outcome",
            "geography_slope_log_odds_per_response_sd",
            "ci_low",
            "ci_high",
            "n_islands",
        },
        "atomic fingerprint",
    )
    _require_columns(
        angles,
        {
            "evidence_scope",
            "stratum",
            "context_a",
            "context_b",
            "n_components",
            "vector_angle_degrees",
            "complete_vector",
        },
        "cross-context geometry",
    )
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

    expected_atomic = {
        (scope, stratum, context, outcome)
        for scope in SCOPES
        for stratum in STRATA
        for context in CONTEXTS
        for outcome in OUTCOME_ORDER
    }
    found_atomic = set(
        atomic[["evidence_scope", "stratum", "context", "outcome"]]
        .astype(str)
        .itertuples(index=False, name=None)
    )
    if found_atomic != expected_atomic:
        missing = sorted(expected_atomic - found_atomic)
        extra = sorted(found_atomic - expected_atomic)
        raise FigureInputError(
            f"atomic fingerprint design cells mismatch; missing={missing[:5]} extra={extra[:5]}"
        )

    expected_angles = {(scope, stratum) for scope in SCOPES for stratum in STRATA}
    found_angles = set(
        angles[["evidence_scope", "stratum"]].astype(str).itertuples(index=False, name=None)
    )
    if found_angles != expected_angles or len(angles) != 4:
        raise FigureInputError("cross-context geometry must contain exactly four scope×stratum cells")
    if not pd.to_numeric(angles["n_components"], errors="coerce").eq(8).all():
        raise FigureInputError("cross-context geometry must use the frozen eight-component vector")

    expected_attenuation = {
        (scope, stratum, source_mode)
        for scope in SCOPES
        for stratum in STRATA
        for source_mode in ("geo_k5", "geo_k10", "geo_k20", "geo50_climate10")
    }
    found_attenuation = set(
        attenuation[["evidence_scope", "stratum", "source_mode"]]
        .astype(str)
        .itertuples(index=False, name=None)
    )
    if found_attenuation != expected_attenuation or len(attenuation) != 16:
        raise FigureInputError("taxonomic attenuation must contain all 16 frozen profiles")
    if set(attenuation["context"].astype(str)) != {"Palearctic"}:
        raise FigureInputError("taxonomic attenuation panel must be the frozen Palearctic response")

    numeric_atomic = ["geography_slope_log_odds_per_response_sd", "ci_low", "ci_high"]
    for column in numeric_atomic:
        atomic[column] = pd.to_numeric(atomic[column], errors="coerce")
    angles["vector_angle_degrees"] = pd.to_numeric(
        angles["vector_angle_degrees"], errors="coerce"
    )
    for column in (
        "observed_score_vector_norm",
        "after_family_residual_vector_norm",
        "after_genus_residual_vector_norm",
        "family_attenuation_fraction",
        "genus_attenuation_fraction",
        "conditional_genus_attenuation_fraction",
    ):
        attenuation[column] = pd.to_numeric(attenuation[column], errors="coerce")

    if atomic[numeric_atomic].isna().any().any():
        raise FigureInputError("atomic fingerprint contains non-numeric plotted values")
    if angles["vector_angle_degrees"].isna().any():
        raise FigureInputError("cross-context geometry contains missing angle values")
    if attenuation[
        [
            "observed_score_vector_norm",
            "after_family_residual_vector_norm",
            "after_genus_residual_vector_norm",
        ]
    ].isna().any().any():
        raise FigureInputError("taxonomic attenuation contains missing vector norms")

    return {"atomic": atomic, "angles": angles, "attenuation": attenuation}


def _range_text(values: pd.Series) -> str:
    x = pd.to_numeric(values, errors="coerce").dropna().to_numpy(float)
    return f"{100 * np.min(x):.1f}–{100 * np.max(x):.1f}%"


def _render_atomic_panel(fig: plt.Figure, slot: Any, atomic: pd.DataFrame) -> list[plt.Axes]:
    nested = GridSpecFromSubplotSpec(2, 2, subplot_spec=slot, hspace=0.18, wspace=0.12)
    axes: list[plt.Axes] = []
    y = np.arange(len(OUTCOME_ORDER))
    for i, context in enumerate(CONTEXTS):
        for j, stratum in enumerate(STRATA):
            ax = fig.add_subplot(nested[i, j])
            axes.append(ax)
            subset = atomic.loc[
                atomic["context"].astype(str).eq(context)
                & atomic["stratum"].astype(str).eq(stratum)
            ]
            for k, scope in enumerate(SCOPES):
                scope_data = (
                    subset.loc[subset["evidence_scope"].astype(str).eq(scope)]
                    .set_index("outcome")
                    .reindex(OUTCOME_ORDER)
                )
                estimate = scope_data["geography_slope_log_odds_per_response_sd"].to_numpy(float)
                low = scope_data["ci_low"].to_numpy(float)
                high = scope_data["ci_high"].to_numpy(float)
                offset = -0.12 if k == 0 else 0.12
                ax.errorbar(
                    estimate,
                    y + offset,
                    xerr=np.vstack([estimate - low, high - estimate]),
                    fmt=SCOPE_MARKERS[scope],
                    markersize=4.5,
                    linewidth=1.1,
                    capsize=2.0,
                    color=SCOPE_COLOURS[scope],
                    label=SCOPE_LABELS[scope],
                    zorder=3,
                )
            ax.axvline(0, color="0.55", linewidth=0.8, linestyle="--", zorder=1)
            for separator in (0.5, 4.5):
                ax.axhline(separator, color="0.87", linewidth=0.7)
            ax.set_yticks(y)
            if j == 0:
                ax.set_yticklabels([OUTCOME_LABELS[outcome] for outcome in OUTCOME_ORDER], fontsize=8.4)
            else:
                ax.set_yticklabels([])
            ax.set_ylim(len(OUTCOME_ORDER) - 0.45, -0.55)
            ax.set_xlim(-0.27, 0.27)
            ax.set_title(
                f"{CONTEXT_LABELS[context]} · {STRATUM_LABELS[stratum]}",
                fontsize=10.2,
                pad=5,
            )
            if i == 1:
                ax.set_xlabel("Isolation slope (log-odds per response SD)", fontsize=9)
            ax.tick_params(axis="x", labelsize=8)
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)
            if i == 0 and j == 1:
                ax.legend(frameon=False, fontsize=8, loc="lower right")

    axes[0].text(-0.33, 1.09, "A", transform=axes[0].transAxes, fontsize=15, fontweight="bold", va="top")
    axes[0].text(
        -0.33,
        1.035,
        "Component response fingerprints",
        transform=axes[0].transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
    )
    for ax in (axes[0], axes[2]):
        ax.text(
            -0.42,
            0.94,
            "COLOUR",
            transform=ax.transAxes,
            fontsize=7.5,
            fontweight="bold",
            rotation=90,
            va="center",
            ha="center",
            color="0.35",
        )
        ax.text(
            -0.42,
            0.61,
            "STRUCTURE",
            transform=ax.transAxes,
            fontsize=7.5,
            fontweight="bold",
            rotation=90,
            va="center",
            ha="center",
            color="0.35",
        )
        ax.text(
            -0.42,
            0.16,
            "REPRODUCTION",
            transform=ax.transAxes,
            fontsize=7.5,
            fontweight="bold",
            rotation=90,
            va="center",
            ha="center",
            color="0.35",
        )
    return axes


def _render_angle_panel(fig: plt.Figure, slot: Any, angles: pd.DataFrame) -> plt.Axes:
    ax = fig.add_subplot(slot)
    rows = [
        ("all_analysis_eligible", "all_native"),
        ("all_analysis_eligible", "native_nonendemic"),
        ("direct_only", "all_native"),
        ("direct_only", "native_nonendemic"),
    ]
    values: list[float] = []
    labels: list[str] = []
    for scope, stratum in rows:
        row = angles.loc[
            angles["evidence_scope"].astype(str).eq(scope)
            & angles["stratum"].astype(str).eq(stratum)
        ].iloc[0]
        values.append(float(row["vector_angle_degrees"]))
        labels.append(f"{SCOPE_LABELS[scope]} · {STRATUM_LABELS[stratum]}")
    y = np.arange(4)
    ax.axvline(90, color="0.65", linestyle="--", linewidth=1)
    ax.hlines(y, 90, values, color="0.75", linewidth=2)
    ax.scatter(values, y, s=38, color="0.2", zorder=3)
    for value, y_value in zip(values, y, strict=True):
        ax.text(value + 1.5, y_value, f"{value:.1f}°", va="center", fontsize=9)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.2)
    ax.set_ylim(3.6, -0.6)
    ax.set_xlim(80, 125)
    ax.set_xlabel("Angle between 8-component vectors", fontsize=9)
    ax.set_title("B  Cross-context orientation", loc="left", fontsize=11, fontweight="bold", pad=8)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="x", labelsize=8)
    ax.text(
        0.02,
        -0.22,
        "90° = orthogonal; >90° = partially opposed",
        transform=ax.transAxes,
        fontsize=8,
        color="0.35",
    )
    return ax


def _render_attenuation_panel(fig: plt.Figure, slot: Any, attenuation: pd.DataFrame) -> plt.Axes:
    ax = fig.add_subplot(slot)
    stages = ["Observed", "After family", "After genus"]
    x = np.arange(3)
    line_styles = {"all_native": "-", "native_nonendemic": "--"}
    normalized_rows: list[np.ndarray] = []
    for row in attenuation.itertuples(index=False):
        observed = float(row.observed_score_vector_norm)
        values = np.asarray(
            [
                1.0,
                float(row.after_family_residual_vector_norm) / observed,
                float(row.after_genus_residual_vector_norm) / observed,
            ]
        )
        normalized_rows.append(values)
        ax.plot(
            x,
            values,
            color=SCOPE_COLOURS[str(row.evidence_scope)],
            linestyle=line_styles[str(row.stratum)],
            linewidth=1.05,
            alpha=0.45,
            marker="o",
            markersize=2.5,
        )
    median = np.median(np.vstack(normalized_rows), axis=0)
    ax.plot(x, median, color="0.1", linewidth=2.6, marker="o", markersize=5, zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels(stages, fontsize=8.5)
    ax.set_ylabel("Response-vector magnitude\n(relative to observed)", fontsize=9)
    ax.set_ylim(0, 1.08)
    ax.set_title("C  Hierarchical attenuation", loc="left", fontsize=11, fontweight="bold", pad=8)
    ax.text(
        0.04,
        0.38,
        f"Family attenuation\n{_range_text(attenuation['family_attenuation_fraction'])}",
        transform=ax.transAxes,
        fontsize=8.5,
    )
    ax.text(
        0.51,
        0.10,
        f"Total genus attenuation\n{_range_text(attenuation['genus_attenuation_fraction'])}",
        transform=ax.transAxes,
        fontsize=8.5,
        fontweight="bold",
    )
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", labelsize=8)
    return ax


def render_figure(
    tables: dict[str, pd.DataFrame],
    output_dir: Path,
    *,
    source_run_id: int,
    source_artifact_id: int,
    source_artifact_digest: str,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    atomic = tables["atomic"]
    angles = tables["angles"]
    attenuation = tables["attenuation"]

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
    _render_atomic_panel(fig, grid[:, 0], atomic)
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

    png_path = output_dir / "chapter1_figure3_response_fingerprint.png"
    svg_path = output_dir / "chapter1_figure3_response_fingerprint.svg"
    pdf_path = output_dir / "chapter1_figure3_response_fingerprint.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    atomic.to_csv(output_dir / "figure3_panel_a_atomic_effects.csv", index=False)
    angles.to_csv(output_dir / "figure3_panel_b_vector_angles.csv", index=False)
    attenuation.to_csv(output_dir / "figure3_panel_c_taxonomic_attenuation.csv", index=False)

    manifest = {
        "contract": "chapter1_figure3_response_fingerprint_v1",
        "status": "figure_rendered_from_frozen_effect_fingerprint",
        "source_workflow_run_id": int(source_run_id),
        "source_artifact_id": int(source_artifact_id),
        "source_artifact_digest": str(source_artifact_digest),
        "panel_a": {
            "contexts": CONTEXTS,
            "strata": STRATA,
            "evidence_scopes": SCOPES,
            "n_atomic_rows": int(len(atomic)),
            "n_outcomes": len(OUTCOME_ORDER),
        },
        "panel_b": {
            "n_vector_angle_rows": int(len(angles)),
            "angle_degree_range": [
                float(angles["vector_angle_degrees"].min()),
                float(angles["vector_angle_degrees"].max()),
            ],
            "inferential_role": "descriptive_geometry_only",
        },
        "panel_c": {
            "n_taxonomic_profiles": int(len(attenuation)),
            "family_attenuation_fraction_range": [
                float(attenuation["family_attenuation_fraction"].min()),
                float(attenuation["family_attenuation_fraction"].max()),
            ],
            "genus_attenuation_fraction_range": [
                float(attenuation["genus_attenuation_fraction"].min()),
                float(attenuation["genus_attenuation_fraction"].max()),
            ],
            "inferential_role": "descriptive_attenuation_not_causal_mediation",
        },
        "new_biological_models_fitted": False,
        "new_p_values_generated": False,
        "claim_boundary": (
            "Figure 3 is a visualization of frozen estimates. Atomic intervals are descriptive "
            "decomposition, vector angles are descriptive geometry, and attenuation is not "
            "causal mediation. Frozen H1-H3 tests remain the inferential basis."
        ),
        "outputs": [png_path.name, svg_path.name, pdf_path.name],
    }
    (output_dir / "chapter1_figure3_response_fingerprint_manifest.json").write_text(
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
    tables = load_inputs(input_dir)
    manifest = render_figure(
        tables,
        output_dir,
        source_run_id=source_run_id,
        source_artifact_id=source_artifact_id,
        source_artifact_digest=source_artifact_digest,
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
