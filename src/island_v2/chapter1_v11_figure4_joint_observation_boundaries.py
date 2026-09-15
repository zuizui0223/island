"""Render P3-defended Chapter 1 Figure 4 from frozen joint observation-bias evidence.

Presentation only. The joint surface has already been fitted and locked. This renderer
does not refit biological models or generate new p-values.
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
    envelope = pd.read_csv(_require(joint_root / "partial_identification_envelope.csv"))
    manifest = json.loads(
        _require(joint_root / "joint_observation_bias_manifest.json").read_text()
    )
    geom = pd.read_csv(_require(geometry_root / "observed_geometry_cross_scope.csv"))
    h5c = json.loads(_require(h5c_root / "h5c_observed_result.json").read_text())
    h5d = pd.read_csv(_require(h5d_root / "h5d_identifiability_summary.csv"))

    if manifest.get("contract") != "chapter1_joint_observation_bias_v1":
        raise FigureInputError("unexpected P3 joint contract")
    if int(manifest.get("n_primary_parameter_surfaces_per_scope", -1)) != 1575:
        raise FigureInputError("unexpected P3 primary surface size")
    if bool(manifest.get("grid_fraction_is_probability", True)):
        raise FigureInputError("grid fractions must not be interpreted as probabilities")

    required_surface = {
        "evidence_scope",
        "target",
        "stratum",
        "trait_resolution_odds_ratio",
        "median_distance_completeness",
        "distance_completeness_odds_ratio",
        "state_recording_odds_ratio",
        "robust_cell",
    }
    required_envelope = {
        "evidence_scope",
        "target",
        "stratum",
        "estimate_lower",
        "estimate_upper",
        "expected_sign_identified",
        "support_identified_across_envelope",
    }
    if missing := required_surface - set(surface.columns):
        raise FigureInputError(f"joint surface missing columns: {sorted(missing)}")
    if missing := required_envelope - set(envelope.columns):
        raise FigureInputError(f"partial envelope missing columns: {sorted(missing)}")

    if len(geom) != 12 or set(geom["classification"].astype(str)) != {
        "monotonic_or_unresolved"
    }:
        raise FigureInputError("geometry result differs from frozen 0/12 promotion")
    if h5c.get("classification") != "no_pollination_mode_specificity_support":
        raise FigureInputError("H5c classification differs from frozen result")
    if len(h5d) != 8 or pd.Series(h5d["qualified"]).astype(bool).any():
        raise FigureInputError("H5d result differs from frozen 0/8 qualification")

    return {
        "surface": surface,
        "envelope": envelope,
        "manifest": manifest,
        "geometry": geom,
        "h5c": h5c,
        "h5d": h5d,
    }


def _robust_matrix(
    surface: pd.DataFrame,
    *,
    target: str,
    scope: str = "direct_only",
    stratum: str = "native_nonendemic",
) -> tuple[np.ndarray, list[float], list[float], tuple[int, int]]:
    sub = surface.loc[
        surface["evidence_scope"].astype(str).eq(scope)
        & surface["target"].astype(str).eq(target)
        & surface["stratum"].astype(str).eq(stratum)
    ].copy()
    if len(sub) != 1575:
        raise FigureInputError(f"{target} expected 1575 cells, found {len(sub)}")
    sub["robust_numeric"] = sub["robust_cell"].astype(bool).astype(float)
    table = sub.pivot_table(
        index="trait_resolution_odds_ratio",
        columns="state_recording_odds_ratio",
        values="robust_numeric",
        aggfunc="mean",
    ).sort_index()
    table = table.reindex(sorted(table.columns, reverse=True), axis=1)
    robust = int(sub["robust_cell"].astype(bool).sum())
    return (
        table.to_numpy(float),
        [float(x) for x in table.columns],
        [float(x) for x in table.index],
        (robust, len(sub)),
    )


def _draw_heatmap(
    ax: plt.Axes,
    matrix: np.ndarray,
    xvals: list[float],
    yvals: list[float],
    *,
    title: str,
    summary: tuple[int, int],
    show_ylabel: bool = True,
) -> Any:
    im = ax.imshow(matrix, aspect="auto", vmin=0, vmax=1, origin="lower")
    ax.set_xticks(np.arange(len(xvals)))
    ax.set_xticklabels(
        [f"{x:g}" for x in xvals], fontsize=7.0, rotation=45, ha="right"
    )
    ax.set_yticks(np.arange(len(yvals)))
    ax.set_yticklabels([f"{y:g}" for y in yvals], fontsize=7.2)
    ax.set_xlabel("species recording OR_D", fontsize=8)
    if show_ylabel:
        ax.set_ylabel("trait-resolution OR_R", fontsize=8)
    else:
        ax.set_ylabel("")
    ax.set_title(title, loc="left", fontsize=9.2, fontweight="bold")
    kept, total = summary
    ax.text(
        0.02,
        0.98,
        f"{kept}/{total} robust",
        transform=ax.transAxes,
        va="top",
        fontsize=8.3,
        fontweight="bold",
    )
    return im


def _panel_formal(ax: plt.Axes, surface: pd.DataFrame) -> tuple[Any, tuple[int, int]]:
    matrix, xvals, yvals, summary = _robust_matrix(
        surface, target="north_tropical_vector_difference"
    )
    im = _draw_heatmap(
        ax,
        matrix,
        xvals,
        yvals,
        title="A  Formal North–Tropical vector contrast",
        summary=summary,
    )
    ax.text(
        0.02,
        0.04,
        "Cells average over the frozen C0 × OR_C grid.\n"
        "Fraction is sensitivity-domain coverage, not probability.",
        transform=ax.transAxes,
        fontsize=7.2,
        va="bottom",
    )
    return im, summary


def _panel_contexts(
    spec: GridSpecFromSubplotSpec,
    fig: plt.Figure,
    surface: pd.DataFrame,
) -> tuple[Any, dict[str, tuple[int, int]]]:
    ax_p = fig.add_subplot(spec[0, 0])
    ax_t = fig.add_subplot(spec[0, 1])
    p = _robust_matrix(surface, target="Palearctic_accessibility")
    t = _robust_matrix(surface, target="tropical_accessibility")
    im = _draw_heatmap(
        ax_p,
        p[0],
        p[1],
        p[2],
        title="B1  Palearctic accessibility",
        summary=p[3],
    )
    _draw_heatmap(
        ax_t,
        t[0],
        t[1],
        t[2],
        title="B2  Tropical accessibility",
        summary=t[3],
        show_ylabel=False,
    )
    ax_p.text(
        -0.18,
        1.08,
        "B  Same bias domain, opposite accessibility robustness",
        transform=ax_p.transAxes,
        fontsize=9.6,
        fontweight="bold",
    )
    return im, {"Palearctic": p[3], "Tropical": t[3]}


def _panel_envelope(ax: plt.Axes, envelope: pd.DataFrame) -> list[dict[str, Any]]:
    wanted = []
    labels = []
    for target, short in (
        ("Palearctic_accessibility", "Palearctic"),
        ("tropical_accessibility", "Tropical"),
    ):
        for scope, scope_label in (
            ("all_analysis_eligible", "all"),
            ("direct_only", "direct"),
        ):
            sub = envelope.loc[
                envelope["target"].astype(str).eq(target)
                & envelope["evidence_scope"].astype(str).eq(scope)
                & envelope["stratum"].astype(str).eq("native_nonendemic")
            ]
            if len(sub) != 1:
                raise FigureInputError(f"missing envelope row {target} {scope}")
            row = sub.iloc[0]
            wanted.append(row)
            labels.append(f"{short} · {scope_label}")

    y = np.arange(len(wanted))[::-1]
    ax.axvline(0, linestyle="--", linewidth=0.9)
    result_rows: list[dict[str, Any]] = []
    for yi, row, label in zip(y, wanted, labels, strict=True):
        lo = float(row["estimate_lower"])
        hi = float(row["estimate_upper"])
        mid = (lo + hi) / 2.0
        ax.errorbar(
            [mid],
            [yi],
            xerr=[[mid - lo], [hi - mid]],
            fmt="o",
            capsize=3,
            markersize=5,
        )
        support = bool(row["support_identified_across_envelope"])
        sign = bool(row["expected_sign_identified"])
        ax.text(
            hi + 0.006,
            yi,
            f"sign {'yes' if sign else 'no'}; support {'yes' if support else 'no'}",
            va="center",
            fontsize=7.4,
        )
        result_rows.append(
            {
                "label": label,
                "lower": lo,
                "upper": hi,
                "sign_identified": sign,
                "support_identified": support,
            }
        )
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Distance-slope partial-identification envelope", fontsize=8.2)
    ax.set_title(
        "C  Deterministic bounds separate identified from fragile axes",
        loc="left",
        fontsize=10.5,
        fontweight="bold",
    )
    ax.text(
        0.0,
        -0.23,
        "Formal North–Tropical vector: support not identified across all deterministic corners.",
        transform=ax.transAxes,
        fontsize=7.2,
        va="top",
    )
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="x", labelsize=7.8)
    return result_rows


def _panel_boundaries(
    ax: plt.Axes,
    geom: pd.DataFrame,
    h5c: dict[str, Any],
    h5d: pd.DataFrame,
) -> dict[str, Any]:
    ax.axis("off")
    geom_promoted = int(
        (geom["classification"].astype(str) != "monotonic_or_unresolved").sum()
    )
    h5d_qualified = int(pd.Series(h5d["qualified"]).astype(bool).sum())
    h5c_p = float(h5c["interaction_p_value"])
    lines = [
        "D  Mechanistic claim boundary remains closed",
        "",
        f"Nonlinear response geometry: {geom_promoted}/12 promoted",
        f"Independent biotic-vs-wind specificity: p = {h5c_p:.3f}",
        f"Distributed-threshold identifiability: {h5d_qualified}/8 qualified",
        "",
        "P3 localizes observation robustness.",
        "It does not estimate true completeness,",
        "identify arbitrary MNAR truth, or rescue H5.",
    ]
    ax.text(
        0.02,
        0.98,
        "\n".join(lines),
        transform=ax.transAxes,
        va="top",
        fontsize=9.2,
        linespacing=1.45,
    )
    return {
        "geometry_promoted": geom_promoted,
        "h5c_p_value": h5c_p,
        "h5d_qualified": h5d_qualified,
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
    data = load_inputs(joint_root, geometry_root, h5c_root, h5d_root)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(13.6, 9.2), constrained_layout=False)
    grid = GridSpec(
        2,
        2,
        figure=fig,
        left=0.08,
        right=0.97,
        bottom=0.09,
        top=0.87,
        hspace=0.42,
        wspace=0.30,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    im_a, formal_summary = _panel_formal(ax_a, data["surface"])
    nested = GridSpecFromSubplotSpec(1, 2, subplot_spec=grid[0, 1], wspace=0.38)
    _im_b, context_summary = _panel_contexts(nested, fig, data["surface"])
    ax_c = fig.add_subplot(grid[1, 0])
    envelope_rows = _panel_envelope(ax_c, data["envelope"])
    ax_d = fig.add_subplot(grid[1, 1])
    boundaries = _panel_boundaries(
        ax_d, data["geometry"], data["h5c"], data["h5d"]
    )

    context_axes = [axis for axis in fig.axes if axis not in (ax_a, ax_c, ax_d)]
    cbar = fig.colorbar(
        im_a,
        ax=[ax_a, *context_axes],
        fraction=0.018,
        pad=0.02,
    )
    cbar.set_label("Fraction robust across C0 × OR_C", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    fig.suptitle(
        "Joint observation bias localizes robust and fragile floral-island claims",
        x=0.08,
        y=0.965,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.08,
        0.925,
        "P3 combines frozen V5 trait-resolution MNAR and V6 species-list detection "
        "assumptions without adding precision for hypothetical species.",
        fontsize=8.8,
    )

    basename = "chapter1_v11_figure4_joint_observation_boundaries"
    outputs = []
    for ext, kwargs in (
        ("png", {"dpi": 300}),
        ("svg", {}),
        ("pdf", {}),
    ):
        path = output_dir / f"{basename}.{ext}"
        fig.savefig(path, bbox_inches="tight", **kwargs)
        outputs.append(path.name)
    plt.close(fig)

    manifest = {
        "contract": "chapter1_v11_figure4_joint_observation_boundaries_v1",
        "status": "rendered_from_frozen_P3_and_prior_claim_boundaries",
        "sources": source_receipts,
        "formal_direct_native_nonendemic": {
            "robust": formal_summary[0],
            "total": formal_summary[1],
        },
        "context_direct_native_nonendemic": {
            key: {"robust": value[0], "total": value[1]}
            for key, value in context_summary.items()
        },
        "partial_identification": envelope_rows,
        "claim_boundaries": boundaries,
        "grid_fraction_is_probability": False,
        "new_biological_models_fitted": False,
        "new_p_values_generated": False,
        "claim_boundary": (
            "Figure 4 distinguishes finite-grid robustness from partial identification. "
            "Palearctic accessibility is the observation-robust core; the formal "
            "North–Tropical contrast is highly finite-grid robust but not identified "
            "across all deterministic corners; tropical accessibility is observation-fragile. "
            "No pollination mechanism or latent true completeness is identified."
        ),
        "outputs": outputs,
    }
    (output_dir / "chapter1_v11_figure4_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
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
    receipts = json.loads(source_receipts_json.read_text())
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
