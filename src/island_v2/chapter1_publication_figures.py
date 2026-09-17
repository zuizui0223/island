from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import typer

from island_v2 import _chapter1_publication_figures_base as base

APP = typer.Typer(add_completion=False)

PUBLICATION_LABELS = base.PUBLICATION_LABELS
CONTEXT_ORDER = base.CONTEXT_ORDER
CONTEXT_LABELS = base.CONTEXT_LABELS
CONTEXT_COLORS = base.CONTEXT_COLORS
CONTEXT_MARKERS = base.CONTEXT_MARKERS
TRAIT_ORDER = base.TRAIT_ORDER
TRAIT_LABELS = base.TRAIT_LABELS
SCOPE_LABELS = base.SCOPE_LABELS
FIGURE_WIDTH_IN = base.FIGURE_WIDTH_IN


def _forest_atomic_clean(
    ax,
    atomic: pd.DataFrame,
    omnibus: pd.DataFrame,
    scope: str,
    xlim: tuple[float, float],
) -> None:
    data = atomic.loc[atomic["evidence_scope"].eq(scope)].copy()
    offsets = dict(zip(CONTEXT_ORDER, (-0.27, -0.09, 0.09, 0.27), strict=True))
    y_base = {trait: len(TRAIT_ORDER) - 1 - i for i, trait in enumerate(TRAIT_ORDER)}
    for context in CONTEXT_ORDER:
        part = data.loc[data["context"].eq(context)].set_index("outcome")
        ys = np.array([y_base[trait] + offsets[context] for trait in TRAIT_ORDER])
        values = np.array([part.loc[trait, "estimate"] for trait in TRAIT_ORDER], dtype=float)
        lows = np.array([part.loc[trait, "ci_low"] for trait in TRAIT_ORDER], dtype=float)
        highs = np.array([part.loc[trait, "ci_high"] for trait in TRAIT_ORDER], dtype=float)
        ax.errorbar(
            values,
            ys,
            xerr=np.vstack([values - lows, highs - values]),
            fmt=CONTEXT_MARKERS[context],
            markersize=3.0,
            color=CONTEXT_COLORS[context],
            ecolor=CONTEXT_COLORS[context],
            elinewidth=0.55,
            capsize=1.3,
            markeredgewidth=0,
        )
    ax.axvline(0, color="#777777", lw=0.5, ls="--")
    ax.set_xlim(*xlim)
    ax.set_yticks([y_base[t] for t in TRAIT_ORDER], [TRAIT_LABELS[t] for t in TRAIT_ORDER])
    ax.set_xlabel("Standardized isolation effect (log-odds scale)")
    ax.grid(axis="x", color="#e6e6e6", lw=0.35)
    ax.spines[["top", "right", "left"]].set_visible(False)

    support = omnibus.loc[omnibus["evidence_scope"].eq(scope)].set_index("context")
    pairs = []
    for context in CONTEXT_ORDER:
        pairs.append(
            f"{CONTEXT_LABELS[context]}: n={int(support.loc[context, 'n_unique_islands']):,}, "
            f"P={support.loc[context, 'p_value']:.2g}"
        )
    support_text = f"{pairs[0]}   |   {pairs[1]}\n{pairs[2]}   |   {pairs[3]}"
    ax.text(
        0.0,
        -0.18,
        support_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=4.15,
        linespacing=1.35,
    )


def _figure2_clean(source_dir: Path, output_dir: Path) -> list[Path]:
    atomic = pd.read_csv(source_dir / "figure2_atomic_coefficients.csv")
    omnibus = pd.read_csv(source_dir / "figure2_omnibus_support.csv")
    summaries = pd.read_csv(source_dir / "figure2_descriptive_summaries.csv")
    max_abs = float(np.nanmax(np.abs(atomic[["ci_low", "ci_high"]].to_numpy())))
    xlim = (-max_abs * 1.08, max_abs * 1.08)

    fig, axes = base.plt.subplots(2, 2, figsize=(FIGURE_WIDTH_IN, 5.65))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()
    _forest_atomic_clean(ax_a, atomic, omnibus, "all_analysis", xlim)
    _forest_atomic_clean(ax_b, atomic, omnibus, "direct_only", xlim)
    ax_a.set_title("Primary evidence", loc="left", fontweight="bold")
    ax_b.set_title("Direct-only sensitivity", loc="left", fontweight="bold")
    base._panel_label(ax_a, "a")
    base._panel_label(ax_b, "b")
    ax_b.set_yticklabels([])

    handles = [
        base.Line2D(
            [0],
            [0],
            marker=CONTEXT_MARKERS[c],
            linestyle="",
            color=CONTEXT_COLORS[c],
            label=CONTEXT_LABELS[c],
            markersize=3.5,
        )
        for c in CONTEXT_ORDER
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.992),
        ncol=4,
        frameon=False,
        columnspacing=0.9,
        handletextpad=0.3,
    )

    classic = summaries.loc[summaries["summary"].eq("classic_orientation")].copy()
    y = np.arange(len(CONTEXT_ORDER))[::-1]
    for scope, marker, fill in (
        ("all_analysis", "o", True),
        ("direct_only", "s", False),
    ):
        part = classic.loc[classic["evidence_scope"].eq(scope)].set_index("context")
        vals = [part.loc[c, "estimate"] for c in CONTEXT_ORDER]
        ax_c.scatter(
            vals,
            y,
            marker=marker,
            s=22,
            facecolors="#0072B2" if fill else "white",
            edgecolors="#0072B2",
            linewidths=0.7,
            label=SCOPE_LABELS[scope],
        )
    ax_c.axvline(0, color="#777777", lw=0.5, ls="--")
    ax_c.set_yticks(y, [CONTEXT_LABELS[c] for c in CONTEXT_ORDER])
    ax_c.set_xlabel("Descriptive mean classic-island direction")
    ax_c.set_title("Recurrent syndrome orientation", loc="left", fontweight="bold")
    ax_c.legend(frameon=False, loc="lower right")
    ax_c.spines[["top", "right", "left"]].set_visible(False)
    ax_c.grid(axis="x", color="#e6e6e6", lw=0.35)
    base._panel_label(ax_c, "c")

    family = summaries.loc[
        summaries["summary"].isin(["reproductive_assurance", "floral_accessibility"])
    ].copy()
    positions = {c: i for i, c in enumerate(CONTEXT_ORDER[::-1])}
    for summary_name, color, offset in (
        ("reproductive_assurance", "#0072B2", -0.12),
        ("floral_accessibility", "#E69F00", 0.12),
    ):
        for scope, marker, open_marker in (
            ("all_analysis", "o", False),
            ("direct_only", "s", True),
        ):
            part = family.loc[
                (family["summary"].eq(summary_name))
                & (family["evidence_scope"].eq(scope))
            ]
            xs = part["estimate"].to_numpy(dtype=float)
            ys = np.array([positions[c] + offset for c in part["context"]])
            ax_d.scatter(
                xs,
                ys,
                marker=marker,
                s=19,
                facecolors="white" if open_marker else color,
                edgecolors=color,
                linewidths=0.7,
            )
    ax_d.axvline(0, color="#777777", lw=0.5, ls="--")
    ax_d.set_yticks(range(4), [CONTEXT_LABELS[c] for c in CONTEXT_ORDER[::-1]])
    ax_d.set_xlabel("Descriptive family mean isolation effect")
    ax_d.set_title("Two positive response families", loc="left", fontweight="bold")
    ax_d.spines[["top", "right", "left"]].set_visible(False)
    ax_d.grid(axis="x", color="#e6e6e6", lw=0.35)
    family_handles = [
        base.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            color="#0072B2",
            label="Reproductive assurance",
            markersize=3.5,
        ),
        base.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            color="#E69F00",
            label="Floral accessibility / generalization",
            markersize=3.5,
        ),
        base.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor="#555555",
            markeredgecolor="#555555",
            color="#555555",
            label="Primary",
            markersize=3.5,
        ),
        base.Line2D(
            [0],
            [0],
            marker="s",
            linestyle="",
            markerfacecolor="white",
            markeredgecolor="#555555",
            color="#555555",
            label="Direct-only",
            markersize=3.5,
        ),
    ]
    ax_d.legend(handles=family_handles, frameon=False, loc="lower right", fontsize=4.9)
    base._panel_label(ax_d, "d")

    fig.subplots_adjust(
        left=0.17,
        right=0.99,
        top=0.91,
        bottom=0.10,
        wspace=0.23,
        hspace=0.53,
    )
    return base._save(fig, output_dir / "main", "figure2_recurrent_plant_response")


def _extended1_clean(source_dir: Path, output_dir: Path) -> list[Path]:
    atomic = pd.read_csv(source_dir / "figure2_atomic_coefficients.csv")
    fig, axes = base.plt.subplots(1, 2, figsize=(FIGURE_WIDTH_IN, 3.45), sharey=True)
    vmax = float(pd.to_numeric(atomic["n_islands"], errors="coerce").max())
    for ax, scope, title in zip(
        axes,
        ("all_analysis", "direct_only"),
        ("Primary", "Direct-only"),
        strict=True,
    ):
        data = atomic.loc[atomic["evidence_scope"].eq(scope)]
        pivot = data.pivot(index="outcome", columns="context", values="n_islands").reindex(
            index=TRAIT_ORDER, columns=CONTEXT_ORDER
        )
        ax.imshow(pivot.to_numpy(), cmap="Greys", aspect="auto", vmin=0, vmax=vmax)
        for yi in range(len(TRAIT_ORDER)):
            for xi in range(len(CONTEXT_ORDER)):
                value = float(pivot.iloc[yi, xi])
                text_color = "white" if value >= 0.58 * vmax else "#111111"
                ax.text(
                    xi,
                    yi,
                    f"{int(value):,}",
                    ha="center",
                    va="center",
                    fontsize=4.8,
                    color=text_color,
                )
        ax.set_xticks(
            range(4),
            [CONTEXT_LABELS[c] for c in CONTEXT_ORDER],
            rotation=35,
            ha="right",
        )
        ax.set_yticks(range(6), [TRAIT_LABELS[t] for t in TRAIT_ORDER])
        ax.set_title(title, loc="left", fontweight="bold")
    base._panel_label(axes[0], "a")
    base._panel_label(axes[1], "b")
    fig.text(
        0.98,
        0.955,
        "Cell labels give analysis-island counts; darker shading indicates larger support.",
        ha="right",
        va="top",
        fontsize=4.4,
        color="#555555",
    )
    fig.subplots_adjust(left=0.20, right=0.98, top=0.90, bottom=0.22, wspace=0.12)
    return base._save(
        fig,
        output_dir / "extended_data",
        "extended_data_figure1_atomic_support",
    )


def _extended4_clean(source_dir: Path, output_dir: Path) -> list[Path]:
    bridge = pd.read_csv(source_dir / "extended_data_table5_functional_bridge.csv")
    traits = [
        "self_compatibility",
        "autonomous_selfing",
        "generalized_form",
        "actinomorphic_symmetry",
    ]
    analyses = [
        "primary",
        "supplemental_only",
        "no_zero_constant",
        "within_publication",
        "within_publication_site",
    ]
    fig, ax = base.plt.subplots(figsize=(FIGURE_WIDTH_IN, 3.65))
    ax.set_xlim(-0.5, len(analyses) - 0.5)
    ax.set_ylim(-0.5, len(traits) - 0.5)
    ax.set_xticks(
        range(len(analyses)),
        [
            "Primary",
            "Supplemental-only",
            "No-zero-constant",
            "Within publication",
            "Within publication × site",
        ],
        rotation=28,
        ha="right",
    )
    ax.set_yticks(range(len(traits)), [TRAIT_LABELS[t] for t in traits])
    reason_labels = {
        "insufficient_rows": "insufficient rows",
        "design_not_full_rank": "design not full rank",
    }
    for yi, trait in enumerate(traits):
        for xi, analysis in enumerate(analyses):
            part = bridge.loc[
                (bridge["trait"].eq(trait)) & (bridge["analysis"].eq(analysis))
            ]
            if part.empty or str(part.iloc[0]["evaluable"]).lower() != "true":
                raw_reason = (
                    "not_evaluable"
                    if part.empty
                    else str(part.iloc[0].get("reason", "not_evaluable"))
                )
                reason = reason_labels.get(raw_reason, raw_reason.replace("_", " "))
                ax.text(
                    xi,
                    yi,
                    "NE",
                    ha="center",
                    va="center",
                    fontsize=5.4,
                    color="#777777",
                )
                ax.text(
                    xi,
                    yi - 0.20,
                    reason,
                    ha="center",
                    va="center",
                    fontsize=3.5,
                    color="#888888",
                )
                continue
            estimate = float(part.iloc[0]["estimate"])
            se = float(part.iloc[0]["se"])
            p = float(part.iloc[0]["two_sided_p"])
            ax.scatter(
                xi,
                yi,
                s=30 + 65 * min(abs(estimate), 0.6) / 0.6,
                color="#009E73",
                edgecolor="white",
                linewidth=0.45,
            )
            ax.text(
                xi,
                yi - 0.20,
                f"β={estimate:.2f}\nP={p:.2g}",
                ha="center",
                va="center",
                fontsize=4.0,
            )
            if np.isfinite(se):
                ax.text(
                    xi,
                    yi + 0.22,
                    f"SE={se:.2f}",
                    ha="center",
                    va="center",
                    fontsize=3.8,
                    color="#555555",
                )
    ax.grid(color="#ededed", lw=0.35)
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax.set_title(
        "Functional compatibility: sensitivity and within-group checks",
        loc="left",
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.20, right=0.98, top=0.92, bottom=0.25)
    return base._save(
        fig,
        output_dir / "extended_data",
        "extended_data_figure4_functional_bridge_sensitivity",
    )


def render_publication_bundle(source_dir: Path, output_dir: Path) -> list[Path]:
    base._style()
    outputs: list[Path] = []
    outputs.extend(base._figure1(source_dir, output_dir))
    outputs.extend(_figure2_clean(source_dir, output_dir))
    outputs.extend(base._figure3(source_dir, output_dir))
    outputs.extend(_extended1_clean(source_dir, output_dir))
    outputs.extend(base._extended2(source_dir, output_dir))
    outputs.extend(base._extended3(source_dir, output_dir))
    outputs.extend(_extended4_clean(source_dir, output_dir))
    return outputs


@APP.command()
def cli(
    source_dir: Path = typer.Option(..., exists=True, file_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    render_publication_bundle(source_dir, output_dir)


if __name__ == "__main__":
    APP()
