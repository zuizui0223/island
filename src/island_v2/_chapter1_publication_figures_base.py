from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import typer  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

APP = typer.Typer(add_completion=False)

MM_PER_INCH = 25.4
FIGURE_WIDTH_IN = 180.0 / MM_PER_INCH

CONTEXT_ORDER = [
    "northern_midlatitude",
    "northern_high_latitude",
    "tropical",
    "southern_extratropical",
]
CONTEXT_LABELS = {
    "northern_midlatitude": "Northern mid-latitude",
    "northern_high_latitude": "Northern high-latitude",
    "tropical": "Tropical",
    "southern_extratropical": "Southern extratropical",
}
CONTEXT_COLORS = {
    "northern_midlatitude": "#0072B2",
    "northern_high_latitude": "#56B4E9",
    "tropical": "#E69F00",
    "southern_extratropical": "#CC79A7",
}
CONTEXT_MARKERS = {
    "northern_midlatitude": "o",
    "northern_high_latitude": "s",
    "tropical": "^",
    "southern_extratropical": "D",
}
TRAIT_ORDER = [
    "generalized_form",
    "actinomorphic_symmetry",
    "shallow_open_tube",
    "self_compatibility",
    "selfing_mating_system",
    "autonomous_selfing",
]
TRAIT_LABELS = {
    "generalized_form": "Generalized floral form",
    "actinomorphic_symmetry": "Actinomorphic symmetry",
    "shallow_open_tube": "Shallow / open tube",
    "self_compatibility": "Self-compatibility",
    "selfing_mating_system": "Selfing mating system",
    "autonomous_selfing": "Autonomous selfing",
}
FAMILY_LABELS = {
    "reproductive_assurance": "Reproductive assurance",
    "floral_accessibility": "Floral accessibility / generalization",
    "floral_architecture": "Floral accessibility / generalization",
}
SCOPE_LABELS = {"all_analysis": "Primary", "direct_only": "Direct-only"}

PUBLICATION_LABELS = (
    "Geographic isolation",
    "Experimental pollen limitation",
    "Recurrent plant response",
    "Reproductive assurance",
    "Floral accessibility / generalization",
    "Trait-informed island inputs",
    "GloPL experimental sites",
    "Primary evidence",
    "Direct-only sensitivity",
    "Current pollen limitation",
    "Post-hoc functional compatibility",
)


def _style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 6.2,
            "axes.labelsize": 6.2,
            "axes.titlesize": 7.0,
            "xtick.labelsize": 5.7,
            "ytick.labelsize": 5.7,
            "legend.fontsize": 5.4,
            "axes.linewidth": 0.45,
            "xtick.major.width": 0.45,
            "ytick.major.width": 0.45,
            "lines.linewidth": 0.65,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def _panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.04,
        label,
        transform=ax.transAxes,
        fontsize=8.3,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def _save(fig: plt.Figure, directory: Path, stem: str) -> list[Path]:
    directory.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    for suffix, kwargs in (
        ("pdf", {}),
        ("svg", {}),
        ("png", {"dpi": 400}),
    ):
        path = directory / f"{stem}.{suffix}"
        fig.savefig(path, bbox_inches="tight", **kwargs)
        outputs.append(path)
    plt.close(fig)
    return outputs


def _world_boundary(ax: plt.Axes, source_dir: Path) -> None:
    geojson = source_dir / "world_land.geojson"
    if not geojson.is_file():
        return
    import geopandas as gpd

    world = gpd.read_file(geojson)
    world.boundary.plot(ax=ax, color="#9a9a9a", linewidth=0.35, zorder=0)


def _figure1(source_dir: Path, output_dir: Path) -> list[Path]:
    islands = pd.read_csv(source_dir / "figure1_islands.csv")
    sites = pd.read_csv(source_dir / "figure1_glopl_sites.csv")
    tested = pd.read_csv(source_dir / "extended_data_table6_glopl_tested_islands.csv")

    fig = plt.figure(figsize=(FIGURE_WIDTH_IN, 3.95))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.45, 1.0], wspace=0.08)
    ax = fig.add_subplot(gs[0, 0])
    design = fig.add_subplot(gs[0, 1])

    _world_boundary(ax, source_dir)
    ax.scatter(
        islands["island_longitude"],
        islands["island_latitude"],
        s=1.4,
        c="#bdbdbd",
        alpha=0.20,
        linewidths=0,
        rasterized=True,
        zorder=1,
    )
    active = islands.loc[islands["trait_informed"].astype(bool)]
    for context in CONTEXT_ORDER:
        part = active.loc[active["analysis_regime"].eq(context)]
        ax.scatter(
            part["island_longitude"],
            part["island_latitude"],
            s=3.0,
            color=CONTEXT_COLORS[context],
            marker=CONTEXT_MARKERS[context],
            alpha=0.72,
            linewidths=0,
            rasterized=True,
            zorder=2,
        )
    ax.scatter(
        sites["Longitude"],
        sites["Latitude"],
        s=4,
        c="#3d3d3d",
        marker="x",
        alpha=0.28,
        linewidths=0.45,
        rasterized=True,
        zorder=3,
    )
    ax.scatter(
        tested["island_longitude"],
        tested["island_latitude"],
        s=28,
        facecolors="white",
        edgecolors="#111111",
        marker="*",
        linewidths=0.65,
        zorder=4,
    )
    ax.set_xlim(-180, 180)
    ax.set_ylim(-58, 84)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)
    _panel_label(ax, "a")
    ax.set_title("Global study coverage", loc="left", fontweight="bold", pad=3)

    handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=3.3, color="#bdbdbd", label="Frozen island universe (8,265)"),
    ]
    for context in CONTEXT_ORDER:
        n = int(active.loc[active["analysis_regime"].eq(context), "island_id"].nunique())
        handles.append(
            Line2D(
                [0],
                [0],
                marker=CONTEXT_MARKERS[context],
                linestyle="",
                markersize=3.4,
                color=CONTEXT_COLORS[context],
                label=f"{CONTEXT_LABELS[context]} ({n:,})",
            )
        )
    handles.extend(
        [
            Line2D([0], [0], marker="x", linestyle="", markersize=3.7, color="#3d3d3d", label="GloPL sites (1,248)"),
            Line2D([0], [0], marker="*", linestyle="", markersize=5.5, markerfacecolor="white", markeredgecolor="#111111", label="GloPL-tested islands (37)"),
        ]
    )
    ax.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(0.0, -0.02),
        ncol=2,
        frameon=False,
        handletextpad=0.35,
        columnspacing=0.9,
    )
    ax.text(
        0.0,
        -0.105,
        "4,453 islands contain trait information; model-specific complete-case n is reported with effect estimates.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        color="#4d4d4d",
        fontsize=5.3,
    )

    design.axis("off")
    _panel_label(design, "b")
    design.set_title("Evidence design", loc="left", fontweight="bold", pad=3)

    def box(y: float, text: str, edge: str = "#4d4d4d") -> None:
        design.text(
            0.5,
            y,
            text,
            ha="center",
            va="center",
            fontsize=6.2,
            bbox={"boxstyle": "round,pad=0.38", "fc": "white", "ec": edge, "lw": 0.65},
            transform=design.transAxes,
        )

    box(0.88, "Geographic isolation")
    box(0.65, "Recurrent plant response", "#0072B2")
    box(0.42, "Experimental pollen limitation", "#E69F00")
    box(0.19, "Reproductive assurance   |   Floral accessibility", "#666666")
    design.annotate("", xy=(0.5, 0.71), xytext=(0.5, 0.84), xycoords="axes fraction", arrowprops={"arrowstyle": "-|>", "lw": 0.65, "color": "#333333"})
    design.annotate("", xy=(0.5, 0.48), xytext=(0.5, 0.84), xycoords="axes fraction", arrowprops={"arrowstyle": "-|>", "lw": 0.65, "color": "#333333"})
    design.annotate("", xy=(0.5, 0.25), xytext=(0.5, 0.60), xycoords="axes fraction", arrowprops={"arrowstyle": "-|>", "lw": 0.65, "color": "#333333"})
    design.annotate("", xy=(0.77, 0.38), xytext=(0.77, 0.23), xycoords="axes fraction", arrowprops={"arrowstyle": "-|>", "lw": 0.6, "linestyle": "--", "color": "#666666"})
    design.text(0.77, 0.305, "functional\ncompatibility", transform=design.transAxes, ha="center", va="center", fontsize=5.0, color="#555555")
    design.text(0.5, 0.06, "Historical trait mediation is not identified.", transform=design.transAxes, ha="center", va="center", fontsize=5.2, color="#555555")

    fig.subplots_adjust(left=0.035, right=0.99, top=0.95, bottom=0.12)
    return _save(fig, output_dir / "main", "figure1_global_coverage")


def _forest_atomic(
    ax: plt.Axes,
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
    text = "   ".join(
        f"{CONTEXT_LABELS[c]}: n={int(support.loc[c, 'n_unique_islands']):,}, P={support.loc[c, 'p_value']:.2g}"
        for c in CONTEXT_ORDER
    )
    ax.text(0.0, -0.20, text, transform=ax.transAxes, ha="left", va="top", fontsize=4.8)


def _figure2(source_dir: Path, output_dir: Path) -> list[Path]:
    atomic = pd.read_csv(source_dir / "figure2_atomic_coefficients.csv")
    omnibus = pd.read_csv(source_dir / "figure2_omnibus_support.csv")
    summaries = pd.read_csv(source_dir / "figure2_descriptive_summaries.csv")
    max_abs = float(np.nanmax(np.abs(atomic[["ci_low", "ci_high"]].to_numpy())))
    xlim = (-max_abs * 1.08, max_abs * 1.08)

    fig, axes = plt.subplots(2, 2, figsize=(FIGURE_WIDTH_IN, 5.55))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()
    _forest_atomic(ax_a, atomic, omnibus, "all_analysis", xlim)
    _forest_atomic(ax_b, atomic, omnibus, "direct_only", xlim)
    ax_a.set_title("Primary evidence", loc="left", fontweight="bold")
    ax_b.set_title("Direct-only sensitivity", loc="left", fontweight="bold")
    _panel_label(ax_a, "a")
    _panel_label(ax_b, "b")
    ax_b.set_yticklabels([])

    handles = [
        Line2D([0], [0], marker=CONTEXT_MARKERS[c], linestyle="", color=CONTEXT_COLORS[c], label=CONTEXT_LABELS[c], markersize=3.5)
        for c in CONTEXT_ORDER
    ]
    ax_a.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, 1.15), ncol=2, frameon=False, columnspacing=0.8, handletextpad=0.3)

    classic = summaries.loc[summaries["summary"].eq("classic_orientation")].copy()
    y = np.arange(len(CONTEXT_ORDER))[::-1]
    for scope, marker, fill in (("all_analysis", "o", True), ("direct_only", "s", False)):
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
    _panel_label(ax_c, "c")

    family = summaries.loc[summaries["summary"].isin(["reproductive_assurance", "floral_accessibility"])].copy()
    positions = {c: i for i, c in enumerate(CONTEXT_ORDER[::-1])}
    for summary_name, color, offset in (
        ("reproductive_assurance", "#0072B2", -0.12),
        ("floral_accessibility", "#E69F00", 0.12),
    ):
        for scope, marker, open_marker in (("all_analysis", "o", False), ("direct_only", "s", True)):
            part = family.loc[(family["summary"].eq(summary_name)) & (family["evidence_scope"].eq(scope))]
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
        Line2D([0], [0], marker="o", linestyle="", color="#0072B2", label="Reproductive assurance", markersize=3.5),
        Line2D([0], [0], marker="o", linestyle="", color="#E69F00", label="Floral accessibility / generalization", markersize=3.5),
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor="#555555", markeredgecolor="#555555", color="#555555", label="Primary", markersize=3.5),
        Line2D([0], [0], marker="s", linestyle="", markerfacecolor="white", markeredgecolor="#555555", color="#555555", label="Direct-only", markersize=3.5),
    ]
    ax_d.legend(handles=family_handles, frameon=False, loc="lower right", fontsize=4.9)
    _panel_label(ax_d, "d")

    fig.subplots_adjust(left=0.17, right=0.99, top=0.94, bottom=0.10, wspace=0.23, hspace=0.48)
    return _save(fig, output_dir / "main", "figure2_recurrent_plant_response")


def _effect_forest(ax: plt.Axes, data: pd.DataFrame, labels: list[str], title: str) -> None:
    data = data.copy().reset_index(drop=True)
    y = np.arange(len(data))[::-1]
    role_style = {
        "primary": ("o", "#0072B2"),
        "sensitivity": ("s", "#56B4E9"),
        "posthoc_shape_diagnostic": ("D", "#E69F00"),
        "posthoc_functional_triangulation": ("o", "#009E73"),
    }
    for i, row in data.iterrows():
        marker, color = role_style.get(str(row.get("inferential_role")), ("o", "#555555"))
        est = float(row["estimate"])
        low = float(row["ci_low"])
        high = float(row["ci_high"])
        ax.errorbar(
            est,
            y[i],
            xerr=np.array([[est - low], [high - est]]),
            fmt=marker,
            markersize=3.5,
            color=color,
            ecolor=color,
            elinewidth=0.65,
            capsize=1.6,
            markeredgewidth=0,
        )
    ax.axvline(0, color="#777777", lw=0.5, ls="--")
    ax.set_yticks(y, labels)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.grid(axis="x", color="#e6e6e6", lw=0.35)
    ax.spines[["top", "right", "left"]].set_visible(False)


def _figure3(source_dir: Path, output_dir: Path) -> list[Path]:
    data = pd.read_csv(source_dir / "figure3_source_data.csv")
    fig, axes = plt.subplots(2, 2, figsize=(FIGURE_WIDTH_IN, 5.45))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    global_row = data.loc[data["row_id"].eq("global_distance")].copy()
    _effect_forest(ax_a, global_row, ["Global distance effect"], "Experimental pollen limitation")
    row = global_row.iloc[0]
    ax_a.set_xlabel("Standardized effect on pollen limitation")
    ax_a.text(
        0.03,
        0.08,
        f"2,969 effect rows | {int(row['n_sites']):,} sites | {int(row['n_publications']):,} publications\n"
        f"estimate = {row['estimate']:.3f} ± {row['se']:.3f}; two-sided P = {row['two_sided_p']:.3g}",
        transform=ax_a.transAxes,
        ha="left",
        va="bottom",
        fontsize=5.3,
    )
    _panel_label(ax_a, "a")

    ids = [
        "global_distance",
        "global_supplemental_only",
        "global_no_zero_constant",
        "mainland_to_offshore_step",
        "within_offshore_gradient",
        "offshore_only_gradient",
    ]
    sens = data.set_index("row_id").loc[ids].reset_index()
    labels = sens["label"].tolist()
    _effect_forest(ax_b, sens, labels, "Sensitivity and shape diagnostics")
    ax_b.set_xlabel("Standardized effect on pollen limitation")
    ax_b.text(0.02, -0.20, "Circles/squares: frozen global model and sensitivities; diamonds: post-hoc shape diagnostics.", transform=ax_b.transAxes, ha="left", va="top", fontsize=4.8)
    _panel_label(ax_b, "b")

    bridge = data.loc[
        data["section"].eq("functional_bridge")
        & data["analysis"].eq("primary")
        & data["evaluable"].astype(str).str.lower().eq("true")
    ].copy()
    trait_order = ["autonomous_selfing", "self_compatibility", "actinomorphic_symmetry", "generalized_form"]
    bridge = bridge.set_index("trait").loc[trait_order].reset_index()
    trait_labels = [TRAIT_LABELS[t] for t in trait_order]
    _effect_forest(ax_c, bridge, trait_labels, "Trait state and current pollen limitation")
    ax_c.set_xlabel("Trait-state effect on pollen limitation")
    support_text = "\n".join(
        f"{TRAIT_LABELS[r.trait]}: {int(r.n_sites):,} sites, {int(r.n_species):,} species"
        for r in bridge.itertuples()
    )
    ax_c.text(0.02, -0.22, support_text, transform=ax_c.transAxes, ha="left", va="top", fontsize=4.5)
    _panel_label(ax_c, "c")

    matrix = data.loc[data["section"].eq("functional_bridge")].copy()
    analyses = ["primary", "supplemental_only", "no_zero_constant", "within_publication", "within_publication_site"]
    traits = ["self_compatibility", "autonomous_selfing", "generalized_form", "actinomorphic_symmetry"]
    ax_d.set_xlim(-0.5, len(analyses) - 0.5)
    ax_d.set_ylim(-0.5, len(traits) - 0.5)
    ax_d.set_xticks(range(len(analyses)), ["Primary", "Supplemental", "No-zero", "Within pub.", "Within pub.×site"], rotation=35, ha="right")
    ax_d.set_yticks(range(len(traits)), [TRAIT_LABELS[t] for t in traits])
    for yi, trait in enumerate(traits):
        for xi, analysis in enumerate(analyses):
            part = matrix.loc[(matrix["trait"].eq(trait)) & (matrix["analysis"].eq(analysis))]
            if part.empty or not str(part.iloc[0]["evaluable"]).lower() == "true":
                ax_d.text(xi, yi, "NE", ha="center", va="center", fontsize=4.8, color="#777777")
                continue
            estimate = float(part.iloc[0]["estimate"])
            ax_d.scatter(xi, yi, s=24 + 55 * min(abs(estimate), 0.6) / 0.6, facecolor="#009E73", edgecolor="white", linewidth=0.45)
            ax_d.text(xi, yi - 0.24, f"{estimate:.2f}", ha="center", va="center", fontsize=4.2)
    ax_d.set_title("Robustness and estimability", loc="left", fontweight="bold")
    ax_d.grid(color="#ededed", lw=0.35)
    ax_d.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax_d.text(0.0, -0.26, "Frozen distance × trait moderation: neither response family was supported; these tests are not reclassified.", transform=ax_d.transAxes, ha="left", va="top", fontsize=4.7)
    _panel_label(ax_d, "d")

    fig.subplots_adjust(left=0.19, right=0.99, top=0.95, bottom=0.12, wspace=0.36, hspace=0.52)
    return _save(fig, output_dir / "main", "figure3_pollen_limitation_bridge")


def _extended1(source_dir: Path, output_dir: Path) -> list[Path]:
    atomic = pd.read_csv(source_dir / "figure2_atomic_coefficients.csv")
    fig, axes = plt.subplots(1, 2, figsize=(FIGURE_WIDTH_IN, 3.45), sharey=True)
    for ax, scope, title in zip(axes, ("all_analysis", "direct_only"), ("Primary", "Direct-only"), strict=True):
        data = atomic.loc[atomic["evidence_scope"].eq(scope)]
        pivot = data.pivot(index="outcome", columns="context", values="n_islands").reindex(index=TRAIT_ORDER, columns=CONTEXT_ORDER)
        image = ax.imshow(pivot.to_numpy(), cmap="Greys", aspect="auto")
        for yi in range(len(TRAIT_ORDER)):
            for xi in range(len(CONTEXT_ORDER)):
                ax.text(xi, yi, f"{int(pivot.iloc[yi, xi]):,}", ha="center", va="center", fontsize=4.8)
        ax.set_xticks(range(4), [CONTEXT_LABELS[c] for c in CONTEXT_ORDER], rotation=35, ha="right")
        ax.set_yticks(range(6), [TRAIT_LABELS[t] for t in TRAIT_ORDER])
        ax.set_title(title, loc="left", fontweight="bold")
        fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02, label="Analysis islands")
    _panel_label(axes[0], "a")
    _panel_label(axes[1], "b")
    fig.subplots_adjust(left=0.20, right=0.98, top=0.93, bottom=0.22, wspace=0.12)
    return _save(fig, output_dir / "extended_data", "extended_data_figure1_atomic_support")


def _extended2(source_dir: Path, output_dir: Path) -> list[Path]:
    islands = pd.read_csv(source_dir / "figure1_islands.csv")
    sites = pd.read_csv(source_dir / "figure1_glopl_sites.csv")
    tested = pd.read_csv(source_dir / "extended_data_table6_glopl_tested_islands.csv")
    fig, axes = plt.subplots(1, 2, figsize=(FIGURE_WIDTH_IN, 3.25))
    labels = ["Frozen islands", "Trait-informed inputs", "GloPL sites", "Island GloPL sites", "GloPL-tested islands"]
    counts = [8265, int(islands["trait_informed"].sum()), sites["site_key"].nunique(), int(sites["is_frozen_island_site"].sum()), tested["island_id"].nunique()]
    axes[0].barh(range(len(labels))[::-1], counts, color="#9a9a9a")
    axes[0].set_yticks(range(len(labels))[::-1], labels)
    axes[0].set_xlabel("Count")
    axes[0].spines[["top", "right", "left"]].set_visible(False)
    axes[0].grid(axis="x", color="#ededed", lw=0.35)
    for y, value in zip(range(len(labels))[::-1], counts, strict=True):
        axes[0].text(value, y, f"  {value:,}", va="center", fontsize=5.1)
    _panel_label(axes[0], "a")
    axes[0].set_title("Sampling and overlap", loc="left", fontweight="bold")

    regime_counts = tested.groupby("analysis_regime")["island_id"].nunique().reindex(CONTEXT_ORDER, fill_value=0)
    axes[1].barh(range(4)[::-1], regime_counts.values, color=[CONTEXT_COLORS[c] for c in CONTEXT_ORDER])
    axes[1].set_yticks(range(4)[::-1], [CONTEXT_LABELS[c] for c in CONTEXT_ORDER])
    axes[1].set_xlabel("GloPL-tested frozen islands")
    axes[1].spines[["top", "right", "left"]].set_visible(False)
    axes[1].grid(axis="x", color="#ededed", lw=0.35)
    for y, value in zip(range(4)[::-1], regime_counts.values, strict=True):
        axes[1].text(value, y, f"  {int(value)}", va="center", fontsize=5.1)
    _panel_label(axes[1], "b")
    axes[1].set_title("Exact-island GloPL overlap", loc="left", fontweight="bold")
    fig.subplots_adjust(left=0.18, right=0.98, top=0.92, bottom=0.16, wspace=0.38)
    return _save(fig, output_dir / "extended_data", "extended_data_figure2_sampling_overlap")


def _extended3(source_dir: Path, output_dir: Path) -> list[Path]:
    data = pd.read_csv(source_dir / "extended_data_table4_glopl_estimates.csv")
    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH_IN, 3.4))
    _effect_forest(ax, data, data["label"].tolist(), "Experimental pollen-limitation model checks")
    ax.set_xlabel("Standardized effect on pollen limitation")
    for i, row in data.reset_index(drop=True).iterrows():
        y = len(data) - 1 - i
        ax.text(ax.get_xlim()[1], y, f"P={row['two_sided_p']:.3g}; sites={int(row['n_sites']):,}", ha="right", va="center", fontsize=4.8)
    fig.subplots_adjust(left=0.25, right=0.97, top=0.91, bottom=0.18)
    return _save(fig, output_dir / "extended_data", "extended_data_figure3_glopl_sensitivity")


def _extended4(source_dir: Path, output_dir: Path) -> list[Path]:
    bridge = pd.read_csv(source_dir / "extended_data_table5_functional_bridge.csv")
    traits = ["self_compatibility", "autonomous_selfing", "generalized_form", "actinomorphic_symmetry"]
    analyses = ["primary", "supplemental_only", "no_zero_constant", "within_publication", "within_publication_site"]
    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH_IN, 3.65))
    ax.set_xlim(-0.5, len(analyses) - 0.5)
    ax.set_ylim(-0.5, len(traits) - 0.5)
    ax.set_xticks(range(len(analyses)), ["Primary", "Supplemental-only", "No-zero-constant", "Within publication", "Within publication × site"], rotation=28, ha="right")
    ax.set_yticks(range(len(traits)), [TRAIT_LABELS[t] for t in traits])
    for yi, trait in enumerate(traits):
        for xi, analysis in enumerate(analyses):
            part = bridge.loc[(bridge["trait"].eq(trait)) & (bridge["analysis"].eq(analysis))]
            if part.empty or not str(part.iloc[0]["evaluable"]).lower() == "true":
                reason = "not evaluable" if part.empty else str(part.iloc[0].get("reason", "not evaluable"))
                ax.text(xi, yi, "NE", ha="center", va="center", fontsize=5.4, color="#777777")
                ax.text(xi, yi - 0.20, reason.replace("_", " ")[:18], ha="center", va="center", fontsize=3.7, color="#888888")
                continue
            estimate = float(part.iloc[0]["estimate"])
            se = float(part.iloc[0]["se"])
            p = float(part.iloc[0]["two_sided_p"])
            ax.scatter(xi, yi, s=30 + 65 * min(abs(estimate), 0.6) / 0.6, color="#009E73", edgecolor="white", linewidth=0.45)
            ax.text(xi, yi - 0.20, f"β={estimate:.2f}\nP={p:.2g}", ha="center", va="center", fontsize=4.0)
            if np.isfinite(se):
                ax.text(xi, yi + 0.22, f"SE={se:.2f}", ha="center", va="center", fontsize=3.8, color="#555555")
    ax.grid(color="#ededed", lw=0.35)
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax.set_title("Functional compatibility: sensitivity and within-group checks", loc="left", fontweight="bold")
    fig.subplots_adjust(left=0.20, right=0.98, top=0.92, bottom=0.25)
    return _save(fig, output_dir / "extended_data", "extended_data_figure4_functional_bridge_sensitivity")


def render_publication_bundle(source_dir: Path, output_dir: Path) -> list[Path]:
    _style()
    outputs: list[Path] = []
    outputs.extend(_figure1(source_dir, output_dir))
    outputs.extend(_figure2(source_dir, output_dir))
    outputs.extend(_figure3(source_dir, output_dir))
    outputs.extend(_extended1(source_dir, output_dir))
    outputs.extend(_extended2(source_dir, output_dir))
    outputs.extend(_extended3(source_dir, output_dir))
    outputs.extend(_extended4(source_dir, output_dir))
    return outputs


@APP.command()
def cli(
    source_dir: Path = typer.Option(..., exists=True, file_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    render_publication_bundle(source_dir, output_dir)


if __name__ == "__main__":
    APP()
