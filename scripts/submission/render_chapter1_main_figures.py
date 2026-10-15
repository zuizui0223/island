from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from xml.etree import ElementTree as ET

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch
from scipy.stats import t as student_t

REGIONS = [
    "northern_midlatitude",
    "northern_high_latitude",
    "tropical",
    "southern_extratropical",
]
REGION_LABELS = {
    "northern_midlatitude": "Northern mid-latitude",
    "northern_high_latitude": "Northern high latitude",
    "tropical": "Tropical",
    "southern_extratropical": "Southern extratropical",
}
H1_OUTCOMES = [
    "self_compatibility",
    "selfing_mating_system",
    "autonomous_selfing",
    "plain_colour",
    "generalized_form",
    "actinomorphic_symmetry",
    "shallow_open_tube",
]
H1_LABELS = {
    "self_compatibility": "Self-compatibility",
    "selfing_mating_system": "Selfing mating system",
    "autonomous_selfing": "Autonomous selfing",
    "plain_colour": "Plain colour",
    "generalized_form": "Generalized form",
    "actinomorphic_symmetry": "Actinomorphic symmetry",
    "shallow_open_tube": "Shallow/open tube",
}


def _style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "svg.hashsalt": "chapter1-main-figures-20260927",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def _save(fig: plt.Figure, out: Path, stem: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(
        out / f"{stem}.pdf",
        bbox_inches="tight",
        metadata={
            "Title": stem,
            "Author": "Chapter 1 submission pipeline",
            "Creator": "matplotlib",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(fig)


def _coord_columns(frame: pd.DataFrame) -> tuple[str, str]:
    lat_candidates = [
        "island_latitude",
        "latitude",
        "lat",
        "centroid_latitude",
        "centroid_lat",
    ]
    lon_candidates = [
        "island_longitude",
        "longitude",
        "lon",
        "long",
        "centroid_longitude",
        "centroid_lon",
    ]
    lat = next((c for c in lat_candidates if c in frame.columns), None)
    lon = next((c for c in lon_candidates if c in frame.columns), None)
    if lat is None or lon is None:
        raise ValueError(
            "Could not locate island latitude/longitude columns in corrected covariates; "
            f"available columns={list(frame.columns)}"
        )
    return lat, lon


def figure1(root: Path, out: Path) -> None:
    cov = pd.read_csv(root / "results/geography_20260924/corrected_geography_covariates.csv")
    sites = pd.read_csv(root / "results/geography_20260924/glopl_corrected_site_distances.csv")
    summary = pd.read_csv(
        root / "submission/chapter1_current/supplement/Table_S1_data_summary.csv"
    )
    former = pd.read_csv(
        root / "results/geography_20260924/formerly_zero_islands_recalculated.csv"
    )
    lat_col, lon_col = _coord_columns(cov)
    if "analysis_regime" not in cov.columns:
        raise ValueError("corrected covariates missing analysis_regime")

    fig = plt.figure(figsize=(12.2, 8.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 0.85])

    ax1 = fig.add_subplot(gs[0, 0])
    for region in REGIONS:
        part = cov.loc[cov["analysis_regime"].astype(str).eq(region)]
        ax1.scatter(
            part[lon_col],
            part[lat_col],
            s=6,
            alpha=0.55,
            label=REGION_LABELS[region],
        )
    ax1.set_xlim(-180, 180)
    ax1.set_ylim(-90, 90)
    ax1.set_xlabel("Longitude")
    ax1.set_ylabel("Latitude")
    ax1.set_title("A | Corrected island analysis universe")
    ax1.legend(loc="lower left", frameon=False, ncol=2)
    ax1.grid(linewidth=0.25, alpha=0.35)

    ax2 = fig.add_subplot(gs[0, 1])
    continent = sites["on_seeded_continent"].astype(str).str.lower().isin({"true", "1"})
    ax2.scatter(
        sites.loc[continent, "longitude"],
        sites.loc[continent, "latitude"],
        s=7,
        alpha=0.35,
        label="Continental",
    )
    ax2.scatter(
        sites.loc[~continent, "longitude"],
        sites.loc[~continent, "latitude"],
        s=12,
        alpha=0.65,
        label="Non-continental",
    )
    ax2.set_xlim(-180, 180)
    ax2.set_ylim(-90, 90)
    ax2.set_xlabel("Longitude")
    ax2.set_ylabel("Latitude")
    ax2.set_title("B | GloPL pollen-supplementation sites")
    ax2.legend(frameon=False)
    ax2.grid(linewidth=0.25, alpha=0.35)

    ax3 = fig.add_subplot(gs[1, 0])
    metrics = {
        row["metric"]: float(row["value"])
        for _, row in summary.iterrows()
        if row["section"] in {"traits", "GloPL"}
    }
    names = ["Species", "Resolved trait cells", "GloPL rows", "GloPL sites"]
    values = [
        metrics["accepted_angiosperm_species"],
        metrics["resolved_species_axis_cells"],
        metrics["effect_rows"],
        metrics["sites"],
    ]
    y = np.arange(len(names))
    ax3.barh(y, values)
    ax3.set_yticks(y, names)
    ax3.invert_yaxis()
    ax3.set_xscale("log")
    ax3.set_xlabel("Count (log scale)")
    ax3.set_title("C | Data coverage")
    for yi, value in zip(y, values, strict=True):
        ax3.text(value * 1.08, yi, f"{int(value):,}", va="center", fontsize=8)

    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis("off")
    median = float(former["spherical_coast_distance_km"].median())
    lines = [
        "D | Geography audit",
        "",
        "8,264 corrected island units",
        "4,379 broad H1 islands",
        "1,113 spurious island zero distances repaired",
        f"median repaired distance = {median:.3f} km",
        "1 continental split component excluded",
        "996 true continental GloPL zeros retained",
        "",
        "Exposure: source-matched GSHHG 2.3.7 coastlines",
        "Metric: minimum minor-great-circle arc separation",
    ]
    ax4.text(0.03, 0.96, "\n".join(lines), va="top", ha="left", fontsize=10)

    fig.suptitle("Figure 1 | Global scope and corrected geographic exposure", y=0.995)
    fig.tight_layout()
    _save(fig, out, "Figure1_global_scope_corrected")


def _box(ax: plt.Axes, x: float, y: float, w: float, h: float, title: str, body: str) -> None:
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        linewidth=1.3,
        fill=False,
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h * 0.68, title, ha="center", va="center", fontweight="bold")
    ax.text(x + w / 2, y + h * 0.34, body, ha="center", va="center", fontsize=8)


def figure2(root: Path, out: Path) -> None:
    del root
    fig, ax = plt.subplots(figsize=(11.5, 6.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    _box(ax, 0.04, 0.68, 0.25, 0.18, "Island geography", "GSHHG 2.3.7\ncorrected isolation")
    _box(ax, 0.375, 0.68, 0.25, 0.18, "Island floras", "GBIF incidence\n106,295 angiosperms")
    _box(ax, 0.71, 0.68, 0.25, 0.18, "Trait evidence", "colour · architecture\nreproductive assurance")
    _box(ax, 0.12, 0.36, 0.30, 0.20, "H1–H2 | Plant response", "multivariate recurrence\nconditional decomposition")
    _box(ax, 0.58, 0.36, 0.30, 0.20, "H3 | Ecological pressure", "independent GloPL\npollen-supplementation data")
    _box(ax, 0.35, 0.08, 0.30, 0.18, "H4 | Functional bridge", "exact-species H2 scores × GloPL\npost-hoc triangulation")

    arrows = [
        ((0.29, 0.77), (0.375, 0.77)),
        ((0.625, 0.77), (0.71, 0.77)),
        ((0.50, 0.68), (0.28, 0.56)),
        ((0.83, 0.68), (0.73, 0.56)),
        ((0.27, 0.36), (0.44, 0.26)),
        ((0.73, 0.36), (0.56, 0.26)),
    ]
    for (x0, y0), (x1, y1) in arrows:
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops={"arrowstyle": "->"})
    ax.text(
        0.5,
        0.94,
        "Figure 2 | Database construction and analytical workflow",
        ha="center",
        fontsize=13,
        fontweight="bold",
    )
    ax.text(
        0.5,
        0.015,
        "H1/H2 use island-trait data; H3 is independent; H4 joins exact species only.",
        ha="center",
        fontsize=8.5,
    )
    _save(fig, out, "Figure2_analytical_workflow")



def figure3(root: Path, out: Path) -> None:
    del root
    fig, ax = plt.subplots(figsize=(9.5, 6.8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    _box(
        ax,
        0.36,
        0.76,
        0.28,
        0.12,
        "Geographic isolation",
        "source separation / connectivity gradient",
    )
    _box(
        ax,
        0.06,
        0.38,
        0.34,
        0.18,
        "Plant response (H1–H2)",
        "reproductive assurance\n+ floral accessibility",
    )
    _box(
        ax,
        0.60,
        0.38,
        0.34,
        0.18,
        "Pollen limitation (H3)",
        "experimental pollen supplementation\nβ = +0.0919",
    )

    for start, end, label in (
        ((0.43, 0.76), (0.28, 0.56), "H1–H2"),
        ((0.57, 0.76), (0.72, 0.56), "H3"),
        ((0.40, 0.47), (0.60, 0.47), "H4: β < 0"),
    ):
        ax.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "->", "lw": 1.8})
        mx = (start[0] + end[0]) / 2
        my = (start[1] + end[1]) / 2
        ax.text(mx, my + 0.035, label, ha="center", va="center", fontweight="bold")

    ax.annotate(
        "",
        xy=(0.27, 0.38),
        xytext=(0.73, 0.38),
        arrowprops={"arrowstyle": "->", "lw": 1.4, "linestyle": "--"},
    )
    ax.text(
        0.50,
        0.30,
        "Historical causal edge not identified",
        ha="center",
        va="center",
        fontweight="bold",
    )
    ax.text(
        0.50,
        0.255,
        "past pollen limitation → selection / sorting / persistence → present traits",
        ha="center",
        va="center",
        fontsize=8.5,
    )
    ax.text(
        0.50,
        0.12,
        "Solid edges are fitted associations, not a mediation model.",
        ha="center",
        va="center",
        fontsize=9,
    )
    ax.text(
        0.50,
        0.95,
        "Figure 3 | Constraint–response triangle and inferential boundary",
        ha="center",
        va="center",
        fontsize=13,
        fontweight="bold",
    )
    _save(fig, out, "Figure3_constraint_response_triangle")



def _q_label(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    if value <= 0:
        return "<1e-300"
    if value < 1e-3:
        return f"{value:.1e}"
    return f"{value:.3f}"


def figure4(root: Path, out: Path) -> None:
    """Promote the frozen final traitwise H1 figure without re-deriving H1."""
    source = root / "results/h1_final_traitwise_t_20261004"
    out.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(
        source / "traitwise_H1.svg",
        out / "Figure4_H1_recurrent_multivariate_response.svg",
    )
    shutil.copyfile(
        source / "traitwise_H1.pdf",
        out / "Figure4_H1_recurrent_multivariate_response.pdf",
    )


def figure5(root: Path, out: Path) -> None:
    h2 = pd.read_csv(
        root / "results/geography_20260924/all/h2_decomposition_models.csv"
    ).copy()
    fit = h2["status"].astype(str).eq("fit")
    primary = fit & h2["response"].isin(
        ["generalized_accessible", "plain_colour"]
    )
    h2["distance_p_finite_cluster_t"] = np.nan
    for idx in h2.index[fit]:
        estimate = float(h2.loc[idx, "distance_estimate"])
        se = float(h2.loc[idx, "distance_se"])
        df = int(h2.loc[idx, "n_clusters"]) - 1
        h2.loc[idx, "distance_p_finite_cluster_t"] = float(
            2.0 * student_t.sf(abs(estimate / se), df=df)
        )
    p = h2.loc[primary, "distance_p_finite_cluster_t"].to_numpy(float)
    order = np.argsort(p)
    ranked = p[order]
    n = len(ranked)
    adj = np.minimum.accumulate(
        (ranked * n / np.arange(1, n + 1))[::-1]
    )[::-1]
    restored = np.empty(n, dtype=float)
    restored[order] = np.clip(adj, 0.0, 1.0)
    h2.loc[primary, "primary_H2b_q_finite_cluster_t"] = restored

    responses = ["selfing_core", "generalized_accessible", "plain_colour"]
    labels = [
        "Reproductive assurance",
        "Accessibility | assurance",
        "Plain colour | assurance",
    ]
    fig, axes = plt.subplots(1, 4, figsize=(13.0, 4.3), sharey=True)
    y = np.arange(len(responses))
    for ax, region in zip(axes, REGIONS, strict=True):
        part = h2.loc[h2["context"].eq(region)].set_index("response")
        est = np.array(
            [float(part.loc[r, "distance_estimate"]) for r in responses]
        )
        se = np.array(
            [float(part.loc[r, "distance_se"]) for r in responses]
        )
        crit = np.array(
            [
                float(
                    student_t.ppf(
                        0.975,
                        df=int(part.loc[r, "n_clusters"]) - 1,
                    )
                )
                for r in responses
            ]
        )
        ax.errorbar(est, y, xerr=crit * se, fmt="o", capsize=3)
        ax.axvline(0, linewidth=0.8)
        access_q = float(
            part.loc[
                "generalized_accessible",
                "primary_H2b_q_finite_cluster_t",
            ]
        )
        ax.set_title(
            f"{REGION_LABELS[region]}\naccessibility q={access_q:.3g}"
        )
        ax.set_xlabel("Isolation coefficient")
        ax.grid(axis="x", linewidth=0.3, alpha=0.4)
    axes[0].set_yticks(y, labels)
    axes[0].invert_yaxis()
    fig.suptitle(
        "Figure 5 | H2: reproductive assurance and additional floral response",
        y=1.02,
    )
    fig.tight_layout()
    _save(fig, out, "Figure5_H2_conditional_decomposition")

def figure6(root: Path, out: Path) -> None:
    h3 = pd.read_csv(
        root / "submission/chapter1_current/supplement/Table_S5_H3_pollen_limitation.csv"
    )
    h4 = pd.read_csv(
        root / "submission/chapter1_current/supplement/Table_S6a_H4_scores.csv"
    )
    atomic = pd.read_csv(root / "results/geography_20260924/h4_atomic_corrected.csv")

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.6))

    order = ["primary", "no_zero_constant", "supplemental_only"]
    labels = ["Primary", "No-zero constant", "Supplemental only"]
    h3p = h3.set_index("analysis").loc[order]
    y = np.arange(3)
    axes[0].errorbar(
        h3p["estimate"].astype(float),
        y,
        xerr=1.96 * h3p["se"].astype(float),
        fmt="o",
        capsize=3,
    )
    axes[0].axvline(0, linewidth=0.8)
    axes[0].set_yticks(y, labels)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Isolation coefficient")
    axes[0].set_title("A | H3 pollen limitation")

    h4p = h4.loc[h4["analysis"].eq("primary")].copy()
    h4_labels = ["Reproductive assurance", "Accessibility"]
    y2 = np.arange(2)
    axes[1].errorbar(
        h4p["estimate"].astype(float),
        y2,
        xerr=1.96 * h4p["se"].astype(float),
        fmt="o",
        capsize=3,
    )
    axes[1].axvline(0, linewidth=0.8)
    axes[1].set_yticks(y2, h4_labels)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Trait-score coefficient")
    axes[1].set_title("B | H4 exact H2 scores")

    ap = atomic.loc[
        atomic["analysis"].eq("primary") & atomic["evaluable"].astype(str).eq("True")
    ].copy()
    ap["estimate"] = pd.to_numeric(ap["trait_state_estimate"], errors="coerce")
    ap["se"] = pd.to_numeric(ap["trait_state_se"], errors="coerce")
    ap = ap.dropna(subset=["estimate", "se"])
    atomic_labels = [H1_LABELS.get(t, t.replace("_", " ").title()) for t in ap["trait"]]
    y3 = np.arange(len(ap))
    axes[2].errorbar(
        ap["estimate"],
        y3,
        xerr=1.96 * ap["se"],
        fmt="o",
        capsize=3,
    )
    axes[2].axvline(0, linewidth=0.8)
    axes[2].set_yticks(y3, atomic_labels)
    axes[2].invert_yaxis()
    axes[2].set_xlabel("Trait-state coefficient")
    axes[2].set_title("C | H4 atomic sensitivities")

    for ax in axes:
        ax.grid(axis="x", linewidth=0.3, alpha=0.4)
    fig.suptitle("Figure 6 | Isolation-associated pollen limitation and functional compatibility", y=1.02)
    fig.tight_layout()
    _save(fig, out, "Figure6_H3_H4_functional_bridge")


def build_manifest(out: Path) -> None:
    expected = {
        1: "Figure1_global_scope_corrected",
        2: "Figure2_analytical_workflow",
        3: "Figure3_constraint_response_triangle",
        4: "Figure4_H1_recurrent_multivariate_response",
        5: "Figure5_H2_conditional_decomposition",
        6: "Figure6_H3_H4_functional_bridge",
    }
    records: list[dict[str, object]] = []
    for number, stem in expected.items():
        svg = out / f"{stem}.svg"
        if not svg.is_file():
            raise FileNotFoundError(svg)
        ET.parse(svg)
        pdf = out / f"{stem}.pdf"
        if not pdf.is_file():
            raise FileNotFoundError(pdf)
        records.append(
            {
                "figure": number,
                "stem": stem,
                "svg": svg.name,
                "svg_bytes": svg.stat().st_size,
                "pdf": pdf.name,
                "pdf_bytes": pdf.stat().st_size,
            }
        )
    (out / "FIGURE_MANIFEST.json").write_text(
        json.dumps(
            {
                "contract": "chapter1_main_figures_corrected_v1",
                "source_surface": "corrected_geography_20260924_plus_final_traitwise_H1_20261004",
                "figures": records,
            },
            indent=2,
        ),
        encoding="utf-8",
    )



# ---- Q1-only five-figure submission display set (2026-10-15) ----
# Original six source figures remain replayable via the default CLI path.
# Do not reuse the source Figure numbers as Q1 Main figure numbers.

Q1_FIGURE_STEMS = {
    1: "Q1_Fig1_distinct_data_streams",
    2: "Q1_Fig2_pollen_limitation_gradient",
    3: "Q1_Fig3_seven_atomic_trait_responses",
    4: "Q1_Fig4_selfing_adjusted_floral_response",
    5: "Q1_Fig5_posthoc_functional_compatibility",
}


def q1_figure1(root: Path, out: Path) -> None:
    """Independent maps and causal boundary at two-column print width (178 mm)."""
    cov = pd.read_csv(root / "results/geography_20260924/corrected_geography_covariates.csv")
    glopl = pd.read_csv(root / "results/geography_20260924/glopl_corrected_site_distances.csv")
    lat, lon = _coord_columns(cov)
    fig = plt.figure(figsize=(7.0, 7.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.35, 0.85, 1.15], hspace=0.5, wspace=0.32)

    ax_island = fig.add_subplot(gs[0, 0])
    for region in REGIONS:
        part = cov.loc[cov["analysis_regime"].astype(str).eq(region)]
        ax_island.scatter(part[lon], part[lat], s=2.2, alpha=0.62,
                          label=REGION_LABELS[region])
    ax_island.set(title="A | Island trait observations", xlim=(-180, 180), ylim=(-90, 90),
                  xlabel="Longitude", ylabel="Latitude")
    ax_island.set_xticks([-180, -90, 0, 90, 180])
    ax_island.set_yticks([-60, 0, 60])
    ax_island.grid(alpha=0.2, linewidth=0.4)

    ax_glopl = fig.add_subplot(gs[0, 1])
    continent = glopl["on_seeded_continent"].astype(str).str.lower().isin({"true", "1"})
    for mask, name in [(continent, "Continental"), (~continent, "Offshore")]:
        ax_glopl.scatter(glopl.loc[mask, "longitude"], glopl.loc[mask, "latitude"],
                         s=3.0, alpha=0.55, label=name)
    ax_glopl.set(title="B | GloPL published sites", xlim=(-180, 180), ylim=(-90, 90),
                 xlabel="Longitude", ylabel="Latitude")
    ax_glopl.set_xticks([-180, -90, 0, 90, 180])
    ax_glopl.set_yticks([-60, 0, 60])
    ax_glopl.legend(loc="lower left", frameon=False, fontsize=7.5)
    ax_glopl.grid(alpha=0.2, linewidth=0.4)

    ax_counts = fig.add_subplot(gs[1, 0]); ax_counts.axis("off")
    ax_counts.text(0.01, 0.99, "C | Island assemblage data", va="top",
                   fontweight="bold", fontsize=10, transform=ax_counts.transAxes)
    ax_counts.text(0.01, 0.81,
                   "8,264 corrected island units\n"
                   "4,379 broad-H1 sampled islands\n"
                   "106,295 accepted species (universe)\n"
                   "222,688/318,885 resolved trait cells\n"
                   "69.83% raw-axis completeness",
                   fontsize=9, va="top", linespacing=1.32,
                   transform=ax_counts.transAxes)

    ax_counts2 = fig.add_subplot(gs[1, 1]); ax_counts2.axis("off")
    ax_counts2.text(0.01, 0.99, "D | GloPL experiments", va="top",
                    fontweight="bold", fontsize=10, transform=ax_counts2.transAxes)
    ax_counts2.text(0.01, 0.81,
                    "2,969 experimental rows\n"
                    "1,248 sites / 919 publications\n"
                    "996 true continental zero sites\n"
                    "Exact H4 matches: 455 / 143 species\n"
                    "(two distinct analysis sets)",
                    fontsize=9, va="top", linespacing=1.32,
                    transform=ax_counts2.transAxes)

    ax_flow = fig.add_subplot(gs[2, :]); ax_flow.axis("off")
    ax_flow.set(xlim=(0, 1), ylim=(0, 1))
    _box(ax_flow, 0.34, 0.72, 0.32, 0.17, "Geographic isolation", "separate exposure data")
    _box(ax_flow, 0.04, 0.28, 0.40, 0.20, "Reproductive / floral states",
         "island-flora H1 + H2")
    _box(ax_flow, 0.56, 0.28, 0.40, 0.20, "Current pollen limitation",
         "GloPL H3 + matched H4")
    for start, stop in [((0.42, 0.72), (0.24, 0.48)),
                        ((0.58, 0.72), (0.76, 0.48))]:
        ax_flow.annotate("", xy=stop, xytext=start,
                         arrowprops={"arrowstyle": "->", "lw": 1.1})
    ax_flow.annotate("", xy=(0.56, 0.37), xytext=(0.44, 0.37),
                     arrowprops={"arrowstyle": "<->", "lw": 1.2})
    ax_flow.annotate("", xy=(0.28, 0.23), xytext=(0.72, 0.23),
                     arrowprops={"arrowstyle": "->", "linestyle": "--", "lw": 1.0})
    ax_flow.text(0.50, 0.11, "Historical pollen-service causation NOT identified",
                 ha="center", fontsize=9, fontweight="bold")
    ax_flow.text(0.50, 0.02, "Solid arrows: fitted associations across separate samples.",
                 ha="center", fontsize=8)

    # Keep the four geographic classes explicitly decoded at actual print width.
    handles, labels = ax_island.get_legend_handles_labels()
    fig.legend(handles, labels, loc="center", bbox_to_anchor=(0.50, 0.568),
               ncol=4, fontsize=7.2, frameon=False)
    fig.suptitle("Figure 1 | Geographic sampling and inference boundary", y=0.995, fontsize=11)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.94, bottom=0.04)
    _save(fig, out, Q1_FIGURE_STEMS[1])


def q1_figure2(root: Path, out: Path) -> None:
    """Frozen H3 coefficients; resized for the 178-mm print-width audit."""
    h3 = pd.read_csv(
        root / "submission/chapter1_current/supplement/Table_S5_H3_pollen_limitation.csv"
    ).set_index("analysis")
    order = ["primary", "no_zero_constant", "supplemental_only"]
    labels = ["Primary", "No-zero constant", "Supplemental-only"]
    work = h3.loc[order]
    b, se = work["estimate"].to_numpy(float), work["se"].to_numpy(float)
    n_pubs = work["n_publications"].to_numpy(int)
    halfwidth = student_t.ppf(0.975, df=n_pubs - 1) * se
    fig, ax = plt.subplots(figsize=(7.0, 3.65))
    y = np.arange(len(order))
    ax.errorbar(b, y, xerr=halfwidth, fmt="o", capsize=3, markersize=5)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(y, [f"{name} (n={n:,})" for name, n in zip(labels, n_pubs, strict=True)])
    ax.invert_yaxis()
    for yi, be, p in zip(y, b, work["finite_publication_t_two_sided_p"].to_numpy(float),
                         strict=True):
        ax.annotate(f"beta={be:+.3f}; p={p:.3f}", xy=(be, yi),
                    xytext=(7, 10 if yi == len(y) - 1 else -12),
                    textcoords="offset points", fontsize=8)
    ax.set(xlabel="Isolation coefficient in GloPL log response ratio",
           title="Figure 2 | Isolation and experimental pollen limitation")
    ax.grid(axis="x", linewidth=0.4, alpha=0.3)
    fig.text(0.10, 0.025, "Publication-t 95% CIs. Supplemental-only result is unsupported.",
             fontsize=8)
    fig.tight_layout(rect=(0, 0.10, 1, 1))
    _save(fig, out, Q1_FIGURE_STEMS[2])


def q1_figure3(root: Path, out: Path) -> None:
    """Refacet the exact 112 frozen H1 fits for readable two-column print width.

    The historical original figure remains byte-for-byte preserved in
    results/h1_final_traitwise_t_20261004 and the six-source display folder.
    No model refit or inference transformation occurs here.
    """
    data = pd.read_csv(
        root / "results/h1_final_traitwise_t_20261004/traitwise_results.csv"
    )
    assert len(data) == 112 and set(data["outcome"]) == set(H1_OUTCOMES)
    scopes = [
        ("all", "broad", "All / broad", "o"),
        ("direct", "broad", "Direct / broad", "s"),
        ("all", "wcvp", "All / WCVP", "D"),
        ("direct", "wcvp", "Direct / WCVP", "^"),
    ]
    colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fig, axes = plt.subplots(4, 1, figsize=(7.0, 9.35), sharex=True)
    rows = np.arange(len(H1_OUTCOMES))
    offsets = [-0.27, -0.09, 0.09, 0.27]
    for ax, region in zip(axes, REGIONS, strict=True):
        r = data.loc[data["context"].eq(region)]
        for i, (evidence, flora, name, marker) in enumerate(scopes):
            part = r.loc[
                r["evidence_scope"].eq(evidence) & r["flora_scope"].eq(flora)
            ].set_index("outcome").loc[H1_OUTCOMES]
            assert len(part) == 7 and (part["status"].astype(str).eq("ok")).all()
            estimate = part["estimate"].to_numpy(float)
            errors = np.vstack((
                estimate - part["ci_low"].to_numpy(float),
                part["ci_high"].to_numpy(float) - estimate,
            ))
            ys = rows + offsets[i]
            supported = part["p_two_sided"].to_numpy(float) < 0.05
            color = colours[i]
            ax.errorbar(estimate, ys, xerr=errors, fmt="none", color=color,
                        linewidth=0.8, capsize=1.5)
            for j, (x, yy) in enumerate(zip(estimate, ys, strict=True)):
                ax.plot(x, yy, marker=marker, markersize=4.1,
                        markeredgecolor=color, markerfacecolor=color if supported[j] else "white",
                        linestyle="none")
        ax.axvline(0, color="black", lw=0.65)
        ax.set_yticks(rows, [H1_LABELS[k] for k in H1_OUTCOMES])
        ax.set_ylim(6.6, -0.62)
        ax.set_title(REGION_LABELS[region], fontsize=10, fontweight="bold", pad=2)
        ax.axhspan(-0.5, 2.5, facecolor="#f2f5f5", alpha=0.8)
        ax.grid(axis="x", linewidth=0.35, alpha=0.25)
        ax.tick_params(labelsize=8.3)
    axes[-1].set_xlabel("Isolation slope (log odds, standardized within trait / region / scope)")
    xmin = float(data["ci_low"].min()) - 0.04
    xmax = float(data["ci_high"].max()) + 0.04
    axes[-1].set_xlim(xmin, xmax)
    handles = [
        plt.Line2D([], [], linestyle="none", marker=m, markersize=5,
                   color=colours[i], label=label)
        for i, (_evidence, _flora, label, m) in enumerate(scopes)
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2,
               bbox_to_anchor=(0.5, 0.051), fontsize=8.3, frameon=False)
    fig.suptitle("Figure 3 | Seven traits in four geographic regions", fontsize=11, y=0.995)
    fig.text(0.51, 0.014,
             "Pointwise 95% cluster-t CIs. Filled = nominal two-sided p < 0.05 (unadjusted).",
             fontsize=8, ha="center")
    fig.tight_layout(rect=(0, 0.085, 1, 0.97), h_pad=0.62)
    _save(fig, out, Q1_FIGURE_STEMS[3])


def q1_figure4(root: Path, out: Path) -> None:
    """Three distinct H2 estimands, stacked rather than squeezing link scales."""
    h2 = pd.read_csv(
        root / "submission/chapter1_current/supplement/Table_S3a_H2_all_finite_cluster.csv"
    )
    responses = [
        ("selfing_core", "A | Selfing core (OLS)"),
        ("generalized_accessible", "B | Floral accessibility | selfing (OLS)"),
        ("plain_colour", "C | Plain colour | selfing (logit)"),
    ]
    fig, axes = plt.subplots(3, 1, figsize=(7.0, 7.25))
    y = np.arange(len(REGIONS))
    ylabels = ["N. midlatitude", "N. high latitude", "Tropical", "S. extratropical"]
    for ax, (response, title) in zip(axes, responses, strict=True):
        part = h2.loc[h2["response"].eq(response)].set_index("context").loc[REGIONS]
        estimate = part["distance_estimate"].to_numpy(float)
        se = part["distance_se"].to_numpy(float)
        halfwidth = student_t.ppf(
            0.975, df=part["n_clusters"].to_numpy(int) - 1
        ) * se
        ax.errorbar(estimate, y, xerr=halfwidth, fmt="o", capsize=3, markersize=5)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_yticks(y, ylabels)
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.2)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Isolation coefficient (response-specific scale)")
        if response == "generalized_accessible":
            qs = part["primary_H2b_q_finite_cluster_t"].to_numpy(float)
            for yi, q in zip(y, qs, strict=True):
                ax.text(0.99, yi, f"BH q={q:.3g}", va="center", ha="right",
                        transform=ax.get_yaxis_transform(), fontsize=8.2,
                        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"})
    fig.suptitle("Figure 4 | Assurance and independent floral associations", y=0.995, fontsize=11)
    fig.text(0.50, 0.018,
             "Response scales differ; H2 is conditional association, not historical mediation.",
             ha="center", fontsize=8.2)
    fig.tight_layout(rect=(0, 0.065, 1, 0.96), h_pad=0.9)
    _save(fig, out, Q1_FIGURE_STEMS[4])


def q1_figure5(root: Path, out: Path) -> None:
    """Non-paired GloPL matched sets with distinct sample sizes."""
    h4 = pd.read_csv(
        root / "submission/chapter1_current/supplement/Table_S6a_H4_scores.csv"
    )
    families = [
        ("selfing_core", "A | Reproductive assurance"),
        ("generalized_accessible", "B | Accessible floral architecture"),
    ]
    order = ["primary", "no_zero_constant", "supplemental_only"]
    fig, axes = plt.subplots(2, 1, figsize=(7.0, 5.45), sharex=True)
    for ax, (syndrome, title) in zip(axes, families, strict=True):
        part = h4.loc[h4["H2_syndrome"].eq(syndrome)].set_index("analysis").loc[order]
        estimate = part["estimate"].to_numpy(float)
        se = part["se"].to_numpy(float)
        n_pub = part["n_publications"].to_numpy(int)
        half = student_t.ppf(0.975, df=n_pub - 1) * se
        y = np.arange(len(order))
        ax.errorbar(estimate, y, xerr=half, fmt="o", capsize=3, markersize=5)
        ax.axvline(0, color="black", lw=0.8)
        labels = ["Primary", "No-zero constant", "Supplemental-only"]
        ax.set_yticks(y, [
            f"{label} (n={n:,} species)" for label, n in zip(
                labels, part["n_species"].to_numpy(int), strict=True
            )
        ])
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.22)
        ax.set_title(title, fontsize=10)
    axes[-1].set_xlabel("Trait-score association with current pollen limitation")
    fig.suptitle("Figure 5 | Exact-species functional consistency (post-hoc)", y=0.995, fontsize=11)
    fig.text(0.5, 0.023,
             "Different matching sets (455 vs 143 primary species); publication-t 95% CIs.",
             ha="center", fontsize=8.2)
    fig.tight_layout(rect=(0, 0.075, 1, 0.95), h_pad=0.7)
    _save(fig, out, Q1_FIGURE_STEMS[5])


def build_q1_manifest(out: Path) -> None:
    records = []
    for number, stem in Q1_FIGURE_STEMS.items():
        svg, pdf = out / f"{stem}.svg", out / f"{stem}.pdf"
        if not svg.is_file() or not pdf.is_file():
            raise FileNotFoundError(f"Missing Q1 Main figure {number}: {svg} / {pdf}")
        ET.parse(svg)
        if pdf.stat().st_size < 1000:
            raise ValueError(f"Suspiciously short Q1 PDF: {pdf}")
        records.append({
            "figure": number, "stem": stem,
            "svg": svg.name, "svg_bytes": svg.stat().st_size,
            "pdf": pdf.name, "pdf_bytes": pdf.stat().st_size,
        })
    (out / "FIGURE_MANIFEST.json").write_text(
        json.dumps({
            "contract": "chapter1_q1_five_main_figures_v1",
            "source_surface": "chapter1_submission_current_corrected_geography_final_traitwise",
            "inferential_boundary": "observational_three_stage_not_historical_mediation",
            "figures": records,
        }, indent=2) + "\n", encoding="utf-8"
    )



def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("submission/chapter1_current/figures"),
    )
    parser.add_argument("--q1-five", action="store_true", help="Render the five Q1 submission figures without changing the six historical source figures.")
    args = parser.parse_args()
    root = args.root.resolve()
    out = (root / args.output_dir).resolve() if not args.output_dir.is_absolute() else args.output_dir
    _style()
    if args.q1_five:
        q1_figure1(root, out)
        q1_figure2(root, out)
        q1_figure3(root, out)
        q1_figure4(root, out)
        q1_figure5(root, out)
        build_q1_manifest(out)
        return
    figure1(root, out)
    figure2(root, out)
    figure3(root, out)
    figure4(root, out)
    figure5(root, out)
    figure6(root, out)
    build_manifest(out)


if __name__ == "__main__":
    main()
