from __future__ import annotations

import argparse
import json
from pathlib import Path
from xml.etree import ElementTree as ET

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch

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



def figure4(root: Path, out: Path) -> None:
    slopes = pd.read_csv(root / "results/geography_20260924/all/beta_binomial_within_slopes.csv")
    omnibus = pd.read_csv(root / "results/geography_20260924/all/beta_binomial_within_omnibus.csv")
    slopes = slopes.loc[slopes["stratum"].eq("all_observed")].copy()
    omnibus = omnibus.loc[omnibus["stratum"].eq("all_observed")].set_index("context")

    fig, axes = plt.subplots(1, 4, figsize=(13.2, 5.3), sharey=True)
    y = np.arange(len(H1_OUTCOMES))
    for ax, region in zip(axes, REGIONS, strict=True):
        part = slopes.loc[slopes["context"].eq(region)].set_index("outcome")
        est = np.array([float(part.loc[o, "geography_slope_log_odds"]) for o in H1_OUTCOMES])
        se = np.array([float(part.loc[o, "cluster_robust_se"]) for o in H1_OUTCOMES])
        ax.errorbar(est, y, xerr=1.96 * se, fmt="o", capsize=2)
        ax.axvline(0, linewidth=0.8)
        q = float(omnibus.loc[region, "q_value"])
        ax.set_title(f"{REGION_LABELS[region]}\njoint q={q:.2g}")
        ax.set_xlabel("Isolation coefficient")
        ax.grid(axis="x", linewidth=0.3, alpha=0.4)
    axes[0].set_yticks(y, [H1_LABELS[o] for o in H1_OUTCOMES])
    axes[0].invert_yaxis()
    fig.suptitle("Figure 4 | H1: recurrent multivariate response with regional expression", y=1.02)
    fig.tight_layout()
    _save(fig, out, "Figure4_H1_recurrent_multivariate_response")


def figure5(root: Path, out: Path) -> None:
    h2 = pd.read_csv(root / "results/geography_20260924/all/h2_decomposition_models.csv")
    responses = ["selfing_core", "generalized_accessible", "plain_colour"]
    labels = ["Reproductive assurance", "Accessibility | assurance", "Plain colour | assurance"]
    fig, axes = plt.subplots(1, 4, figsize=(13.0, 4.3), sharey=True)
    y = np.arange(len(responses))
    for ax, region in zip(axes, REGIONS, strict=True):
        part = h2.loc[h2["context"].eq(region)].set_index("response")
        est = np.array([float(part.loc[r, "distance_estimate"]) for r in responses])
        se = np.array([float(part.loc[r, "distance_se"]) for r in responses])
        ax.errorbar(est, y, xerr=1.96 * se, fmt="o", capsize=3)
        ax.axvline(0, linewidth=0.8)
        access_q = float(part.loc["generalized_accessible", "primary_H2b_q"])
        ax.set_title(f"{REGION_LABELS[region]}\naccessibility q={access_q:.3g}")
        ax.set_xlabel("Isolation coefficient")
        ax.grid(axis="x", linewidth=0.3, alpha=0.4)
    axes[0].set_yticks(y, labels)
    axes[0].invert_yaxis()
    fig.suptitle("Figure 5 | H2: reproductive assurance and additional floral response", y=1.02)
    fig.tight_layout()
    _save(fig, out, "Figure5_H2_conditional_decomposition")


def figure6(root: Path, out: Path) -> None:
    h3 = pd.read_csv(
        root / "submission/chapter1_current/supplement/Table_S5_H3_pollen_limitation.csv"
    )
    h4 = pd.read_csv(root / "results/geography_20260924/h4_exact_corrected.csv")
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
                "source_surface": "corrected_geography_20260924",
                "figures": records,
            },
            indent=2,
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("submission/chapter1_current/figures"),
    )
    args = parser.parse_args()
    root = args.root.resolve()
    out = (root / args.output_dir).resolve() if not args.output_dir.is_absolute() else args.output_dir
    _style()
    figure1(root, out)
    figure2(root, out)
    figure3(root, out)
    figure4(root, out)
    figure5(root, out)
    figure6(root, out)
    build_manifest(out)


if __name__ == "__main__":
    main()
