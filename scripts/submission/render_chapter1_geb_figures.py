"""GEB-only figure rendering from the frozen Chapter 1 source tables.

Replots global point locations on a Mollweide equal-area map and changes panel
title lettering from A/B/C to (a)/(b)/(c) for the other live figures.
No statistical re-fitting. Original Ecology Letters figure assets are untouched.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
from pathlib import Path

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import render_chapter1_main_figures as base  # noqa: E402

PROJECTION = "Mollweide"
EQUATOR_5KKM_DEGREES = 5000.0 / 111.195


def _rad_lon(series: pd.Series) -> np.ndarray:
    v = pd.to_numeric(series, errors="coerce").to_numpy(float)
    v = ((v + 180.0) % 360.0) - 180.0
    return np.deg2rad(v)


def _rad_lat(series: pd.Series) -> np.ndarray:
    v = pd.to_numeric(series, errors="coerce").to_numpy(float)
    return np.deg2rad(np.clip(v, -89.999, 89.999))


def _equatorial_scale(ax: matplotlib.axes.Axes) -> None:
    # Great-circle length near equator on the mean-radius sphere, not at 30/60 deg.
    start = -158.0
    end = start + EQUATOR_5KKM_DEGREES
    ax.plot(np.deg2rad([start, end]), [0.0, 0.0], lw=2.5, solid_capstyle="butt")
    ax.text(np.deg2rad(start + 3.0), np.deg2rad(-12.0), "5,000 km at equator",
            fontsize=7, ha="left", va="top")


def render_equal_area_figure1(root: Path, out: Path) -> None:
    cov = pd.read_csv(root / "results/geography_20260924/corrected_geography_covariates.csv")
    sites = pd.read_csv(root / "results/geography_20260924/glopl_corrected_site_distances.csv")
    summary = pd.read_csv(root / "submission/chapter1_current/supplement/Table_S1_data_summary.csv")
    formerly_zero = pd.read_csv(root / "results/geography_20260924/formerly_zero_islands_recalculated.csv")
    lat_col, lon_col = base._coord_columns(cov)
    assert "analysis_regime" in cov, "Geographic regime missing from frozen covariates"

    fig = plt.figure(figsize=(12.9, 8.7))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 0.85])

    ax1 = fig.add_subplot(gs[0, 0], projection="mollweide")
    for region in base.REGIONS:
        part = cov.loc[cov["analysis_regime"].astype(str).eq(region)]
        ax1.scatter(_rad_lon(part[lon_col]), _rad_lat(part[lat_col]),
                    s=6, alpha=0.6, label=base.REGION_LABELS[region])
    ax1.grid(linewidth=0.3, alpha=0.4)
    ax1.set_title("(a) Global island analysis units", fontsize=10)
    _equatorial_scale(ax1)
    ax1.legend(loc="lower center", bbox_to_anchor=(0.5, -0.25),
               frameon=False, ncol=2, fontsize=7)

    ax2 = fig.add_subplot(gs[0, 1], projection="mollweide")
    cont = sites["on_seeded_continent"].astype(str).str.lower().isin({"true", "1"})
    for part, lab, marker_size in [
        (sites.loc[cont], "Continental GloPL sites", 8),
        (sites.loc[~cont], "Offshore GloPL sites", 13),
    ]:
        ax2.scatter(_rad_lon(part["longitude"]), _rad_lat(part["latitude"]),
                    s=marker_size, alpha=0.58, label=lab)
    ax2.grid(linewidth=0.3, alpha=0.4)
    ax2.set_title("(b) Experimental pollen-supplementation sites", fontsize=10)
    _equatorial_scale(ax2)
    ax2.legend(loc="lower center", bbox_to_anchor=(0.5, -0.22),
               frameon=False, ncol=2, fontsize=7)

    ax3 = fig.add_subplot(gs[1, 0])
    metrics = {
        row["metric"]: float(row["value"]) for _, row in summary.iterrows()
        if row["section"] in {"traits", "GloPL"}
    }
    labs = ["Species", "Resolved trait cells", "GloPL rows", "GloPL sites"]
    vals = [
        metrics["accepted_angiosperm_species"],
        metrics["resolved_species_axis_cells"],
        metrics["effect_rows"],
        metrics["sites"],
    ]
    y = np.arange(len(labs))
    ax3.barh(y, vals)
    ax3.set_yticks(y, labs)
    ax3.invert_yaxis()
    ax3.set_xscale("log")
    ax3.set_xlabel("Count (log scale)")
    ax3.set_title("(c) Dataset coverage")
    for yi, value in zip(y, vals, strict=True):
        ax3.text(value * 1.08, yi, f"{int(value):,}", va="center", fontsize=8)

    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis("off")
    med = float(formerly_zero["spherical_coast_distance_km"].median())
    lines = [
        "(d) Geographic exposure audit",
        "",
        "8,264 corrected island units",
        "4,379 sampled H1 islands",
        "1,113 spurious island zero distances repaired",
        f"median repaired distance: {med:.3f} km",
        "one continental split component excluded",
        "996 actual continental GloPL zeros retained",
        "",
        "Source: GSHHG 2.3.7 matched coastlines",
        "Distance: great-circle coastline separation",
        "Map projection: Mollweide (equal area)",
        "Equatorial scale: 5,000 km",
    ]
    ax4.text(0.03, 0.98, "\n".join(lines), fontsize=9, va="top")
    fig.suptitle("Figure 1 | Worldwide scope and corrected isolation", y=1.01)
    fig.tight_layout()
    base._save(fig, out, "Figure1_global_scope_corrected")


def render(root: Path, out: Path) -> dict:
    root = root.resolve()
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    base._style()
    render_equal_area_figure1(root, out)

    old_title = matplotlib.axes.Axes.set_title

    def panel_title(self, label, *args, **kwargs):
        if isinstance(label, str):
            match = re.match(r"^([A-D])\s*\|\s*(.+)$", label)
            if match:
                label = f"({match.group(1).lower()}) {match.group(2)}"
        return old_title(self, label, *args, **kwargs)

    matplotlib.axes.Axes.set_title = panel_title
    try:
        base.figure2(root, out)
        base.figure3(root, out)
        # Figure 4 is frozen H1 output: copied unchanged, without refitting.
        for ext in ("svg", "pdf"):
            source = root / "results/h1_final_traitwise_t_20261004" / ("traitwise_H1." + ext)
            assert source.is_file(), source
            shutil.copyfile(source, out / ("Figure4_H1_recurrent_multivariate_response." + ext))
        base.figure5(root, out)
        base.figure6(root, out)
    finally:
        matplotlib.axes.Axes.set_title = old_title
    base.build_manifest(out)

    audit = {
        "contract": "chapter1_geb_equal_area_figures_v1",
        "scientific_refit": False,
        "input_analysis_selector": "config/chapter1_submission_current.json",
        "figure1_projection": PROJECTION,
        "figure1_equal_area": True,
        "equator_scale_km": 5000,
        "figure4_frozen_copy": True,
        "n_figures": 6,
        "figure6_offshore_posthoc": True,
        "source_original_figure_assets_unchanged": True,
    }
    (out / "GEB_FIGURE_AUDIT.json").write_text(json.dumps(audit, indent=2) + "\n")
    return audit


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=Path("."))
    p.add_argument("--output-dir", type=Path, required=True)
    opts = p.parse_args()
    print(json.dumps(render(opts.root, opts.output_dir), indent=2))


if __name__ == "__main__":
    main()
