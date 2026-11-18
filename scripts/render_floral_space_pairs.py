#!/usr/bin/env python3
"""Render aggregate pair-sharing results and their regional permutation references.

Only permutation_summary.csv and manifest.json are read. Missing comparisons
remain visible, and no observed/null ratio is required (a zero null is valid).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd


REGIONS = (
    "northern_high_latitude", "northern_midlatitude", "tropical", "southern_extratropical",
)
REGION_LABELS = (
    "Northern\nhigh latitudes", "Northern\nmidlatitudes", "Tropics", "Southern\nextratropics",
)
SCOPES = ("all", "direct", "all_single_state", "direct_single_state")
SETTINGS = (("all", "All", -0.16, "#287E8D", "o"),
            ("direct", "Direct", 0.16, "#BC7042", "s"))
STATISTICS = (
    "observed_median_jaccard", "null_median", "null_q025", "null_q975",
    "p_two_sided", "q_bh_16", "q_by_16",
)


def load_results(root: Path) -> tuple[pd.DataFrame, dict]:
    """Reject ambiguous/incomplete summaries instead of silently plotting zeros."""
    table = pd.read_csv(root / "permutation_summary.csv")
    parameters = json.loads((root / "manifest.json").read_text())["parameters"]
    required = {"evidence_scope", "analysis_regime", "status", "n_pairs", "n_islands",
                "n_permutations", *STATISTICS}
    if missing := required - set(table):
        raise ValueError(f"Missing summary columns: {sorted(missing)}")
    expected = {(scope, regime) for scope in SCOPES for regime in REGIONS}
    keys = list(zip(table.evidence_scope, table.analysis_regime, strict=True))
    if len(keys) != len(expected) or set(keys) != expected:
        raise ValueError("Summary must contain each of the 16 scope/regime comparisons exactly once")
    for key in ("permutations", "min_species", "min_unique", "max_pairs", "seed"):
        value = parameters[key]
        minimum = 0 if key == "seed" else 1
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise ValueError(f"Manifest parameter {key} must be an integer >= {minimum}")
    for column in ("n_pairs", "n_islands", "n_permutations"):
        values = pd.to_numeric(table[column], errors="raise")
        if not np.isfinite(values).all() or (values < 0).any() or (values % 1 != 0).any():
            raise ValueError(f"{column} must contain nonnegative integer counts")
        table[column] = values.astype(int)
    for column in STATISTICS:
        table[column] = pd.to_numeric(table[column], errors="raise")
    if not table.status.isin(["ok", "no_eligible_pairs"]).all():
        raise ValueError("Only completed ('ok') or 'no_eligible_pairs' comparisons can be plotted")
    complete = table.loc[table.status.eq("ok")]
    values = complete[list(STATISTICS)].to_numpy()
    if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
        raise ValueError("Completed comparisons require finite Jaccard and p/q values in [0, 1]")
    if ((complete.n_pairs < 1) | (complete.n_islands < 2)
            | (complete.n_pairs > parameters["max_pairs"])
            | (complete.n_permutations != parameters["permutations"])).any():
        raise ValueError("Completed support/permutation counts contradict manifest parameters")
    if ((complete.null_q025 > complete.null_median)
            | (complete.null_median > complete.null_q975)).any():
        raise ValueError("Null median must lie within the reported 95% null range")
    empty = table.loc[table.status.eq("no_eligible_pairs")]
    if (empty.n_pairs != 0).any() or (empty.n_permutations != 0).any():
        raise ValueError("Missing comparisons must have zero pairs and zero permutations")
    if empty[list(STATISTICS[:4])].notna().any().any():
        raise ValueError("Missing comparisons must not contain invented Jaccard statistics")
    for suffix in ("", "_single_state"):
        all_rows = table.loc[table.evidence_scope.eq("all" + suffix)].set_index("analysis_regime")
        direct = table.loc[table.evidence_scope.eq("direct" + suffix)].set_index("analysis_regime")
        columns = ["n_pairs", "n_islands", "status"]
        if not all_rows[columns].sort_index().equals(direct[columns].sort_index()):
            raise ValueError("All and Direct must share pair support within each encoding")
    return table.set_index(["evidence_scope", "analysis_regime"]), parameters


def q_label(value: float) -> str:
    """Show small adjusted values without rounding them to zero."""
    if not np.isfinite(value):
        return "NA"
    return f"{value:.3g}"


def render(root: Path) -> None:
    table, parameters = load_results(root)
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11, "axes.titlesize": 13,
        "axes.titleweight": "bold", "axes.spines.top": False, "axes.spines.right": False,
        "svg.fonttype": "none", "svg.hashsalt": "floral-space-pairs", "savefig.facecolor": "white",
    })
    plotted = table.loc[table.status.eq("ok"), list(STATISTICS[:4])].to_numpy()
    if plotted.size:
        low, high = float(plotted.min()), float(plotted.max())
        padding = max(0.025, (high - low) * 0.10)
        limits = (max(0, low - padding), min(1, high + padding))
    else:
        limits = (0, 1)

    fig = plt.figure(figsize=(19, 9.4))
    panels = fig.add_gridspec(1, 2, left=0.105, right=0.985, top=0.77, bottom=0.30, wspace=0.13)
    fig.suptitle("Floral-profile sharing between islands after removing shared taxa",
                 fontsize=19, weight="bold", y=0.97)
    fig.text(0.5, 0.918,
             "Observed median Jaccard similarity and a regional, within-genus permutation reference",
             ha="center", fontsize=12, color="#4A5056")
    fig.text(0.5, 0.882,
             f"{parameters['permutations']:,} profile-label permutations per completed comparison; "
             "one assignment is reused across all sampled pairs in its regime",
             ha="center", fontsize=10.5, color="#4A5056")

    for panel, suffix, title in zip(panels, ("", "_single_state"),
                                    ("A  All recorded state sets", "B  One recorded state per trait"),
                                    strict=True):
        cells = panel.subgridspec(1, 2, width_ratios=(1.25, 1), wspace=0.06)
        ax, annotation = fig.add_subplot(cells[0]), fig.add_subplot(cells[1])
        ax.set_title(title, loc="left", pad=30)
        annotation.set_title("Sampled support · BY q", loc="left", fontsize=10.5, pad=13)
        for scope, _, offset, color, marker in SETTINGS:
            for index, regime in enumerate(REGIONS):
                row = table.loc[(scope + suffix, regime)]
                y = index + offset
                if row.status == "ok":
                    ax.plot([row.null_q025, row.null_q975], [y, y], color=color,
                            linewidth=6, alpha=0.25, solid_capstyle="round", zorder=1)
                    ax.plot(row.null_median, y, marker="|", color=color, markersize=15,
                            markeredgewidth=1.7, linestyle="", zorder=2, clip_on=False)
                    ax.plot(row.observed_median_jaccard, y, marker=marker, linestyle="",
                            color=color, markerfacecolor=color if row.q_by_16 < 0.05 else "white",
                            markersize=8, markeredgewidth=1.8, zorder=3, clip_on=False)
                    label = (f"{row.n_pairs:,} pairs · {row.n_islands:,} islands"
                             f" · q={q_label(row.q_by_16)}")
                else:
                    ax.text(0.03, y, "No eligible pairs", transform=ax.get_yaxis_transform(),
                            color=color, fontsize=10, va="center")
                    label = f"0 pairs · {row.n_islands:,} islands · q=NA"
                annotation.text(0, y, label, color=color, fontsize=9.5, va="center")
        for boundary in (0.5, 1.5, 2.5):
            ax.axhline(boundary, color="#EEEEEE", linewidth=0.8, zorder=0)
            annotation.axhline(boundary, color="#EEEEEE", linewidth=0.8, zorder=0)
        ax.set(xlim=limits, ylim=(3.5, -0.5), xlabel="Median profile Jaccard similarity")
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=3))
        ax.set_yticks(range(4), REGION_LABELS if not suffix else [""] * 4)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0, pad=10, labelsize=10)
        ax.tick_params(axis="x", labelsize=10)
        annotation.set(ylim=(3.5, -0.5), xlim=(0, 1))
        annotation.set_axis_off()

    handles = [Line2D([0], [0], color=color, marker=marker, linestyle="", markersize=8,
                      markerfacecolor="white", label=f"Observed: {label}")
               for _, label, _, color, marker in SETTINGS]
    handles.extend([
        Line2D([0], [0], color="#7A838A", linewidth=6, alpha=0.3, label="Null 95% range"),
        Line2D([0], [0], color="#7A838A", marker="|", markersize=15, linestyle="",
               markeredgewidth=1.7, label="Null median"),
        Line2D([0], [0], color="#525B62", marker="o", linestyle="", markersize=8,
               label="Filled observed: BY q < 0.05"),
    ])
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.205),
               ncol=5, frameon=False, fontsize=10, columnspacing=2)
    notes = [
        "Horizontal bands are central 95% permutation null ranges, not effect confidence intervals. "
        "Outlined observed points: BY q ≥ 0.05.",
        "q values use BY adjustment across the fixed 16-comparison family; BH values are retained "
        "in the results table. Missing comparisons remain in that family.",
        f"Support: both islands have ≥{parameters['min_species']:,} complete taxon-name records "
        f"and retain ≥{parameters['min_unique']:,} nonshared names in both scopes; "
        f"at most {parameters['max_pairs']:,} pairs per regime and encoding.",
        "All and Direct share sampled pair IDs within each encoding. Encodings can use different "
        "pairs; overlapping scopes are not independent replication.",
        "Each sampled pair has equal weight. The regional genus reference fixes profile "
        "frequencies per taxon; differences do not establish an evolutionary mechanism.",
        f"Permutation and pair-selection seed: {parameters['seed']}. "
        "Both panels use the same raw Jaccard axis scale.",
    ]
    for y, note in zip(np.linspace(0.18, 0.035, len(notes)), notes, strict=True):
        fig.text(0.105, y, note, fontsize=10, color="#4A5056", va="top")
    fig.savefig(root / "pair_sharing_results.png", dpi=200)
    svg = root / "pair_sharing_results.svg"
    fig.savefig(svg, metadata={"Date": None})
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path("results/floral_space_pairs_20261106"))
    render(parser.parse_args().results)
