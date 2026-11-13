"""Render reviewed aggregate outputs; all plotted points are observed or derived data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter, ScalarFormatter
import numpy as np
import pandas as pd


REGIONS = ["northern_high_latitude", "northern_midlatitude", "tropical", "southern_extratropical"]
REGION_LABELS = ["Northern\nhigh latitudes", "Northern\nmidlatitudes", "Tropics", "Southern\nextratropics"]
SCOPES = ["all", "direct", "all_single_state", "direct_single_state"]


def check_report_parameters(root: Path) -> None:
    """This dated report has fixed comparison cohorts and labels."""
    for scope in SCOPES:
        summary = json.loads((root / scope / "summary.json").read_text())
        if summary["min_species"] != 50 or summary["rarefaction_size"] != 50:
            raise ValueError("This report renderer requires min_species=50 and rarefaction_size=50")
    if (root / "permutation_summary.csv").exists():
        manifest = json.loads((root / "permutation_manifest.json").read_text())
        if (manifest["minimum_complete_name_units_per_island"] != 50
                or manifest["permutations_requested_per_comparison"] != 999):
            raise ValueError("This report renderer requires 50-name eligibility and 999 permutations")


def save_figure(fig: plt.Figure, root: Path, name: str) -> None:
    fig.savefig(root / f"{name}.png", dpi=180)
    svg_path = root / f"{name}.svg"
    fig.savefig(svg_path)
    svg_path.write_text("\n".join(line.rstrip() for line in svg_path.read_text().splitlines()) + "\n")


def render(root: Path) -> None:
    check_report_parameters(root)
    tables = {s: pd.read_csv(root / s / "island_metrics.csv") for s in SCOPES}
    summaries = {s: json.loads((root / s / "summary.json").read_text()) for s in SCOPES}
    a = tables["all"].query("n_species >= 50").copy()
    paired = tables["all"].merge(tables["direct"], on="island_id", suffixes=("_all", "_direct"))
    paired = paired.loc[paired.n_species_all.ge(50) & paired.n_species_direct.ge(50)].copy()
    paired.to_csv(root / "paired_all_direct_islands.csv.gz", index=False)
    support = set.intersection(*[set(t.loc[t.n_species.ge(50), "island_id"]) for t in tables.values()])
    scope_rows = []
    for scope, table in tables.items():
        eligible = table.loc[table.n_species.ge(50)]
        common = eligible.loc[eligible.island_id.isin(support)]
        for support_name, part in [("scope_eligible", eligible), ("four_scope_common", common)]:
            scope_rows.append({
                "scope": scope, "support": support_name, "n_islands": len(part),
                "n_complete_taxon_names": summaries[scope]["n_complete_species"],
                "n_observed_profiles": summaries[scope]["n_profiles"],
                **{f"median_{col}": float(part[col].median()) for col in [
                    "n_species", "n_profiles", "singleton_fraction", "rarefied_profiles",
                    "profile_ratio_to_genus_pool", "joint_coverage",
                ]},
            })
    pd.DataFrame(scope_rows).to_csv(root / "scope_comparison.csv", index=False)
    comparison = {
        "paired_all_direct_islands": len(paired), "four_scope_common_islands": len(support),
        "rarefaction_size": 50,
        "paired_rarefied_richness_spearman": float(paired.rarefied_profiles_all.corr(
            paired.rarefied_profiles_direct, method="spearman")),
        "paired_median_rarefied_richness_difference_all_minus_direct": float(
            (paired.rarefied_profiles_all - paired.rarefied_profiles_direct).median()),
        "paired_singleton_fraction_median_all": float(paired.singleton_fraction_all.median()),
        "paired_singleton_fraction_median_direct": float(paired.singleton_fraction_direct.median()),
        "interpretation": "Matched islands and equal taxon sample size; missingness and imputation remain observational differences, not a causal treatment.",
    }
    (root / "matched_scope_summary.json").write_text(json.dumps(comparison, indent=2) + "\n")
    a["richness_group"] = pd.cut(a.n_species, [49, 99, 249, 499, 999, np.inf],
                                 labels=["50-99", "100-249", "250-499", "500-999", "1000+"])
    stratified = a.groupby("richness_group", observed=True).agg(
        n_islands=("island_id", "size"), median_taxon_names=("n_species", "median"),
        median_profiles=("n_profiles", "median"),
        median_singleton_fraction=("singleton_fraction", "median"),
        median_mean_redundancy=("mean_redundancy", "median"),
    ).reset_index()
    stratified.to_csv(root / "richness_strata.csv", index=False)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titleweight": "bold", "axes.titlesize": 12,
                         "svg.fonttype": "none", "savefig.facecolor": "white"})
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), layout="constrained")
    fig.suptitle("Island floral profiles: variety, shared types and their carriers", fontsize=19, weight="bold")
    fig.get_layout_engine().set(h_pad=0.16, w_pad=0.13, hspace=0.09, wspace=0.08)

    ax = axes[0, 0]
    scatter = ax.scatter(a.n_species, a.n_profiles, c=a.singleton_fraction, cmap="viridis",
                         vmin=0, vmax=1, s=17, alpha=0.8, linewidths=0, rasterized=True)
    ax.set(xscale="log", xlabel="Taxon names with complete colour and form", ylabel="Distinct recorded profiles")
    ax.set_title(f"A  {len(a):,} islands with at least 50 complete records", loc="left", pad=12)
    ax.set_xticks([50, 100, 250, 500, 1000, 2500])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_ylim(bottom=0)
    colorbar = fig.colorbar(scatter, ax=ax, fraction=0.043, pad=0.02)
    colorbar.set_label("Profiles with only one carrier")
    colorbar.ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.text(0.04, 0.94, f"Median singleton share: {a.singleton_fraction.median():.1%}",
            transform=ax.transAxes, va="top", bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9})

    ax = axes[0, 1]
    pair = pd.read_csv(root / "all/pairwise_comparisons.csv")
    ax.scatter(pair.species_jaccard, pair.profile_jaccard, s=10, alpha=0.20, color="#AAB2BA", linewidths=0,
               label="All carriers", rasterized=True)
    ax.scatter(pair.species_jaccard, pair.profile_jaccard_without_shared_species, s=10, alpha=0.35,
               color="#287E8D", linewidths=0, label="Shared taxa removed", rasterized=True)
    ax.plot([0, 1], [0, 1], color="#7A7A7A", linewidth=0.9, linestyle="--")
    ax.set(xlim=(-0.01, 1), ylim=(-0.01, 1), xlabel="Taxonomic overlap before removal (Jaccard)", ylabel="Profile overlap (Jaccard)")
    ax.set_title("B  Profile sharing persists after removing shared taxa", loc="left", pad=12)
    ax.legend(frameon=False, loc="upper left", markerscale=2)
    ax.text(0.97, 0.04,
            f"{len(pair):,} sampled island pairs\nMedian after removal: {pair.profile_jaccard_without_shared_species.median():.1%}\nDescriptive comparison; coarse categories increase similarity",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.92})

    ax = axes[1, 0]
    profiles = pd.read_csv(root / "all/profile_coordinates.csv")
    diagnostic = pd.read_csv(root / "all/ordination_diagnostics.csv")
    fractions = np.diff(np.r_[0, diagnostic.positive_inertia_fraction.to_numpy()])
    palette = {"white": "#C8C8C3", "green_brown_inconspicuous": "#708143",
               "yellow_orange": "#DDA432", "red_pink": "#C65378", "blue_purple": "#6A5EB3"}
    colour_labels = {"white": "White", "green_brown_inconspicuous": "Green / brown",
                     "yellow_orange": "Yellow / orange", "red_pink": "Red / pink", "blue_purple": "Blue / purple"}
    profiles["colour_group"] = profiles.flower_primary_color.where(
        ~profiles.flower_primary_color.str.contains("|", regex=False), "multiple")
    for group in ["multiple", *palette]:
        part = profiles.loc[profiles.colour_group.eq(group)]
        ax.scatter(part.PCoA1, part.PCoA2, s=14 + 3*np.sqrt(part.n_species_scope),
                   c=palette.get(group, "#62727C"), alpha=0.68, edgecolors="white", linewidths=0.35,
                   label=colour_labels.get(group, "Multiple recorded colours"))
    ax.axhline(0, color="#EEEEEE", zorder=0)
    ax.axvline(0, color="#EEEEEE", zorder=0)
    ax.set_xlabel(f"PCoA 1 ({fractions[0]:.1%})")
    ax.set_ylabel(f"PCoA 2 ({fractions[1]:.1%})")
    ax.set_title("C  Shared PCoA coordinates: a two-axis projection", loc="left", pad=12)
    ax.text(0.03, 0.97, f"{len(profiles):,} recorded profiles\nTwo axes represent {fractions[:2].sum():.1%} of inertia",
            transform=ax.transAxes, va="top", fontsize=10,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9})
    handles = [Line2D([0], [0], marker="o", linestyle="", color=palette.get(g, "#62727C"),
                      label=colour_labels.get(g, "Multiple colours"), markersize=6)
               for g in [*palette, "multiple"]]
    ax.legend(handles=handles, loc="lower left", frameon=True, facecolor="white",
              edgecolor="none", framealpha=0.92, ncol=2, fontsize=8)

    ax = axes[1, 1]
    for scope, offset, color, marker in [("all", -0.11, "#287E8D", "o"), ("direct", 0.11, "#BC7042", "s")]:
        for i, region in enumerate(REGIONS):
            values = paired.loc[paired.analysis_regime_all.eq(region), f"profile_ratio_to_genus_pool_{scope}"]
            q25, median, q75 = values.quantile([0.25, 0.5, 0.75])
            ax.errorbar(i + offset, median, yerr=[[median-q25], [q75-median]], fmt=marker,
                        color=color, ecolor=color, capsize=4, markersize=6, linewidth=1.6)
    ax.axhline(1, color="#777777", linewidth=1, linestyle="--")
    ax.set_xticks(range(4), REGION_LABELS, fontsize=10)
    ax.set_ylabel("Observed / genus-fixed expected profile count")
    ax.set_ylim(0.80, 1.20)
    ax.set_title("D  Richness relative to a genus-fixed finite pool", loc="left", pad=12)
    ax.legend(handles=[Line2D([0], [0], marker="o", color="#287E8D", linestyle="", label="All"),
                       Line2D([0], [0], marker="s", color="#BC7042", linestyle="", label="Direct")],
              frameon=False, loc="upper right")
    ax.text(0.03, 0.04, f"Same {len(paired):,} islands · median and interquartile range\nPools: all recorded taxa in each broad latitudinal regime",
            transform=ax.transAxes, fontsize=9, va="bottom")
    for ax in axes.flat:
        ax.tick_params(labelsize=10)
    save_figure(fig, root, "floral_space_overview")
    plt.close(fig)

    if (root / "permutation_summary.csv").exists():
        render_permutations(root)


def render_permutations(root: Path) -> None:
    """Compare observed statistics with their own randomization distributions."""
    check_report_parameters(root)
    results = pd.read_csv(root / "permutation_summary.csv")
    if len(results) != 16 or not results.status.eq("ok").all():
        raise ValueError("Permutation figure requires all 16 completed comparisons")
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.8), sharey=True)
    fig.subplots_adjust(left=0.115, right=0.985, top=0.77, bottom=0.24, wspace=0.12)
    fig.suptitle("Single-state profiles are concentrated in northern midlatitude islands",
                 fontsize=18, weight="bold", y=0.98)
    fig.text(0.5, 0.916,
             "Observed richness relative to a genus-fixed reference; 999 shared-taxon permutations per comparison",
             ha="center", fontsize=11, color="#4A5056")
    settings = [("all", "All", -0.15, "#287E8D", "o"),
                ("direct", "Direct", 0.15, "#BC7042", "s")]
    for ax, suffix, title in zip(axes, ["", "_single_state"],
                                 ["A  All recorded state sets", "B  One recorded state per trait"], strict=True):
        for scope, label, offset, color, marker in settings:
            for i, region in enumerate(REGIONS):
                row = results.loc[results.evidence_scope.eq(scope + suffix)
                                  & results.analysis_regime.eq(region)].iloc[0]
                y = i + offset
                ax.plot([row.null_q025, row.null_q975], [y, y], color=color,
                        linewidth=5, alpha=0.22, solid_capstyle="round", zorder=1)
                ax.plot(row.null_median, y, marker="|", color=color, markersize=12,
                        markeredgewidth=1.5, zorder=2)
                ax.plot(row.observed_ratio, y, marker=marker, linestyle="",
                        color=color, markerfacecolor=color if row.q_bh_16 < 0.05 else "white",
                        markersize=8, markeredgewidth=1.7, zorder=3)
                ax.text(1.077, y, f"n={int(row.n_islands):,}  q={row.q_bh_16:.3f}",
                        va="center", fontsize=8.7, color=color)
        ax.axvline(1, color="#A3A3A3", linewidth=0.8, linestyle="--", zorder=0)
        ax.set(xlim=(0.855, 1.15), ylim=(3.45, -0.45))
        ax.set_xticks([0.9, 0.95, 1, 1.05], ["0.90", "0.95", "1.00", "1.05"])
        ax.set_title(title, loc="left", fontsize=13, pad=18)
        ax.set_xlabel("Median observed / expected profile count", fontsize=11, labelpad=10)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0, pad=10)
        for boundary in [0.5, 1.5, 2.5]:
            ax.axhline(boundary, color="#EEEEEE", linewidth=0.8, zorder=0)
    axes[0].set_yticks(range(4), REGION_LABELS, fontsize=10)
    handles = [Line2D([0], [0], color=color, marker=marker, linestyle="", markersize=7,
                      label=label) for _, label, _, color, marker in settings]
    handles += [Line2D([0], [0], color="#A5A9AC", linewidth=5, alpha=0.5,
                       label="95% permutation range"),
                Line2D([0], [0], color="#70777C", marker="|", linestyle="", markersize=12,
                       label="Permutation median")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.1),
               ncol=4, frameon=False, fontsize=10)
    fig.text(0.115, 0.07,
             "Filled points: BH q < 0.05 across 16 exploratory comparisons. Bands describe the null distribution, not effect confidence intervals.",
             fontsize=9.5, color="#4A5056")
    fig.text(0.115, 0.035,
             "n = islands with at least 50 complete taxon-name records. All and Direct overlap; single-state filtering changes the analyzed records.",
             fontsize=9.5, color="#4A5056")
    save_figure(fig, root, "permutation_results")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path("results/floral_space_20261106"))
    render(parser.parse_args().results)
