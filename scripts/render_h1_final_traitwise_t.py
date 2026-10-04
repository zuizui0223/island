"""Render all traitwise estimates directly from the committed replay table."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

root = Path(__file__).resolve().parents[1]
p = root / "results/h1_final_traitwise_t_20261004"
d = pd.read_csv(p / "traitwise_results.csv")
labels = [
    "Self-compatible",
    "Selfing mating system",
    "Autonomous / delayed selfing",
    "Plain colour",
    "Generalized floral form",
    "Radial symmetry",
    "Shallow / open tube",
]
traits = d.outcome.unique()
regions = d.context.unique()
scopes = [
    ("all", "broad", "#156477", "o", -0.24),
    ("direct", "broad", "#156477", "s", -0.08),
    ("all", "wcvp", "#b35b15", "o", 0.08),
    ("direct", "wcvp", "#b35b15", "s", 0.24),
]
fig, axes = plt.subplots(1, 4, figsize=(16, 7), sharex=True, sharey=True)
for ax, region, title in zip(
    axes,
    regions,
    ["Northern mid-latitudes", "Northern high latitudes", "Tropics", "Southern extratropics"],
):
    ax.axvline(0, color="#777777", lw=0.8)
    for low, high in [(-0.5, 2.5), (3.5, 6.5)]:
        ax.axhspan(low, high, color="#f2f5f5", zorder=0)
    for ev, flora, col, marker, offset in scopes:
        for y, trait in enumerate(traits):
            r = d.loc[
                (d.evidence_scope == ev)
                & (d.flora_scope == flora)
                & (d.context == region)
                & (d.outcome == trait)
            ].iloc[0]
            ax.plot([r.ci_low, r.ci_high], [y + offset, y + offset], color=col, lw=1, alpha=0.75)
            ax.plot(
                r.estimate,
                y + offset,
                marker,
                color=col,
                mfc=col if r.p_two_sided < 0.05 else "white",
                ms=5,
                zorder=3,
            )
    ax.set_title(title, fontsize=12, pad=13)
    ax.set_ylim(6.6, -0.6)
    ax.set_yticks(range(7), labels)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0, labelsize=11)
    ax.set_xlabel("Isolation slope (log odds)", fontsize=10)
fig.suptitle("Seven traits, four regions — no pooled syndrome score", fontsize=18, y=0.98)
legend = [
    Line2D([0], [0], marker=m, color=c, lw=1, label=f"{ev.title()} · {fl}")
    for ev, fl, c, m, _ in scopes
]
fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.56, 0.095), ncol=4, frameon=False)
fig.text(
    0.31,
    0.04,
    "Bars: pointwise 95% cluster-t intervals. Filled: two-sided unadjusted P < 0.05.\nIsolation is standardized within each trait × region × scope; WCVP indicates regional native compatibility.",
    fontsize=10,
)
fig.subplots_adjust(left=0.19, right=0.985, bottom=0.23, top=0.88, wspace=0.13)
for ext in ["png", "pdf", "svg"]:
    fig.savefig(p / f"traitwise_H1.{ext}", dpi=180)
