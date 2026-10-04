"""Rebuild poster H1-H4 from final repository inference and original plot observations."""

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import yaml
from scipy.special import expit
from scipy.stats import t

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from island_v2.chapter1_all_data_probability import _fit_single_beta_binomial, _standardize

parser = argparse.ArgumentParser()
parser.add_argument("--workspace", type=Path, required=True)
args = parser.parse_args()
repo = Path(__file__).resolve().parents[1]
w = args.workspace
s = w / "outputs/corrected_figures/source_data"
out = repo / "results/poster_q1_20261004"
out.mkdir(parents=True, exist_ok=True)
cfg = yaml.safe_load((repo / "config/chapter1_h1_final_traitwise_t_20261004.yml").read_text())
fin = pd.read_csv(repo / "results/h1_final_traitwise_t_20261004/all_broad.csv")
obs = pd.read_csv(s / "h1_observed_islands.csv.gz")
cs = cfg["contexts"]
traits = cfg["model_outcomes"]
cols = ["#3E6DA1", "#669347", "#CD8638", "#8E639C"]
labs = ["N. mid-latitude", "N. high-latitude", "Tropical", "S. extratropical"]
tlabs = [
    "Self-compatible",
    "Selfing mating system",
    "Autonomous selfing",
    "Plain colour",
    "Generalized form",
    "Radial symmetry",
    "Shallow / open tube",
]
plt.rcParams.update(
    {
        "font.family": "Arial",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "text.color": "#163849",
        "axes.labelcolor": "#163849",
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def save(fig, name):
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(out / f"{name}.{ext}", dpi=250)
    plt.close(fig)


rows = []
checks = []
for c in cs:
    for trait in traits:
        a = obs.loc[obs.analysis_regime.eq(c) & obs.outcome.eq(trait)]
        predictors = [*cfg["baseline_covariates"], cfg["geography_column"]]
        X = np.column_stack([np.ones(len(a)), *[_standardize(a[v]) for v in predictors]])
        fit = _fit_single_beta_binomial(
            a.successes.to_numpy(float),
            a.trials.to_numpy(float),
            X,
            ["intercept", *predictors],
            max_iter=5000,
        )
        assert fit["success"]
        labels = a.spatial_block.to_numpy(str)
        g = len(np.unique(labels))
        k = len(fit["theta"])
        scores = np.array([fit["score"][labels == b].sum(axis=0) for b in np.unique(labels)])
        V = (
            fit["bread"]
            @ scores.T
            @ scores
            @ fit["bread"].T
            * g
            / (g - 1)
            * (len(a) - 1)
            / (len(a) - k)
        )
        ref = fin.loc[fin.context.eq(c) & fin.outcome.eq(trait)].iloc[0]
        j = len(predictors)
        assert (
            abs(fit["theta"][j] - ref.estimate) < 1e-6 and abs(np.sqrt(V[j, j]) - ref.se) < 1e-5
        ), (
            c,
            trait,
            len(a),
            ref.n_islands,
            fit["theta"][j],
            ref.estimate,
            np.sqrt(V[j, j]),
            ref.se,
        )
        checks.append(
            {
                "context": c,
                "trait": trait,
                "slope_delta": float(fit["theta"][j] - ref.estimate),
                "se_delta": float(np.sqrt(V[j, j]) - ref.se),
            }
        )
        z = X[:, -1]
        grid = np.linspace(z.min(), z.max(), 81)
        Z = np.zeros((81, k))
        Z[:, 0] = 1
        Z[:, j] = grid
        eta = Z @ fit["theta"]
        err = np.sqrt(np.einsum("ij,jk,ik->i", Z, V, Z))
        crit = t.ppf(0.975, g - 1)
        rows.extend(
            {
                "context": c,
                "outcome": trait,
                "z_isolation": x,
                "prediction": y,
                "ci_low": lo,
                "ci_high": hi,
            }
            for x, y, lo, hi in zip(
                grid, expit(eta), expit(eta - crit * err), expit(eta + crit * err)
            )
        )
pred = pd.DataFrame(rows)
pred.to_csv(out / "h1_predictions.csv", index=False)
fig, axs = plt.subplots(2, 4, figsize=(14.2, 5.14), sharex=True, sharey=True)
fig.subplots_adjust(left=0.06, right=0.99, top=0.87, bottom=0.16, hspace=0.63, wspace=0.17)
for trait, title, ax in zip(traits, tlabs, axs.flat):
    for idx, (c, col) in enumerate(zip(cs, cols)):
        a = obs.loc[obs.analysis_regime.eq(c) & obs.outcome.eq(trait)]
        q = pred.loc[pred.context.eq(c) & pred.outcome.eq(trait)]
        ax.scatter(
            a.z_isolation,
            a.observed_proportion,
            s=4,
            alpha=0.095,
            color=col,
            edgecolors="none",
            rasterized=True,
        )
        ax.fill_between(q.z_isolation, q.ci_low, q.ci_high, color=col, alpha=0.16, lw=0)
        ax.plot(q.z_isolation, q.prediction, color=col, lw=1.8)
        p = fin.loc[fin.context.eq(c) & fin.outcome.eq(trait), "p_two_sided"].iloc[0]
        star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        ax.text(
            0.1 + idx * 0.265,
            1.025,
            star,
            color=col,
            transform=ax.transAxes,
            fontsize=15,
            ha="center",
            fontweight="bold",
        )
    ax.set_title(title, fontsize=17, fontweight="bold", pad=23)
    ax.set_ylim(-0.03, 1.03)
    ax.set_yticks([0, 0.5, 1])
    ax.set_xticks([-1, 0, 1, 2])
    ax.tick_params(labelsize=14)
axs.flat[-1].axis("off")
axs.flat[-1].legend(
    handles=[Line2D([], [], color=c, lw=2, label=l) for c, l in zip(cols, labs)],
    loc="center",
    frameon=False,
    fontsize=16,
)
fig.supxlabel("Isolation (SD within each region and trait)", fontsize=17, y=0.025)
fig.supylabel("Trait frequency", fontsize=17, x=0.007)
save(fig, "poster_h1")
audit = repo / "results/h1_final_directional_20261003"
h2 = pd.read_csv(audit / "h2_all_finite_cluster_audit.csv")
fig, axs = plt.subplots(1, 2, figsize=(14.2, 2.5), sharey=True)
fig.subplots_adjust(left=0.23, right=0.86, top=0.80, bottom=0.32, wspace=0.62)
for ax, response, col, title in zip(
    axs,
    ["generalized_accessible", "plain_colour"],
    ["#407847", "#7866A4"],
    ["STRUCTURE: floral access", "COLOUR: plainness"],
):
    ax.axvline(0, c="#687580", ls="--", lw=1)
    for i, c in enumerate(cs):
        r = h2.loc[h2.context.eq(c) & h2.response.eq(response)].iloc[0]
        q = r.primary_H2b_q_finite_cluster_t
        co = col if q < 0.05 else "#687580"
        ax.errorbar(
            r.distance_estimate,
            3 - i,
            xerr=t.ppf(0.975, r.n_clusters - 1) * r.distance_se,
            fmt="o",
            color=co,
            mfc=co if q < 0.05 else "white",
            ms=8,
            capsize=2.5,
        )
        ax.text(
            1.035,
            3 - i,
            "q" + ("<.001" if q < 0.001 else f"={q:.3f}"),
            transform=ax.get_yaxis_transform(),
            fontsize=18,
            va="center",
            color=co,
        )
    ax.set_ylim(-0.6, 3.6)
    ax.set_xlim(-0.05, 0.21)
    ax.set_xticks([0, 0.1, 0.2])
    ax.tick_params(labelsize=18)
    ax.set_title(title, fontsize=21, color=col, fontweight="bold")
axs[0].set_yticks(range(4), labs[::-1], fontsize=20)
fig.supxlabel("Isolation effect after selfing adjustment (95% t CI)", fontsize=19, y=0.02)
save(fig, "poster_h2")
h3 = pd.read_csv(audit / "h3_finite_publication_audit.csv").query("analysis=='primary'").iloc[0]
a = pd.read_csv(s / "h3_observed_measurement_cells.csv.gz")
pr = pd.read_csv(s / "h3_frozen_partial_prediction.csv")
fig, axs = plt.subplots(1, 2, figsize=(13.4, 3.1))
fig.subplots_adjust(left=0.085, right=0.985, bottom=0.27, top=0.79, wspace=0.34)
axs[0].scatter(
    a.z_log1p_distance_to_major_continent_km,
    a.adjusted_partial_response,
    s=10,
    alpha=0.2,
    c="#087E83",
    edgecolors="none",
    rasterized=True,
)
axs[0].plot(pr.z_isolation, pr.prediction, c="#087E83")
axs[0].set_title("All 1,408 adjusted cells", fontsize=20)
axs[0].set_ylabel("Adjusted pollen limitation", fontsize=18)
x = pr.z_isolation.to_numpy()
dx = x - x.min()
b = h3.estimate
ci = t.ppf(0.975, h3.n_publications - 1) * h3.se
axs[1].fill_between(x, (b - ci) * dx, (b + ci) * dx, color="#087E83", alpha=0.22)
axs[1].plot(x, b * dx, c="#087E83", lw=2.6)
axs[1].set_title(
    f"Model contrast + 95% t CI; p={h3.finite_publication_t_two_sided_p:.4f}", fontsize=18
)
axs[1].set_ylabel("Predicted difference\nfrom minimum isolation", fontsize=17)
for ax in axs:
    ax.set_xlabel("Isolation (SD of log distance)", fontsize=18)
    ax.tick_params(labelsize=16)
save(fig, "poster_h3")
h4 = pd.read_csv(audit / "h4_finite_publication_audit.csv").query("analysis=='primary'")
fig, axs = plt.subplots(1, 2, figsize=(14.2, 2.8))
fig.subplots_adjust(left=0.10, right=0.99, bottom=0.25, top=0.77, wspace=0.25)
for ax, (_, r), col, label in zip(
    axs, h4.iterrows(), ["#087E83", "#407847"], ["Selfing", "Floral access"]
):
    x = np.linspace(0, 1, 101)
    ci = t.ppf(0.975, r.n_publications - 1) * r.se
    ax.axhline(0, c="#687580", ls="--", lw=0.7)
    ax.fill_between(x, (r.estimate - ci) * x, (r.estimate + ci) * x, color=col, alpha=0.23)
    ax.plot(x, r.estimate * x, c=col, lw=2.5)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.61, 0.055)
    ax.set_xticks([0, 0.5, 1])
    ax.set_yticks([-0.6, -0.3, 0])
    ax.tick_params(labelsize=17)
    ax.set_xlabel(label + " score (0–1)", fontsize=19)
    ax.set_title(
        f"{label}: {int(r.n_species)} species, {int(r.n_cells)} cells\nβ={r.estimate:.3f}; p={r.finite_publication_t_two_sided_p:.4f}",
        fontsize=18,
        fontweight="bold",
    )
axs[0].set_ylabel("Predicted difference\nfrom score 0 (95% t CI)", fontsize=18)
save(fig, "poster_h4")
(out / "curve_verification.json").write_text(json.dumps(checks, indent=2), encoding="utf-8")
print(
    "All 28 H1 curves reproduce final slopes within 1e-6 and SE within 1e-5 (numerical Hessian tolerance). H2-H4 use current audit tables."
)
