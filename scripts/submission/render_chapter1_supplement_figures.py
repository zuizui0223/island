from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REGION_LABELS = {
    "northern_midlatitude": "Northern mid-latitude",
    "northern_high_latitude": "Northern high latitude",
    "tropical": "Tropical",
    "southern_extratropical": "Southern extratropical",
}
OUTCOME_LABELS = {
    "self_compatibility": "Self-compatibility",
    "selfing_mating_system": "Selfing mating system",
    "autonomous_selfing": "Autonomous selfing",
    "plain_colour": "Plain colour",
    "generalized_form": "Generalized form",
    "actinomorphic_symmetry": "Actinomorphic symmetry",
    "shallow_open_tube": "Shallow/open tube",
}
COLOUR_LABELS = {
    "white": "White",
    "red_pink": "Red/pink",
    "yellow_orange": "Yellow/orange",
    "blue_purple": "Blue/purple",
    "green_brown_inconspicuous": "Green/brown/inconspicuous",
}


def _style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "svg.hashsalt": "chapter1-supplement-20260927",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def _save(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(
        output_dir / f"{stem}.pdf",
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


def figure_s1(root: Path, output_dir: Path) -> None:
    frame = pd.read_csv(root / "results/geography_20260924/formerly_zero_islands_recalculated.csv")
    values = frame["spherical_coast_distance_km"].astype(float)
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    bins = np.logspace(
        np.log10(max(values.min(), 1e-4)),
        np.log10(values.max()),
        35,
    )
    ax.hist(values, bins=bins, edgecolor="black", linewidth=0.3)
    ax.set_xscale("log")
    ax.set_xlabel("Corrected coastline distance (km; log scale)")
    ax.set_ylabel("Number of formerly zero-distance islands")
    ax.set_title("Figure S1 | Corrected distances for 1,113 formerly zero islands")
    median = float(values.median())
    ax.axvline(median, linestyle="--", linewidth=1.2)
    ax.text(
        median,
        ax.get_ylim()[1] * 0.92,
        f"median = {median:.3f} km",
        rotation=90,
        va="top",
        ha="right",
    )
    ax.text(
        0.01,
        0.98,
        "Previous exposure: 0 km for every island shown",
        transform=ax.transAxes,
        va="top",
    )
    _save(fig, output_dir, "Figure_S1_corrected_distance_distribution")


def figure_s2(root: Path, output_dir: Path) -> None:
    frame = pd.read_csv(
        root / "submission/chapter1_current/supplement/Table_S2a_H1_atomic.csv"
    )
    frame = frame.loc[
        frame["evidence_scope"].eq("direct_only")
        & frame["stratum"].eq("all_observed")
    ].copy()
    outcomes = list(OUTCOME_LABELS)
    contexts = list(REGION_LABELS)
    fig, axes = plt.subplots(1, 4, figsize=(12.5, 5.2), sharey=True)
    y = np.arange(len(outcomes))
    for ax, context in zip(axes, contexts, strict=True):
        part = frame.loc[frame["context"].eq(context)].set_index("outcome")
        est = np.array([float(part.loc[x, "geography_slope_log_odds"]) for x in outcomes])
        se = np.array([float(part.loc[x, "cluster_robust_se"]) for x in outcomes])
        ax.errorbar(est, y, xerr=1.96 * se, fmt="o", capsize=2)
        ax.axvline(0, linewidth=0.8)
        ax.set_title(REGION_LABELS[context])
        ax.set_xlabel("Isolation coefficient")
        ax.grid(axis="x", linewidth=0.3, alpha=0.4)
    axes[0].set_yticks(y, [OUTCOME_LABELS[x] for x in outcomes])
    axes[0].invert_yaxis()
    fig.suptitle("Figure S2 | Direct-only H1 atomic responses by region", y=1.02)
    fig.tight_layout()
    _save(fig, output_dir, "Figure_S2_H1_direct_forest")


def _matrix_plot(
    matrix: pd.DataFrame,
    q_matrix: pd.DataFrame,
    *,
    title: str,
    xlabel: str,
    ylabel: str,
    output_dir: Path,
    stem: str,
    row_labels: list[str] | None = None,
) -> None:
    arr = matrix.to_numpy(dtype=float)
    vmax = float(np.nanmax(np.abs(arr))) if np.isfinite(arr).any() else 1.0
    vmax = max(vmax, 1e-6)
    fig_height = max(4.0, 0.3 * matrix.shape[0] + 1.8)
    fig, ax = plt.subplots(figsize=(8.5, fig_height))
    im = ax.imshow(arr, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(matrix.shape[1]), matrix.columns, rotation=30, ha="right")
    labels = row_labels if row_labels is not None else list(matrix.index)
    ax.set_yticks(np.arange(matrix.shape[0]), labels)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = arr[i, j]
            if not np.isfinite(value):
                continue
            q = q_matrix.iloc[i, j]
            star = ""
            if pd.notna(q):
                qf = float(q)
                if qf < 0.001:
                    star = "***"
                elif qf < 0.01:
                    star = "**"
                elif qf < 0.05:
                    star = "*"
            ax.text(j, i, f"{value:.2f}{star}", ha="center", va="center", fontsize=7)
    cbar = fig.colorbar(im, ax=ax, shrink=0.75)
    cbar.set_label("Isolation coefficient")
    fig.tight_layout()
    _save(fig, output_dir, stem)


def figure_s3(root: Path, output_dir: Path) -> None:
    frame = pd.read_csv(
        root / "results/geography_20260924/all/raw_patterns/raw_colour_model_results.csv"
    )
    part = frame.loc[
        frame["stratum"].eq("all_observed")
        & frame["support_tier"].eq("confirmatory")
        & frame["model"].eq("conditional_selfing")
        & frame["status"].eq("fit")
    ].copy()
    contexts = list(REGION_LABELS)
    colours = [x for x in COLOUR_LABELS if x in set(part["colour"])]
    matrix = (
        part.pivot(index="colour", columns="context", values="distance_estimate")
        .reindex(index=colours, columns=contexts)
        .rename(columns=REGION_LABELS)
    )
    q_matrix = (
        part.pivot(index="colour", columns="context", values="distance_q")
        .reindex(index=colours, columns=contexts)
        .rename(columns=REGION_LABELS)
    )
    _matrix_plot(
        matrix,
        q_matrix,
        title="Figure S3 | Selfing-adjusted raw colour responses",
        xlabel="Region",
        ylabel="Colour state",
        output_dir=output_dir,
        stem="Figure_S3_H2_raw_colour_matrix",
        row_labels=[COLOUR_LABELS[x] for x in colours],
    )


def _compact_combination_label(value: str) -> str:
    label = value.replace("__", " | ").replace("_given_colour", "")
    label = label.replace("_", " ")
    return label


def figure_s4(root: Path, output_dir: Path) -> None:
    frame = pd.read_csv(
        root
        / "results/geography_20260924/all/raw_patterns/"
        "raw_colour_conditioned_architecture_model_results.csv"
    )
    part = frame.loc[
        frame["stratum"].eq("all_observed")
        & frame["support_tier"].eq("confirmatory")
        & frame["model"].eq("conditional_selfing")
        & frame["status"].eq("fit")
    ].copy()
    contexts = list(REGION_LABELS)
    support = (
        part.groupby("combination", as_index=False)
        .agg(
            min_q=("distance_q", "min"),
            max_abs=("distance_estimate", lambda x: float(np.nanmax(np.abs(x)))),
        )
        .sort_values(["min_q", "max_abs"], ascending=[True, False])
    )
    keep = support["combination"].head(24).tolist()
    matrix = (
        part.loc[part["combination"].isin(keep)]
        .pivot(index="combination", columns="context", values="distance_estimate")
        .reindex(index=keep, columns=contexts)
        .rename(columns=REGION_LABELS)
    )
    q_matrix = (
        part.loc[part["combination"].isin(keep)]
        .pivot(index="combination", columns="context", values="distance_q")
        .reindex(index=keep, columns=contexts)
        .rename(columns=REGION_LABELS)
    )
    _matrix_plot(
        matrix,
        q_matrix,
        title="Figure S4 | Colour-conditioned architecture responses (top 24 by q)",
        xlabel="Region",
        ylabel="Colour | architecture combination",
        output_dir=output_dir,
        stem="Figure_S4_H2_colour_conditioned_architecture",
        row_labels=[_compact_combination_label(x) for x in keep],
    )


def figure_s5(root: Path, output_dir: Path) -> None:
    frame = pd.read_csv(
        root / "submission/chapter1_current/supplement/Table_S5_H3_pollen_limitation.csv"
    )
    order = ["primary", "no_zero_constant", "supplemental_only"]
    part = frame.set_index("analysis").loc[order].reset_index()
    y = np.arange(len(part))
    est = part["estimate"].astype(float).to_numpy()
    se = part["se"].astype(float).to_numpy()
    labels = ["Primary", "No-zero constant", "Supplemental only"]
    fig, ax = plt.subplots(figsize=(6.6, 3.3))
    ax.errorbar(est, y, xerr=1.96 * se, fmt="o", capsize=3)
    ax.axvline(0, linewidth=0.8)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("Isolation coefficient for pollen limitation")
    ax.set_title("Figure S5 | H3 sensitivity estimates")
    ax.grid(axis="x", linewidth=0.3, alpha=0.4)
    fig.tight_layout()
    _save(fig, output_dir, "Figure_S5_H3_sensitivity_forest")


def figure_s6(root: Path, output_dir: Path) -> None:
    frame = pd.read_csv(
        root / "submission/chapter1_current/supplement/Table_S6b_H4_atomic.csv"
    )
    part = frame.loc[frame["analysis"].eq("primary") & frame["evaluable"].astype(str).eq("True")].copy()
    part["estimate"] = pd.to_numeric(part["trait_state_estimate"], errors="coerce")
    part["se"] = pd.to_numeric(part["trait_state_se"], errors="coerce")
    part = part.dropna(subset=["estimate", "se"])
    labels = [OUTCOME_LABELS.get(x, x.replace("_", " ").title()) for x in part["trait"]]
    y = np.arange(len(part))
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    ax.errorbar(
        part["estimate"].to_numpy(float),
        y,
        xerr=1.96 * part["se"].to_numpy(float),
        fmt="o",
        capsize=3,
    )
    ax.axvline(0, linewidth=0.8)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("Trait-state coefficient on current pollen limitation")
    ax.set_title("Figure S6 | H4 primary atomic-trait sensitivities")
    ax.grid(axis="x", linewidth=0.3, alpha=0.4)
    fig.tight_layout()
    _save(fig, output_dir, "Figure_S6_H4_atomic_forest")


def figure_s7(root: Path, output_dir: Path) -> None:
    frame = pd.read_csv(
        root / "submission/chapter1_current/supplement/Table_S2a_H1_atomic.csv"
    )
    part = frame.loc[frame["stratum"].eq("all_observed")].copy()
    all_rows = part.loc[part["evidence_scope"].eq("all_analysis_eligible")].copy()
    direct_rows = part.loc[part["evidence_scope"].eq("direct_only")].copy()
    merged = all_rows.merge(
        direct_rows,
        on=["context", "outcome"],
        suffixes=("_all", "_direct"),
        validate="one_to_one",
    )
    merged["support_ratio"] = (
        merged["n_islands_direct"].astype(float) / merged["n_islands_all"].astype(float)
    )
    contexts = list(REGION_LABELS)
    outcomes = list(OUTCOME_LABELS)
    matrix = (
        merged.pivot(index="outcome", columns="context", values="support_ratio")
        .reindex(index=outcomes, columns=contexts)
        .rename(columns=REGION_LABELS)
    )
    arr = matrix.to_numpy(float)
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    im = ax.imshow(arr, aspect="auto", cmap="viridis", vmin=0.0, vmax=1.0)
    ax.set_xticks(np.arange(matrix.shape[1]), matrix.columns, rotation=30, ha="right")
    ax.set_yticks(np.arange(matrix.shape[0]), [OUTCOME_LABELS[x] for x in outcomes])
    ax.set_title("Figure S7 | Direct-only island support relative to all-analysis scope")
    ax.set_xlabel("Region")
    ax.set_ylabel("H1 atomic response")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = arr[i, j]
            if np.isfinite(value):
                ax.text(j, i, f"{100 * value:.0f}%", ha="center", va="center", fontsize=7)
    cbar = fig.colorbar(im, ax=ax, shrink=0.75)
    cbar.set_label("Direct-only / all-analysis island support")
    fig.tight_layout()
    _save(fig, output_dir, "Figure_S7_evidence_scope_support")


def build_manifest(output_dir: Path) -> None:
    records = []
    for path in sorted(output_dir.iterdir()):
        if path.suffix.lower() not in {".svg", ".pdf"}:
            continue
        records.append(
            {
                "file": path.name,
                "bytes": path.stat().st_size,
            }
        )
    (output_dir / "FIGURE_MANIFEST.json").write_text(
        json.dumps(
            {
                "contract": "chapter1_supplement_figures_v1",
                "n_figures": 7,
                "formats": ["svg", "pdf"],
                "files": records,
                "source_surface": "corrected_geography_20260924",
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
        default=Path("submission/chapter1_current/supplement/figures"),
    )
    args = parser.parse_args()
    root = args.root.resolve()
    output_dir = (root / args.output_dir).resolve() if not args.output_dir.is_absolute() else args.output_dir

    _style()
    figure_s1(root, output_dir)
    figure_s2(root, output_dir)
    figure_s3(root, output_dir)
    figure_s4(root, output_dir)
    figure_s5(root, output_dir)
    figure_s6(root, output_dir)
    figure_s7(root, output_dir)
    build_manifest(output_dir)


if __name__ == "__main__":
    main()
