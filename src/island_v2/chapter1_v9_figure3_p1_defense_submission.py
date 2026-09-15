"""P1-defended submission renderer for Chapter 1 Figure 3.

The figure combines the frozen taxonomic-depth point estimates with two post-baseline
robustness tests: the matched-complexity pseudo-genus null (P1c) and paired spatial-
block uncertainty (P1d). No biological model is refit here.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import typer  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

from island_v2.chapter1_v8_figure3_assembly_depth import (  # noqa: E402
    _attenuation_panel,
    load_inputs,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


class Figure3P1Error(ValueError):
    """Raised when a frozen P1 input violates the publication contract."""


def _load_p1_inputs(p1c_root: Path, p1d_root: Path) -> dict[str, Any]:
    p1c_result_path = p1c_root / "chapter1_p1c_matched_genus_null_result.json"
    p1c_perm_path = p1c_root / "p1c_matched_genus_null_permutations.csv.gz"
    p1d_path = p1d_root / "p1b_paired_bootstrap_summary.csv"
    if not p1c_result_path.exists() or not p1c_perm_path.exists() or not p1d_path.exists():
        raise Figure3P1Error("missing frozen P1c/P1d input")

    p1c = json.loads(p1c_result_path.read_text(encoding="utf-8"))
    if p1c.get("status") != "true_genus_exceeds_matched_complexity_null":
        raise Figure3P1Error("P1c did not pass the frozen matched-complexity null")
    if int(p1c.get("valid_permutations", -1)) != 2000:
        raise Figure3P1Error("P1c must contain all 2000 valid permutations")
    if bool(p1c.get("permutations_regenerated", True)):
        raise Figure3P1Error("P1c aggregation must not regenerate permutations")

    permutations = pd.read_csv(p1c_perm_path)
    if len(permutations) != 2000 or permutations["permutation_id"].nunique() != 2000:
        raise Figure3P1Error("unexpected P1c permutation table")
    if not permutations["valid"].astype(bool).all():
        raise Figure3P1Error("invalid P1c permutation present")

    p1d = pd.read_csv(p1d_path)
    direct = p1d.loc[p1d["evidence_scope"].astype(str).eq("direct_only")].copy()
    if len(direct) != 8:
        raise Figure3P1Error("expected eight direct-only P1d profiles")
    if direct["extra_attenuation_ci_above_zero"].astype(bool).any():
        raise Figure3P1Error("P1d publication contract expects 0/8 extra-attenuation CIs above zero")
    return {"p1c": p1c, "permutations": permutations, "p1d": direct}


def _matched_null_panel(ax: plt.Axes, p1c: dict[str, Any], permutations: pd.DataFrame) -> None:
    values = pd.to_numeric(
        permutations["primary_median_conditional_attenuation"], errors="coerce"
    ).dropna()
    observed = float(p1c["observed_primary_median_conditional_attenuation"])
    ax.hist(values.to_numpy(float), bins=36, edgecolor="white", linewidth=0.4)
    ax.axvline(observed, linestyle="--", linewidth=1.6)
    ax.axvline(float(np.median(values)), linestyle=":", linewidth=1.2)
    ax.set_xlabel("Conditional genus attenuation")
    ax.set_ylabel("Matched pseudo-genus permutations")
    ax.set_title("B  True genus vs matched-complexity null", loc="left", fontweight="bold")
    ax.text(
        0.98,
        0.96,
        f"true genus = {observed:.3f}\nnull median = {np.median(values):.3f}\n57/2000 ≥ true\np = {float(p1c['one_sided_randomization_p_value']):.3f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.6,
    )


def _paired_uncertainty_panel(ax: plt.Axes, p1d: pd.DataFrame) -> None:
    order = p1d.sort_values(["stratum", "source_mode"]).reset_index(drop=True)
    y = np.arange(len(order))
    point = pd.to_numeric(order["observed_family_to_genus_extra_attenuation"]).to_numpy(float)
    low = pd.to_numeric(order["family_to_genus_extra_attenuation_ci_low"]).to_numpy(float)
    high = pd.to_numeric(order["family_to_genus_extra_attenuation_ci_high"]).to_numpy(float)
    xerr = np.vstack([point - low, high - point])
    ax.errorbar(point, y, xerr=xerr, fmt="o", capsize=2.5, linewidth=1.0)
    ax.axvline(0.0, linestyle="--", linewidth=1.0)
    labels = [
        f"{'native' if s == 'all_native' else 'NNE'} · {m}"
        for s, m in zip(order["stratum"].astype(str), order["source_mode"].astype(str))
    ]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.2)
    ax.invert_yaxis()
    ax.set_xlabel("Additional family→genus attenuation\n(point estimate and paired spatial-block 95% CI)")
    ax.set_title("C  Exact increment is spatially imprecise", loc="left", fontweight="bold")
    ax.text(
        0.98,
        0.04,
        "0/8 CIs exclude zero",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        fontweight="bold",
    )


def _decision_panel(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.set_title("D  What survives the cross-examination?", loc="left", fontweight="bold")
    rows = [
        ("Same observations and weights at all stages", "PASS"),
        ("True genus > matched arbitrary fine grouping", "PASS"),
        ("Total genus attenuation remains large", "PASS"),
        ("Exact family→genus increment precisely estimated", "NO"),
        ("Causal assembly mechanism identified", "NO"),
    ]
    y = 0.88
    for label, status in rows:
        ax.text(0.02, y, label, transform=ax.transAxes, fontsize=8.7, va="center")
        ax.text(0.97, y, status, transform=ax.transAxes, fontsize=9.2, fontweight="bold", ha="right", va="center")
        y -= 0.15
    ax.text(
        0.02,
        0.06,
        "Defensible claim:\nThe Palearctic floral-island response is strongly genus-structured,\nbut the size of the family→genus increment is not precisely localized.",
        transform=ax.transAxes,
        fontsize=8.6,
        va="bottom",
        fontweight="bold",
    )


def render_submission(
    *,
    effect_root: Path,
    pr142_root: Path,
    p1c_root: Path,
    p1d_root: Path,
    output_dir: Path,
    provenance: dict[str, Any],
) -> dict[str, Any]:
    baseline = load_inputs(effect_root, pr142_root)
    p1 = _load_p1_inputs(p1c_root, p1d_root)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(14.4, 9.5), constrained_layout=False)
    grid = GridSpec(2, 2, figure=fig, width_ratios=[1.22, 1.0], height_ratios=[1.0, 1.0],
                    hspace=0.38, wspace=0.30, top=0.88, bottom=0.10, left=0.08, right=0.97)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    _attenuation_panel(ax_a, baseline["attenuation"])
    ax_a.set_title("A  Frozen taxonomic attenuation profiles", loc="left", fontweight="bold")
    _matched_null_panel(ax_b, p1["p1c"], p1["permutations"])
    _paired_uncertainty_panel(ax_c, p1["p1d"])
    _decision_panel(ax_d)

    fig.suptitle(
        "The Palearctic floral-island response is specifically genus-structured",
        x=0.08,
        y=0.97,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.08,
        0.935,
        "True genus boundaries outperform matched arbitrary within-family partitions, while the exact family→genus increment remains spatially imprecise.",
        fontsize=9,
    )

    basename = "chapter1_v9_figure3_p1_defense_submission"
    paths = [output_dir / f"{basename}.{ext}" for ext in ("png", "svg", "pdf")]
    fig.savefig(paths[0], dpi=300, bbox_inches="tight")
    fig.savefig(paths[1], bbox_inches="tight")
    fig.savefig(paths[2], bbox_inches="tight")
    plt.close(fig)

    p1["p1d"].to_csv(output_dir / "figure3_p1d_direct_only_paired_uncertainty.csv", index=False)
    p1["permutations"].to_csv(output_dir / "figure3_p1c_matched_genus_null.csv.gz", index=False)
    baseline["attenuation"].to_csv(output_dir / "figure3_frozen_attenuation_profiles.csv", index=False)

    manifest = {
        "contract": "chapter1_v9_figure3_p1_defense_submission_v1",
        "status": "rendered_from_frozen_baseline_and_p1_robustness_results",
        "provenance": provenance,
        "p1c": {
            "valid_permutations": 2000,
            "observed_statistic": float(p1["p1c"]["observed_primary_median_conditional_attenuation"]),
            "one_sided_p": float(p1["p1c"]["one_sided_randomization_p_value"]),
        },
        "p1d": {
            "direct_only_profiles": 8,
            "extra_attenuation_ci_above_zero_profiles": 0,
        },
        "new_biological_models_fitted": False,
        "permutations_regenerated": False,
        "claim_boundary": (
            "Figure 3 supports genus-specific taxonomic structure beyond matched grouping complexity, "
            "but not a precisely estimated incremental family-to-genus attenuation or a causal assembly mechanism."
        ),
        "outputs": [p.name for p in paths],
    }
    (output_dir / "chapter1_v9_figure3_p1_defense_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("render")
def render_command(
    effect_root: Path = typer.Option(..., exists=True, file_okay=False),
    pr142_root: Path = typer.Option(..., exists=True, file_okay=False),
    p1c_root: Path = typer.Option(..., exists=True, file_okay=False),
    p1d_root: Path = typer.Option(..., exists=True, file_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    provenance = {
        "effect_fingerprint": {"run_id": 34796763611, "artifact_id": 10330230959, "digest": "sha256:a2db683b6d46ca09d4eae6d9a3bfcaf85a4f172b5589227958bc388df8edba60"},
        "primary_analysis": {"run_id": 34232450884, "artifact_id": 10058653212, "digest": "sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999"},
        "p1c": {"run_id": 34941827774, "artifact_id": 10385820775, "digest": "sha256:861d18fe87b9190619f925a2446be5fd4d460b818825578930883257c0a6ed16", "source_permutation_run_id": 34936193944},
        "p1d": {"run_id": 34939113182, "artifact_id": 10383754545, "digest": "sha256:822c69e9a1be4eb2c340391c1fb84dfba91b45b8a765549bb8ad1ec6379ef480"},
    }
    typer.echo(json.dumps(render_submission(effect_root=effect_root, pr142_root=pr142_root, p1c_root=p1c_root, p1d_root=p1d_root, output_dir=output_dir, provenance=provenance), indent=2))


if __name__ == "__main__":
    app()
