"""Submission Figure 2 for the P2-defended Chapter 1 surface.

This figure uses only the frozen P2 artifact.  It shows the formal same-layer
northern-midlatitude versus tropical H2 contrast and keeps the separate
Palearctic--Neotropical realm comparison as a claim-boundary panel.
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

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _read(root: Path, name: str) -> pd.DataFrame:
    path = root / name
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def _profile_label(scope: str, stratum: str) -> str:
    scope_label = "All" if scope == "all_analysis_eligible" else "Direct"
    stratum_label = "native" if stratum == "all_native" else "NNE"
    return f"{scope_label} · {stratum_label}"


def _draw_vectors(ax: plt.Axes, frame: pd.DataFrame, *, title: str) -> None:
    ax.axhline(0, linewidth=0.8, color="0.75")
    ax.axvline(0, linewidth=0.8, color="0.75")
    for idx, row in frame.reset_index(drop=True).iterrows():
        north = np.array([
            float(row["northern_midlatitude_accessibility_generalization_slope_joint"]),
            float(row["northern_midlatitude_reproductive_assurance_slope_joint"]),
        ])
        trop = np.array([
            float(row["tropical_accessibility_generalization_slope_joint"]),
            float(row["tropical_reproductive_assurance_slope_joint"]),
        ])
        label = _profile_label(str(row["evidence_scope"]), str(row["stratum"]))
        line = ax.plot([north[0], trop[0]], [north[1], trop[1]], marker="o", linewidth=1.4, alpha=0.8, label=label)[0]
        ax.scatter([north[0]], [north[1]], marker="o", s=55, color=line.get_color())
        ax.scatter([trop[0]], [trop[1]], marker="^", s=65, color=line.get_color())
        ax.annotate("N", tuple(north), xytext=(4, 3), textcoords="offset points", fontsize=7)
        ax.annotate("T", tuple(trop), xytext=(4, 3), textcoords="offset points", fontsize=7)
    ax.set_xlabel("Accessibility/generalization slope")
    ax.set_ylabel("Reproductive-assurance slope")
    ax.set_title(title, loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7, loc="best")


def render(
    *,
    input_dir: Path,
    output_dir: Path,
    source_run_id: int,
    source_artifact_id: int,
    source_artifact_digest: str,
) -> dict[str, Any]:
    p2a = _read(input_dir, "p2a_common_island_vector_results.csv")
    boot = _read(input_dir, "p2a_paired_block_bootstrap_summary.csv")
    p2b = _read(input_dir, "p2b_coobserved_species_vector_results.csv")
    realm = _read(input_dir, "p2c_realm_boundary_audit.csv")
    result = json.loads((input_dir / "chapter1_p2_component_nonconcordance_result.json").read_text())

    output_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.3))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    _draw_vectors(ax_a, p2a, title="A  Same-island North–Tropical response vectors")

    primary = p2a.loc[
        p2a["evidence_scope"].eq("direct_only")
        & p2a["stratum"].eq("native_nonendemic")
    ].iloc[0]
    bsum = boot.loc[
        boot["evidence_scope"].eq("direct_only")
        & boot["stratum"].eq("native_nonendemic")
    ].iloc[0]
    ax_b.axvline(0, linewidth=0.8, color="0.75")
    ax_b.errorbar(
        [float(primary["determinant"])],
        [0],
        xerr=[[
            float(primary["determinant"]) - float(bsum["determinant_ci_low"])
        ], [
            float(bsum["determinant_ci_high"]) - float(primary["determinant"])
        ]],
        fmt="o",
        capsize=4,
    )
    ax_b.set_yticks([0])
    ax_b.set_yticklabels(["det(N,T)"])
    ax_b.set_xlabel("2D determinant (95% paired-block interval)")
    ax_b.set_title("B  Joint difference survives; strong non-collinearity does not", loc="left", fontweight="bold")
    ax_b.text(
        0.02,
        0.82,
        f"joint vector p = {float(primary['p_value']):.4g}\n"
        f"angle = {float(primary['angle_degrees']):.1f}°\n"
        f"angle 95% = [{float(bsum['angle_ci_low']):.1f}, {float(bsum['angle_ci_high']):.1f}]°",
        transform=ax_b.transAxes,
        va="top",
        fontsize=9,
    )

    p2b_primary = p2b.loc[
        p2b["evidence_scope"].eq("direct_only")
        & p2b["stratum"].eq("native_nonendemic")
    ].copy()
    _draw_vectors(ax_c, p2b_primary, title="C  Same co-observed species denominator")
    row = p2b_primary.iloc[0]
    ax_c.text(
        0.02,
        0.98,
        f"853 co-observed species\n280 NNE islands\njoint vector p = {float(row['p_value']):.4g}",
        transform=ax_c.transAxes,
        va="top",
        fontsize=8.5,
    )

    ax_d.axvline(0.05, linewidth=1.0, linestyle="--", color="0.55")
    realm = realm.copy()
    realm["label"] = [
        _profile_label(str(s), str(t)) for s, t in zip(realm["evidence_scope"], realm["stratum"], strict=True)
    ]
    y = np.arange(len(realm))
    ax_d.scatter(realm["p_value"].astype(float), y, s=55)
    ax_d.set_yticks(y)
    ax_d.set_yticklabels(realm["label"], fontsize=8)
    ax_d.set_xlim(0, max(0.45, float(realm["p_value"].max()) + 0.03))
    ax_d.set_xlabel("Frozen Palearctic–Neotropical direct-test p")
    ax_d.set_title("D  Realm contrast remains a claim boundary", loc="left", fontweight="bold")
    ax_d.text(
        0.02,
        0.08,
        "No realm profile passes p < 0.05.\nDo not relabel Palearctic vs tropical\nas the formal H2 direct contrast.",
        transform=ax_d.transAxes,
        fontsize=8.5,
        va="bottom",
    )

    fig.suptitle(
        "The formal North–Tropical floral/reproductive response difference survives common support",
        x=0.06,
        y=0.98,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.06,
        0.945,
        "P2 post-baseline audit: same islands, paired spatial blocks, and a stricter same-species sensitivity; no new trait definitions.",
        fontsize=9,
        color="0.35",
    )
    fig.tight_layout(rect=[0.04, 0.04, 0.99, 0.92])

    base = "chapter1_v10_figure2_p2_submission"
    png = output_dir / f"{base}.png"
    svg = output_dir / f"{base}.svg"
    pdf = output_dir / f"{base}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    manifest = {
        "contract": "chapter1_v10_figure2_p2_submission_v1",
        "source_p2_run_id": int(source_run_id),
        "source_p2_artifact_id": int(source_artifact_id),
        "source_p2_artifact_digest": str(source_artifact_digest),
        "p2_status": result["status"],
        "primary_profile": result["primary_profile"],
        "primary_joint_vector_p": result["p2a_primary_p_value"],
        "primary_determinant": result["p2a_primary_determinant"],
        "primary_determinant_ci": result["p2a_primary_determinant_ci"],
        "common_species_executable": result["p2b_common_species_sensitivity_executable"],
        "palearctic_neotropical_used_as_primary_direct_contrast": False,
        "new_biological_models_fitted": False,
        "outputs": [png.name, svg.name, pdf.name],
    }
    (output_dir / "chapter1_v10_figure2_p2_submission_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command()
def main(
    input_dir: Path = typer.Option(..., exists=True, file_okay=False),
    output_dir: Path = typer.Option(...),
    source_run_id: int = typer.Option(...),
    source_artifact_id: int = typer.Option(...),
    source_artifact_digest: str = typer.Option(...),
) -> None:
    typer.echo(json.dumps(render(
        input_dir=input_dir,
        output_dir=output_dir,
        source_run_id=source_run_id,
        source_artifact_id=source_artifact_id,
        source_artifact_digest=source_artifact_digest,
    ), indent=2))


if __name__ == "__main__":
    app()
