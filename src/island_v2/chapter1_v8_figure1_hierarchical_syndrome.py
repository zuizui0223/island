"""Render Chapter 1 v8 Figure 1 as an inference map from frozen evidence locks.

Presentation only: no biological model is refit and no adverse or non-identified result
is converted into mechanism support. Quantitative annotations come from frozen v8 locks.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import typer  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

app = typer.Typer(add_completion=False, no_args_is_help=True)

DEFAULT_FIGURE3_LOCK = Path("config/chapter1_v8_figure3_submission_result_lock.json")
DEFAULT_FIGURE4_LOCK = Path("config/chapter1_v8_figure4_result_lock.json")
DEFAULT_N1_LOCK = Path("config/chapter1_nee_n1_result_lock.json")
DEFAULT_SUBMISSION_FREEZE = Path("docs/chapter1_submission_freeze_20260914.md")
DEFAULT_MAXIMAL_INTEGRATION = Path("docs/chapter1_v8_maximal_data_integration_20260914.md")


class FigureInputError(ValueError):
    """Raised when a frozen Figure 1 source no longer matches its contract."""


def _load_json(path: Path, contract: str) -> dict[str, Any]:
    if not path.is_file():
        raise FigureInputError(f"missing frozen Figure 1 source: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("contract") != contract:
        raise FigureInputError(
            f"unexpected contract for {path}: {payload.get('contract')!r} != {contract!r}"
        )
    return payload


def load_inputs(
    *,
    figure3_lock: Path = DEFAULT_FIGURE3_LOCK,
    figure4_lock: Path = DEFAULT_FIGURE4_LOCK,
    n1_lock: Path = DEFAULT_N1_LOCK,
    submission_freeze: Path = DEFAULT_SUBMISSION_FREEZE,
    maximal_integration: Path = DEFAULT_MAXIMAL_INTEGRATION,
) -> dict[str, Any]:
    """Read and fail-closed validate the frozen quantities shown in Figure 1."""
    fig3 = _load_json(figure3_lock, "chapter1_v8_figure3_submission_result_lock_v1")
    fig4 = _load_json(figure4_lock, "chapter1_v8_figure4_result_lock_v1")
    n1 = _load_json(n1_lock, "chapter1_nee_n1_result_lock_v1")
    if not submission_freeze.is_file() or not maximal_integration.is_file():
        raise FigureInputError("canonical v8 narrative source is missing")

    freeze_text = submission_freeze.read_text(encoding="utf-8")
    maximal_text = maximal_integration.read_text(encoding="utf-8")
    f3 = fig3.get("frozen_results", {})
    f4 = fig4.get("frozen_results", {})
    gate = n1.get("n1", {}).get("gate", {})

    if f3.get("support_ladder") != "4/4 -> 4/4 -> 0/4":
        raise FigureInputError("Figure 3 support ladder changed")
    genus = list(f3.get("genus_attenuation_fraction_range", []))
    conditional = list(f3.get("conditional_family_to_genus_attenuation_range", []))
    if genus != [0.788, 0.859] or conditional != [0.706, 0.791]:
        raise FigureInputError("Figure 3 attenuation range changed")

    expected_f4 = {
        "v6_palearctic_survival": "99/100",
        "v6_tropical_survival": "35/75",
        "geometry_promoted": "0/12",
        "h5d_qualified": "0/8",
    }
    for key, expected in expected_f4.items():
        if f4.get(key) != expected:
            raise FigureInputError(f"Figure 4 frozen result changed: {key}")

    h5c_p = float(f4.get("h5c_p_value"))
    h5c_est = float(f4.get("h5c_interaction_estimate"))
    if abs(h5c_p - 0.4122068222871858) > 1e-12:
        raise FigureInputError("H5c p-value changed")
    if gate.get("N1_pass") is not False:
        raise FigureInputError("Figure 1 requires the frozen failed N1 gate")
    n1_p = float(gate.get("global_Wald_p_value"))
    if abs(n1_p - 0.6551581817640355) > 1e-12:
        raise FigureInputError("N1 global p-value changed")
    if gate.get("failure_action") != "stop_before_N2_and_keep_frozen_Chapter1":
        raise FigureInputError("N1 stop rule changed")

    required = [
        "H4 promotion is 0/16",
        "The broad Palearctic primary response survives the finite predeclared MNAR trait-resolution grid",
        "Chapter 2 / `izu-core` is reserved for **how and why functionally**",
    ]
    for phrase in required:
        if phrase not in freeze_text:
            raise FigureInputError(f"canonical submission freeze no longer contains: {phrase}")
    if "- GloBI breadth: 0/4 promoted." not in maximal_text:
        raise FigureInputError("canonical GloBI breadth boundary changed")

    return {
        "support_ladder": f3["support_ladder"],
        "genus_attenuation_pct": [100 * float(x) for x in genus],
        "family_to_genus_attenuation_pct": [100 * float(x) for x in conditional],
        "v6_palearctic_survival": f4["v6_palearctic_survival"],
        "v6_tropical_survival": f4["v6_tropical_survival"],
        "geometry_promoted": f4["geometry_promoted"],
        "h5c_estimate": h5c_est,
        "h5c_p_value": h5c_p,
        "h5d_qualified": f4["h5d_qualified"],
        "h4_promoted": "0/16",
        "globi_breadth_promoted": "0/4",
        "n1_p_value": n1_p,
        "n1_pass": False,
        "n1_failure_action": gate["failure_action"],
        "source_locks": {
            "figure3": {
                "workflow_run_id": fig3["workflow_run_id"],
                "artifact_id": fig3["artifact_id"],
                "artifact_digest": fig3["artifact_digest"],
            },
            "figure4": {
                "workflow_run_id": fig4["workflow_run_id"],
                "artifact_id": fig4["artifact_id"],
                "artifact_digest": fig4["artifact_digest"],
            },
            "n1": {
                "workflow_run_id": n1["n1"]["workflow_run_id"],
                "artifact_id": n1["n1"]["artifact_id"],
                "artifact_digest": n1["n1"]["artifact_digest"],
            },
        },
    }


def _box(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    text: str,
    *,
    facecolor: str = "white",
    edgecolor: str = "0.35",
    fontsize: float = 8.1,
    weight: str = "normal",
    linewidth: float = 1.0,
) -> None:
    x, y = xy
    ax.add_patch(
        FancyBboxPatch(
            (x, y), width, height,
            boxstyle="round,pad=0.012,rounding_size=0.012",
            facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth,
        )
    )
    ax.text(
        x + width / 2, y + height / 2, text,
        ha="center", va="center", fontsize=fontsize, fontweight=weight,
    )


def _arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    dashed: bool = False,
    color: str = "0.3",
    linewidth: float = 1.2,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start, end, arrowstyle="-|>", mutation_scale=10,
            linewidth=linewidth, linestyle="--" if dashed else "-",
            color=color, shrinkA=2, shrinkB=2,
        )
    )


def _setup(ax: plt.Axes, letter: str, title: str) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(
        0.0, 1.02, f"{letter}  {title}", transform=ax.transAxes,
        ha="left", va="bottom", fontsize=11.2, fontweight="bold",
    )


def _panel_a(ax: plt.Axes) -> None:
    _setup(ax, "A", "Why can an island floral syndrome appear?")
    _box(ax, (0.03, 0.72), 0.30, 0.18, "SOURCE POOL\nA   B   C   D   E", facecolor="#f5f5f5", weight="bold", fontsize=9)
    _box(ax, (0.67, 0.72), 0.30, 0.18, "REMOTE ISLAND\nA       C       E", facecolor="#f5f5f5", weight="bold", fontsize=9)
    _arrow(ax, (0.34, 0.81), (0.66, 0.81))
    ax.text(0.50, 0.84, "geographic isolation", ha="center", va="bottom", fontsize=8)

    cards = [
        (0.05, "Repeated within-lineage response\n\nThe same lineage changes phenotype\nafter colonization."),
        (0.36, "Hierarchical lineage sorting\n\nLineages differ in persistence,\nchanging assemblage trait composition."),
        (0.67, "Mixed generator\n\nMembership change and within-lineage\nresponse can coexist."),
    ]
    for x, text in cards:
        _box(ax, (x, 0.27), 0.28, 0.28, text, facecolor="#fbfbfb", fontsize=7.8)
    ax.text(
        0.5, 0.11,
        "An assemblage mean alone cannot distinguish adaptation from membership change.",
        ha="center", fontsize=8.8, fontweight="bold",
    )
    ax.text(
        0.5, 0.04,
        "Pollinator icons are deliberately absent: Chapter 1 did not identify that causal link.",
        ha="center", fontsize=7.6, color="0.35",
    )


def _panel_b(ax: plt.Axes, values: dict[str, Any]) -> None:
    _setup(ax, "B", "The response branches, then collapses at an assembly depth")
    _box(ax, (0.05, 0.80), 0.90, 0.12, "ONE UNIVERSAL RESPONSE?  →  NO", facecolor="#f7f7f7", weight="bold", fontsize=9.2)
    _arrow(ax, (0.50, 0.79), (0.50, 0.69))
    _box(
        ax, (0.05, 0.53), 0.90, 0.15,
        "BIOGEOGRAPHIC BRANCHING\nPalearctic: accessibility/generalization ↑  +  assurance ↑\nTropical: accessibility/generalization ↓  +  assurance ↑",
        facecolor="#f7f7f7", fontsize=8.4, weight="bold",
    )
    _arrow(ax, (0.50, 0.52), (0.50, 0.43))
    ladder = values["support_ladder"].replace(" -> ", "  →  ")
    g0, g1 = values["genus_attenuation_pct"]
    c0, c1 = values["family_to_genus_attenuation_pct"]
    _box(
        ax, (0.05, 0.16), 0.90, 0.26,
        (
            "WHERE DOES THE PALEARCTIC RESPONSE LIVE?\n\n"
            f"observed        family-adjusted        genus-adjusted\n{ladder}\n\n"
            f"genus attenuation: {g0:.1f}–{g1:.1f}% of observed vector\n"
            f"family → genus conditional attenuation: {c0:.1f}–{c1:.1f}%"
        ),
        facecolor="#eef3f8", edgecolor="#496b8a", fontsize=8.5,
        weight="bold", linewidth=1.2,
    )
    ax.text(
        0.50, 0.08, "ASSEMBLY DEPTH: strongest attenuation at FAMILY → GENUS",
        ha="center", fontsize=9.2, fontweight="bold", color="#294b69",
    )


def _panel_c(ax: plt.Axes, values: dict[str, Any]) -> None:
    _setup(ax, "C", "Cross-examination: what survives, what is not promoted?")
    cards = [
        ("Trait missingness\n(V5)", "Palearctic core survives\nfinite MNAR grid", "#edf4ed"),
        ("Species-list bias\n(V6)", f"Palearctic {values['v6_palearctic_survival']}\nTropical {values['v6_tropical_survival']}", "#edf4ed"),
        ("Area mechanism\n(H4)", f"{values['h4_promoted']} promoted", "#f4f1ed"),
        ("Common nonlinear\nbreakpoint", f"{values['geometry_promoted']} promoted", "#f4f1ed"),
        ("Channel heterogeneity\n(N1)", f"p = {values['n1_p_value']:.3f}\nstop before N2", "#f4f1ed"),
        ("GloBI source\nbreadth", f"{values['globi_breadth_promoted']} promoted", "#f4f1ed"),
        ("Biotic vs wind\nspecificity", f"β = {values['h5c_estimate']:+.3f}\np = {values['h5c_p_value']:.3f}", "#f4f1ed"),
        ("Distributed\nthresholds", f"{values['h5d_qualified']} design cells qualified", "#f2eef5"),
    ]
    xs = [0.03, 0.27, 0.51, 0.75]
    ys = [0.58, 0.28]
    for i, (title, body, face) in enumerate(cards):
        _box(
            ax, (xs[i % 4], ys[i // 4]), 0.21, 0.22,
            f"{title}\n\n{body}", facecolor=face, fontsize=7.0,
        )
    _box(
        ax, (0.07, 0.045), 0.86, 0.15,
        "STRONG PLANT-SIDE PATTERN\nUPSTREAM POLLINATOR MECHANISM REMAINS UNIDENTIFIED",
        facecolor="#f7f7f7", fontsize=7.8, weight="bold",
    )


def _panel_d(ax: plt.Axes) -> None:
    _setup(ax, "D", "Global assembly depth ≠ local response geometry")
    _box(
        ax, (0.04, 0.24), 0.39, 0.58,
        (
            "CHAPTER 1 — GLOBAL\n\n"
            "many lineages + many island histories\n↓\naveraging\n↓\n"
            "smooth / context-dependent\nassemblage gradient\n\n"
            "answers:\nWHERE?\nWHICH COMPONENTS?\nAT WHAT ASSEMBLY DEPTH?"
        ),
        facecolor="#eef3f8", edgecolor="#496b8a", fontsize=7.9, weight="bold",
    )
    _box(
        ax, (0.57, 0.24), 0.39, 0.58,
        (
            "CHAPTER 2 — izu-core\n\n"
            "one deeply resolved island system\ninteraction state\n↓\n"
            "effective service\n↓\nreproductive outcome\n↓\nphenotype\n\n"
            "compare cline / step / shared breakpoint /\nchannel-specific geometry"
        ),
        facecolor="#f6f0e8", edgecolor="#9a7244", fontsize=7.7, weight="bold",
    )
    _arrow(ax, (0.44, 0.50), (0.56, 0.50), dashed=True, color="0.35", linewidth=1.3)
    ax.text(
        0.50, 0.13,
        "A local threshold can average into a smooth global gradient,\n"
        "but Chapter 1 cannot identify that generator against heterogeneous smooth clines.",
        ha="center", va="center", fontsize=7.7, color="0.30",
    )


def render_figure(
    *,
    output_dir: Path,
    figure3_lock: Path = DEFAULT_FIGURE3_LOCK,
    figure4_lock: Path = DEFAULT_FIGURE4_LOCK,
    n1_lock: Path = DEFAULT_N1_LOCK,
    submission_freeze: Path = DEFAULT_SUBMISSION_FREEZE,
    maximal_integration: Path = DEFAULT_MAXIMAL_INTEGRATION,
) -> dict[str, Any]:
    values = load_inputs(
        figure3_lock=figure3_lock,
        figure4_lock=figure4_lock,
        n1_lock=n1_lock,
        submission_freeze=submission_freeze,
        maximal_integration=maximal_integration,
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(15.0, 10.8), constrained_layout=False)
    grid = fig.add_gridspec(
        2, 2, left=0.045, right=0.975, bottom=0.055, top=0.90,
        hspace=0.27, wspace=0.18,
    )
    _panel_a(fig.add_subplot(grid[0, 0]))
    _panel_b(fig.add_subplot(grid[0, 1]), values)
    _panel_c(fig.add_subplot(grid[1, 0]), values)
    _panel_d(fig.add_subplot(grid[1, 1]))

    fig.suptitle(
        "From an island floral syndrome to a hierarchy-of-assembly test",
        x=0.05, y=0.975, ha="left", fontsize=16, fontweight="bold",
    )
    fig.text(
        0.05, 0.935,
        "Isolation-associated trait composition can arise from repeated lineage response, hierarchical sorting, or both; mechanism is named only after assembly depth and rival explanations are tested.",
        fontsize=9.2, color="0.32",
    )

    basename = "chapter1_v8_figure1_hierarchical_syndrome"
    png = output_dir / f"{basename}.png"
    svg = output_dir / f"{basename}.svg"
    pdf = output_dir / f"{basename}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    manifest = {
        "contract": "chapter1_v8_figure1_hierarchical_syndrome_v1",
        "status": "rendered_from_frozen_v8_evidence_locks",
        "source_locks": values["source_locks"],
        "frozen_results": {
            "support_ladder": values["support_ladder"],
            "genus_attenuation_pct": values["genus_attenuation_pct"],
            "family_to_genus_attenuation_pct": values["family_to_genus_attenuation_pct"],
            "v6_palearctic_survival": values["v6_palearctic_survival"],
            "v6_tropical_survival": values["v6_tropical_survival"],
            "h4_promoted": values["h4_promoted"],
            "geometry_promoted": values["geometry_promoted"],
            "n1_p_value": values["n1_p_value"],
            "n1_pass": values["n1_pass"],
            "globi_breadth_promoted": values["globi_breadth_promoted"],
            "h5c_estimate": values["h5c_estimate"],
            "h5c_p_value": values["h5c_p_value"],
            "h5d_qualified": values["h5d_qualified"],
        },
        "new_biological_models_fitted": False,
        "new_p_values_generated": False,
        "n2_opened": False,
        "claim_boundary": (
            "Figure 1 is an inference map. Family-to-genus attenuation localizes hierarchical "
            "expression but is not causal mediation; failed/non-identified H4, N1, H5c, H5d, "
            "geometry and GloBI alternatives cannot be inverted into evidence that those "
            "mechanisms are absent. Chapter 2 remains the prospective functional-mechanism test."
        ),
        "outputs": [png.name, svg.name, pdf.name],
    }
    (output_dir / "chapter1_v8_figure1_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("render")
def render_command(
    output_dir: Path = typer.Option(...),
    figure3_lock: Path = typer.Option(DEFAULT_FIGURE3_LOCK, exists=True),
    figure4_lock: Path = typer.Option(DEFAULT_FIGURE4_LOCK, exists=True),
    n1_lock: Path = typer.Option(DEFAULT_N1_LOCK, exists=True),
    submission_freeze: Path = typer.Option(DEFAULT_SUBMISSION_FREEZE, exists=True),
    maximal_integration: Path = typer.Option(DEFAULT_MAXIMAL_INTEGRATION, exists=True),
) -> None:
    manifest = render_figure(
        output_dir=output_dir,
        figure3_lock=figure3_lock,
        figure4_lock=figure4_lock,
        n1_lock=n1_lock,
        submission_freeze=submission_freeze,
        maximal_integration=maximal_integration,
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
