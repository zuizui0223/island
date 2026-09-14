"""Render Chapter 1 v8 Figure 1 from the frozen Figure 2–4 result locks.

Figure 1 is an inference map, not a new analysis. It validates the already frozen
result locks, then visualizes the biological question, rival generators, assembly
depth, falsification boundaries, and the scale bridge to izu-core.
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


class FigureInputError(ValueError):
    """Raised when a frozen result lock is missing or inconsistent."""


def _read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FigureInputError(f"missing frozen result lock: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def load_locks(repo_root: Path) -> dict[str, dict[str, Any]]:
    fig2 = _read_json(repo_root / "config/chapter1_v8_figure2_result_lock.json")
    fig3 = _read_json(repo_root / "config/chapter1_v8_figure3_submission_result_lock.json")
    fig4 = _read_json(repo_root / "config/chapter1_v8_figure4_result_lock.json")

    if fig2.get("contract") != "chapter1_v8_figure2_result_lock_v1":
        raise FigureInputError("unexpected Figure 2 lock contract")
    if fig3.get("contract") != "chapter1_v8_figure3_submission_result_lock_v1":
        raise FigureInputError("unexpected Figure 3 lock contract")
    if fig4.get("contract") != "chapter1_v8_figure4_result_lock_v1":
        raise FigureInputError("unexpected Figure 4 lock contract")

    f3 = fig3.get("frozen_results", {})
    if f3.get("support_ladder") != "4/4 -> 4/4 -> 0/4":
        raise FigureInputError("Figure 3 support ladder differs from frozen manuscript claim")
    if f3.get("family_attenuation_fraction_range") != [0.196, 0.334]:
        raise FigureInputError("Figure 3 family attenuation range changed")
    if f3.get("genus_attenuation_fraction_range") != [0.788, 0.859]:
        raise FigureInputError("Figure 3 genus attenuation range changed")
    if f3.get("conditional_family_to_genus_attenuation_range") != [0.706, 0.791]:
        raise FigureInputError("Figure 3 conditional attenuation range changed")

    f4 = fig4.get("frozen_results", {})
    expected = {
        "v6_palearctic_survival": "99/100",
        "v6_north_tropical_vector_survival": "70/75",
        "v6_tropical_survival": "35/75",
        "geometry_promoted": "0/12",
        "h5d_qualified": "0/8",
    }
    for key, value in expected.items():
        if f4.get(key) != value:
            raise FigureInputError(f"Figure 4 frozen result changed: {key}")
    if abs(float(f4.get("h5c_p_value", float("nan"))) - 0.4122068222871858) > 1e-12:
        raise FigureInputError("Figure 4 H5c p-value changed")

    return {"figure2": fig2, "figure3": fig3, "figure4": fig4}


def _box(
    ax: plt.Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    text: str,
    *,
    fontsize: float = 8.2,
    linewidth: float = 1.0,
    fill: str = "white",
    edge: str = "0.35",
    weight: str = "normal",
) -> FancyBboxPatch:
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.015",
        linewidth=linewidth,
        edgecolor=edge,
        facecolor=fill,
        transform=ax.transAxes,
    )
    ax.add_patch(patch)
    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=fontsize,
        fontweight=weight,
        linespacing=1.25,
    )
    return patch


def _arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    dashed: bool = False,
    linewidth: float = 1.2,
) -> None:
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=9,
        linewidth=linewidth,
        color="0.35",
        linestyle="--" if dashed else "-",
        transform=ax.transAxes,
        connectionstyle="arc3,rad=0",
    )
    ax.add_patch(arrow)


def _panel_a(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_title(
        "A  Why do floral traits change with isolation?",
        loc="left",
        fontsize=11,
        fontweight="bold",
        pad=8,
    )
    ax.text(
        0.02,
        0.91,
        "The same assemblage pattern can arise from different biological generators.",
        transform=ax.transAxes,
        fontsize=8.3,
        color="0.35",
    )

    _box(ax, 0.04, 0.66, 0.24, 0.14, "SOURCE POOL\nA   B   C   D   E", fontsize=9, weight="bold")
    _box(ax, 0.70, 0.66, 0.24, 0.14, "REMOTE ISLAND\nA   C   E", fontsize=9, weight="bold")
    _arrow(ax, (0.29, 0.73), (0.69, 0.73))
    ax.text(0.49, 0.77, "geographic filter", transform=ax.transAxes, ha="center", fontsize=8, color="0.35")

    _box(ax, 0.03, 0.37, 0.28, 0.17, "WITHIN-LINEAGE RESPONSE\nlineage persists\ntrait value changes", fontsize=8.2)
    _box(ax, 0.36, 0.37, 0.28, 0.17, "HIERARCHICAL SORTING\nlineages differ in\narrival / establishment / persistence", fontsize=8.2, linewidth=1.5, weight="bold")
    _box(ax, 0.69, 0.37, 0.28, 0.17, "MIXED GENERATOR\nsorting +\nwithin-lineage response", fontsize=8.2)
    _arrow(ax, (0.82, 0.65), (0.83, 0.55))
    ax.text(
        0.50,
        0.20,
        "An assemblage mean alone cannot distinguish adaptation from membership change.",
        transform=ax.transAxes,
        ha="center",
        fontsize=9,
        fontweight="bold",
    )
    ax.text(
        0.50,
        0.09,
        "Primary channels: attraction / colour   ·   structural accessibility   ·   reproductive assurance",
        transform=ax.transAxes,
        ha="center",
        fontsize=8.1,
        color="0.35",
    )


def _panel_b(ax: plt.Axes, fig3: dict[str, Any]) -> None:
    ax.set_axis_off()
    ax.set_title(
        "B  The response branches, then attenuates at a taxonomic depth",
        loc="left",
        fontsize=11,
        fontweight="bold",
        pad=8,
    )
    f3 = fig3["frozen_results"]
    _box(ax, 0.04, 0.76, 0.26, 0.12, "ONE UNIVERSAL\nRESPONSE?", fontsize=9, weight="bold")
    _arrow(ax, (0.31, 0.82), (0.40, 0.82))
    ax.text(0.355, 0.855, "NO", transform=ax.transAxes, ha="center", fontsize=8.5, fontweight="bold")
    _box(
        ax,
        0.41,
        0.71,
        0.54,
        0.22,
        "BIOGEOGRAPHIC BRANCHING\nPalearctic: accessibility ↑  +  assurance ↑\nTropical: accessibility ↓  +  assurance ↑\n(tropical accessibility is less observation-robust)",
        fontsize=8.5,
        weight="bold",
    )
    _arrow(ax, (0.68, 0.70), (0.68, 0.61))

    ax.text(0.05, 0.58, "ASSEMBLY DEPTH", transform=ax.transAxes, fontsize=9.5, fontweight="bold")
    stages = [("Observed", "4/4"), ("Family-adjusted", "4/4"), ("Genus-adjusted", "0/4")]
    xs = [0.08, 0.39, 0.70]
    for i, ((label, score), x) in enumerate(zip(stages, xs, strict=True)):
        _box(ax, x, 0.36, 0.22, 0.14, f"{label}\n{score}", fontsize=8.7, weight="bold")
        if i < 2:
            _arrow(ax, (x + 0.23, 0.43), (xs[i + 1] - 0.01, 0.43))
    ax.plot([0.49, 0.80], [0.32, 0.32], transform=ax.transAxes, color="0.3", linewidth=2.4)
    ax.text(
        0.645,
        0.255,
        "largest attenuation: family → genus",
        transform=ax.transAxes,
        ha="center",
        fontsize=9,
        fontweight="bold",
    )
    ax.text(
        0.645,
        0.17,
        f"genus attenuation {100*f3['genus_attenuation_fraction_range'][0]:.1f}–{100*f3['genus_attenuation_fraction_range'][1]:.1f}%\n"
        f"of observed vector; conditional family→genus {100*f3['conditional_family_to_genus_attenuation_range'][0]:.1f}–{100*f3['conditional_family_to_genus_attenuation_range'][1]:.1f}%",
        transform=ax.transAxes,
        ha="center",
        fontsize=8.2,
    )
    ax.text(
        0.05,
        0.07,
        "This localizes hierarchical expression; it does not prove dispersal alone or absence of within-lineage change.",
        transform=ax.transAxes,
        fontsize=7.8,
        color="0.35",
    )


def _panel_c(ax: plt.Axes, fig4: dict[str, Any]) -> None:
    ax.set_axis_off()
    ax.set_title(
        "C  Cross-examination: what survives, what does not identify a mechanism?",
        loc="left",
        fontsize=11,
        fontweight="bold",
        pad=8,
    )
    f4 = fig4["frozen_results"]
    cards = [
        ("Trait missingness (V5)", "Palearctic core survives finite MNAR stress"),
        ("Species-list bias (V6)", f"Palearctic {f4['v6_palearctic_survival']} survives; tropical {f4['v6_tropical_survival']}"),
        ("Area mechanism (H4)", "0/16 promoted"),
        ("Common nonlinear geometry", f"{f4['geometry_promoted']} promoted"),
        ("Pollinator channel heterogeneity", "N1 p = 0.655; modest effects underpowered"),
        ("Sampled GloBI source breadth", "0/4 promoted"),
        ("Biotic vs wind specificity", f"interaction +{f4['h5c_interaction_estimate']:.3f}; p = {f4['h5c_p_value']:.3f}"),
        ("Distributed thresholds", f"{f4['h5d_qualified']} identifiable; false threshold calls 19–25.5%"),
    ]
    cols = [0.03, 0.51]
    top = 0.82
    h = 0.145
    gap = 0.018
    for idx, (label, outcome) in enumerate(cards):
        col = idx // 4
        row = idx % 4
        y = top - row * (h + gap)
        _box(ax, cols[col], y, 0.45, h, f"{label}\n{outcome}", fontsize=7.9, linewidth=1.0)
    _box(
        ax,
        0.15,
        0.05,
        0.70,
        0.10,
        "STRONG PLANT-SIDE PATTERN; UPSTREAM POLLINATOR MECHANISM REMAINS UNIDENTIFIED",
        fontsize=8.7,
        linewidth=1.6,
        weight="bold",
    )


def _panel_d(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_title(
        "D  Global assembly depth ≠ local response geometry",
        loc="left",
        fontsize=11,
        fontweight="bold",
        pad=8,
    )
    _box(
        ax,
        0.03,
        0.48,
        0.39,
        0.35,
        "CHAPTER 1 — GLOBAL\nmany lineages + many island histories\n↓ averaging\ncontext-dependent assemblage gradient\n\nanswers:\nWHERE?  WHICH COMPONENTS?\nAT WHAT TAXONOMIC DEPTH?",
        fontsize=8.3,
        weight="bold",
    )
    _box(
        ax,
        0.58,
        0.43,
        0.39,
        0.45,
        "CHAPTER 2 — IZU-CORE\none deeply resolved island system\n\ninteraction state\n↓\neffective service\n↓\nreproductive outcome\n↓\nphenotype\n\ncompare cline / step / shared breakpoint / channel geometry",
        fontsize=8.1,
        weight="bold",
    )
    _arrow(ax, (0.43, 0.66), (0.57, 0.66), dashed=True, linewidth=1.5)
    ax.text(
        0.50,
        0.30,
        "A local threshold can average into a smooth global gradient,\nbut Chapter 1 cannot identify that generator against heterogeneous smooth clines.",
        transform=ax.transAxes,
        ha="center",
        fontsize=8.4,
        fontweight="bold",
    )
    ax.text(
        0.50,
        0.12,
        "MECHANISM REMAINS TO BE RESOLVED LOCALLY",
        transform=ax.transAxes,
        ha="center",
        fontsize=10,
        fontweight="bold",
    )


def render_figure(*, repo_root: Path, output_dir: Path) -> dict[str, Any]:
    locks = load_locks(repo_root)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(16.0, 10.8), constrained_layout=False)
    grid = fig.add_gridspec(
        2,
        2,
        left=0.04,
        right=0.985,
        bottom=0.055,
        top=0.90,
        hspace=0.28,
        wspace=0.16,
    )
    axes = [fig.add_subplot(grid[i, j]) for i in range(2) for j in range(2)]
    _panel_a(axes[0])
    _panel_b(axes[1], locks["figure3"])
    _panel_c(axes[2], locks["figure4"])
    _panel_d(axes[3])

    fig.suptitle(
        "From ‘why are island flowers dull?’ to a hierarchy-of-assembly test",
        x=0.045,
        y=0.975,
        ha="left",
        fontsize=16,
        fontweight="bold",
    )
    fig.text(
        0.045,
        0.94,
        "A visible floral syndrome can reflect repeated lineage response, hierarchical sorting, or both; mechanism is promoted only with independent evidence.",
        fontsize=9.5,
        color="0.30",
    )

    basename = "chapter1_v8_figure1_hierarchical_syndrome_inference"
    png = output_dir / f"{basename}.png"
    svg = output_dir / f"{basename}.svg"
    pdf = output_dir / f"{basename}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    manifest = {
        "contract": "chapter1_v8_figure1_inference_map_v1",
        "status": "rendered_from_frozen_figure2_figure3_figure4_locks",
        "source_locks": {
            "figure2": locks["figure2"]["contract"],
            "figure3": locks["figure3"]["contract"],
            "figure4": locks["figure4"]["contract"],
        },
        "frozen_support_ladder": locks["figure3"]["frozen_results"]["support_ladder"],
        "genus_attenuation_fraction_range": locks["figure3"]["frozen_results"]["genus_attenuation_fraction_range"],
        "v6_palearctic_survival": locks["figure4"]["frozen_results"]["v6_palearctic_survival"],
        "v6_tropical_survival": locks["figure4"]["frozen_results"]["v6_tropical_survival"],
        "geometry_promoted": locks["figure4"]["frozen_results"]["geometry_promoted"],
        "h5c_p_value": locks["figure4"]["frozen_results"]["h5c_p_value"],
        "h5d_qualified": locks["figure4"]["frozen_results"]["h5d_qualified"],
        "new_biological_models_fitted": False,
        "new_p_values_generated": False,
        "claim_boundary": (
            "Figure 1 is a conceptual inference map backed by existing frozen result locks. "
            "It does not identify pollinator causation, prove evolution absent, or infer a global threshold."
        ),
        "outputs": [png.name, svg.name, pdf.name],
    }
    (output_dir / "chapter1_v8_figure1_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command()
def main(
    repo_root: Path = typer.Option(Path("."), exists=True, file_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    typer.echo(json.dumps(render_figure(repo_root=repo_root, output_dir=output_dir), indent=2))


if __name__ == "__main__":
    app()
