"""Render Chapter 1 v8 Figure 1 as a claim-safe inference map.

Figure 1 is presentation-only. It visualizes the already frozen Chapter 1 inference
architecture and a small set of numerical anchors; it does not refit a biological
model, recompute a test statistic, or promote a mechanism.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import typer  # noqa: E402
from matplotlib.patches import Circle, FancyBboxPatch, Polygon  # noqa: E402

app = typer.Typer(add_completion=False, no_args_is_help=True)

SPEC_PATH = Path("docs/chapter1_figure1_hierarchical_syndrome_spec_20260914.md")
ANCHORS = {
    "support_ladder": "4/4 -> 4/4 -> 0/4",
    "genus_attenuation": "78.8–85.9%",
    "family_to_genus_attenuation": "70.6–79.1%",
    "h5c_p": "p = 0.412",
    "h5d_identifiable": "0/8 identifiable",
}
SPEC_TOKENS = [
    "4/4   ────────>   4/4   ───────────>   0/4",
    "78.8–85.9%",
    "70.6–79.1%",
    "p = 0.412",
    "0/8 identifiable",
    "Do not show a pollinator icon as the cause",
    "family→genus attenuation visually central",
]

INK = "#252525"
MUTED = "#666666"
LINE = "#b8b8b8"
BLUE = "#456f9e"
PURPLE = "#8d6397"
ORANGE = "#c8783d"
TEAL = "#4f8579"
PANEL_BG = "#fafafa"
SOFT_BLUE = "#eef4f9"
SOFT_ORANGE = "#fbf1e9"
SOFT_PURPLE = "#f5eef7"
SOFT_GREY = "#f2f2f2"


class FigureInputError(ValueError):
    """Raised when the canonical Figure 1 design contract has drifted."""


def validate_spec(spec_path: Path) -> str:
    """Require the canonical claim-safe Figure 1 specification and frozen anchors."""
    if not spec_path.is_file():
        raise FigureInputError(f"missing Figure 1 specification: {spec_path}")
    text = spec_path.read_text(encoding="utf-8")
    missing = [token for token in SPEC_TOKENS if token not in text]
    if missing:
        raise FigureInputError(f"Figure 1 specification drifted; missing tokens: {missing}")
    return text


def _panel(ax: plt.Axes, letter: str, title: str) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.add_patch(
        FancyBboxPatch(
            (0.005, 0.01),
            0.99,
            0.98,
            boxstyle="round,pad=0.008,rounding_size=0.02",
            facecolor=PANEL_BG,
            edgecolor="#dddddd",
            linewidth=0.9,
        )
    )
    ax.text(0.025, 0.965, letter, va="top", fontsize=13, fontweight="bold", color=INK)
    ax.text(0.075, 0.965, title, va="top", fontsize=10.8, fontweight="bold", color=INK)


def _box(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    text: str,
    *,
    facecolor: str = "white",
    edgecolor: str = LINE,
    fontsize: float = 8.0,
    weight: str = "normal",
    ha: str = "center",
) -> None:
    x, y = xy
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.008,rounding_size=0.014",
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=1.0,
        )
    )
    tx = x + width / 2 if ha == "center" else x + 0.02
    ax.text(
        tx,
        y + height / 2,
        text,
        ha=ha,
        va="center",
        fontsize=fontsize,
        color=INK,
        fontweight=weight,
        linespacing=1.22,
    )


def _arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    dashed: bool = False,
    color: str = MUTED,
    width: float = 1.0,
) -> None:
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops={
            "arrowstyle": "->",
            "color": color,
            "lw": width,
            "linestyle": "--" if dashed else "-",
            "shrinkA": 1,
            "shrinkB": 1,
        },
    )


def _island(ax: plt.Axes, center: tuple[float, float], scale: float, color: str) -> None:
    cx, cy = center
    pts = [
        (cx - 0.9 * scale, cy - 0.05 * scale),
        (cx - 0.55 * scale, cy + 0.35 * scale),
        (cx - 0.15 * scale, cy + 0.18 * scale),
        (cx + 0.20 * scale, cy + 0.48 * scale),
        (cx + 0.78 * scale, cy + 0.12 * scale),
        (cx + 0.55 * scale, cy - 0.35 * scale),
        (cx - 0.10 * scale, cy - 0.44 * scale),
    ]
    ax.add_patch(Polygon(pts, closed=True, facecolor=color, edgecolor="white", linewidth=0.8))


def _panel_a(ax: plt.Axes) -> None:
    _panel(ax, "A", "Why do floral traits change with isolation?")
    ax.text(
        0.50,
        0.885,
        "ADAPTATION OR HIERARCHICAL ASSEMBLY?",
        ha="center",
        fontsize=11.4,
        fontweight="bold",
        color=INK,
    )

    _island(ax, (0.13, 0.72), 0.07, TEAL)
    _island(ax, (0.87, 0.72), 0.06, ORANGE)
    ax.text(0.13, 0.79, "near source", ha="center", fontsize=7.5, color=MUTED)
    ax.text(0.87, 0.79, "remote island", ha="center", fontsize=7.5, color=MUTED)
    _arrow(ax, (0.22, 0.72), (0.78, 0.72), color=LINE)
    ax.text(0.50, 0.745, "increasing source separation", ha="center", fontsize=7.2, color=MUTED)

    channel_y = 0.625
    for x, label, color in (
        (0.25, "attraction / colour", BLUE),
        (0.50, "accessibility", ORANGE),
        (0.75, "reproductive assurance", PURPLE),
    ):
        ax.add_patch(Circle((x, channel_y), 0.016, facecolor=color, edgecolor="white", linewidth=0.5))
        ax.text(x, channel_y - 0.042, label, ha="center", va="top", fontsize=7.2, color=INK)

    _box(
        ax,
        (0.04, 0.27),
        0.28,
        0.20,
        "Repeated within-lineage response\n\nsame lineage, trait changes",
        facecolor=SOFT_BLUE,
        edgecolor=BLUE,
        fontsize=7.7,
        weight="bold",
    )
    _box(
        ax,
        (0.36, 0.27),
        0.28,
        0.20,
        "Hierarchical lineage sorting\n\nsource pool A B C D E\nisland pool A   C   E",
        facecolor=SOFT_ORANGE,
        edgecolor=ORANGE,
        fontsize=7.4,
        weight="bold",
    )
    _box(
        ax,
        (0.68, 0.27),
        0.28,
        0.20,
        "Mixed generator\n\nsorting + within-lineage\nresponse",
        facecolor=SOFT_PURPLE,
        edgecolor=PURPLE,
        fontsize=7.7,
        weight="bold",
    )
    ax.text(
        0.50,
        0.135,
        "Assemblage means alone cannot distinguish these generators.",
        ha="center",
        fontsize=8.2,
        color=INK,
        fontweight="bold",
    )
    ax.text(
        0.50,
        0.075,
        "No pollinator-specific causal chain is assumed here.",
        ha="center",
        fontsize=7.4,
        color=MUTED,
    )


def _panel_b(ax: plt.Axes) -> None:
    _panel(ax, "B", "The response branches, then collapses at a taxonomic depth")
    ax.text(0.08, 0.86, "ONE UNIVERSAL RESPONSE?", fontsize=8.4, fontweight="bold", color=INK)
    _box(ax, (0.70, 0.815), 0.18, 0.09, "NO", facecolor=SOFT_GREY, fontsize=9, weight="bold")
    _arrow(ax, (0.64, 0.86), (0.69, 0.86))
    _arrow(ax, (0.79, 0.81), (0.79, 0.73))

    _box(
        ax,
        (0.08, 0.58),
        0.80,
        0.14,
        "CONTEXT BRANCHING\nPalearctic: accessibility ↑ + assurance ↑\nTropical: accessibility ↓ + assurance ↑",
        facecolor=SOFT_BLUE,
        edgecolor=BLUE,
        fontsize=8.0,
        weight="bold",
    )
    _arrow(ax, (0.48, 0.575), (0.48, 0.50))
    ax.text(0.08, 0.485, "WHERE DOES THE STRONGEST RESPONSE LIVE?", fontsize=8.3, fontweight="bold")

    xs = [0.19, 0.49, 0.79]
    labels = ["observed", "family-adjusted", "genus-adjusted"]
    support = ["4/4", "4/4", "0/4"]
    for x, label, value in zip(xs, labels, support, strict=True):
        _box(
            ax,
            (x - 0.105, 0.285),
            0.21,
            0.13,
            f"{label}\n{value}",
            facecolor="white" if value != "0/4" else SOFT_ORANGE,
            edgecolor=BLUE if value != "0/4" else ORANGE,
            fontsize=8.0,
            weight="bold",
        )
    _arrow(ax, (0.30, 0.35), (0.375, 0.35), width=1.3)
    _arrow(ax, (0.60, 0.35), (0.675, 0.35), width=2.3, color=ORANGE)
    ax.text(0.64, 0.235, "ASSEMBLY DEPTH", ha="center", fontsize=9.5, color=ORANGE, fontweight="bold")
    ax.text(
        0.64,
        0.185,
        "family → genus",
        ha="center",
        fontsize=8.3,
        color=INK,
        fontweight="bold",
    )
    ax.text(
        0.64,
        0.125,
        "genus attenuation 78.8–85.9%\nconditional drop 70.6–79.1%",
        ha="center",
        fontsize=7.5,
        color=MUTED,
    )


def _evidence_row(
    ax: plt.Axes,
    y: float,
    badge: str,
    label: str,
    detail: str,
    facecolor: str,
    edgecolor: str,
) -> None:
    _box(ax, (0.05, y - 0.035), 0.22, 0.07, badge, facecolor=facecolor, edgecolor=edgecolor, fontsize=7.0, weight="bold")
    ax.text(0.30, y + 0.012, label, va="center", fontsize=7.6, fontweight="bold", color=INK)
    ax.text(0.30, y - 0.025, detail, va="center", fontsize=6.9, color=MUTED)


def _panel_c(ax: plt.Axes) -> None:
    _panel(ax, "C", "Cross-examination: what did not explain the pattern?")
    rows = [
        ("SURVIVES", "Trait missingness (V5)", "Palearctic core survives finite MNAR grid", SOFT_BLUE, BLUE),
        ("SURVIVES", "Species-list bias (V6)", "Palearctic 99/100; remote survey 80/80", SOFT_BLUE, BLUE),
        ("NOT PROMOTED", "Area mechanism (H4)", "0/16 promoted", SOFT_GREY, MUTED),
        ("NOT PROMOTED", "Common breakpoint", "0/12 nonlinear shapes promoted", SOFT_GREY, MUTED),
        ("NOT PROMOTED", "Channel heterogeneity (N1)", "joint p = 0.655; modest effects underpowered", SOFT_GREY, MUTED),
        ("NOT PROMOTED", "GloBI source breadth", "0/4 promoted", SOFT_GREY, MUTED),
        ("NOT PROMOTED", "Biotic vs wind specificity", "interaction +0.065, p = 0.412", SOFT_GREY, MUTED),
        ("NON-IDENTIFIABLE", "Distributed thresholds", "0/8 identifiable; false calls 19–25.5%", SOFT_PURPLE, PURPLE),
    ]
    for y, row in zip([0.84, 0.75, 0.66, 0.57, 0.48, 0.39, 0.30, 0.21], rows, strict=True):
        _evidence_row(ax, y, *row)
    _box(
        ax,
        (0.07, 0.055),
        0.86,
        0.085,
        "Strong plant-side pattern; upstream pollinator mechanism still unidentified.",
        facecolor=SOFT_ORANGE,
        edgecolor=ORANGE,
        fontsize=8.0,
        weight="bold",
    )


def _panel_d(ax: plt.Axes) -> None:
    _panel(ax, "D", "Global assembly depth ≠ local response geometry")
    _box(
        ax,
        (0.05, 0.50),
        0.37,
        0.33,
        "CHAPTER 1 · GLOBAL\n\nmany lineages + histories\n↓\nassemblage gradient\n\nWHERE?\nWHICH COMPONENTS?\nAT WHAT ASSEMBLY DEPTH?",
        facecolor=SOFT_BLUE,
        edgecolor=BLUE,
        fontsize=7.7,
        weight="bold",
    )
    _box(
        ax,
        (0.58, 0.50),
        0.37,
        0.33,
        "CHAPTER 2 · LOCAL\n\ninteraction state\n↓\neffective service\n↓\nreproductive outcome\n↓\nphenotype",
        facecolor=SOFT_PURPLE,
        edgecolor=PURPLE,
        fontsize=7.7,
        weight="bold",
    )
    _arrow(ax, (0.43, 0.66), (0.57, 0.66), dashed=True, color=ORANGE, width=1.4)
    ax.text(0.50, 0.70, "scale bridge", ha="center", fontsize=7.0, color=ORANGE, fontweight="bold")
    _box(
        ax,
        (0.08, 0.25),
        0.84,
        0.14,
        "A local threshold can average into a smooth global gradient,\nbut Chapter 1 cannot identify that generator against heterogeneous smooth clines.",
        facecolor=SOFT_ORANGE,
        edgecolor=ORANGE,
        fontsize=7.6,
        weight="bold",
    )
    ax.text(
        0.50,
        0.12,
        "MECHANISM REMAINS TO BE RESOLVED LOCALLY",
        ha="center",
        fontsize=9.6,
        color=INK,
        fontweight="bold",
    )


def render_figure(*, spec_path: Path, output_dir: Path) -> dict[str, Any]:
    validate_spec(spec_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(14.0, 10.0))
    fig.subplots_adjust(left=0.035, right=0.985, top=0.91, bottom=0.045, hspace=0.16, wspace=0.08)
    _panel_a(axes[0, 0])
    _panel_b(axes[0, 1])
    _panel_c(axes[1, 0])
    _panel_d(axes[1, 1])

    fig.suptitle(
        "From an island floral syndrome to a hierarchy-of-assembly test",
        x=0.04,
        y=0.975,
        ha="left",
        fontsize=16,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.04,
        0.94,
        "Direction and assembly depth are identified before any named mechanism is promoted.",
        fontsize=9.5,
        color=MUTED,
    )

    basename = "chapter1_v8_figure1_hierarchical_inference"
    png_path = output_dir / f"{basename}.png"
    svg_path = output_dir / f"{basename}.svg"
    pdf_path = output_dir / f"{basename}.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    manifest = {
        "contract": "chapter1_v8_figure1_hierarchical_inference_v1",
        "status": "rendered_from_frozen_inference_spec",
        "design_spec": str(spec_path),
        "new_biological_models_fitted": False,
        "new_p_values_generated": False,
        "numerical_anchors": ANCHORS,
        "claim_boundary": (
            "Figure 1 is an inference map. It localizes the strongest broad response at the "
            "family-to-genus transition and summarizes adverse tests without promoting an "
            "upstream pollinator-specific causal mechanism or treating genus attenuation as "
            "absence of evolution."
        ),
        "solid_arrow_meaning": "relationship directly established or organizational flow",
        "dashed_arrow_meaning": "mechanistic or cross-scale link remains unresolved/prospective",
        "outputs": [png_path.name, svg_path.name, pdf_path.name],
    }
    (output_dir / "chapter1_v8_figure1_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("render")
def render_command(
    output_dir: Path = typer.Option(...),
    spec_path: Path = typer.Option(SPEC_PATH, exists=True),
) -> None:
    manifest = render_figure(spec_path=spec_path, output_dir=output_dir)
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
