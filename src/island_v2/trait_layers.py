"""Assign trait candidates to M0/M1/M2 layers and an analysis-eligibility track.

This governs analysis eligibility only. It never decides, alters, or emits a
trait value; it reads a reviewed trait-candidate table and annotates each row
with its trait layer (M0/M1/M2) and analysis track:

- ``main_direct``          — species-direct evidence (plus species-indirect for
  M0) eligible for the primary analysis;
- ``expanded_sensitivity`` — hierarchical (genus/family) inference, or
  species-level evidence below a layer's main bar, reported only as a labelled
  sensitivity track;
- ``excluded``             — unresolved / below-track evidence.

Reproductive-assurance and pollinator-guild traits are restricted to
species-direct evidence in the main analysis; their genus/family inference is
always isolated to the sensitivity track.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

TRACK_MAIN = "main_direct"
TRACK_SENSITIVITY = "expanded_sensitivity"
TRACK_EXCLUDED = "excluded"
UNCLASSIFIED_LAYER = "unclassified_trait"

_SPECIES_LEVEL_SCOPES = {"species_direct", "species_indirect"}
_HIERARCHICAL_SCOPES = {"genus_inference", "family_inference", "genus_direct", "family_direct"}


@app.callback()
def main() -> None:
    """Annotate reviewed trait candidates with layer and analysis track."""


def load_layers(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or "layers" not in config:
        raise typer.BadParameter("trait layers config must contain a layers mapping")
    return config


def trait_to_layer(config: dict[str, Any]) -> dict[str, str]:
    """Map every trait name to its layer key."""
    mapping: dict[str, str] = {}
    for layer_key, layer in config["layers"].items():
        for trait in layer.get("traits", []) or []:
            mapping[str(trait)] = layer_key
    return mapping


def analysis_track(
    trait: str,
    candidate_kind: str,
    evidence_scope: str,
    config: dict[str, Any],
) -> tuple[str, str]:
    """Return ``(trait_layer, analysis_track)`` for a candidate.

    Fail-closed: hierarchical inference is never main; anything not clearly
    species-direct/indirect at a layer's accepted scope drops out of the main
    analysis rather than being admitted by default.
    """
    layer_map = trait_to_layer(config)
    layer = layer_map.get(str(trait).strip())
    if layer is None:
        return UNCLASSIFIED_LAYER, TRACK_EXCLUDED

    kind = str(candidate_kind or "").strip().lower()
    scope = str(evidence_scope or "").strip().lower()
    accepted = {str(s).lower() for s in config["layers"][layer].get("main_analysis_accepted_scopes", [])}

    if kind == "hierarchical_inference" or scope in _HIERARCHICAL_SCOPES:
        return layer, TRACK_SENSITIVITY
    if scope in accepted:
        return layer, TRACK_MAIN
    if scope in _SPECIES_LEVEL_SCOPES:
        # Species-level but below this layer's main bar (e.g. species_indirect for M1).
        return layer, TRACK_SENSITIVITY
    return layer, TRACK_EXCLUDED


def annotate_candidates(candidates: pd.DataFrame, config: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Add trait_layer and analysis_track columns; return the frame and a summary."""
    trait_col = next((c for c in ("trait_name", "trait") if c in candidates.columns), None)
    if trait_col is None:
        raise typer.BadParameter("candidate table must have a 'trait_name' or 'trait' column")
    table = candidates.copy()
    kind = table["candidate_kind"] if "candidate_kind" in table.columns else pd.Series([""] * len(table))
    scope = table["evidence_scope"] if "evidence_scope" in table.columns else pd.Series([""] * len(table))

    layers: list[str] = []
    tracks: list[str] = []
    for trait_value, kind_value, scope_value in zip(
        table[trait_col].astype(str), kind.astype(str), scope.astype(str), strict=False
    ):
        layer, track = analysis_track(trait_value, kind_value, scope_value, config)
        layers.append(layer)
        tracks.append(track)
    table["trait_layer"] = layers
    table["analysis_track"] = tracks

    def _counts(series: pd.Series) -> dict[str, int]:
        return {str(k): int(v) for k, v in series.value_counts().sort_index().items()}

    summary = {
        "note": (
            "Analysis eligibility only; no trait value is decided or altered. "
            "Hierarchical inference is isolated to expanded_sensitivity and is never "
            "a primary result. Reproductive-assurance and pollinator-guild traits "
            "require species-direct evidence for the main analysis."
        ),
        "n_candidates": int(len(table)),
        "by_layer": _counts(table["trait_layer"]),
        "by_track": _counts(table["analysis_track"]),
        "n_main_direct": int((table["analysis_track"] == TRACK_MAIN).sum()),
        "n_expanded_sensitivity": int((table["analysis_track"] == TRACK_SENSITIVITY).sum()),
        "n_excluded": int((table["analysis_track"] == TRACK_EXCLUDED).sum()),
    }
    return table, summary


@app.command("annotate")
def annotate(
    candidates_csv: Path = typer.Option(..., exists=True, help="Reviewed trait-candidate CSV (trait_name/trait)."),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/trait_layers.yml")),
) -> None:
    """Write layer- and track-annotated candidates plus a summary."""
    config = load_layers(config_path)
    candidates = pd.read_csv(candidates_csv, dtype=str).fillna("")
    table, summary = annotate_candidates(candidates, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_dir / "trait_candidates_layered.csv", index=False)
    (output_dir / "trait_layer_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    typer.echo(
        f"Annotated {summary['n_candidates']} candidate(s): "
        f"{summary['n_main_direct']} main_direct, "
        f"{summary['n_expanded_sensitivity']} expanded_sensitivity, "
        f"{summary['n_excluded']} excluded."
    )


if __name__ == "__main__":
    app()
