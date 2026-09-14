"""Scope adapter for the frozen V6 species-detection tipping analysis.

The final PR142 syndrome artifact contains additional floristic strata that V6 does not
analyse. The core V6 implementation intentionally fails closed if an affected score lacks
its matching recorded-richness denominator. This adapter materializes an ephemeral view
of the pinned artifact containing only the V6-contracted strata for syndrome scores while
symlinking every other immutable input unchanged.

No biological value, sensitivity grid, model, or claim rule is changed.
"""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2.chapter1_species_detection_tipping import load_config, run_sensitivity

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _link(source: Path, target: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    os.symlink(source.resolve(), target)


def build_scoped_artifact_view(
    *,
    artifact_root: Path,
    config: dict[str, Any],
    view_root: Path,
) -> dict[str, int]:
    """Create a minimal immutable-input view with syndrome rows limited to V6 strata."""
    static_paths = [
        "fixed/canonical/input/chapter1_status_flora.csv.gz",
        "fixed/isolation/results/purpose_shortest_island_data.csv",
        "fixed/realm/realm/island_biogeographic_realm_assignment.csv",
        "fixed/canonical/observation_bias/observation_selection_coefficients.csv",
    ]
    for relative in static_paths:
        _link(artifact_root / relative, view_root / relative)

    strata = {str(x) for x in config["strata"]}
    counts: dict[str, int] = {}
    for short_scope in ("all", "direct"):
        source = artifact_root / f"syndrome/{short_scope}/island_syndrome_scores.csv.gz"
        if not source.is_file():
            raise FileNotFoundError(source)
        frame = pd.read_csv(source)
        if "stratum" not in frame.columns:
            raise ValueError(f"syndrome score table lacks stratum: {source}")
        scoped = frame.loc[frame["stratum"].astype(str).isin(strata)].copy()
        if scoped.empty:
            raise ValueError(f"no V6 strata remain in {source}")
        if set(scoped["stratum"].astype(str)) != strata:
            raise ValueError(
                f"not every frozen V6 stratum is present in {source}: "
                f"{sorted(set(scoped['stratum'].astype(str)))}"
            )
        target = view_root / f"syndrome/{short_scope}/island_syndrome_scores.csv.gz"
        target.parent.mkdir(parents=True, exist_ok=True)
        scoped.to_csv(target, index=False, compression="gzip")
        counts[short_scope] = int(len(scoped))
    return counts


def run_scoped(
    *,
    artifact_root: Path,
    config_path: Path,
    pattern_config_path: Path,
    branching_config_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    config = load_config(config_path)
    with tempfile.TemporaryDirectory(prefix="chapter1-v6-") as tmp:
        view_root = Path(tmp)
        scope_counts = build_scoped_artifact_view(
            artifact_root=artifact_root,
            config=config,
            view_root=view_root,
        )
        manifest = run_sensitivity(
            artifact_root=view_root,
            config_path=config_path,
            pattern_config_path=pattern_config_path,
            branching_config_path=branching_config_path,
            output_dir=output_dir,
        )
    manifest["scoped_syndrome_rows"] = scope_counts
    manifest["scope_adapter"] = "all_native_and_native_nonendemic_only"
    (output_dir / "chapter1_species_detection_tipping_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("run")
def run_command(
    artifact_root: Path = typer.Option(..., exists=True, file_okay=False),
    config_path: Path = typer.Option(
        Path("config/chapter1_species_detection_tipping.yml"), exists=True, dir_okay=False
    ),
    pattern_config_path: Path = typer.Option(
        Path("config/chapter1_pr136_biogeographic_pattern.yml"), exists=True, dir_okay=False
    ),
    branching_config_path: Path = typer.Option(
        Path("config/chapter1_global_branching.yml"), exists=True, dir_okay=False
    ),
    output_dir: Path = typer.Option(...),
) -> None:
    typer.echo(
        json.dumps(
            run_scoped(
                artifact_root=artifact_root,
                config_path=config_path,
                pattern_config_path=pattern_config_path,
                branching_config_path=branching_config_path,
                output_dir=output_dir,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    app()
