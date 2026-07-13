"""Run exploratory M0-M4 path models on the already-materialized integrated dataset.

This runner intentionally performs no new acquisition and applies no Bombus applicability
filter. It is a path-analysis snapshot using the existing island trait compositions and
Bombus environmental compatibility. Production observation-process refinement remains a
separate later analysis.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

from island_v2.m0_m4_sem import _load_config, build_sem_island_data, fit_model_ladder

app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.command("run")
def run(
    trait_composition_csv: Path = typer.Option(..., exists=True),
    island_metrics_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/m0_m4_sem.yml"), exists=True),
    analysis_tier: str = typer.Option("primary"),
) -> None:
    """Fit M0-M4 paths to all islands available in the existing integrated dataset."""
    config = _load_config(config_path)
    composition = pd.read_csv(trait_composition_csv)
    metrics = pd.read_csv(island_metrics_csv)

    island_data, metadata = build_sem_island_data(composition, metrics, config, covariates=None)
    baseline = metadata["available_baseline_covariates"]
    paths, fit, indirect = fit_model_ladder(island_data, baseline)

    output_dir.mkdir(parents=True, exist_ok=True)
    island_data.to_csv(output_dir / "path_island_data.csv", index=False)
    paths.to_csv(output_dir / "path_coefficients.csv", index=False)
    fit.to_csv(output_dir / "path_model_fit.csv", index=False)
    indirect.to_csv(output_dir / "path_indirect_effects.csv", index=False)

    manifest = {
        "contract": "m0_m4_existing_data_paths_v1",
        "analysis_tier": analysis_tier,
        "n_islands": int(island_data["island_id"].nunique()),
        "baseline_covariates_used": baseline,
        "missing_baseline_covariates": metadata["missing_baseline_covariates"],
        "bombus_predictor": config["bombus_channel"]["column"],
        "scope": "all islands present in the existing integrated dataset; no applicability gate",
        "interpretation": (
            "Exploratory path analysis using existing Bombus environmental compatibility and "
            "trait-composition proxies. No new Apidae acquisition or observation-process correction."
        ),
    }
    (output_dir / "path_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
