"""Efficient runner for trait cascade LOO validation.

Taxonomic LOO predictions are computed once per species-direct cell. The full
support-threshold grid is then evaluated by switching among cached genus and
family predictions, leaving unqualified targets unresolved rather than re-running
the LOO calculation for
every threshold pair.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2.trait_inference_validation import (
    _classification_metrics,
    _load_yaml,
    focused_binary_metrics,
    leave_one_species_out_predictions,
    load_species_direct_cells,
    prepare_direct_cells,
    summarize_prediction_layers,
    support_stratified_metrics,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def threshold_grid_from_cached_predictions(
    predictions: pd.DataFrame,
    genus_thresholds: list[int],
    family_thresholds: list[int],
    genus_min_consensus: float = 1.0,
    family_min_consensus: float = 0.8,
    family_min_supporting_genera: int = 2,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for genus_min in genus_thresholds:
        for family_min in family_thresholds:
            work = predictions.copy()
            use_genus = (
                work["genus_prediction"].astype(str).ne("")
                & pd.to_numeric(work["genus_support_species"], errors="coerce").ge(genus_min)
                & pd.to_numeric(work["genus_top_tie_n"], errors="coerce").eq(1)
                & pd.to_numeric(work["genus_winner_share"], errors="coerce").ge(
                    genus_min_consensus
                )
            )
            use_family = (
                ~use_genus
                & work["family_prediction"].astype(str).ne("")
                & pd.to_numeric(work["family_support_species"], errors="coerce").ge(family_min)
                & pd.to_numeric(work["family_top_tie_n"], errors="coerce").eq(1)
                & pd.to_numeric(work["family_winner_share"], errors="coerce").ge(
                    family_min_consensus
                )
                & pd.to_numeric(
                    work["family_supporting_genera"], errors="coerce"
                ).ge(family_min_supporting_genera)
            )
            work["grid_prediction_tier"] = "unresolved_no_evidence"
            work["grid_prediction"] = ""
            work.loc[use_family, "grid_prediction_tier"] = "family_inference"
            work.loc[use_family, "grid_prediction"] = work.loc[
                use_family, "family_prediction"
            ]
            work.loc[use_genus, "grid_prediction_tier"] = "genus_inference"
            work.loc[use_genus, "grid_prediction"] = work.loc[
                use_genus, "genus_prediction"
            ]

            for trait_name, trait in work.groupby("trait_name", sort=True):
                metrics = _classification_metrics(
                    trait, "actual_value", "grid_prediction"
                )
                tier_counts = trait["grid_prediction_tier"].value_counts().to_dict()
                rows.append(
                    {
                        "trait_name": trait_name,
                        "genus_min_support": int(genus_min),
                        "family_min_support": int(family_min),
                        **metrics,
                        "n_genus_predictions": int(
                            tier_counts.get("genus_inference", 0)
                        ),
                        "n_family_predictions": int(
                            tier_counts.get("family_inference", 0)
                        ),
                        "n_unresolved_predictions": int(
                            tier_counts.get("unresolved_no_evidence", 0)
                        ),
                    }
                )
    return pd.DataFrame(rows)


@app.command("run")
def run(
    cascade_root: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    validation_config: Path = typer.Option(
        Path("config/trait_inference_validation.yml"), exists=True
    ),
    cascade_config: Path = typer.Option(
        Path("config/trait_fill_cascade.yml"), exists=True
    ),
) -> None:
    validation = _load_yaml(validation_config)
    cascade = _load_yaml(cascade_config)
    traits = {str(value) for value in validation["validate_traits"]}
    direct = prepare_direct_cells(load_species_direct_cells(cascade_root, traits))
    genus_min = int(cascade.get("min_genus_support", 1))
    family_min = int(cascade.get("min_family_support", 3))
    family_genera_min = int(cascade.get("min_family_supporting_genera", 2))
    genus_consensus = float(cascade.get("min_genus_consensus", 1.0))
    family_consensus = float(cascade.get("min_family_consensus", 0.8))
    inference_policies = dict(cascade.get("inference_policies") or {})
    predictions = leave_one_species_out_predictions(
        direct,
        genus_min_support=genus_min,
        family_min_support=family_min,
        genus_min_consensus=genus_consensus,
        family_min_consensus=family_consensus,
        family_min_supporting_genera=family_genera_min,
        inference_policies=inference_policies,
    )
    layer_metrics = summarize_prediction_layers(predictions)
    grid = threshold_grid_from_cached_predictions(
        predictions,
        [int(value) for value in validation["cascade_policy"]["genus_support_grid"]],
        [int(value) for value in validation["cascade_policy"]["family_support_grid"]],
        genus_min_consensus=genus_consensus,
        family_min_consensus=family_consensus,
        family_min_supporting_genera=family_genera_min,
    )
    support_metrics = support_stratified_metrics(
        predictions,
        [[int(value) for value in pair] for pair in validation["reporting"]["support_bins"]],
    )
    focused = validation["focused_binary_diagnostic"]
    binary = focused_binary_metrics(
        predictions,
        str(focused["trait_name"]),
        str(focused["positive_value"]),
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    direct.drop(columns=["raw_distribution"]).to_csv(
        output_dir / "species_direct_cells.csv", index=False
    )
    predictions.to_csv(output_dir / "trait_inference_loo_predictions.csv", index=False)
    layer_metrics.to_csv(output_dir / "trait_inference_loo_metrics.csv", index=False)
    grid.to_csv(output_dir / "trait_inference_threshold_grid.csv", index=False)
    support_metrics.to_csv(output_dir / "trait_inference_support_metrics.csv", index=False)
    binary.to_csv(output_dir / "focused_binary_validation.csv", index=False)

    manifest = {
        "contract": str(validation["contract"]),
        "runner": "cached_single_pass_loo",
        "n_species_direct_cells": int(len(direct)),
        "n_direct_species": int(direct["accepted_species"].nunique()),
        "trait_direct_counts": {
            str(key): int(value)
            for key, value in direct["trait_name"].value_counts().sort_index().items()
        },
        "operational_cascade_thresholds": {
            "min_genus_support": genus_min,
            "min_family_support": family_min,
            "min_family_supporting_genera": family_genera_min,
            "min_genus_consensus": genus_consensus,
            "min_family_consensus": family_consensus,
            "trait_specific_policy_count": len(inference_policies),
        },
        "operational_layer_metrics": layer_metrics.to_dict("records"),
        "focused_binary_metrics": binary.to_dict("records"),
        "threshold_grid_rows": int(len(grid)),
        "interpretation": (
            "Leave-one-species-out validation removes each target species from all taxonomic support. "
            "The full support-threshold grid reuses cached LOO predictions and is reported without cherry-picking."
        ),
    }
    (output_dir / "trait_inference_validation_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
