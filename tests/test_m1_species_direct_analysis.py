from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.m1_species_direct_analysis import (
    build_species_direct_data,
    fit_inference_grid,
    spatial_cv,
)


def _config() -> dict:
    return {
        "analysis_contract": "test_species_direct",
        "outcome": {
            "trait_name": "pollination_functional_guild",
            "required_fill_tier": "species_direct",
            "positive_value": "bumblebees",
        },
        "bombus_predictor": {"column": "bombus_channel_state"},
        "baseline_full": [
            "log_island_area_km2",
            "abs_latitude",
            "climate_pc1",
            "climate_pc2",
            "climate_pc3",
            "climate_pc4",
            "log_trait_species_richness",
            "log_gbif_flora_records",
        ],
        "baseline_sensitivity_no_abs_latitude": [
            "log_island_area_km2",
            "climate_pc1",
            "climate_pc2",
            "climate_pc3",
            "climate_pc4",
            "log_trait_species_richness",
            "log_gbif_flora_records",
        ],
        "spatial_cv": {
            "folds": 5,
            "repeats": 2,
            "seed": 20260713,
            "weight_schemes": ["row_weighted", "equal_species_weighted"],
        },
        "reporting": {
            "nominal_z_threshold": 1.96,
            "minimum_rows": 100,
            "minimum_species_clusters": 10,
            "minimum_family_clusters": 5,
            "minimum_spatial_blocks": 5,
        },
    }


def _synthetic_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(17)
    species = [f"Species {i:03d}" for i in range(80)]
    species_family = {name: f"Family {i % 16:02d}" for i, name in enumerate(species)}
    species_trait = {name: int(rng.random() < 0.42) for name in species}

    reference_rows = []
    master_rows = []
    for island_index in range(100):
        island_id = f"island_{island_index:03d}"
        block = f"block_{island_index // 5:02d}"
        bombus = float(rng.uniform(0, 1))
        reference_rows.append(
            {
                "island_id": island_id,
                "spatial_block": block,
                "bombus_channel_state": bombus,
                "log_island_area_km2": np.log1p(10 + island_index),
                "abs_latitude": 5 + island_index % 70,
                "climate_pc1": np.sin(island_index / 9),
                "climate_pc2": np.cos(island_index / 11),
                "climate_pc3": np.sin(island_index / 6),
                "climate_pc4": np.cos(island_index / 8),
                "log_trait_species_richness": np.log1p(50 + island_index % 30),
            }
        )

        # Bombus-compatible islands preferentially contain bumblebee-guild species.
        scores = np.array(
            [
                0.2 + 1.8 * bombus * species_trait[name] + rng.normal(0, 0.15)
                for name in species
            ]
        )
        chosen = np.argsort(scores)[-35:]
        for species_index in chosen:
            name = species[species_index]
            # One direct outcome row.
            master_rows.append(
                {
                    "island_id": island_id,
                    "accepted_species": name,
                    "family": species_family[name],
                    "n_unique_gbif_ids": int(rng.integers(1, 8)),
                    "trait_name": "pollination_functional_guild",
                    "filled_value": "bumblebees" if species_trait[name] else "other_bees",
                    "fill_tier": "species_direct",
                }
            )
            # An additional trait row repeats occurrence metadata and tests effort de-duplication.
            master_rows.append(
                {
                    "island_id": island_id,
                    "accepted_species": name,
                    "family": species_family[name],
                    "n_unique_gbif_ids": master_rows[-1]["n_unique_gbif_ids"],
                    "trait_name": "mating_system",
                    "filled_value": "mixed_mating",
                    "fill_tier": "species_direct",
                }
            )
    return pd.DataFrame(master_rows), pd.DataFrame(reference_rows)


def test_species_direct_builder_and_equal_species_weights() -> None:
    master, reference = _synthetic_data()
    data, metadata = build_species_direct_data(master, reference, _config())
    assert metadata["n_rows"] == 3500
    assert metadata["n_islands"] == 100
    assert metadata["n_species"] >= 20
    assert metadata["n_families"] >= 10
    species_weight_totals = data.groupby("accepted_species")["equal_species_weighted"].sum()
    assert np.allclose(species_weight_totals, species_weight_totals.iloc[0])
    assert np.isfinite(data["log_gbif_flora_records"]).all()


def test_two_way_clustered_inference_recovers_positive_bombus_composition_signal() -> None:
    master, reference = _synthetic_data()
    data, _ = build_species_direct_data(master, reference, _config())
    coefficients, fits = fit_inference_grid(data, _config())
    assert fits["converged"].all()
    primary = coefficients.loc[
        coefficients["model"].eq("M1")
        & coefficients["baseline_spec"].eq("full")
        & coefficients["weight_scheme"].eq("row_weighted")
        & coefficients["covariance_scheme"].eq("species_x_spatial_block")
        & coefficients["predictor"].eq("bombus_channel_state")
    ]
    assert len(primary) == 1
    assert float(primary.iloc[0]["estimate_log_odds"]) > 0
    assert int(primary.iloc[0]["n_cluster_a"]) >= 20
    assert int(primary.iloc[0]["n_cluster_b"]) == 20


def test_species_direct_spatial_cv_uses_whole_blocks_and_detects_incremental_signal() -> None:
    master, reference = _synthetic_data()
    data, _ = build_species_direct_data(master, reference, _config())
    folds, summary = spatial_cv(data, _config())
    assert folds["repeat"].nunique() == 2
    assert folds["fold"].nunique() == 5
    assert set(summary["weight_scheme"]) == {
        "row_weighted",
        "equal_species_weighted",
    }
    row_full = summary.loc[
        summary["baseline_spec"].eq("full")
        & summary["weight_scheme"].eq("row_weighted")
    ].iloc[0]
    assert row_full["proportion_repeats_m1_better_loglik"] >= 0.5
