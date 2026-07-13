from __future__ import annotations

import json

import numpy as np
import pandas as pd

from island_v2.m1_probabilistic_trait_uncertainty import (
    aggregate_species_states,
    build_probability_mapping,
    run_imputations,
)


def _config() -> dict:
    return {
        "analysis_contract": "test_probabilistic_m1",
        "analysis_role": "sensitivity",
        "outcome": {
            "trait_name": "pollination_functional_guild",
            "positive_value": "bumblebees",
            "included_fill_tiers": [
                "species_direct",
                "genus_inference",
                "family_inference",
            ],
        },
        "imputation": {
            "n_imputations": 20,
            "random_seed": 20260713,
        },
        "baseline_specs": {
            "full": [
                "log_island_area_km2",
                "abs_latitude",
                "climate_pc1",
                "climate_pc2",
                "climate_pc3",
                "climate_pc4",
                "log_trait_species_richness",
            ],
            "no_abs_latitude": [
                "log_island_area_km2",
                "climate_pc1",
                "climate_pc2",
                "climate_pc3",
                "climate_pc4",
                "log_trait_species_richness",
            ],
        },
        "bombus_predictor": "bombus_channel_state",
        "spatial_cluster": "spatial_block",
        "nominal_z_threshold": 1.96,
        "reporting": {"minimum_imputations": 20},
    }


def _synthetic_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(404)
    species = [f"Species {index:03d}" for index in range(60)]
    species_probabilities = {
        name: float(rng.beta(1.5, 2.5)) for name in species
    }
    rows = []
    reference = []
    for island_index in range(80):
        island_id = f"island_{island_index:03d}"
        bombus = float(rng.uniform(0, 1))
        reference.append(
            {
                "island_id": island_id,
                "spatial_block": f"block_{island_index // 4:02d}",
                "bombus_channel_state": bombus,
                "log_island_area_km2": np.log1p(10 + island_index),
                "abs_latitude": 5 + island_index % 70,
                "climate_pc1": np.sin(island_index / 8),
                "climate_pc2": np.cos(island_index / 10),
                "climate_pc3": np.sin(island_index / 6),
                "climate_pc4": np.cos(island_index / 7),
                "log_trait_species_richness": np.log1p(50 + island_index % 30),
            }
        )
        # Higher Bombus compatibility makes higher-p species more likely to occur.
        scores = np.array(
            [species_probabilities[name] * (0.5 + 2.0 * bombus) + rng.normal(0, 0.08) for name in species]
        )
        selected = np.argsort(scores)[-30:]
        for species_index in selected:
            name = species[species_index]
            probability = species_probabilities[name]
            if species_index < 10:
                fill_tier = "species_direct"
                filled = "bumblebees" if probability >= 0.5 else "other_bees"
                distribution = {filled: 1.0}
            elif species_index < 35:
                fill_tier = "genus_inference"
                filled = "bumblebees" if probability >= 0.5 else "other_bees"
                distribution = {
                    "bumblebees": probability,
                    "other_bees": 1.0 - probability,
                }
            else:
                fill_tier = "family_inference"
                filled = "bumblebees" if probability >= 0.5 else "other_bees"
                distribution = {
                    "bumblebees": probability,
                    "other_bees": 1.0 - probability,
                }
            rows.append(
                {
                    "island_id": island_id,
                    "accepted_species": name,
                    "trait_name": "pollination_functional_guild",
                    "filled_value": filled,
                    "fill_tier": fill_tier,
                    "value_distribution": json.dumps(distribution),
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(reference)


def test_probability_mapping_draws_one_state_per_species_across_islands() -> None:
    master, reference = _synthetic_tables()
    mapping, metadata = build_probability_mapping(master, reference, _config())
    assert metadata["n_species"] > 20
    assert metadata["n_islands"] == 80
    state = np.arange(len(mapping.species_table)) % 2
    successes = aggregate_species_states(mapping, state.astype(float))
    manual = pd.DataFrame(
        {
            "island_code": mapping.island_codes,
            "species_code": mapping.species_codes,
        }
    )
    manual["state"] = manual["species_code"].map(
        dict(enumerate(state.astype(float)))
    )
    expected = (
        manual.groupby("island_code")["state"]
        .sum()
        .reindex(range(len(mapping.island_table)), fill_value=0)
        .to_numpy(dtype=float)
    )
    assert np.array_equal(successes, expected)


def test_probabilistic_imputation_returns_pooled_and_comparator_results() -> None:
    master, reference = _synthetic_tables()
    mapping, _ = build_probability_mapping(master, reference, _config())
    draws, pooled, comparators = run_imputations(mapping, _config())
    assert draws["imputation"].nunique() == 20
    assert set(pooled["baseline_spec"]) == {"full", "no_abs_latitude"}
    assert set(comparators["comparator"]) == {
        "expected_probability_counts",
        "hard_modal_broad_without_global_fallback",
    }
    assert np.isfinite(pooled["pooled_estimate_log_odds"]).all()
    assert (pooled["n_imputations"] == 20).all()
