from __future__ import annotations

import pandas as pd

from island_v2.m0_m4_sem import build_sem_island_data, fit_model_ladder


def _config() -> dict[str, object]:
    return {
        "proxies": {
            "floral_phenotype_proxy": {
                "trait_name": "pollination_functional_guild",
                "positive_values": ["bumblebees"],
            },
            "reproductive_assurance_proxy": {
                "components": [
                    {"trait_name": "autonomous_selfing_capacity", "positive_values": ["autonomous"]},
                    {"trait_name": "self_incompatibility", "positive_values": ["SC", "likely_SC"]},
                ]
            },
            "alternative_functional_replacement_proxy": {
                "trait_name": "pollination_functional_guild",
                "positive_values": ["other_bees", "flies", "birds", "bats", "wind", "mixed"],
            },
        },
        "bombus_channel": {"column": "environmental_compatibility_max"},
        "baseline_candidates": ["n_trait_species", "log_island_area_km2"],
    }


def test_build_sem_island_data_uses_declared_trait_values() -> None:
    composition = pd.DataFrame(
        {
            "island_id": ["i1", "i1", "i1", "i1", "i1", "i2", "i2", "i2", "i2"],
            "trait_name": [
                "pollination_functional_guild",
                "pollination_functional_guild",
                "autonomous_selfing_capacity",
                "self_incompatibility",
                "self_incompatibility",
                "pollination_functional_guild",
                "autonomous_selfing_capacity",
                "self_incompatibility",
                "self_incompatibility",
            ],
            "filled_value": [
                "bumblebees",
                "other_bees",
                "autonomous",
                "SC",
                "SI",
                "flies",
                "autonomous",
                "likely_SC",
                "SI",
            ],
            "species_share_within_trait": [0.4, 0.6, 0.5, 0.25, 0.75, 1.0, 0.2, 0.4, 0.6],
        }
    )
    metrics = pd.DataFrame(
        {
            "island_id": ["i1", "i2"],
            "environmental_compatibility_max": [0.8, 0.2],
            "n_trait_species": [10, 20],
        }
    )
    result, metadata = build_sem_island_data(composition, metrics, _config())
    indexed = result.set_index("island_id")
    assert indexed.loc["i1", "floral_phenotype_proxy"] == 0.4
    assert indexed.loc["i1", "alternative_functional_replacement_proxy"] == 0.6
    assert indexed.loc["i1", "reproductive_assurance_proxy"] == 0.375
    assert indexed.loc["i2", "floral_phenotype_proxy"] == 0.0
    assert metadata["available_baseline_covariates"] == ["n_trait_species"]
    assert metadata["missing_baseline_covariates"] == ["log_island_area_km2"]


def test_model_ladder_emits_mediation_and_interaction_models() -> None:
    rows = []
    for index in range(40):
        bombus = index / 39
        mediator = 0.2 + 0.5 * bombus
        alternative = 1.0 - bombus
        outcome = 0.1 + 0.3 * bombus + 0.4 * mediator - 0.2 * alternative
        rows.append(
            {
                "island_id": f"i{index}",
                "n_trait_species": 50 + index,
                "bombus_channel_state": bombus,
                "reproductive_assurance_proxy": mediator,
                "alternative_functional_replacement_proxy": alternative,
                "floral_phenotype_proxy": outcome,
            }
        )
    data = pd.DataFrame(rows)
    paths, fit, indirect = fit_model_ladder(data, ["n_trait_species"])
    assert set(fit["model"]) == {"M0", "M1", "M2", "M3", "M4"}
    assert "bombus_x_alternative" in set(paths["predictor"])
    assert set(indirect["model"]) == {"M2", "M3"}
    assert indirect["indirect_effect"].notna().all()
