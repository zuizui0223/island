import numpy as np
import pandas as pd

from island_v2.chapter1_pr138_regional_lineage_decomposition import (
    AXIS_COMPONENTS,
    build_regional_lineage_responses,
)


def test_genus_decomposition_is_zero_when_species_equal_their_genus_means():
    species_scores = pd.DataFrame(
        [
            {
                "accepted_species": species,
                "syndrome": syndrome,
                "syndrome_concordance": score,
            }
            for species, genus_scores in {
                "Alpha one": {"generalized_accessible": 0.6, "selfing_core": -0.2},
                "Alpha two": {"generalized_accessible": 0.6, "selfing_core": -0.2},
                "Beta one": {"generalized_accessible": -0.4, "selfing_core": 0.4},
                "Beta two": {"generalized_accessible": -0.4, "selfing_core": 0.4},
            }.items()
            for syndrome, score in genus_scores.items()
        ]
    )
    status_flora = pd.DataFrame(
        [
            {
                "island_id": island,
                "accepted_species": species,
                "origin_status": "native",
                "endemic_status": "nonendemic",
                "floristic_status": "native_nonendemic",
            }
            for island in ["i1", "i2"]
            for species in ["Alpha one", "Alpha two", "Beta one", "Beta two"]
        ]
    )
    gift_flora = pd.DataFrame(
        {
            "entity_ID": [1, 1, 1, 1],
            "work_species": ["Alpha one", "Alpha two", "Beta one", "Beta two"],
        }
    )
    source_assignments = pd.DataFrame(
        {
            "island_id": ["i1", "i2"],
            "source_mode": ["geo_k1", "geo_k1"],
            "source_rank": [1, 1],
            "entity_ID": [1, 1],
            "source_distance_km": [10.0, 20.0],
        }
    )
    # Island and source floras are identical, so both the observed source-adjusted
    # score and the deterministic genus-mean source-adjusted score are exactly zero.
    source_adjusted_scores = pd.DataFrame(
        [
            {
                "source_mode": "geo_k1",
                "island_id": island,
                "syndrome": component,
                "syndrome_score": 0.0,
                "n_species": 4,
                "stratum": stratum,
            }
            for island in ["i1", "i2"]
            for component in AXIS_COMPONENTS.values()
            for stratum in ["all_native", "native_nonendemic"]
        ]
    )
    source_config = {
        "source_region_scores": {
            "minimum_trait_scored_species_per_region_syndrome": 1,
        },
        "response": {"source_expectation_requires_minimum_source_regions": 1},
        "source_assignment": {"primary_modes": {"geo_k1": {"k": 1}}},
    }

    responses, audit = build_regional_lineage_responses(
        source_adjusted_scores,
        species_scores,
        status_flora,
        gift_flora,
        source_assignments,
        source_config,
    )

    assert not responses.empty
    assert set(responses["axis"]) == set(AXIS_COMPONENTS)
    assert np.allclose(responses["observed"].to_numpy(float), 0.0)
    assert np.allclose(responses["genus_expected"].to_numpy(float), 0.0)
    assert np.allclose(responses["residual"].to_numpy(float), 0.0)
    assert not audit.empty


def test_genus_decomposition_preserves_frozen_observed_scores_exactly():
    species_scores = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "syndrome": syndrome,
                "syndrome_concordance": score,
            }
            for syndrome, score in [
                ("generalized_accessible", 0.25),
                ("selfing_core", -0.10),
            ]
        ]
    )
    status_flora = pd.DataFrame(
        {
            "island_id": ["i1"],
            "accepted_species": ["Alpha one"],
            "origin_status": ["native"],
            "endemic_status": ["nonendemic"],
            "floristic_status": ["native_nonendemic"],
        }
    )
    gift_flora = pd.DataFrame({"entity_ID": [1], "work_species": ["Alpha one"]})
    source_assignments = pd.DataFrame(
        {
            "island_id": ["i1"],
            "source_mode": ["geo_k1"],
            "source_rank": [1],
            "entity_ID": [1],
            "source_distance_km": [10.0],
        }
    )
    frozen = {
        "generalized_accessible": 0.137,
        "selfing_core": -0.221,
    }
    source_adjusted_scores = pd.DataFrame(
        [
            {
                "source_mode": "geo_k1",
                "island_id": "i1",
                "syndrome": component,
                "syndrome_score": value,
                "n_species": 1,
                "stratum": stratum,
            }
            for component, value in frozen.items()
            for stratum in ["all_native", "native_nonendemic"]
        ]
    )
    source_config = {
        "source_region_scores": {
            "minimum_trait_scored_species_per_region_syndrome": 1,
        },
        "response": {"source_expectation_requires_minimum_source_regions": 1},
        "source_assignment": {"primary_modes": {"geo_k1": {"k": 1}}},
    }

    responses, _ = build_regional_lineage_responses(
        source_adjusted_scores,
        species_scores,
        status_flora,
        gift_flora,
        source_assignments,
        source_config,
    )

    for axis, component in AXIS_COMPONENTS.items():
        values = responses.loc[responses["axis"].eq(axis), "observed"].to_numpy(float)
        assert np.allclose(values, frozen[component])
