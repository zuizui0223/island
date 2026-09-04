import numpy as np
import pandas as pd

from island_v2.chapter1_pr138_si_genus_source_null import run_si_genus_source_null


def _configs():
    pattern = {
        "geography_column": "log_distance_to_continent_km",
        "cluster_column": "spatial_block",
        "baseline_covariates": ["log_island_area_km2", "climate_pc1"],
    }
    source = {
        "source_region_scores": {
            "minimum_trait_scored_species_per_region_syndrome": 1,
        },
        "source_assignment": {"primary_modes": {"geo_k5": {}}},
        "response": {"source_expectation_requires_minimum_source_regions": 1},
    }
    return pattern, source


def _common_tables(species, island_species):
    ledger = pd.DataFrame(
        {
            "accepted_species": species,
            "trait_name": "self_incompatibility",
            "normalized_value": "SI",
        }
    )
    status_rows = []
    for island_id, members in island_species.items():
        for sp in members:
            status_rows.append(
                {
                    "island_id": island_id,
                    "accepted_species": sp,
                    "origin_status": "native",
                    "floristic_status": "native_nonendemic",
                }
            )
    status = pd.DataFrame(status_rows)
    n = len(island_species)
    ids = list(island_species)
    distance = np.linspace(-1.0, 1.0, n)
    covariates = pd.DataFrame(
        {
            "island_id": ids,
            "log_distance_to_continent_km": distance,
            "log_island_area_km2": np.sin(np.arange(n) / 7),
            "climate_pc1": np.cos(np.arange(n) / 9),
            "spatial_block": [f"b{i // 6}" for i in range(n)],
        }
    )
    realm = pd.DataFrame(
        {"island_id": ids, "biogeographic_realm": "Palearctic"}
    )
    gift_rows = []
    for entity in [1, 2]:
        for sp in species:
            gift_rows.append({"entity_ID": entity, "work_species": sp})
    gift = pd.DataFrame(gift_rows)
    assignments = pd.DataFrame(
        [
            {
                "island_id": island_id,
                "source_mode": "geo_k5",
                "entity_ID": entity,
            }
            for island_id in ids
            for entity in [1, 2]
        ]
    )
    return ledger, status, gift, assignments, covariates, realm


def _score_table(attraction_by_species):
    rows = []
    for species, attraction in attraction_by_species.items():
        rows.extend(
            [
                {
                    "accepted_species": species,
                    "syndrome": "large_bee_like",
                    "syndrome_concordance": -attraction,
                },
                {
                    "accepted_species": species,
                    "syndrome": "generalized_accessible",
                    "syndrome_concordance": attraction,
                },
            ]
        )
    return pd.DataFrame(rows)


def test_genus_composition_can_absorb_si_restricted_slope():
    species = ["Alpha one", "Alpha two", "Beta one", "Beta two"]
    attraction = {
        "Alpha one": 0.8,
        "Alpha two": 0.8,
        "Beta one": -0.8,
        "Beta two": -0.8,
    }
    island_species = {}
    n = 60
    for i in range(n):
        island_species[f"i{i}"] = (
            ["Beta one", "Beta two"] if i < n // 2 else ["Alpha one", "Alpha two"]
        )
    ledger, status, gift, assignments, covariates, realm = _common_tables(
        species, island_species
    )
    scores = _score_table(attraction)
    pattern, source = _configs()
    outputs = run_si_genus_source_null(
        ledger,
        scores,
        ledger,
        scores,
        status,
        gift,
        assignments,
        covariates,
        realm,
        pattern,
        source,
        n_draws=200,
        seed=13801,
    )
    result = outputs["results"]
    assert (result["observed_distance_estimate"] > 0).all()
    assert np.allclose(result["residual_distance_estimate"], 0.0, atol=1e-10)
    assert np.allclose(
        result["observed_distance_estimate"], result["null_distance_mean"], atol=1e-10
    )
    assert outputs["manifest"]["source_pool_recomputed_each_draw"] is True


def test_within_genus_sorting_survives_genus_null():
    species = [
        "Alpha low", "Alpha high",
        "Beta low", "Beta high",
        "Gamma low", "Gamma high",
    ]
    attraction = {
        "Alpha low": -0.8,
        "Alpha high": 0.8,
        "Beta low": -0.6,
        "Beta high": 0.6,
        "Gamma low": -0.7,
        "Gamma high": 0.7,
    }
    island_species = {}
    n = 60
    for i in range(n):
        island_species[f"i{i}"] = (
            ["Alpha low", "Beta low", "Gamma low"]
            if i < n // 2
            else ["Alpha high", "Beta high", "Gamma high"]
        )
    ledger, status, gift, assignments, covariates, realm = _common_tables(
        species, island_species
    )
    scores = _score_table(attraction)
    pattern, source = _configs()
    outputs = run_si_genus_source_null(
        ledger,
        scores,
        ledger,
        scores,
        status,
        gift,
        assignments,
        covariates,
        realm,
        pattern,
        source,
        n_draws=400,
        seed=13802,
    )
    result = outputs["results"]
    assert (result["observed_distance_estimate"] > 0).all()
    assert (result["residual_distance_estimate"] > 0.5).all()
    assert (result["empirical_one_sided_p_observed_gt_null"] < 0.05).all()
