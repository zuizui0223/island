import numpy as np
import pandas as pd

from island_v2.chapter1_pr138_si_genus_influence import (
    FULL_SENTINEL,
    run_evidence_mode,
)


def _synthetic_inputs():
    genus_values = {"A": -1.0, "B": -0.5, "C": 0.5, "D": 1.0}
    species = []
    ledger = []
    scores = []
    for genus, value in genus_values.items():
        for index in range(5):
            species_name = f"{genus} species{index}"
            species.append(species_name)
            ledger.append(
                {
                    "accepted_species": species_name,
                    "trait_name": "self_incompatibility",
                    "normalized_value": "SI",
                }
            )
            scores.extend(
                [
                    {
                        "accepted_species": species_name,
                        "syndrome": "large_bee_like",
                        "syndrome_concordance": -value,
                    },
                    {
                        "accepted_species": species_name,
                        "syndrome": "generalized_accessible",
                        "syndrome_concordance": value,
                    },
                ]
            )

    flora = []
    covariates = []
    realms = []
    assignments = []
    for index in range(80):
        island_id = f"i{index}"
        distance = -1 + 2 * index / 79
        genera = ["A", "B"] if distance < 0 else ["C", "D"]
        if index % 4 == 0:
            genera += ["C"] if distance < 0 else ["B"]
        for genus in genera:
            for species_index in range(5):
                flora.append(
                    {
                        "island_id": island_id,
                        "accepted_species": f"{genus} species{species_index}",
                        "origin_status": "native",
                        "floristic_status": "native_nonendemic",
                    }
                )
        covariates.append(
            {
                "island_id": island_id,
                "log_distance_to_continent_km": distance,
                "log_island_area_km2": np.sin(index / 7),
                "spatial_block": f"b{index // 5}",
            }
        )
        realms.append(
            {"island_id": island_id, "biogeographic_realm": "Palearctic"}
        )
        assignments.append(
            {"island_id": island_id, "source_mode": "geo_k5", "entity_ID": 1}
        )

    gift_flora = [{"entity_ID": 1, "work_species": name} for name in species]
    return tuple(
        map(
            pd.DataFrame,
            [ledger, scores, flora, gift_flora, assignments, covariates, realms],
        )
    )


def test_leave_one_genus_preserves_distributed_positive_signal():
    ledger, scores, flora, gift, assignments, covariates, realms = _synthetic_inputs()
    pattern = {
        "geography_column": "log_distance_to_continent_km",
        "baseline_covariates": ["log_island_area_km2"],
        "cluster_column": "spatial_block",
    }
    source = {
        "source_region_scores": {
            "minimum_trait_scored_species_per_region_syndrome": 1,
        },
        "response": {"source_expectation_requires_minimum_source_regions": 1},
        "source_assignment": {"primary_modes": {"geo_k5": {}}},
    }
    outputs = run_evidence_mode(
        ledger,
        scores,
        flora,
        gift,
        assignments,
        covariates,
        realms,
        pattern,
        source,
        evidence_mode="all_analysis_eligible",
    )

    results = outputs["results"]
    ranking = outputs["ranking"]
    manifest = outputs["manifest"]
    assert FULL_SENTINEL in set(results["omitted_genus"])
    assert manifest["n_scored_si_genera_tested"] == 4

    deletions = results.loc[
        ~results["omitted_genus"].eq(FULL_SENTINEL) & results["status"].eq("fit")
    ]
    assert len(deletions) == 8
    assert (deletions["distance_estimate"] > 0).all()
    assert set(ranking["omitted_genus"]) == {"A", "B", "C", "D"}
    assert (ranking["n_positive"] == 2).all()
