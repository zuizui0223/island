import numpy as np
import pandas as pd

from island_v2.chapter1_area_capacity_moderation import (
    _fit_between,
    _fit_within,
    build_families,
)

AREA_CONFIG = {
    "distance_column": "distance",
    "area_column": "area",
    "cluster_column": "cluster",
    "control_columns": ["climate"],
}


def _long_family(
    *,
    interactions: dict[str, float],
    n_per_context: int = 150,
) -> pd.DataFrame:
    rng = np.random.default_rng(13927)
    rows = []
    for context_index, (context, interaction) in enumerate(interactions.items()):
        for i in range(n_per_context):
            island = f"{context}_{i}"
            distance = -1.5 + 3.0 * i / (n_per_context - 1)
            area = np.sin(i / 8.0) + 0.05 * distance
            climate = np.cos(i / 11.0)
            cluster = f"{context}_b{i // 5}"
            for response, offset in [("access", 0.0), ("assurance", 0.08)]:
                score = (
                    (0.35 + offset) * distance
                    + 0.05 * area
                    + interaction * distance * area
                    + 0.03 * climate
                    + rng.normal(0, 0.025)
                )
                rows.append(
                    {
                        "family": "plant",
                        "source_mode": "not_applicable",
                        "island_id": island,
                        "stratum": "all_native",
                        "response": response,
                        "response_score": score,
                        "context": context,
                        "distance": distance,
                        "area": area,
                        "climate": climate,
                        "cluster": cluster,
                    }
                )
    return pd.DataFrame(rows)


def test_negative_distance_by_area_is_classified_as_small_island_amplification_signal():
    data = _long_family(interactions={"north": -0.55})
    coefficients, omnibus = _fit_within(
        data,
        AREA_CONFIG,
        family="plant",
        source_mode="not_applicable",
        stratum="all_native",
        context="north",
        context_column="context",
        context_layer="analysis_regime",
        support_tier="confirmatory",
        threshold=50,
    )
    assert omnibus["status"] == "fit"
    assert omnibus["p_value"] < 0.05
    assert coefficients["distance_estimate_at_mean_area"].gt(0).all()
    assert coefficients["distance_x_area_estimate"].lt(0).all()
    assert (
        coefficients["distance_slope_at_small_area_z_minus1"]
        > coefficients["distance_slope_at_large_area_z_plus1"]
    ).all()


def test_direct_context_test_recovers_difference_in_area_moderation():
    data = _long_family(interactions={"north": -0.55, "tropical": 0.20})
    differences, omnibus = _fit_between(
        data,
        AREA_CONFIG,
        family="plant",
        source_mode="not_applicable",
        stratum="all_native",
        context_a="north",
        context_b="tropical",
        context_column="context",
        context_layer="analysis_regime",
        support_tier="confirmatory",
        threshold=50,
    )
    assert omnibus["status"] == "fit"
    assert omnibus["p_value"] < 0.05
    assert differences["distance_x_area_difference_b_minus_a"].gt(0).all()
    assert differences["difference_p"].lt(0.05).all()


def test_support_gate_fails_closed_without_two_supported_responses():
    data = _long_family(interactions={"north": -0.4}, n_per_context=40)
    coefficients, omnibus = _fit_within(
        data,
        AREA_CONFIG,
        family="plant",
        source_mode="not_applicable",
        stratum="all_native",
        context="north",
        context_column="context",
        context_layer="analysis_regime",
        support_tier="confirmatory",
        threshold=50,
    )
    assert coefficients.empty
    assert omnibus["status"] == "not_testable"


def test_lineage_family_uses_only_predeclared_scope_and_matching():
    lineage = pd.DataFrame(
        [
            {
                "island_id": "i1",
                "stratum": "all_native",
                "source_mode": "geo_k5",
                "evidence_scope": scope,
                "source_matching": matching,
                "minimum_represented_genera": minimum,
                "entry_enrichment": 0.2,
                "loading_increment": 0.1,
            }
            for scope, matching, minimum in [
                ("broad", "prevalence_richness", 5),
                ("broad", "prevalence_only", 5),
                ("exact_si_direct", "prevalence_richness", 3),
            ]
        ]
    )
    area_config = {
        "families": {
            "source_lineage_broad": {
                "kind": "lineage_scores",
                "evidence_scope": "broad",
                "source_matching": "prevalence_richness",
                "minimum_represented_genera": 5,
                "responses": ["entry_enrichment", "loading_increment"],
            }
        }
    }
    families = build_families(
        pd.DataFrame(), lineage, {"branch_axes": {}}, area_config
    )
    result = families["source_lineage_broad"]
    assert len(result) == 2
    assert set(result["response"]) == {"entry_enrichment", "loading_increment"}
