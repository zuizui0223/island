import pandas as pd

from island_v2.chapter1_pr138_source_pool_sensitivity import (
    build_mainland_syndrome_scores,
    build_source_assignments,
    build_source_expectations,
    match_gift_species_to_frozen_scores,
)


def _config():
    return {
        "source_assignment": {
            "primary_modes": {
                "geo_k5": {"rule": "five_nearest_source_region_representative_points", "k": 5},
                "geo_k10": {"rule": "ten_nearest_source_region_representative_points", "k": 10},
                "geo_k20": {"rule": "twenty_nearest_source_region_representative_points", "k": 20},
                "geo50_climate10": {
                    "rule": "ten_climate_nearest_among_fifty_geographically_nearest_source_regions",
                    "geographic_candidate_k": 50,
                    "k": 10,
                },
            }
        }
    }


def test_species_matching_uses_exact_then_unique_binomial_and_fails_closed_on_ambiguity():
    scores = pd.DataFrame(
        {
            "accepted_species": ["Alpha beta", "Gamma delta subsp. one", "Gamma delta subsp. two"],
            "syndrome": ["x", "x", "x"],
            "syndrome_concordance": [1.0, 0.5, -0.5],
        }
    )
    flora = pd.DataFrame(
        {
            "entity_ID": [1, 1, 1],
            "work_species": ["Alpha beta", "Gamma delta", "Missing species"],
        }
    )
    expanded, audit = match_gift_species_to_frozen_scores(flora, scores)
    assert set(expanded["accepted_species"]) == {"Alpha beta"}
    counts = audit.set_index("match_method")["n_gift_entity_species_rows"].to_dict()
    assert counts["exact_normalized"] == 1
    assert counts["unmatched"] == 2


def test_mainland_scores_and_expectations_are_equal_weighted_across_selected_regions():
    matched = pd.DataFrame(
        {
            "entity_ID": [1, 1, 2, 2, 3, 3],
            "accepted_species": ["A a", "B b", "C c", "D d", "E e", "F f"],
            "syndrome": ["x"] * 6,
            "syndrome_concordance": [1.0, 1.0, 0.0, 0.0, -1.0, -1.0],
        }
    )
    mainland = build_mainland_syndrome_scores(matched, min_species=2)
    assignments = pd.DataFrame(
        {
            "island_id": ["i1", "i1", "i1"],
            "source_mode": ["m", "m", "m"],
            "source_rank": [1, 2, 3],
            "entity_ID": [1, 2, 3],
            "source_distance_km": [10.0, 20.0, 30.0],
            "climate_distance": [0.1, 0.2, 0.3],
        }
    )
    expected = build_source_expectations(assignments, mainland, min_source_regions=3)
    row = expected.iloc[0]
    assert row["n_source_regions"] == 3
    assert abs(row["source_expectation"] - 0.0) < 1e-12
    assert bool(row["source_expectation_eligible"])


def test_source_assignment_depends_only_on_geography_and_climate_columns():
    cov = pd.DataFrame(
        {
            "island_id": ["i1"],
            "island_latitude": [0.0],
            "island_longitude": [0.0],
            "climate_pc1": [0.0],
            "climate_pc2": [0.0],
            "climate_pc3": [0.0],
            "climate_pc4": [0.0],
            "syndrome_score": [999.0],
        }
    )
    source = pd.DataFrame(
        {
            "entity_ID": list(range(1, 31)),
            "source_latitude": [0.0] * 30,
            "source_longitude": [float(x) for x in range(1, 31)],
            "climate_pc1": [float(x) / 10 for x in range(30)],
            "climate_pc2": [0.0] * 30,
            "climate_pc3": [0.0] * 30,
            "climate_pc4": [0.0] * 30,
            "mainland_syndrome_score": list(reversed(range(30))),
        }
    )
    first = build_source_assignments(cov, source, _config())
    cov.loc[0, "syndrome_score"] = -999.0
    source["mainland_syndrome_score"] = range(30)
    second = build_source_assignments(cov, source, _config())
    pd.testing.assert_frame_equal(first.reset_index(drop=True), second.reset_index(drop=True))
    geo5 = first.loc[first["source_mode"].eq("geo_k5")]
    assert geo5["entity_ID"].tolist() == [1, 2, 3, 4, 5]
