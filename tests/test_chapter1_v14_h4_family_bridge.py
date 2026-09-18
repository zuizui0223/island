import numpy as np
import pandas as pd

from island_v2.chapter1_v14_h4_family_bridge import (
    build_family_scores,
    run_family_bridge,
)


MEASUREMENT = {
    "PL_Effect_Size_Type1": "FS",
    "PL_Effect_Size_Type2": "Sup",
    "Constant_added": "False",
    "Level_of_Supplementation": "flower",
}


def _preflight(components: list[str], n: int = 40) -> pd.DataFrame:
    rows = []
    for i in range(n):
        state = float(i % 2)
        for trait in components:
            rows.append(
                {
                    "species_key": f"sp{i}",
                    "trait": trait,
                    "trait_state": state,
                }
            )
    return pd.DataFrame(rows)


def _effects(n: int = 40) -> pd.DataFrame:
    rows = []
    distance = np.linspace(-1.0, 1.0, n)
    for i in range(n):
        score = float(i % 2)
        rows.append(
            {
                "row_id": i,
                "species_key": f"sp{i}",
                "site_key": f"site{i}",
                "study_key": f"study{i}",
                "analysis_regime": "northern_midlatitude",
                "z_distance": distance[i],
                "PL_Effect_Size": 0.2 * distance[i] - 0.8 * score,
                **MEASUREMENT,
            }
        )
    return pd.DataFrame(rows)


def test_build_family_scores_requires_two_components():
    frame = pd.DataFrame(
        [
            {"species_key": "a", "trait": "x", "trait_state": 1},
            {"species_key": "a", "trait": "y", "trait_state": 0},
            {"species_key": "b", "trait": "x", "trait_state": 1},
        ]
    )
    score = build_family_scores(
        frame,
        components=["x", "y", "z"],
        score_name="score",
        minimum_components=2,
    )
    assert score["species_key"].tolist() == ["a"]
    assert score.loc[0, "n_components"] == 2
    assert score.loc[0, "score"] == 0.5


def test_family_bridge_recovers_negative_family_effects():
    ra_components = [
        "self_compatibility",
        "selfing_mating_system",
        "autonomous_selfing",
    ]
    acc_components = [
        "generalized_form",
        "actinomorphic_symmetry",
        "shallow_open_tube",
    ]
    config = {
        "analysis": {
            "publication_total_weight": 1.0,
            "minimum_nonmissing_components_per_species": 2,
        },
        "families": {
            "reproductive_assurance": {
                "score_name": "reproductive_assurance_score",
                "parent": "reproductive_assurance",
                "components": ra_components,
            },
            "accessibility_generalization": {
                "score_name": "accessibility_generalization_score",
                "parent": "floral_architecture",
                "components": acc_components,
            },
        },
    }
    results, manifest = run_family_bridge(
        _preflight(ra_components),
        _effects(),
        _preflight(acc_components),
        _effects(),
        config,
    )
    primary = results.loc[results["analysis"].eq("primary")]
    assert len(primary) == 2
    assert primary["evaluable"].astype(bool).all()
    assert primary["estimate"].lt(0).all()
    assert manifest["all_primary_estimates_negative"] is True
