import numpy as np
import pandas as pd

from island_v2.chapter1_v14_h4_exact_h2_score_bridge import (
    exact_h2_scores,
    run_exact_h2_bridge,
)


MEASUREMENT = {
    "PL_Effect_Size_Type1": "FS",
    "PL_Effect_Size_Type2": "Sup",
    "Constant_added": "False",
    "Level_of_Supplementation": "flower",
}


def _scores(n: int = 40) -> pd.DataFrame:
    rows = []
    for i in range(n):
        value = float(i % 2)
        for syndrome in ("selfing_core", "generalized_accessible"):
            rows.append(
                {
                    "accepted_species": f"Sp_{i}",
                    "syndrome": syndrome,
                    "soft_membership": value,
                }
            )
    return pd.DataFrame(rows)


def _effects(n: int = 40) -> pd.DataFrame:
    distance = np.linspace(-1.0, 1.0, n)
    rows = []
    for i in range(n):
        score = float(i % 2)
        rows.append(
            {
                "row_id": i,
                "species_key": f"sp {i}",
                "site_key": f"site{i}",
                "study_key": f"study{i}",
                "analysis_regime": "northern_midlatitude",
                "z_distance": distance[i],
                "PL_Effect_Size": 0.2 * distance[i] - 0.8 * score,
                **MEASUREMENT,
            }
        )
    return pd.DataFrame(rows)


def test_exact_h2_scores_normalize_species_names():
    result = exact_h2_scores(
        _scores(2),
        syndrome="selfing_core",
        score_name="score",
    )
    assert result["species_key"].tolist() == ["sp 0", "sp 1"]
    assert result["score"].tolist() == [0.0, 1.0]


def test_exact_h2_bridge_recovers_negative_family_effects():
    config = {
        "analysis": {"publication_total_weight": 1.0},
        "families": {
            "reproductive_assurance": {
                "H2_syndrome": "selfing_core",
                "score_name": "selfing_core_score",
                "parent": "reproductive_assurance",
            },
            "accessibility_generalization": {
                "H2_syndrome": "generalized_accessible",
                "score_name": "generalized_accessible_score",
                "parent": "floral_architecture",
            },
        },
    }
    results, manifest = run_exact_h2_bridge(
        _scores(),
        _effects(),
        _effects(),
        config,
    )
    primary = results.loc[results["analysis"].eq("primary")]
    assert len(primary) == 2
    assert primary["evaluable"].astype(bool).all()
    assert primary["estimate"].lt(0).all()
    assert manifest["all_primary_estimates_negative"] is True
