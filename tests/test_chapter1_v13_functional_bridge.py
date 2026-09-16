from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.chapter1_v13_functional_bridge import (
    aggregate_species_measurement_cells,
    fit_global_trait_level,
    fit_within_group_trait_level,
)

MEASUREMENT_COLUMNS = (
    "PL_Effect_Size_Type1",
    "PL_Effect_Size_Type2",
    "Constant_added",
    "Level_of_Supplementation",
)


def _raw_row(
    *,
    study: str,
    site: str,
    species: str,
    trait_state: int,
    effect: float,
    context: str = "northern_midlatitude",
    distance: float = 0.0,
) -> dict[str, object]:
    return {
        "study_key": study,
        "site_key": site,
        "species_key": species,
        "analysis_regime": context,
        "z_distance": distance,
        "trait": "autonomous_selfing",
        "trait_state": trait_state,
        "PL_Effect_Size": effect,
        "PL_Effect_Size_Type1": "logRR",
        "PL_Effect_Size_Type2": "Sup",
        "Constant_added": "False",
        "Constant_added_bool": False,
        "Level_of_Supplementation": "full",
    }


def test_aggregate_cells_gives_each_publication_total_weight_one() -> None:
    rows = pd.DataFrame(
        [
            _raw_row(study="p1", site="a", species="s1", trait_state=0, effect=1.0),
            _raw_row(study="p1", site="b", species="s2", trait_state=1, effect=0.5),
            _raw_row(study="p1", site="c", species="s3", trait_state=0, effect=1.1),
            _raw_row(study="p2", site="d", species="s4", trait_state=1, effect=0.4),
        ]
    )
    cells = aggregate_species_measurement_cells(rows, publication_total_weight=1.0)
    totals = cells.groupby("study_key")["analysis_weight"].sum()
    assert np.allclose(totals.to_numpy(), 1.0)


def _synthetic_cells(effect: float = -0.5) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    contexts = ["northern_midlatitude", "tropical"]
    for i in range(40):
        study = f"p{i:02d}"
        context = contexts[i % 2]
        distance = (i - 19.5) / 10.0
        study_offset = ((i % 7) - 3) * 0.015
        for state in (0, 1):
            state_noise = ((i + state) % 5 - 2) * 0.01
            rows.append(
                {
                    "study_key": study,
                    "site_key": f"site{i:02d}",
                    "species_key": f"sp{i:02d}_{state}",
                    "analysis_regime": context,
                    "z_distance": distance,
                    "trait": "autonomous_selfing",
                    "trait_state": state,
                    "PL_Effect_Size": (
                        1.0
                        + 0.12 * distance
                        + (0.15 if context == "tropical" else 0.0)
                        + effect * state
                        + study_offset
                        + state_noise
                    ),
                    "PL_Effect_Size_Type1": "logRR",
                    "PL_Effect_Size_Type2": "Sup" if i % 3 else "Full",
                    "Constant_added": "False",
                    "Constant_added_bool": False,
                    "Level_of_Supplementation": "full",
                    "analysis_weight": 0.5,
                }
            )
    return pd.DataFrame(rows)


def test_global_trait_level_recovers_negative_trait_effect() -> None:
    result = fit_global_trait_level(_synthetic_cells(effect=-0.5))
    assert result["evaluable"] is True
    assert result["trait_state_estimate"] < -0.4
    assert result["trait_state_one_sided_negative_p"] < 0.05


def test_within_publication_requires_both_trait_states() -> None:
    cells = _synthetic_cells().drop_duplicates("study_key").copy()
    result = fit_within_group_trait_level(cells, ["study_key"])
    assert result["evaluable"] is False
    assert result["reason"] == "no_groups_with_both_trait_states"


def _paired_site_cells() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for i in range(24):
        study = f"pair{i:02d}"
        site = f"site{i:02d}"
        local_effect = -0.35 - 0.03 * (i % 4)
        for state in (0, 1):
            rows.append(
                {
                    "study_key": study,
                    "site_key": site,
                    "species_key": f"pairsp{i:02d}_{state}",
                    "analysis_regime": "tropical" if i % 2 else "northern_midlatitude",
                    "z_distance": (i - 11.5) / 8.0,
                    "trait": "autonomous_selfing",
                    "trait_state": state,
                    "PL_Effect_Size": 1.2 + local_effect * state + 0.01 * (i % 3),
                    "PL_Effect_Size_Type1": "logRR",
                    "PL_Effect_Size_Type2": "Sup",
                    "Constant_added": "False",
                    "Constant_added_bool": False,
                    "Level_of_Supplementation": "full",
                    "analysis_weight": 0.5,
                }
            )
    return pd.DataFrame(rows)


def test_within_site_recovers_paired_negative_difference() -> None:
    result = fit_within_group_trait_level(
        _paired_site_cells(), ["study_key", "site_key"]
    )
    assert result["evaluable"] is True
    assert result["n_fixed_effect_groups"] == 24
    assert result["trait_state_estimate"] < -0.3
    assert result["trait_state_one_sided_negative_p"] < 0.05
