import copy

import pandas as pd

from island_v2.chapter1_wcvp_final_directional_h1 import (
    STRATUM,
    build_replay_config,
    relabel_stratum,
)


def test_replay_config_only_changes_floristic_sampling_frame():
    parent = {
        "strata": ["all_observed", "all_native", "native_nonendemic"],
        "flora_roles": {
            "broad_primary": "all_observed",
            "status_sensitivities": ["all_native", "native_nonendemic"],
        },
        "model_outcomes": ["a", "b"],
        "response_families": {"x": ["a"], "y": ["b"]},
        "directional_score": {"outcome_weights": {"a": 0.5, "b": 0.5}},
        "contexts": ["north", "south"],
        "baseline_covariates": ["area"],
    }
    frozen = copy.deepcopy(parent)
    replay = build_replay_config(parent)

    assert parent == frozen
    assert replay["strata"] == ["all_native"]
    assert replay["flora_roles"] == {
        "broad_primary": "all_native",
        "status_sensitivities": [],
    }
    for key in (
        "model_outcomes",
        "response_families",
        "directional_score",
        "contexts",
        "baseline_covariates",
    ):
        assert replay[key] == frozen[key]


def test_relabel_only_changes_output_label():
    frame = pd.DataFrame(
        {
            "stratum": ["all_native", "other"],
            "estimate": [0.1, 0.2],
        }
    )
    out = relabel_stratum(frame)
    assert out["stratum"].tolist() == [STRATUM, "other"]
    assert out["estimate"].tolist() == [0.1, 0.2]
    assert frame["stratum"].tolist() == ["all_native", "other"]
