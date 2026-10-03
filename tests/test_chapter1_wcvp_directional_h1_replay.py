from pathlib import Path

import yaml

from island_v2.chapter1_h1_final_directional import validate_score_weights
from island_v2.chapter1_wcvp_directional_h1_replay import (
    regional_native_h1_config,
)


def test_wcvp_replay_changes_only_floristic_stratum():
    parent = yaml.safe_load(
        Path("config/chapter1_h1_final_directional.yml").read_text(
            encoding="utf-8"
        )
    )
    replay = regional_native_h1_config(parent)

    assert replay["model_outcomes"] == parent["model_outcomes"]
    assert replay["response_families"] == parent["response_families"]
    assert replay["directional_score"] == parent["directional_score"]
    assert replay["broad_outcomes"] == parent["broad_outcomes"]
    assert replay["baseline_covariates"] == parent["baseline_covariates"]
    assert replay["contexts"] == parent["contexts"]
    assert replay["global_synthesis"] == parent["global_synthesis"]

    assert replay["strata"] == ["all_native"]
    assert replay["flora_roles"]["broad_primary"] == "all_native"
    assert replay["flora_roles"]["status_sensitivities"] == []


def test_wcvp_replay_keeps_exact_frozen_seven_trait_weights():
    parent = yaml.safe_load(
        Path("config/chapter1_h1_final_directional.yml").read_text(
            encoding="utf-8"
        )
    )
    replay = regional_native_h1_config(parent)
    weights = validate_score_weights(replay)

    assert set(weights) == {
        "self_compatibility",
        "selfing_mating_system",
        "autonomous_selfing",
        "plain_colour",
        "generalized_form",
        "actinomorphic_symmetry",
        "shallow_open_tube",
    }
    assert "flower_size_reduction" not in weights
