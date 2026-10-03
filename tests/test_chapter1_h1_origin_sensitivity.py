from pathlib import Path

import yaml


def _load(path: str):
    return yaml.safe_load(Path(path).read_text(encoding="utf-8"))


def test_origin_sensitivity_keeps_frozen_seven_trait_h1_unchanged():
    parent = _load("config/chapter1_h1_final_directional.yml")
    sens = _load("config/chapter1_h1_origin_sensitivity.yml")

    assert sens["model_outcomes"] == parent["model_outcomes"]
    assert sens["response_families"] == parent["response_families"]
    assert sens["directional_score"]["outcome_weights"] == parent[
        "directional_score"
    ]["outcome_weights"]
    assert sens["broad_outcomes"] == parent["broad_outcomes"]
    assert sens["baseline_covariates"] == parent["baseline_covariates"]
    assert sens["contexts"] == parent["contexts"]
    assert sens["global_synthesis"] == parent["global_synthesis"]
    assert "flower_size_reduction" not in sens["model_outcomes"]


def test_origin_sensitivity_is_source_native_only():
    sens = _load("config/chapter1_h1_origin_sensitivity.yml")
    assert sens["flora_roles"]["broad_primary"] == "all_observed"
    assert sens["strata"] == ["all_observed", "all_native"]
    assert sens["flora_roles"]["status_sensitivities"] == ["all_native"]
    assert "all_introduced" not in sens["strata"]
