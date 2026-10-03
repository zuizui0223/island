from pathlib import Path

import pandas as pd
import yaml

from island_v2.chapter1_all_data_probability import _stratum_mask


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


def test_introduced_stratum_uses_source_backed_origin_status_only():
    frame = pd.DataFrame(
        {
            "origin_status": [
                "native",
                "introduced",
                "unresolved",
                "introduced_or_uncertain",
            ],
            "floristic_status": [
                "native_nonendemic",
                "introduced",
                "unresolved",
                "introduced",
            ],
        }
    )
    assert _stratum_mask(frame, "all_introduced").tolist() == [
        False,
        True,
        False,
        False,
    ]


def test_origin_sensitivity_only_swaps_status_strata():
    parent = _load("config/chapter1_h1_final_directional.yml")
    sens = _load("config/chapter1_h1_origin_sensitivity.yml")
    assert parent["flora_roles"]["broad_primary"] == "all_observed"
    assert sens["flora_roles"]["broad_primary"] == "all_observed"
    assert sens["strata"] == [
        "all_observed",
        "all_native",
        "all_introduced",
    ]
