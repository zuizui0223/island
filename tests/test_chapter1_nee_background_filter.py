from __future__ import annotations

import pandas as pd
import pytest
import yaml

from island_v2.chapter1_nee_background_filter import filter_frozen_background


def _policy() -> dict:
    return yaml.safe_load(
        open("config/chapter1_nee_channel_observation_policy.yml", encoding="utf-8")
    )


def test_non_bombus_bee_background_keeps_only_frozen_seven_families() -> None:
    records = pd.DataFrame(
        {
            "family": ["Apidae", "Halictidae", "Vespidae", "Crabronidae", "Megachilidae"],
            "species": ["a", "b", "c", "d", "e"],
        }
    )
    result, receipt = filter_frozen_background(records, "non_bombus_bees", _policy())
    assert set(result["family"]) == {"Apidae", "Halictidae", "Megachilidae"}
    assert receipt["n_input_exact_island_records"] == 5
    assert receipt["n_frozen_background_records"] == 3
    assert receipt["n_records_excluded_by_background_definition"] == 2
    assert receipt["uses_functional_target_catalog"] is False


def test_bombus_background_is_apidae() -> None:
    records = pd.DataFrame(
        {"family": ["Apidae", "Halictidae"], "species": ["Bombus x", "Other bee"]}
    )
    result, receipt = filter_frozen_background(records, "bombus", _policy())
    assert list(result["family"]) == ["Apidae"]
    assert receipt["filter_rule"] == "family_exact_Apidae"


def test_exact_broad_taxon_channels_pass_through_without_functional_filtering() -> None:
    records = pd.DataFrame({"species": ["x", "y"], "family": ["f1", "f2"]})
    for channel in ["lepidoptera", "flower_visiting_birds", "diptera"]:
        result, receipt = filter_frozen_background(records, channel, _policy())
        assert len(result) == 2
        assert receipt["filter_rule"] == "exact_acquisition_taxon_is_frozen_background"
        assert receipt["uses_functional_target_catalog"] is False


def test_non_bombus_filter_fails_without_family_column() -> None:
    with pytest.raises(ValueError, match="requires family"):
        filter_frozen_background(pd.DataFrame({"species": ["x"]}), "non_bombus_bees", _policy())


def test_unregistered_channel_fails_closed() -> None:
    with pytest.raises(ValueError, match="unregistered"):
        filter_frozen_background(pd.DataFrame(), "unknown", _policy())
