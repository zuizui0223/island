from __future__ import annotations

import yaml


def _config() -> dict:
    return yaml.safe_load(
        open("config/chapter1_nee_gbif_background_acquisition.yml", encoding="utf-8")
    )


def test_four_new_background_campaigns_are_frozen() -> None:
    config = _config()
    assert set(config["campaigns"]) == {
        "non_bombus_bees",
        "lepidoptera",
        "flower_visiting_birds",
        "diptera",
    }
    assert config["campaigns"]["non_bombus_bees"]["gbif_acquisition_taxon"] == "Hymenoptera"
    assert config["campaigns"]["non_bombus_bees"]["analytical_background_filter"] == "frozen_seven_bee_families"
    assert config["campaigns"]["lepidoptera"]["gbif_rank"] == "ORDER"
    assert config["campaigns"]["flower_visiting_birds"]["gbif_rank"] == "CLASS"
    assert config["campaigns"]["diptera"]["gbif_rank"] == "ORDER"


def test_bombus_remains_existing_apidae_route() -> None:
    config = _config()
    assert config["bombus"]["not_rerun_in_new_four_channel_campaign"] is True
    assert "Apidae" in config["bombus"]["acquisition_route"]


def test_spatial_assignment_is_frozen_exact_island() -> None:
    config = _config()
    assert config["spatial_scope"]["exact_island_assignment_required"] is True
    assert config["spatial_scope"]["acquisition_buffer_m"] == 2000
    assert config["spatial_scope"]["max_islands_per_block"] == 125
