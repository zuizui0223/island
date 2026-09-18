from pathlib import Path

import yaml


def test_v14_hypothesis_order_and_roles():
    config = yaml.safe_load(
        Path("config/chapter1_v14_hypothesis_architecture.yml").read_text(encoding="utf-8")
    )
    assert list(config["hypotheses"]) == [
        "H1_global_island_syndrome",
        "H2_floral_shift_decomposition",
        "H3_global_pollen_limitation",
        "H4_functional_bridge",
    ]
    h1 = config["hypotheses"]["H1_global_island_syndrome"]["primary_response_vector"]
    assert h1["colour_dulling"] == ["plain_colour"]
    assert len(h1["reproductive_assurance"]) == 3
    assert len(h1["accessibility_generalization"]) == 3
    assert config["claim_ceiling"]["conditional_H2_is_causal_mediation"] is False


def test_v14_H1_probability_contains_seven_atoms_and_three_domains():
    config = yaml.safe_load(
        Path("config/chapter1_v14_all_data_probability.yml").read_text(encoding="utf-8")
    )
    assert len(config["model_outcomes"]) == 7
    assert "plain_colour" in config["model_outcomes"]
    assert set(config["response_families"]) == {
        "reproductive_assurance",
        "colour_dulling",
        "accessibility_generalization",
    }
