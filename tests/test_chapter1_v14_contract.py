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


def test_v14_manuscript_has_reordered_hypothesis_sections():
    text = Path(
        "docs/chapter1_manuscript_full_v14_reordered_hypotheses_20260918.md"
    ).read_text(encoding="utf-8")
    methods = [
        text.index("### H1: recurrent global floral/reproductive island syndrome"),
        text.index("### H2: selfing-syndrome versus pollinator-facing floral decomposition"),
        text.index("### H3: global experimental pollen limitation"),
        text.index("### H4: post-hoc functional triangulation"),
    ]
    assert methods == sorted(methods)
    assert text.count("### H4: post-hoc functional triangulation") == 1
    assert text.count(
        "### H4: current pollen limitation is lower in key island-syndrome trait states"
    ) == 1
    assert "## Claim ceiling for v14" in text
    assert "## Claim ceiling for v13" not in text
    assert "seven-response" in text
    assert "plain_colour" in text
