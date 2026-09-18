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
    h4_result_heading = (
        "### H4: H2 response scores are associated with lower current pollen limitation"
    )
    assert text.count(h4_result_heading) == 1
    assert "### Supplementary prospective validation audit" in text
    assert "## Claim ceiling for v14" in text
    assert "## Claim ceiling for v13" not in text
    assert "seven-response" in text
    assert "plain_colour" in text


def test_v14_H2_uses_raw_trait_concordance_not_weighted_guild_scores():
    config = yaml.safe_load(
        Path("config/chapter1_v14_h2_decomposition.yml").read_text(encoding="utf-8")
    )
    raw = config["raw_pollination_syndrome_concordance"]
    assert raw["role"] == "primary_pollination_syndrome_concordance_layer"
    assert raw["named_weighted_scores_role"].startswith("historical_secondary")
    assert "large_bee_like" not in config["continuous_models"]
    assert "raw_architecture_given_raw_colour_conditional_on_selfing_core" in raw["estimands"]


def test_v14_H4_keeps_validation_audit_outside_result():
    config = yaml.safe_load(
        Path("config/chapter1_v14_h4_evidence_hierarchy.yml").read_text(encoding="utf-8")
    )
    layers = config["evidence_layers"]
    assert set(layers) == {
        "v13_exact_species_discovery",
        "v14_exact_H2_score_discovery",
        "v14_atomic_reconstruction_sensitivity",
    }
    exact = layers["v14_exact_H2_score_discovery"]
    assert exact["reproductive_assurance"]["two_sided_p"] < 0.05
    assert exact["accessibility_generalization"]["two_sided_p"] < 0.05

    audit = config["validation_audit"]
    assert audit["part_of_H4_result"] is False
    assert audit["role"] == "supplementary_nonpromoted_validation_audit"
    assert audit["wild_temporal"]["outcomes_unblinded"] is False
    assert audit["pollimcrop"]["outcomes_unblinded"] is False
    assert audit["pollimcrop"]["admitted_hypotheses"] == []

def test_v14_H4_exact_bridge_uses_literal_H2_species_scores():
    config = yaml.safe_load(
        Path("config/chapter1_v14_h4_exact_h2_score_bridge.yml").read_text(
            encoding="utf-8"
        )
    )
    assert config["H2_score_source"]["score_column"] == "soft_membership"
    assert (
        config["families"]["reproductive_assurance"]["H2_syndrome"]
        == "selfing_core"
    )
    assert (
        config["families"]["accessibility_generalization"]["H2_syndrome"]
        == "generalized_accessible"
    )


def test_v14_H4_atomic_reconstruction_is_sensitivity_only():
    config = yaml.safe_load(
        Path("config/chapter1_v14_h4_family_bridge.yml").read_text(encoding="utf-8")
    )
    assert (
        config["inferential_role"]
        == "posthoc_atomic_family_reconstruction_sensitivity"
    )
    access = config["families"]["accessibility_generalization"]
    assert access["role"] == "atomic_reconstruction_sensitivity"
    assert access["weights"] == {
        "generalized_form": 1.0,
        "actinomorphic_symmetry": 0.75,
        "shallow_open_tube": 1.0,
    }
