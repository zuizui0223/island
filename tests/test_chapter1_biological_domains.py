from pathlib import Path

import pandas as pd
import yaml

from island_v2.chapter1_biological_domains import (
    aggregate_common_species_scores,
    audit_h3,
    domain_specs,
)


def _config():
    return yaml.safe_load(
        Path("config/chapter1_biological_domains.yml").read_text(
            encoding="utf-8"
        )
    )


def _syndrome():
    return yaml.safe_load(
        Path("config/chapter1_pr138_pollination_syndromes.yml").read_text(
            encoding="utf-8"
        )
    )


def test_domains_reuse_frozen_assurance_and_accessibility_definitions():
    specs = domain_specs(_syndrome(), _config())
    assert specs["reproductive_assurance"] == _syndrome()["syndromes"][
        "selfing_core"
    ]
    assert specs["accessibility_specialization"] == _syndrome()["syndromes"][
        "generalized_accessible"
    ]


def test_signal_display_direction_only_reuses_frozen_flower_size_recode():
    specs = domain_specs(_syndrome(), _config())
    size = specs["reduced_flower_size"]["traits"]
    assert set(size) == {"flower_size_class"}
    assert set(size["flower_size_class"]["preferred"]) == {
        "very_small",
        "small",
    }
    assert set(size["flower_size_class"]["opposed"]) == {
        "large",
        "very_large",
    }
    text = str(specs).casefold()
    assert "flower_primary_color" not in text
    assert "inflorescence_display" not in text


def test_common_species_support_uses_same_species_for_all_three_domains():
    flora = pd.DataFrame(
        {
            "island_id": ["i1", "i1", "i1", "i2", "i2"],
            "accepted_species": ["a", "b", "c", "a", "d"],
        }
    )
    scores = pd.DataFrame(
        {
            "accepted_species": ["a", "b", "c", "d"],
            "reproductive_assurance": [1.0, 1.0, 1.0, 1.0],
            "accessibility_specialization": [0.5, 0.5, None, 0.5],
            "reduced_flower_size": [0.2, None, 0.2, 0.2],
        }
    )
    out = aggregate_common_species_scores(
        flora,
        scores,
        minimum_species=1,
    )
    assert set(out["island_id"]) == {"i1", "i2"}
    assert out.set_index("island_id").loc["i1", "common_n_species"] == 1
    assert out.set_index("island_id").loc["i2", "common_n_species"] == 2


def test_h3_is_replayed_unchanged_with_finite_publication_reference():
    main = {
        "corrected": {
            "global_gradient": {
                "distance_slope": 0.1,
                "distance_slope_se": 0.04,
                "n_publications": 100,
            },
            "sensitivities": {
                "supplemental_only": {
                    "distance_slope": 0.03,
                    "distance_slope_se": 0.05,
                    "n_publications": 60,
                },
                "no_zero_constant": {
                    "distance_slope": 0.09,
                    "distance_slope_se": 0.04,
                    "n_publications": 95,
                },
            },
        }
    }
    offshore = {
        "offshore_continuous_gradient": {
            "estimate": 0.2,
            "se": 0.09,
            "n_publications": 80,
        }
    }
    out = audit_h3(main, offshore)
    primary = out.loc[out["analysis"].eq("corrected_global")].iloc[0]
    assert primary["estimate"] == 0.1
    assert primary["finite_publication_p_two_sided"] < 0.05
    assert len(out) == 4


def test_contract_marks_reclassification_posthoc_and_blocks_hard_attraction_order():
    config = _config()
    assert (
        config["inferential_role"]
        == "posthoc_biological_domain_reclassification_sensitivity"
    )
    guards = set(config["reviewer_guards"])
    assert "colour_hue_is_not_attraction_intensity" in guards
    assert "inflorescence_type_is_not_ordinal_display_investment" in guards
    assert "this_reclassification_is_posthoc_and_cannot_replace_confirmatory_H1" in guards
