import importlib
from pathlib import Path

import pandas as pd
import pytest
import typer

clade = importlib.import_module("island_v2.anemophilous_clade_rule")


DECLARED = {"family": {"Poaceae"}, "genus": {"Plantago"}}

LOOSE = clade.CladeThresholds(
    min_species=3,
    min_dominance=0.9,
    min_species_loo_accuracy=0.9,
    min_lineage_loo_accuracy=0.9,
    min_source_lineages=3,
)


def _evidence(records):
    return pd.DataFrame(
        records,
        columns=[
            "accepted_species",
            "family",
            "trait_name",
            "normalized_value",
            "evidence_scope",
            "source_lineage",
        ],
    )


INVARIANT = _evidence(
    [
        ("Poa alpha", "Poaceae", "flower_primary_color", "green", "species_direct", "doi:a"),
        ("Poa beta", "Poaceae", "flower_primary_color", "green", "species_direct", "doi:b"),
        ("Poa gamma", "Poaceae", "flower_primary_color", "green", "species_direct", "doi:c"),
        ("Poa delta", "Poaceae", "flower_primary_color", "green", "synonym_direct", "doi:d"),
    ]
)


def test_invariant_declared_clade_passes_every_gate():
    rules = clade.build_clade_rules(INVARIANT, DECLARED, LOOSE)

    assert len(rules) == 1
    rule = rules.iloc[0]
    assert rule["clade"] == "Poaceae"
    assert rule["clade_rank"] == "family"
    assert rule["inferred_value"] == "green"
    assert rule["n_direct_species"] == 4
    assert rule["dominance"] == pytest.approx(1.0)
    assert rule["species_loo_accuracy"] == pytest.approx(1.0)
    assert rule["lineage_loo_accuracy"] == pytest.approx(1.0)
    assert bool(rule["eligible"])
    assert rule["tier"] == clade.TIER


def test_undeclared_clade_is_never_scored():
    evidence = _evidence(
        [
            ("Rubia alpha", "Rubiaceae", "flower_primary_color", "white", "species_direct", "d:a"),
            ("Rubia beta", "Rubiaceae", "flower_primary_color", "white", "species_direct", "d:b"),
            ("Rubia gamma", "Rubiaceae", "flower_primary_color", "white", "species_direct", "d:c"),
        ]
    )
    assert clade.build_clade_rules(evidence, DECLARED, LOOSE).empty


def test_a_real_counterexample_fails_the_rule_closed():
    evidence = pd.concat(
        [
            INVARIANT,
            _evidence(
                [
                    (
                        "Poa epsilon",
                        "Poaceae",
                        "flower_primary_color",
                        "purple",
                        "species_direct",
                        "doi:e",
                    )
                ]
            ),
        ],
        ignore_index=True,
    )
    rule = clade.build_clade_rules(evidence, DECLARED, LOOSE).iloc[0]

    assert rule["counterexample_species"] == 1
    assert rule["dominance"] == pytest.approx(0.8)
    assert not bool(rule["eligible"])
    assert "dominance" in rule["ineligible_reason"]


def test_single_source_lineage_cannot_satisfy_support_alone():
    evidence = _evidence(
        [
            ("Poa alpha", "Poaceae", "flower_primary_color", "green", "species_direct", "doi:only"),
            ("Poa beta", "Poaceae", "flower_primary_color", "green", "species_direct", "doi:only"),
            ("Poa gamma", "Poaceae", "flower_primary_color", "green", "species_direct", "doi:only"),
        ]
    )
    rule = clade.build_clade_rules(evidence, DECLARED, LOOSE).iloc[0]

    assert rule["n_source_lineages"] == 1
    assert not bool(rule["eligible"])
    assert "lineages" in rule["ineligible_reason"]


def test_low_support_fails_even_when_perfectly_dominant():
    evidence = INVARIANT.iloc[:2]
    rule = clade.build_clade_rules(evidence, DECLARED, LOOSE).iloc[0]

    assert rule["dominance"] == pytest.approx(1.0)
    assert not bool(rule["eligible"])
    assert "support" in rule["ineligible_reason"]


def test_non_direct_scopes_are_excluded_from_support():
    evidence = _evidence(
        [
            ("Poa alpha", "Poaceae", "flower_primary_color", "green", "species_direct", "doi:a"),
            ("Poa beta", "Poaceae", "flower_primary_color", "green", "validated_low", "doi:b"),
            ("Poa gamma", "Poaceae", "flower_primary_color", "green", "family_inference", "doi:c"),
            ("Poa delta", "Poaceae", "flower_primary_color", "green", "global_fallback", "doi:d"),
        ]
    )
    rule = clade.build_clade_rules(evidence, DECLARED, LOOSE).iloc[0]

    # only the one species_direct row counts toward support
    assert rule["n_direct_species"] == 1
    assert not bool(rule["eligible"])


def test_species_with_self_contradictory_evidence_is_dropped():
    evidence = pd.concat(
        [
            INVARIANT,
            _evidence(
                [
                    (
                        "Poa zeta",
                        "Poaceae",
                        "flower_primary_color",
                        "green",
                        "species_direct",
                        "doi:z",
                    ),
                    (
                        "Poa zeta",
                        "Poaceae",
                        "flower_primary_color",
                        "purple",
                        "species_direct",
                        "doi:z2",
                    ),
                ]
            ),
        ],
        ignore_index=True,
    )
    rule = clade.build_clade_rules(evidence, DECLARED, LOOSE).iloc[0]

    # Poa zeta disagrees with itself and is neither counted nor treated as a
    # counterexample
    assert rule["n_direct_species"] == 4
    assert rule["counterexample_species"] == 0


def test_declared_genus_rule_is_scored_at_genus_rank():
    evidence = _evidence(
        [
            ("Plantago a", "Plantaginaceae", "floral_form", "reduced", "species_direct", "d:a"),
            ("Plantago b", "Plantaginaceae", "floral_form", "reduced", "species_direct", "d:b"),
            ("Plantago c", "Plantaginaceae", "floral_form", "reduced", "species_direct", "d:c"),
            ("Veronica a", "Plantaginaceae", "floral_form", "zygomorphic", "species_direct", "d:d"),
        ]
    )
    rules = clade.build_clade_rules(evidence, DECLARED, LOOSE)

    assert len(rules) == 1
    rule = rules.iloc[0]
    assert rule["clade"] == "Plantago"
    assert rule["clade_rank"] == "genus"
    assert bool(rule["eligible"])
    # Veronica is in the same family but is not a declared clade, so it neither
    # contributes to nor breaks the Plantago rule
    assert rule["n_direct_species"] == 3


def test_apply_emits_tier_rows_only_for_eligible_rules_and_never_overwrites_direct():
    rules = clade.build_clade_rules(INVARIANT, DECLARED, LOOSE)
    unresolved = pd.DataFrame(
        {
            "accepted_species": ["Poa omega", "Poa alpha", "Rubia alpha"],
            "family": ["Poaceae", "Poaceae", "Rubiaceae"],
            "trait_name": ["flower_primary_color"] * 3,
        }
    )
    applied = clade.apply_clade_rules(rules, unresolved, INVARIANT)

    assert applied["accepted_species"].tolist() == ["Poa omega"]
    assert applied.iloc[0]["inferred_value"] == "green"
    assert applied.iloc[0]["tier"] == clade.TIER


def test_apply_emits_nothing_when_no_rule_is_eligible():
    rules = clade.build_clade_rules(INVARIANT.iloc[:2], DECLARED, LOOSE)
    unresolved = pd.DataFrame(
        {
            "accepted_species": ["Poa omega"],
            "family": ["Poaceae"],
            "trait_name": ["flower_primary_color"],
        }
    )
    assert clade.apply_clade_rules(rules, unresolved, INVARIANT).empty


def test_default_thresholds_are_stricter_than_the_genus_track():
    validated_low = importlib.import_module("island_v2.validated_low_inference")
    genus = validated_low.DEFAULT_THRESHOLDS

    assert clade.DEFAULT_THRESHOLDS.min_species > genus.min_species
    assert clade.DEFAULT_THRESHOLDS.min_dominance >= genus.colour_min_dominance


def test_malformed_evidence_fails_closed():
    with pytest.raises(typer.BadParameter):
        clade.build_clade_rules(INVARIANT.drop(columns=["normalized_value"]), DECLARED, LOOSE)


def test_repository_config_declares_the_clades():
    declared = clade.load_declared_clades(Path("config/island_weighted_acquisition.yml"))

    assert "Poaceae" in declared["family"]
    assert "Cyperaceae" in declared["family"]
    assert "Plantago" in declared["genus"]
    # Plantaginaceae sensu APG is not uniformly anemophilous
    assert "Plantaginaceae" not in declared["family"]
