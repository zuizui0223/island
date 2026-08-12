import importlib

import pandas as pd
import pytest
import typer

value = importlib.import_module("island_v2.island_weighted_acquisition_value")


SCOPE_CLASSES = {
    "anemophilous_families": ["Poaceae"],
    "anemophilous_genera": ["Plantago"],
    "non_angiosperm_families": ["Equisetaceae"],
}

MASTER = pd.DataFrame(
    {
        "accepted_species": [
            "Wide grassa",
            "Plantago testa",
            "Equisetum testa",
            "Narrow enda",
            "Narrow endb",
        ],
        "genus": ["Wide", "Plantago", "Equisetum", "Narrow", "Narrow"],
        "family": ["Poaceae", "Plantaginaceae", "Equisetaceae", "Rubiaceae", "Rubiaceae"],
        "n_islands": [4, 3, 2, 1, 1],
    }
)

ISLAND_SPECIES = pd.DataFrame(
    {
        "island_id": [
            "i1", "i1", "i1",
            "i2", "i2",
            "i3", "i3",
            "i4",
        ],
        "species": [
            "Wide grassa", "Plantago testa", "Narrow enda",
            "Wide grassa", "Plantago testa",
            "Wide grassa", "Equisetum testa",
            "Wide grassa",
        ],
    }
)


def test_species_value_ranks_by_island_footprint_and_shares_sum_to_one():
    frame = value.build_species_value(MASTER, SCOPE_CLASSES)

    assert frame["accepted_species"].tolist()[0] == "Wide grassa"
    assert frame["priority_rank"].tolist() == [1, 2, 3, 4, 5]
    assert frame["flora_mass_share"].sum() == pytest.approx(1.0)
    assert frame["cumulative_flora_mass_share"].iloc[-1] == pytest.approx(1.0)
    # ties are broken by name so the ranking is stable across runs
    assert frame.loc[frame["n_islands"] == 1, "accepted_species"].tolist() == [
        "Narrow enda",
        "Narrow endb",
    ]


def test_scope_class_separates_planning_lanes_without_asserting_trait_values():
    frame = value.build_species_value(MASTER, SCOPE_CLASSES).set_index("accepted_species")

    assert frame.loc["Wide grassa", "scope_class"] == value.ANEMOPHILOUS
    assert frame.loc["Plantago testa", "scope_class"] == value.ANEMOPHILOUS
    assert frame.loc["Equisetum testa", "scope_class"] == value.NON_ANGIOSPERM
    assert frame.loc["Narrow enda", "scope_class"] == value.BIOTIC


def test_plantaginaceae_outside_plantago_stays_biotic():
    master = pd.DataFrame(
        {
            "accepted_species": ["Veronica testa"],
            "genus": ["Veronica"],
            "family": ["Plantaginaceae"],
            "n_islands": [5],
        }
    )
    frame = value.build_species_value(master, SCOPE_CLASSES)
    assert frame.loc[0, "scope_class"] == value.BIOTIC


def test_island_readiness_applies_both_gates():
    readiness = value.island_readiness(
        ISLAND_SPECIES, {"Wide grassa"}, min_covered_species=1, min_covered_fraction=0.6
    ).set_index("island_id")

    assert readiness.loc["i1", "n_flora_species"] == 3
    assert readiness.loc["i1", "covered_fraction"] == pytest.approx(1 / 3)
    assert not bool(readiness.loc["i1", "analysis_ready"])
    assert bool(readiness.loc["i4", "analysis_ready"])

    # the count gate alone must not pass an unrepresentative fragment
    strict = value.island_readiness(
        ISLAND_SPECIES, {"Wide grassa"}, min_covered_species=2, min_covered_fraction=0.0
    ).set_index("island_id")
    assert not bool(strict.loc["i4", "analysis_ready"])


def test_frontier_prefers_island_weighted_targeting_at_equal_budget():
    species_value = value.build_species_value(MASTER, SCOPE_CLASSES)
    frontier, summary = build = value.build_frontier(
        species_value,
        ISLAND_SPECIES,
        targets=[1],
        min_covered_species=1,
        min_covered_fraction=0.5,
    )
    assert build is not None

    weighted = frontier.loc[frontier["strategy"] == "island_weighted"].iloc[0]
    plain = frontier.loc[frontier["strategy"] == "unweighted_species_count"].iloc[0]

    # unweighted name order reaches "Equisetum testa" first, the widest-ranging
    # species is reached only by island-weighted targeting
    assert weighted["flora_mass_share"] > plain["flora_mass_share"]
    assert weighted["n_islands_analysis_ready"] >= plain["n_islands_analysis_ready"]
    assert summary["total_flora_mass"] == 11
    assert summary["n_islands"] == 4


def test_scope_breakdown_splits_a_budget_by_evidence_lane():
    species_value = value.build_species_value(MASTER, SCOPE_CLASSES)
    breakdown = value.scope_breakdown(species_value, budget=3)

    assert breakdown["n_species"].sum() == 3
    assert set(breakdown["target_budget"]) == {3}
    assert set(breakdown["scope_class"]) == {value.ANEMOPHILOUS, value.NON_ANGIOSPERM}


def test_missing_columns_fail_closed():
    with pytest.raises(typer.BadParameter):
        value.build_species_value(MASTER.drop(columns=["n_islands"]), SCOPE_CLASSES)
    with pytest.raises(typer.BadParameter):
        value.island_readiness(ISLAND_SPECIES.drop(columns=["species"]), set(), 1, 0.5)
    with pytest.raises(typer.BadParameter):
        value.build_species_value(MASTER.assign(n_islands=0), SCOPE_CLASSES)


AXES = {
    "flower_colour": ["flower_primary_color"],
    "reproductive_assurance": ["self_incompatibility"],
}

TIERS = {
    "direct_only": {"include_validated_low": False, "role": "primary"},
    "direct_plus_validated_low": {"include_validated_low": True, "role": "sensitivity"},
}

EVIDENCE = pd.DataFrame(
    {
        "accepted_species": [
            "Wide grassa",
            "Plantago testa",
            "Narrow enda",
            "Equisetum testa",
            "Wide grassa",
        ],
        "trait_name": [
            "flower_primary_color",
            "flower_primary_color",
            "flower_primary_color",
            "flower_primary_color",
            "self_incompatibility",
        ],
        "evidence_scope": [
            "species_direct",
            "validated_low",
            "family_inference",
            "global_fallback",
            "synonym_direct",
        ],
    }
)


def test_direct_tier_excludes_validated_low_and_prohibited_scopes():
    direct = value.covered_species_by_axis(EVIDENCE, AXES, include_validated_low=False)

    assert direct["flower_colour"] == {"Wide grassa"}
    assert direct["reproductive_assurance"] == {"Wide grassa"}
    # family inference and global fallback are never admitted by either tier
    assert "Narrow enda" not in direct["flower_colour"]
    assert "Equisetum testa" not in direct["flower_colour"]


def test_sensitivity_tier_admits_validated_low_only():
    sensitivity = value.covered_species_by_axis(EVIDENCE, AXES, include_validated_low=True)

    assert sensitivity["flower_colour"] == {"Wide grassa", "Plantago testa"}
    assert "Narrow enda" not in sensitivity["flower_colour"]
    assert "Equisetum testa" not in sensitivity["flower_colour"]


def test_traits_outside_the_declared_axes_are_ignored():
    evidence = pd.DataFrame(
        {
            "accepted_species": ["Wide grassa", "Wide grassa"],
            "trait_name": ["reward_type", "pollen_vector_mode"],
            "evidence_scope": ["species_direct", "species_direct"],
        }
    )
    assert value.covered_species_by_axis(evidence, AXES, include_validated_low=False) == {}


def test_evaluate_reports_primary_and_sensitivity_separately():
    table, summary = value.evaluate_tiers(
        EVIDENCE,
        ISLAND_SPECIES,
        AXES,
        TIERS,
        min_covered_species=1,
        min_covered_fraction=0.5,
    )

    assert set(table["tier"]) == {"direct_only", "direct_plus_validated_low"}
    assert set(table["axis"]) == {"flower_colour", "reproductive_assurance"}

    primary = table.loc[
        (table["tier"] == "direct_only") & (table["axis"] == "flower_colour")
    ].iloc[0]
    sensitivity = table.loc[
        (table["tier"] == "direct_plus_validated_low") & (table["axis"] == "flower_colour")
    ].iloc[0]

    assert primary["role"] == "primary"
    assert sensitivity["role"] == "sensitivity"
    # the genus-inference tier can only ever add species, never remove them
    assert sensitivity["n_covered_species"] > primary["n_covered_species"]
    assert summary["primary_metric"].startswith("n_islands_analysis_ready")
    assert summary["sensitivity_min_axis_ready"] >= summary["primary_min_axis_ready"]


def test_evaluate_fails_closed_on_a_malformed_ledger():
    with pytest.raises(typer.BadParameter):
        value.covered_species_by_axis(
            EVIDENCE.drop(columns=["evidence_scope"]), AXES, include_validated_low=False
        )


def test_source_yield_scores_a_source_by_islands_not_species_count():
    species_value = value.build_species_value(MASTER, SCOPE_CLASSES)

    # a source that only supplies the narrowest endemic: many species, no islands
    narrow = value.source_yield(
        baseline=set(),
        candidate={"Narrow enda", "Narrow endb"},
        island_species=ISLAND_SPECIES,
        species_value=species_value,
        budget=5,
        min_covered_species=1,
        min_covered_fraction=0.6,
    )
    # a source that only supplies the widest species: fewer species, more islands
    wide = value.source_yield(
        baseline=set(),
        candidate={"Wide grassa"},
        island_species=ISLAND_SPECIES,
        species_value=species_value,
        budget=5,
        min_covered_species=1,
        min_covered_fraction=0.6,
    )

    assert narrow["n_new_species"] > wide["n_new_species"]
    assert wide["n_islands_analysis_ready_gained"] > narrow["n_islands_analysis_ready_gained"]
    assert wide["new_flora_mass_share"] > narrow["new_flora_mass_share"]


def test_source_yield_excludes_species_already_in_the_baseline():
    species_value = value.build_species_value(MASTER, SCOPE_CLASSES)
    summary = value.source_yield(
        baseline={"Wide grassa"},
        candidate={"Wide grassa", "Plantago testa"},
        island_species=ISLAND_SPECIES,
        species_value=species_value,
        budget=5,
        min_covered_species=1,
        min_covered_fraction=0.5,
    )

    assert summary["n_candidate_species"] == 2
    assert summary["n_new_species"] == 1
    assert summary["n_islands_analysis_ready_lost"] == 0


def test_source_yield_reports_priority_head_hit_rate():
    species_value = value.build_species_value(MASTER, SCOPE_CLASSES)
    summary = value.source_yield(
        baseline=set(),
        candidate={"Wide grassa", "Narrow endb"},
        island_species=ISLAND_SPECIES,
        species_value=species_value,
        budget=1,
        min_covered_species=1,
        min_covered_fraction=0.5,
    )

    # only Wide grassa is inside the top-1 priority head
    assert summary["n_new_species_in_priority_head"] == 1
    assert summary["priority_head_hit_rate"] == pytest.approx(0.5)


def test_repository_config_declares_axes_and_tiers():
    from pathlib import Path

    config = value.load_config(Path("config/island_weighted_acquisition.yml"))

    assert set(config["axes"]) == {
        "flower_colour",
        "floral_structural_complexity",
        "reproductive_assurance",
    }
    # reward_type and pollen_vector_mode are separate outcomes and must not fill
    # the floral-structure axis
    structure = config["axes"]["floral_structural_complexity"]
    assert "reward_type" not in structure
    assert "pollen_vector_mode" not in structure

    tiers = config["evaluation_tiers"]
    assert tiers["direct_only"]["role"] == "primary"
    assert tiers["direct_only"]["include_validated_low"] is False
    assert tiers["direct_plus_validated_low"]["role"] == "sensitivity"


def test_repository_config_is_loadable_and_declares_gates(tmp_path):
    from pathlib import Path

    config = value.load_config(Path("config/island_weighted_acquisition.yml"))

    assert config["readiness"]["min_covered_species"] > 0
    assert 0 < config["readiness"]["min_covered_fraction"] <= 1
    assert "Poaceae" in config["scope_classes"]["anemophilous_families"]
    # Plantaginaceae sensu APG is not uniformly anemophilous and must stay out
    assert "Plantaginaceae" not in config["scope_classes"]["anemophilous_families"]
    assert "Plantago" in config["scope_classes"]["anemophilous_genera"]
