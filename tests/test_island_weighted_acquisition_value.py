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


def test_repository_config_is_loadable_and_declares_gates(tmp_path):
    from pathlib import Path

    config = value.load_config(Path("config/island_weighted_acquisition.yml"))

    assert config["readiness"]["min_covered_species"] > 0
    assert 0 < config["readiness"]["min_covered_fraction"] <= 1
    assert "Poaceae" in config["scope_classes"]["anemophilous_families"]
    # Plantaginaceae sensu APG is not uniformly anemophilous and must stay out
    assert "Plantaginaceae" not in config["scope_classes"]["anemophilous_families"]
    assert "Plantago" in config["scope_classes"]["anemophilous_genera"]
