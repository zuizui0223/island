import pandas as pd
import pytest
import typer

from island_v2.core_pilot_nomination import nominate_core_pilots


def _applicability(rows):
    cols = ["island_id", "source_region_id", "applicability", "regional_analysis_only", "native_bombus_status"]
    return pd.DataFrame(rows, columns=cols)


def test_eligible_requires_applicable_core_plus_collected_flora():
    appl = _applicability(
        [
            ["isl_jp", "temperate_japan", "applicable", "false", "native_present"],
            ["isl_jp_sparse", "temperate_japan", "applicable", "false", "native_present"],
            ["isl_cape", "cape", "structurally_not_applicable", "false", "native_absent"],
            ["isl_conting", "montane", "applicable", "true", "native_present"],  # regional-only, not Core
        ]
    )
    effort = pd.DataFrame(
        {
            "island_id": ["isl_jp", "isl_jp_sparse", "isl_cape", "isl_conting"],
            "n_species": [40, 2, 200, 50],
        }
    )
    table, report = nominate_core_pilots(appl, effort, min_flora_species=5)
    by = table.set_index("island_id")

    assert bool(by.loc["isl_jp", "core_pilot_eligible"]) is True
    # applicable Core but <5 collected flora species -> not eligible, flagged as a gap
    assert bool(by.loc["isl_jp_sparse", "core_pilot_eligible"]) is False
    assert bool(by.loc["isl_jp_sparse", "is_core_applicable"]) is True
    # out-of-domain island with rich flora is never eligible
    assert bool(by.loc["isl_cape", "core_pilot_eligible"]) is False
    # regional-contingency (regional_analysis_only) is not Core
    assert bool(by.loc["isl_conting", "is_core_applicable"]) is False

    assert report["n_core_pilot_eligible"] == 1
    assert report["eligible_island_ids"] == ["isl_jp"]
    assert [g["island_id"] for g in report["applicable_but_no_collected_flora"]] == ["isl_jp_sparse"]


def test_no_eligible_island_reports_gap_without_inventing_one():
    # Only out-of-domain islands have flora (the real current state).
    appl = _applicability(
        [
            ["isl_cape", "cape", "structurally_not_applicable", "false", "native_absent"],
            ["isl_annobon", "gog", "structurally_not_applicable", "false", "native_absent"],
        ]
    )
    effort = pd.DataFrame({"island_id": ["isl_cape", "isl_annobon"], "n_species": [153, 128]})
    _table, report = nominate_core_pilots(appl, effort, min_flora_species=5)
    assert report["n_core_pilot_eligible"] == 0
    assert report["n_core_applicable_islands"] == 0
    assert report["eligible_island_ids"] == []


def test_prohibited_outcome_columns_are_rejected():
    appl = _applicability([["isl", "r", "applicable", "false", "native_present"]])
    appl["flower_primary_color"] = "white"  # outcome trait must not influence nomination
    effort = pd.DataFrame({"island_id": ["isl"], "n_species": [10]})
    with pytest.raises(typer.BadParameter):
        nominate_core_pilots(appl, effort, min_flora_species=5)
