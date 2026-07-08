from pathlib import Path

import pandas as pd

from island_v2.trait_colour import (
    annotate_admissibility,
    load_schema,
    validate_colour,
)

SCHEMA = load_schema(Path("config/colour_schema.yml"))


def _valid_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "accepted_species": ["A a", "B b", "C c"],
            "raw_colour_description": ["white flowers", "red cultivar", "opens yellow, fades orange"],
            "primary_colour": ["white", "red", "yellow"],
            "secondary_colour": ["", "", "orange"],
            "within_species_variation": ["none", "none", "anthesis_colour_change"],
            "geographic_scope": ["species_range", "regional", "species_range"],
            "wild_or_cultivated": ["wild", "cultivated", "wild"],
            "colour_evidence_type": ["flora_or_description", "image", "field_observation"],
        }
    )


def test_valid_table_has_no_issues():
    assert validate_colour(_valid_table(), SCHEMA) == []


def test_missing_required_column_is_flagged():
    table = _valid_table().drop(columns=["primary_colour"])
    issues = validate_colour(table, SCHEMA)
    assert any("missing required columns" in i for i in issues)


def test_invalid_vocabulary_is_flagged():
    table = _valid_table()
    table.loc[0, "primary_colour"] = "chartreuse"
    table.loc[1, "within_species_variation"] = "sometimes"
    issues = validate_colour(table, SCHEMA)
    assert any("invalid primary_colour" in i for i in issues)
    assert any("invalid within_species_variation" in i for i in issues)


def test_cultivated_row_is_non_admissible():
    result, summary = annotate_admissibility(_valid_table(), SCHEMA)
    # row 1 is a cultivated plant -> not a species-level primary colour
    assert not result.loc[1, "species_level_primary_colour_admissible"]
    # row 0 is wild, no variation artefact -> admissible
    assert result.loc[0, "species_level_primary_colour_admissible"]


def test_anthesis_colour_change_row_is_non_admissible():
    result, _ = annotate_admissibility(_valid_table(), SCHEMA)
    # row 2 is an anthesis colour change -> variation artefact, not admissible
    assert not result.loc[2, "species_level_primary_colour_admissible"]


def test_annotate_retains_all_rows_and_counts():
    table = _valid_table()
    result, summary = annotate_admissibility(table, SCHEMA)
    assert len(result) == len(table)  # nothing dropped
    assert summary["n_rows"] == 3
    assert summary["n_species_level_admissible"] == 1
    assert summary["n_flagged_non_admissible"] == 2
