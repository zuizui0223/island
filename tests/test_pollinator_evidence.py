from pathlib import Path

import pandas as pd

from island_v2.pollinator_evidence import (
    functional_replacement_index,
    load_schema,
    validate_long_format,
)

SCHEMA = load_schema(Path("config/pollinator_evidence_schema.yml"))


def _rows(records):
    cols = SCHEMA["required_columns"]
    return pd.DataFrame(records, columns=cols)


def _row(species, guild, evidence, effectiveness, review="accepted"):
    return {
        "accepted_species": species,
        "pollinator_guild": guild,
        "evidence_type": evidence,
        "effectiveness_class": effectiveness,
        "source_citation": "Author 2020",
        "source_url": "https://example.org/x",
        "study_location": "Island A",
        "season": "summer",
        "wild_or_cultivated": "wild",
        "review_status": review,
    }


def test_valid_long_format_passes():
    table = _rows([_row("A a", "bumblebees", "pollen_deposition", "confirmed_effective_pollinator")])
    assert validate_long_format(table, SCHEMA) == []


def test_invalid_vocabulary_is_reported():
    table = _rows([_row("A a", "bees", "seen_flying_near", "maybe")])
    issues = validate_long_format(table, SCHEMA)
    assert any("evidence_type" in i for i in issues)
    assert any("effectiveness_class" in i for i in issues)


def test_collapsed_guild_column_is_rejected():
    table = _rows([_row("A a", "bees", "floral_visitor", "visitor_only")])
    table["mixed_or_generalist"] = "yes"
    issues = validate_long_format(table, SCHEMA)
    assert any("collapse" in i or "forbidden" in i for i in issues)


def test_index_flags_functional_replacement_from_two_confirmed_guilds():
    table = _rows(
        [
            _row("Multi sp", "bumblebees", "exclusion_experiment", "confirmed_effective_pollinator"),
            _row("Multi sp", "hoverflies", "fruit_or_seed_set_contribution", "confirmed_effective_pollinator"),
            _row("Single sp", "bumblebees", "floral_visitor", "visitor_only"),
            _row("Ignored sp", "bees", "pollen_contact", "confirmed_effective_pollinator", review="unreviewed"),
        ]
    )
    index = functional_replacement_index(table, SCHEMA).set_index("accepted_species")
    assert bool(index.loc["Multi sp", "functional_replacement_candidate"]) is True
    assert int(index.loc["Multi sp", "n_confirmed_effective_guilds"]) == 2
    assert index.loc["Multi sp", "strongest_evidence_type"] == "exclusion_experiment"
    # A single visitor-only guild is not a replacement candidate.
    assert bool(index.loc["Single sp", "functional_replacement_candidate"]) is False
    # Unreviewed rows never enter the derived index.
    assert "Ignored sp" not in index.index


def test_index_never_collapses_raw_table():
    table = _rows(
        [
            _row("A a", "bumblebees", "floral_visitor", "visitor_only"),
            _row("A a", "moths", "pollen_deposition", "probable_pollinator"),
        ]
    )
    before = table.copy()
    functional_replacement_index(table, SCHEMA)
    # deriving the index does not mutate the raw long-format rows.
    pd.testing.assert_frame_equal(table, before)
