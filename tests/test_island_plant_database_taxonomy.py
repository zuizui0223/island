from __future__ import annotations

import pandas as pd

from island_v2.island_plant_database_taxonomy import (
    COL_XR_CHECKLIST_KEY,
    deterministic_pilot,
    normalize_match,
    query_for_row,
    summarize,
)


def _row(name: str = "Alpha one") -> pd.Series:
    return pd.Series(
        {
            "accepted_species": name,
            "genus": name.split()[0],
            "family": "Alphaaceae",
        }
    )


def _classification() -> list[dict[str, str]]:
    return [
        {"key": "P", "name": "Plantae", "rank": "KINGDOM"},
        {"key": "F", "name": "Alphaaceae", "rank": "FAMILY"},
        {"key": "G", "name": "Alpha", "rank": "GENUS"},
        {"key": "S", "name": "Alpha one", "rank": "SPECIES"},
    ]


def test_query_is_species_strict_and_preserves_context() -> None:
    query = query_for_row(_row())
    assert query == {
        "scientificName": "Alpha one",
        "taxonRank": "SPECIES",
        "kingdom": "Plantae",
        "strict": True,
        "verbose": False,
        "genus": "Alpha",
        "family": "Alphaaceae",
    }


def test_exact_accepted_species_is_only_candidate_not_promoted() -> None:
    result = {
        "usage": {
            "key": "COL1",
            "name": "Alpha one Author",
            "canonicalName": "Alpha one",
            "authorship": "Author",
            "rank": "SPECIES",
            "status": "ACCEPTED",
        },
        "classification": _classification(),
        "diagnostics": {"matchType": "EXACT", "confidence": 99},
        "synonym": False,
    }
    row = normalize_match(_row(), result)
    assert row["automatic_resolution_candidate"] == "true"
    assert row["resolution_status"] == "exact_accepted_candidate"
    assert row["review_status"] == "automatic_candidate_not_promoted"
    assert row["candidate_accepted_key"] == "COL1"
    assert row["candidate_kingdom"] == "Plantae"
    assert row["checklist_key"] == COL_XR_CHECKLIST_KEY


def test_exact_synonym_maps_to_accepted_candidate() -> None:
    result = {
        "usage": {
            "key": "SYN1",
            "name": "Alpha old",
            "canonicalName": "Alpha old",
            "rank": "SPECIES",
            "status": "SYNONYM",
        },
        "acceptedUsage": {
            "key": "ACC1",
            "name": "Alpha new Author",
            "canonicalName": "Alpha new",
            "authorship": "Author",
            "rank": "SPECIES",
            "status": "ACCEPTED",
        },
        "classification": _classification(),
        "diagnostics": {"matchType": "EXACT", "confidence": 100},
        "synonym": True,
    }
    row = normalize_match(_row("Alpha old"), result)
    assert row["automatic_resolution_candidate"] == "true"
    assert row["resolution_status"] == "exact_synonym_to_accepted_candidate"
    assert row["matched_usage_key"] == "SYN1"
    assert row["candidate_accepted_key"] == "ACC1"
    assert row["candidate_accepted_name"] == "Alpha new"


def test_variant_or_low_confidence_stays_review_required() -> None:
    result = {
        "usage": {
            "key": "COL1",
            "canonicalName": "Alpha one",
            "rank": "SPECIES",
            "status": "ACCEPTED",
        },
        "classification": _classification(),
        "diagnostics": {
            "matchType": "VARIANT",
            "confidence": 92,
            "processingFlags": ["LOW_CONFIDENCE"],
        },
        "synonym": False,
    }
    row = normalize_match(_row(), result)
    assert row["automatic_resolution_candidate"] == "false"
    assert row["resolution_status"] == "review_required"
    assert row["processing_flags"] == "LOW_CONFIDENCE"


def test_none_is_unmatched_not_higher_rank_promotion() -> None:
    row = normalize_match(
        _row(),
        {"diagnostics": {"matchType": "NONE", "confidence": 0, "processingFlags": ["NO_MATCH"]}},
    )
    assert row["automatic_resolution_candidate"] == "false"
    assert row["resolution_status"] == "unmatched"
    assert row["candidate_accepted_key"] == ""


def test_deterministic_pilot_spans_full_sorted_table() -> None:
    source = pd.DataFrame(
        {
            "accepted_species": [f"Species {i:03d}" for i in range(100)],
            "genus": ["Species"] * 100,
            "family": ["Testaceae"] * 100,
        }
    )
    sample = deterministic_pilot(source, 10)
    assert sample.index.tolist() == list(range(10))
    assert sample["accepted_species"].tolist() == [f"Species {i:03d}" for i in range(0, 100, 10)]


def test_summary_separates_automatic_review_and_unmatched() -> None:
    frame = pd.DataFrame(
        [
            {"automatic_resolution_candidate": "true", "resolution_status": "exact_accepted_candidate", "match_type": "EXACT"},
            {"automatic_resolution_candidate": "false", "resolution_status": "review_required", "match_type": "VARIANT"},
            {"automatic_resolution_candidate": "false", "resolution_status": "unmatched", "match_type": "NONE"},
        ]
    )
    report = summarize(frame)
    assert report["n_input"] == 3
    assert report["n_automatic_resolution_candidates"] == 1
    assert report["n_review_required"] == 1
    assert report["n_unmatched"] == 1
