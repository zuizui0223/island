from __future__ import annotations

import json

from island_v2.fetch_gift_traits import map_gift_trait_value, prepare_gift_rows


def test_gift_mappings_are_conservative_and_explicit() -> None:
    assert map_gift_trait_value("3.21.1", "orange/yellow") == [
        ("flower_primary_color", "yellow_orange")
    ]
    assert map_gift_trait_value("3.21.1", "pink/white") == [
        ("flower_primary_color", "multicolored_variable")
    ]
    assert map_gift_trait_value("3.4.1", "bisexual") == [("sex_system", "hermaphroditic")]
    assert map_gift_trait_value("3.4.1", "dioecious/monoecious") == []
    assert map_gift_trait_value("3.6.2", "insect/wind") == [
        ("pollen_vector_mode", "mixed"),
        ("pollination_functional_guild", "mixed_or_generalist"),
    ]
    assert map_gift_trait_value("3.6.3", "fly") == [
        ("pollen_vector_mode", "biotic"),
        ("pollination_functional_guild", "flies"),
    ]
    assert map_gift_trait_value("3.6.3", "water") == [
        ("pollination_functional_guild", "wind_water")
    ]
    assert map_gift_trait_value("3.6.3", "other") == []


def _raw(**updates: str) -> dict[str, str]:
    row = {
        "trait_derived_ID": "6805406",
        "ref_ID": "23",
        "orig_ID": "13642",
        "trait_ID": "3.21.1",
        "trait_value": "yellow",
        "derived": "0",
        "bias_deriv": "0",
        "cf_genus": "0",
        "genus": "Acacia",
        "cf_species": "0",
        "aff_species": "0",
        "species_epithet": "spirorbis",
        "subtaxon": "",
        "author": "Labill.",
        "matched": "1",
        "resolved": "1",
        "work_ID": "216",
        "_expected_ref_ID": "23",
        "_expected_trait_ID": "3.21.1",
        "_query_url": "https://gift.example/index3.2.php?query=traits_raw&refid=23",
    }
    row.update(updates)
    return row


def test_prepare_gift_rows_preserves_raw_evidence_and_audits_api_leakage() -> None:
    references = {
        "23": {
            "ref_ID": "23",
            "restricted": "0",
            "ref_long": "Wagner & Lorence (2002) Flora of the Marquesas Islands.",
        }
    }
    species = {"216": {"work_ID": "216", "work_species": "Acacia spirorbis"}}
    rows = [
        _raw(),
        _raw(trait_derived_ID="orphan", ref_ID="", work_ID=""),
        _raw(trait_derived_ID="infra", subtaxon="var. spirorbis"),
    ]

    prepared, audit = prepare_gift_rows(
        rows,
        references_by_id=references,
        species_by_id=species,
        version="3.2",
    )

    assert len(prepared) == 1
    assert prepared.loc[0, "scientific_name"] == "Acacia spirorbis"
    assert prepared.loc[0, "trait_value"] == "yellow_orange"
    assert prepared.loc[0, "source_record_id"].startswith("gift:3.2:6805406:")
    evidence = json.loads(prepared.loc[0, "evidence_excerpt"])
    assert evidence["trait_value"] == "yellow"
    assert evidence["work_species"] == "Acacia spirorbis"
    assert set(audit["reason"]) == {
        "api_query_scope_mismatch",
        "infraspecific_source_scope",
    }


def test_prepare_gift_rows_excludes_restricted_reference() -> None:
    prepared, audit = prepare_gift_rows(
        [_raw()],
        references_by_id={"23": {"ref_ID": "23", "restricted": "1", "ref_long": "x"}},
        species_by_id={"216": {"work_ID": "216", "work_species": "Acacia spirorbis"}},
        version="3.2",
    )

    assert prepared.empty
    assert audit.loc[0, "reason"] == "missing_or_restricted_reference"
