from __future__ import annotations

import zipfile
from pathlib import Path

import pandas as pd
import pytest

from island_v2.floraweb_synonym_recovery import (
    _collapse_recovered_evidence,
    _wfo_accepted_species,
    audit_two_backbone_responses,
    build_wfo_local_prefilter,
)
from island_v2.floraweb_trait_source import STRICT_EVIDENCE_COLUMNS


def test_local_wfo_prefilter_joins_synonym_taxon_id_to_taxon_name_id(
    tmp_path: Path,
) -> None:
    archive_path = tmp_path / "wfo.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr(
            "name.tsv",
            "ID\tscientificName\trank\n"
            "name-syn\tOldus name\tspecies\n"
            "name-accepted\tNewus name\tspecies\n",
        )
        archive.writestr(
            "synonym.tsv",
            "nameID\ttaxonID\nname-syn\ttaxon-accepted\n",
        )
        archive.writestr(
            "taxon.tsv",
            "ID\tnameID\ntaxon-accepted\tname-accepted\n",
        )
    audit = pd.DataFrame([{"source_name": "Oldus name", "accepted": False}])
    master = pd.DataFrame([{"accepted_species": "Newus name", "family": "Newaceae"}])

    result = build_wfo_local_prefilter(
        wfo_zip=archive_path,
        floraweb_name_audit=audit,
        master=master,
        universe_species={"Newus name"},
    )

    assert result.loc[0, "local_candidate_accepted_species"] == "Newus name"
    assert "taxon-accepted" in result.loc[0, "wfo_routes_json"]
    assert "name-accepted" in result.loc[0, "wfo_routes_json"]


def _response(*, gbif_species: str = "Newus name", placement_family: str = "Newaceae"):
    return {
        "source_name": "Oldus name",
        "local_candidate_accepted_species": "Newus name",
        "retrieved_at_utc": "2026-08-28T00:00:00Z",
        "wfo_status": 200,
        "wfo": {
            "parsedName": {"rank": "species", "canonical_form": "Oldus name"},
            "params": {
                "fuzzyNameParts": 0,
                "checkHomonyms": True,
                "checkRank": True,
                "classificationVersion": "2026-06",
            },
            "match": {
                "wfo_id": "wfo-new",
                "placement": (
                    f"Code/Plantae/{placement_family}/Newus/name"
                    "$Code/Plantae/Oldaceae/Oldus/name"
                ),
            },
        },
        "gbif_status": 200,
        "gbif": {
            "matchType": "EXACT",
            "rank": "SPECIES",
            "kingdom": "Plantae",
            "species": gbif_species,
            "family": "Newaceae",
            "usageKey": 1,
            "acceptedUsageKey": 2,
        },
    }


def test_wfo_accepted_name_uses_placement_before_synonym_separator() -> None:
    assert _wfo_accepted_species(_response()) == "Newus name"


@pytest.mark.parametrize(
    ("record", "reason"),
    [
        (_response(gbif_species="Otherus name"), "backbones_disagree"),
        (_response(placement_family="Wrongaceae"), "family_conflict"),
    ],
)
def test_two_backbone_gate_rejects_disagreement_and_family_conflict(
    record: dict[str, object], reason: str
) -> None:
    master = pd.DataFrame([{"accepted_species": "Newus name", "family": "Newaceae"}])
    result = audit_two_backbone_responses([record], master, {"Newus name"})
    assert result.loc[0, "mapping_reason"] == reason
    assert result.loc[0, "accepted_species"] == ""


def test_two_backbone_gate_accepts_only_strict_consensus() -> None:
    master = pd.DataFrame([{"accepted_species": "Newus name", "family": "Newaceae"}])
    result = audit_two_backbone_responses([_response()], master, {"Newus name"})
    assert result.loc[0, "mapping_reason"] == "accepted_strict_two_backbone"
    assert result.loc[0, "accepted_species"] == "Newus name"


def test_recovered_provider_rows_collapse_to_one_underlying_lineage() -> None:
    rows = []
    for usage_id in ("one", "two"):
        row = {column: "" for column in STRICT_EVIDENCE_COLUMNS}
        row.update(
            {
                "accepted_species": "Newus name",
                "trait_name": "flower_primary_color",
                "normalized_value": "white",
                "source_record_id": f"floraweb:name-use-id:{usage_id}:colour",
                "source_url": f"https://example.test/{usage_id}",
                "source_excerpt": "Blüten weiß.",
                "quality": "medium",
            }
        )
        rows.append(row)
    evidence, provenance, conflicts = _collapse_recovered_evidence(
        pd.DataFrame(rows),
        original_names_by_usage={"one": {"Oldus name"}, "two": {"Oldus name"}},
        mapping_by_name={
            "Oldus name": {
                "wfo_classification_version": "2026-06",
                "wfo_match_id": "wfo-new",
                "gbif_accepted_usage_key": "2",
            }
        },
        exact_provider_keys=set(),
        independent=False,
    )
    assert len(evidence) == 1
    assert evidence.loc[0, "source_lineage"].startswith(
        "rothmaler:floraweb:accepted-species:"
    )
    assert provenance.loc[0, "supporting_rows"] == 2
    assert conflicts.empty

    excluded, excluded_provenance, _ = _collapse_recovered_evidence(
        pd.DataFrame(rows),
        original_names_by_usage={"one": {"Oldus name"}, "two": {"Oldus name"}},
        mapping_by_name={},
        exact_provider_keys={("Newus name", "flower_primary_color")},
        independent=False,
    )
    assert excluded.empty
    assert excluded_provenance.loc[0, "decision"] == (
        "excluded_existing_exact_floraweb_species_trait"
    )
