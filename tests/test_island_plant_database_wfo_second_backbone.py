from __future__ import annotations

import csv
import io
import zipfile
from pathlib import Path

import pandas as pd

import island_v2.island_plant_database_wfo_second_backbone as wfo


def _write_backbone(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "taxonID",
        "scientificName",
        "family",
        "taxonRank",
        "taxonomicStatus",
        "acceptedNameUsageID",
        "kingdom",
    ]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("classification.csv", buffer.getvalue())


def _queue() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "submitted_name": "Alpha one",
                "submitted_family": "Alphaaceae",
                "triage_class": "colxr_none_no_match",
                "candidate_accepted_name": "",
                "candidate_accepted_rank": "",
                "second_backbone_eligible": "true",
            },
            {
                "submitted_name": "Beta old",
                "submitted_family": "Betaceae",
                "triage_class": "colxr_variant_species_name",
                "candidate_accepted_name": "Beta new",
                "candidate_accepted_rank": "SPECIES",
                "second_backbone_eligible": "true",
            },
            {
                "submitted_name": "Gamma three",
                "submitted_family": "Gammaaceae",
                "triage_class": "colxr_none_no_match",
                "candidate_accepted_name": "",
                "candidate_accepted_rank": "",
                "second_backbone_eligible": "true",
            },
            {
                "submitted_name": "Delta four",
                "submitted_family": "Deltaaceae",
                "triage_class": "colxr_none_no_match",
                "candidate_accepted_name": "",
                "candidate_accepted_rank": "",
                "second_backbone_eligible": "true",
            },
        ]
    )


def test_exact_and_synonym_rescue_candidates_are_fail_closed(tmp_path: Path, monkeypatch) -> None:
    backbone = tmp_path / "wfo.zip"
    _write_backbone(
        backbone,
        [
            {
                "taxonID": "wfo-1",
                "scientificName": "Alpha one",
                "family": "Alphaaceae",
                "taxonRank": "species",
                "taxonomicStatus": "Accepted",
                "acceptedNameUsageID": "",
                "kingdom": "Plantae",
            },
            {
                "taxonID": "wfo-2",
                "scientificName": "Beta old",
                "family": "Betaceae",
                "taxonRank": "species",
                "taxonomicStatus": "Synonym",
                "acceptedNameUsageID": "wfo-3",
                "kingdom": "Plantae",
            },
            {
                "taxonID": "wfo-3",
                "scientificName": "Beta new",
                "family": "Betaceae",
                "taxonRank": "species",
                "taxonomicStatus": "Accepted",
                "acceptedNameUsageID": "",
                "kingdom": "Plantae",
            },
            {
                "taxonID": "wfo-4",
                "scientificName": "Gamma three",
                "family": "Wrongaceae",
                "taxonRank": "species",
                "taxonomicStatus": "Accepted",
                "acceptedNameUsageID": "",
                "kingdom": "Plantae",
            },
        ],
    )
    monkeypatch.setattr(wfo, "EXPECTED_QUEUE_ROWS", 4)
    audit, summary = wfo.audit_queue(_queue(), backbone)

    assert audit["wfo_resolution_status"].tolist() == [
        "exact_species_rescue_candidate",
        "exact_species_rescue_candidate",
        "exact_species_family_conflict",
        "no_exact_accepted_species_concept",
    ]
    assert audit["wfo_rescue_candidate"].tolist() == ["true", "true", "false", "false"]
    assert audit.loc[1, "wfo_accepted_name"] == "Beta new"
    assert audit.loc[1, "wfo_colxr_target_concordance"] == "same_target"
    assert summary["n_wfo_rescue_candidates"] == 2
    assert summary["n_manual_review"] == 2


def test_different_colxr_target_and_multiple_wfo_targets_do_not_rescue(tmp_path: Path, monkeypatch) -> None:
    backbone = tmp_path / "wfo.zip"
    _write_backbone(
        backbone,
        [
            {
                "taxonID": "wfo-10",
                "scientificName": "Alpha one",
                "family": "Alphaaceae",
                "taxonRank": "species",
                "taxonomicStatus": "Accepted",
                "acceptedNameUsageID": "",
                "kingdom": "Plantae",
            },
            {
                "taxonID": "wfo-20",
                "scientificName": "Beta old",
                "family": "Betaceae",
                "taxonRank": "species",
                "taxonomicStatus": "Synonym",
                "acceptedNameUsageID": "wfo-21",
                "kingdom": "Plantae",
            },
            {
                "taxonID": "wfo-22",
                "scientificName": "Beta old",
                "family": "Betaceae",
                "taxonRank": "species",
                "taxonomicStatus": "Synonym",
                "acceptedNameUsageID": "wfo-23",
                "kingdom": "Plantae",
            },
            {
                "taxonID": "wfo-21",
                "scientificName": "Beta first",
                "family": "Betaceae",
                "taxonRank": "species",
                "taxonomicStatus": "Accepted",
                "acceptedNameUsageID": "",
                "kingdom": "Plantae",
            },
            {
                "taxonID": "wfo-23",
                "scientificName": "Beta second",
                "family": "Betaceae",
                "taxonRank": "species",
                "taxonomicStatus": "Accepted",
                "acceptedNameUsageID": "",
                "kingdom": "Plantae",
            },
        ],
    )
    queue = _queue().iloc[:2].copy()
    queue.loc[0, "candidate_accepted_name"] = "Alpha other"
    monkeypatch.setattr(wfo, "EXPECTED_QUEUE_ROWS", 2)
    audit, _ = wfo.audit_queue(queue, backbone)

    assert audit.loc[0, "wfo_resolution_status"] == "exact_species_colxr_target_conflict"
    assert audit.loc[0, "wfo_rescue_candidate"] == "false"
    assert audit.loc[1, "wfo_resolution_status"] == "ambiguous_exact_name_multiple_species_targets"
    assert audit.loc[1, "wfo_unique_accepted_species_targets"] == 2
