from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from island_v2.island_plant_database_gbif import (
    build_gbif_candidate_core,
    build_island_download_map,
)


def _write_inputs(tmp_path: Path) -> dict[str, Path]:
    pairs = pd.DataFrame(
        [
            {
                "island_id": "isl-a",
                "species": "Alpha one",
                "n_records": "3",
                "n_unique_gbif_ids": "3",
                "basis_of_record_set": "PRESERVED_SPECIMEN|HUMAN_OBSERVATION",
                "review_status": "unresolved_taxonomy",
            },
            {
                "island_id": "isl-a",
                "species": "Beta two",
                "n_records": "2",
                "n_unique_gbif_ids": "2",
                "basis_of_record_set": "HUMAN_OBSERVATION",
                "review_status": "unresolved_taxonomy",
            },
            {
                "island_id": "isl-b",
                "species": "Alpha one",
                "n_records": "1",
                "n_unique_gbif_ids": "1",
                "basis_of_record_set": "PRESERVED_SPECIMEN",
                "review_status": "unresolved_taxonomy",
            },
        ]
    )
    taxa = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "genus": "Alpha",
                "family": "Alphaaceae",
                "n_islands": "2",
                "n_records": "4",
            },
            {
                "accepted_species": "Beta two",
                "genus": "Beta",
                "family": "Betaceae",
                "n_islands": "1",
                "n_records": "2",
            },
        ]
    )
    # isl-a is deliberately represented by two query-side parts to mimic an
    # antimeridian split. Both must collapse back to one analysis island.
    members = pd.DataFrame(
        [
            {
                "block_id": "block-a-west",
                "island_id": "isl-a__antimeridian_part1",
                "analysis_island_id": "isl-a",
            },
            {
                "block_id": "block-a-east",
                "island_id": "isl-a__antimeridian_part2",
                "analysis_island_id": "isl-a",
            },
            {
                "block_id": "block-b",
                "island_id": "isl-b",
                "analysis_island_id": "isl-b",
            },
        ]
    )
    islands = pd.DataFrame([{"island_id": "isl-a"}, {"island_id": "isl-b"}])
    campaign = {
        "campaign_id": "test-campaign",
        "ledger": [
            {
                "block_id": "block-a-west",
                "request_status": "succeeded",
                "download_key": "DL-AW",
                "doi": "10.15468/dl.aw",
            },
            {
                "block_id": "block-a-east",
                "request_status": "succeeded",
                "download_key": "DL-AE",
                "doi": "10.15468/dl.ae",
            },
            {
                "block_id": "block-b",
                "request_status": "succeeded",
                "download_key": "DL-B",
                "doi": "10.15468/dl.b",
            },
        ],
    }

    paths = {
        "pairs": tmp_path / "pairs.csv.gz",
        "taxa": tmp_path / "taxa.csv",
        "members": tmp_path / "members.csv",
        "islands": tmp_path / "islands.csv",
        "campaign": tmp_path / "campaign.json",
    }
    pairs.to_csv(paths["pairs"], index=False, compression="gzip")
    taxa.to_csv(paths["taxa"], index=False)
    members.to_csv(paths["members"], index=False)
    islands.to_csv(paths["islands"], index=False)
    paths["campaign"].write_text(json.dumps(campaign), encoding="utf-8")
    return paths


def test_analysis_island_provenance_collapses_antimeridian_parts(tmp_path: Path) -> None:
    paths = _write_inputs(tmp_path)
    members = pd.read_csv(paths["members"], dtype=str)
    campaign = json.loads(paths["campaign"].read_text(encoding="utf-8"))

    mapping = build_island_download_map(members, campaign)

    assert mapping["island_id"].tolist() == ["isl-a", "isl-b"]
    a = mapping.set_index("island_id").loc["isl-a"]
    assert set(a["download_keys"].split("|")) == {"DL-AE", "DL-AW"}
    assert set(a["download_dois"].split("|")) == {"10.15468/dl.ae", "10.15468/dl.aw"}


def test_build_candidate_core_keeps_presence_provisional(tmp_path: Path) -> None:
    paths = _write_inputs(tmp_path)

    taxa, island_taxa, evidence, rights, manifest = build_gbif_candidate_core(
        paths["pairs"],
        paths["taxa"],
        paths["members"],
        paths["campaign"],
        paths["islands"],
    )

    assert len(taxa) == 2
    assert len(island_taxa) == 3
    assert len(evidence) == 3
    assert len(rights) == 3
    assert island_taxa["membership_status"].eq("candidate").all()
    assert island_taxa["establishment_status"].eq("unknown").all()
    assert island_taxa["endemic_status"].eq("unknown").all()
    assert island_taxa["release_status"].eq("review_required").all()
    assert evidence["rights_status"].eq("review_required").all()
    assert rights["rights_status"].eq("review_required").all()
    assert manifest["scientific_interpretation"]["absence_inferred_from_missing_records"] is False
    assert manifest["counts"]["multi_download_analysis_islands"] == 1
    assert manifest["counts"]["max_downloads_per_analysis_island"] == 2

    a_taxon = taxa.set_index("accepted_name").loc["Alpha one", "taxon_id"]
    row = evidence[(evidence["island_id"] == "isl-a") & (evidence["taxon_id"] == a_taxon)].iloc[0]
    assert set(row["source_record_id"].split("|")) == {"DL-AE", "DL-AW"}
    assert set(row["dataset_key_or_doi"].split("|")) == {
        "10.15468/dl.ae",
        "10.15468/dl.aw",
    }


def test_candidate_core_rejects_unknown_database_island(tmp_path: Path) -> None:
    paths = _write_inputs(tmp_path)
    islands = pd.DataFrame([{"island_id": "isl-a"}])
    islands.to_csv(paths["islands"], index=False)

    with pytest.raises(ValueError, match="outside Database 2.0"):
        build_gbif_candidate_core(
            paths["pairs"],
            paths["taxa"],
            paths["members"],
            paths["campaign"],
            paths["islands"],
        )


def test_block_members_without_analysis_id_remain_backward_compatible(tmp_path: Path) -> None:
    paths = _write_inputs(tmp_path)
    members = pd.DataFrame(
        [
            {"block_id": "block-a-west", "island_id": "isl-a"},
            {"block_id": "block-b", "island_id": "isl-b"},
        ]
    )
    campaign = json.loads(paths["campaign"].read_text(encoding="utf-8"))

    mapping = build_island_download_map(members, campaign)

    assert set(mapping["island_id"]) == {"isl-a", "isl-b"}
