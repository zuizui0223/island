import csv
import hashlib
from pathlib import Path

import pytest

from island_v2.fetch_austraits_bulk import (
    _verify_declared_checksum,
    discover_long_trait_table,
    select_latest_release,
    select_plain_text_archive,
    standardize_dataset,
)


def test_selects_latest_release_and_zip_archive() -> None:
    payload = {
        "hits": {
            "hits": [
                {
                    "metadata": {"publication_date": "2024-01-01", "version": "5.0"},
                    "files": [
                        {"key": "austraits-5.0.rds", "links": {"self": "https://example/rds"}}
                    ],
                },
                {
                    "metadata": {"publication_date": "2025-01-01", "version": "6.0"},
                    "files": [
                        {"key": "austraits-6.0.rds", "links": {"self": "https://example/rds2"}},
                        {"key": "austraits-6.0.zip", "links": {"self": "https://example/zip"}},
                    ],
                },
            ]
        }
    }
    release = select_latest_release(payload)
    assert release["metadata"]["version"] == "6.0"
    assert select_plain_text_archive(release)["key"] == "austraits-6.0.zip"


def test_discovers_and_standardizes_requested_dataset(tmp_path: Path) -> None:
    irrelevant = tmp_path / "taxa.csv"
    irrelevant.write_text("taxon_name,family\nA a,Aaceae\n", encoding="utf-8")
    traits = tmp_path / "traits.csv"
    with traits.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["taxon_name", "trait_name", "trait_value", "unit", "dataset_id"],
        )
        writer.writeheader()
        writer.writerows(
            [
                {
                    "taxon_name": "Campanula  punctata",
                    "trait_name": "flower_length",
                    "trait_value": "42",
                    "unit": "mm",
                    "dataset_id": "eFLOWER_Dun_2022",
                },
                {
                    "taxon_name": "Other species",
                    "trait_name": "leaf_area",
                    "trait_value": "3",
                    "unit": "mm2",
                    "dataset_id": "other_dataset",
                },
            ]
        )

    path, columns = discover_long_trait_table(tmp_path)
    assert path == traits
    output = tmp_path / "standardized.csv"
    stats = standardize_dataset(
        path,
        output,
        columns=columns,
        dataset_key="eFLOWER_Dun_2022",
    )

    assert stats == {
        "input_rows_in_discovered_table": 2,
        "standardized_rows": 1,
        "unique_taxa": 1,
        "unique_traits": 1,
        "unique_nonempty_units": 1,
    }
    rows = list(csv.DictReader(output.open(encoding="utf-8")))
    assert rows[0]["taxon_name"] == "Campanula punctata"
    assert rows[0]["trait_name"] == "flower_length"
    assert rows[0]["trait_unit"] == "mm"
    assert rows[0]["source_record_id"] == "austraits:eFLOWER_Dun_2022:1"


def test_standardize_can_filter_traits_across_all_datasets(tmp_path: Path) -> None:
    source = tmp_path / "traits.csv"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["taxon_name", "trait_name", "value", "dataset_id"],
        )
        writer.writeheader()
        writer.writerows(
            [
                {
                    "taxon_name": "Poa annua",
                    "trait_name": "sex_type",
                    "value": "hermaphrodite",
                    "dataset_id": "Dataset_A",
                },
                {
                    "taxon_name": "Poa annua",
                    "trait_name": "leaf_area",
                    "value": "4",
                    "dataset_id": "Dataset_A",
                },
                {
                    "taxon_name": "Acacia dealbata",
                    "trait_name": "pollination_syndrome",
                    "value": "insect",
                    "dataset_id": "Dataset_B",
                },
            ]
        )

    output = tmp_path / "selected.csv"
    stats = standardize_dataset(
        source,
        output,
        columns={
            "taxon": "taxon_name",
            "trait": "trait_name",
            "value": "value",
            "unit": None,
            "dataset": "dataset_id",
        },
        dataset_key=None,
        trait_names={"sex_type", "pollination_syndrome"},
    )

    rows = list(csv.DictReader(output.open(encoding="utf-8")))
    assert stats["standardized_rows"] == 2
    assert [row["trait_name"] for row in rows] == ["sex_type", "pollination_syndrome"]
    assert [row["dataset_key"] for row in rows] == ["Dataset_A", "Dataset_B"]
    assert [row["source_record_id"] for row in rows] == [
        "austraits:Dataset_A:1",
        "austraits:Dataset_B:3",
    ]


def test_verifies_declared_archive_checksum(tmp_path: Path) -> None:
    archive = tmp_path / "austraits.zip"
    archive.write_bytes(b"verified public archive")
    expected = hashlib.md5(archive.read_bytes()).hexdigest()  # noqa: S324 - source checksum contract

    assert _verify_declared_checksum(archive, f"md5:{expected}") is True
    with pytest.raises(ValueError, match="checksum mismatch"):
        _verify_declared_checksum(archive, "md5:00000000000000000000000000000000")
