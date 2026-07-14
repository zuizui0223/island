import csv
from pathlib import Path

from island_v2.fetch_austraits_bulk import (
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
                    "files": [{"key": "austraits-5.0.rds", "links": {"self": "https://example/rds"}}],
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
            fieldnames=["taxon_name", "trait_name", "trait_value", "dataset_id"],
        )
        writer.writeheader()
        writer.writerows(
            [
                {
                    "taxon_name": "Campanula  punctata",
                    "trait_name": "flower_length",
                    "trait_value": "42",
                    "dataset_id": "eFLOWER_Dun_2022",
                },
                {
                    "taxon_name": "Other species",
                    "trait_name": "leaf_area",
                    "trait_value": "3",
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
    }
    rows = list(csv.DictReader(output.open(encoding="utf-8")))
    assert rows[0]["taxon_name"] == "Campanula punctata"
    assert rows[0]["trait_name"] == "flower_length"
