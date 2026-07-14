import csv
import json
from pathlib import Path

from island_v2.free_bulk_trait_pilot import audit_bulk_source, main


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_audit_reports_exact_matches_and_unmapped_traits(tmp_path: Path) -> None:
    source = tmp_path / "source.csv"
    _write_csv(
        source,
        ["taxon_name", "trait_name", "value"],
        [
            {"taxon_name": "Campanula punctata", "trait_name": "flower_length", "value": "4.2"},
            {"taxon_name": "Campanula  punctata", "trait_name": "flower_diameter", "value": "3.1"},
            {"taxon_name": "Campanula punctata", "trait_name": "flower_length", "value": "4.5"},
            {"taxon_name": "Other species", "trait_name": "unknown_trait", "value": "x"},
        ],
    )

    report = audit_bulk_source(
        source,
        {"Campanula punctata"},
        taxon_column="taxon_name",
        trait_column="trait_name",
    )

    assert report["input_records"] == 4
    assert report["unique_source_taxa"] == 2
    assert report["exact_master_matches"] == 1
    assert report["mapped_records"] == 3
    assert report["mapped_records_in_master"] == 3
    assert report["candidate_fallback_replacement_cells"] == 2
    assert report["candidate_fallback_replacement_cells_by_target_trait"] == {
        "flower_diameter_mm": 1,
        "flower_length_mm": 1,
    }
    assert report["unmapped_trait_records"] == 1
    assert report["policy"]["fuzzy_matching"] is False


def test_cli_writes_json_report(tmp_path: Path) -> None:
    master = tmp_path / "master.csv"
    source = tmp_path / "source.csv"
    output = tmp_path / "report.json"
    _write_csv(master, ["scientific_name"], [{"scientific_name": "Campanula punctata"}])
    _write_csv(
        source,
        ["taxon_name", "trait_name"],
        [{"taxon_name": "Campanula punctata", "trait_name": "flower_perianth_fusion"}],
    )

    assert main([
        "--source-csv", str(source),
        "--master-csv", str(master),
        "--output-json", str(output),
    ]) == 0

    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["exact_master_matches"] == 1
    assert payload["mapped_records_by_target_trait"] == {
        "flower_perianth_fusion_fraction": 1
    }
    assert payload["candidate_fallback_replacement_cells"] == 1
