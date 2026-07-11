from __future__ import annotations

import json
import zipfile
from pathlib import Path

import pandas as pd

from island_v2.eol_traitbank_inventory import scan_eol_traitbank_inventory


def _write_archive(path: Path) -> None:
    terms = pd.DataFrame(
        [
            {
                "uri": "http://rs.tdwg.org/dwc/terms/habitat",
                "name": "habitat",
                "type": "measurement",
            },
            {
                "uri": "http://eol.org/schema/terms/PollinationSyndrome",
                "name": "pollination syndrome",
                "type": "measurement",
            },
            {
                "uri": "http://example.org/flowerTubeLength",
                "name": "flower tube length",
                "type": "measurement",
            },
            {"uri": "eol:value:wind", "name": "wind pollination", "type": "value"},
            {"uri": "eol:value:bee", "name": "bee pollination", "type": "value"},
        ]
    )
    traits = pd.DataFrame(
        [
            {
                "eol_pk": "t1",
                "page_id": "1",
                "predicate": "http://rs.tdwg.org/dwc/terms/habitat",
                "object_term": "marine",
            },
            {
                "eol_pk": "t2",
                "page_id": "2",
                "predicate": "http://eol.org/schema/terms/PollinationSyndrome",
                "object_term": "eol:value:wind",
            },
            {
                "eol_pk": "t3",
                "page_id": "3",
                "predicate": "http://eol.org/schema/terms/PollinationSyndrome",
                "object_term": "eol:value:bee",
            },
            {
                "eol_pk": "t4",
                "page_id": "4",
                "predicate": "http://example.org/flowerTubeLength",
                "normal_measurement": "24",
                "normal_units": "mm",
            },
        ]
    )
    pages = pd.DataFrame(
        [
            {"page_id": "1", "canonical": "Marine thing"},
            {"page_id": "2", "canonical": "Poa annua"},
            {"page_id": "3", "canonical": "Campanula microdonta"},
            {"page_id": "4", "canonical": "Campanula punctata"},
        ]
    )
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("trait_bank/terms.csv", terms.to_csv(index=False))
        archive.writestr("trait_bank/traits.csv", traits.to_csv(index=False))
        archive.writestr("trait_bank/pages.csv", pages.to_csv(index=False))


def test_inventory_scans_all_predicates_and_only_relevant_values(tmp_path: Path) -> None:
    archive = tmp_path / "traits_all.zip"
    output = tmp_path / "inventory"
    _write_archive(archive)

    report = scan_eol_traitbank_inventory(
        eol_archive=archive,
        output_dir=output,
        chunksize=2,
    )

    predicates = pd.read_csv(report["predicate_inventory_csv"], dtype=str).fillna("")
    values = pd.read_csv(report["relevant_value_inventory_csv"], dtype=str).fillna("")
    terms = pd.read_csv(report["relevant_terms_csv"], dtype=str).fillna("")
    manifest = json.loads(Path(report["manifest_json"]).read_text(encoding="utf-8"))

    pollination = predicates.loc[
        predicates["predicate"].eq("http://eol.org/schema/terms/PollinationSyndrome")
    ].iloc[0]
    assert int(pollination["n_rows"]) == 2
    assert pollination["currently_mapped"].lower() == "true"

    assert set(values["value_label"]) == {"wind pollination", "bee pollination", "24 mm"}
    assert "marine" not in set(values["value"])
    assert "flower tube length" in set(terms["name"])
    assert manifest["n_trait_rows_scanned"] == 4
    assert manifest["n_observed_predicates"] == 3
    assert manifest["n_currently_mapped_predicates_observed"] == 1


def test_inventory_honours_smoke_row_cap(tmp_path: Path) -> None:
    archive = tmp_path / "traits_all.zip"
    output = tmp_path / "inventory"
    _write_archive(archive)

    report = scan_eol_traitbank_inventory(
        eol_archive=archive,
        output_dir=output,
        chunksize=3,
        max_trait_rows=2,
    )

    assert report["n_trait_rows_scanned"] == 2
    predicates = pd.read_csv(report["predicate_inventory_csv"], dtype=str).fillna("")
    assert int(predicates["n_rows"].astype(int).sum()) == 2
