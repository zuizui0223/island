from __future__ import annotations

import zipfile
from pathlib import Path

import pandas as pd

from island_v2.eol_traitbank_inventory_fast import (
    scan_fast_inventory,
    selected_trait_columns,
)


def _archive(path: Path) -> None:
    terms = pd.DataFrame(
        [
            {
                "uri": "http://eol.org/schema/terms/PollinationSyndrome",
                "name": "pollination syndrome",
                "type": "measurement",
            },
            {
                "uri": "eol:value:bee",
                "name": "bee pollination",
                "type": "value",
            },
        ]
    )
    traits = pd.DataFrame(
        [
            {
                "eol_pk": "t1",
                "page_id": "1",
                "predicate": "http://eol.org/schema/terms/PollinationSyndrome",
                "object_term": "eol:value:bee",
                "source": "unused provenance column",
                "measurement_method": "unused method column",
                "metadata": "unused metadata column",
            },
            {
                "eol_pk": "t2",
                "page_id": "2",
                "predicate": "http://rs.tdwg.org/dwc/terms/habitat",
                "object_term": "forest",
                "source": "unused provenance column",
                "measurement_method": "unused method column",
                "metadata": "unused metadata column",
            },
        ]
    )
    pages = pd.DataFrame(
        [
            {"page_id": "1", "canonical": "Campanula microdonta"},
            {"page_id": "2", "canonical": "Examplea habitatensis"},
        ]
    )
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("trait_bank/terms.csv", terms.to_csv(index=False))
        archive.writestr("trait_bank/traits.csv", traits.to_csv(index=False))
        archive.writestr("trait_bank/pages.csv", pages.to_csv(index=False))


def test_selected_trait_columns_omit_irrelevant_fields(tmp_path: Path) -> None:
    archive = tmp_path / "traits_all.zip"
    _archive(archive)

    all_columns, selected = selected_trait_columns(archive)

    assert set(all_columns) == {
        "eol_pk",
        "page_id",
        "predicate",
        "object_term",
        "source",
        "measurement_method",
        "metadata",
    }
    assert selected == ["predicate", "object_term"]


def test_fast_inventory_matches_relevant_predicate_without_all_columns(
    tmp_path: Path,
) -> None:
    archive = tmp_path / "traits_all.zip"
    output = tmp_path / "inventory"
    _archive(archive)

    report = scan_fast_inventory(
        eol_archive=archive,
        output_dir=output,
        chunksize=1,
    )

    predicates = pd.read_csv(report["predicate_inventory_csv"])
    values = pd.read_csv(report["relevant_value_inventory_csv"])
    pollination = predicates.loc[
        predicates["predicate"].eq(
            "http://eol.org/schema/terms/PollinationSyndrome"
        )
    ].iloc[0]

    assert report["n_trait_rows_scanned"] == 2
    assert report["n_columns_in_traits_table"] == 7
    assert report["n_columns_read"] == 2
    assert report["columns_read"] == ["predicate", "object_term"]
    assert int(pollination["n_rows"]) == 1
    assert pollination["currently_mapped"]
    assert values["value_label"].tolist() == ["bee pollination"]
