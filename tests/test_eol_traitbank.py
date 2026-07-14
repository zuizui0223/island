from __future__ import annotations

import json
import shutil
import zipfile
from pathlib import Path

import pandas as pd

from island_v2.bulk_trait_ingest import ingest_bulk_traits
from island_v2.eol_traitbank import prepare_eol_traitbank_long


def _write_eol_archive(path: Path) -> None:
    pages = pd.DataFrame(
        [
            {"page_id": "10", "canonical": "Poa annua"},
            {"page_id": "20", "canonical": "Campanula microdonta"},
            {"page_id": "30", "canonical": "Unmatched plant"},
        ]
    )
    terms = pd.DataFrame(
        [
            {
                "uri": "http://eol.org/schema/terms/PollinationSyndrome",
                "name": "pollination syndrome",
                "type": "measurement",
            },
            {
                "uri": "http://eol.org/schema/terms/FlowerColor",
                "name": "flower color",
                "type": "measurement",
            },
            {"uri": "eol:value:wind", "name": "wind pollination", "type": "value"},
            {"uri": "eol:value:bee", "name": "bee pollination", "type": "value"},
            {"uri": "eol:value:white", "name": "white", "type": "value"},
            {"uri": "eol:value:unsupported", "name": "chartreuse", "type": "value"},
            {"uri": "eol:trait:unsupported", "name": "leaf size", "type": "measurement"},
        ]
    )
    traits = pd.DataFrame(
        [
            {
                "eol_pk": "t1",
                "resource_pk": "r1",
                "page_id": "10",
                "scientific_name": "Poa annua L.",
                "predicate": "http://eol.org/schema/terms/PollinationSyndrome",
                "value_uri": "eol:value:wind",
                "citation": "EOL test fixture",
            },
            {
                "eol_pk": "t2",
                "resource_pk": "r2",
                "page_id": "20",
                "scientific_name": "",
                "predicate": "http://eol.org/schema/terms/FlowerColor",
                # The real 2024 EOL export calls this column value_uri.
                "value_uri": "eol:value:white",
                "citation": "EOL test fixture",
            },
            {
                "eol_pk": "t3",
                "resource_pk": "r3",
                "page_id": "20",
                "scientific_name": "Campanula microdonta",
                "predicate": "http://eol.org/schema/terms/PollinationSyndrome",
                "value_uri": "eol:value:bee",
                "citation": "EOL test fixture",
            },
            {
                "eol_pk": "t4",
                "resource_pk": "r4",
                "page_id": "20",
                "scientific_name": "Campanula microdonta",
                "predicate": "http://eol.org/schema/terms/FlowerColor",
                "value_uri": "eol:value:unsupported",
                "citation": "EOL test fixture",
            },
            {
                "eol_pk": "t5",
                "resource_pk": "r5",
                "page_id": "30",
                "scientific_name": "Unmatched plant",
                "predicate": "eol:trait:unsupported",
                "value_uri": "eol:value:white",
                "citation": "EOL test fixture",
            },
        ]
    )
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("trait_bank/pages.csv", pages.to_csv(index=False))
        archive.writestr("trait_bank/terms.csv", terms.to_csv(index=False))
        archive.writestr("trait_bank/traits.csv", traits.to_csv(index=False))
        archive.writestr("trait_bank/metadata.csv", "eol_pk,predicate,object_term\n")
        archive.writestr("trait_bank/inferred.csv", "page_id,eol_pk\n")


def _write_master(path: Path) -> None:
    pd.DataFrame(
        [
            {
                "accepted_species": "Poa annua",
                "genus": "Poa",
                "family": "Poaceae",
                "n_islands": 3,
                "n_records": 20,
            },
            {
                "accepted_species": "Campanula microdonta",
                "genus": "Campanula",
                "family": "Campanulaceae",
                "n_islands": 2,
                "n_records": 4,
            },
        ]
    ).to_csv(path, index=False)


def test_prepare_eol_traitbank_and_join_scientific_names(tmp_path: Path) -> None:
    archive = tmp_path / "traits_all.zip"
    master = tmp_path / "master.csv"
    prepared = tmp_path / "prepared"
    ingested = tmp_path / "ingested"
    _write_eol_archive(archive)
    _write_master(master)

    prepare_report = prepare_eol_traitbank_long(
        eol_archive=archive,
        output_dir=prepared,
        config_path=Path("config/eol_traitbank_terms.yml"),
        source_name="eol_test",
        chunksize=2,
    )

    long = pd.read_csv(prepare_report["long_csv"], dtype=str).fillna("")
    assert set(zip(long["trait_name"], long["trait_value"], strict=False)) == {
        ("pollen_vector_mode", "abiotic_wind"),
        ("flower_primary_color", "white"),
        ("pollination_functional_guild", "other_bees"),
    }
    assert (
        long.loc[long["source_record_id"].eq("t2"), "scientific_name"].item()
        == "Campanula microdonta"
    )

    unmapped = pd.read_csv(prepare_report["unmapped_terms_csv"], dtype=str).fillna("")
    assert set(unmapped["reason"]) == {"unmapped_value", "unmapped_predicate"}
    assert prepare_report["n_trait_rows_scanned"] == 5
    assert prepare_report["n_long_rows"] == 3

    ingest_report = ingest_bulk_traits(
        master_csv=master,
        source_csv=Path(prepare_report["long_csv"]),
        output_dir=ingested,
        profiles_path=Path("config/bulk_trait_sources.yml"),
        profile_name="v2_canonical_long_csv",
        source_name="eol_test",
    )
    candidates = pd.read_csv(ingest_report["candidate_csv"], dtype=str).fillna("")
    assert len(candidates) == 3
    assert candidates.loc[candidates["source_record_id"].eq("t1"), "name_match_method"].item() == (
        "accepted_name_author_stripped"
    )
    assert candidates.loc[candidates["source_record_id"].eq("t2"), "name_match_method"].item() == (
        "accepted_name_exact"
    )
    manifest = json.loads(Path(prepare_report["manifest_json"]).read_text(encoding="utf-8"))
    assert manifest["inferred_csv_policy"].startswith("EOL inferred.csv is not used")


def test_prepare_eol_can_read_extracted_content_and_keep_archive_checksum(
    tmp_path: Path,
) -> None:
    archive = tmp_path / "traits_all.zip"
    extracted = tmp_path / "extracted"
    prepared = tmp_path / "prepared"
    _write_eol_archive(archive)
    shutil.unpack_archive(archive, extracted)

    report = prepare_eol_traitbank_long(
        eol_archive=archive,
        eol_content_dir=extracted,
        output_dir=prepared,
        config_path=Path("config/eol_traitbank_terms.yml"),
        source_name="eol_extracted_test",
        chunksize=2,
    )

    assert report["eol_content_path"] == str(extracted)
    assert len(report["eol_archive_sha256"]) == 64
    assert report["n_trait_rows_scanned"] == 5
    assert report["n_long_rows"] == 3
