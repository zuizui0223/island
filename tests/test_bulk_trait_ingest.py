from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from island_v2.bulk_trait_ingest import ingest_bulk_traits


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


def _write_source(path: Path) -> None:
    pd.DataFrame(
        [
            {
                "scientific_name": "Poa annua L.",
                "trait_name": "pollen_vector_mode",
                "trait_value": "biotic",
                "source_record_id": "r1",
                "source_citation": "Example curated dataset",
                "source_url": "https://example.org/r1",
                "evidence_excerpt": "Biotic pollen vector recorded.",
            },
            {
                "scientific_name": "Campanula microdonta",
                "trait_name": "flower_primary_color",
                "trait_value": "white",
                "source_record_id": "r2",
                "source_citation": "Example curated dataset",
                "source_url": "https://example.org/r2",
                "evidence_excerpt": "White flower recorded.",
            },
            {
                "scientific_name": "Poa annua subsp. annua",
                "trait_name": "pollen_vector_mode",
                "trait_value": "biotic",
                "source_record_id": "r3",
                "source_citation": "Example curated dataset",
                "source_url": "https://example.org/r3",
                "evidence_excerpt": "An infraspecific record must not be collapsed to the species.",
            },
            {
                "scientific_name": "Campanula microdonta",
                "trait_name": "mating_system",
                "trait_value": "not_recorded",
                "source_record_id": "r4",
                "source_citation": "Example curated dataset",
                "source_url": "https://example.org/r4",
                "evidence_excerpt": "An unsupported source code remains unmapped.",
            },
        ]
    ).to_csv(path, index=False)


def test_bulk_ingest_keeps_direct_evidence_and_audits_non_matches(tmp_path: Path) -> None:
    master = tmp_path / "master.csv"
    source = tmp_path / "source.csv"
    output = tmp_path / "output"
    _write_master(master)
    _write_source(source)

    report = ingest_bulk_traits(
        master_csv=master,
        source_csv=source,
        output_dir=output,
        profiles_path=Path("config/bulk_trait_sources.yml"),
        profile_name="v2_canonical_long_csv",
        source_name="example_bulk_v1",
        source_citation_override="",
        source_url_override="",
    )

    candidates = pd.read_csv(report["candidate_csv"], dtype=str).fillna("")
    assert len(candidates) == 2
    assert set(candidates["accepted_species"]) == {"Poa annua", "Campanula microdonta"}
    assert set(candidates["candidate_kind"]) == {"source_backed"}
    assert set(candidates["evidence_scope"]) == {"species_direct"}
    assert set(candidates["review_status"]) == {"pending"}
    assert candidates.loc[
        candidates["accepted_species"].eq("Poa annua"), "name_match_method"
    ].item() == ("accepted_name_author_stripped")

    unmatched = pd.read_csv(report["unmatched_taxa_csv"], dtype=str).fillna("")
    assert unmatched["source_scientific_name"].tolist() == ["Poa annua subsp. annua"]

    unmapped = pd.read_csv(report["unmapped_values_csv"], dtype=str).fillna("")
    assert unmapped["reason"].tolist() == ["unmapped_or_invalid_trait_value"]
    assert report["n_source_rows"] == 4
    assert report["n_candidate_rows"] == 2
    assert report["n_unmatched_taxon_rows"] == 1
    assert report["n_unmapped_value_rows"] == 1

    coverage = pd.read_csv(report["coverage_by_trait_csv"], dtype=str).fillna("")
    pollen = coverage.loc[coverage["trait_name"].eq("pollen_vector_mode")].iloc[0]
    assert pollen["n_species_with_direct_evidence"] == "1"
    assert pollen["n_island_species_pairs_covered"] == "3"
    assert pollen["analysis_role"] == "global_bulk_priority"

    manifest = json.loads(Path(report["manifest_json"]).read_text(encoding="utf-8"))
    assert manifest["non_inference_rule"].startswith("Unmatched taxa")


def test_bulk_ingest_rejects_unknown_source_profile(tmp_path: Path) -> None:
    master = tmp_path / "master.csv"
    source = tmp_path / "source.csv"
    _write_master(master)
    _write_source(source)

    try:
        ingest_bulk_traits(
            master_csv=master,
            source_csv=source,
            output_dir=tmp_path / "output",
            profiles_path=Path("config/bulk_trait_sources.yml"),
            profile_name="not_a_profile",
            source_name="example",
        )
    except Exception as error:
        assert "Unknown source profile" in str(error)
    else:
        raise AssertionError("Expected an unknown profile to fail explicitly")


def test_austraits_eflower_maps_only_directly_corresponding_categories(tmp_path: Path) -> None:
    master = tmp_path / "master.csv"
    source = tmp_path / "austraits.csv"
    _write_master(master)
    pd.DataFrame(
        [
            {
                "taxon_name": "Campanula microdonta",
                "trait_name": "flower_colour",
                "trait_value": "blue_purple",
                "source_record_id": "a1",
            },
            {
                "taxon_name": "Campanula microdonta",
                "trait_name": "flower_perianth_symmetry",
                "trait_value": "actinomorphic_general",
                "source_record_id": "a2",
            },
            {
                "taxon_name": "Campanula microdonta",
                "trait_name": "flower_perianth_symmetry",
                "trait_value": "disymmetric",
                "source_record_id": "a3",
            },
            {
                "taxon_name": "Campanula microdonta",
                "trait_name": "flower_structural_sex_type",
                "trait_value": "bisexual",
                "source_record_id": "a4",
            },
        ]
    ).to_csv(source, index=False)

    report = ingest_bulk_traits(
        master_csv=master,
        source_csv=source,
        output_dir=tmp_path / "output",
        profiles_path=Path("config/bulk_trait_sources.yml"),
        profile_name="austraits_eflower",
        source_name="austraits_eflower_v7.0.0",
        source_citation_override="AusTraits v7.0.0, eFLOWER_Dun_2022",
        source_url_override="https://doi.org/10.5281/zenodo.15718081",
    )

    candidates = pd.read_csv(report["candidate_csv"], dtype=str).fillna("")
    assert candidates[["trait_name", "standardized_value"]].to_dict("records") == [
        {"trait_name": "flower_primary_color", "standardized_value": "blue_purple"},
        {"trait_name": "floral_symmetry", "standardized_value": "actinomorphic"},
    ]
    assert candidates["source_url"].eq("https://doi.org/10.5281/zenodo.15718081").all()
    assert candidates["evidence_excerpt"].str.contains("austraits_eflower").all()

    unmapped = pd.read_csv(report["unmapped_values_csv"], dtype=str).fillna("")
    assert set(unmapped["source_record_id"]) == {"a3", "a4"}
    assert set(unmapped["reason"]) == {
        "unmapped_trait_name",
        "unmapped_or_invalid_trait_value",
    }


def test_austraits_all_profiles_separate_vector_and_functional_guild(tmp_path: Path) -> None:
    master = tmp_path / "master.csv"
    source = tmp_path / "austraits_all.csv"
    _write_master(master)
    pd.DataFrame(
        [
            {
                "taxon_name": "Campanula microdonta",
                "trait_name": "flower_colour",
                "trait_value": "blue_purple white_cream",
                "source_record_id": "all1",
            },
            {
                "taxon_name": "Campanula microdonta",
                "trait_name": "sex_type",
                "trait_value": "dioecious",
                "source_record_id": "all2",
            },
            {
                "taxon_name": "Campanula microdonta",
                "trait_name": "pollination_syndrome",
                "trait_value": "insect wind",
                "source_record_id": "all3",
            },
            {
                "taxon_name": "Poa annua",
                "trait_name": "pollination_system",
                "trait_value": "bee",
                "source_record_id": "all4",
            },
            {
                "taxon_name": "Campanula microdonta",
                "trait_name": "sex_type",
                "trait_value": "unknown",
                "source_record_id": "all5",
            },
        ]
    ).to_csv(source, index=False)

    common = {
        "master_csv": master,
        "source_csv": source,
        "profiles_path": Path("config/bulk_trait_sources.yml"),
        "source_name": "austraits_all_v7.0.0",
        "source_citation_override": "AusTraits v7.0.0",
        "source_url_override": "https://doi.org/10.5281/zenodo.15718081",
    }
    primary_report = ingest_bulk_traits(
        **common,
        output_dir=tmp_path / "primary",
        profile_name="austraits_all_primary",
    )
    guild_report = ingest_bulk_traits(
        **common,
        output_dir=tmp_path / "guild",
        profile_name="austraits_all_pollinator_guild",
    )

    primary = pd.read_csv(primary_report["candidate_csv"], dtype=str).fillna("")
    assert set(zip(primary["trait_name"], primary["standardized_value"], strict=True)) == {
        ("flower_primary_color", "multicolored_variable"),
        ("sex_system", "dioecious"),
        ("pollen_vector_mode", "mixed"),
        ("pollen_vector_mode", "biotic"),
    }
    primary_unmapped = pd.read_csv(primary_report["unmapped_values_csv"], dtype=str).fillna("")
    assert primary_unmapped.loc[
        primary_unmapped["source_record_id"].eq("all5"), "reason"
    ].item() == ("unmapped_or_invalid_trait_value")

    guild = pd.read_csv(guild_report["candidate_csv"], dtype=str).fillna("")
    assert set(zip(guild["trait_name"], guild["standardized_value"], strict=True)) == {
        ("pollination_functional_guild", "mixed_or_generalist"),
        ("pollination_functional_guild", "other_bees"),
    }
