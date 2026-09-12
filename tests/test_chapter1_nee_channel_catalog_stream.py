from __future__ import annotations

from pathlib import Path

import pandas as pd

from island_v2.chapter1_nee_channel_catalog import build_catalog, load_config
from island_v2.chapter1_nee_channel_catalog_stream import CatalogAccumulator


CONFIG = Path("config/chapter1_nee_channel_taxon_catalog.yml")


def _evidence() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "channel_id": "lepidoptera",
                "pollinator_species": "Danaus plexippus",
                "pollinator_genus": "Danaus",
                "pollinator_family": "Nymphalidae",
                "pollinator_order": "Lepidoptera",
                "pollinator_class": "Insecta",
                "plant_taxon": "Plant alpha",
                "interaction_type_id": "RO:visit",
                "interaction_type_name": "visitsFlowersOf",
                "evidence_strength": "flower_visit",
                "reference_key": "citation:A",
                "reference_citation": "A",
                "reference_doi": "",
                "reference_url": "",
                "source_namespace": "mock",
                "source_archive_uri": "",
                "source_doi": "",
            },
            {
                "channel_id": "lepidoptera",
                "pollinator_species": "Danaus plexippus",
                "pollinator_genus": "Danaus",
                "pollinator_family": "Nymphalidae",
                "pollinator_order": "Lepidoptera",
                "pollinator_class": "Insecta",
                "plant_taxon": "Plant beta",
                "interaction_type_id": "RO:visit",
                "interaction_type_name": "visitsFlowersOf",
                "evidence_strength": "flower_visit",
                "reference_key": "citation:B",
                "reference_citation": "B",
                "reference_doi": "",
                "reference_url": "",
                "source_namespace": "mock",
                "source_archive_uri": "",
                "source_doi": "",
            },
            {
                "channel_id": "bombus",
                "pollinator_species": "Bombus terrestris",
                "pollinator_genus": "Bombus",
                "pollinator_family": "Apidae",
                "pollinator_order": "Hymenoptera",
                "pollinator_class": "Insecta",
                "plant_taxon": "Plant gamma",
                "interaction_type_id": "RO:pollinates",
                "interaction_type_name": "pollinates",
                "evidence_strength": "pollination",
                "reference_key": "citation:C",
                "reference_citation": "C",
                "reference_doi": "",
                "reference_url": "",
                "source_namespace": "mock",
                "source_archive_uri": "",
                "source_doi": "",
            },
        ]
    )


def test_streaming_accumulator_matches_in_memory_catalog() -> None:
    config = load_config(CONFIG)
    evidence = _evidence()
    expected = build_catalog(
        evidence,
        config,
        source_version_doi="10.5281/zenodo.20546682",
        source_sha256="abc",
    ).sort_values(["channel_id", "pollinator_species"]).reset_index(drop=True)

    accumulator = CatalogAccumulator()
    accumulator.add(evidence.iloc[:1])
    accumulator.add(evidence.iloc[1:])
    actual = accumulator.to_frame(
        config,
        source_version_doi="10.5281/zenodo.20546682",
        source_sha256="abc",
    ).sort_values(["channel_id", "pollinator_species"]).reset_index(drop=True)

    pd.testing.assert_frame_equal(actual, expected)


def test_streaming_accumulator_keeps_visit_only_singletons_sensitivity_only() -> None:
    config = load_config(CONFIG)
    accumulator = CatalogAccumulator()
    accumulator.add(_evidence().iloc[:1])
    catalog = accumulator.to_frame(
        config,
        source_version_doi="10.5281/zenodo.20546682",
        source_sha256="abc",
    )
    assert catalog.loc[0, "catalog_tier"] == "sensitivity"
