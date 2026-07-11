from __future__ import annotations

import pandas as pd

from island_v2.trait_media_scout import collect_media


def test_collect_media_deduplicates_images_and_preserves_provenance() -> None:
    def getter(url: str, params: dict[str, object]) -> dict[str, object]:
        if url.endswith("/species/match"):
            return {"usageKey": 123}
        if url.endswith("/occurrence/search"):
            return {
                "results": [
                    {
                        "key": 1,
                        "scientificName": "Plantus alba",
                        "basisOfRecord": "HUMAN_OBSERVATION",
                        "countryCode": "JP",
                        "eventDate": "2026-01-01",
                        "media": [
                            {
                                "identifier": "https://example.org/a.jpg",
                                "type": "StillImage",
                                "creator": ["A. Collector", "B. Photographer"],
                                "license": {"name": "CC-BY"},
                            },
                            {"identifier": "https://example.org/a.jpg", "type": "StillImage"},
                            {"identifier": "https://example.org/b.jpg", "type": "StillImage"},
                        ],
                    }
                ]
            }
        raise AssertionError(url)

    frame, report = collect_media(
        pd.DataFrame({"accepted_species": ["Plantus alba"]}),
        getter,
        max_taxa=10,
        max_occurrences=50,
        max_images_per_species=5,
    )

    assert frame["media_identifier"].tolist() == [
        "https://example.org/a.jpg",
        "https://example.org/b.jpg",
    ]
    assert frame["gbif_taxon_key"].tolist() == ["123", "123"]
    assert frame.iloc[0]["media_creator"] == "A. Collector | B. Photographer"
    assert '"name": "CC-BY"' in frame.iloc[0]["media_license"]
    assert report["n_species_with_media"] == 1
    assert report["n_media_rows"] == 2


def test_collect_media_prefers_distinct_occurrences_before_second_views() -> None:
    def getter(url: str, params: dict[str, object]) -> dict[str, object]:
        if url.endswith("/species/match"):
            return {"usageKey": 123}
        if url.endswith("/occurrence/search"):
            return {
                "results": [
                    {
                        "key": 1,
                        "media": [
                            {"identifier": "https://example.org/1a.jpg", "type": "StillImage"},
                            {"identifier": "https://example.org/1b.jpg", "type": "StillImage"},
                        ],
                    },
                    {
                        "key": 2,
                        "media": [
                            {"identifier": "https://example.org/2a.jpg", "type": "StillImage"},
                            {"identifier": "https://example.org/2b.jpg", "type": "StillImage"},
                        ],
                    },
                    {
                        "key": 3,
                        "media": [
                            {"identifier": "https://example.org/3a.jpg", "type": "StillImage"}
                        ],
                    },
                ]
            }
        raise AssertionError(url)

    frame, report = collect_media(
        pd.DataFrame({"accepted_species": ["Plantus alba"]}),
        getter,
        max_taxa=10,
        max_occurrences=50,
        max_images_per_species=4,
    )

    assert frame["occurrence_key"].tolist() == ["1", "2", "3", "1"]
    assert frame["media_identifier"].tolist() == [
        "https://example.org/1a.jpg",
        "https://example.org/2a.jpg",
        "https://example.org/3a.jpg",
        "https://example.org/1b.jpg",
    ]
    assert report["selection_policy"] == (
        "round_robin_distinct_occurrences_before_additional_images"
    )
