from __future__ import annotations

import pandas as pd

from island_v2.trait_media_scout import collect_media


def test_collect_media_deduplicates_images_and_preserves_provenance() -> None:
    def getter(url: str, params: dict[str, object]) -> dict[str, object]:
        if url.endswith('/species/match'):
            return {'usageKey': 123}
        if url.endswith('/occurrence/search'):
            return {
                'results': [
                    {
                        'key': 1,
                        'scientificName': 'Plantus alba',
                        'basisOfRecord': 'HUMAN_OBSERVATION',
                        'countryCode': 'JP',
                        'eventDate': '2026-01-01',
                        'media': [
                            {'identifier': 'https://example.org/a.jpg', 'type': 'StillImage'},
                            {'identifier': 'https://example.org/a.jpg', 'type': 'StillImage'},
                            {'identifier': 'https://example.org/b.jpg', 'type': 'StillImage'},
                        ],
                    }
                ]
            }
        raise AssertionError(url)

    frame, report = collect_media(
        pd.DataFrame({'accepted_species': ['Plantus alba']}),
        getter,
        max_taxa=10,
        max_occurrences=50,
        max_images_per_species=5,
    )

    assert frame['media_identifier'].tolist() == [
        'https://example.org/a.jpg',
        'https://example.org/b.jpg',
    ]
    assert frame['gbif_taxon_key'].tolist() == ['123', '123']
    assert report['n_species_with_media'] == 1
    assert report['n_media_rows'] == 2
