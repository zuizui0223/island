from __future__ import annotations

from io import BytesIO

import pandas as pd
from PIL import Image

from island_v2.trait_visual_packet import build_contact_sheets, freeze_images


def _jpeg_bytes() -> bytes:
    buffer = BytesIO()
    Image.new("RGB", (40, 30), "white").save(buffer, format="JPEG")
    return buffer.getvalue()


def test_freeze_images_hashes_and_builds_contact_sheet(tmp_path) -> None:
    media = pd.DataFrame(
        [
            {
                "accepted_species": "Plantus alba",
                "media_identifier": "https://example.org/a.jpg",
                "occurrence_key": "1",
                "media_license": "CC-BY",
            },
            {
                "accepted_species": "Plantus alba",
                "media_identifier": "https://example.org/b.jpg",
                "occurrence_key": "2",
                "media_license": "CC-BY",
            },
        ]
    )

    frozen, report = freeze_images(
        media,
        tmp_path,
        images_per_species=5,
        downloader=lambda _: _jpeg_bytes(),
    )

    assert report["n_species_with_downloaded_image"] == 1
    assert report["n_images_downloaded"] == 2
    assert frozen["image_sha256"].str.len().eq(64).all()
    sheets = build_contact_sheets(frozen, tmp_path, species_per_sheet=5)
    assert len(sheets) == 1
    assert sheets[0].exists()
