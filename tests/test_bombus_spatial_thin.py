from __future__ import annotations

import pandas as pd
import pytest

from island_v2.bombus_spatial_thin import spatial_thin


def test_spatial_thin_keeps_one_record_per_species_cell() -> None:
    frame = pd.DataFrame(
        {
            "bombus_species": ["Bombus alpha", "Bombus alpha", "Bombus alpha", "Bombus beta"],
            "decimal_latitude": [35.01, 35.09, 35.40, 35.02],
            "decimal_longitude": [139.01, 139.08, 139.40, 139.02],
            "year": [2020, 2021, 2019, 2022],
            "gbif_id": ["2", "1", "3", "4"],
        }
    )

    result = spatial_thin(frame, cell_size_degrees=0.25)

    assert len(result) == 3
    assert result.groupby("bombus_species").size().to_dict() == {
        "Bombus alpha": 2,
        "Bombus beta": 1,
    }
    assert not result.duplicated(
        ["bombus_species", "spatial_thin_lat_cell", "spatial_thin_lon_cell"]
    ).any()


def test_spatial_thin_rejects_nonpositive_cell_size() -> None:
    frame = pd.DataFrame(
        {
            "bombus_species": ["Bombus alpha"],
            "decimal_latitude": [0.0],
            "decimal_longitude": [0.0],
        }
    )
    with pytest.raises(Exception, match="must be positive"):
        spatial_thin(frame, cell_size_degrees=0)
