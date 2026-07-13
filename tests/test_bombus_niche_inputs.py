from __future__ import annotations

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point

from island_v2.bombus_niche_inputs import ENVIRONMENT_COLUMNS, build_global_targets


def test_environment_set_covers_mean_seasonality_extremes_and_elevation() -> None:
    assert ENVIRONMENT_COLUMNS == [
        "bio1",
        "bio4",
        "bio5",
        "bio6",
        "bio12",
        "bio15",
        "bio18",
        "bio19",
        "elevation",
    ]


def test_global_targets_cross_every_island_with_every_species() -> None:
    islands = gpd.GeoDataFrame(
        {
            "island_id": ["island_b", "island_a", "island_a"],
            "geometry": [Point(1, 1), Point(0, 0), Point(0, 0)],
        },
        crs=4326,
    )
    species = pd.DataFrame(
        {
            "bombus_species": ["Bombus beta", "Bombus alpha"],
            "species_key": [2, 1],
        }
    )

    targets = build_global_targets(islands, species, "island_id")

    assert len(targets) == 4
    assert set(zip(targets["island_id"], targets["bombus_species"], strict=True)) == {
        ("island_a", "Bombus alpha"),
        ("island_a", "Bombus beta"),
        ("island_b", "Bombus alpha"),
        ("island_b", "Bombus beta"),
    }
    assert targets["projection_scope"].eq(
        "global_environmental_screen_not_source_pool"
    ).all()
