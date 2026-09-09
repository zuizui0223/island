from __future__ import annotations

import geopandas as gpd
import pandas as pd
from shapely.geometry import box

from island_v2.chapter1_pr138_realm_sensitivity import assign_biogeographic_realms


def test_realm_assignment_is_point_intersection_only_and_fail_closed():
    covariates = pd.DataFrame(
        [
            {"island_id": "a", "island_latitude": 0.0, "island_longitude": 0.0},
            {"island_id": "b", "island_latitude": 0.0, "island_longitude": 10.0},
            {"island_id": "c", "island_latitude": 50.0, "island_longitude": 50.0},
        ]
    )
    realms = gpd.GeoDataFrame(
        {
            "REALM": ["Nearctic", "Neotropical"],
            "geometry": [box(-2, -2, 2, 2), box(8, -2, 12, 2)],
        },
        crs="EPSG:4326",
    )
    out = assign_biogeographic_realms(covariates, realms).set_index("island_id")
    assert out.loc["a", "biogeographic_realm"] == "Nearctic"
    assert out.loc["b", "biogeographic_realm"] == "Neotropical"
    assert out.loc["c", "biogeographic_realm"] == ""
    assert out.loc["c", "realm_assignment_status"] == "unresolved_no_intersection"


def test_realm_aliases_are_normalized():
    covariates = pd.DataFrame(
        [{"island_id": "x", "island_latitude": 0.0, "island_longitude": 0.0}]
    )
    realms = gpd.GeoDataFrame(
        {"realm": ["Palaearctic"], "geometry": [box(-2, -2, 2, 2)]},
        crs="EPSG:4326",
    )
    out = assign_biogeographic_realms(covariates, realms)
    assert out.iloc[0]["biogeographic_realm"] == "Palearctic"
