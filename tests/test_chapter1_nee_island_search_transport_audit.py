from __future__ import annotations

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import box

from island_v2.chapter1_nee_island_search_transport_audit import (
    _error_class,
    audit_transport,
    prepare_audit_frame,
)


def _config() -> dict:
    return {
        "contract": "chapter1_nee_island_search_transport_audit_v1",
        "status": "reporting_only_after_canonical_search",
        "inputs": {"channels": ["bombus", "diptera"]},
    }


def _observations() -> pd.DataFrame:
    rows = []
    for channel in ["bombus", "diptera"]:
        for index in range(10):
            error = "RuntimeError: 400 Bad Request" if index >= 8 else ""
            rows.append(
                {
                    "island_id": f"i{index}",
                    "channel_id": channel,
                    "observation_state": "unresolved" if error else "insufficient_effort",
                    "search_complete": "false" if error else "true",
                    "search_truncated": "false",
                    "search_error": error,
                }
            )
    return pd.DataFrame(rows)


def _covariates() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "island_id": [f"i{i}" for i in range(10)],
            "area_km2": [float(i + 1) for i in range(10)],
            "distance_to_continent_km": [float((i + 1) * 2) for i in range(10)],
            "log_island_area_km2": [float(i) / 10 for i in range(10)],
            "log_distance_to_continent_km": [float(i) / 20 for i in range(10)],
        }
    )


def _islands() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        {
            "island_id": [f"i{i}" for i in range(10)],
            "geometry": [box(i, 0, i + 0.5, 0.5) for i in range(10)],
        },
        geometry="geometry",
        crs="EPSG:4326",
    )


def test_transport_audit_is_reporting_only_and_counts_errors() -> None:
    summary, deciles, receipt = audit_transport(
        _observations(), _covariates(), _islands(), _config()
    )
    assert len(summary) == 2
    assert set(summary["channel_id"]) == {"bombus", "diptera"}
    assert summary.set_index("channel_id").loc["bombus", "n_search_errors"] == 2
    assert summary.set_index("channel_id").loc["bombus", "search_error_fraction"] == pytest.approx(0.2)
    assert summary.set_index("channel_id").loc["bombus", "n_error_http_400"] == 2
    assert receipt["reporting_only"] is True
    assert receipt["transport_audit_cannot_change_N1_pass"] is True
    assert receipt["rerun_or_recode_permitted"] is False
    assert set(deciles["diagnostic"]) == {
        "exact_polygon_wkt_length",
        "log1p_area_km2",
        "log1p_distance_to_continent_km",
    }


def test_error_rows_are_not_recoded_by_audit() -> None:
    obs = _observations()
    original = obs.copy(deep=True)
    prepare_audit_frame(obs, _covariates(), _islands())
    pd.testing.assert_frame_equal(obs, original)


def test_duplicate_observation_rows_are_rejected() -> None:
    obs = pd.concat([_observations(), _observations().iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="unique by island_id x channel_id"):
        prepare_audit_frame(obs, _covariates(), _islands())


def test_error_classification() -> None:
    assert _error_class("") == "none"
    assert _error_class("RuntimeError: 400 Bad Request") == "http_400"
    assert _error_class("HTTPStatusError: 429 Too Many Requests") == "http_429"
    assert _error_class("URL query is too long") == "query_too_long"
    assert _error_class("RemoteProtocolError: Server disconnected") == "server_disconnect"
    assert _error_class("ReadTimeout: timed out") == "timeout"
