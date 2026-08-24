import importlib

import geopandas as gpd
import pandas as pd
from shapely.geometry import box

mod = importlib.import_module("island_v2.gift_floristic_status_match")


def test_spatial_match_accepts_strong_unambiguous_exact_island():
    gift = gpd.GeoDataFrame(
        {
            "entity_ID": ["101"],
            "entity_class": ["Island"],
            "geometry": [box(0, 0, 1, 1)],
        },
        crs="EPSG:4326",
    )
    gshhg = gpd.GeoDataFrame(
        {
            "island_id": ["i1", "i2"],
            "geometry": [box(0, 0, 1, 1), box(3, 3, 4, 4)],
        },
        crs="EPSG:4326",
    )
    matches, candidates = mod.match_gift_islands_to_gshhg(gift, gshhg)
    assert matches.iloc[0]["spatial_match_status"] == "accepted"
    assert matches.iloc[0]["island_id"] == "i1"
    assert len(candidates) == 1


def test_spatial_match_rejects_island_group_and_weak_overlap():
    gift = gpd.GeoDataFrame(
        {
            "entity_ID": ["group", "weak"],
            "entity_class": ["Island Group", "Island"],
            "geometry": [box(0, 0, 2, 1), box(10, 10, 12, 12)],
        },
        crs="EPSG:4326",
    )
    gshhg = gpd.GeoDataFrame(
        {
            "island_id": ["i1", "i2"],
            "geometry": [box(0, 0, 1, 1), box(10, 10, 10.5, 10.5)],
        },
        crs="EPSG:4326",
    )
    matches, _ = mod.match_gift_islands_to_gshhg(gift, gshhg)
    assert set(matches["entity_ID"]) == {"weak"}
    assert matches.iloc[0]["spatial_match_status"] == "insufficient_overlap"


def test_status_classification_requires_status_metadata_and_list_endemism():
    native_endemic = pd.Series(
        {
            "native_indicated": "1",
            "natural_indicated": "1",
            "native": "1",
            "naturalized": "0",
            "quest_native": "0",
            "end_list": "1",
            "endemic_list": "1",
            "quest_end_list": "0",
            "matched": "1",
            "resolved": "1",
            "questionable": "0",
            "subtaxon": "",
        }
    )
    assert mod.classify_gift_status(native_endemic)[:2] == ("native", "endemic")

    nonendemic = native_endemic.copy()
    nonendemic["endemic_list"] = "0"
    assert mod.classify_gift_status(nonendemic)[:2] == ("native", "nonendemic")

    no_endemism_scope = native_endemic.copy()
    no_endemism_scope["end_list"] = "0"
    assert mod.classify_gift_status(no_endemism_scope)[:2] == ("native", "unresolved")


def test_native_naturalized_conflict_fails_closed():
    row = pd.Series(
        {
            "native_indicated": "1",
            "natural_indicated": "1",
            "native": "1",
            "naturalized": "1",
            "quest_native": "0",
            "end_list": "1",
            "endemic_list": "0",
            "quest_end_list": "0",
            "matched": "1",
            "resolved": "1",
        }
    )
    assert mod.classify_gift_status(row)[:2] == ("unresolved", "unresolved")


def test_build_status_ledger_keeps_only_frozen_flora_memberships():
    checklists = pd.DataFrame(
        {
            "list_ID": ["11", "11", "11"],
            "work_species": ["Alpha one", "Beta two", "Outside three"],
            "native": ["1", "0", "1"],
            "quest_native": ["0", "0", "0"],
            "naturalized": ["0", "1", "0"],
            "endemic_list": ["1", "0", "0"],
            "quest_end_list": ["0", "0", "0"],
            "matched": ["1", "1", "1"],
            "resolved": ["1", "1", "1"],
            "questionable": ["0", "0", "0"],
            "subtaxon": ["", "", ""],
        }
    )
    meta = pd.DataFrame(
        {
            "list_ID": ["11"],
            "ref_ID": ["7"],
            "entity_ID": ["101"],
            "entity_class": ["Island"],
            "suit_geo": ["1"],
            "restricted": ["0"],
            "native_indicated": ["1"],
            "natural_indicated": ["1"],
            "end_list": ["1"],
            "end_ref": ["1"],
            "geo_entity": ["Test Island"],
        }
    )
    matches = pd.DataFrame(
        {
            "entity_ID": ["101"],
            "island_id": ["i1"],
            "spatial_match_status": ["accepted"],
        }
    )
    flora = pd.DataFrame(
        {"island_id": ["i1", "i1"], "species": ["Alpha one", "Beta two"]}
    )
    ledger, audit = mod.build_status_ledger(
        checklists, meta, matches, flora, gift_version="3.2"
    )
    assert set(ledger["accepted_species"]) == {"Alpha one", "Beta two"}
    alpha = ledger.loc[ledger["accepted_species"].eq("Alpha one")].iloc[0]
    beta = ledger.loc[ledger["accepted_species"].eq("Beta two")].iloc[0]
    assert (alpha["origin_status"], alpha["endemic_status"]) == ("native", "endemic")
    assert beta["origin_status"] == "introduced"
    assert not audit.loc[audit["accepted_species"].eq("Outside three"), "in_frozen_flora"].iloc[0]
