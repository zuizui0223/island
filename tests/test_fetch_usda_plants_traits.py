from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from island_v2.fetch_usda_plants_traits import (
    FILTERS_URL,
    PROFILE_ROOT,
    RESULTS_URL,
    fetch_usda_plants_traits,
    map_usda_flower_colors,
    prepare_usda_rows,
)


def _plants() -> list[dict[str, object]]:
    return [
        {
            "id": 10,
            "symbol": "ALON",
            "scientificNameWithoutAuthor": "Alpha one",
            "isSynonym": False,
        },
        {
            "id": 20,
            "symbol": "BETW",
            "scientificNameWithoutAuthor": "Beta two",
            "isSynonym": False,
        },
        {
            "id": 30,
            "symbol": "GAMV",
            "scientificNameWithoutAuthor": "Gamma three var. minor",
            "isSynonym": False,
        },
    ]


def _filters() -> list[dict[str, object]]:
    return [
        {
            "name": "Flower Color",
            "type": "Morphology",
            "filters": [
                {"id": 1, "display": "Blue", "plantIds": [10, 11]},
                {"id": 2, "display": "White", "plantIds": [11, 20]},
                {"id": 3, "display": "Yellow", "plantIds": [30, 40]},
            ],
        },
        {
            "name": "Group",
            "type": "Basic",
            "filters": [
                {"id": 1, "display": "Dicot", "plantIds": [10, 11, 30]},
                {"id": 2, "display": "Monocot", "plantIds": [20]},
                {"id": 3, "display": "Gymnosperm", "plantIds": [40]},
            ],
        },
    ]


def _resolution() -> list[dict[str, object]]:
    return [
        {
            "requested_id": 11,
            "url": f"{PROFILE_ROOT}/11",
            "response": {"Id": 10, "Rank": "Species", "Group": "Dicot"},
        }
    ]


def test_map_usda_flower_colors_is_conservative() -> None:
    assert map_usda_flower_colors(["Blue", "Purple"]) == "blue_purple"
    assert map_usda_flower_colors(["Blue", "White"]) == "multicolored_variable"
    assert map_usda_flower_colors(["Yellow"]) == "yellow_orange"
    assert map_usda_flower_colors(["Chartreuse"]) is None


def test_prepare_usda_rows_resolves_aliases_and_audits_exclusions() -> None:
    prepared, audit = prepare_usda_rows(
        _plants(),
        _filters(),
        profile_resolutions=_resolution(),
        accessed_on="2026-07-14",
    )

    assert prepared[["scientific_name", "trait_value"]].to_dict("records") == [
        {"scientific_name": "Alpha one", "trait_value": "multicolored_variable"},
        {"scientific_name": "Beta two", "trait_value": "white"},
    ]
    alpha = prepared.loc[prepared["scientific_name"].eq("Alpha one")].iloc[0]
    evidence = json.loads(alpha["evidence_excerpt"])
    assert evidence == {
        "characteristic_name": "Flower Color",
        "raw_values": ["Blue", "White"],
        "source_plant_ids": [10, 11],
        "resolved_plant_id": 10,
        "symbol": "ALON",
        "scientificNameWithoutAuthor": "Alpha one",
        "source_groups": ["Dicot"],
        "filter_url": FILTERS_URL,
        "profile_resolution_urls": [f"{PROFILE_ROOT}/11"],
    }
    assert "Accessed 2026-07-14" in alpha["source_citation"]
    assert alpha["source_url"] == FILTERS_URL
    assert set(audit["reason"]) == {
        "non_angiosperm_source_group",
        "non_binomial_or_hybrid_taxon",
    }


def test_fetch_usda_plants_traits_snapshots_all_api_inputs(tmp_path: Path) -> None:
    payloads = {
        RESULTS_URL: _plants(),
        FILTERS_URL: _filters(),
        f"{PROFILE_ROOT}/11": _resolution()[0]["response"],
        f"{PROFILE_ROOT}/40": {"Id": 40, "Rank": "Species", "Group": "Gymnosperm"},
    }
    calls: list[str] = []

    def getter(url: str) -> object:
        calls.append(url)
        return payloads[url]

    manifest = fetch_usda_plants_traits(
        output_dir=tmp_path,
        max_workers=1,
        getter=getter,
        accessed_on="2026-07-14",
    )

    assert calls[:2] == [RESULTS_URL, FILTERS_URL]
    assert set(calls[2:]) == {f"{PROFILE_ROOT}/11", f"{PROFILE_ROOT}/40"}
    assert manifest["n_flower_color_source_ids"] == 5
    assert manifest["n_missing_ids_resolved_by_profile"] == 2
    assert manifest["n_prepared_canonical_rows"] == 2
    assert (tmp_path / "metadata" / "characteristicSearchResults.json").exists()
    assert (tmp_path / "metadata" / "characteristicSearchFilterResults.json").exists()
    assert (tmp_path / "metadata" / "missingPlantProfileResolutions.json").exists()
    assert (tmp_path / "usda_plants_fetch_manifest.json").exists()
    prepared = pd.read_csv(tmp_path / "usda_plants_prepared_canonical_long.csv.gz")
    assert prepared["scientific_name"].tolist() == ["Alpha one", "Beta two"]
