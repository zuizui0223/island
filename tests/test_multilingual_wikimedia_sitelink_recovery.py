from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from island_v2.multilingual_wikimedia_recovery import RECOVERY_TASK
from island_v2.multilingual_wikimedia_sitelink_recovery import (
    WAVE_TASK,
    fetch_english_wikidata_ids,
    fetch_resolved_language_extracts,
    fetch_wikidata_sitelinks,
    run_sitelink_recovery_wave,
)


def test_fetch_english_qids_follows_redirects() -> None:
    def getter(url: str, params: dict) -> dict:
        assert "en.wikipedia.org" in url
        assert params["prop"] == "pageprops"
        return {
            "query": {
                "redirects": [
                    {"from": "Campanula microdonta", "to": "Campanula punctata"}
                ],
                "pages": {
                    "1": {
                        "title": "Campanula punctata",
                        "pageprops": {"wikibase_item": "Q123"},
                    },
                    "-1": {"title": "Missing plant", "missing": ""},
                },
            }
        }

    qids, completed, errors = fetch_english_wikidata_ids(
        getter,
        ["Campanula microdonta", "Missing plant"],
    )

    assert qids == {"Campanula microdonta": "Q123"}
    assert completed == {"Campanula microdonta", "Missing plant"}
    assert errors == []


def test_fetch_wikidata_sitelinks_returns_local_titles() -> None:
    def getter(url: str, params: dict) -> dict:
        assert "wikidata.org" in url
        assert params["sitefilter"] == "jawiki|zhwiki|ruwiki"
        return {
            "entities": {
                "Q123": {
                    "sitelinks": {
                        "jawiki": {"title": "シマホタルブクロ"},
                        "zhwiki": {"title": "紫斑風鈴草"},
                    }
                }
            }
        }

    resolved, errors = fetch_wikidata_sitelinks(
        getter,
        {"Campanula microdonta": "Q123"},
        ("ja", "zh", "ru"),
    )

    assert resolved == {
        "Campanula microdonta": {
            "ja": "シマホタルブクロ",
            "zh": "紫斑風鈴草",
        }
    }
    assert errors == []


def test_fetch_resolved_extract_retains_scientific_name() -> None:
    def getter(url: str, params: dict) -> dict:
        assert "ja.wikipedia.org" in url
        assert params["titles"] == "シマホタルブクロ"
        return {
            "query": {
                "pages": {
                    "1": {
                        "title": "シマホタルブクロ",
                        "extract": "送粉者はマルハナバチである。",
                    }
                }
            }
        }

    sources, completed, errors = fetch_resolved_language_extracts(
        getter,
        {"Campanula microdonta": "シマホタルブクロ"},
        "ja",
    )

    assert completed == {"Campanula microdonta"}
    assert errors == []
    assert sources["accepted_species"].tolist() == ["Campanula microdonta"]
    assert sources["source_type"].tolist() == ["wikipedia_sitelink_extract_ja"]
    assert "シマホタルブクロ" in sources["source_url"].item()


def test_sitelink_wave_recovers_vernacular_pages_and_reopens_downstream(
    tmp_path: Path,
) -> None:
    campaign_dir = tmp_path / "campaign"
    campaign_dir.mkdir()
    pd.DataFrame(
        [
            {
                "accepted_species": "Campanula microdonta",
                "genus": "Campanula",
                "family": "Campanulaceae",
                "n_islands": "2",
                "n_records": "10",
                "machine_biotic_candidate": "False",
                "reproductive_wikimedia_status": "processed",
                "floral_access_wikimedia_status": "not_applicable",
                "alternative_guild_wikimedia_status": "not_applicable",
                "alternative_guild_openalex_status": "not_applicable",
            },
            {
                "accepted_species": "Poa annua",
                "genus": "Poa",
                "family": "Poaceae",
                "n_islands": "4",
                "n_records": "20",
                "machine_biotic_candidate": "False",
                "reproductive_wikimedia_status": "pending",
                "floral_access_wikimedia_status": "pending",
                "alternative_guild_wikimedia_status": "pending",
                "alternative_guild_openalex_status": "pending",
            },
        ]
    ).to_csv(campaign_dir / "campaign_ledger.csv.gz", index=False, compression="gzip")
    (campaign_dir / "campaign_status.json").write_text(
        json.dumps({"n_machine_biotic_candidates": 0}) + "\n",
        encoding="utf-8",
    )

    def getter(url: str, params: dict) -> dict:
        if "en.wikipedia.org" in url and params.get("prop") == "pageprops":
            return {
                "query": {
                    "pages": {
                        "1": {
                            "title": "Campanula microdonta",
                            "pageprops": {"wikibase_item": "Q123"},
                        }
                    }
                }
            }
        if "wikidata.org" in url:
            return {
                "entities": {
                    "Q123": {
                        "sitelinks": {
                            "jawiki": {"title": "シマホタルブクロ"},
                            "zhwiki": {"title": "紫斑風鈴草"},
                        }
                    }
                }
            }
        if "ja.wikipedia.org" in url and params["titles"] == "シマホタルブクロ":
            return {
                "query": {
                    "pages": {
                        "1": {
                            "title": "シマホタルブクロ",
                            "extract": "本種の送粉者はマルハナバチである。",
                        }
                    }
                }
            }
        if "zh.wikipedia.org" in url and params["titles"] == "紫斑風鈴草":
            return {
                "query": {
                    "pages": {
                        "1": {
                            "title": "紫斑風鈴草",
                            "extract": "该植物为雌雄同株。",
                        }
                    }
                }
            }
        requested = params.get("titles", "")
        return {"query": {"pages": {"-1": {"title": requested, "missing": ""}}}}

    summary = run_sitelink_recovery_wave(
        campaign_dir=campaign_dir,
        batch_size=100,
        languages=("ja", "zh", "ru"),
        getter=getter,
    )

    assert summary["task"] == WAVE_TASK
    assert summary["n_species_attempted"] == 1
    assert summary["n_species_with_english_qid"] == 1
    assert summary["n_species_with_any_target_sitelink"] == 1
    assert summary["n_species_with_sitelink_source_text"] == 1
    assert summary["n_candidates_total"] == 3
    assert summary["n_new_biotic_species"] == 1

    ledger = pd.read_csv(campaign_dir / "campaign_ledger.csv.gz", dtype=str).fillna("")
    campanula = ledger.loc[
        ledger["accepted_species"].eq("Campanula microdonta")
    ].iloc[0]
    assert campanula["machine_biotic_candidate"].lower() == "true"
    assert campanula["floral_access_wikimedia_status"] == "pending"
    assert campanula["alternative_guild_wikimedia_status"] == "pending"
    assert campanula[f"{WAVE_TASK}_status"] == "processed"

    wave_dir = next(
        (campaign_dir / "waves").glob(
            "wave_*_reproductive_multilingual_sitelink_recovery"
        )
    )
    candidates = pd.read_csv(
        wave_dir / "machine_candidates.csv.gz",
        dtype=str,
    ).fillna("")
    assert set(candidates["campaign_task"]) == {RECOVERY_TASK}
    observed = set(
        candidates[["trait_name", "candidate_value"]].itertuples(
            index=False,
            name=None,
        )
    )
    assert ("pollination_functional_guild", "bumblebees") in observed
    assert ("pollen_vector_mode", "biotic") in observed
    assert ("sex_system", "monoecious") in observed
