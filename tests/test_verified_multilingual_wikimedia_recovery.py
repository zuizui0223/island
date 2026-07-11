from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from island_v2.verified_multilingual_wikimedia_recovery import (
    WAVE_TASK,
    fetch_verified_sitelinks,
    filter_uncertain_candidates,
    run_verified_recovery_wave,
)


def test_fetch_verified_sitelinks_accepts_only_exact_p225() -> None:
    def getter(url: str, params: dict) -> dict:
        assert "wikidata.org" in url
        assert params["props"] == "claims|sitelinks"
        return {
            "entities": {
                "Q1": {
                    "claims": {
                        "P225": [
                            {
                                "mainsnak": {
                                    "datavalue": {
                                        "value": "Campanula microdonta"
                                    }
                                }
                            }
                        ]
                    },
                    "sitelinks": {
                        "jawiki": {"title": "シマホタルブクロ"},
                    },
                },
                "Q2": {
                    "claims": {
                        "P225": [
                            {
                                "mainsnak": {
                                    "datavalue": {"value": "Glossopteris"}
                                }
                            }
                        ]
                    },
                    "sitelinks": {
                        "ruwiki": {"title": "Глоссоптерис"},
                    },
                },
                "Q3": {
                    "claims": {},
                    "sitelinks": {
                        "zhwiki": {"title": "未知植物"},
                    },
                },
            }
        }

    resolved, audit, errors = fetch_verified_sitelinks(
        getter,
        {
            "Campanula microdonta": "Q1",
            "Glossopteris indica": "Q2",
            "Unknowna example": "Q3",
        },
        ("ja", "zh", "ru"),
    )

    assert resolved == {
        "Campanula microdonta": {"ja": "シマホタルブクロ"}
    }
    assert errors == []
    statuses = dict(
        audit[["accepted_species", "verification_status"]].itertuples(
            index=False,
            name=None,
        )
    )
    assert statuses == {
        "Campanula microdonta": "p225_verified",
        "Glossopteris indica": "p225_mismatch",
        "Unknowna example": "missing_p225",
    }


def test_filter_uncertain_candidates_rejects_unknown_russian_state() -> None:
    candidates = pd.DataFrame(
        [
            {
                "accepted_species": "Glossopteris indica",
                "trait_name": "sex_system",
                "candidate_value": "monoecious",
                "source_type": "wikipedia_sitelink_extract_ru",
                "source_excerpt": (
                    "Неизвестно, были ли растения однодомными или двудомными."
                ),
            },
            {
                "accepted_species": "Aucuba japonica",
                "trait_name": "sex_system",
                "candidate_value": "dioecious",
                "source_type": "wikipedia_sitelink_extract_ja",
                "source_excerpt": "雌雄異株である。",
            },
        ]
    )

    retained, rejected = filter_uncertain_candidates(candidates)

    assert retained["accepted_species"].tolist() == ["Aucuba japonica"]
    assert rejected["accepted_species"].tolist() == ["Glossopteris indica"]
    assert rejected["audit_reason"].tolist() == [
        "explicit_uncertainty_context"
    ]


def test_verified_wave_uses_twenty_title_extract_batches(tmp_path: Path) -> None:
    campaign_dir = tmp_path / "campaign"
    campaign_dir.mkdir()
    species = [f"Examplea species{i}" for i in range(45)]
    ledger = pd.DataFrame(
        [
            {
                "accepted_species": name,
                "genus": "Examplea",
                "family": "Exampleaceae",
                "n_islands": "1",
                "n_records": str(100 - index),
                "machine_biotic_candidate": "False",
                "reproductive_wikimedia_status": "processed",
                "floral_access_wikimedia_status": "not_applicable",
                "alternative_guild_wikimedia_status": "not_applicable",
                "alternative_guild_openalex_status": "not_applicable",
            }
            for index, name in enumerate(species)
        ]
    )
    ledger.to_csv(
        campaign_dir / "campaign_ledger.csv.gz",
        index=False,
        compression="gzip",
    )
    (campaign_dir / "campaign_status.json").write_text(
        json.dumps({"n_machine_biotic_candidates": 0}) + "\n",
        encoding="utf-8",
    )

    extract_batch_sizes: list[int] = []

    def getter(url: str, params: dict) -> dict:
        if "en.wikipedia.org" in url:
            pages = {
                str(index): {
                    "title": name,
                    "pageprops": {"wikibase_item": f"Q{index}"},
                }
                for index, name in enumerate(species)
            }
            return {"query": {"pages": pages}}
        if "wikidata.org" in url:
            entities = {}
            for qid in params["ids"].split("|"):
                index = int(qid[1:])
                name = species[index]
                entities[qid] = {
                    "claims": {
                        "P225": [
                            {
                                "mainsnak": {
                                    "datavalue": {"value": name}
                                }
                            }
                        ]
                    },
                    "sitelinks": {
                        "jawiki": {"title": f"現地名{index}"},
                    },
                }
            return {"entities": entities}
        if "ja.wikipedia.org" in url:
            requested = params["titles"].split("|")
            extract_batch_sizes.append(len(requested))
            pages = {
                str(index): {
                    "title": title,
                    "extract": "雌雄異株である。",
                }
                for index, title in enumerate(requested)
            }
            return {"query": {"pages": pages}}
        return {"query": {"pages": {}}}

    summary = run_verified_recovery_wave(
        campaign_dir=campaign_dir,
        batch_size=45,
        languages=("ja",),
        getter=getter,
    )

    assert extract_batch_sizes == [20, 20, 5]
    assert summary["n_species_attempted"] == 45
    assert summary["wikidata_identity_status_counts"] == {
        "p225_verified": 45
    }
    assert summary["n_species_with_verified_target_sitelink"] == 45
    assert summary["n_species_with_source_text"] == 45
    assert summary["n_candidates_total"] == 45
    assert summary["n_uncertain_candidates_rejected"] == 0

    wave_dir = next(
        (campaign_dir / "waves").glob(
            "wave_*_reproductive_multilingual_verified_sitelink_recovery"
        )
    )
    assert wave_dir.name.endswith(WAVE_TASK)
    identity = pd.read_csv(
        wave_dir / "wikidata_identity_audit.csv.gz",
        dtype=str,
    ).fillna("")
    assert set(identity["verification_status"]) == {"p225_verified"}
