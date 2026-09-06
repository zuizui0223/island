from __future__ import annotations

import pandas as pd

from island_v2.useful_plants_reproductive_source import build_package


def _row(provider: str, species: str = "Alpha one") -> dict[str, object]:
    status = (
        "accepted_wild_self_fertile"
        if provider == "ken_fern"
        else "accepted_self_fertile_statement"
    )
    return {
        "accepted_species": species,
        "source_url": f"https://{provider}.example/{species.replace(' ', '+')}",
        "supporting_excerpt": "Self-fertile: Yes; Pollinators: Bees; Cultivation Status: Wild",
        "normalized_value": "SC",
        "trait_name": "self_incompatibility",
        "status": status,
        "page_sha256": provider * 8,
        "content_fingerprint": "same-excerpt",
        "nonsexual_qualified": False,
        "cultivar_qualified": False,
        "dioecious_conflict": False,
        "retrieved_at_utc": "2026-08-27T00:00:00Z",
        "cache_file": f"{provider}.html",
    }


def test_provider_redistributions_share_one_lineage() -> None:
    evidence, audit = build_package(
        pd.DataFrame([_row("ken_fern")]), pd.DataFrame([_row("pfaf")])
    )
    assert len(evidence) == 2
    assert evidence["source_lineage"].nunique() == 1
    assert evidence["source_provider"].nunique() == 2
    assert set(audit["decision"]) == {"accept"}


def test_self_fertile_maps_to_sc_not_autonomous_selfing() -> None:
    empty_pfaf = pd.DataFrame(columns=_row("pfaf").keys())
    evidence, _ = build_package(pd.DataFrame([_row("ken_fern")]), empty_pfaf)
    assert set(evidence["trait_name"]) == {"self_incompatibility"}
    assert set(evidence["normalized_value"]) == {"SC"}


def test_nonsexual_row_is_fail_closed() -> None:
    row = _row("ken_fern")
    row["nonsexual_qualified"] = True
    empty_pfaf = pd.DataFrame(columns=_row("pfaf").keys())
    evidence, audit = build_package(pd.DataFrame([row]), empty_pfaf)
    assert evidence.empty
    assert audit.iloc[0]["reason"] == "nonsexual_reproduction_category"
