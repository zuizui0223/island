from __future__ import annotations

import pandas as pd

from island_v2.useful_plants_reproduction_remine import build_package


def _candidate(provider: str, species: str = "Alpha beta") -> dict[str, str]:
    return {
        "accepted_species": species,
        "trait_name": "self_incompatibility",
        "normalized_value": "SI",
        "provider": provider,
        "source_url": f"https://{provider}.example/{species}",
        "supporting_excerpt": "Plants are self-sterile and require cross-fertilization [280].",
        "page_sha256": provider * 16,
        "excerpt_sha256": provider * 12,
        "content_fingerprint": "same-claim",
        "cache_file": f"{provider}.html",
    }


def _review(**overrides: str) -> dict[str, str]:
    row = {
        "accepted_species": "Alpha beta",
        "trait_name": "self_incompatibility",
        "normalized_value": "SI",
        "decision": "accept",
        "reason": "explicit self-sterility",
        "underlying_citation_ref": "280",
        "underlying_citation": "Stearn 2002.",
        "reviewer": "reviewer",
        "reviewed_at_utc": "2026-08-27T15:00:00Z",
    }
    row.update(overrides)
    return row


def test_provider_mirrors_and_species_share_underlying_citation_lineage() -> None:
    candidates = [
        _candidate("ken_fern"),
        _candidate("pfaf"),
        _candidate("ken_fern", "Alpha gamma"),
    ]
    review = pd.DataFrame([_review(), _review(accepted_species="Alpha gamma")])
    evidence, _ = build_package(candidates, review)
    assert len(evidence) == 3
    assert evidence.source_lineage.nunique() == 1
    assert set(evidence.trait_name) == {"self_incompatibility"}


def test_apomixis_review_remains_excluded() -> None:
    candidate = _candidate("ken_fern")
    review = pd.DataFrame(
        [
            _review(
                decision="exclude",
                reason="apomixis does not establish sexual self-compatibility",
            )
        ]
    )
    evidence, audit = build_package([candidate], review)
    assert evidence.empty
    assert audit.iloc[0].package_decision == "exclude"


def test_self_fertility_never_maps_to_autonomous_selfing() -> None:
    evidence, _ = build_package([_candidate("ken_fern")], pd.DataFrame([_review()]))
    assert set(evidence.trait_name) == {"self_incompatibility"}
    assert "autonomous_selfing_capacity" not in set(evidence.trait_name)
