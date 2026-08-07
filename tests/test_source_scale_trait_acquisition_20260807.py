from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.fetch_europe_pmc_reproduction_scale_20260807 import VALUE_PATTERNS
from analysis.integrate_europe_pmc_reproduction_scale_20260807 import (
    INVALIDATED_PRIOR,
    PROMOTIONS,
    common_evidence,
    review_candidates,
)

ROOT = Path(__file__).resolve().parents[1]
CANDIDATES = (
    ROOT
    / "data/v2/staging/traits/direct_llm_pilot"
    / "20260807_europe_pmc_reproduction_scale"
    / "fulltext_reproductive_candidates.csv.gz"
)


def test_every_fulltext_candidate_has_an_explicit_review_decision() -> None:
    candidates = pd.read_csv(CANDIDATES, dtype=str).fillna("")
    review = review_candidates(candidates)

    assert len(review) == 43
    assert review["candidate_id"].nunique() == 43
    assert int(review["accepted_correct"].sum()) == 17
    assert int(review["promotion_accepted"].sum()) == 8
    assert review["review_decision"].value_counts().to_dict() == {
        "rejected": 26,
        "corroborating_not_counted": 9,
        "accepted": 8,
    }


def test_promoted_rows_are_unique_direct_species_traits() -> None:
    candidates = pd.read_csv(CANDIDATES, dtype=str).fillna("")
    common = common_evidence(review_candidates(candidates))

    assert len(common) == len(PROMOTIONS) == 8
    assert not common.duplicated(["accepted_species", "trait_name"]).any()
    assert set(common["quality"]) == {"high", "medium"}
    assert set(common["evidence_scope"]) == {"species_direct"}
    assert common["source_excerpt"].str.len().gt(0).all()
    assert common["source_lineage"].str.startswith(("doi:", "pmcid:")).all()


def test_known_panicum_comparison_misattribution_is_invalidated() -> None:
    assert (
        "Panicum hallii",
        "self_incompatibility",
        "SI",
        "url:https://europepmc.org/article/MED/28449656",
    ) in INVALIDATED_PRIOR


def test_self_fertile_is_not_an_autonomous_selfing_pattern() -> None:
    autonomous_patterns = VALUE_PATTERNS["autonomous_selfing_capacity"]
    assert not any(
        pattern.search("The plant is self-fertile.") for _, pattern in autonomous_patterns
    )
