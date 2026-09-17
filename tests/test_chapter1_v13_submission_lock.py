from __future__ import annotations

import json
from pathlib import Path

import pytest

from island_v2.chapter1_v13_submission_lock import validate_v13_lock

LOCK = Path("config/chapter1_v13_unified_island_syndrome_result_lock.json")


def _valid_lock() -> dict:
    return {
        "contract": "chapter1_v13_unified_island_syndrome_result_lock_v2",
        "parents": {
            "all_data": {"contract": "chapter1_all_data_route_result_lock_v2"},
            "v12_two_panel": {"contract": "chapter1_v12_two_panel_result_lock_v7"},
            "v12_glopl": {"contract": "chapter1_v12_h5_glopl_extension_result_lock_v4"},
            "functional_bridge": {
                "contract": "chapter1_v13_functional_bridge_result_lock_v1",
                "artifact_id": 10465048981,
                "artifact_digest": "sha256:50b6ca693ec989690854a9e9aadb5cd8a672b5fbcf3ff84e5aabbab06a6abc8e",
            },
        },
        "architecture": {
            "H1": "global_recurrent_island_syndrome",
            "H2": "global_pollination_constraint",
            "H3": "dual_plant_response_pathways",
            "H4": "posthoc_functional_triangulation",
        },
        "H1_global_recurrent_island_syndrome": {
            "between_region_contrasts_used_in_submission": False,
        },
        "H4_functional_bridge": {
            "inferential_role": "posthoc_functional_triangulation",
            "artifact_id": 10465048981,
            "historical_selection_identified": False,
        },
        "supplementary": {"GloBI_role": "non_promoted_sampling_sensitive_context_evidence"},
        "claim_ceiling": {
            "historical_pollen_limitation_selected_traits": False,
            "pollen_limitation_mediates_global_syndrome": False,
            "GloBI_identifies_pollinator_mechanism": False,
            "equal_response_vectors_across_regions": False,
            "between_region_difference_is_primary_claim": False,
        },
    }


def test_valid_lock_passes() -> None:
    report = validate_v13_lock(_valid_lock())
    assert report["verified"] is True


def test_h5_reintroduction_fails() -> None:
    lock = _valid_lock()
    lock["architecture"]["H5"] = "taxonomic_realization"
    with pytest.raises(ValueError, match="H1-H4"):
        validate_v13_lock(lock)


def test_between_region_promotion_fails() -> None:
    lock = _valid_lock()
    lock["H1_global_recurrent_island_syndrome"]["between_region_contrasts_used_in_submission"] = True
    with pytest.raises(ValueError, match="between-region"):
        validate_v13_lock(lock)


def test_confirmatory_relabel_fails() -> None:
    lock = _valid_lock()
    lock["H4_functional_bridge"]["inferential_role"] = "confirmatory"
    with pytest.raises(ValueError, match="posthoc"):
        validate_v13_lock(lock)


def test_historical_selection_claim_fails() -> None:
    lock = _valid_lock()
    lock["claim_ceiling"]["historical_pollen_limitation_selected_traits"] = True
    with pytest.raises(ValueError, match="historical"):
        validate_v13_lock(lock)


def test_globi_mechanism_promotion_fails() -> None:
    lock = _valid_lock()
    lock["claim_ceiling"]["GloBI_identifies_pollinator_mechanism"] = True
    with pytest.raises(ValueError, match="GloBI"):
        validate_v13_lock(lock)


def test_missing_functional_bridge_provenance_fails() -> None:
    lock = _valid_lock()
    del lock["parents"]["functional_bridge"]["artifact_digest"]
    with pytest.raises(ValueError, match="functional bridge provenance"):
        validate_v13_lock(lock)


def test_repository_lock_passes_validator() -> None:
    if not LOCK.exists():
        pytest.skip("v13 result lock not created yet")
    report = validate_v13_lock(json.loads(LOCK.read_text(encoding="utf-8")))
    assert report["verified"] is True
