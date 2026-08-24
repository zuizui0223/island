from __future__ import annotations

from pathlib import Path

import pandas as pd

from island_v2.review_inventory_artifact import build_inventory_audit


def _candidate(candidate_id: str, trait: str, value: str, raw: str) -> dict[str, object]:
    return {
        "candidate_id": candidate_id,
        "accepted_species": f"Alpha {candidate_id}",
        "trait_name": trait,
        "normalized_value": value,
        "raw_value": raw,
        "source_url": f"https://flora.example/plants/alpha-{candidate_id}",
        "page_title": f"Alpha {candidate_id}",
        "domain": "flora.example",
        "language": "en",
        "supporting_excerpt": f"Field: {raw}",
        "source_lineage": f"url:https://flora.example/plants/alpha-{candidate_id}",
        "name_match_method": "exact_accepted_name",
        "matched_page_name": f"Alpha {candidate_id}",
        "source_tier": "B",
        "source_type": "university_extension",
        "domain_review_status": "approved_registry_trait",
        "content_sha256": candidate_id * 64,
    }


def _prior() -> pd.DataFrame:
    from island_v2.review_inventory_artifact import AUDIT_COLUMNS

    return pd.DataFrame(columns=AUDIT_COLUMNS)


def test_fail_closed_filter_keeps_only_ontology_valid_unambiguous_rows() -> None:
    candidates = pd.DataFrame(
        [
            _candidate("a", "flower_primary_color", "white", "White"),
            _candidate("b", "flower_size_class", "medium", "< 1 inch"),
            _candidate(
                "c",
                "inflorescence_display",
                "other_described",
                "Panicle; Umbel",
            ),
            _candidate("d", "floral_form", "other_described", "Bell; Tubular"),
        ]
    )
    audit, exclusions, report = build_inventory_audit(
        candidates,
        _prior(),
        ontology_yaml=Path("config/trait_ontology.yml"),
        reviewer="test",
        reviewed_at_utc="2026-08-08T00:00:00Z",
        source_run_id="1",
        source_artifact="artifact-1",
    )
    assert audit["candidate_id"].tolist() == ["a"]
    assert set(exclusions["candidate_id"]) == {"b", "c", "d"}
    reasons = exclusions.set_index("candidate_id")["pre_audit_exclusion_reasons"]
    assert "ambiguous_upper_bound_measurement" in reasons["b"]
    assert "ontology_value_not_allowed" in reasons["c"]
    assert "multistate_collapsed_to_other_described" in reasons["d"]
    assert report["new_audit_precision"] == 1.0
    assert report["network_refetches"] == 0
