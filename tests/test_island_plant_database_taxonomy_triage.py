from __future__ import annotations

import pandas as pd

from island_v2.island_plant_database_taxonomy_triage import build_triage


def _row(
    name: str,
    *,
    match_type: str = "NONE",
    matched_rank: str = "",
    matched_status: str = "",
    accepted_rank: str = "",
    accepted_status: str = "",
    kingdom: str = "",
    flags: str = "NO_MATCH",
    resolution_status: str = "unmatched",
) -> dict[str, str]:
    return {
        "submitted_name": name,
        "match_type": match_type,
        "matched_rank": matched_rank,
        "matched_status": matched_status,
        "candidate_accepted_rank": accepted_rank,
        "candidate_accepted_status": accepted_status,
        "candidate_kingdom": kingdom,
        "processing_flags": flags,
        "issues": "",
        "resolution_status": resolution_status,
    }


def test_triage_separates_source_scope_and_second_backbone_candidates() -> None:
    frame = pd.DataFrame(
        [
            _row("Alpha beta"),
            _row("Alpha"),
            _row("Alpha beta x gamma"),
            _row("Acacia Island-phrase"),
        ]
    )
    out, summary = build_triage(frame)
    assert out["source_name_rank_class"].tolist() == [
        "binomial_candidate",
        "genus_or_higher_rank",
        "infra_or_hybrid_candidate",
        "unresolved_name_rank",
    ]
    assert out["second_backbone_eligible"].tolist() == ["true", "false", "false", "false"]
    assert summary["n_second_backbone_eligible"] == 1
    assert summary["n_source_scope_review"] == 3


def test_non_plantae_match_is_not_second_backbone_species_candidate() -> None:
    frame = pd.DataFrame(
        [
            _row(
                "Gonatopus example",
                match_type="EXACT",
                matched_rank="SPECIES",
                matched_status="ACCEPTED",
                accepted_rank="SPECIES",
                accepted_status="ACCEPTED",
                kingdom="Animalia",
                flags="",
                resolution_status="review_required",
            )
        ]
    )
    out, summary = build_triage(frame)
    assert out.loc[0, "triage_class"] == "kingdom_conflict_non_plantae"
    assert out.loc[0, "source_scope_status"] == "kingdom_conflict_review"
    assert out.loc[0, "second_backbone_eligible"] == "false"
    assert summary["n_kingdom_conflict_review"] == 1


def test_colxr_failure_classes_remain_fail_closed() -> None:
    frame = pd.DataFrame(
        [
            _row("Alpha beta", flags="MULTIPLE_MATCHES_SAME_CONFIDENCE"),
            _row("Gamma delta", flags="LOW_CONFIDENCE"),
            _row(
                "Epsilon zeta",
                match_type="EXACT",
                matched_rank="SPECIES",
                matched_status="PROVISIONALLY_ACCEPTED",
                accepted_rank="SPECIES",
                accepted_status="PROVISIONALLY_ACCEPTED",
                kingdom="Plantae",
                flags="",
                resolution_status="review_required",
            ),
            _row(
                "Eta theta",
                match_type="EXACT",
                matched_rank="SPECIES",
                matched_status="SYNONYM",
                accepted_rank="SUBSPECIES",
                kingdom="Plantae",
                flags="",
                resolution_status="review_required",
            ),
            _row(
                "Iota kappa",
                match_type="VARIANT",
                matched_rank="SPECIES",
                matched_status="ACCEPTED",
                accepted_rank="SPECIES",
                accepted_status="ACCEPTED",
                kingdom="Plantae",
                flags="",
                resolution_status="review_required",
            ),
            _row(
                "Lambda mu",
                match_type="HIGHERRANK",
                matched_rank="GENUS",
                matched_status="ACCEPTED",
                accepted_rank="GENUS",
                accepted_status="ACCEPTED",
                kingdom="Plantae",
                flags="",
                resolution_status="review_required",
            ),
        ]
    )
    out, summary = build_triage(frame)
    assert out["triage_class"].tolist() == [
        "colxr_none_ambiguous",
        "colxr_none_low_confidence",
        "colxr_exact_provisionally_accepted",
        "colxr_exact_target_infraspecific",
        "colxr_variant_species_name",
        "colxr_higher_rank_only",
    ]
    assert out["second_backbone_eligible"].eq("true").all()
    assert summary["n_second_backbone_eligible"] == 6
    assert summary["n_source_scope_review"] == 0
    assert summary["n_kingdom_conflict_review"] == 0
