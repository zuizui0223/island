from __future__ import annotations

import pandas as pd
import pytest

from island_v2.efloras_medium_checkpoint import apply_reviewed_audit, high_precision_candidate


def candidate(trait: str, value: str, excerpt: str) -> dict[str, str]:
    return {
        "trait_name": trait,
        "candidate_value": value,
        "source_excerpt": excerpt,
    }


def test_size_gate_keeps_whole_flower_measurements_and_rejects_parts() -> None:
    assert high_precision_candidate(
        candidate("flower_size_class", "medium", "Flowers 16–18 mm diam.")
    )[0]
    assert not high_precision_candidate(
        candidate("flower_size_class", "small", "Flowers: stamens 3–4 mm.")
    )[0]
    assert not high_precision_candidate(
        candidate("flower_size_class", "small", "Corolla yellow; petals 5–10 mm.")
    )[0]


def test_form_gate_requires_the_flower_corolla_or_perianth_subject() -> None:
    assert high_precision_candidate(
        candidate("floral_form", "funnelform", "Corolla white, funnelform, 50–80 mm.")
    )[0]
    assert not high_precision_candidate(
        candidate("floral_form", "campanulate", "Flowers: calyx campanulate, 4 mm.")
    )[0]
    assert not high_precision_candidate(
        candidate("floral_form", "rotate", "Flowers with stellate-pubescent bracts.")
    )[0]


def test_tube_gate_rejects_calyx_width_fruit_and_subject_switches() -> None:
    assert high_precision_candidate(
        candidate("tube_depth_class", "intermediate", "Corolla tube 7–10 mm long.")
    )[0]
    assert not high_precision_candidate(
        candidate("tube_depth_class", "shallow", "Calyx tube 2–3 mm.")
    )[0]
    assert not high_precision_candidate(
        candidate("tube_depth_class", "shallow", "Fruiting perianth tube 2 mm.")
    )[0]


def test_display_gate_excludes_comparisons_and_cones() -> None:
    assert high_precision_candidate(
        candidate(
            "inflorescence_display",
            "raceme_spike_panicle",
            "Inflorescences terminal, paniculate cymes, 3–7 cm.",
        )
    )[0]
    assert not high_precision_candidate(
        candidate(
            "inflorescence_display",
            "solitary",
            "In contrast, Other species has solitary flowers.",
        )
    )[0]
    assert not high_precision_candidate(
        candidate(
            "inflorescence_display",
            "composite_display",
            "Pollen-cone capitula of 6–14 cones.",
        )
    )[0]


def test_cultivar_context_is_fail_closed() -> None:
    accepted, reason = high_precision_candidate(
        candidate(
            "flower_size_class",
            "large",
            "Flowers 8 cm diam. in cultivars.",
        )
    )
    assert not accepted
    assert reason == "cultivar_or_horticultural_context"


def test_size_gate_rejects_measurements_of_surrounding_organs() -> None:
    assert not high_precision_candidate(
        candidate(
            "flower_size_class",
            "medium",
            "Inflorescence a raceme with more than 10 flowers, 5-8 cm.",
        )
    )[0]
    assert not high_precision_candidate(
        candidate(
            "flower_size_class",
            "medium",
            "Filaments 3 cm, extending to corolla mouth, 1 shorter, 1-2 cm.",
        )
    )[0]


def test_form_gate_rejects_comparative_corolla_reference() -> None:
    assert not high_precision_candidate(
        candidate(
            "floral_form",
            "tubular",
            "Corona segments paler than corolla, stipitate, tubular.",
        )
    )[0]


def test_display_gate_rejects_incomplete_solitary_multistate() -> None:
    accepted, reason = high_precision_candidate(
        candidate(
            "inflorescence_display",
            "solitary",
            "Flowers axillary, solitary or 2-7-fascicled.",
        )
    )
    assert not accepted
    assert reason == "incomplete_multistate_context"


def test_tube_gate_rejects_nearby_scale_and_filament_measurements() -> None:
    assert not high_precision_candidate(
        candidate(
            "tube_depth_class",
            "absent_or_open",
            "Scales shorter than corolla tube length, bridged at 0.1-0.2 mm.",
        )
    )[0]
    assert not high_precision_candidate(
        candidate(
            "tube_depth_class",
            "deep",
            "Filaments inserted at mid perianth tube, erect, 4-5 cm.",
        )
    )[0]


def test_taxonomic_comparison_note_is_not_direct_evidence() -> None:
    accepted, reason = high_precision_candidate(
        candidate(
            "tube_depth_class",
            "deep",
            "Author treated P. sinica as a species; its corolla tube was 17 mm, "
            "which does not agree with P. arenosa, but is included in the "
            "circumscription of P. arenosa.",
        )
    )
    assert not accepted
    assert reason == "taxonomic_comparative_note"


def test_size_gate_rejects_in_flower_and_spike_measurements() -> None:
    assert not high_precision_candidate(
        candidate(
            "flower_size_class",
            "medium",
            "Pedicels 0.2-1.5 cm in flower, to 3 cm in fruit.",
        )
    )[0]
    assert not high_precision_candidate(
        candidate(
            "flower_size_class",
            "large",
            "Terminal spike with a few female flowers, linear, 4-8 cm.",
        )
    )[0]


def test_size_gate_uses_long_dimension_for_two_dimensional_measurement() -> None:
    assert high_precision_candidate(
        candidate(
            "flower_size_class",
            "large",
            "Corollas yellow, laminae 48-69 x 9-16 mm.",
        )
    )[0]
    assert not high_precision_candidate(
        candidate(
            "flower_size_class",
            "medium",
            "Corollas yellow, laminae 48-69 x 9-16 mm.",
        )
    )[0]


def test_display_gate_rejects_unmapped_solitary_alternatives_and_cyathium() -> None:
    for excerpt in (
        "Flowers axillary, solitary or paired.",
        "Inflorescences solitary flowers, rarely verticils.",
        "Inflorescence a terminal solitary cyathium.",
    ):
        assert not high_precision_candidate(
            candidate("inflorescence_display", "solitary", excerpt)
        )[0]


def test_size_gate_rejects_post_measurement_peduncle_and_tepal_lengths() -> None:
    for excerpt in (
        "Female flowers solitary or 2-5 on 6-8 mm peduncle.",
        "Flowers: tepals white, lanceolate, 3-5 mm, apex acuminate.",
    ):
        assert not high_precision_candidate(
            candidate("flower_size_class", "small", excerpt)
        )[0]


def test_size_gate_rejects_unmapped_range_and_other_species_note() -> None:
    assert not high_precision_candidate(
        candidate(
            "flower_size_class",
            "small",
            "Corolla 3-10 mm or much smaller when cleistogamous.",
        )
    )[0]
    assert not high_precision_candidate(
        candidate(
            "flower_size_class",
            "very_large",
            "Papaver bracteatum, which some authors have included in P. orientale, "
            "is similar but has flowers to 20 cm diam.",
        )
    )[0]


def test_size_gate_rejects_bud_and_spur_measurements() -> None:
    for excerpt in (
        "Female flower ellipsoid or ovoid in bud, 5-7 mm.",
        "Nodes with 2-5 flowers clustered on a 1-3 mm spur.",
    ):
        assert not high_precision_candidate(
            candidate("flower_size_class", "small", excerpt)
        )[0]


def test_reviewed_audit_can_only_add_decisions() -> None:
    row = {
        "candidate_id": "candidate-1",
        "trait_name": "floral_form",
        "accepted_species": "Example species",
        "normalized_value": "tubular",
        "source_url": "https://example.org/species",
        "source_excerpt": "Corolla tubular.",
        "accepted_correct": "",
        "cultivar_status": "",
        "reviewer": "",
        "reviewed_at_utc": "",
        "audit_reason": "",
        "audit_hash": "hash-1",
    }
    reviewed_row = row | {
        "accepted_correct": "true",
        "cultivar_status": "wild_species_description",
        "reviewer": "Reviewer",
        "reviewed_at_utc": "2026-08-10T00:00:00Z",
        "audit_reason": "quote supports value",
    }
    result = apply_reviewed_audit(pd.DataFrame([row]), pd.DataFrame([reviewed_row]))
    assert result.loc[0, "accepted_correct"] == "true"

    changed = reviewed_row | {"normalized_value": "bell_campanulate"}
    with pytest.raises(ValueError, match="immutable column"):
        apply_reviewed_audit(pd.DataFrame([row]), pd.DataFrame([changed]))
