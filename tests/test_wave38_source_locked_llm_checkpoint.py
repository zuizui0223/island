from __future__ import annotations

from island_v2.wave38_source_locked_llm_checkpoint import candidate_supported


def _row(trait: str, value: str, quote: str) -> dict[str, str]:
    return {
        "trait_name": trait,
        "provisional_candidate_value": value,
        "evidence_quote": quote,
        "source_url": "https://example.test/species",
    }


def test_rejects_sex_system_as_mating_system() -> None:
    accepted, reason, _ = candidate_supported(
        _row(
            "mating_system",
            "predominantly_outcrossing",
            "A dioecious species, both male and female forms need to be grown.",
        )
    )
    assert not accepted
    assert reason == "sex_system_or_pollination_not_mating_system"


def test_rejects_self_fertile_no_as_si() -> None:
    accepted, reason, _ = candidate_supported(
        _row("self_incompatibility", "SI", "Self-fertile. No.")
    )
    assert not accepted
    assert reason == "self_fertility_category_error"


def test_accepts_explicit_sc_and_autonomous_selfing_separately() -> None:
    assert candidate_supported(
        _row("self_incompatibility", "SC", "Self-fertile: Yes")
    )[0]
    assert candidate_supported(
        _row(
            "autonomous_selfing_capacity",
            "autonomous",
            "The flowers are self-pollinating.",
        )
    )[0]


def test_rejects_calyx_shape_as_whole_flower_form() -> None:
    accepted, reason, _ = candidate_supported(
        _row("floral_form", "bell_campanulate", "Calyx campanulate")
    )
    assert not accepted
    assert reason == "form_not_explicit"


def test_accepts_explicit_flower_form_and_symmetry() -> None:
    assert candidate_supported(
        _row("floral_form", "tubular", "The plant bears dark red tubular flowers.")
    )[0]
    assert candidate_supported(
        _row(
            "floral_symmetry",
            "zygomorphic",
            "Each flower is bilaterally symmetrical.",
        )
    )[0]


def test_preserves_explicit_multicolour_states() -> None:
    accepted, _, normalized = candidate_supported(
        _row(
            "flower_primary_color",
            "multicolored_variable",
            "The corolla is white or violet.",
        )
    )
    assert accepted
    assert normalized == "blue_purple|white|multicolored_variable"


def test_ignores_fruit_colour_after_flower_clause() -> None:
    accepted, _, _ = candidate_supported(
        _row(
            "flower_primary_color",
            "red_pink",
            "The flowers are white; fruits turn red when ripe.",
        )
    )
    assert not accepted


def test_numeric_flower_and_tube_thresholds_match_shared_contract() -> None:
    assert candidate_supported(
        _row("flower_size_class", "small", "Flowers are 8-14 mm wide.")
    )[0]
    assert candidate_supported(
        _row("tube_depth_class", "intermediate", "The corolla tube is 6-8 mm long.")
    )[0]


def test_rejects_other_described_and_cultivar_evidence() -> None:
    assert not candidate_supported(
        _row("floral_form", "other_described", "Flowers with spreading petals.")
    )[0]
    assert not candidate_supported(
        _row("flower_primary_color", "red_pink", "Cultivar flowers are red.")
    )[0]
