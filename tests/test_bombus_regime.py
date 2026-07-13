from island_v2.bombus_regime import classify_regime


def test_structural_absence_is_not_treated_as_arrival_failure() -> None:
    assert classify_regime("structurally_not_applicable", 0.9) == "structural_bombus_absence_comparison"


def test_applicable_region_uses_environment_only_after_applicability_freeze() -> None:
    assert classify_regime("applicable", 0.8) == "applicable_environment_compatible"
    assert classify_regime("applicable", 0.2) == "applicable_environment_mismatch"


def test_unresolved_applicability_stays_unresolved() -> None:
    assert classify_regime("unresolved", 0.9) == "applicability_unresolved"
