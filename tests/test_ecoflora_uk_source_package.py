from __future__ import annotations

import importlib.util
from pathlib import Path

MODULE_PATH = (
    Path(__file__).parents[1]
    / "analysis"
    / "prepare_ecoflora_uk_source_package_20260808.py"
)
SPEC = importlib.util.spec_from_file_location("ecoflora_source_package", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_flower_colour_preserves_multiple_states() -> None:
    assert MODULE._states("white, yellow", MODULE.COLOUR_PATTERNS) == "white|yellow_orange"
    assert (
        MODULE._states("dark purple brown", MODULE.COLOUR_PATTERNS)
        == "blue_purple|green_brown_inconspicuous"
    )


def test_reproductive_traits_are_not_used_as_cross_trait_proxies() -> None:
    assert MODULE._mating_system("cross or automatic self") == "mixed_mating"
    assert MODULE._self_incompatibility("none") == "SC"
    assert MODULE._self_incompatibility("present, mechanism unknown") == "SI"
    assert MODULE._cleistogamy("some cleistogamous flowers") == "facultative"
    assert MODULE._cleistogamy("cleistogamy not recorded") == ""


def test_ambiguous_or_non_mating_reproductive_values_fail_closed() -> None:
    assert MODULE._mating_system("normally cross, normally self") == ""
    assert MODULE._mating_system("apomictic, normally cross") == ""
    assert MODULE._cleistogamy("pseudo-cleistogamous") == ""
