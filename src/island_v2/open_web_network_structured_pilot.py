"""CLI wrapper adding explicit key-value fields to broad Web extraction."""

from __future__ import annotations

from island_v2 import open_web_network_pilot as network
from island_v2.open_web_structured_traits import extract_structured_characteristics

_RULE_EXTRACTOR = network.extract_rules


def _combined_extract_rules(
    page,
    *,
    trait_name: str,
    config,
    treatment_end_pattern: str = "",
):
    sentence_rows = _RULE_EXTRACTOR(
        page,
        trait_name=trait_name,
        config=config,
        treatment_end_pattern=treatment_end_pattern,
    )
    structured_rows = extract_structured_characteristics(page, trait_name=trait_name)
    return list(dict.fromkeys([*sentence_rows, *structured_rows]))


network.extract_rules = _combined_extract_rules


if __name__ == "__main__":
    network.app()
