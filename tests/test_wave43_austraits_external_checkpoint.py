from __future__ import annotations

import pandas as pd
import pytest

from island_v2.wave43_austraits_external_checkpoint import validate_mapping


def _row(
    source: str,
    *,
    accepted: str = "Alpha outside",
    reason: str = "accepted_strict_two_backbone_external",
    source_family: str = "Testaceae",
    wfo_family: str = "Testaceae",
    gbif_family: str = "Testaceae",
) -> dict[str, str]:
    return {
        "source_name": source,
        "source_family": source_family,
        "accepted_species": accepted if reason.startswith("accepted_") else "",
        "mapping_reason": reason,
        "wfo_accepted_species": accepted,
        "wfo_family": wfo_family,
        "wfo_release": "2026-06",
        "gbif_accepted_species": accepted,
        "gbif_family": gbif_family,
        "gbif_match_type": "EXACT",
        "gbif_rank": "SPECIES",
        "gbif_kingdom": "Plantae",
    }


def test_wave43_mapping_accepts_only_external_exact_two_backbone_species() -> None:
    validate_mapping(
        pd.DataFrame([_row("Alpha old")]),
        target_species={"Alpha target"},
        expected_names=1,
    )


@pytest.mark.parametrize(
    ("row", "target"),
    [
        (_row("Alpha old", accepted="Alpha target"), {"Alpha target"}),
        (_row("Alpha old", source_family="Otheraceae"), set()),
        (_row("Alpha old", gbif_family="Otheraceae"), set()),
        (_row("Alpha old", accepted="Alpha outside subsp"), set()),
    ],
)
def test_wave43_mapping_fails_closed_on_identity_gate_violation(
    row: dict[str, str], target: set[str]
) -> None:
    with pytest.raises(ValueError, match="strict identity gate"):
        validate_mapping(pd.DataFrame([row]), target_species=target, expected_names=1)


def test_wave43_rejected_mapping_cannot_retain_accepted_species() -> None:
    row = _row("Alpha old", reason="backbones_disagree")
    row["accepted_species"] = "Alpha outside"

    with pytest.raises(ValueError, match="rejected mappings"):
        validate_mapping(pd.DataFrame([row]), target_species=set(), expected_names=1)
