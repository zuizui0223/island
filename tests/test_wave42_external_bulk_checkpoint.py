from __future__ import annotations

import pandas as pd
import pytest

from island_v2.wave42_external_bulk_checkpoint import validate_mapping


def _row(
    source: str,
    *,
    accepted: str = "Alpha outside",
    reason: str = "accepted_strict_two_backbone_external",
    wfo_family: str = "Testaceae",
    gbif_family: str = "Testaceae",
) -> dict[str, str]:
    return {
        "source_name": source,
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


def test_wave42_mapping_accepts_only_external_exact_two_backbone_species() -> None:
    mapping = pd.DataFrame([_row("Alpha old")])

    validate_mapping(mapping, target_species={"Alpha target"}, expected_names=1)


@pytest.mark.parametrize(
    ("row", "target"),
    [
        (_row("Alpha old", accepted="Alpha target"), {"Alpha target"}),
        (_row("Alpha old", gbif_family="Otheraceae"), set()),
        (_row("Alpha old", accepted="Alpha outside subsp"), set()),
    ],
)
def test_wave42_mapping_fails_closed_on_identity_gate_violation(
    row: dict[str, str], target: set[str]
) -> None:
    with pytest.raises(ValueError, match="strict identity gate"):
        validate_mapping(pd.DataFrame([row]), target_species=target, expected_names=1)


def test_wave42_rejected_mapping_cannot_retain_accepted_species() -> None:
    row = _row("Alpha old", reason="backbones_disagree")
    row["accepted_species"] = "Alpha outside"

    with pytest.raises(ValueError, match="rejected mappings"):
        validate_mapping(pd.DataFrame([row]), target_species=set(), expected_names=1)
