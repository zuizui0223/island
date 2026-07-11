from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import typer

from island_v2.global_candidate_index_audited import (
    apply_candidate_exclusions,
    load_candidate_exclusions,
)


def _raw_candidates() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "accepted_species": "Glossopteris indica",
                "trait_name": "sex_system",
                "candidate_value": "dioecious",
                "source_url": "https://ru.wikipedia.org/wiki/Glossopteris",
                "source_excerpt": "Неизвестно, были ли растения однодомными или двудомными.",
                "wave_id": "wave_00038_test",
            },
            {
                "accepted_species": "Glossopteris indica",
                "trait_name": "sex_system",
                "candidate_value": "monoecious",
                "source_url": "https://ru.wikipedia.org/wiki/Glossopteris",
                "source_excerpt": "Неизвестно, были ли растения однодомными или двудомными.",
                "wave_id": "wave_00038_test",
            },
            {
                "accepted_species": "Aucuba japonica",
                "trait_name": "sex_system",
                "candidate_value": "dioecious",
                "source_url": "https://ja.wikipedia.org/wiki/アオキ_(植物)",
                "source_excerpt": "雌雄異株である。",
                "wave_id": "wave_00038_test",
            },
        ]
    )


def test_apply_candidate_exclusions_retains_raw_audit() -> None:
    exclusions = pd.DataFrame(
        [
            {
                "accepted_species": "Glossopteris indica",
                "trait_name": "sex_system",
                "candidate_value": "dioecious",
                "source_url": "https://ru.wikipedia.org/wiki/Glossopteris",
                "exclusion_reason": "genus scope and uncertainty",
                "exclusion_rule": "scope_v1",
                "added_at": "2026-07-11",
            },
            {
                "accepted_species": "Glossopteris indica",
                "trait_name": "sex_system",
                "candidate_value": "monoecious",
                "source_url": "https://ru.wikipedia.org/wiki/Glossopteris",
                "exclusion_reason": "genus scope and uncertainty",
                "exclusion_rule": "scope_v1",
                "added_at": "2026-07-11",
            },
        ]
    )

    retained, matches = apply_candidate_exclusions(_raw_candidates(), exclusions)

    assert retained["accepted_species"].tolist() == ["Aucuba japonica"]
    assert len(matches) == 2
    assert set(matches["candidate_value"]) == {"dioecious", "monoecious"}
    assert set(matches["exclusion_rule"]) == {"scope_v1"}
    assert "source_excerpt" in matches.columns


def test_load_candidate_exclusions_rejects_duplicate_exact_keys(tmp_path: Path) -> None:
    path = tmp_path / "machine_candidate_exclusions.csv"
    row = {
        "accepted_species": "Glossopteris indica",
        "trait_name": "sex_system",
        "candidate_value": "dioecious",
        "source_url": "https://ru.wikipedia.org/wiki/Glossopteris",
        "exclusion_reason": "duplicate",
        "exclusion_rule": "scope_v1",
        "added_at": "2026-07-11",
    }
    pd.DataFrame([row, row]).to_csv(path, index=False)

    with pytest.raises(typer.BadParameter, match="duplicate exact keys"):
        load_candidate_exclusions(tmp_path)
