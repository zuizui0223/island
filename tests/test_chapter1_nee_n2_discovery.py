from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from island_v2.chapter1_nee_n2_discovery import (
    build_ranked_universe,
    load_config,
    select_wave,
    source_genus_token,
)

CONFIG = Path("config/chapter1_nee_n2_dependency_discovery.yml")


def _config() -> dict[str, object]:
    return load_config(CONFIG)


def _assignments() -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for island, entities in {
        "i1": ["e1", "e2", "e3", "e4", "e5"],
        "i2": ["e1", "e2", "e6", "e7", "e8"],
        "i3": ["e2", "e3", "e6", "e9", "e10"],
    }.items():
        for rank, entity in enumerate(entities, start=1):
            rows.append(
                {
                    "island_id": island,
                    "source_mode": "geo_k5",
                    "source_rank": str(rank),
                    "entity_ID": entity,
                }
            )
    rows.append(
        {
            "island_id": "i1",
            "source_mode": "geo_k10",
            "source_rank": "1",
            "entity_ID": "ignored",
        }
    )
    return pd.DataFrame(rows)


def _flora() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"entity_ID": "e1", "work_species": "Alpha one"},
            {"entity_ID": "e2", "work_species": "Alpha two"},
            {"entity_ID": "e2", "work_species": "Beta one"},
            {"entity_ID": "e3", "work_species": "Beta two"},
            {"entity_ID": "e4", "work_species": "Gamma one"},
            {"entity_ID": "e5", "work_species": "Delta one"},
            {"entity_ID": "e6", "work_species": "Beta three"},
            {"entity_ID": "e7", "work_species": "Epsilon one"},
            {"entity_ID": "e8", "work_species": "Zeta one"},
            {"entity_ID": "e9", "work_species": "Eta one"},
            {"entity_ID": "e10", "work_species": "Theta one"},
            {"entity_ID": "e1", "work_species": "× Hybridus one"},
            {"entity_ID": "e3", "work_species": "×Cephalorhiza hybrida"},
            {"entity_ID": "ignored", "work_species": "Outside one"},
        ]
    )


def test_literal_genus_token_holds_out_hybrid_markers() -> None:
    assert source_genus_token("Campanula punctata") == ("Campanula", None)
    assert source_genus_token("× Hybridus one") == ("×", "hybrid_or_nothogenus_marker")
    assert source_genus_token("×Cephalorhiza hybrida") == (
        "×Cephalorhiza",
        "hybrid_or_nothogenus_marker",
    )


def test_ranked_universe_uses_source_opportunity_only_and_is_deterministic() -> None:
    ranked, holdouts, receipt = build_ranked_universe(_assignments(), _flora(), _config())
    # Alpha and Beta are available to the same three islands; Beta wins only because
    # it occurs in three selected mainland source entities versus two for Alpha.
    assert ranked.iloc[0]["accepted_genus"] == "Beta"
    assert ranked.iloc[0]["n_islands_with_source_available_genus"] == 3
    assert ranked.iloc[0]["n_selected_mainland_source_entities_containing_genus"] == 3
    assert ranked.iloc[1]["accepted_genus"] == "Alpha"
    assert ranked.iloc[1]["n_islands_with_source_available_genus"] == 3
    assert ranked.iloc[1]["n_selected_mainland_source_entities_containing_genus"] == 2
    assert set(holdouts["raw_genus_token"]) == {"×", "×Cephalorhiza"}
    assert receipt["n_islands"] == 3
    assert receipt["uses_island_genus_entry"] is False
    assert receipt["uses_N2_effect_direction"] is False
    assert receipt["dependency_evidence_opened"] is False


def test_tie_break_is_lexicographic_after_both_source_scores() -> None:
    assignments = _assignments()
    flora = _flora()
    extra = pd.DataFrame(
        [
            {"entity_ID": "e4", "work_species": "Aardvark one"},
            {"entity_ID": "e4", "work_species": "Beacon one"},
        ]
    )
    ranked, _, _ = build_ranked_universe(assignments, pd.concat([flora, extra], ignore_index=True), _config())
    pair = ranked.loc[ranked["accepted_genus"].isin(["Aardvark", "Beacon"])].reset_index(drop=True)
    assert pair["n_islands_with_source_available_genus"].nunique() == 1
    assert pair["n_selected_mainland_source_entities_containing_genus"].nunique() == 1
    assert pair["accepted_genus"].tolist() == ["Aardvark", "Beacon"]


def test_missing_or_duplicate_geo_k5_rank_is_hard_stop() -> None:
    missing = _assignments().loc[lambda t: ~((t.island_id == "i1") & (t.source_mode == "geo_k5") & (t.source_rank == "5"))]
    with pytest.raises(ValueError, match="exactly ranks"):
        build_ranked_universe(missing, _flora(), _config())

    duplicated = pd.concat([_assignments(), _assignments().iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate island x source_rank"):
        build_ranked_universe(duplicated, _flora(), _config())


def test_wave_partition_is_fixed_at_500_and_select_wave_is_pure() -> None:
    rows = []
    assignments = _assignments()
    selected_entities = sorted(set(assignments.loc[assignments.source_mode.eq("geo_k5"), "entity_ID"]))
    for index in range(1001):
        letters = "".join(chr(65 + ((index // (26 ** p)) % 26)) for p in [2, 1, 0])
        genus = "Genus" + letters
        rows.append({"entity_ID": selected_entities[index % len(selected_entities)], "work_species": f"{genus} species"})
    ranked, _, receipt = build_ranked_universe(assignments, pd.DataFrame(rows), _config())
    assert receipt["n_waves"] == 3
    assert len(select_wave(ranked, 1)) == 500
    assert len(select_wave(ranked, 2)) == 500
    assert len(select_wave(ranked, 3)) == 1
    assert select_wave(ranked, 2)["global_rank"].min() == 501
    assert select_wave(ranked, 2)["global_rank"].max() == 1000
