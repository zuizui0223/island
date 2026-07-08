from pathlib import Path

import pandas as pd
import pytest
import typer

from island_v2.trait_lead_packet import compress_review_packet, load_config

CONFIG = load_config(Path("config/lead_review_packet.yml"))


def _leads(n_species, per_cell, groups=("A", "B", "C"), taxon_score=3, start=0):
    rows = []
    for s in range(start, start + n_species):
        for g in groups:
            for k in range(per_cell):
                rows.append(
                    {
                        "query_taxon": f"Genus sp{s:03d}",
                        "query_trait_group": g,
                        "title": f"paper {s}-{g}-{k}",
                        "taxon_relevance_score": str(taxon_score),
                        "title_relevance_score": str(per_cell - k),  # first is strongest
                        "abstract_relevance_score": "0",
                        "priority_rank": "1",
                    }
                )
    return pd.DataFrame(rows)


def test_tiers_are_split_S_packet_A_rescue_B_background():
    s = _leads(2, 2, taxon_score=4, start=0)   # Tier S
    a = _leads(2, 2, taxon_score=1, start=10)  # Tier A (genus only)
    b = _leads(2, 3, taxon_score=0, start=20)  # Tier B (background)
    leads = pd.concat([s, a, b], ignore_index=True)
    packet, rescue, background, summary = compress_review_packet(leads, CONFIG)
    assert summary["n_tier_S_species_confident"] == len(s)
    assert summary["n_tier_A_genus_flora"] == len(a)
    assert summary["n_tier_B_background"] == len(b)
    # packet holds only Tier S; rescue only Tier A; background only Tier B
    assert (packet["taxon_relevance_score"] >= 3).all()
    assert rescue["taxon_relevance_score"].isin([1, 2]).all()
    assert (background["taxon_relevance_score"] == 0).all()


def test_tier_is_derived_when_column_absent():
    # A lead table without taxon_relevance_tier (older discovery output).
    leads = _leads(1, 1, taxon_score=4)
    assert "taxon_relevance_tier" not in leads.columns
    packet, rescue, background, summary = compress_review_packet(leads, CONFIG)
    assert summary["n_tier_S_species_confident"] == 3  # 1 species x 3 groups, all score 4


def test_small_tier_s_is_kept_whole():
    packet, rescue, background, summary = compress_review_packet(_leads(2, 2, taxon_score=3), CONFIG)
    assert summary["n_selected_packet"] == 12  # 2*3*2, below target_max
    assert len(rescue) == 0 and len(background) == 0


def test_large_tier_s_is_compressed_into_band_and_spread():
    leads = _leads(50, 6, taxon_score=4)  # 900 Tier-S leads
    packet, rescue, background, summary = compress_review_packet(leads, CONFIG)
    assert summary["n_tier_S_species_confident"] == 900
    assert CONFIG["target_min"] <= summary["n_selected_packet"] <= CONFIG["target_max"]
    assert summary["n_species_in_packet"] == 50
    assert packet["query_taxon"].value_counts().max() <= (CONFIG["target_max"] // 50) + 1


def test_packet_keeps_strongest_lead_per_cell_first():
    leads = _leads(1, 5, groups=("A",), taxon_score=3)
    tiny = {**CONFIG, "target_max": 2, "target_min": 1}
    packet, *_ , summary = compress_review_packet(leads, tiny)
    assert summary["n_selected_packet"] == 2
    assert packet["title_relevance_score"].max() == 5


def test_prohibited_column_is_rejected():
    leads = _leads(1, 1).assign(native_status="unknown")
    with pytest.raises(typer.BadParameter, match="prohibited"):
        compress_review_packet(leads, CONFIG)


def test_missing_required_column_is_rejected():
    leads = _leads(1, 1).drop(columns=["taxon_relevance_score"])
    with pytest.raises(typer.BadParameter, match="missing columns"):
        compress_review_packet(leads, CONFIG)
