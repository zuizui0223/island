from pathlib import Path

import pandas as pd
import pytest
import typer

from island_v2.trait_lead_packet import compress_review_packet, load_config

CONFIG = load_config(Path("config/lead_review_packet.yml"))


def _leads(n_species: int, per_cell: int, groups=("A", "B", "C"), taxon_score=3) -> pd.DataFrame:
    rows = []
    for s in range(n_species):
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


def test_taxon_gate_drops_zero_taxon_relevance_leads():
    leads = _leads(2, 2)
    # add off-species leads that matched a trait but not the taxon
    off = _leads(2, 3, taxon_score=0)
    combined = pd.concat([leads, off], ignore_index=True)
    packet, summary = compress_review_packet(combined, CONFIG)
    assert (packet["taxon_relevance_score"] > 0).all()
    assert summary["n_gated_out_zero_taxon"] == len(off)
    assert summary["n_after_taxon_gate"] == len(leads)


def test_small_gated_set_is_kept_whole():
    leads = _leads(2, 2)  # 2*3*2 = 12 leads, below target_max
    packet, summary = compress_review_packet(leads, CONFIG)
    assert summary["n_selected"] == 12
    assert summary["n_after_taxon_gate"] == 12


def test_large_set_is_compressed_into_target_band_and_spread():
    # 50 species x 3 groups x 6 = 900 gated leads -> must compress to <= target_max
    leads = _leads(50, 6)
    packet, summary = compress_review_packet(leads, CONFIG)
    assert summary["n_after_taxon_gate"] == 900
    assert summary["n_selected"] <= CONFIG["target_max"]
    assert summary["n_selected"] >= CONFIG["target_min"]
    # every species represented (round-robin spread), none monopolises
    assert summary["n_species_in_packet"] == 50
    top_species_count = packet["query_taxon"].value_counts().max()
    assert top_species_count <= (CONFIG["target_max"] // 50) + 1


def test_compression_keeps_strongest_lead_per_cell_first():
    # one species, one group, 5 leads of descending title relevance; cap forces top picks
    leads = _leads(1, 5, groups=("A",))
    tiny_cfg = {**CONFIG, "target_max": 2, "target_min": 1}
    packet, summary = compress_review_packet(leads, tiny_cfg)
    assert summary["n_selected"] == 2
    # strongest title_relevance_score (5) must be kept
    assert packet["title_relevance_score"].max() == 5


def test_prohibited_column_is_rejected():
    leads = _leads(1, 1).assign(native_status="unknown")
    with pytest.raises(typer.BadParameter, match="prohibited"):
        compress_review_packet(leads, CONFIG)


def test_missing_required_column_is_rejected():
    leads = _leads(1, 1).drop(columns=["taxon_relevance_score"])
    with pytest.raises(typer.BadParameter, match="missing columns"):
        compress_review_packet(leads, CONFIG)
