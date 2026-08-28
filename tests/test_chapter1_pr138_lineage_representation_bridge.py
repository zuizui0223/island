import numpy as np
import pandas as pd

from island_v2.chapter1_pr138_lineage_representation_bridge import (
    compute_island_enrichment,
    evidence_definitions,
    exact_si_species,
    gift_source_functional_counts,
    source_genus_positions,
)


def test_broad_direct_scope_reuses_frozen_broad_source_contract() -> None:
    marker = pd.DataFrame({"accepted_species": ["Alpha one"]})
    direct = pd.DataFrame({"accepted_species": ["Beta one"]})
    definitions = evidence_definitions(marker, direct, marker, direct)
    broad_direct = next(
        item for item in definitions if item["evidence_scope"] == "broad_direct"
    )
    assert broad_direct["species_scores"] is direct
    assert broad_direct["trait_ledger"] is None
    assert broad_direct["minimum_source_scored_species"] == 1
    assert broad_direct["minimum_represented_genera"] == [5]
    assert broad_direct["restricted_to_exact_si"] is False
    assert broad_direct["raw_gift_availability"] is True


def test_loading_increment_separates_species_loading_from_genus_entry() -> None:
    prevalence = np.array([1, 1, 1, 1], dtype=int)
    richness = np.ones(4, dtype=float)
    positions = np.array([-1.0, -0.5, 0.5, 1.0])

    neutral = compute_island_enrichment(
        prevalence,
        richness,
        np.array([1.0, 0.0, 0.0, 1.0]),
        positions,
        matching="prevalence_richness",
        minimum_represented_genera=2,
    )
    assert neutral is not None
    assert neutral["entry_enrichment"] == 0.0
    assert neutral["species_enrichment"] == 0.0
    assert neutral["loading_increment"] == 0.0

    loaded = compute_island_enrichment(
        prevalence,
        richness,
        np.array([1.0, 0.0, 0.0, 3.0]),
        positions,
        matching="prevalence_richness",
        minimum_represented_genera=2,
    )
    assert loaded is not None
    assert loaded["entry_enrichment"] == 0.0
    assert loaded["species_enrichment"] == 0.5
    assert loaded["loading_increment"] == 0.5


def test_source_prevalence_matching_absorbs_supply_only_difference() -> None:
    result = compute_island_enrichment(
        np.array([1, 2], dtype=int),
        np.array([1.0, 1.0]),
        np.array([1.0, 1.0]),
        np.array([-1.0, 1.0]),
        matching="prevalence_only",
        minimum_represented_genera=2,
    )
    assert result is not None
    assert result["entry_enrichment"] == 0.0
    assert result["species_enrichment"] == 0.0


def test_source_position_uses_gift_mainland_species_only() -> None:
    scores = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "syndrome": "large_bee_like",
                "syndrome_concordance": 1.0,
            },
            {
                "accepted_species": "Alpha one",
                "syndrome": "generalized_accessible",
                "syndrome_concordance": -1.0,
            },
            {
                "accepted_species": "Alpha two",
                "syndrome": "large_bee_like",
                "syndrome_concordance": -1.0,
            },
            {
                "accepted_species": "Alpha two",
                "syndrome": "generalized_accessible",
                "syndrome_concordance": 1.0,
            },
        ]
    )
    gift = pd.DataFrame([{"entity_ID": 1, "work_species": "Alpha two"}])
    positions, _, _ = source_genus_positions(gift, scores)
    assert len(positions) == 1
    assert positions.loc[0, "genus"] == "Alpha"
    assert positions.loc[0, "source_functional_position"] == 1.0
    assert positions.loc[0, "n_source_scored_species"] == 1


def test_gift_source_counts_separate_memberships_from_unique_species() -> None:
    matched = pd.DataFrame(
        [
            {"entity_ID": 1, "accepted_species": "Alpha one"},
            {"entity_ID": 2, "accepted_species": "Alpha one"},
            {"entity_ID": 2, "accepted_species": "Beta one"},
            {"entity_ID": 3, "accepted_species": "Gamma one"},
        ]
    )
    functional = pd.DataFrame(
        [
            {"accepted_species": "Alpha one"},
            {"accepted_species": "Beta one"},
        ]
    )

    assert gift_source_functional_counts(matched, functional) == {
        "n_gift_source_functional_memberships": 3,
        "n_unique_gift_source_functional_species": 2,
    }


def test_exact_si_species_fails_closed_on_exact_state() -> None:
    ledger = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "trait_name": "self_incompatibility",
                "normalized_value": "SI",
            },
            {
                "accepted_species": "Beta one",
                "trait_name": "self_incompatibility",
                "normalized_value": "SC",
            },
            {
                "accepted_species": "Gamma one",
                "trait_name": "self_incompatibility",
                "normalized_value": "SC|SI",
            },
        ]
    )
    assert exact_si_species(ledger) == {"Alpha one"}
