from __future__ import annotations

import pandas as pd
import pytest

from island_v2.chapter1_three_level_filter_evolution import audit_ledger


def _row(
    *,
    system: str,
    taxon: str,
    direction: str,
    mechanism: str,
    persistent: str = "no",
    post: str = "no",
) -> dict[str, str]:
    return {
        "system": system,
        "taxon": taxon,
        "study": "study",
        "doi": "10.test/example",
        "island_mainland_design": "design",
        "population_trait_measurement": "measurement",
        "persistent_under_common_environment_or_genetic_marker": persistent,
        "observed_direction": direction,
        "mechanism_identified": mechanism,
        "post_colonization_evolution_identified": post,
        "role": "role",
    }


def test_three_level_audit_counts_null_and_directional_population_rows() -> None:
    table = pd.DataFrame(
        [
            _row(
                system="current_repo",
                taxon="global",
                direction="genus composition positive",
                mechanism="assemblage_filtering_no_demographic_subprocess",
            ),
            _row(
                system="population_pilot",
                taxon="A",
                direction="island mean differs numerically",
                mechanism="no_uniform_island_shift",
            ),
            _row(
                system="external_validation",
                taxon="B",
                direction="island SC; mainland SI",
                mechanism="within_species_genetic_divergence_consistent_with_colonization_filter",
                persistent="genetic_marker",
            ),
        ]
    )
    result = audit_ledger(table)
    assert result["n_rows"] == 3
    assert result["n_null_or_nonuniform_within_species_rows"] == 1
    assert result["n_rows_with_directional_same_species_divergence"] == 1
    assert result["n_rows_with_common_environment_or_genetic_support"] == 1
    assert result["n_post_colonization_evolution_identified"] == 0
    assert result["three_level_status"]["assemblage_filtering_identified"] is True
    assert (
        result["three_level_status"][
            "founder_filtering_vs_post_colonization_evolution_separated"
        ]
        is False
    )


def test_post_colonization_claim_fails_without_persistent_evidence() -> None:
    table = pd.DataFrame(
        [
            _row(
                system="external_validation",
                taxon="A",
                direction="island shift",
                mechanism="post_colonization_change_identified",
                persistent="no",
                post="yes",
            )
        ]
    )
    with pytest.raises(ValueError, match="lacks persistent/genetic evidence"):
        audit_ledger(table)


def test_post_colonization_claim_fails_without_explicit_mechanism() -> None:
    table = pd.DataFrame(
        [
            _row(
                system="external_validation",
                taxon="A",
                direction="island shift",
                mechanism="within_species_divergence",
                persistent="genetic_marker",
                post="yes",
            )
        ]
    )
    with pytest.raises(ValueError, match="lacks explicit mechanism"):
        audit_ledger(table)
