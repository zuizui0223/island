from __future__ import annotations

import pandas as pd
import pytest
from scipy import sparse

from island_v2.chapter1_species_sorting_identifiability import (
    _decompose_row,
    audit_identifiability,
)


def test_species_sorting_identity_and_prevalence_matching() -> None:
    # Source candidates:
    # A(state 0), B(state 1) occur in one source each;
    # C(state 1) occurs in two sources.
    # Island contains B and C. Prevalence-matched expectation is
    # mean([mean(A,B), mean(C)]) = mean([0.5, 1.0]) = 0.75.
    states = pd.Series([0.0, 1.0, 1.0]).to_numpy()
    prevalence = sparse.csr_matrix([[1, 1, 2]])
    observed = sparse.csr_matrix([[0, 1, 1]])
    result = _decompose_row(
        prevalence_row=prevalence,
        observed_row=observed,
        states=states,
        n_species=3,
    )
    assert result is not None
    assert result["raw_h1_mean"] == pytest.approx(1.0)
    assert result["source_species_expectation"] == pytest.approx(0.75)
    assert result["species_sorting_enrichment"] == pytest.approx(0.25)
    assert result["identity_error"] == pytest.approx(0.0)
    assert result["source_overlap_fraction"] == pytest.approx(1.0)


def test_source_unavailable_observed_species_is_not_silently_used() -> None:
    states = pd.Series([0.0, 1.0, 1.0]).to_numpy()
    prevalence = sparse.csr_matrix([[1, 1, 0]])
    observed = sparse.csr_matrix([[0, 1, 1]])
    result = _decompose_row(
        prevalence_row=prevalence,
        observed_row=observed,
        states=states,
        n_species=3,
    )
    assert result is not None
    assert result["n_observed_trait_species"] == 2
    assert result["n_observed_source_candidate_species"] == 1
    assert result["source_overlap_fraction"] == pytest.approx(0.5)
    assert result["raw_h1_mean"] == pytest.approx(1.0)


def test_identifiability_fails_closed_without_population_trait_axis() -> None:
    audit = pd.DataFrame(
        {
            "accepted_species": ["A a", "B b"],
            "trait_name": ["self_incompatibility", "self_incompatibility"],
            "resolved_for_primary": [True, True],
            "canonical_signature": ["SC", "SI"],
        }
    )
    ledger = pd.DataFrame(
        {
            "accepted_species": ["A a", "B b"],
            "trait_name": ["self_incompatibility", "self_incompatibility"],
            "normalized_value": ["SC", "SI"],
        }
    )
    flora = pd.DataFrame(
        {
            "island_id": ["i1", "i2", "i1", "i2"],
            "accepted_species": ["A a", "A a", "B b", "B b"],
            "floristic_status": ["native_nonendemic"] * 4,
        }
    )
    gift = pd.DataFrame(
        {
            "entity_ID": [1, 2],
            "accepted_species": ["A a", "B b"],
        }
    )
    config = {
        "scope": {"floristic_status": "native_nonendemic"},
        "identifiability": {
            "current_trait_unit": "accepted_species",
            "required_for_within_species_change": [
                "population_or_locality_identifier",
                "island_vs_source_population_role",
                "population_specific_trait_state",
            ],
        },
    }
    result, candidates = audit_identifiability(
        state_audit_all=audit,
        state_audit_direct=audit,
        trait_ledger_all=ledger,
        trait_ledger_direct=ledger,
        status_flora=flora,
        matched_gift=gift,
        config=config,
    )
    assert result["within_species_change_estimable"] is False
    assert result["within_species_change_zero"] is False
    assert result["n_candidate_repeated_lineages_for_future_population_sampling"] == 2
    assert set(candidates["accepted_species"]) == {"A a", "B b"}
