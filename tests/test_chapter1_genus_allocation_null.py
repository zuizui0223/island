from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from island_v2.chapter1_genus_allocation_null import (
    _allocation_share,
    _csr_slots,
    _estimability,
    _linear_slope_coefficients,
    _match_observed_to_candidate,
    _model_slope_weights,
    _shuffle_groups,
)
from island_v2.chapter1_species_sorting_identifiability import _decompose_row


def test_linear_coefficients_match_brute_decomposition() -> None:
    states = np.array([0.0, 1.0, 1.0, 0.0])
    genus_codes = np.array([0, 0, 1, 1], dtype=np.int32)
    prevalence = sparse.csr_matrix(
        np.array(
            [
                [1, 1, 2, 0],
                [1, 2, 1, 1],
                [0, 1, 0, 2],
            ],
            dtype=np.int16,
        )
    )
    observed = sparse.csr_matrix(
        np.array(
            [
                [0, 1, 1, 0],
                [1, 0, 1, 1],
                [0, 1, 0, 0],
            ],
            dtype=np.int8,
        )
    )
    slope_weights = np.array([0.2, -0.3, 0.1])

    ci, cs, cp = _csr_slots(prevalence)
    oi_all, os_all, _ = _csr_slots(observed)
    oi, os, op = _match_observed_to_candidate(
        candidate_island=ci,
        candidate_species=cs,
        candidate_prevalence=cp.astype(np.int16),
        observed_island=oi_all,
        observed_species=os_all,
        n_species=len(states),
    )
    coeff = _linear_slope_coefficients(
        n_species=len(states),
        n_islands=3,
        candidate_island=ci,
        candidate_species=cs,
        candidate_prevalence=cp.astype(np.int16),
        observed_island=oi,
        observed_species=os,
        observed_prevalence=op,
        genus_codes=genus_codes,
        slope_weights=slope_weights,
    )

    brute = []
    for island in range(3):
        result = _decompose_row(
            prevalence_row=prevalence.getrow(island),
            observed_row=observed.getrow(island),
            states=states,
            genus_codes=genus_codes,
        )
        assert result is not None
        brute.append(result)

    for component, brute_name in (
        ("total_sorting", "total_species_sorting_enrichment"),
        ("genus_structure", "genus_species_structure_enrichment"),
        ("within_genus", "within_genus_species_sorting_enrichment"),
    ):
        expected = sum(
            slope_weights[i] * float(brute[i][brute_name])
            for i in range(3)
        )
        observed_beta = float(coeff[component] @ states)
        assert observed_beta == pytest.approx(expected, abs=1e-12)

    assert float(coeff["total_sorting"] @ states) == pytest.approx(
        float(coeff["genus_structure"] @ states)
        + float(coeff["within_genus"] @ states),
        abs=1e-12,
    )


def test_estimability_reports_only_variable_multi_species_cells() -> None:
    states = np.array([0.0, 1.0, 1.0, 0.0])
    genus_codes = np.array([0, 0, 1, 1], dtype=np.int32)
    prevalence = sparse.csr_matrix(
        np.array(
            [
                [1, 1, 2, 0],
                [1, 2, 1, 1],
                [0, 1, 0, 2],
            ],
            dtype=np.int16,
        )
    )
    observed = sparse.csr_matrix(
        np.array(
            [
                [0, 1, 1, 0],
                [1, 0, 1, 1],
                [0, 1, 0, 0],
            ],
            dtype=np.int8,
        )
    )

    ci, cs, cp = _csr_slots(prevalence)
    oi_all, os_all, _ = _csr_slots(observed)
    oi, os, op = _match_observed_to_candidate(
        candidate_island=ci,
        candidate_species=cs,
        candidate_prevalence=cp.astype(np.int16),
        observed_island=oi_all,
        observed_species=os_all,
        n_species=len(states),
    )
    coeff = _linear_slope_coefficients(
        n_species=len(states),
        n_islands=3,
        candidate_island=ci,
        candidate_species=cs,
        candidate_prevalence=cp.astype(np.int16),
        observed_island=oi,
        observed_species=os,
        observed_prevalence=op,
        genus_codes=genus_codes,
        slope_weights=np.array([0.2, -0.3, 0.1]),
    )
    audit = _estimability(
        states=states,
        candidate_species=cs,
        genus_inverse=coeff["genus_inverse"],
        genus_count=coeff["genus_count"],
        observed_genus_group=coeff["observed_genus_group"],
    )

    assert audit["observed_source_candidate_slots"] == 6
    assert audit["estimable_slots"] == 3
    assert audit["estimable_slot_fraction"] == pytest.approx(0.5)
    assert audit["observed_genus_prevalence_cells"] == 5
    assert audit["estimable_cells"] == 2
    assert audit["estimable_cell_fraction"] == pytest.approx(0.4)


def test_shuffle_groups_preserve_genus_positive_counts() -> None:
    states = np.array([0.0, 1.0, 1.0, 0.0, 1.0])
    genus_codes = np.array([0, 0, 1, 1, 1], dtype=np.int32)
    candidate_species = np.array([0, 1, 2, 3, 4], dtype=np.int32)
    groups = _shuffle_groups(
        candidate_species=candidate_species,
        genus_codes=genus_codes,
        states=states,
    )
    rng = np.random.default_rng(7)
    shuffled = states.copy()
    before = {
        int(g): float(states[genus_codes == g].sum())
        for g in np.unique(genus_codes)
    }
    for group in groups:
        shuffled[group] = rng.permutation(states[group])
    after = {
        int(g): float(shuffled[genus_codes == g].sum())
        for g in np.unique(genus_codes)
    }
    assert before == after


def test_model_weights_equal_full_ols_distance_coefficient() -> None:
    islands = [f"i{i}" for i in range(12)]
    x = np.arange(12, dtype=float)
    covariates = pd.DataFrame(
        {
            "island_id": islands,
            "log_island_area_km2": x + 1.0,
            "climate_pc1": np.sin(x),
            "climate_pc2": np.cos(x),
            "climate_pc3": (x % 3) - 1.0,
            "climate_pc4": (x**2) / 10.0,
            "log_distance_to_continent_km": np.log1p(x + 0.5),
        }
    )
    model = _model_slope_weights(
        covariates=covariates,
        islands=islands,
        n_observed=np.ones(12, dtype=int),
        minimum_islands=5,
    )
    assert model is not None
    weights, n = model
    assert n == 12

    response = 0.3 * x + np.sin(x / 2.0)

    def standardize(a: np.ndarray) -> np.ndarray:
        return (a - np.mean(a)) / np.std(a, ddof=0)

    design = np.column_stack(
        [
            np.ones(12),
            standardize(covariates["log_island_area_km2"].to_numpy()),
            standardize(covariates["climate_pc1"].to_numpy()),
            standardize(covariates["climate_pc2"].to_numpy()),
            standardize(covariates["climate_pc3"].to_numpy()),
            standardize(covariates["climate_pc4"].to_numpy()),
            standardize(covariates["log_distance_to_continent_km"].to_numpy()),
        ]
    )
    expected = float(np.linalg.lstsq(design, response, rcond=None)[0][-1])
    assert float(weights @ response) == pytest.approx(expected, abs=1e-12)


def test_absolute_allocation_share_is_bounded() -> None:
    assert _allocation_share(8.0, 2.0) == pytest.approx(0.8)
    assert _allocation_share(-8.0, 2.0) == pytest.approx(0.8)
    assert np.isnan(_allocation_share(0.0, 0.0))
