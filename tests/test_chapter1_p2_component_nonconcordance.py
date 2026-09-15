from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.chapter1_p2_component_nonconcordance import (
    _angle_degrees,
    _build_coobserved_species_scores,
    _common_island_scores,
    _determinant,
)


def test_common_island_scores_keeps_only_islands_with_both_axes() -> None:
    frame = pd.DataFrame(
        [
            {"island_id": "a", "stratum": "all_native", "syndrome": "accessibility_generalization", "syndrome_score": 1.0, "n_species": 5},
            {"island_id": "a", "stratum": "all_native", "syndrome": "reproductive_assurance", "syndrome_score": 2.0, "n_species": 4},
            {"island_id": "b", "stratum": "all_native", "syndrome": "accessibility_generalization", "syndrome_score": 1.0, "n_species": 5},
        ]
    )
    out = _common_island_scores(
        frame,
        axes=["accessibility_generalization", "reproductive_assurance"],
        stratum="all_native",
    )
    assert set(out["island_id"]) == {"a"}
    assert out["syndrome"].nunique() == 2


def test_vector_geometry_detects_noncollinearity() -> None:
    a = np.array([1.0, 1.0])
    b = np.array([-1.0, 1.0])
    assert _determinant(a, b) == 2.0
    assert np.isclose(_angle_degrees(a, b), 90.0)
    c = 2.0 * a
    assert np.isclose(_determinant(a, c), 0.0)
    assert abs(_angle_degrees(a, c)) < 1e-5


def test_coobserved_species_uses_same_species_denominator() -> None:
    species = pd.DataFrame(
        [
            {"accepted_species": "sp1", "syndrome": "generalized_accessible", "syndrome_concordance": 1.0},
            {"accepted_species": "sp1", "syndrome": "selfing_core", "syndrome_concordance": 0.0},
            {"accepted_species": "sp2", "syndrome": "generalized_accessible", "syndrome_concordance": 0.0},
            {"accepted_species": "sp2", "syndrome": "selfing_core", "syndrome_concordance": 1.0},
            {"accepted_species": "sp3", "syndrome": "generalized_accessible", "syndrome_concordance": 1.0},
        ]
    )
    flora = pd.DataFrame(
        [
            {"island_id": "i1", "accepted_species": "sp1", "origin_status": "native", "endemic_status": "nonendemic", "floristic_status": "native_nonendemic"},
            {"island_id": "i1", "accepted_species": "sp2", "origin_status": "native", "endemic_status": "nonendemic", "floristic_status": "native_nonendemic"},
            {"island_id": "i1", "accepted_species": "sp3", "origin_status": "native", "endemic_status": "nonendemic", "floristic_status": "native_nonendemic"},
        ]
    )
    out, diag = _build_coobserved_species_scores(
        species,
        flora,
        source_axes={
            "accessibility_generalization": "generalized_accessible",
            "reproductive_assurance": "selfing_core",
        },
        strata=["all_native", "native_nonendemic"],
        minimum_species=2,
    )
    assert diag["n_coobserved_species"] == 2
    subset = out.loc[out["stratum"].eq("native_nonendemic")]
    assert set(subset["n_species"]) == {2}
    assert subset["syndrome"].nunique() == 2
