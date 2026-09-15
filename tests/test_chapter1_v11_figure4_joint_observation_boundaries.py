import json
from pathlib import Path

import pandas as pd

from island_v2.chapter1_v11_figure4_joint_observation_boundaries import (
    _primary_surface_matrix,
    load_inputs,
)


def _joint_surface() -> pd.DataFrame:
    rows = []
    for target in ("Palearctic_vector", "tropical_vector"):
        for or_r in (0.5, 1.0, 2.0):
            for or_d in (0.5, 1.0):
                for c0 in (0.5, 0.95):
                    for or_c in (0.25, 1.0):
                        rows.append(
                            {
                                "evidence_scope": "direct_only",
                                "target": target,
                                "stratum": "native_nonendemic",
                                "surface_type": "joint_selection_grid",
                                "status": "fit",
                                "trait_resolution_odds_ratio": or_r,
                                "median_distance_completeness": c0,
                                "distance_completeness_odds_ratio": or_c,
                                "state_recording_odds_ratio": or_d,
                                "robust_cell": not (target == "tropical_vector" and or_d == 0.5),
                            }
                        )
    return pd.DataFrame(rows)


def test_primary_surface_matrix_collapses_only_c0_orc() -> None:
    surface = _joint_surface()
    pal, ors_r, ors_d = _primary_surface_matrix(surface, "Palearctic_vector")
    trop, _, _ = _primary_surface_matrix(surface, "tropical_vector")
    assert pal.shape == (3, 2)
    assert ors_r == [0.5, 1.0, 2.0]
    assert ors_d == [0.5, 1.0]
    assert (pal == 1.0).all()
    assert (trop[:, 0] == 0.0).all()
    assert (trop[:, 1] == 1.0).all()


def test_load_inputs_accepts_frozen_joint_and_boundary_shapes(tmp_path: Path) -> None:
    joint = tmp_path / "joint"
    geom = tmp_path / "geom"
    h5c = tmp_path / "h5c"
    h5d = tmp_path / "h5d"
    for path in (joint, geom, h5c, h5d):
        path.mkdir()

    _joint_surface().to_csv(
        joint / "joint_surface_classification.csv.gz", index=False, compression="gzip"
    )
    pd.DataFrame(
        [
            {
                "evidence_scope": "direct_only",
                "target": "Palearctic_vector",
                "stratum": "native_nonendemic",
                "n_fit_cells": 24,
                "n_robust_cells": 24,
                "robust_fraction_of_fit_grid": 1.0,
            }
        ]
    ).to_csv(joint / "joint_robustness_summary.csv", index=False)
    envelope_rows = []
    for target, lo, hi, robust in (
        ("Palearctic_accessibility", 0.01, 0.05, True),
        ("tropical_accessibility", -0.2, 0.02, False),
        ("Palearctic_vector", None, None, True),
        ("tropical_vector", None, None, False),
        ("north_tropical_vector_difference", None, None, False),
    ):
        envelope_rows.append(
            {
                "evidence_scope": "direct_only",
                "target": target,
                "stratum": "native_nonendemic",
                "estimate_lower": lo,
                "estimate_upper": hi,
                "expected_sign_identified": robust,
                "support_identified_across_envelope": robust,
                "envelope_robust": robust,
            }
        )
    pd.DataFrame(envelope_rows).to_csv(
        joint / "partial_identification_envelope.csv", index=False
    )
    (joint / "joint_observation_bias_manifest.json").write_text(
        json.dumps(
            {
                "contract": "chapter1_joint_observation_bias_v1",
                "n_primary_parameter_surfaces_per_scope": 1575,
                "grid_fraction_is_probability": False,
            }
        )
    )

    pd.DataFrame(
        [
            {
                "all_observed_D": 0,
                "all_critical_D": 1,
                "direct_observed_D": 0,
                "direct_critical_D": 1,
                "classification": "monotonic_or_unresolved",
            }
            for _ in range(12)
        ]
    ).to_csv(geom / "observed_geometry_cross_scope.csv", index=False)
    (h5c / "h5c_observed_result.json").write_text(
        json.dumps(
            {
                "classification": "no_pollination_mode_specificity_support",
                "interaction_estimate": 0.06,
                "interaction_ci_low": -0.09,
                "interaction_ci_high": 0.22,
                "interaction_p_value": 0.41,
            }
        )
    )
    pd.DataFrame(
        [
            {
                "classification_accuracy": 0.75,
                "false_distributed_selection_under_smooth": 0.2,
                "qualified": False,
            }
            for _ in range(8)
        ]
    ).to_csv(h5d / "h5d_identifiability_summary.csv", index=False)

    loaded = load_inputs(joint, geom, h5c, h5d)
    assert len(loaded["geometry"]) == 12
    assert len(loaded["h5d"]) == 8
    assert loaded["manifest"]["grid_fraction_is_probability"] is False


# P3 renderer validation is intentionally driven by frozen artifacts in CI.
