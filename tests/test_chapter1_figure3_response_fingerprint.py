from pathlib import Path

import pandas as pd

from island_v2.chapter1_figure3_response_fingerprint import (
    CONTEXTS,
    OUTCOME_ORDER,
    SCOPES,
    STRATA,
    load_inputs,
    render_figure,
)


def _write_inputs(root: Path) -> None:
    atomic_rows = []
    for scope_index, scope in enumerate(SCOPES):
        for stratum in STRATA:
            for context_index, context in enumerate(CONTEXTS):
                for outcome_index, outcome in enumerate(OUTCOME_ORDER):
                    estimate = (0.01 + 0.005 * outcome_index) * (
                        1 if context_index == 0 else -1
                    )
                    estimate += 0.002 * scope_index
                    atomic_rows.append(
                        {
                            "evidence_scope": scope,
                            "stratum": stratum,
                            "context": context,
                            "outcome": outcome,
                            "geography_slope_log_odds_per_response_sd": estimate,
                            "ci_low": estimate - 0.02,
                            "ci_high": estimate + 0.02,
                            "n_islands": 80,
                        }
                    )
    pd.DataFrame(atomic_rows).to_csv(root / "atomic_response_fingerprint.csv", index=False)

    angle_rows = []
    for scope_index, scope in enumerate(SCOPES):
        for stratum_index, stratum in enumerate(STRATA):
            angle_rows.append(
                {
                    "evidence_scope": scope,
                    "stratum": stratum,
                    "context_a": CONTEXTS[0],
                    "context_b": CONTEXTS[1],
                    "n_components": 8,
                    "cosine_similarity": -0.1,
                    "vector_angle_degrees": 95.0 + 5 * scope_index + stratum_index,
                    "complete_vector": True,
                }
            )
    pd.DataFrame(angle_rows).to_csv(
        root / "atomic_cross_context_vector_geometry.csv", index=False
    )

    attenuation_rows = []
    for scope in SCOPES:
        for stratum in STRATA:
            for source_mode in ("geo_k5", "geo_k10", "geo_k20", "geo50_climate10"):
                attenuation_rows.append(
                    {
                        "evidence_scope": scope,
                        "source_mode": source_mode,
                        "stratum": stratum,
                        "context_layer": "biogeographic_realm",
                        "context": "Palearctic",
                        "observed_score_vector_norm": 0.05,
                        "after_family_residual_vector_norm": 0.035,
                        "after_genus_residual_vector_norm": 0.01,
                        "family_attenuation_fraction": 0.30,
                        "genus_attenuation_fraction": 0.80,
                        "conditional_genus_attenuation_fraction": 0.7142857,
                    }
                )
    pd.DataFrame(attenuation_rows).to_csv(
        root / "taxonomic_vector_attenuation.csv", index=False
    )


def test_figure3_renders_from_complete_frozen_design(tmp_path: Path) -> None:
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()
    _write_inputs(input_dir)

    tables = load_inputs(input_dir)
    manifest = render_figure(
        tables,
        output_dir,
        source_run_id=10,
        source_artifact_id=20,
        source_artifact_digest="sha256:test",
    )

    assert manifest["new_biological_models_fitted"] is False
    assert manifest["new_p_values_generated"] is False
    assert manifest["panel_a"]["n_atomic_rows"] == 64
    assert manifest["panel_b"]["n_vector_angle_rows"] == 4
    assert manifest["panel_c"]["n_taxonomic_profiles"] == 16
    for name in manifest["outputs"]:
        path = output_dir / name
        assert path.is_file()
        assert path.stat().st_size > 1_000


def test_figure3_fails_closed_when_one_atomic_cell_is_missing(tmp_path: Path) -> None:
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    _write_inputs(input_dir)
    atomic_path = input_dir / "atomic_response_fingerprint.csv"
    atomic = pd.read_csv(atomic_path).iloc[:-1].copy()
    atomic.to_csv(atomic_path, index=False)

    try:
        load_inputs(input_dir)
    except ValueError as exc:
        assert "design cells mismatch" in str(exc)
    else:
        raise AssertionError("missing frozen atomic design cell must fail closed")
