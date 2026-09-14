from pathlib import Path

import pandas as pd
import yaml

from island_v2.chapter1_effect_fingerprint import run_synthesis


OUTCOMES = [
    "plain_colour",
    "generalized_form",
    "actinomorphic_symmetry",
    "shallow_open_tube",
    "small_flower",
    "self_compatibility",
    "selfing_mating_system",
    "autonomous_selfing",
]


def _config() -> dict:
    return {
        "contract": "chapter1_effect_fingerprint_v1",
        "pinned_input": {
            "workflow_run_id": 1,
            "artifact_id": 2,
            "digest": "sha256:test",
        },
        "primary_taxonomic_axes": ["generalized_accessible", "selfing_core"],
        "primary_taxonomic_context": {
            "context_layer": "biogeographic_realm",
            "context": "Palearctic",
        },
        "stages": ["observed_score", "after_family_residual", "after_genus_residual"],
        "strata": ["all_native", "native_nonendemic"],
        "source_modes": ["geo_k5"],
        "taxonomic_attenuation": {
            "do_not_compute_ratio_when_observed_norm_below": 1e-12,
        },
        "atomic_fingerprint": {
            "contexts": ["northern_midlatitude", "tropical"],
            "support_tier": "confirmatory",
            "outcome_order": OUTCOMES,
            "outcome_domains": {
                "flower_colour": ["plain_colour"],
                "floral_structural_complexity": [
                    "generalized_form",
                    "actinomorphic_symmetry",
                    "shallow_open_tube",
                    "small_flower",
                ],
                "reproductive_assurance": [
                    "self_compatibility",
                    "selfing_mating_system",
                    "autonomous_selfing",
                ],
            },
        },
        "claim_boundary": "test boundary",
    }


def _write_taxonomic(root: Path, scope: str) -> None:
    out = root / "taxonomic-depth" / scope
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    stage_slopes = {
        "observed_score": (0.3, 0.4),
        "after_family_residual": (0.18, 0.24),
        "after_genus_residual": (0.03, 0.04),
    }
    for stratum in ("all_native", "native_nonendemic"):
        for stage, slopes in stage_slopes.items():
            for axis, slope in zip(("generalized_accessible", "selfing_core"), slopes, strict=True):
                rows.append(
                    {
                        "source_mode": "geo_k5",
                        "context_layer": "biogeographic_realm",
                        "axis_set": f"taxonomic_stage__{stage}",
                        "stratum": stratum,
                        "support_tier": "confirmatory",
                        "context": "Palearctic",
                        "syndrome": f"{stage}__{axis}",
                        "distance_slope": slope,
                        "cluster_robust_se": 0.01,
                        "p_value": 0.01,
                        "n_islands": 100,
                    }
                )
    pd.DataFrame(rows).to_csv(out / "slopes.csv", index=False)


def _write_atomic(root: Path, short_scope: str) -> None:
    out = root / "atomic" / short_scope
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for stratum in ("all_native", "native_nonendemic"):
        for context in ("northern_midlatitude", "tropical"):
            for index, outcome in enumerate(OUTCOMES):
                base = 0.05 if index % 2 == 0 else -0.05
                slope = base if context == "northern_midlatitude" else -base
                rows.append(
                    {
                        "stratum": stratum,
                        "support_tier": "confirmatory",
                        "context": context,
                        "outcome": outcome,
                        "geography_slope_log_odds_per_response_sd": slope,
                        "cluster_robust_se": 0.01,
                        "p_value": 0.01,
                        "n_islands": 80,
                    }
                )
    pd.DataFrame(rows).to_csv(out / "observed_within_outcome_slopes.csv", index=False)


def test_effect_fingerprint_builds_attenuation_and_domains(tmp_path: Path) -> None:
    root = tmp_path / "artifact"
    _write_taxonomic(root, "all_analysis_eligible")
    _write_taxonomic(root, "direct_only")
    _write_atomic(root, "all")
    _write_atomic(root, "direct")
    config_path = tmp_path / "config.yml"
    config_path.write_text(yaml.safe_dump(_config()), encoding="utf-8")
    output = tmp_path / "output"

    manifest = run_synthesis(
        artifact_root=root,
        config_path=config_path,
        output_dir=output,
    )

    attenuation = pd.read_csv(output / "taxonomic_vector_attenuation.csv")
    assert len(attenuation) == 4
    assert attenuation["genus_attenuation_fraction"].between(0.89, 0.91).all()
    assert attenuation["conditional_genus_attenuation_fraction"].between(0.83, 0.84).all()

    fingerprint = pd.read_csv(output / "atomic_response_fingerprint.csv")
    assert set(fingerprint["domain"]) == {
        "flower_colour",
        "floral_structural_complexity",
        "reproductive_assurance",
    }
    assert set(fingerprint["direction"]) == {"positive", "negative"}

    summary = pd.read_csv(output / "atomic_domain_fingerprint_summary.csv")
    assert len(summary) == 2 * 2 * 2 * 3

    geometry = pd.read_csv(output / "atomic_cross_context_vector_geometry.csv")
    assert len(geometry) == 4
    assert geometry["complete_vector"].all()
    assert geometry["vector_angle_degrees"].between(179.9, 180.0).all()

    assert manifest["new_p_values_generated"] is False
    assert manifest["mechanism_promoted"] is False
