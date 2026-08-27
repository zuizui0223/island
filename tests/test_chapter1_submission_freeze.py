import json
from pathlib import Path

import pandas as pd
import pytest

from island_v2.chapter1_submission_freeze import validate_submission_freeze


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


def _build_fixture(tmp_path: Path) -> tuple[Path, Path]:
    root = tmp_path / "artifact"
    result = root / "result"
    result.mkdir(parents=True)
    common = {
        "source_mode": "not_applicable",
        "support_tier": "confirmatory",
        "threshold": 50,
    }
    coefficients = []
    for stratum in ["all_native", "native_nonendemic"]:
        for source_mode in ["geo_k5", "geo_k10", "geo_k20", "geo50_climate10"]:
            coefficients.append(
                {
                    **common,
                    "family": "source_lineage_broad",
                    "source_mode": source_mode,
                    "context_layer": "biogeographic_realm",
                    "stratum": stratum,
                    "context": "Palearctic",
                    "response": "entry_enrichment",
                    "n_islands": 141,
                    "distance_q": 0.01,
                    "distance_x_area_q": 0.01,
                    "area_moderation_state": "distance_effect_stronger_on_smaller_islands",
                }
            )
        coefficients.append(
            {
                **common,
                "family": "plant_side_branches",
                "context_layer": "biogeographic_realm",
                "stratum": stratum,
                "context": "Neotropical",
                "response": "reproductive_assurance",
                "n_islands": 59,
                "distance_q": 0.01,
                "distance_x_area_q": 0.01,
                "area_moderation_state": "distance_effect_stronger_on_smaller_islands",
            }
        )
    coefficients = pd.DataFrame(coefficients)
    coefficients.to_csv(result / "area_moderation_coefficients.csv", index=False)

    within = pd.DataFrame(
        [
            {
                **common,
                "family": "plant_side_branches",
                "context_layer": "biogeographic_realm",
                "stratum": stratum,
                "context": "Neotropical",
                "n_clusters": 14,
                "q_vector_family": 0.01,
            }
            for stratum in ["all_native", "native_nonendemic"]
        ]
    )
    within.to_csv(result / "area_moderation_within_context.csv", index=False)

    between_coefficients = []
    for stratum in ["all_native", "native_nonendemic"]:
        between_coefficients.append(
            {
                **common,
                "family": "plant_side_branches",
                "context_layer": "analysis_regime",
                "stratum": stratum,
                "context_a": "northern_midlatitude",
                "context_b": "tropical",
                "response": "accessibility_generalization",
                "difference_q": 0.01,
                "difference_axis_supported": True,
            }
        )
        for context_layer, supported in [
            ("analysis_regime", True),
            ("biogeographic_realm", stratum == "all_native"),
        ]:
            between_coefficients.append(
                {
                    **common,
                    "family": "source_lineage_broad",
                    "source_mode": f"{stratum}_{context_layer}",
                    "context_layer": context_layer,
                    "stratum": stratum,
                    "context_a": "north",
                    "context_b": "tropical",
                    "response": "entry_enrichment",
                    "difference_q": 0.01 if supported else 0.5,
                    "difference_axis_supported": supported,
                }
            )
    between_coefficients = pd.DataFrame(between_coefficients)
    between_coefficients.to_csv(
        result / "area_moderation_between_coefficients.csv", index=False
    )
    between = between_coefficients.drop_duplicates(
        [
            "family",
            "source_mode",
            "context_layer",
            "stratum",
            "support_tier",
            "threshold",
            "context_a",
            "context_b",
        ]
    ).copy()
    between["q_between_family"] = 0.01
    between.to_csv(result / "area_moderation_between_context.csv", index=False)

    _write_json(
        result / "area_moderation_manifest.json",
        {
            "contract": "chapter1_area_capacity_moderation_v1",
            "continuous_area_no_threshold": True,
            "joint_vector_gate_before_axis_classification": True,
            "n_coefficient_rows": len(coefficients),
            "n_within_rows": len(within),
            "n_between_coefficient_rows": len(between_coefficients),
            "n_between_rows": len(between),
            "claim_ceiling": {
                "pollinator_mobility_observed": False,
                "effective_pollination_service_observed": False,
            },
        },
    )
    _write_json(
        root / "frozen/result/lineage_representation_manifest.json",
        {
            "contract": "chapter1_pr138_lineage_representation_bridge_v1",
            "claim_boundary": (
                "No pollinator loss, historical source ancestry, or in-situ evolution."
            ),
        },
    )
    for relative in [
        "frozen/result/lineage_representation_island_scores.csv.gz",
        "frozen/frozen/input/canonical/chapter1_checkpoint.json",
        "frozen/frozen/input/input/pathway/syndrome/island_syndrome_scores.csv.gz",
        "frozen/frozen/input/input/isolation/results/purpose_shortest_island_data.csv",
    ]:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"fixture")

    lock = tmp_path / "lock.json"
    _write_json(
        lock,
        {
            "contract": "chapter1_submission_freeze_lock_v1",
            "workflow_run_id": 1,
            "artifact_id": 2,
            "artifact_digest": "sha256:test",
            "expires_at": "2026-11-25T00:00:00Z",
            "verified_summary": {
                "area_moderation_contract": "chapter1_area_capacity_moderation_v1",
                "n_coefficient_rows": len(coefficients),
                "n_within_rows": len(within),
                "n_between_coefficient_rows": len(between_coefficients),
                "n_between_rows": len(between),
                "palearctic_broad_entry_smaller_island_scenarios": 8,
                "palearctic_broad_entry_total_scenarios": 8,
                "broad_regime_entry_difference_supported_scenarios": 2,
                "formal_realm_entry_difference_supported_scenarios": 1,
                "neotropical_reproductive_assurance_smaller_island_strata": 2,
                "north_tropical_accessibility_difference_supported_strata": 2,
                "neotropical_reproductive_assurance_min_islands": 59,
                "neotropical_reproductive_assurance_spatial_blocks": 14,
            },
            "claim_ceiling": {
                "pollinator_mobility_observed": False,
                "effective_pollination_service_observed": False,
            },
        },
    )
    return lock, root


def test_submission_freeze_fixture_passes(tmp_path: Path) -> None:
    lock, root = _build_fixture(tmp_path)
    report = validate_submission_freeze(lock, root)
    assert report["verified"] is True
    assert report["headline"]["palearctic_broad_entry_smaller_island_scenarios"] == 8
    assert report["archive_status"] == "temporary_actions_retention_not_durable_archive"


def test_submission_freeze_fails_if_causal_ceiling_is_raised(tmp_path: Path) -> None:
    lock, root = _build_fixture(tmp_path)
    manifest_path = root / "result/area_moderation_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["claim_ceiling"]["pollinator_mobility_observed"] = True
    _write_json(manifest_path, manifest)
    with pytest.raises(ValueError, match="claim ceiling was raised"):
        validate_submission_freeze(lock, root)


def test_submission_freeze_fails_if_artifact_file_is_missing(tmp_path: Path) -> None:
    lock, root = _build_fixture(tmp_path)
    (root / "result/area_moderation_between_context.csv").unlink()
    with pytest.raises(ValueError, match="missing files"):
        validate_submission_freeze(lock, root)
