"""Validate the Chapter 1 submission-candidate artifact and claim ceiling."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

AREA_FILES = {
    "coefficients": "result/area_moderation_coefficients.csv",
    "within": "result/area_moderation_within_context.csv",
    "between_coefficients": "result/area_moderation_between_coefficients.csv",
    "between": "result/area_moderation_between_context.csv",
    "manifest": "result/area_moderation_manifest.json",
}

EMBEDDED_FILES = {
    "lineage_manifest": "frozen/result/lineage_representation_manifest.json",
    "lineage_scores": "frozen/result/lineage_representation_island_scores.csv.gz",
    "canonical_checkpoint": "frozen/frozen/input/canonical/chapter1_checkpoint.json",
    "island_scores": (
        "frozen/frozen/input/input/pathway/syndrome/island_syndrome_scores.csv.gz"
    ),
    "covariates": (
        "frozen/frozen/input/input/isolation/results/purpose_shortest_island_data.csv"
    ),
}


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise TypeError(f"expected JSON object: {path}")
    return value


def _as_bool(series: pd.Series) -> pd.Series:
    values = series.astype(str).str.casefold()
    unknown = set(values.unique()) - {"true", "false"}
    if unknown:
        raise ValueError(f"invalid boolean values: {sorted(unknown)}")
    return values.eq("true")


def _require_unique(table: pd.DataFrame, keys: list[str], label: str) -> None:
    missing = set(keys) - set(table.columns)
    if missing:
        raise ValueError(f"{label} missing key columns: {sorted(missing)}")
    if table.duplicated(keys).any():
        raise ValueError(f"{label} has duplicate keys: {keys}")


def _require_q_values(table: pd.DataFrame, columns: list[str], label: str) -> None:
    for column in columns:
        values = pd.to_numeric(table[column], errors="coerce").dropna()
        if not values.between(0, 1).all():
            raise ValueError(f"{label}.{column} contains values outside [0, 1]")


def _require_files(root: Path) -> dict[str, Path]:
    paths = {
        name: root / relative
        for name, relative in {**AREA_FILES, **EMBEDDED_FILES}.items()
    }
    missing = [str(path.relative_to(root)) for path in paths.values() if not path.is_file()]
    if missing:
        raise ValueError(f"submission artifact missing files: {missing}")
    return paths


def validate_submission_freeze(lock_path: Path, artifact_root: Path) -> dict[str, Any]:
    """Fail closed unless the downloaded artifact matches the frozen scientific contract."""
    lock = _load_json(lock_path)
    if lock.get("contract") != "chapter1_submission_freeze_lock_v1":
        raise ValueError("unexpected submission-freeze lock contract")
    expected = lock["verified_summary"]
    paths = _require_files(artifact_root)

    manifest = _load_json(paths["manifest"])
    if manifest.get("contract") != expected["area_moderation_contract"]:
        raise ValueError("area-moderation contract changed")
    if not manifest.get("continuous_area_no_threshold"):
        raise ValueError("continuous-area no-threshold contract was not retained")
    if not manifest.get("joint_vector_gate_before_axis_classification"):
        raise ValueError("axis classification no longer follows the joint-vector gate")

    for key, observed in manifest.get("claim_ceiling", {}).items():
        if key.endswith(("_observed", "_identified")) and observed is not False:
            raise ValueError(f"artifact claim ceiling was raised: {key}")
    if any(value is not False for value in lock["claim_ceiling"].values()):
        raise ValueError("submission lock contains an affirmative causal claim")

    coefficients = pd.read_csv(paths["coefficients"])
    within = pd.read_csv(paths["within"])
    between_coefficients = pd.read_csv(paths["between_coefficients"])
    between = pd.read_csv(paths["between"])
    tables = {
        "n_coefficient_rows": coefficients,
        "n_within_rows": within,
        "n_between_coefficient_rows": between_coefficients,
        "n_between_rows": between,
    }
    for count_key, table in tables.items():
        if len(table) != int(expected[count_key]):
            raise ValueError(
                f"{count_key} changed: observed {len(table)}, expected {expected[count_key]}"
            )
        if int(manifest[count_key]) != len(table):
            raise ValueError(f"manifest {count_key} does not match its CSV")

    common = [
        "family",
        "source_mode",
        "context_layer",
        "stratum",
        "support_tier",
        "threshold",
    ]
    _require_unique(
        coefficients,
        [*common, "context", "response"],
        "area_moderation_coefficients",
    )
    _require_unique(
        within,
        [*common, "context"],
        "area_moderation_within_context",
    )
    _require_unique(
        between_coefficients,
        [*common, "context_a", "context_b", "response"],
        "area_moderation_between_coefficients",
    )
    _require_unique(
        between,
        [*common, "context_a", "context_b"],
        "area_moderation_between_context",
    )
    _require_q_values(coefficients, ["distance_q", "distance_x_area_q"], "coefficients")
    _require_q_values(within, ["q_vector_family"], "within")
    _require_q_values(between_coefficients, ["difference_q"], "between_coefficients")
    _require_q_values(between, ["q_between_family"], "between")

    confirmatory = coefficients.loc[coefficients["support_tier"].eq("confirmatory")]
    pal_entry = confirmatory.loc[
        confirmatory["family"].eq("source_lineage_broad")
        & confirmatory["context_layer"].eq("biogeographic_realm")
        & confirmatory["context"].eq("Palearctic")
        & confirmatory["response"].eq("entry_enrichment")
    ]
    pal_small = pal_entry["area_moderation_state"].eq(
        "distance_effect_stronger_on_smaller_islands"
    )
    if len(pal_entry) != int(expected["palearctic_broad_entry_total_scenarios"]):
        raise ValueError("Palearctic broad genus-entry scenario count changed")
    if int(pal_small.sum()) != int(
        expected["palearctic_broad_entry_smaller_island_scenarios"]
    ):
        raise ValueError("Palearctic broad genus-entry area pattern changed")

    neo_reproductive = confirmatory.loc[
        confirmatory["family"].eq("plant_side_branches")
        & confirmatory["context_layer"].eq("biogeographic_realm")
        & confirmatory["context"].eq("Neotropical")
        & confirmatory["response"].eq("reproductive_assurance")
    ]
    neo_small = neo_reproductive["area_moderation_state"].eq(
        "distance_effect_stronger_on_smaller_islands"
    )
    if int(neo_small.sum()) != int(
        expected["neotropical_reproductive_assurance_smaller_island_strata"]
    ):
        raise ValueError("Neotropical reproductive-assurance area pattern changed")
    if int(neo_reproductive["n_islands"].min()) != int(
        expected["neotropical_reproductive_assurance_min_islands"]
    ):
        raise ValueError("Neotropical reproductive-assurance support changed")

    confirmatory_between = between_coefficients.loc[
        between_coefficients["support_tier"].eq("confirmatory")
    ]
    access_difference = confirmatory_between.loc[
        confirmatory_between["family"].eq("plant_side_branches")
        & confirmatory_between["context_layer"].eq("analysis_regime")
        & confirmatory_between["response"].eq("accessibility_generalization")
    ]
    access_supported = _as_bool(access_difference["difference_axis_supported"])
    if int(access_supported.sum()) != int(
        expected["north_tropical_accessibility_difference_supported_strata"]
    ):
        raise ValueError("north-tropical accessibility moderation difference changed")

    broad_entry_between = confirmatory_between.loc[
        confirmatory_between["family"].eq("source_lineage_broad")
        & confirmatory_between["response"].eq("entry_enrichment")
    ].copy()
    broad_entry_between["supported"] = _as_bool(
        broad_entry_between["difference_axis_supported"]
    )
    supported_by_layer = broad_entry_between.groupby("context_layer")["supported"].sum()
    if int(supported_by_layer.get("analysis_regime", 0)) != int(
        expected["broad_regime_entry_difference_supported_scenarios"]
    ):
        raise ValueError("broad-regime lineage-entry contrast changed")
    if int(supported_by_layer.get("biogeographic_realm", 0)) != int(
        expected["formal_realm_entry_difference_supported_scenarios"]
    ):
        raise ValueError("formal-realm lineage-entry contrast changed")

    neo_within = within.loc[
        within["support_tier"].eq("confirmatory")
        & within["family"].eq("plant_side_branches")
        & within["context_layer"].eq("biogeographic_realm")
        & within["context"].eq("Neotropical")
    ]
    if set(pd.to_numeric(neo_within["n_clusters"], errors="raise")) != {
        int(expected["neotropical_reproductive_assurance_spatial_blocks"])
    }:
        raise ValueError("Neotropical spatial-block support changed")

    lineage_manifest = _load_json(paths["lineage_manifest"])
    if lineage_manifest.get("contract") != "chapter1_pr138_lineage_representation_bridge_v1":
        raise ValueError("embedded lineage-representation contract changed")
    boundary = str(lineage_manifest.get("claim_boundary", "")).casefold()
    for forbidden in ["pollinator loss", "historical source ancestry", "in-situ evolution"]:
        if forbidden not in boundary:
            raise ValueError(f"embedded lineage claim boundary no longer excludes {forbidden}")

    return {
        "contract": "chapter1_submission_freeze_verification_v1",
        "verified": True,
        "workflow_run_id": int(lock["workflow_run_id"]),
        "artifact_id": int(lock["artifact_id"]),
        "artifact_digest": str(lock["artifact_digest"]),
        "artifact_expires_at": str(lock["expires_at"]),
        "archive_status": "temporary_actions_retention_not_durable_archive",
        "row_counts": {key: len(table) for key, table in tables.items()},
        "headline": {
            "palearctic_broad_entry_smaller_island_scenarios": int(pal_small.sum()),
            "neotropical_reproductive_assurance_smaller_island_strata": int(
                neo_small.sum()
            ),
            "north_tropical_accessibility_difference_supported_strata": int(
                access_supported.sum()
            ),
            "broad_regime_entry_difference_supported_scenarios": int(
                supported_by_layer.get("analysis_regime", 0)
            ),
            "formal_realm_entry_difference_supported_scenarios": int(
                supported_by_layer.get("biogeographic_realm", 0)
            ),
        },
        "claim_ceiling": lock["claim_ceiling"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lock-path", type=Path, required=True)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--report-path", type=Path, required=True)
    args = parser.parse_args()
    report = validate_submission_freeze(args.lock_path, args.artifact_root)
    args.report_path.parent.mkdir(parents=True, exist_ok=True)
    args.report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
