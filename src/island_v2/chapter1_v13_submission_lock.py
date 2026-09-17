"""Fail-closed validator for the Chapter 1 v13 unified submission lock."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


CONTRACT = "chapter1_v13_unified_island_syndrome_result_lock_v1"
EXPECTED_PARENTS = {
    "all_data": "chapter1_all_data_route_result_lock_v2",
    "v12_two_panel": "chapter1_v12_two_panel_result_lock_v7",
    "v12_glopl": "chapter1_v12_h5_glopl_extension_result_lock_v4",
    "functional_bridge": "chapter1_v13_functional_bridge_result_lock_v1",
}
EXPECTED_ARCHITECTURE = {
    "H1": "global_recurrent_island_syndrome",
    "H2": "global_pollination_constraint",
    "H3": "dual_plant_response_pathways",
    "H4": "posthoc_functional_triangulation",
    "H5": "taxonomic_realization",
}
FALSE_CLAIMS = (
    "historical_pollen_limitation_selected_traits",
    "pollen_limitation_mediates_global_syndrome",
    "GloBI_identifies_pollinator_mechanism",
    "all_observed_equals_native_assembly",
    "taxonomic_attenuation_identifies_cause",
)


def _require_mapping(value: object, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    return value


def validate_v13_lock(lock: dict[str, Any]) -> dict[str, Any]:
    """Validate provenance, architecture, inferential role, and claim ceilings."""

    if lock.get("contract") != CONTRACT:
        raise ValueError("unexpected v13 submission-lock contract")

    parents = _require_mapping(lock.get("parents"), "parents")
    for name, expected_contract in EXPECTED_PARENTS.items():
        parent = _require_mapping(parents.get(name), f"parents.{name}")
        if parent.get("contract") != expected_contract:
            raise ValueError(f"parent contract changed: {name}")

    bridge = _require_mapping(parents.get("functional_bridge"), "parents.functional_bridge")
    artifact_id = bridge.get("artifact_id")
    digest = str(bridge.get("artifact_digest", ""))
    if not artifact_id or not digest.startswith("sha256:"):
        raise ValueError("functional bridge provenance is incomplete")

    architecture = _require_mapping(lock.get("architecture"), "architecture")
    if architecture != EXPECTED_ARCHITECTURE:
        raise ValueError("v13 H1-H5 architecture changed")

    h4 = _require_mapping(lock.get("H4_functional_bridge"), "H4_functional_bridge")
    if h4.get("inferential_role") != "posthoc_functional_triangulation":
        raise ValueError("functional bridge must remain posthoc functional triangulation")
    if h4.get("artifact_id") != artifact_id:
        raise ValueError("functional bridge artifact provenance mismatch")
    if h4.get("historical_selection_identified") is not False:
        raise ValueError("historical selection cannot be identified by H4")

    supplementary = _require_mapping(lock.get("supplementary"), "supplementary")
    if supplementary.get("GloBI_role") != "non_promoted_sampling_sensitive_context_evidence":
        raise ValueError("GloBI must remain non-promoted supplementary evidence")

    ceiling = _require_mapping(lock.get("claim_ceiling"), "claim_ceiling")
    for claim in FALSE_CLAIMS:
        if ceiling.get(claim) is not False:
            if claim.startswith("GloBI"):
                raise ValueError("GloBI mechanism promotion is prohibited")
            if claim.startswith("historical"):
                raise ValueError("historical pollen-limitation selection is not identified")
            raise ValueError(f"claim ceiling was raised: {claim}")

    h5 = lock.get("H5_taxonomic_realization")
    if h5 is not None:
        h5 = _require_mapping(h5, "H5_taxonomic_realization")
        if h5.get("causal_mediation_identified") is not False:
            raise ValueError("taxonomic realization cannot be relabelled causal mediation")

    return {
        "contract": "chapter1_v13_submission_lock_verification_v1",
        "verified": True,
        "functional_bridge_artifact_id": int(artifact_id),
        "functional_bridge_artifact_digest": digest,
        "inferential_role": str(h4["inferential_role"]),
        "claim_ceiling_verified": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lock-path", type=Path, required=True)
    parser.add_argument("--report-path", type=Path, required=True)
    args = parser.parse_args()
    lock = json.loads(args.lock_path.read_text(encoding="utf-8"))
    report = validate_v13_lock(lock)
    args.report_path.parent.mkdir(parents=True, exist_ok=True)
    args.report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
