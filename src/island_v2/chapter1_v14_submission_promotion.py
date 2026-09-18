"""Fail-closed promotion gate for the Chapter 1 v14 candidate result surface."""

from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
from typing import Any

RESULT_CONTRACT = "chapter1_v14_reordered_hypotheses_result_v1"
PROMOTED_CONTRACT = "chapter1_v14_canonical_result_lock_v1"
TOL = 1e-9


def _load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def _close(observed: object, expected: object, *, label: str) -> None:
    a = float(observed)
    b = float(expected)
    if not math.isfinite(a) or not math.isfinite(b):
        raise ValueError(f"{label}: non-finite value")
    if not math.isclose(a, b, rel_tol=TOL, abs_tol=TOL):
        raise ValueError(f"{label}: observed {a} != expected {b}")


def _rows_by(rows: list[dict[str, Any]], *keys: str) -> dict[tuple[str, ...], dict[str, Any]]:
    out: dict[tuple[str, ...], dict[str, Any]] = {}
    for row in rows:
        key = tuple(str(row[k]) for k in keys)
        if key in out:
            raise ValueError(f"duplicate result row for {key}")
        out[key] = row
    return out


def _validate_h1(
    result: dict[str, Any],
    expected: dict[str, Any],
) -> None:
    scope_map = {"all": "all_analysis", "direct": "direct_only"}
    for observed_scope, expected_scope in scope_map.items():
        payload = result["scopes"][observed_scope]
        omnibus = _rows_by(payload["H1_all_observed_omnibus"], "context")
        orientation = _rows_by(payload["H1_equal_domain_orientation"], "context")
        for context, target in expected["H1"][expected_scope].items():
            o = omnibus[(context,)]
            d = orientation[(context,)]
            _close(o["p_value"], target["joint_p"], label=f"H1 {observed_scope} {context} p")
            _close(o["q_value"], target["joint_q"], label=f"H1 {observed_scope} {context} q")
            if not bool(o["vector_supported"]):
                raise ValueError(f"H1 {observed_scope} {context}: vector no longer supported")
            _close(
                d["equal_domain_orientation"],
                target["equal_domain_orientation"],
                label=f"H1 {observed_scope} {context} orientation",
            )


def _validate_h2(
    result: dict[str, Any],
    expected: dict[str, Any],
) -> None:
    scope_map = {"all": "all_analysis", "direct": "direct_only"}
    response_map = {
        "selfing_core": "selfing_core",
        "generalized_accessible": "generalized_accessible_given_selfing",
        "plain_colour": "plain_colour_given_selfing",
    }
    for observed_scope, expected_scope in scope_map.items():
        rows = _rows_by(result["scopes"][observed_scope]["H2_models"], "context", "response")
        for context, target_context in expected["H2"][expected_scope].items():
            for response, target_key in response_map.items():
                row = rows[(context, response)]
                target = target_context[target_key]
                _close(
                    row["distance_estimate"],
                    target["beta"],
                    label=f"H2 {observed_scope} {context} {response} beta",
                )
                _close(
                    row["distance_p"],
                    target["p"],
                    label=f"H2 {observed_scope} {context} {response} p",
                )
                if "q" in target:
                    _close(
                        row["primary_H2b_q"],
                        target["q"],
                        label=f"H2 {observed_scope} {context} {response} q",
                    )


def _validate_bridge_rows(
    observed_rows: list[dict[str, Any]],
    expected: dict[str, Any],
    *,
    label: str,
) -> None:
    rows = _rows_by(observed_rows, "family", "analysis")
    for family, payload in expected["families"].items():
        for analysis in ("primary", "supplemental_only", "no_zero_constant"):
            target = payload[analysis]
            row = rows[(family, analysis)]
            _close(row["estimate"], target["estimate"], label=f"{label} {family} {analysis} beta")
            _close(row["se"], target["se"], label=f"{label} {family} {analysis} se")
            if int(row["n_species"]) != int(target["n_species"]):
                raise ValueError(f"{label} {family} {analysis}: n_species changed")
            if int(row["n_publications"]) != int(target["n_publications"]):
                raise ValueError(f"{label} {family} {analysis}: n_publications changed")


def validate_reproduction(
    result: dict[str, Any],
    preflight: dict[str, Any],
    exact_bridge: dict[str, Any],
    reconstruction: dict[str, Any],
) -> dict[str, Any]:
    if result.get("contract") != RESULT_CONTRACT:
        raise ValueError("unexpected v14 workflow result contract")
    if preflight.get("canonical") is not False:
        raise ValueError("preflight lock must remain non-canonical before promotion")
    _validate_h1(result, preflight)
    _validate_h2(result, preflight)
    _validate_bridge_rows(
        result["H4_exact_H2_score_bridge_results"],
        exact_bridge,
        label="H4 exact",
    )
    _validate_bridge_rows(
        result["H4_atomic_reconstruction_results"],
        reconstruction,
        label="H4 reconstruction",
    )
    return {
        "verified": True,
        "result_contract": RESULT_CONTRACT,
        "H1_reproduced": True,
        "H2_reproduced": True,
        "H4_exact_H2_score_bridge_reproduced": True,
        "H4_atomic_reconstruction_reproduced": True,
    }


def promote(
    result: dict[str, Any],
    preflight: dict[str, Any],
    verification: dict[str, Any],
    *,
    ci_run_id: int,
    ci_artifact_id: int,
    ci_artifact_digest: str,
    head_sha: str,
) -> dict[str, Any]:
    if verification.get("verified") is not True:
        raise ValueError("cannot promote an unverified result")
    if not ci_artifact_digest.startswith("sha256:"):
        raise ValueError("CI artifact digest must be a sha256 digest")
    if len(head_sha) != 40:
        raise ValueError("head SHA must be a full 40-character Git SHA")

    out = copy.deepcopy(preflight)
    out["contract"] = PROMOTED_CONTRACT
    out["status"] = "canonical_v14_reproduced"
    out["canonical"] = True
    out["promotion_provenance"] = {
        "head_sha": head_sha,
        "ci_run_id": int(ci_run_id),
        "ci_artifact_id": int(ci_artifact_id),
        "ci_artifact_digest": ci_artifact_digest,
        "workflow_result_contract": str(result["contract"]),
        "verification": verification,
    }
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)

    verify = sub.add_parser("verify")
    verify.add_argument("--result-summary", type=Path, required=True)
    verify.add_argument("--preflight-lock", type=Path, required=True)
    verify.add_argument("--exact-bridge-lock", type=Path, required=True)
    verify.add_argument("--reconstruction-lock", type=Path, required=True)
    verify.add_argument("--report", type=Path, required=True)

    promote_parser = sub.add_parser("promote")
    promote_parser.add_argument("--result-summary", type=Path, required=True)
    promote_parser.add_argument("--preflight-lock", type=Path, required=True)
    promote_parser.add_argument("--exact-bridge-lock", type=Path, required=True)
    promote_parser.add_argument("--reconstruction-lock", type=Path, required=True)
    promote_parser.add_argument("--ci-run-id", type=int, required=True)
    promote_parser.add_argument("--ci-artifact-id", type=int, required=True)
    promote_parser.add_argument("--ci-artifact-digest", required=True)
    promote_parser.add_argument("--head-sha", required=True)
    promote_parser.add_argument("--output", type=Path, required=True)

    args = parser.parse_args()
    result = _load(args.result_summary)
    preflight = _load(args.preflight_lock)
    exact = _load(args.exact_bridge_lock)
    reconstruction = _load(args.reconstruction_lock)
    verification = validate_reproduction(result, preflight, exact, reconstruction)

    if args.command == "verify":
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(verification, indent=2) + "\n", encoding="utf-8")
        return

    promoted = promote(
        result,
        preflight,
        verification,
        ci_run_id=args.ci_run_id,
        ci_artifact_id=args.ci_artifact_id,
        ci_artifact_digest=args.ci_artifact_digest,
        head_sha=args.head_sha,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(promoted, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
