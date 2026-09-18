"""Freeze the secondary PolLimCrop H4 transportability result without reinterpretation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

RESULT_CONTRACT = "chapter1_h4_pollimcrop_response_mapping_v1"
PREFLIGHT_LOCK_CONTRACT = "chapter1_h4_pollimcrop_transportability_preflight_result_lock_v1"
HYPOTHESES = (
    "H4a_reproductive_assurance",
    "H4b_accessibility_generalization",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise typer.BadParameter(f"invalid JSON object: {path}")
    return value


def _load_yaml(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise typer.BadParameter(f"invalid YAML object: {path}")
    return value


def build_result_lock(
    result: dict[str, Any],
    preflight_lock: dict[str, Any],
    *,
    result_json_path: Path,
    preflight_lock_path: Path,
    response_mapping_path: Path,
    workflow_run_id: int,
    workflow_head_sha: str,
    artifact_id: int,
    artifact_digest: str,
) -> dict[str, Any]:
    if result.get("contract") != RESULT_CONTRACT:
        raise typer.BadParameter("unexpected PolLimCrop result contract")
    if preflight_lock.get("lock_contract") != PREFLIGHT_LOCK_CONTRACT:
        raise typer.BadParameter("unexpected PolLimCrop preflight result lock")
    decision = preflight_lock.get("decision", {})
    if not bool(decision.get("outcome_extraction_authorized")):
        raise typer.BadParameter(
            "preflight lock did not authorize PolLimCrop outcome extraction"
        )

    mapping = _load_yaml(response_mapping_path)
    if mapping.get("contract") != RESULT_CONTRACT:
        raise typer.BadParameter("unexpected frozen response mapping")
    frozen_mapping_sha = str(
        preflight_lock.get("frozen_inputs", {}).get("response_mapping_sha256", "")
    )
    observed_mapping_sha = _sha256(response_mapping_path)
    if frozen_mapping_sha != observed_mapping_sha:
        raise typer.BadParameter("response mapping changed after support lock")

    admitted = [str(x) for x in decision.get("admitted_hypotheses", [])]
    failed = [str(x) for x in decision.get("support_failed_hypotheses", [])]
    if not admitted:
        raise typer.BadParameter("no admitted PolLimCrop hypothesis in support lock")

    primary = result.get("primary_results", {})
    for name in HYPOTHESES:
        if name not in primary:
            raise typer.BadParameter(f"result missing hypothesis: {name}")
    for name in failed:
        item = primary[name]
        if bool(item.get("outcomes_used")):
            raise typer.BadParameter(
                f"support-failed hypothesis used outcomes: {name}"
            )
        if item.get("reason") != "frozen_preflight_support_gate_failed":
            raise typer.BadParameter(
                f"support-failed hypothesis has unexpected result state: {name}"
            )

    supported = [str(x) for x in result.get("supported_primary_hypotheses", [])]
    if any(name not in admitted for name in supported):
        raise typer.BadParameter(
            "supported result includes hypothesis not admitted by support lock"
        )
    if int(result.get("n_supported_primary_hypotheses", -1)) != len(supported):
        raise typer.BadParameter("supported-hypothesis count mismatch")

    classifications: dict[str, str] = {}
    for name in HYPOTHESES:
        if name not in admitted:
            classifications[name] = "not_evaluable_support_gate_failed"
            continue
        item = primary[name]
        if not bool(item.get("evaluable")):
            classifications[name] = "admitted_but_model_not_evaluable"
        elif bool(item.get("supported")):
            classifications[name] = "transportability_supported"
        else:
            classifications[name] = "transportability_not_supported"

    return {
        "contract": "chapter1_h4_pollimcrop_transportability_result_lock_v1",
        "date": "2026-09-18",
        "status": "completed_secondary_external_domain_transportability",
        "inferential_role": "secondary_external_domain_transportability",
        "workflow": {
            "run_id": int(workflow_run_id),
            "head_sha": str(workflow_head_sha),
            "artifact_id": int(artifact_id),
            "artifact_digest": str(artifact_digest),
        },
        "frozen_inputs": {
            "preflight_result_lock": str(preflight_lock_path),
            "preflight_result_lock_sha256": _sha256(preflight_lock_path),
            "response_mapping": str(response_mapping_path),
            "response_mapping_sha256": observed_mapping_sha,
            "result_json_sha256": _sha256(result_json_path),
        },
        "support_decision": {
            "admitted_hypotheses": admitted,
            "support_failed_hypotheses": failed,
        },
        "classification": classifications,
        "primary_results": primary,
        "supported_primary_hypotheses": supported,
        "n_supported_primary_hypotheses": len(supported),
        "claim_ceiling": result.get("claim_ceiling", {}),
        "boundary": {
            "wild_temporal_replication_repaired": False,
            "existing_v13_H4_reclassified_confirmatory_for_wild_flora": False,
            "historical_selection_identified": False,
            "mediation_identified": False,
        },
    }


@app.command("freeze")
def freeze(
    result_json: Path = typer.Option(..., exists=True, dir_okay=False),
    preflight_lock_path: Path = typer.Option(..., exists=True, dir_okay=False),
    response_mapping_path: Path = typer.Option(..., exists=True, dir_okay=False),
    workflow_run_id: int = typer.Option(..., min=1),
    workflow_head_sha: str = typer.Option(...),
    artifact_id: int = typer.Option(..., min=1),
    artifact_digest: str = typer.Option(...),
    output_json: Path = typer.Option(...),
) -> None:
    result = _load_json(result_json)
    preflight_lock = _load_json(preflight_lock_path)
    lock = build_result_lock(
        result,
        preflight_lock,
        result_json_path=result_json,
        preflight_lock_path=preflight_lock_path,
        response_mapping_path=response_mapping_path,
        workflow_run_id=workflow_run_id,
        workflow_head_sha=workflow_head_sha,
        artifact_id=artifact_id,
        artifact_digest=artifact_digest,
    )
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(lock, indent=2) + "\n", encoding="utf-8")
    typer.echo(json.dumps(lock, indent=2))


if __name__ == "__main__":
    app()
