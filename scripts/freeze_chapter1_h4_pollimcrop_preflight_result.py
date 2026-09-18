"""Freeze an outcome-blind PolLimCrop H4 preflight result.

This utility consumes only the already-created metadata-only PREFLIGHT.json plus the
frozen preflight and response-mapping contracts. It never reads the PolLimCrop dataset
or any pollen-limitation outcome value.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

PREFLIGHT_CONTRACT = "chapter1_h4_pollimcrop_transportability_preflight_v1"
MAPPING_CONTRACT = "chapter1_h4_pollimcrop_response_mapping_v1"
HYPOTHESES = (
    "H4a_reproductive_assurance",
    "H4b_accessibility_generalization",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_yaml(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise typer.BadParameter(f"invalid YAML object: {path}")
    return value


def build_result_lock(
    preflight: dict[str, Any],
    *,
    preflight_config_path: Path,
    response_mapping_path: Path,
    workflow_run_id: int,
    workflow_head_sha: str,
    artifact_id: int,
    artifact_digest: str,
) -> dict[str, Any]:
    if preflight.get("contract") != PREFLIGHT_CONTRACT:
        raise typer.BadParameter("unexpected PolLimCrop preflight result contract")
    if preflight.get("outcomes_read") is not False:
        raise typer.BadParameter("PolLimCrop preflight result must remain outcome-blind")

    config = _load_yaml(preflight_config_path)
    mapping = _load_yaml(response_mapping_path)
    if config.get("contract") != PREFLIGHT_CONTRACT:
        raise typer.BadParameter("unexpected PolLimCrop preflight config")
    if mapping.get("contract") != MAPPING_CONTRACT:
        raise typer.BadParameter("unexpected PolLimCrop response mapping")
    if mapping.get("status") != "frozen_before_row_level_PolLimCrop_outcome_read":
        raise typer.BadParameter("PolLimCrop response mapping is not frozen")

    mapping_result = preflight.get("response_mapping", {})
    observed_mapping_sha = str(mapping_result.get("sha256", ""))
    expected_mapping_sha = _sha256(response_mapping_path)
    if observed_mapping_sha != expected_mapping_sha:
        raise typer.BadParameter("preflight response-mapping SHA-256 mismatch")
    if str(mapping_result.get("response_column", "")) != str(
        mapping["source_method_basis"]["response_column"]
    ):
        raise typer.BadParameter("preflight response-column token mismatch")
    if mapping_result.get("outcome_values_read") is not False:
        raise typer.BadParameter("response mapping preflight must not read outcomes")

    expected_scope = {
        key: int(value)
        for key, value in config["source"]["expected_scope"].items()
    }
    observed_scope = {
        key: int(value)
        for key, value in preflight.get("validated_expected_scope", {}).items()
    }
    if observed_scope != expected_scope:
        raise typer.BadParameter(
            f"PolLimCrop validated scope mismatch: expected {expected_scope}, "
            f"observed {observed_scope}"
        )

    support = preflight.get("support", {})
    missing = [name for name in HYPOTHESES if name not in support]
    if missing:
        raise typer.BadParameter(f"preflight missing support results: {missing}")

    admitted = [
        name for name in HYPOTHESES if bool(support[name].get("evaluable"))
    ]
    rejected = [name for name in HYPOTHESES if name not in admitted]
    authorized = bool(admitted)

    return {
        "contract": PREFLIGHT_CONTRACT,
        "lock_contract": "chapter1_h4_pollimcrop_transportability_preflight_result_lock_v1",
        "date": "2026-09-18",
        "status": (
            "outcome_blind_support_gate_passed"
            if authorized
            else "outcome_blind_support_gate_not_met"
        ),
        "inferential_role": "secondary_external_domain_transportability_support_decision",
        "workflow": {
            "run_id": int(workflow_run_id),
            "head_sha": str(workflow_head_sha),
            "artifact_id": int(artifact_id),
            "artifact_digest": str(artifact_digest),
        },
        "frozen_inputs": {
            "preflight_config": str(preflight_config_path),
            "preflight_config_sha256": _sha256(preflight_config_path),
            "response_mapping": str(response_mapping_path),
            "response_mapping_sha256": expected_mapping_sha,
            "response_column": str(
                mapping["source_method_basis"]["response_column"]
            ),
        },
        "source_format": preflight.get("source_format", {}),
        "figshare": preflight.get("figshare", {}),
        "validated_expected_scope": observed_scope,
        "metadata_columns_read": preflight.get("metadata_columns_read", []),
        "support": support,
        "decision": {
            "admitted_hypotheses": admitted,
            "support_failed_hypotheses": rejected,
            "outcome_extraction_authorized": authorized,
            "outcomes_may_be_read_only_for_admitted_hypotheses": True,
            "threshold_relaxation_authorized": False,
            "wild_temporal_replication_repaired": False,
        },
        "blinding": {
            "preflight_outcomes_read": False,
            "row_level_PL_effectsize_read_before_lock": False,
            "support_thresholds_changed": False,
            "trait_definitions_changed": False,
        },
        "claim_ceiling": config["claim_ceiling"],
    }


@app.command("freeze")
def freeze(
    preflight_json: Path = typer.Option(..., exists=True, dir_okay=False),
    preflight_config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    response_mapping_path: Path = typer.Option(..., exists=True, dir_okay=False),
    workflow_run_id: int = typer.Option(..., min=1),
    workflow_head_sha: str = typer.Option(...),
    artifact_id: int = typer.Option(..., min=1),
    artifact_digest: str = typer.Option(...),
    output_json: Path = typer.Option(...),
) -> None:
    preflight = json.loads(preflight_json.read_text(encoding="utf-8"))
    result = build_result_lock(
        preflight,
        preflight_config_path=preflight_config_path,
        response_mapping_path=response_mapping_path,
        workflow_run_id=workflow_run_id,
        workflow_head_sha=workflow_head_sha,
        artifact_id=artifact_id,
        artifact_digest=artifact_digest,
    )
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    typer.echo(json.dumps(result, indent=2))


if __name__ == "__main__":
    app()
