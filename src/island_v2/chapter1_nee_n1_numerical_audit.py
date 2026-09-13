"""Monotone numerical-validity audit for the frozen Chapter 1 N1 fit.

This module does not define an alternative estimator. It replays only the two already-
frozen primary N1 fits and asks whether their numerical outputs are valid enough to be
considered by the original gate. A failed numerical audit can only block N1; it cannot
turn an original N1 failure into a pass or modify the frozen model after outcomes are
opened.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_nee_n1_model import (
    build_primary_frame,
    design_matrix,
    fit_clustered_logit,
    load_config as load_n1_config,
)
from island_v2.chapter1_nee_n1_run import prepare_canonical_covariates

app = typer.Typer(add_completion=False, no_args_is_help=True)


def load_audit_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("N1 numerical audit config must be a mapping")
    if payload.get("contract") != "chapter1_nee_n1_numerical_validity_audit_v1":
        raise typer.BadParameter("unexpected N1 numerical audit contract")
    if payload.get("monotonicity", {}).get("can_only_block_original_N1_gate") is not True:
        raise typer.BadParameter("numerical audit must remain monotone and block-only")
    return payload


def fit_diagnostic(fit: dict[str, Any]) -> dict[str, Any]:
    beta = np.asarray(fit.get("beta", []), dtype=float)
    covariance = np.asarray(fit.get("covariance", []), dtype=float)
    standard_errors = np.asarray(fit.get("standard_errors", []), dtype=float)
    ll = float(fit.get("log_likelihood", float("nan")))
    finite_beta = bool(beta.size > 0 and np.isfinite(beta).all())
    finite_covariance = bool(covariance.size > 0 and np.isfinite(covariance).all())
    finite_se = bool(standard_errors.size > 0 and np.isfinite(standard_errors).all())
    finite_ll = bool(math.isfinite(ll))
    converged = bool(fit.get("converged", False))
    valid = converged and finite_beta and finite_covariance and finite_se and finite_ll
    return {
        "converged": converged,
        "iterations": int(fit.get("iterations", 0)),
        "n_rows": int(fit.get("n_rows", 0)),
        "n_clusters": int(fit.get("n_clusters", 0)),
        "max_absolute_coefficient": float(np.max(np.abs(beta))) if beta.size else None,
        "finite_coefficients": finite_beta,
        "finite_cluster_robust_covariance": finite_covariance,
        "finite_cluster_robust_standard_errors": finite_se,
        "finite_log_likelihood": finite_ll,
        "log_likelihood": ll if finite_ll else None,
        "numerically_valid": valid,
    }


def audit_n1_numerical_validity(
    projected: pd.DataFrame,
    qualification: pd.DataFrame,
    canonical_covariates_raw: pd.DataFrame,
    n1_config: dict[str, Any],
    audit_config: dict[str, Any],
    *,
    expected_islands: int = 8265,
) -> dict[str, Any]:
    covariates = prepare_canonical_covariates(
        canonical_covariates_raw,
        expected_islands=expected_islands,
    )
    frame, metadata = build_primary_frame(projected, qualification, covariates, n1_config)
    if metadata.get("status") != "qualified":
        return {
            "contract": audit_config["contract"],
            "status": "N1_not_evaluable_for_original_gate",
            "numerical_validity_pass": False,
            "original_gate_may_be_considered": False,
            "audit_itself_can_promote_N1": False,
            "failure_action": "stop_before_N2_and_keep_frozen_Chapter1",
            "claim_ceiling": "N1_not_promoted",
            "primary_support": metadata,
            "common_slope_fit": None,
            "heterogeneous_slope_fit": None,
        }

    channels = list(metadata["confirmatory_channels"])
    reference = str(metadata["reference_channel"])
    X0, names0, _ = design_matrix(frame, channels, reference, heterogeneous=False)
    X1, names1, _ = design_matrix(frame, channels, reference, heterogeneous=True)
    common = fit_clustered_logit(frame, X0, names0)
    heterogeneous = fit_clustered_logit(frame, X1, names1)
    common_diag = fit_diagnostic(common)
    heterogeneous_diag = fit_diagnostic(heterogeneous)
    valid = bool(common_diag["numerically_valid"] and heterogeneous_diag["numerically_valid"])

    if valid:
        status = str(audit_config["pass"]["status"])
        action = str(audit_config["pass"]["action"])
        failure_action = None
        claim_ceiling = None
    else:
        status = str(audit_config["failure"]["status"])
        action = None
        failure_action = str(audit_config["failure"]["failure_action"])
        claim_ceiling = str(audit_config["failure"]["claim_ceiling"])

    return {
        "contract": audit_config["contract"],
        "status": status,
        "numerical_validity_pass": valid,
        "original_gate_may_be_considered": valid,
        "audit_itself_can_promote_N1": False,
        "action": action,
        "failure_action": failure_action,
        "claim_ceiling": claim_ceiling,
        "inferential_p_values_promotable": valid,
        "model_comparison_promotable": valid,
        "primary_support": {
            "confirmatory_channels": channels,
            "reference_channel": reference,
            "n_rows": int(len(frame)),
            "n_islands": int(frame["island_id"].nunique()),
            "n_spatial_blocks": int(frame["spatial_block"].astype(str).nunique()),
            "n_source_regions_raw": int(frame["source_region_id"].astype(str).nunique()),
            "n_source_regions_model": int(frame["source_region_model"].astype(str).nunique()),
        },
        "common_slope_fit": common_diag,
        "heterogeneous_slope_fit": heterogeneous_diag,
        "prohibited_rescue": list(audit_config["prohibited"]),
    }


def write_audit(result: dict[str, Any], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "N1_numerical_validity_audit.json").write_text(
        json.dumps(result, indent=2) + "\n",
        encoding="utf-8",
    )


@app.command("run")
def run_command(
    projected_csv: Path = typer.Option(..., exists=True),
    qualification_csv: Path = typer.Option(..., exists=True),
    canonical_covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    n1_config_path: Path = typer.Option(Path("config/chapter1_nee_n1_model.yml"), exists=True),
    audit_config_path: Path = typer.Option(
        Path("config/chapter1_nee_n1_numerical_audit.yml"), exists=True
    ),
    expected_islands: int = typer.Option(8265, min=1),
) -> None:
    n1_config = load_n1_config(n1_config_path)
    audit_config = load_audit_config(audit_config_path)
    result = audit_n1_numerical_validity(
        pd.read_csv(projected_csv, dtype=str).fillna(""),
        pd.read_csv(qualification_csv, dtype=str).fillna(""),
        pd.read_csv(canonical_covariates_csv),
        n1_config,
        audit_config,
        expected_islands=expected_islands,
    )
    write_audit(result, output_dir)
    typer.echo(json.dumps(result, indent=2))


if __name__ == "__main__":
    app()
