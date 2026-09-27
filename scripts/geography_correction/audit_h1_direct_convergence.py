from __future__ import annotations

import json
import math
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import typer
import yaml
from scipy.optimize import minimize

from island_v2 import chapter1_all_data_probability as h1

app = typer.Typer(add_completion=False, no_args_is_help=True)

TARGET_STRATUM = "all_observed"
TARGET_CONTEXT = "northern_high_latitude"
TARGET_OUTCOME = "shallow_open_tube"


def _retry_fit(
    y: np.ndarray,
    n: np.ndarray,
    design: np.ndarray,
    names: list[str],
    *,
    initial_theta: np.ndarray,
) -> dict[str, object]:
    p = design.shape[1]

    def value_grad(theta: np.ndarray) -> tuple[float, np.ndarray]:
        ll, score = h1._beta_binomial_value_score(theta, y, n, design)
        return -float(np.sum(ll)), -np.sum(score, axis=0)

    result = minimize(
        lambda t: value_grad(t)[0],
        np.asarray(initial_theta, dtype=float),
        jac=lambda t: value_grad(t)[1],
        method="L-BFGS-B",
        bounds=[(None, None)] * p + [(-6.0, 12.0)],
        options={
            "maxiter": 10000,
            "maxls": 200,
            "ftol": 1e-12,
            "gtol": 1e-7,
        },
    )
    theta = np.asarray(result.x, dtype=float)
    dimension = len(theta)
    hessian = np.zeros((dimension, dimension), dtype=float)
    for j in range(dimension):
        step = 1e-5 * (1.0 + abs(theta[j]))
        plus = theta.copy()
        minus = theta.copy()
        plus[j] += step
        minus[j] -= step
        grad_plus = value_grad(plus)[1]
        grad_minus = value_grad(minus)[1]
        hessian[:, j] = (grad_plus - grad_minus) / (2.0 * step)
    hessian = (hessian + hessian.T) / 2.0
    bread = np.linalg.pinv(hessian, rcond=1e-10)
    loglik, score = h1._beta_binomial_value_score(theta, y, n, design)
    grad = value_grad(theta)[1]
    return {
        "success": bool(result.success),
        "message": str(result.message),
        "theta": theta,
        "bread": bread,
        "score": score,
        "names": [*names, "log_kappa"],
        "log_likelihood": float(np.sum(loglik)),
        "kappa": float(np.exp(np.clip(theta[-1], -6.0, 12.0))),
        "gradient_norm": float(np.linalg.norm(grad)),
        "n_iterations": int(getattr(result, "nit", -1)),
        "n_function_evals": int(getattr(result, "nfev", -1)),
    }


def _target_design(
    data: pd.DataFrame,
    config: dict[str, object],
) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    context_col = str(config["context_column"])
    geography = str(config["geography_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    work = data.loc[
        data["stratum"].eq(TARGET_STRATUM)
        & data[context_col].eq(TARGET_CONTEXT)
        & data["outcome"].eq(TARGET_OUTCOME)
    ].copy()
    columns = [np.ones(len(work), dtype=float)]
    names = [f"{TARGET_OUTCOME}:intercept"]
    for predictor in baseline:
        columns.append(h1._standardize(work[predictor]))
        names.append(f"{TARGET_OUTCOME}:z_{predictor}")
    columns.append(h1._standardize(work[geography]))
    names.append(f"{TARGET_OUTCOME}:z_{geography}")
    return work, np.column_stack(columns), names


def _extract_target(slopes: pd.DataFrame) -> dict[str, object]:
    row = slopes.loc[
        slopes["stratum"].eq(TARGET_STRATUM)
        & slopes["context"].eq(TARGET_CONTEXT)
        & slopes["outcome"].eq(TARGET_OUTCOME)
    ].iloc[0]
    return {
        "estimate": float(row["geography_slope_log_odds"]),
        "se": float(row["cluster_robust_se"]),
        "p_value": float(row["p_value"]),
        "optimizer_success": bool(row["optimizer_success"]),
        "n_islands": int(row["n_islands"]),
        "n_species_trials": int(row["n_species_trials"]),
        "kappa": float(row["kappa"]),
    }


def _extract_omnibus(omnibus: pd.DataFrame) -> dict[str, object]:
    row = omnibus.loc[
        omnibus["stratum"].eq(TARGET_STRATUM)
        & omnibus["context"].eq(TARGET_CONTEXT)
    ].iloc[0]
    return {
        "n_retained_outcomes": int(row["n_retained_outcomes"]),
        "retained_outcomes": str(row["retained_outcomes"]),
        "joint_wald_chisq": float(row["joint_wald_chisq"]),
        "joint_df": int(row["joint_df"]),
        "p_value": float(row["p_value"]),
        "q_value": float(row["q_value"]),
        "all_optimizers_converged": bool(row["all_optimizers_converged"]),
        "vector_supported": bool(row["vector_supported"]),
    }


@app.command()
def main(
    counts_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    covariates_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    frozen_slopes_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    frozen_omnibus_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    output: Path = typer.Option(...),
) -> None:
    counts = pd.read_csv(counts_csv)
    covariates = pd.read_csv(covariates_csv)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))

    frozen_slopes = pd.read_csv(frozen_slopes_csv)
    frozen_omnibus = pd.read_csv(frozen_omnibus_csv)
    frozen_target = _extract_target(frozen_slopes)
    frozen_joint = _extract_omnibus(frozen_omnibus)

    default_tables = h1.run_probability_analysis(counts, covariates, config)
    default_target = _extract_target(default_tables[0])
    default_joint = _extract_omnibus(default_tables[2])

    prepared = h1._prepare(counts, covariates, config)
    target_work, target_design, target_names = _target_design(prepared, config)
    default_fit = h1._fit_single_beta_binomial(
        target_work["successes"].to_numpy(float),
        target_work["trials"].to_numpy(float),
        target_design,
        target_names,
        max_iter=int(config.get("max_iter", 1000)),
    )
    retry_fit = _retry_fit(
        target_work["successes"].to_numpy(float),
        target_work["trials"].to_numpy(float),
        target_design,
        target_names,
        initial_theta=np.asarray(default_fit["theta"], dtype=float),
    )
    geo_name = f"{TARGET_OUTCOME}:z_{config['geography_column']}"
    geo_idx = retry_fit["names"].index(geo_name)
    retry_estimate = float(np.asarray(retry_fit["theta"])[geo_idx])

    # Re-run the seven-response vector with a fail-then-retry fitter.
    original_fitter = h1._fit_single_beta_binomial

    def robust_fitter(
        y: np.ndarray,
        n: np.ndarray,
        design: np.ndarray,
        names: list[str],
        *,
        max_iter: int,
    ) -> dict[str, object]:
        first = original_fitter(y, n, design, names, max_iter=max_iter)
        if bool(first["success"]):
            first["retry_used"] = False
            return first
        second = _retry_fit(
            y,
            n,
            design,
            names,
            initial_theta=np.asarray(first["theta"], dtype=float),
        )
        second["retry_used"] = True
        second["first_message"] = str(first["message"])
        return second

    h1._fit_single_beta_binomial = robust_fitter
    try:
        robust_tables = h1.run_probability_analysis(counts, covariates, config)
    finally:
        h1._fit_single_beta_binomial = original_fitter

    robust_target = _extract_target(robust_tables[0])
    robust_joint = _extract_omnibus(robust_tables[2])

    # Fully converged sensitivity that removes only the warned component.
    six_config = deepcopy(config)
    six_config["model_outcomes"] = [
        x for x in config["model_outcomes"] if x != TARGET_OUTCOME
    ]
    six_tables = h1.run_probability_analysis(counts, covariates, six_config)
    six_joint = _extract_omnibus(six_tables[2])

    estimate_delta = abs(robust_target["estimate"] - frozen_target["estimate"])
    retry_delta = abs(retry_estimate - frozen_target["estimate"])
    audit_pass = bool(
        retry_fit["success"]
        and robust_target["optimizer_success"]
        and robust_joint["all_optimizers_converged"]
        and robust_joint["q_value"] < 0.05
        and six_joint["all_optimizers_converged"]
        and six_joint["q_value"] < 0.05
        and estimate_delta < 1e-4
        and retry_delta < 1e-4
    )

    report = {
        "contract": "chapter1_h1_direct_northern_high_optimizer_audit_v1",
        "target": {
            "stratum": TARGET_STRATUM,
            "context": TARGET_CONTEXT,
            "outcome": TARGET_OUTCOME,
        },
        "frozen_corrected": {
            "target_fit": frozen_target,
            "joint_vector": frozen_joint,
        },
        "default_replay": {
            "target_fit": default_target,
            "joint_vector": default_joint,
            "target_optimizer_message": str(default_fit["message"]),
        },
        "enhanced_retry": {
            "success": bool(retry_fit["success"]),
            "message": str(retry_fit["message"]),
            "gradient_norm": float(retry_fit["gradient_norm"]),
            "n_iterations": int(retry_fit["n_iterations"]),
            "n_function_evals": int(retry_fit["n_function_evals"]),
            "target_estimate": retry_estimate,
            "absolute_delta_from_frozen": retry_delta,
        },
        "robust_seven_response_replay": {
            "target_fit": robust_target,
            "joint_vector": robust_joint,
            "absolute_target_estimate_delta_from_frozen": estimate_delta,
        },
        "six_response_sensitivity": {
            "excluded_outcome": TARGET_OUTCOME,
            "joint_vector": six_joint,
        },
        "decision": {
            "audit_pass": audit_pass,
            "criterion": (
                "enhanced target retry converges; robust seven-response replay "
                "converges and remains FDR-supported; six-response sensitivity "
                "converges and remains FDR-supported; target estimate drift <1e-4"
            ),
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    typer.echo(json.dumps(report, indent=2))
    if not audit_pass:
        raise typer.Exit(code=2)


if __name__ == "__main__":
    app()
