"""Audit the one non-converged corrected H1 Direct-only component.

This is a diagnostic only. It does not overwrite the corrected submission tables.
It reruns the all-observed northern-high-latitude Direct-only beta-binomial fits
from exact frozen inputs, retries each outcome from multiple dispersion starts,
and reports both the full seven-response omnibus and a six-response sensitivity
that excludes the originally non-converged shallow/open-tube component.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.optimize import minimize

from island_v2 import chapter1_all_data_probability as h1


def _fit_multistart(
    y: np.ndarray,
    n: np.ndarray,
    design: np.ndarray,
    names: list[str],
    *,
    max_iter: int,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    y = np.asarray(y, dtype=float)
    n = np.asarray(n, dtype=float)
    design = np.asarray(design, dtype=float)
    p = design.shape[1]

    def value_grad(theta: np.ndarray) -> tuple[float, np.ndarray]:
        ll, score = h1._beta_binomial_value_score(theta, y, n, design)
        return -float(np.sum(ll)), -np.sum(score, axis=0)

    attempts: list[dict[str, object]] = []
    candidates: list[tuple[object, float, float]] = []
    for kappa_start in (5.0, 20.0, 50.0, 200.0, 1000.0, 10000.0):
        start = np.zeros(p + 1, dtype=float)
        start[p] = math.log(kappa_start)
        result = minimize(
            lambda t: value_grad(t)[0],
            start,
            jac=lambda t: value_grad(t)[1],
            method="L-BFGS-B",
            bounds=[(None, None)] * p + [(-6.0, 12.0)],
            options={"maxiter": max(5000, int(max_iter)), "ftol": 1e-11, "gtol": 1e-6},
        )
        theta = np.asarray(result.x, dtype=float)
        objective, grad = value_grad(theta)
        grad_inf = float(np.max(np.abs(grad)))
        loglik = -float(objective)
        attempts.append(
            {
                "kappa_start": kappa_start,
                "success": bool(result.success),
                "message": str(result.message),
                "nit": int(getattr(result, "nit", -1)),
                "log_likelihood": loglik,
                "gradient_inf_norm": grad_inf,
                "kappa": float(np.exp(np.clip(theta[-1], -6.0, 12.0))),
            }
        )
        if np.isfinite(loglik) and np.all(np.isfinite(theta)):
            candidates.append((result, loglik, grad_inf))

    if not candidates:
        raise RuntimeError("all multistart attempts were non-finite")

    successful = [x for x in candidates if bool(x[0].success)]
    pool = successful if successful else candidates
    result, loglik, grad_inf = max(pool, key=lambda x: x[1])
    theta = np.asarray(result.x, dtype=float)
    dim = len(theta)
    hessian = np.zeros((dim, dim), dtype=float)
    for j in range(dim):
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
    _, score = h1._beta_binomial_value_score(theta, y, n, design)
    fit = {
        "success": bool(result.success),
        "message": str(result.message),
        "theta": theta,
        "bread": bread,
        "score": score,
        "names": [*names, "log_kappa"],
        "log_likelihood": loglik,
        "gradient_inf_norm": grad_inf,
        "kappa": float(np.exp(np.clip(theta[-1], -6.0, 12.0))),
    }
    return fit, attempts


def _omnibus(
    fits: list[dict[str, object]],
    clusters: list[np.ndarray],
    slope_indices: list[int],
) -> dict[str, float | int]:
    covariance, _, theta = h1._assemble_cluster_covariance(fits, clusters)
    slopes = theta[slope_indices]
    slope_cov = covariance[np.ix_(slope_indices, slope_indices)]
    rank = int(np.linalg.matrix_rank(slope_cov))
    statistic = float(slopes @ np.linalg.pinv(slope_cov) @ slopes)
    return {
        "joint_wald_chisq": statistic,
        "joint_df": rank,
        "p_value": h1._chi_square_sf_integer_df(statistic, rank),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--geometry", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    config = yaml.safe_load(
        (args.sources / "island--chapter1_v14_all_data_probability.yml").read_text(
            encoding="utf-8"
        )
    )
    counts = pd.read_csv(
        args.sources / "v14-artifact/direct/h1/all_data_probability_counts.csv.gz"
    )
    covariates = pd.read_csv(args.geometry / "corrected_geography_covariates.csv")
    data = h1._prepare(counts, covariates, config)

    stratum = "all_observed"
    context_value = "northern_high_latitude"
    geography = str(config["geography_column"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    threshold = int(config["minimum_islands_per_outcome"])
    work = data.loc[
        data["stratum"].eq(stratum) & data[context].eq(context_value)
    ].copy()
    support = work.groupby("outcome")["island_id"].nunique()
    retained = [
        str(x)
        for x in config["model_outcomes"]
        if int(support.get(str(x), 0)) >= threshold
    ]

    fits: list[dict[str, object]] = []
    clusters: list[np.ndarray] = []
    slope_indices: list[int] = []
    audit_rows: list[dict[str, object]] = []
    all_attempts: dict[str, object] = {}
    offset = 0

    frozen_slopes = pd.read_csv(
        args.geometry / "direct/beta_binomial_within_slopes.csv"
    )
    frozen_slopes = frozen_slopes.loc[
        frozen_slopes["stratum"].eq(stratum)
        & frozen_slopes["context"].eq(context_value)
    ].set_index("outcome")

    for outcome in retained:
        part = work.loc[work["outcome"].eq(outcome)].copy()
        columns = [np.ones(len(part), dtype=float)]
        names = [f"{outcome}:intercept"]
        for predictor in baseline:
            columns.append(h1._standardize(part[predictor]))
            names.append(f"{outcome}:z_{predictor}")
        columns.append(h1._standardize(part[geography]))
        slope_name = f"{outcome}:z_{geography}"
        names.append(slope_name)
        design = np.column_stack(columns)

        default = h1._fit_single_beta_binomial(
            part["successes"].to_numpy(float),
            part["trials"].to_numpy(float),
            design,
            names,
            max_iter=int(config.get("max_iter", 800)),
        )
        retry, attempts = _fit_multistart(
            part["successes"].to_numpy(float),
            part["trials"].to_numpy(float),
            design,
            names,
            max_iter=int(config.get("max_iter", 800)),
        )
        fits.append(retry)
        clusters.append(part[cluster].to_numpy(str))
        local_slope = names.index(slope_name)
        slope_indices.append(offset + local_slope)
        frozen = frozen_slopes.loc[outcome]
        retry_estimate = float(np.asarray(retry["theta"])[local_slope])
        default_estimate = float(np.asarray(default["theta"])[local_slope])
        audit_rows.append(
            {
                "outcome": outcome,
                "n_islands": int(part["island_id"].nunique()),
                "frozen_optimizer_success": bool(frozen["optimizer_success"]),
                "default_success": bool(default["success"]),
                "default_message": str(default["message"]),
                "retry_success": bool(retry["success"]),
                "retry_message": str(retry["message"]),
                "frozen_slope": float(frozen["geography_slope_log_odds"]),
                "default_slope": default_estimate,
                "retry_slope": retry_estimate,
                "retry_minus_frozen": retry_estimate
                - float(frozen["geography_slope_log_odds"]),
                "default_log_likelihood": float(default["log_likelihood"]),
                "retry_log_likelihood": float(retry["log_likelihood"]),
                "retry_gradient_inf_norm": float(retry["gradient_inf_norm"]),
                "retry_kappa": float(retry["kappa"]),
            }
        )
        all_attempts[outcome] = attempts
        offset += len(retry["names"])

    full = _omnibus(fits, clusters, slope_indices)
    drop_name = "shallow_open_tube"
    keep_positions = [i for i, name in enumerate(retained) if name != drop_name]
    fits6 = [fits[i] for i in keep_positions]
    clusters6 = [clusters[i] for i in keep_positions]
    slope_indices6: list[int] = []
    offset6 = 0
    for i in keep_positions:
        fit = fits[i]
        outcome = retained[i]
        slope_name = f"{outcome}:z_{geography}"
        local = list(fit["names"]).index(slope_name)
        slope_indices6.append(offset6 + local)
        offset6 += len(fit["names"])
    six = _omnibus(fits6, clusters6, slope_indices6)

    frozen_omnibus = pd.read_csv(
        args.geometry / "direct/beta_binomial_within_omnibus.csv"
    )
    frozen_row = frozen_omnibus.loc[
        frozen_omnibus["stratum"].eq(stratum)
        & frozen_omnibus["context"].eq(context_value)
    ].iloc[0]

    summary = {
        "status": "pass"
        if all(bool(f["success"]) for f in fits)
        and float(full["p_value"]) < 0.05
        and float(six["p_value"]) < 0.05
        else "review",
        "scope": "direct_only",
        "stratum": stratum,
        "context": context_value,
        "original_nonconverged_outcome": drop_name,
        "frozen": {
            "p_value": float(frozen_row["p_value"]),
            "q_value": float(frozen_row["q_value"]),
            "all_optimizers_converged": bool(frozen_row["all_optimizers_converged"]),
        },
        "multistart_seven_response": full,
        "drop_nonconverged_six_response_sensitivity": six,
        "all_multistart_fits_successful": all(bool(f["success"]) for f in fits),
        "max_abs_retry_minus_frozen_slope": max(
            abs(float(row["retry_minus_frozen"])) for row in audit_rows
        ),
        "claim_boundary": (
            "Diagnostic only. The six-response result is a sensitivity, not a "
            "replacement predeclared H1 estimand."
        ),
    }

    pd.DataFrame(audit_rows).to_csv(args.output / "fit_audit.csv", index=False)
    (args.output / "attempts.json").write_text(
        json.dumps(all_attempts, indent=2), encoding="utf-8"
    )
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
