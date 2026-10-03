"""Reviewer-robust final inference for Chapter 1 H1-H2.

This module does not invent a new response after seeing the October raw-state results.
H1 reuses the seven v14 outcomes whose positive classes were already frozen in
`config/chapter1_v14_all_data_probability.yml`, and reuses the already-frozen
equal-domain orientation (three biological domains receive equal total weight).

The change is inferential:
- H1 is reduced from a high-dimensional direction-free Wald test to a one-degree-of-
  freedom directional projection along the frozen classic-island orientation.
- uncertainty for H1 and the continuous H2 branches is estimated by exact delete-one-
  spatial-cluster jackknife refits, with t reference df = G-1.
- the 97-state three-axis analysis remains a descriptive realization/diagnostic layer.

This is a robustness re-estimation of a frozen biological contrast, not a prospective
confirmation and not evidence for within-lineage evolution or causal mediation.
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
from scipy.stats import t as student_t

from island_v2.chapter1_all_data_probability import (
    _fit_single_beta_binomial,
    build_broad_counts,
)
from island_v2.chapter1_v14_h2_decomposition import build_plant_scores

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _holm(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=values.index, dtype=float)
    ok = p.notna()
    if not ok.any():
        return out
    x = p.loc[ok].to_numpy(float)
    order = np.argsort(x)
    ranked = x[order]
    n = len(ranked)
    adjusted = np.maximum.accumulate(ranked * (n - np.arange(n)))
    adjusted = np.clip(adjusted, 0.0, 1.0)
    restored = np.empty(n, dtype=float)
    restored[order] = adjusted
    out.loc[ok] = restored
    return out


def _scale_parameters(frame: pd.DataFrame, columns: list[str]) -> dict[str, tuple[float, float]]:
    result: dict[str, tuple[float, float]] = {}
    for column in columns:
        x = pd.to_numeric(frame[column], errors="coerce").to_numpy(float)
        mean = float(np.mean(x))
        sd = float(np.std(x, ddof=0))
        if not np.isfinite(sd) or sd <= 0:
            raise ValueError(f"constant or invalid predictor: {column}")
        result[column] = (mean, sd)
    return result


def _scaled(series: pd.Series, parameters: tuple[float, float]) -> np.ndarray:
    mean, sd = parameters
    x = pd.to_numeric(series, errors="coerce").to_numpy(float)
    return (x - mean) / sd


def _fit_beta_slope(
    frame: pd.DataFrame,
    *,
    geography: str,
    baseline: list[str],
    scaling: dict[str, tuple[float, float]],
    max_iter: int,
) -> tuple[float, bool]:
    names = ["intercept"]
    columns: list[np.ndarray] = [np.ones(len(frame), dtype=float)]
    for predictor in baseline:
        columns.append(_scaled(frame[predictor], scaling[predictor]))
        names.append(f"z_{predictor}")
    columns.append(_scaled(frame[geography], scaling[geography]))
    names.append(f"z_{geography}")
    fit = _fit_single_beta_binomial(
        frame["successes"].to_numpy(float),
        frame["trials"].to_numpy(float),
        np.column_stack(columns),
        names,
        max_iter=max_iter,
    )
    index = names.index(f"z_{geography}")
    estimate = float(fit["theta"][index])
    return estimate, bool(fit["success"])


def _jackknife_summary(full: float, leave_one_out: np.ndarray) -> dict[str, float | int]:
    x = np.asarray(leave_one_out, dtype=float)
    x = x[np.isfinite(x)]
    g = int(len(x))
    if g < 3:
        return {
            "jackknife_n_clusters": g,
            "jackknife_se": float("nan"),
            "t_value": float("nan"),
            "df": float("nan"),
            "p_two_sided": float("nan"),
            "p_one_sided_positive": float("nan"),
            "jackknife_bias_corrected_estimate": float("nan"),
            "loo_min": float("nan"),
            "loo_max": float("nan"),
            "loo_nonpositive_count": 0,
            "loo_nonpositive_fraction": float("nan"),
            "max_abs_loo_shift": float("nan"),
        }
    mean_loo = float(np.mean(x))
    variance = float((g - 1.0) / g * np.sum((x - mean_loo) ** 2))
    se = float(math.sqrt(max(variance, 0.0)))
    t_value = float(full / se) if se > 0 else float("nan")
    df = g - 1
    two = (
        float(2.0 * student_t.sf(abs(t_value), df))
        if math.isfinite(t_value)
        else float("nan")
    )
    one = (
        float(student_t.sf(t_value, df))
        if math.isfinite(t_value)
        else float("nan")
    )
    bias_corrected = float(g * full - (g - 1.0) * mean_loo)
    nonpositive = int(np.sum(x <= 0))
    return {
        "jackknife_n_clusters": g,
        "jackknife_se": se,
        "t_value": t_value,
        "df": int(df),
        "p_two_sided": two,
        "p_one_sided_positive": one,
        "jackknife_bias_corrected_estimate": bias_corrected,
        "loo_min": float(np.min(x)),
        "loo_max": float(np.max(x)),
        "loo_nonpositive_count": nonpositive,
        "loo_nonpositive_fraction": float(nonpositive / g),
        "max_abs_loo_shift": float(np.max(np.abs(x - full))),
    }


def _cluster_balance(parts: dict[str, pd.DataFrame], cluster: str) -> dict[str, float | int]:
    frames = [
        p[["island_id", cluster]].drop_duplicates()
        for p in parts.values()
        if not p.empty
    ]
    if not frames:
        return {}
    island_clusters = pd.concat(frames, ignore_index=True).drop_duplicates("island_id")
    sizes = island_clusters.groupby(cluster)["island_id"].nunique().to_numpy(float)
    total = float(np.sum(sizes))
    effective = (
        float(total * total / np.sum(sizes * sizes))
        if total > 0 and np.sum(sizes * sizes) > 0
        else float("nan")
    )
    mean = float(np.mean(sizes))
    cv = float(np.std(sizes, ddof=0) / mean) if mean > 0 else float("nan")
    return {
        "n_clusters_union": int(len(sizes)),
        "cluster_size_min_islands": int(np.min(sizes)),
        "cluster_size_median_islands": float(np.median(sizes)),
        "cluster_size_max_islands": int(np.max(sizes)),
        "cluster_size_cv": cv,
        "effective_cluster_count": effective,
    }


def _prepare_h1_parts(
    counts: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    stratum: str,
    context_value: str,
    outcomes: list[str],
    geography: str,
    context: str,
    cluster: str,
    baseline: list[str],
    minimum_islands: int,
    minimum_positive_islands: int,
    minimum_negative_islands: int,
) -> tuple[dict[str, pd.DataFrame], dict[str, dict[str, tuple[float, float]]], list[dict[str, Any]]]:
    cov_cols = ["island_id", geography, context, cluster, *baseline]
    merged = counts.loc[
        counts["stratum"].astype(str).eq(stratum)
        & counts["outcome"].astype(str).isin(outcomes)
    ].merge(
        covariates[cov_cols].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    for column in ["successes", "trials", geography, *baseline]:
        merged[column] = pd.to_numeric(merged[column], errors="coerce")
    merged[context] = merged[context].fillna("").astype(str)
    merged[cluster] = merged[cluster].fillna("").astype(str)
    merged = merged.loc[merged[context].eq(context_value)].dropna(
        subset=["successes", "trials", geography, *baseline]
    )
    merged = merged.loc[
        merged["trials"].gt(0)
        & merged["successes"].ge(0)
        & merged["successes"].le(merged["trials"])
        & merged[cluster].ne("")
    ].copy()

    parts: dict[str, pd.DataFrame] = {}
    scaling: dict[str, dict[str, tuple[float, float]]] = {}
    support: list[dict[str, Any]] = []
    for outcome in outcomes:
        part = merged.loc[merged["outcome"].astype(str).eq(outcome)].copy()
        n_islands = int(part["island_id"].nunique())
        positive = int(part.loc[part["successes"].gt(0), "island_id"].nunique())
        negative = int(
            part.loc[part["successes"].lt(part["trials"]), "island_id"].nunique()
        )
        n_clusters = int(part[cluster].nunique())
        eligible = (
            n_islands >= minimum_islands
            and positive >= minimum_positive_islands
            and negative >= minimum_negative_islands
            and n_clusters >= 3
        )
        support.append(
            {
                "stratum": stratum,
                "context": context_value,
                "outcome": outcome,
                "n_islands": n_islands,
                "positive_islands": positive,
                "negative_islands": negative,
                "n_clusters": n_clusters,
                "eligible": eligible,
            }
        )
        if not eligible:
            continue
        try:
            params = _scale_parameters(part, [*baseline, geography])
        except ValueError:
            continue
        parts[outcome] = part
        scaling[outcome] = params
    return parts, scaling, support


def _contrast_value(vector: dict[str, float], weights: dict[str, float]) -> float:
    return float(sum(float(weights[name]) * float(vector[name]) for name in weights))


def _fit_h1_vector(
    parts: dict[str, pd.DataFrame],
    scaling: dict[str, dict[str, tuple[float, float]]],
    *,
    outcomes: list[str],
    geography: str,
    baseline: list[str],
    cluster: str,
    max_iter: int,
    omit_cluster: str | None = None,
) -> tuple[dict[str, float], bool]:
    vector: dict[str, float] = {}
    converged = True
    for outcome in outcomes:
        part = parts[outcome]
        if omit_cluster is not None:
            part = part.loc[part[cluster].astype(str).ne(str(omit_cluster))].copy()
        if len(part) < 20:
            return {}, False
        estimate, success = _fit_beta_slope(
            part,
            geography=geography,
            baseline=baseline,
            scaling=scaling[outcome],
            max_iter=max_iter,
        )
        if not success:
            estimate_retry, success_retry = _fit_beta_slope(
                part,
                geography=geography,
                baseline=baseline,
                scaling=scaling[outcome],
                max_iter=max_iter * 3,
            )
            estimate = estimate_retry
            success = success_retry
        vector[outcome] = estimate
        converged = converged and success
    return vector, converged


def run_h1(
    counts: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
    *,
    evidence_scope: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    h1 = config["h1"]
    geography = str(config["geography_column"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    outcomes = [str(x) for x in h1["outcomes"]]
    contrasts = {
        str(name): {str(k): float(v) for k, v in weights.items()}
        for name, weights in h1["contrasts"].items()
    }
    rows: list[dict[str, Any]] = []
    support_rows: list[dict[str, Any]] = []
    slope_rows: list[dict[str, Any]] = []

    for stratum in [str(x) for x in h1["strata"]]:
        for context_value in [str(x) for x in config["contexts"]]:
            parts, scaling, support = _prepare_h1_parts(
                counts,
                covariates,
                stratum=stratum,
                context_value=context_value,
                outcomes=outcomes,
                geography=geography,
                context=context,
                cluster=cluster,
                baseline=baseline,
                minimum_islands=int(config["minimum_islands"]),
                minimum_positive_islands=int(config["minimum_positive_islands"]),
                minimum_negative_islands=int(config["minimum_negative_islands"]),
            )
            for item in support:
                item["evidence_scope"] = evidence_scope
                support_rows.append(item)

            if set(parts) != set(outcomes):
                for contrast_name in contrasts:
                    rows.append(
                        {
                            "evidence_scope": evidence_scope,
                            "stratum": stratum,
                            "context": context_value,
                            "contrast": contrast_name,
                            "status": "not_testable_all_frozen_outcomes_required",
                        }
                    )
                continue

            full_vector, full_converged = _fit_h1_vector(
                parts,
                scaling,
                outcomes=outcomes,
                geography=geography,
                baseline=baseline,
                cluster=cluster,
                max_iter=int(config["max_iter"]),
            )
            if set(full_vector) != set(outcomes):
                continue
            for outcome in outcomes:
                slope_rows.append(
                    {
                        "evidence_scope": evidence_scope,
                        "stratum": stratum,
                        "context": context_value,
                        "outcome": outcome,
                        "estimate": full_vector[outcome],
                        "optimizer_converged": full_converged,
                    }
                )

            clusters = sorted(
                set().union(
                    *[
                        set(part[cluster].astype(str).unique())
                        for part in parts.values()
                    ]
                )
            )
            loo_vectors: list[dict[str, float]] = []
            loo_clusters: list[str] = []
            loo_failures: list[str] = []
            for cluster_value in clusters:
                vector, converged = _fit_h1_vector(
                    parts,
                    scaling,
                    outcomes=outcomes,
                    geography=geography,
                    baseline=baseline,
                    cluster=cluster,
                    max_iter=int(config["max_iter"]),
                    omit_cluster=cluster_value,
                )
                if converged and set(vector) == set(outcomes):
                    loo_vectors.append(vector)
                    loo_clusters.append(cluster_value)
                else:
                    loo_failures.append(cluster_value)

            balance = _cluster_balance(parts, cluster)
            for contrast_name, weights in contrasts.items():
                required = set(weights)
                if not required.issubset(full_vector):
                    continue
                full_value = _contrast_value(full_vector, weights)
                loo = np.array(
                    [
                        _contrast_value(vector, weights)
                        for vector in loo_vectors
                        if required.issubset(vector)
                    ],
                    dtype=float,
                )
                summary = _jackknife_summary(full_value, loo)
                rows.append(
                    {
                        "evidence_scope": evidence_scope,
                        "stratum": stratum,
                        "context": context_value,
                        "contrast": contrast_name,
                        "status": (
                            "fit"
                            if len(loo_failures) == 0
                            else "fit_with_loo_failures"
                        ),
                        "estimate": full_value,
                        "full_vector_all_optimizers_converged": full_converged,
                        "n_expected_clusters": len(clusters),
                        "n_successful_loo_clusters": len(loo_clusters),
                        "n_failed_loo_clusters": len(loo_failures),
                        "failed_loo_clusters": "|".join(loo_failures),
                        **balance,
                        **summary,
                    }
                )

    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame["holm_p_one_sided_within_contrast"] = np.nan
        frame["holm_p_two_sided_within_contrast"] = np.nan
        fit = frame["status"].astype(str).str.startswith("fit")
        for (_, stratum, contrast_name), index in frame.loc[fit].groupby(
            ["evidence_scope", "stratum", "contrast"]
        ).groups.items():
            idx = list(index)
            frame.loc[idx, "holm_p_one_sided_within_contrast"] = _holm(
                frame.loc[idx, "p_one_sided_positive"]
            )
            frame.loc[idx, "holm_p_two_sided_within_contrast"] = _holm(
                frame.loc[idx, "p_two_sided"]
            )
        frame["directional_supported"] = (
            frame["estimate"].gt(0)
            & frame["holm_p_one_sided_within_contrast"].le(
                float(config["alpha"])
            )
        ).fillna(False)
    return frame, pd.DataFrame(support_rows), pd.DataFrame(slope_rows)


def _ols_fit(
    frame: pd.DataFrame,
    *,
    response: str,
    predictors: list[str],
    scaling: dict[str, tuple[float, float]],
) -> tuple[np.ndarray, list[str]]:
    names = ["intercept", *[f"z_{p}" for p in predictors]]
    columns = [np.ones(len(frame), dtype=float)]
    for predictor in predictors:
        columns.append(_scaled(frame[predictor], scaling[predictor]))
    x = np.column_stack(columns)
    y = pd.to_numeric(frame[response], errors="coerce").to_numpy(float)
    beta = np.linalg.pinv(x.T @ x) @ x.T @ y
    return beta, names


def _run_h2_model(
    frame: pd.DataFrame,
    *,
    response: str,
    predictors: list[str],
    geography: str,
    cluster: str,
    minimum_clusters: int,
) -> dict[str, Any]:
    required = [response, *predictors, cluster, "island_id"]
    work = frame[required].copy()
    for column in [response, *predictors]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work[cluster] = work[cluster].fillna("").astype(str)
    work = work.dropna(subset=[response, *predictors])
    work = work.loc[work[cluster].ne("")].copy()
    clusters = sorted(work[cluster].unique())
    if len(clusters) < minimum_clusters or len(work) < 20:
        return {
            "status": "not_testable",
            "n_islands": int(work["island_id"].nunique()),
            "n_clusters": int(len(clusters)),
        }
    try:
        scaling = _scale_parameters(work, predictors)
    except ValueError:
        return {"status": "not_testable", "n_islands": int(len(work)), "n_clusters": len(clusters)}
    beta, names = _ols_fit(
        work, response=response, predictors=predictors, scaling=scaling
    )
    index = names.index(f"z_{geography}")
    full = float(beta[index])
    loo: list[float] = []
    failed: list[str] = []
    for value in clusters:
        subset = work.loc[work[cluster].ne(value)]
        try:
            beta_loo, names_loo = _ols_fit(
                subset, response=response, predictors=predictors, scaling=scaling
            )
            loo.append(float(beta_loo[names_loo.index(f"z_{geography}")]))
        except (ValueError, np.linalg.LinAlgError):
            failed.append(str(value))
    summary = _jackknife_summary(full, np.asarray(loo, dtype=float))
    sizes = (
        work[["island_id", cluster]]
        .drop_duplicates()
        .groupby(cluster)["island_id"]
        .nunique()
        .to_numpy(float)
    )
    total = float(np.sum(sizes))
    effective = float(total * total / np.sum(sizes * sizes))
    return {
        "status": "fit" if not failed else "fit_with_loo_failures",
        "n_islands": int(work["island_id"].nunique()),
        "n_clusters": int(len(clusters)),
        "estimate": full,
        "n_failed_loo_clusters": len(failed),
        "failed_loo_clusters": "|".join(failed),
        "cluster_size_cv": float(np.std(sizes, ddof=0) / np.mean(sizes)),
        "effective_cluster_count": effective,
        **summary,
    }


def run_h2(
    syndrome_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
    *,
    evidence_scope: str,
) -> pd.DataFrame:
    h2 = config["h2"]
    geography = str(config["geography_column"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    plant = build_plant_scores(syndrome_scores, str(h2["stratum"]))
    cov_cols = ["island_id", geography, context, cluster, *baseline]
    merged = plant.merge(
        covariates[cov_cols].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="one_to_one",
    )
    merged[context] = merged[context].fillna("").astype(str)
    rows: list[dict[str, Any]] = []
    models = [
        (
            "H2a_reproductive_assurance",
            "selfing_core",
            [geography, *baseline],
        ),
        (
            "H2b_accessibility_given_assurance",
            "generalized_accessible",
            [geography, "selfing_core", *baseline],
        ),
    ]
    for context_value in [str(x) for x in config["contexts"]]:
        part = merged.loc[merged[context].eq(context_value)].copy()
        for role, response, predictors in models:
            result = _run_h2_model(
                part,
                response=response,
                predictors=predictors,
                geography=geography,
                cluster=cluster,
                minimum_clusters=int(h2["minimum_clusters"]),
            )
            rows.append(
                {
                    "evidence_scope": evidence_scope,
                    "context": context_value,
                    "analysis_role": role,
                    "response": response,
                    **result,
                }
            )
    frame = pd.DataFrame(rows)
    frame["holm_p_one_sided_within_pathway"] = np.nan
    frame["holm_p_two_sided_within_pathway"] = np.nan
    fit = frame["status"].astype(str).str.startswith("fit")
    for (_, role), index in frame.loc[fit].groupby(
        ["evidence_scope", "analysis_role"]
    ).groups.items():
        idx = list(index)
        frame.loc[idx, "holm_p_one_sided_within_pathway"] = _holm(
            frame.loc[idx, "p_one_sided_positive"]
        )
        frame.loc[idx, "holm_p_two_sided_within_pathway"] = _holm(
            frame.loc[idx, "p_two_sided"]
        )
    frame["directional_supported"] = (
        frame["estimate"].gt(0)
        & frame["holm_p_one_sided_within_pathway"].le(float(config["alpha"]))
    ).fillna(False)
    return frame


def _load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if (
        not isinstance(config, dict)
        or config.get("contract") != "chapter1_final_reviewer_robust_v1"
    ):
        raise typer.BadParameter("unexpected final reviewer-robust contract")
    return config


@app.command("run")
def run(
    status_flora_csv: Path = typer.Option(..., exists=True),
    state_audit_csv: Path = typer.Option(..., exists=True),
    syndrome_scores_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    v14_h1_config_path: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
) -> None:
    config = _load_config(config_path)
    v14_h1 = yaml.safe_load(v14_h1_config_path.read_text(encoding="utf-8"))
    status_flora = pd.read_csv(status_flora_csv)
    state_audit = pd.read_csv(state_audit_csv)
    covariates = pd.read_csv(covariates_csv)
    syndrome_scores = pd.read_csv(syndrome_scores_csv)

    counts_config = dict(v14_h1)
    counts_config["strata"] = [str(x) for x in config["h1"]["strata"]]
    counts = build_broad_counts(status_flora, state_audit, counts_config)

    h1, h1_support, h1_slopes = run_h1(
        counts, covariates, config, evidence_scope=evidence_scope
    )
    h2 = run_h2(
        syndrome_scores, covariates, config, evidence_scope=evidence_scope
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    h1.to_csv(output_dir / "h1_directional_cluster_jackknife.csv", index=False)
    h1_support.to_csv(output_dir / "h1_outcome_support.csv", index=False)
    h1_slopes.to_csv(output_dir / "h1_atomic_slopes.csv", index=False)
    h2.to_csv(output_dir / "h2_cluster_jackknife.csv", index=False)
    counts.to_csv(
        output_dir / "h1_counts.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    primary = h1.loc[
        h1["stratum"].eq("all_observed")
        & h1["contrast"].eq("classic_equal_domain")
    ].copy()
    summary = {
        "contract": config["contract"],
        "evidence_scope": evidence_scope,
        "h1_primary": {
            "n_contexts_tested": int(primary["status"].astype(str).str.startswith("fit").sum()),
            "n_contexts_positive": int(primary["estimate"].gt(0).sum()),
            "n_contexts_directional_supported_after_holm": int(
                primary["directional_supported"].sum()
            ),
            "all_four_positive": bool(len(primary) == 4 and primary["estimate"].gt(0).all()),
            "all_four_directional_supported_after_holm": bool(
                len(primary) == 4 and primary["directional_supported"].all()
            ),
        },
        "claim_boundary": config["claim_boundary"],
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(summary, indent=2))


if __name__ == "__main__":
    app()
