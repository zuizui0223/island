"""Finite-cluster directional H1 inference for Chapter 1.

The confirmatory estimand is the pre-oriented seven-indicator classic island-syndrome
score already frozen in v14. Indicators are first averaged within their three biological
domains and the three domain means are then weighted equally. This prevents domains
with more measured indicators from receiving more inferential weight.

The module deliberately separates that directional H1 test from the newer three-axis
raw-state analysis. Raw-state models remain useful for describing how each measurement
domain reorganizes, but their high-dimensional direction-free Wald tests are not used as
confirmatory evidence for a classic island syndrome.

Inference is one-dimensional. We report a cluster-robust standard error, a finite-cluster
t test with G-1 degrees of freedom, and a deterministic linearized wild-cluster sign-flip
sensitivity. The four-region recurrence claim is an intersection-union test: all four
predeclared regional directional effects must be positive and individually supported.
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
    _prepare,
    _standardize,
    build_broad_counts,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def validate_score_weights(config: dict[str, Any]) -> dict[str, float]:
    weights = {
        str(k): float(v)
        for k, v in config["directional_score"]["outcome_weights"].items()
    }
    expected = [str(x) for x in config["model_outcomes"]]
    if set(weights) != set(expected):
        missing = sorted(set(expected) - set(weights))
        extra = sorted(set(weights) - set(expected))
        raise typer.BadParameter(
            "directional score weights do not match model outcomes; "
            f"missing={missing}, extra={extra}"
        )
    total = float(sum(weights.values()))
    if not math.isclose(total, 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise typer.BadParameter(
            f"directional score weights must sum to 1, got {total}"
        )

    families = config["response_families"]
    target = 1.0 / len(families)
    for family, outcomes in families.items():
        family_weight = float(sum(weights[str(x)] for x in outcomes))
        if not math.isclose(
            family_weight,
            target,
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise typer.BadParameter(
                f"family {family} has weight {family_weight}, expected {target}"
            )
    return weights


def _stack_cluster_influences(
    fits: list[dict[str, Any]],
    clusters: list[np.ndarray],
    slope_global_indices: list[int],
    slope_weights: np.ndarray,
) -> tuple[np.ndarray, float, int, int, list[str]]:
    offsets: list[int] = []
    total_size = 0
    for fit in fits:
        offsets.append(total_size)
        total_size += len(fit["names"])

    bread = np.zeros((total_size, total_size), dtype=float)
    cluster_scores: dict[str, np.ndarray] = {}
    total_rows = 0
    for fit, labels, offset in zip(fits, clusters, offsets, strict=True):
        dim = len(fit["names"])
        bread[offset : offset + dim, offset : offset + dim] = fit["bread"]
        labels = np.asarray(labels).astype(str)
        total_rows += len(labels)
        for cluster in np.unique(labels):
            if cluster not in cluster_scores:
                cluster_scores[cluster] = np.zeros(total_size, dtype=float)
            cluster_scores[cluster][offset : offset + dim] += fit["score"][
                labels == cluster
            ].sum(axis=0)

    contrast = np.zeros(total_size, dtype=float)
    for index, weight in zip(
        slope_global_indices,
        slope_weights,
        strict=True,
    ):
        contrast[index] = float(weight)

    cluster_labels = list(cluster_scores)
    influences = np.asarray(
        [
            float(contrast @ bread @ cluster_scores[label])
            for label in cluster_labels
        ],
        dtype=float,
    )
    g = int(len(influences))
    correction = 1.0
    if g > 1 and total_rows > total_size:
        correction = (g / (g - 1.0)) * (
            (total_rows - 1.0) / (total_rows - total_size)
        )
    return (
        influences,
        float(correction),
        g,
        int(total_size),
        cluster_labels,
    )


def _wild_signflip_p(
    estimate: float,
    se: float,
    influences: np.ndarray,
    *,
    replications: int,
    seed: int,
) -> float:
    if not math.isfinite(estimate) or not math.isfinite(se) or se <= 0:
        return float("nan")
    if len(influences) < 2 or replications < 1:
        return float("nan")
    rng = np.random.default_rng(int(seed))
    observed = float(estimate / se)
    batch = 2000
    exceed = 0
    done = 0
    while done < int(replications):
        n = min(batch, int(replications) - done)
        multipliers = rng.choice(
            np.array([-1.0, 1.0]),
            size=(n, len(influences)),
        )
        t_star = (multipliers @ influences) / se
        exceed += int(np.count_nonzero(t_star >= observed))
        done += n
    return float((exceed + 1.0) / (int(replications) + 1.0))


def _fit_directional_score(
    prepared: pd.DataFrame,
    *,
    stratum: str,
    context_value: str,
    config: dict[str, Any],
    weights: dict[str, float],
    seed: int,
) -> tuple[dict[str, Any], pd.DataFrame]:
    geography = str(config["geography_column"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    outcomes = [str(x) for x in config["model_outcomes"]]
    threshold = int(config["minimum_islands_per_outcome"])

    work = prepared.loc[
        prepared["stratum"].eq(stratum)
        & prepared[context].eq(context_value)
    ].copy()
    support = work.groupby("outcome")["island_id"].nunique()
    missing = [x for x in outcomes if int(support.get(x, 0)) < threshold]
    if missing:
        return {
            "stratum": stratum,
            "context": context_value,
            "status": "not_testable",
            "missing_or_under_supported_outcomes": "|".join(missing),
            "n_clusters": int(work[cluster].nunique()),
        }, pd.DataFrame()

    fits: list[dict[str, Any]] = []
    cluster_parts: list[np.ndarray] = []
    slope_indices: list[int] = []
    rows: list[dict[str, Any]] = []
    offset = 0
    max_iter = int(config.get("max_iter", 1000))
    retry_max_iter = int(config.get("retry_max_iter", max_iter))

    for outcome in outcomes:
        part = work.loc[work["outcome"].eq(outcome)].copy()
        columns = [np.ones(len(part), dtype=float)]
        names = [f"{outcome}:intercept"]
        for predictor in baseline:
            columns.append(_standardize(part[predictor]))
            names.append(f"{outcome}:z_{predictor}")
        columns.append(_standardize(part[geography]))
        slope_name = f"{outcome}:z_{geography}"
        names.append(slope_name)
        design = np.column_stack(columns)
        fit = _fit_single_beta_binomial(
            part["successes"].to_numpy(float),
            part["trials"].to_numpy(float),
            design,
            names,
            max_iter=max_iter,
        )
        retry_used = False
        if not bool(fit["success"]) and retry_max_iter > max_iter:
            fit = _fit_single_beta_binomial(
                part["successes"].to_numpy(float),
                part["trials"].to_numpy(float),
                design,
                names,
                max_iter=retry_max_iter,
            )
            retry_used = True
        fits.append(fit)
        cluster_parts.append(part[cluster].to_numpy(str))
        local_slope = names.index(slope_name)
        slope_indices.append(offset + local_slope)
        estimate = float(fit["theta"][local_slope])
        domain = next(
            family
            for family, members in config["response_families"].items()
            if outcome in members
        )
        rows.append(
            {
                "stratum": stratum,
                "context": context_value,
                "outcome": outcome,
                "domain": domain,
                "directional_weight": float(weights[outcome]),
                "geography_slope_log_odds": estimate,
                "weighted_contribution": float(weights[outcome])
                * estimate,
                "n_islands": int(part["island_id"].nunique()),
                "optimizer_success": bool(fit["success"]),
                "retry_used": retry_used,
            }
        )
        offset += len(fit["names"])

    if not all(bool(fit["success"]) for fit in fits):
        failed = [
            outcome
            for outcome, fit in zip(outcomes, fits, strict=True)
            if not fit["success"]
        ]
        return {
            "stratum": stratum,
            "context": context_value,
            "status": "optimizer_failure",
            "failed_outcomes": "|".join(failed),
            "n_clusters": int(work[cluster].nunique()),
        }, pd.DataFrame(rows)

    slopes = np.asarray(
        [
            float(
                fit["theta"][
                    fit["names"].index(f"{outcome}:z_{geography}")
                ]
            )
            for outcome, fit in zip(outcomes, fits, strict=True)
        ],
        dtype=float,
    )
    weight_vector = np.asarray(
        [weights[x] for x in outcomes],
        dtype=float,
    )
    estimate = float(weight_vector @ slopes)
    (
        influences,
        correction,
        g,
        parameter_dimension,
        cluster_labels,
    ) = _stack_cluster_influences(
        fits,
        cluster_parts,
        slope_indices,
        weight_vector,
    )
    variance = float(correction * np.sum(np.square(influences)))
    se = float(math.sqrt(max(variance, 0.0)))
    t_value = float(estimate / se) if se > 0 else float("nan")
    df = int(g - 1) if g > 1 else 0
    p_one_t = (
        float(student_t.sf(t_value, df=df))
        if df > 0 and math.isfinite(t_value)
        else float("nan")
    )
    p_two_t = (
        float(2.0 * student_t.sf(abs(t_value), df=df))
        if df > 0 and math.isfinite(t_value)
        else float("nan")
    )
    squared = np.square(influences)
    total_sq = float(np.sum(squared))
    fourth_sum = float(np.sum(np.square(squared)))
    effective_clusters = (
        float(total_sq * total_sq / fourth_sum)
        if total_sq > 0 and fourth_sum > 0
        else float("nan")
    )
    max_cluster_variance_share = (
        float(np.max(squared) / total_sq)
        if total_sq > 0
        else float("nan")
    )
    top_cluster = (
        str(cluster_labels[int(np.argmax(squared))])
        if total_sq > 0 and cluster_labels
        else ""
    )
    p_wild = _wild_signflip_p(
        estimate,
        se,
        influences * math.sqrt(correction),
        replications=int(
            config["directional_score"]["wild_cluster_replications"]
        ),
        seed=int(seed),
    )
    alpha = float(config["alpha"])
    return {
        "stratum": stratum,
        "context": context_value,
        "status": "fit",
        "estimate": estimate,
        "cluster_robust_se": se,
        "t_value": t_value,
        "finite_cluster_df": df,
        "p_one_sided_t": p_one_t,
        "p_two_sided_t": p_two_t,
        "p_one_sided_wild": p_wild,
        "n_unique_islands": int(
            work.loc[
                work["outcome"].isin(outcomes),
                "island_id",
            ].nunique()
        ),
        "n_clusters": g,
        "effective_clusters": effective_clusters,
        "max_cluster_variance_share": max_cluster_variance_share,
        "top_influence_cluster": top_cluster,
        "sandwich_finite_sample_correction": correction,
        "stacked_parameter_dimension": parameter_dimension,
        "all_optimizers_converged": True,
        "positive_direction": bool(estimate > 0),
        "t_supported": bool(estimate > 0 and p_one_t <= alpha),
        "wild_supported": bool(estimate > 0 and p_wild <= alpha),
        "robust_supported": bool(
            estimate > 0
            and p_one_t <= alpha
            and p_wild <= alpha
        ),
    }, pd.DataFrame(rows)


def meta_directional_summary(
    results: pd.DataFrame,
    *,
    contexts: list[str],
    alpha: float,
) -> dict[str, Any]:
    """Random-effects summary of the four predeclared regional directional scores.

    The regional estimates are based on disjoint geographic strata. A
    DerSimonian-Laird tau-squared estimate is paired with a modified
    Knapp-Hartung variance inflation (never below the conventional random-effects
    variance) and a t reference with k-1 degrees of freedom. This summary asks
    whether there is a positive global average directional component while
    explicitly quantifying between-region heterogeneity. It is not a substitute
    for the stricter intersection-union recurrence test.
    """
    fit = results.loc[
        results["status"].eq("fit")
        & results["context"].isin(contexts)
    ].copy()
    if len(fit) != len(contexts) or set(fit["context"]) != set(contexts):
        return {
            "status": "not_testable",
            "n_required_contexts": len(contexts),
            "n_fitted_contexts": int(len(fit)),
        }
    y = pd.to_numeric(fit["estimate"], errors="coerce").to_numpy(float)
    se = pd.to_numeric(
        fit["cluster_robust_se"], errors="coerce"
    ).to_numpy(float)
    if (
        not np.isfinite(y).all()
        or not np.isfinite(se).all()
        or np.any(se <= 0)
    ):
        return {
            "status": "not_testable",
            "reason": "nonfinite_regional_estimate_or_se",
        }
    k = len(y)
    w = 1.0 / np.square(se)
    fixed = float(np.sum(w * y) / np.sum(w))
    fixed_se = float(math.sqrt(1.0 / np.sum(w)))
    q = float(np.sum(w * np.square(y - fixed)))
    q_df = int(k - 1)
    c_value = float(np.sum(w) - np.sum(np.square(w)) / np.sum(w))
    tau2 = float(max(0.0, (q - q_df) / c_value)) if c_value > 0 else 0.0
    wr = 1.0 / (np.square(se) + tau2)
    random_estimate = float(np.sum(wr * y) / np.sum(wr))
    hk_scale = float(
        np.sum(wr * np.square(y - random_estimate)) / q_df
    ) if q_df > 0 else float("nan")
    modified_hk_scale = float(max(1.0, hk_scale))
    random_se = float(
        math.sqrt(modified_hk_scale / np.sum(wr))
    )
    t_value = (
        float(random_estimate / random_se)
        if random_se > 0
        else float("nan")
    )
    p_one = (
        float(student_t.sf(t_value, df=q_df))
        if q_df > 0 and math.isfinite(t_value)
        else float("nan")
    )
    p_two = (
        float(2.0 * student_t.sf(abs(t_value), df=q_df))
        if q_df > 0 and math.isfinite(t_value)
        else float("nan")
    )
    i2 = (
        float(max(0.0, (q - q_df) / q))
        if q > 0
        else 0.0
    )
    return {
        "status": "fit",
        "n_contexts": int(k),
        "all_context_estimates_positive": bool(np.all(y > 0)),
        "fixed_effect_estimate": fixed,
        "fixed_effect_se": fixed_se,
        "cochran_q": q,
        "heterogeneity_df": q_df,
        "i2": i2,
        "tau2_dl": tau2,
        "random_effects_estimate": random_estimate,
        "modified_knapp_hartung_se": random_se,
        "modified_knapp_hartung_scale": modified_hk_scale,
        "random_effects_t": t_value,
        "random_effects_p_one_sided": p_one,
        "random_effects_p_two_sided": p_two,
        "positive_global_average_supported": bool(
            random_estimate > 0 and p_one <= alpha
        ),
        "interpretation": (
            "positive global average component with explicit regional heterogeneity; "
            "not evidence that every region independently supports the syndrome"
        ),
    }


def intersection_union_summary(
    results: pd.DataFrame,
    *,
    contexts: list[str],
    alpha: float,
) -> dict[str, Any]:
    fit = results.loc[
        results["status"].eq("fit")
        & results["context"].isin(contexts)
    ].copy()
    if len(fit) != len(contexts) or set(fit["context"]) != set(contexts):
        return {
            "status": "not_testable",
            "n_required_contexts": len(contexts),
            "n_fitted_contexts": int(len(fit)),
        }
    p_t = float(fit["p_one_sided_t"].max())
    p_wild = float(fit["p_one_sided_wild"].max())
    all_positive = bool(fit["positive_direction"].all())
    return {
        "status": "fit",
        "n_required_contexts": len(contexts),
        "n_fitted_contexts": int(len(fit)),
        "all_context_estimates_positive": all_positive,
        "iut_p_one_sided_t": p_t,
        "iut_p_one_sided_wild": p_wild,
        "recurrent_supported_t": bool(
            all_positive and p_t <= alpha
        ),
        "recurrent_supported_wild": bool(
            all_positive and p_wild <= alpha
        ),
        "recurrent_robust_supported": bool(
            all_positive
            and p_t <= alpha
            and p_wild <= alpha
        ),
    }


def run_final_directional_h1(
    status_flora: pd.DataFrame,
    state_audit: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
    *,
    evidence_scope: str,
) -> dict[str, pd.DataFrame | dict[str, Any]]:
    weights = validate_score_weights(config)
    counts = build_broad_counts(status_flora, state_audit, config)
    prepared = _prepare(counts, covariates, config)
    rows: list[dict[str, Any]] = []
    component_parts: list[pd.DataFrame] = []
    contexts = [str(x) for x in config["contexts"]]
    strata = [str(x) for x in config["strata"]]
    base_seed = int(config["directional_score"]["wild_cluster_seed"])
    counter = 0
    for stratum in strata:
        for context in contexts:
            row, components = _fit_directional_score(
                prepared,
                stratum=stratum,
                context_value=context,
                config=config,
                weights=weights,
                seed=base_seed + counter,
            )
            row["evidence_scope"] = evidence_scope
            rows.append(row)
            if not components.empty:
                components.insert(0, "evidence_scope", evidence_scope)
                component_parts.append(components)
            counter += 1
    results = pd.DataFrame(rows)
    primary_stratum = str(config["flora_roles"]["broad_primary"])
    primary = results.loc[
        results["stratum"].eq(primary_stratum)
    ].copy()
    iut = intersection_union_summary(
        primary,
        contexts=contexts,
        alpha=float(config["alpha"]),
    )
    iut["evidence_scope"] = evidence_scope
    iut["stratum"] = primary_stratum
    meta = meta_directional_summary(
        primary,
        contexts=contexts,
        alpha=float(config["alpha"]),
    )
    meta["evidence_scope"] = evidence_scope
    meta["stratum"] = primary_stratum

    cluster_column = str(config["cluster_column"])
    context_column = str(config["context_column"])
    leaveout_rows: list[dict[str, Any]] = []
    for index, parent in primary.reset_index(drop=True).iterrows():
        if parent.get("status") != "fit":
            continue
        top_cluster = str(parent.get("top_influence_cluster", ""))
        context_value = str(parent["context"])
        if not top_cluster:
            continue
        filtered = prepared.loc[
            ~(
                prepared[context_column].astype(str).eq(context_value)
                & prepared[cluster_column].astype(str).eq(top_cluster)
            )
        ].copy()
        row, _ = _fit_directional_score(
            filtered,
            stratum=primary_stratum,
            context_value=context_value,
            config=config,
            weights=weights,
            seed=base_seed + 1000 + index,
        )
        row["evidence_scope"] = evidence_scope
        row["sensitivity"] = "drop_top_variance_cluster"
        row["excluded_cluster"] = top_cluster
        row["parent_max_cluster_variance_share"] = float(
            parent.get("max_cluster_variance_share", float("nan"))
        )
        leaveout_rows.append(row)
    leaveout = pd.DataFrame(leaveout_rows)
    leaveout_iut = intersection_union_summary(
        leaveout,
        contexts=contexts,
        alpha=float(config["alpha"]),
    )
    leaveout_iut["evidence_scope"] = evidence_scope
    leaveout_iut["stratum"] = primary_stratum
    leaveout_iut["sensitivity"] = "drop_top_variance_cluster"
    leaveout_meta = meta_directional_summary(
        leaveout,
        contexts=contexts,
        alpha=float(config["alpha"]),
    )
    leaveout_meta["evidence_scope"] = evidence_scope
    leaveout_meta["stratum"] = primary_stratum
    leaveout_meta["sensitivity"] = "drop_top_variance_cluster"

    return {
        "results": results,
        "components": (
            pd.concat(component_parts, ignore_index=True)
            if component_parts
            else pd.DataFrame()
        ),
        "counts": counts,
        "iut": iut,
        "meta": meta,
        "top_cluster_sensitivity": leaveout,
        "top_cluster_iut": leaveout_iut,
        "top_cluster_meta": leaveout_meta,
    }


@app.command("run")
def run(
    status_flora_csv: Path = typer.Option(..., exists=True),
    state_audit_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
) -> None:
    config = yaml.safe_load(
        config_path.read_text(encoding="utf-8")
    )
    result = run_final_directional_h1(
        pd.read_csv(status_flora_csv, dtype=str).fillna(""),
        pd.read_csv(state_audit_csv, dtype=str).fillna(""),
        pd.read_csv(covariates_csv),
        config,
        evidence_scope=evidence_scope,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    result["results"].to_csv(
        output_dir / "directional_score_results.csv",
        index=False,
    )
    result["components"].to_csv(
        output_dir / "directional_score_components.csv",
        index=False,
    )
    result["counts"].to_csv(
        output_dir / "directional_score_counts.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    result["top_cluster_sensitivity"].to_csv(
        output_dir / "top_cluster_leaveout_results.csv",
        index=False,
    )
    (output_dir / "directional_score_iut.json").write_text(
        json.dumps(result["iut"], indent=2) + "\n",
        encoding="utf-8",
    )
    (output_dir / "directional_score_meta.json").write_text(
        json.dumps(result["meta"], indent=2) + "\n",
        encoding="utf-8",
    )
    (output_dir / "top_cluster_leaveout_iut.json").write_text(
        json.dumps(result["top_cluster_iut"], indent=2) + "\n",
        encoding="utf-8",
    )
    (output_dir / "top_cluster_leaveout_meta.json").write_text(
        json.dumps(result["top_cluster_meta"], indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(
        json.dumps(
            {
                "primary_iut": result["iut"],
                "primary_meta": result["meta"],
                "top_cluster_leaveout_iut": result["top_cluster_iut"],
                "top_cluster_leaveout_meta": result["top_cluster_meta"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    app()
