"""All-data probability analysis for Chapter 1 island floral composition.

This module is an additive analysis surface. It does not overwrite the frozen v10
submission results. The intended hierarchy is:

1. all-analysis-eligible trait evidence (High/Medium + validated Low) is the primary
   evidence scope;
2. direct-only High/Medium evidence is a sensitivity analysis;
3. all observed island floras are analysed as the broad sampling-frame result;
4. all-native and native-nonendemic floras are status-resolved sensitivities;
5. beta-binomial regression is the primary probability model for island trait counts;
6. grouped-binomial regression remains a model-form sensitivity elsewhere in the repo.

The all-observed result describes composition of the observed island flora. It cannot,
by itself, identify native colonisation, in-situ evolution, or historical assembly because
introduced and unresolved-status records are retained by construction.
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
from scipy.optimize import minimize
from scipy.special import betaln, digamma, gammaln

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _truthy(series: pd.Series) -> pd.Series:
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def _bh(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=values.index, dtype=float)
    ok = p.notna()
    if not ok.any():
        return out
    x = p.loc[ok].to_numpy(float)
    order = np.argsort(x)
    ranked = x[order]
    n = len(ranked)
    adjusted = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    restored = np.empty(n, dtype=float)
    restored[order] = np.clip(adjusted, 0.0, 1.0)
    out.loc[ok] = restored
    return out


def _chi_square_sf_integer_df(statistic: float, df: int) -> float:
    if df <= 0 or not math.isfinite(statistic) or statistic < 0:
        return float("nan")
    x = statistic / 2.0
    if df % 2 == 0:
        s = 1.0
        q = math.exp(-x)
    else:
        s = 0.5
        q = math.erfc(math.sqrt(x))
    target = df / 2.0
    while s < target - 1e-12:
        q += math.exp(s * math.log(x) - x - math.lgamma(s + 1.0)) if x > 0 else 0.0
        s += 1.0
    return float(min(max(q, 0.0), 1.0))


def _normal_two_sided_p(z: float) -> float:
    return math.erfc(abs(float(z)) / math.sqrt(2.0)) if math.isfinite(z) else float("nan")


def _classify_signature(value: object, positive: set[str], negative: set[str]) -> float:
    tokens = {x.strip() for x in str(value or "").split("|") if x.strip()}
    if not tokens:
        return float("nan")
    if tokens <= positive:
        return 1.0
    if tokens <= negative:
        return 0.0
    return float("nan")


def _stratum_mask(frame: pd.DataFrame, stratum: str) -> pd.Series:
    if stratum == "all_observed":
        return pd.Series(True, index=frame.index)
    if stratum == "all_native":
        return frame["origin_status"].astype(str).eq("native")
    if stratum == "native_nonendemic":
        return frame["floristic_status"].astype(str).eq("native_nonendemic")
    if stratum == "endemic":
        return frame["floristic_status"].astype(str).eq("endemic")
    raise typer.BadParameter(f"unknown flora stratum: {stratum}")


def build_broad_counts(
    status_flora: pd.DataFrame,
    state_audit: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    required_flora = {
        "island_id",
        "accepted_species",
        "origin_status",
        "floristic_status",
    }
    required_audit = {
        "accepted_species",
        "trait_name",
        "resolved_for_primary",
        "canonical_signature",
    }
    missing = required_flora - set(status_flora.columns)
    if missing:
        raise typer.BadParameter(f"status flora missing columns: {sorted(missing)}")
    missing = required_audit - set(state_audit.columns)
    if missing:
        raise typer.BadParameter(f"state audit missing columns: {sorted(missing)}")

    flora = status_flora.copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    audit = state_audit.loc[_truthy(state_audit["resolved_for_primary"])].copy()
    audit["accepted_species"] = audit["accepted_species"].astype(str)
    audit["trait_name"] = audit["trait_name"].astype(str)

    rows: list[pd.DataFrame] = []
    for outcome, spec in config["broad_outcomes"].items():
        states = audit.loc[
            audit["trait_name"].eq(str(spec["trait_name"])),
            ["accepted_species", "canonical_signature"],
        ].drop_duplicates("accepted_species")
        positive = {str(x) for x in spec["positive_states"]}
        negative = {str(x) for x in spec["negative_states"]}
        states = states.copy()
        states["state"] = [
            _classify_signature(v, positive, negative)
            for v in states["canonical_signature"]
        ]
        states = states.dropna(subset=["state"])[["accepted_species", "state"]]
        joined = flora.merge(states, on="accepted_species", how="inner", validate="many_to_one")
        for stratum in [str(x) for x in config["strata"]]:
            subset = joined.loc[
                _stratum_mask(joined, stratum),
                ["island_id", "accepted_species", "state"],
            ].drop_duplicates(["island_id", "accepted_species"])
            if subset.empty:
                continue
            summary = (
                subset.groupby("island_id", as_index=False)
                .agg(successes=("state", "sum"), trials=("state", "size"))
            )
            summary["outcome"] = str(outcome)
            summary["stratum"] = stratum
            rows.append(summary)
    if not rows:
        return pd.DataFrame(
            columns=["island_id", "successes", "trials", "outcome", "stratum"]
        )
    return pd.concat(rows, ignore_index=True)


def _standardize(series: pd.Series) -> np.ndarray:
    x = pd.to_numeric(series, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def _beta_binomial_value_score(
    theta: np.ndarray,
    y: np.ndarray,
    n: np.ndarray,
    design: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    p = design.shape[1]
    beta = theta[:p]
    kappa = float(np.exp(np.clip(theta[p], -6.0, 12.0)))
    eta = np.clip(design @ beta, -30.0, 30.0)
    prob = 1.0 / (1.0 + np.exp(-eta))
    alpha = np.clip(prob * kappa, 1e-10, None)
    beta_shape = np.clip((1.0 - prob) * kappa, 1e-10, None)
    loglik = (
        gammaln(n + 1.0)
        - gammaln(y + 1.0)
        - gammaln(n - y + 1.0)
        + betaln(y + alpha, n - y + beta_shape)
        - betaln(alpha, beta_shape)
    )
    dldp = kappa * (
        digamma(y + alpha)
        - digamma(alpha)
        - digamma(n - y + beta_shape)
        + digamma(beta_shape)
    )
    score_eta = dldp * prob * (1.0 - prob)
    score_beta = design * score_eta[:, None]
    dldk = (
        prob * (digamma(y + alpha) - digamma(alpha))
        + (1.0 - prob) * (digamma(n - y + beta_shape) - digamma(beta_shape))
        - digamma(n + kappa)
        + digamma(kappa)
    )
    score_logk = (kappa * dldk)[:, None]
    score = np.concatenate([score_beta, score_logk], axis=1)
    return loglik, score


def _fit_single_beta_binomial(
    y: np.ndarray,
    n: np.ndarray,
    design: np.ndarray,
    names: list[str],
    *,
    max_iter: int,
) -> dict[str, Any]:
    y = np.asarray(y, dtype=float)
    n = np.asarray(n, dtype=float)
    design = np.asarray(design, dtype=float)
    p = design.shape[1]
    start = np.zeros(p + 1, dtype=float)
    start[p] = math.log(50.0)

    def value_grad(theta: np.ndarray) -> tuple[float, np.ndarray]:
        ll, score = _beta_binomial_value_score(theta, y, n, design)
        return -float(np.sum(ll)), -np.sum(score, axis=0)

    result = minimize(
        lambda t: value_grad(t)[0],
        start,
        jac=lambda t: value_grad(t)[1],
        method="L-BFGS-B",
        bounds=[(None, None)] * p + [(-6.0, 12.0)],
        options={"maxiter": int(max_iter), "ftol": 1e-10, "gtol": 1e-6},
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
    loglik, score = _beta_binomial_value_score(theta, y, n, design)
    return {
        "success": bool(result.success),
        "message": str(result.message),
        "theta": theta,
        "bread": bread,
        "score": score,
        "names": [*names, "log_kappa"],
        "log_likelihood": float(np.sum(loglik)),
        "kappa": float(np.exp(np.clip(theta[-1], -6.0, 12.0))),
    }


def _assemble_cluster_covariance(
    fits: list[dict[str, Any]],
    clusters: list[np.ndarray],
) -> tuple[np.ndarray, list[str], np.ndarray]:
    offsets: list[int] = []
    size = 0
    names: list[str] = []
    for fit in fits:
        offsets.append(size)
        names.extend(fit["names"])
        size += len(fit["names"])
    bread = np.zeros((size, size), dtype=float)
    cluster_scores: dict[str, np.ndarray] = {}
    total_rows = 0
    for fit, labels, offset in zip(fits, clusters, offsets, strict=True):
        dim = len(fit["names"])
        bread[offset : offset + dim, offset : offset + dim] = fit["bread"]
        labels = np.asarray(labels).astype(str)
        total_rows += len(labels)
        for cluster in np.unique(labels):
            if cluster not in cluster_scores:
                cluster_scores[cluster] = np.zeros(size, dtype=float)
            cluster_scores[cluster][offset : offset + dim] += fit["score"][
                labels == cluster
            ].sum(axis=0)
    meat = np.zeros((size, size), dtype=float)
    for score in cluster_scores.values():
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    g = len(cluster_scores)
    if g > 1 and total_rows > size:
        covariance *= (g / (g - 1.0)) * ((total_rows - 1.0) / (total_rows - size))
    theta = np.concatenate([fit["theta"] for fit in fits])
    return covariance, names, theta


def _prepare(
    counts: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    geography = str(config["geography_column"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    required = {"island_id", geography, context, cluster, *baseline}
    missing = required - set(covariates.columns)
    if missing:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")
    data = counts.merge(
        covariates[list(required)].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    for column in ["successes", "trials", geography, *baseline]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data[context] = data[context].fillna("").astype(str)
    data[cluster] = data[cluster].fillna("").astype(str)
    data = data.dropna(subset=["successes", "trials", geography, *baseline])
    return data.loc[
        data["trials"].gt(0)
        & data["successes"].ge(0)
        & data["successes"].le(data["trials"])
        & data[context].ne("")
        & data[cluster].ne("")
    ].copy()


def _fit_within(
    data: pd.DataFrame,
    *,
    stratum: str,
    context_value: str,
    threshold: int,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    geography = str(config["geography_column"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    work = data.loc[
        data["stratum"].eq(stratum) & data[context].eq(context_value)
    ].copy()
    counts = work.groupby("outcome")["island_id"].nunique()
    retained = [
        str(x)
        for x in config["model_outcomes"]
        if int(counts.get(str(x), 0)) >= int(threshold)
    ]
    if len(retained) < int(config["minimum_outcomes_per_vector"]):
        return pd.DataFrame(), {
            "stratum": stratum,
            "context": context_value,
            "status": "not_testable",
            "n_retained_outcomes": len(retained),
            "retained_outcomes": "|".join(retained),
        }

    fits: list[dict[str, Any]] = []
    cluster_parts: list[np.ndarray] = []
    slope_local_index: list[int] = []
    outcome_rows: list[dict[str, Any]] = []
    offset = 0
    for outcome in retained:
        part = work.loc[work["outcome"].eq(outcome)].copy()
        columns = [np.ones(len(part), dtype=float)]
        names = [f"{outcome}:intercept"]
        for predictor in baseline:
            columns.append(_standardize(part[predictor]))
            names.append(f"{outcome}:z_{predictor}")
        columns.append(_standardize(part[geography]))
        names.append(f"{outcome}:z_{geography}")
        design = np.column_stack(columns)
        fit = _fit_single_beta_binomial(
            part["successes"].to_numpy(float),
            part["trials"].to_numpy(float),
            design,
            names,
            max_iter=int(config.get("max_iter", 800)),
        )
        fits.append(fit)
        cluster_parts.append(part[cluster].to_numpy(str))
        slope_local_index.append(offset + names.index(f"{outcome}:z_{geography}"))
        outcome_rows.append(
            {
                "stratum": stratum,
                "context": context_value,
                "outcome": outcome,
                "n_islands": int(part["island_id"].nunique()),
                "n_species_trials": int(part["trials"].sum()),
                "kappa": float(fit["kappa"]),
                "optimizer_success": bool(fit["success"]),
            }
        )
        offset += len(fit["names"])

    covariance, _, theta = _assemble_cluster_covariance(fits, cluster_parts)
    slopes = theta[slope_local_index]
    slope_cov = covariance[np.ix_(slope_local_index, slope_local_index)]
    slope_se = np.sqrt(np.clip(np.diag(slope_cov), 0.0, None))
    for row, estimate, stderr in zip(outcome_rows, slopes, slope_se, strict=True):
        z = float(estimate / stderr) if stderr > 0 else float("nan")
        row.update(
            {
                "geography_slope_log_odds": float(estimate),
                "cluster_robust_se": float(stderr),
                "p_value": _normal_two_sided_p(z),
            }
        )
    rank = int(np.linalg.matrix_rank(slope_cov))
    statistic = (
        float(slopes @ np.linalg.pinv(slope_cov) @ slopes) if rank > 0 else float("nan")
    )
    omnibus = {
        "stratum": stratum,
        "context": context_value,
        "status": "fit",
        "n_retained_outcomes": len(retained),
        "retained_outcomes": "|".join(retained),
        "n_unique_islands": int(
            work.loc[work["outcome"].isin(retained), "island_id"].nunique()
        ),
        "n_clusters": int(work.loc[work["outcome"].isin(retained), cluster].nunique()),
        "joint_wald_chisq": statistic,
        "joint_df": rank,
        "p_value": _chi_square_sf_integer_df(statistic, rank),
        "all_optimizers_converged": all(bool(f["success"]) for f in fits),
    }
    return pd.DataFrame(outcome_rows), omnibus


def _fit_between(
    data: pd.DataFrame,
    *,
    stratum: str,
    context_a: str,
    context_b: str,
    threshold: int,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    geography = str(config["geography_column"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    work = data.loc[
        data["stratum"].eq(stratum) & data[context].isin([context_a, context_b])
    ].copy()
    support = work.groupby(["outcome", context])["island_id"].nunique().unstack(fill_value=0)
    for value in (context_a, context_b):
        if value not in support.columns:
            support[value] = 0
    retained = [
        str(x)
        for x in config["model_outcomes"]
        if str(x) in support.index
        if int(support.loc[str(x), context_a]) >= int(threshold)
        and int(support.loc[str(x), context_b]) >= int(threshold)
    ]
    if len(retained) < int(config["minimum_outcomes_per_vector"]):
        return pd.DataFrame(), {
            "stratum": stratum,
            "context_a": context_a,
            "context_b": context_b,
            "status": "not_testable",
            "n_retained_outcomes": len(retained),
            "retained_outcomes": "|".join(retained),
        }

    fits: list[dict[str, Any]] = []
    cluster_parts: list[np.ndarray] = []
    interaction_indices: list[int] = []
    rows: list[dict[str, Any]] = []
    offset = 0
    for outcome in retained:
        part = work.loc[work["outcome"].eq(outcome)].copy()
        context_b_indicator = part[context].eq(context_b).to_numpy(float)
        columns = [np.ones(len(part), dtype=float), context_b_indicator]
        names = [f"{outcome}:intercept", f"{outcome}:context[{context_b}]"]
        for predictor in baseline:
            z = _standardize(part[predictor])
            columns.extend([z, z * context_b_indicator])
            names.extend(
                [
                    f"{outcome}:z_{predictor}",
                    f"{outcome}:z_{predictor}:context[{context_b}]",
                ]
            )
        z_geo = _standardize(part[geography])
        columns.extend([z_geo, z_geo * context_b_indicator])
        interaction_name = f"{outcome}:z_{geography}:context[{context_b}]"
        names.extend([f"{outcome}:z_{geography}", interaction_name])
        fit = _fit_single_beta_binomial(
            part["successes"].to_numpy(float),
            part["trials"].to_numpy(float),
            np.column_stack(columns),
            names,
            max_iter=int(config.get("max_iter", 800)),
        )
        fits.append(fit)
        cluster_parts.append(part[cluster].to_numpy(str))
        interaction_indices.append(offset + names.index(interaction_name))
        rows.append(
            {
                "stratum": stratum,
                "context_a": context_a,
                "context_b": context_b,
                "outcome": outcome,
                "n_islands_context_a": int(support.loc[outcome, context_a]),
                "n_islands_context_b": int(support.loc[outcome, context_b]),
                "kappa": float(fit["kappa"]),
                "optimizer_success": bool(fit["success"]),
            }
        )
        offset += len(fit["names"])

    covariance, _, theta = _assemble_cluster_covariance(fits, cluster_parts)
    vector = theta[interaction_indices]
    vector_cov = covariance[np.ix_(interaction_indices, interaction_indices)]
    vector_se = np.sqrt(np.clip(np.diag(vector_cov), 0.0, None))
    for row, estimate, stderr in zip(rows, vector, vector_se, strict=True):
        z = float(estimate / stderr) if stderr > 0 else float("nan")
        row.update(
            {
                "slope_difference_b_minus_a": float(estimate),
                "cluster_robust_se": float(stderr),
                "p_value": _normal_two_sided_p(z),
            }
        )
    rank = int(np.linalg.matrix_rank(vector_cov))
    statistic = (
        float(vector @ np.linalg.pinv(vector_cov) @ vector)
        if rank > 0
        else float("nan")
    )
    omnibus = {
        "stratum": stratum,
        "context_a": context_a,
        "context_b": context_b,
        "status": "fit",
        "n_retained_outcomes": len(retained),
        "retained_outcomes": "|".join(retained),
        "n_unique_islands": int(
            work.loc[work["outcome"].isin(retained), "island_id"].nunique()
        ),
        "n_clusters": int(work.loc[work["outcome"].isin(retained), cluster].nunique()),
        "joint_wald_chisq": statistic,
        "joint_df": rank,
        "p_value": _chi_square_sf_integer_df(statistic, rank),
        "all_optimizers_converged": all(bool(f["success"]) for f in fits),
    }
    return pd.DataFrame(rows), omnibus


def run_probability_analysis(
    counts: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    data = _prepare(counts, covariates, config)
    contexts = [str(x) for x in config["contexts"]]
    strata = [str(x) for x in config["strata"]]
    threshold = int(config["minimum_islands_per_outcome"])
    within_parts: list[pd.DataFrame] = []
    between_parts: list[pd.DataFrame] = []
    within_rows: list[dict[str, Any]] = []
    between_rows: list[dict[str, Any]] = []
    for stratum in strata:
        for context_value in contexts:
            slopes, result = _fit_within(
                data,
                stratum=stratum,
                context_value=context_value,
                threshold=threshold,
                config=config,
            )
            if not slopes.empty:
                within_parts.append(slopes)
            within_rows.append(result)
        for context_a, context_b in config["between_contexts"]:
            slopes, result = _fit_between(
                data,
                stratum=stratum,
                context_a=str(context_a),
                context_b=str(context_b),
                threshold=threshold,
                config=config,
            )
            if not slopes.empty:
                between_parts.append(slopes)
            between_rows.append(result)
    within = pd.DataFrame(within_rows)
    between = pd.DataFrame(between_rows)
    if "p_value" in within.columns:
        within["q_value"] = within.groupby("stratum", group_keys=False)["p_value"].transform(_bh)
        within["vector_supported"] = within["q_value"].le(float(config["alpha"])).fillna(False)
    if "p_value" in between.columns:
        between["q_value"] = between.groupby("stratum", group_keys=False)["p_value"].transform(_bh)
        between["difference_supported"] = between["q_value"].le(float(config["alpha"])).fillna(False)
    return (
        pd.concat(within_parts, ignore_index=True)
        if within_parts
        else pd.DataFrame(),
        pd.concat(between_parts, ignore_index=True)
        if between_parts
        else pd.DataFrame(),
        within,
        between,
    )


@app.command("run")
def run(
    status_flora_csv: Path = typer.Option(..., exists=True),
    state_audit_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    counts = build_broad_counts(
        pd.read_csv(status_flora_csv),
        pd.read_csv(state_audit_csv),
        config,
    )
    within_slopes, between_slopes, within, between = run_probability_analysis(
        counts,
        pd.read_csv(covariates_csv),
        config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(output_dir / "all_data_probability_counts.csv.gz", index=False)
    within_slopes.to_csv(output_dir / "beta_binomial_within_slopes.csv", index=False)
    between_slopes.to_csv(output_dir / "beta_binomial_between_slopes.csv", index=False)
    within.to_csv(output_dir / "beta_binomial_within_omnibus.csv", index=False)
    between.to_csv(output_dir / "beta_binomial_between_omnibus.csv", index=False)
    manifest = {
        "contract": str(config["contract"]),
        "evidence_scope": str(evidence_scope),
        "model_family": "beta_binomial_logit",
        "primary_evidence_scope": str(config["evidence_roles"]["primary"]),
        "direct_sensitivity_scope": str(config["evidence_roles"]["direct_sensitivity"]),
        "broad_flora_scope": "all_observed",
        "status_sensitivities": ["all_native", "native_nonendemic"],
        "minimum_species_per_island_cutoff": None,
        "claim_boundary": (
            "all_observed estimates describe the observed island flora and may include "
            "introduced or unresolved-status species; native assembly/evolution claims "
            "require status-resolved sensitivities"
        ),
    }
    (output_dir / "all_data_probability_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
