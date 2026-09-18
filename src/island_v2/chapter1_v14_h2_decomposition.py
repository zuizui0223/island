"""Chapter 1 v14 H2: selfing versus pollinator-facing floral decomposition.

H2 asks whether the floral portion of the island syndrome is fully reducible to
measured reproductive assurance.  It deliberately separates:
- H2a: selfing_core ~ isolation + controls
- H2b: plain_colour / generalized_accessible ~ isolation + selfing_core + controls

The conditional H2b coefficient is decomposition, not causal mediation or proof
of direct selection by a named pollinator.
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

from island_v2.chapter1_all_data_probability import (
    _assemble_cluster_covariance,
    _fit_single_beta_binomial,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _normal_two_sided_p(z: float) -> float:
    return math.erfc(abs(float(z)) / math.sqrt(2.0)) if math.isfinite(z) else float("nan")


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


def _z(series: pd.Series) -> np.ndarray:
    x = pd.to_numeric(series, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != "chapter1_v14_h2_decomposition_v1":
        raise typer.BadParameter("unexpected v14 H2 decomposition contract")
    return config


def build_plant_scores(scores: pd.DataFrame, stratum: str) -> pd.DataFrame:
    required = {"island_id", "stratum", "syndrome", "syndrome_score"}
    if missing := required - set(scores.columns):
        raise typer.BadParameter(f"syndrome score table missing columns: {sorted(missing)}")
    part = scores.loc[scores["stratum"].astype(str).eq(stratum)].copy()
    wide = part.pivot_table(
        index="island_id", columns="syndrome", values="syndrome_score", aggfunc="first"
    )
    needed = {"selfing_core", "generalized_accessible"}
    if missing := needed - set(wide.columns):
        raise typer.BadParameter(f"syndrome score table lacks required axes: {sorted(missing)}")
    return wide.reset_index()


def _clustered_ols(
    data: pd.DataFrame,
    *,
    response: str,
    predictors: list[str],
    cluster_column: str,
) -> dict[str, Any]:
    required = [response, *predictors, cluster_column]
    work = data[required].copy()
    for column in [response, *predictors]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work[cluster_column] = work[cluster_column].fillna("").astype(str)
    work = work.dropna(subset=[response, *predictors])
    work = work.loc[work[cluster_column].ne("")].copy()
    if len(work) < 20:
        return {"status": "not_testable", "n_islands": int(len(work))}
    names = ["intercept"]
    columns = [np.ones(len(work), dtype=float)]
    for predictor in predictors:
        try:
            columns.append(_z(work[predictor]))
        except ValueError:
            continue
        names.append(f"z_{predictor}")
    if len(columns) < 2:
        return {"status": "not_testable", "n_islands": int(len(work))}
    X = np.column_stack(columns)
    y = pd.to_numeric(work[response], errors="coerce").to_numpy(float)
    bread = np.linalg.pinv(X.T @ X)
    beta = bread @ X.T @ y
    residual = y - X @ beta
    labels = work[cluster_column].astype(str).to_numpy()
    unique = np.unique(labels)
    meat = np.zeros((X.shape[1], X.shape[1]), dtype=float)
    for label in unique:
        mask = labels == label
        score = X[mask].T @ residual[mask]
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    n = len(work)
    k = X.shape[1]
    g = len(unique)
    if g > 1 and n > k:
        covariance *= (g / (g - 1.0)) * ((n - 1.0) / (n - k))
    coefficients: dict[str, dict[str, float]] = {}
    for index, name in enumerate(names):
        se = float(math.sqrt(max(float(covariance[index, index]), 0.0)))
        estimate = float(beta[index])
        z_value = estimate / se if se > 0 else float("nan")
        coefficients[name] = {
            "estimate": estimate,
            "se": se,
            "p": _normal_two_sided_p(z_value),
        }
    return {
        "status": "fit",
        "n_islands": int(n),
        "n_clusters": int(g),
        "coefficients": coefficients,
    }


def _fit_plain_colour(
    counts: pd.DataFrame,
    plant_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    context_value: str,
    config: dict[str, Any],
) -> dict[str, Any]:
    geography = str(config["geography_column"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    stratum = str(config["primary_stratum"])
    part = counts.loc[
        counts["stratum"].astype(str).eq(stratum)
        & counts["outcome"].astype(str).eq("plain_colour")
    ].copy()
    needed_cov = ["island_id", geography, context, cluster, *baseline]
    work = (
        part.merge(plant_scores[["island_id", "selfing_core"]], on="island_id", how="left")
        .merge(covariates[needed_cov].drop_duplicates("island_id"), on="island_id", how="left")
    )
    numeric = ["successes", "trials", geography, "selfing_core", *baseline]
    for column in numeric:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.loc[work[context].astype(str).eq(context_value)].dropna(
        subset=[*numeric, cluster]
    )
    work = work.loc[work["trials"].gt(0)]
    if work["island_id"].nunique() < int(config["minimum_islands"]):
        return {"status": "not_testable", "n_islands": int(work["island_id"].nunique())}

    predictor_order = [*baseline, "selfing_core", geography]
    names = ["intercept"]
    columns = [np.ones(len(work), dtype=float)]
    for predictor in predictor_order:
        columns.append(_z(work[predictor]))
        names.append(f"z_{predictor}")
    fit = _fit_single_beta_binomial(
        work["successes"].to_numpy(float),
        work["trials"].to_numpy(float),
        np.column_stack(columns),
        names,
        max_iter=1000,
    )
    covariance, assembled_names, theta = _assemble_cluster_covariance(
        [fit], [work[cluster].astype(str).to_numpy()]
    )
    coefficients: dict[str, dict[str, float]] = {}
    for predictor in (geography, "selfing_core"):
        name = f"z_{predictor}"
        index = assembled_names.index(name)
        se = float(math.sqrt(max(float(covariance[index, index]), 0.0)))
        estimate = float(theta[index])
        z_value = estimate / se if se > 0 else float("nan")
        coefficients[predictor] = {
            "estimate": estimate,
            "se": se,
            "p": _normal_two_sided_p(z_value),
        }
    return {
        "status": "fit",
        "n_islands": int(work["island_id"].nunique()),
        "n_clusters": int(work[cluster].nunique()),
        "kappa": float(fit["kappa"]),
        "optimizer_success": bool(fit["success"]),
        "coefficients": coefficients,
    }


def run_decomposition(
    counts: pd.DataFrame,
    syndrome_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
    *,
    evidence_scope: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    geography = str(config["geography_column"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    contexts = [str(x) for x in config["contexts"]]
    plant = build_plant_scores(syndrome_scores, str(config["primary_stratum"]))
    needed_cov = ["island_id", geography, context, cluster, *baseline]
    continuous = plant.merge(
        covariates[needed_cov].drop_duplicates("island_id"), on="island_id", how="left"
    )
    numeric_columns = [
        geography,
        *baseline,
        "selfing_core",
        "generalized_accessible",
    ]
    for column in numeric_columns:
        continuous[column] = pd.to_numeric(continuous[column], errors="coerce")
    continuous[context] = continuous[context].fillna("").astype(str)
    continuous[cluster] = continuous[cluster].fillna("").astype(str)

    rows: list[dict[str, Any]] = []
    for context_value in contexts:
        part = continuous.loc[continuous[context].eq(context_value)].copy()
        for response, conditional, role in (
            ("selfing_core", False, "H2a_reproductive_assurance"),
            ("generalized_accessible", True, "H2b_accessibility_conditional_on_selfing"),
        ):
            predictors = [geography, *baseline]
            if conditional:
                predictors.insert(1, "selfing_core")
            result = _clustered_ols(
                part, response=response, predictors=predictors, cluster_column=cluster
            )
            row: dict[str, Any] = {
                "evidence_scope": evidence_scope,
                "context": context_value,
                "response": response,
                "analysis_role": role,
                "model_family": "clustered_ols",
                "status": result["status"],
                "n_islands": int(result.get("n_islands", 0)),
                "n_clusters": int(result.get("n_clusters", 0)),
            }
            if result["status"] == "fit":
                distance = result["coefficients"][f"z_{geography}"]
                row.update(
                    distance_estimate=distance["estimate"],
                    distance_se=distance["se"],
                    distance_p=distance["p"],
                )
                if conditional and "z_selfing_core" in result["coefficients"]:
                    selfing = result["coefficients"]["z_selfing_core"]
                    row.update(
                        selfing_core_estimate=selfing["estimate"],
                        selfing_core_se=selfing["se"],
                        selfing_core_p=selfing["p"],
                    )
            rows.append(row)

        colour = _fit_plain_colour(
            counts, plant, covariates, context_value=context_value, config=config
        )
        row = {
            "evidence_scope": evidence_scope,
            "context": context_value,
            "response": "plain_colour",
            "analysis_role": "H2b_colour_dulling_conditional_on_selfing",
            "model_family": "beta_binomial_logit",
            "status": colour["status"],
            "n_islands": int(colour.get("n_islands", 0)),
            "n_clusters": int(colour.get("n_clusters", 0)),
        }
        if colour["status"] == "fit":
            distance = colour["coefficients"][geography]
            selfing = colour["coefficients"]["selfing_core"]
            row.update(
                distance_estimate=distance["estimate"],
                distance_se=distance["se"],
                distance_p=distance["p"],
                selfing_core_estimate=selfing["estimate"],
                selfing_core_se=selfing["se"],
                selfing_core_p=selfing["p"],
                kappa=colour["kappa"],
                optimizer_success=colour["optimizer_success"],
            )
        rows.append(row)

    results = pd.DataFrame(rows)
    results["primary_H2b_q"] = np.nan
    primary = (
        results["status"].eq("fit")
        & results["response"].isin(["plain_colour", "generalized_accessible"])
    )
    results.loc[primary, "primary_H2b_q"] = _bh(results.loc[primary, "distance_p"])

    summary = {
        "contract": str(config["contract"]),
        "evidence_scope": evidence_scope,
        "interpretation": (
            "Conditional H2b persistence means the floral response is not reducible "
            "to measured selfing_core; it is not causal mediation or direct-selection proof."
        ),
        "n_contexts_selfing_core_positive": int(
            results.loc[
                results["response"].eq("selfing_core") & results["status"].eq("fit"),
                "distance_estimate",
            ].gt(0).sum()
        ),
        "n_contexts_plain_colour_positive_after_selfing": int(
            results.loc[
                results["response"].eq("plain_colour") & results["status"].eq("fit"),
                "distance_estimate",
            ].gt(0).sum()
        ),
        "n_contexts_accessibility_positive_after_selfing": int(
            results.loc[
                results["response"].eq("generalized_accessible") & results["status"].eq("fit"),
                "distance_estimate",
            ].gt(0).sum()
        ),
    }
    return results, summary


@app.command("run")
def run(
    counts_csv: Path = typer.Option(..., exists=True),
    syndrome_scores_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    results, summary = run_decomposition(
        pd.read_csv(counts_csv),
        pd.read_csv(syndrome_scores_csv),
        pd.read_csv(covariates_csv),
        config,
        evidence_scope=evidence_scope,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "h2_decomposition_models.csv", index=False)
    (output_dir / "h2_decomposition_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(summary, indent=2))


if __name__ == "__main__":
    app()
