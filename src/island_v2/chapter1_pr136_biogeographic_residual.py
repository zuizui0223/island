"""SECONDARY PR138 analysis: genus-fixed decomposition of observed syndrome patterns.

This module is not the primary test of whether a floral island syndrome exists. The
primary layer is the observed status-stratified assemblage composition in
``chapter1_pr138_biogeographic_pattern.py``.

Here the response is observed direct trait share minus the same-genus null expectation.
The purpose is to ask whether an observed biogeographic pattern persists beyond genus
composition, or whether it is compatible with floristic / lineage assembly. Absence of
a residual signal does not erase an observed assemblage-level pattern.

No pollinator predictor enters this model. Pollination-syndrome concordance is a later
Discussion-level interpretation.
"""

from __future__ import annotations

import itertools
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_context_analysis import _chi_square_sf_integer_df

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _normal_two_sided_p(z: float) -> float:
    return math.erfc(abs(float(z)) / math.sqrt(2.0))


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
    adjusted = np.clip(adjusted, 0.0, 1.0)
    restored = np.empty(n, dtype=float)
    restored[order] = adjusted
    out.loc[ok] = restored
    return out


def _standardized_masked(work: pd.DataFrame, mask: np.ndarray, column: str) -> np.ndarray:
    x = pd.to_numeric(work.loc[mask, column], errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError(f"constant or invalid predictor: {column}")
    out = np.zeros(len(work), dtype=float)
    out[mask] = (x - mean) / sd
    return out


def _joint_wald(beta: np.ndarray, covariance: np.ndarray) -> tuple[float, int, float]:
    rank = int(np.linalg.matrix_rank(covariance))
    if rank <= 0:
        return float("nan"), 0, float("nan")
    statistic = float(beta @ np.linalg.pinv(covariance) @ beta)
    return statistic, rank, _chi_square_sf_integer_df(statistic, rank)


def _fit_weighted_clustered_design(
    y: np.ndarray,
    weights: np.ndarray,
    design: np.ndarray,
    names: list[str],
    clusters: np.ndarray,
) -> tuple[pd.DataFrame, np.ndarray, dict[str, Any]]:
    n, p = design.shape
    unique_clusters = np.unique(clusters)
    n_clusters = int(len(unique_clusters))
    if n < max(10, p + 3) or n_clusters < 2:
        return pd.DataFrame(), np.empty((0, 0)), {
            "status": "insufficient_complete_rows",
            "n_rows": n,
            "n_clusters": n_clusters,
        }
    xtwx = design.T @ (weights[:, None] * design)
    bread = np.linalg.pinv(xtwx)
    beta = bread @ (design.T @ (weights * y))
    residual = y - design @ beta
    meat = np.zeros((p, p), dtype=float)
    for cluster in unique_clusters:
        mask = clusters == cluster
        score = design[mask].T @ (weights[mask] * residual[mask])
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    if n_clusters > 1 and n > p:
        covariance *= (n_clusters / (n_clusters - 1.0)) * ((n - 1.0) / (n - p))
    se = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
    rows = []
    for name, estimate, stderr in zip(names, beta, se, strict=True):
        z = float(estimate / stderr) if stderr > 0 else float("nan")
        rows.append({
            "predictor": name,
            "estimate": float(estimate),
            "cluster_robust_se": float(stderr),
            "z_value": z,
            "p_value": _normal_two_sided_p(z) if math.isfinite(z) else float("nan"),
        })
    return pd.DataFrame(rows), covariance, {
        "status": "fit",
        "n_rows": n,
        "n_clusters": n_clusters,
    }


def _prepare(genus_null: pd.DataFrame, covariates: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    geography = str(config["geography_column"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    required_null = {
        "island_id", "outcome", "stratum", "observed_n_species",
        "deviation_observed_minus_null",
    }
    required_cov = {"island_id", geography, context, cluster, *baseline}
    missing = required_null - set(genus_null.columns)
    if missing:
        raise typer.BadParameter(f"genus-null table missing columns: {sorted(missing)}")
    missing = required_cov - set(covariates.columns)
    if missing:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")
    data = genus_null.merge(
        covariates[["island_id", geography, context, cluster, *baseline]].drop_duplicates("island_id"),
        on="island_id", how="left", validate="many_to_one"
    )
    numeric = ["observed_n_species", "deviation_observed_minus_null", geography, *baseline]
    for column in numeric:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data[context] = data[context].fillna("").astype(str)
    data[cluster] = data[cluster].fillna("").astype(str)
    data = data.dropna(subset=numeric)
    return data.loc[
        data["observed_n_species"].gt(0) & data[context].ne("") & data[cluster].ne("")
    ].copy()


def _fit_within(
    data: pd.DataFrame,
    *,
    stratum: str,
    context_value: str,
    support_tier: str,
    threshold: int,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    geography = str(config["geography_column"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    min_outcomes = int(config["minimum_outcomes_per_vector"])
    work = data.loc[
        data["stratum"].astype(str).eq(stratum) & data[context].eq(context_value)
    ].copy()
    counts = work.groupby("outcome")["island_id"].nunique()
    retained = sorted(counts.loc[counts.ge(threshold)].index.astype(str))
    if len(retained) < min_outcomes:
        return pd.DataFrame(), {
            "stratum": stratum,
            "support_tier": support_tier,
            "threshold": threshold,
            "context": context_value,
            "status": "not_testable",
            "n_retained_outcomes": len(retained),
            "retained_outcomes": "|".join(retained),
        }
    work = work.loc[work["outcome"].astype(str).isin(retained)].copy()
    names: list[str] = []
    columns: list[np.ndarray] = []
    geography_names: list[str] = []
    for outcome in retained:
        mask = work["outcome"].astype(str).eq(outcome).to_numpy()
        indicator = mask.astype(float)
        names.append(f"outcome[{outcome}]")
        columns.append(indicator)
        for predictor in baseline:
            z = _standardized_masked(work, mask, predictor)
            names.append(f"outcome[{outcome}]:z_{predictor}")
            columns.append(z)
        z_geo = _standardized_masked(work, mask, geography)
        geo_name = f"outcome[{outcome}]:z_{geography}"
        names.append(geo_name)
        columns.append(z_geo)
        geography_names.append(geo_name)
    coef, covariance, fit = _fit_weighted_clustered_design(
        work["deviation_observed_minus_null"].to_numpy(float),
        work["observed_n_species"].to_numpy(float),
        np.column_stack(columns), names, work[cluster].to_numpy(str)
    )
    if coef.empty:
        return coef, {
            "stratum": stratum,
            "support_tier": support_tier,
            "threshold": threshold,
            "context": context_value,
            "status": fit["status"],
            "n_retained_outcomes": len(retained),
            "retained_outcomes": "|".join(retained),
        }
    beta = coef.set_index("predictor")["estimate"]
    indices = [names.index(name) for name in geography_names]
    slopes = np.array([float(beta[name]) for name in geography_names])
    cov = covariance[np.ix_(indices, indices)]
    stat, df, p_value = _joint_wald(slopes, cov)
    rows = []
    ci = coef.set_index("predictor")
    for outcome, name in zip(retained, geography_names, strict=True):
        r = ci.loc[name]
        rows.append({
            "stratum": stratum,
            "support_tier": support_tier,
            "context": context_value,
            "outcome": outcome,
            "geography_slope": float(r["estimate"]),
            "cluster_robust_se": float(r["cluster_robust_se"]),
            "p_value": float(r["p_value"]),
            "n_islands": int(counts.loc[outcome]),
        })
    return pd.DataFrame(rows), {
        "stratum": stratum,
        "support_tier": support_tier,
        "threshold": threshold,
        "context": context_value,
        "status": "fit",
        "n_retained_outcomes": len(retained),
        "retained_outcomes": "|".join(retained),
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "joint_wald_chisq": stat,
        "joint_df": df,
        "p_value": p_value,
    }


def _fit_between(
    data: pd.DataFrame,
    *,
    stratum: str,
    context_a: str,
    context_b: str,
    support_tier: str,
    threshold: int,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    geography = str(config["geography_column"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    min_outcomes = int(config["minimum_outcomes_per_vector"])
    work = data.loc[
        data["stratum"].astype(str).eq(stratum) & data[context].isin([context_a, context_b])
    ].copy()
    support = work.groupby(["outcome", context])["island_id"].nunique().unstack(fill_value=0)
    for value in (context_a, context_b):
        if value not in support.columns:
            support[value] = 0
    retained = sorted(support.index[
        support[context_a].ge(threshold) & support[context_b].ge(threshold)
    ].astype(str))
    if len(retained) < min_outcomes:
        return pd.DataFrame(), {
            "stratum": stratum,
            "support_tier": support_tier,
            "threshold": threshold,
            "context_a": context_a,
            "context_b": context_b,
            "status": "not_testable",
            "n_retained_outcomes": len(retained),
            "retained_outcomes": "|".join(retained),
        }
    work = work.loc[work["outcome"].astype(str).isin(retained)].copy()
    b_indicator = work[context].eq(context_b).to_numpy(float)
    names: list[str] = []
    columns: list[np.ndarray] = []
    interaction_names: list[str] = []
    for outcome in retained:
        mask = work["outcome"].astype(str).eq(outcome).to_numpy()
        indicator = mask.astype(float)
        names.append(f"outcome[{outcome}]")
        columns.append(indicator)
        names.append(f"outcome[{outcome}]:context[{context_b}]")
        columns.append(indicator * b_indicator)
        for predictor in baseline:
            z = _standardized_masked(work, mask, predictor)
            base_name = f"outcome[{outcome}]:z_{predictor}"
            names.append(base_name)
            columns.append(z)
            names.append(f"{base_name}:context[{context_b}]")
            columns.append(z * b_indicator)
        z_geo = _standardized_masked(work, mask, geography)
        geo_name = f"outcome[{outcome}]:z_{geography}"
        names.append(geo_name)
        columns.append(z_geo)
        interaction_name = f"{geo_name}:context[{context_b}]"
        names.append(interaction_name)
        columns.append(z_geo * b_indicator)
        interaction_names.append(interaction_name)
    coef, covariance, fit = _fit_weighted_clustered_design(
        work["deviation_observed_minus_null"].to_numpy(float),
        work["observed_n_species"].to_numpy(float),
        np.column_stack(columns), names, work[cluster].to_numpy(str)
    )
    if coef.empty:
        return coef, {
            "stratum": stratum,
            "support_tier": support_tier,
            "threshold": threshold,
            "context_a": context_a,
            "context_b": context_b,
            "status": fit["status"],
            "n_retained_outcomes": len(retained),
            "retained_outcomes": "|".join(retained),
        }
    beta = coef.set_index("predictor")["estimate"]
    indices = [names.index(name) for name in interaction_names]
    vector = np.array([float(beta[name]) for name in interaction_names])
    cov = covariance[np.ix_(indices, indices)]
    stat, df, p_value = _joint_wald(vector, cov)
    rows = []
    ci = coef.set_index("predictor")
    for outcome, name in zip(retained, interaction_names, strict=True):
        r = ci.loc[name]
        rows.append({
            "stratum": stratum,
            "support_tier": support_tier,
            "context_a": context_a,
            "context_b": context_b,
            "outcome": outcome,
            "slope_difference_b_minus_a": float(r["estimate"]),
            "cluster_robust_se": float(r["cluster_robust_se"]),
            "p_value": float(r["p_value"]),
            "n_islands_context_a": int(support.loc[outcome, context_a]),
            "n_islands_context_b": int(support.loc[outcome, context_b]),
        })
    return pd.DataFrame(rows), {
        "stratum": stratum,
        "support_tier": support_tier,
        "threshold": threshold,
        "context_a": context_a,
        "context_b": context_b,
        "status": "fit",
        "n_retained_outcomes": len(retained),
        "retained_outcomes": "|".join(retained),
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "joint_wald_chisq": stat,
        "joint_df": df,
        "p_value": p_value,
    }


def run_biogeographic_residual(
    genus_null: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    data = _prepare(genus_null, covariates, config)
    contexts = [str(x) for x in config["contexts"]]
    strata = [str(x) for x in config["strata"]]
    tiers = {str(k): int(v) for k, v in config["support_tiers"].items()}
    within_slope_parts = []
    between_slope_parts = []
    within_rows = []
    between_rows = []
    for stratum in strata:
        for tier, threshold in tiers.items():
            for context_value in contexts:
                slopes, result = _fit_within(
                    data, stratum=stratum, context_value=context_value,
                    support_tier=tier, threshold=threshold, config=config
                )
                if not slopes.empty:
                    within_slope_parts.append(slopes)
                within_rows.append(result)
            for context_a, context_b in itertools.combinations(contexts, 2):
                slopes, result = _fit_between(
                    data, stratum=stratum, context_a=context_a, context_b=context_b,
                    support_tier=tier, threshold=threshold, config=config
                )
                if not slopes.empty:
                    between_slope_parts.append(slopes)
                between_rows.append(result)
    within = pd.DataFrame(within_rows)
    between = pd.DataFrame(between_rows)
    if "p_value" in within.columns:
        within["q_within_stratum_tier"] = within.groupby(
            ["stratum", "support_tier"], group_keys=False
        )["p_value"].transform(_bh)
        within["residual_filtering_supported"] = within["q_within_stratum_tier"].le(
            float(config["alpha"])
        )
    if "p_value" in between.columns:
        between["q_within_stratum_tier"] = between.groupby(
            ["stratum", "support_tier"], group_keys=False
        )["p_value"].transform(_bh)
        between["biogeographic_difference_supported"] = between[
            "q_within_stratum_tier"
        ].le(float(config["alpha"]))
    return (
        pd.concat(within_slope_parts, ignore_index=True) if within_slope_parts else pd.DataFrame(),
        pd.concat(between_slope_parts, ignore_index=True) if between_slope_parts else pd.DataFrame(),
        within,
        between,
    )


@app.command("run")
def run(
    genus_null_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(
        Path("config/chapter1_pr136_biogeographic_residual.yml"), exists=True
    ),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    within_slopes, between_slopes, within, between = run_biogeographic_residual(
        pd.read_csv(genus_null_csv), pd.read_csv(covariates_csv), config
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    within_slopes.to_csv(output_dir / "biogeographic_residual_within_slopes.csv", index=False)
    between_slopes.to_csv(output_dir / "biogeographic_residual_between_slopes.csv", index=False)
    within.to_csv(output_dir / "biogeographic_residual_within_omnibus.csv", index=False)
    between.to_csv(output_dir / "biogeographic_residual_between_omnibus.csv", index=False)
    manifest = {
        "contract": str(config["contract"]),
        "role": "secondary genus-fixed decomposition",
        "response": "observed direct trait share minus same-genus null expectation",
        "pollinator_predictors_in_primary_model": False,
        "interpretation": "Residual null cannot erase an observed assemblage-level pattern.",
    }
    (output_dir / "biogeographic_residual_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
