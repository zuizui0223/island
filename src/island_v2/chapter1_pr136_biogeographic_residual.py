"""PR136 primary pattern test: biogeographic contingency of genus-fixed residual traits.

The Chapter 1 primary analysis deliberately does *not* put Bombus or any other
pollinator into the model. It first establishes the pattern: whether the classic
floral island syndrome strengthens along the mainland-distance gradient in northern
mid-latitude island floras, and whether tropical and southern island floras show a
weaker, absent, heterogeneous, or different residual trait vector after floristic
status and genus composition are represented.

Pollination syndromes are interpreted only after the residual vector is established.
Bombus absence, butterfly mobility, bird pollination, Baker's law, and reduced need
for attraction are candidate explanations, not primary predictors.

A significant result in one region and a nonsignificant result in another is not
evidence of regional heterogeneity; the between-context vector contrast is tested
directly.
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


def _normal_two_sided_p(z: float) -> float:
    return math.erfc(abs(float(z)) / math.sqrt(2.0))


def _standardize(values: pd.Series) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


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
        return (
            pd.DataFrame(),
            np.empty((0, 0)),
            {"status": "insufficient_complete_rows", "n_rows": n, "n_clusters": n_clusters},
        )

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
    rows: list[dict[str, Any]] = []
    for name, estimate, stderr in zip(names, beta, se, strict=True):
        z = float(estimate / stderr) if stderr > 0 else float("nan")
        rows.append(
            {
                "predictor": name,
                "estimate": float(estimate),
                "cluster_robust_se": float(stderr),
                "z_value": z,
                "p_value": _normal_two_sided_p(z) if math.isfinite(z) else float("nan"),
            }
        )
    return (
        pd.DataFrame(rows),
        covariance,
        {"status": "fit", "n_rows": n, "n_clusters": n_clusters},
    )


def _prepare(
    genus_null: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    geography = str(config["geography_column"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    required_null = {
        "island_id",
        "outcome",
        "stratum",
        "observed_n_species",
        "deviation_observed_minus_null",
    }
    required_cov = {"island_id", geography, context, cluster, *baseline}
    missing = required_null - set(genus_null.columns)
    if missing:
        raise typer.BadParameter(f"genus-null table missing columns: {sorted(missing)}")
    missing = required_cov - set(covariates.columns)
    if missing:
        raise typer.BadParameter(f"covariate table missing columns: {sorted(missing)}")

    data = genus_null.merge(
        covariates[list(required_cov)].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    numeric = [
        "observed_n_species",
        "deviation_observed_minus_null",
        geography,
        *baseline,
    ]
    for column in numeric:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data[context] = data[context].fillna("").astype(str)
    data[cluster] = data[cluster].fillna("").astype(str)
    data = data.dropna(subset=numeric)
    data = data.loc[
        data["observed_n_species"].gt(0)
        & data[context].ne("")
        & data[cluster].ne("")
    ].copy()
    allowed = {str(x) for x in config["contexts"]}
    return data.loc[data[context].isin(allowed)].copy()


def _within_context(
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
    minimum_outcomes = int(config.get("minimum_outcomes_per_vector", 2))

    work = data.loc[
        data["stratum"].astype(str).eq(stratum)
        & data[context].eq(context_value)
    ].copy()
    counts = work.groupby("outcome")["island_id"].nunique()
    outcomes = sorted(counts.loc[counts.ge(threshold)].index.astype(str))
    work = work.loc[work["outcome"].astype(str).isin(outcomes)].copy()
    base_result = {
        "stratum": stratum,
        "support_tier": support_tier,
        "threshold": threshold,
        "context": context_value,
        "n_retained_outcomes": len(outcomes),
        "retained_outcomes": "|".join(outcomes),
    }
    if len(outcomes) < minimum_outcomes:
        return pd.DataFrame(), {**base_result, "status": "not_testable"}

    names: list[str] = []
    columns: list[np.ndarray] = []
    geography_names: list[str] = []
    for outcome in outcomes:
        mask = work["outcome"].astype(str).eq(outcome).to_numpy()
        indicator = mask.astype(float)
        names.append(f"outcome[{outcome}]")
        columns.append(indicator)
        for predictor in baseline:
            vector = np.zeros(len(work), dtype=float)
            vector[mask] = _standardize(work.loc[mask, predictor])
            names.append(f"outcome[{outcome}]:z_{predictor}")
            columns.append(vector)
        vector = np.zeros(len(work), dtype=float)
        vector[mask] = _standardize(work.loc[mask, geography])
        name = f"outcome[{outcome}]:z_{geography}"
        names.append(name)
        columns.append(vector)
        geography_names.append(name)

    coefficients, covariance, fit = _fit_weighted_clustered_design(
        work["deviation_observed_minus_null"].to_numpy(float),
        work["observed_n_species"].to_numpy(float),
        np.column_stack(columns),
        names,
        work[cluster].to_numpy(str),
    )
    if coefficients.empty:
        return coefficients, {**base_result, "status": fit["status"]}
    beta = coefficients.set_index("predictor")["estimate"]
    indices = [names.index(x) for x in geography_names]
    vector = np.array([float(beta[x]) for x in geography_names])
    cov = covariance[np.ix_(indices, indices)]
    statistic, df, p_value = _joint_wald(vector, cov)
    coefficients.insert(0, "context", context_value)
    coefficients.insert(0, "support_tier", support_tier)
    coefficients.insert(0, "stratum", stratum)
    return coefficients, {
        **base_result,
        "status": "fit",
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "joint_residual_vector_chisq": statistic,
        "joint_residual_vector_df": df,
        "p_value": p_value,
    }


def _between_contexts(
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
    minimum_outcomes = int(config.get("minimum_outcomes_per_vector", 2))

    work = data.loc[
        data["stratum"].astype(str).eq(stratum)
        & data[context].isin([context_a, context_b])
    ].copy()
    counts = work.groupby(["outcome", context])["island_id"].nunique().unstack(fill_value=0)
    if context_a not in counts.columns or context_b not in counts.columns:
        return pd.DataFrame(), {
            "stratum": stratum,
            "support_tier": support_tier,
            "threshold": threshold,
            "context_a": context_a,
            "context_b": context_b,
            "status": "not_testable",
            "n_retained_outcomes": 0,
            "retained_outcomes": "",
        }
    outcomes = sorted(
        counts.index[
            counts[context_a].ge(threshold) & counts[context_b].ge(threshold)
        ].astype(str)
    )
    base_result = {
        "stratum": stratum,
        "support_tier": support_tier,
        "threshold": threshold,
        "context_a": context_a,
        "context_b": context_b,
        "n_retained_outcomes": len(outcomes),
        "retained_outcomes": "|".join(outcomes),
    }
    if len(outcomes) < minimum_outcomes:
        return pd.DataFrame(), {**base_result, "status": "not_testable"}

    work = work.loc[work["outcome"].astype(str).isin(outcomes)].copy()
    context_b_indicator = work[context].eq(context_b).to_numpy(float)
    names: list[str] = []
    columns: list[np.ndarray] = []
    interaction_names: list[str] = []
    for outcome in outcomes:
        mask = work["outcome"].astype(str).eq(outcome).to_numpy()
        indicator = mask.astype(float)
        names.append(f"outcome[{outcome}]")
        columns.append(indicator)
        for predictor in baseline:
            vector = np.zeros(len(work), dtype=float)
            vector[mask] = _standardize(work.loc[mask, predictor])
            names.append(f"outcome[{outcome}]:z_{predictor}")
            columns.append(vector)
        names.append(f"outcome[{outcome}]:context[{context_b}]")
        columns.append(indicator * context_b_indicator)
        geography_vector = np.zeros(len(work), dtype=float)
        geography_vector[mask] = _standardize(work.loc[mask, geography])
        names.append(f"outcome[{outcome}]:z_{geography}")
        columns.append(geography_vector)
        interaction = f"outcome[{outcome}]:z_{geography}:context[{context_b}]"
        names.append(interaction)
        columns.append(geography_vector * context_b_indicator)
        interaction_names.append(interaction)

    coefficients, covariance, fit = _fit_weighted_clustered_design(
        work["deviation_observed_minus_null"].to_numpy(float),
        work["observed_n_species"].to_numpy(float),
        np.column_stack(columns),
        names,
        work[cluster].to_numpy(str),
    )
    if coefficients.empty:
        return coefficients, {**base_result, "status": fit["status"]}
    beta = coefficients.set_index("predictor")["estimate"]
    indices = [names.index(x) for x in interaction_names]
    vector = np.array([float(beta[x]) for x in interaction_names])
    cov = covariance[np.ix_(indices, indices)]
    statistic, df, p_value = _joint_wald(vector, cov)
    coefficients.insert(0, "context_b", context_b)
    coefficients.insert(0, "context_a", context_a)
    coefficients.insert(0, "support_tier", support_tier)
    coefficients.insert(0, "stratum", stratum)
    return coefficients, {
        **base_result,
        "status": "fit",
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "joint_context_difference_chisq": statistic,
        "joint_context_difference_df": df,
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
    within_coef: list[pd.DataFrame] = []
    between_coef: list[pd.DataFrame] = []
    within_rows: list[dict[str, Any]] = []
    between_rows: list[dict[str, Any]] = []

    for stratum in strata:
        for tier, threshold in tiers.items():
            for context_value in contexts:
                coefficients, result = _within_context(
                    data,
                    stratum=stratum,
                    context_value=context_value,
                    support_tier=tier,
                    threshold=threshold,
                    config=config,
                )
                if not coefficients.empty:
                    within_coef.append(coefficients)
                within_rows.append(result)
            for context_a, context_b in itertools.combinations(contexts, 2):
                coefficients, result = _between_contexts(
                    data,
                    stratum=stratum,
                    context_a=context_a,
                    context_b=context_b,
                    support_tier=tier,
                    threshold=threshold,
                    config=config,
                )
                if not coefficients.empty:
                    between_coef.append(coefficients)
                between_rows.append(result)

    within = pd.DataFrame(within_rows)
    between = pd.DataFrame(between_rows)
    if not within.empty:
        within["q_within_stratum_tier"] = within.groupby(
            ["stratum", "support_tier"], group_keys=False
        )["p_value"].apply(_bh)
        within["residual_filtering_supported"] = within[
            "q_within_stratum_tier"
        ].le(float(config.get("alpha", 0.05)))
    if not between.empty:
        between["q_within_stratum_tier"] = between.groupby(
            ["stratum", "support_tier"], group_keys=False
        )["p_value"].apply(_bh)
        between["biogeographic_difference_supported"] = between[
            "q_within_stratum_tier"
        ].le(float(config.get("alpha", 0.05)))

    return (
        pd.concat(within_coef, ignore_index=True) if within_coef else pd.DataFrame(),
        pd.concat(between_coef, ignore_index=True) if between_coef else pd.DataFrame(),
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
    within_coef, between_coef, within, between = run_biogeographic_residual(
        pd.read_csv(genus_null_csv), pd.read_csv(covariates_csv), config
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    within_coef.to_csv(output_dir / "biogeographic_residual_within_coefficients.csv", index=False)
    between_coef.to_csv(output_dir / "biogeographic_residual_between_coefficients.csv", index=False)
    within.to_csv(output_dir / "biogeographic_residual_within_omnibus.csv", index=False)
    between.to_csv(output_dir / "biogeographic_residual_between_omnibus.csv", index=False)
    manifest = {
        "contract": "chapter1_pr136_biogeographic_residual_v2",
        "primary_response": "observed direct trait share minus same-genus null expectation",
        "primary_modifier": config["context_column"],
        "pollinator_predictors_in_primary_model": False,
        "pollination_syndrome_role": "post-analysis concordance and discussion only",
        "invalid_inference": "significant in one region and nonsignificant in another is not a regional difference",
    }
    (output_dir / "biogeographic_residual_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
