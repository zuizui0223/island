"""Redesigned all-data Chapter 1 H1-H3 analysis.

H1 asks whether one distance-response vector can describe all predeclared contexts.
H2 retains the planned direct North--Tropical contrast from the all-data probability
analysis. H3 tests continuous distance-by-area moderation within contexts and direct
context differences in that moderation.

H4-H5 are intentionally not fitted here because they require native-status/source-pool
or independent mechanism evidence. Their conditional role is frozen in
``config/chapter1_hypothesis_ladder_v2.yml``.
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
    _chi_square_sf_integer_df,
    _fit_single_beta_binomial,
    _normal_two_sided_p,
    _standardize,
    build_broad_counts,
    run_probability_analysis,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _joint_test(vector: np.ndarray, covariance: np.ndarray) -> tuple[float, int, float]:
    vector = np.asarray(vector, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    rank = int(np.linalg.matrix_rank(covariance))
    if rank <= 0:
        return float("nan"), 0, float("nan")
    statistic = float(vector @ np.linalg.pinv(covariance) @ vector)
    return statistic, rank, _chi_square_sf_integer_df(statistic, rank)


def _prepare(
    counts: pd.DataFrame,
    covariates: pd.DataFrame,
    probability_config: dict[str, Any],
    ladder_config: dict[str, Any],
) -> pd.DataFrame:
    geography = str(ladder_config["geography_column"])
    area = str(ladder_config["area_column"])
    context = str(ladder_config["context_column"])
    cluster = str(ladder_config["cluster_column"])
    controls = [str(x) for x in ladder_config["control_columns"]]
    required = {"island_id", geography, area, context, cluster, *controls}
    missing = required - set(covariates.columns)
    if missing:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")
    data = counts.merge(
        covariates[list(required)].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    for column in ["successes", "trials", geography, area, *controls]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data[context] = data[context].fillna("").astype(str)
    data[cluster] = data[cluster].fillna("").astype(str)
    data = data.dropna(subset=["successes", "trials", geography, area, *controls])
    data = data.loc[
        data["trials"].gt(0)
        & data["successes"].ge(0)
        & data["successes"].le(data["trials"])
        & data[context].ne("")
        & data[cluster].ne("")
    ].copy()
    allowed = {str(x) for x in ladder_config["contexts"]}
    data = data.loc[data[context].isin(allowed)].copy()
    outcomes = {str(x) for x in probability_config["model_outcomes"]}
    return data.loc[data["outcome"].astype(str).isin(outcomes)].copy()


def _retained_outcomes_all_contexts(
    work: pd.DataFrame,
    contexts: list[str],
    outcomes: list[str],
    context_column: str,
    threshold: int,
) -> list[str]:
    support = work.groupby(["outcome", context_column])["island_id"].nunique().unstack(fill_value=0)
    for value in contexts:
        if value not in support.columns:
            support[value] = 0
    return [
        outcome
        for outcome in outcomes
        if outcome in support.index
        and all(int(support.loc[outcome, value]) >= threshold for value in contexts)
    ]


def fit_h1_global_context_heterogeneity(
    data: pd.DataFrame,
    probability_config: dict[str, Any],
    ladder_config: dict[str, Any],
    *,
    stratum: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Fit a global H1 test of equality of distance slopes across all contexts."""
    geography = str(ladder_config["geography_column"])
    area = str(ladder_config["area_column"])
    context = str(ladder_config["context_column"])
    cluster = str(ladder_config["cluster_column"])
    controls = [str(x) for x in ladder_config["control_columns"]]
    contexts = [str(x) for x in ladder_config["contexts"]]
    reference = contexts[0]
    alternatives = contexts[1:]
    outcomes = [str(x) for x in probability_config["model_outcomes"]]
    threshold = int(ladder_config["minimum_islands_per_outcome"])
    work = data.loc[data["stratum"].eq(stratum)].copy()
    retained = _retained_outcomes_all_contexts(work, contexts, outcomes, context, threshold)
    if len(retained) < int(ladder_config["minimum_outcomes_per_vector"]):
        return pd.DataFrame(), {
            "hypothesis": "H1_universal_response",
            "stratum": stratum,
            "status": "not_testable",
            "n_retained_outcomes": len(retained),
            "retained_outcomes": "|".join(retained),
        }

    fits: list[dict[str, Any]] = []
    clusters: list[np.ndarray] = []
    target_indices: list[int] = []
    target_meta: list[tuple[str, str]] = []
    offset = 0
    for outcome in retained:
        part = work.loc[work["outcome"].eq(outcome)].copy()
        indicators = {value: part[context].eq(value).to_numpy(float) for value in alternatives}
        columns = [np.ones(len(part), dtype=float)]
        names = [f"{outcome}:intercept"]
        for value in alternatives:
            columns.append(indicators[value])
            names.append(f"{outcome}:context[{value}]")

        for predictor in [area, *controls]:
            z = _standardize(part[predictor])
            columns.append(z)
            names.append(f"{outcome}:z_{predictor}")
            for value in alternatives:
                columns.append(z * indicators[value])
                names.append(f"{outcome}:z_{predictor}:context[{value}]")

        z_distance = _standardize(part[geography])
        columns.append(z_distance)
        names.append(f"{outcome}:z_{geography}")
        for value in alternatives:
            interaction = f"{outcome}:z_{geography}:context[{value}]"
            columns.append(z_distance * indicators[value])
            names.append(interaction)
            target_indices.append(offset + names.index(interaction))
            target_meta.append((outcome, value))

        fit = _fit_single_beta_binomial(
            part["successes"].to_numpy(float),
            part["trials"].to_numpy(float),
            np.column_stack(columns),
            names,
            max_iter=int(probability_config.get("max_iter", 1000)),
        )
        fits.append(fit)
        clusters.append(part[cluster].to_numpy(str))
        offset += len(fit["names"])

    covariance, _, theta = _assemble_cluster_covariance(fits, clusters)
    vector = theta[target_indices]
    cov = covariance[np.ix_(target_indices, target_indices)]
    se = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    rows: list[dict[str, Any]] = []
    for (outcome, value), estimate, stderr in zip(target_meta, vector, se, strict=True):
        z = estimate / stderr if stderr > 0 else float("nan")
        rows.append(
            {
                "hypothesis": "H1_universal_response",
                "stratum": stratum,
                "reference_context": reference,
                "comparison_context": value,
                "outcome": outcome,
                "distance_slope_difference_log_odds": float(estimate),
                "cluster_robust_se": float(stderr),
                "p_value": _normal_two_sided_p(float(z)),
            }
        )
    statistic, df, p_value = _joint_test(vector, cov)
    omnibus = {
        "hypothesis": "H1_universal_response",
        "stratum": stratum,
        "status": "fit",
        "reference_context": reference,
        "contexts": "|".join(contexts),
        "n_retained_outcomes": len(retained),
        "retained_outcomes": "|".join(retained),
        "n_unique_islands": int(work.loc[work["outcome"].isin(retained), "island_id"].nunique()),
        "n_clusters": int(work.loc[work["outcome"].isin(retained), cluster].nunique()),
        "joint_wald_chisq": statistic,
        "joint_df": df,
        "p_value": p_value,
        "universal_response_rejected": bool(math.isfinite(p_value) and p_value < float(ladder_config["alpha"])),
        "all_optimizers_converged": all(bool(fit["success"]) for fit in fits),
    }
    return pd.DataFrame(rows), omnibus


def _fit_h3_within_context(
    data: pd.DataFrame,
    probability_config: dict[str, Any],
    ladder_config: dict[str, Any],
    *,
    stratum: str,
    context_value: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    geography = str(ladder_config["geography_column"])
    area = str(ladder_config["area_column"])
    context = str(ladder_config["context_column"])
    cluster = str(ladder_config["cluster_column"])
    controls = [str(x) for x in ladder_config["control_columns"]]
    threshold = int(ladder_config["minimum_islands_per_outcome"])
    outcomes = [str(x) for x in probability_config["model_outcomes"]]
    work = data.loc[data["stratum"].eq(stratum) & data[context].eq(context_value)].copy()
    support = work.groupby("outcome")["island_id"].nunique()
    retained = [outcome for outcome in outcomes if int(support.get(outcome, 0)) >= threshold]
    if len(retained) < int(ladder_config["minimum_outcomes_per_vector"]):
        return pd.DataFrame(), {
            "hypothesis": "H3_island_capacity_moderation",
            "stratum": stratum,
            "context": context_value,
            "status": "not_testable",
            "n_retained_outcomes": len(retained),
            "retained_outcomes": "|".join(retained),
        }

    fits: list[dict[str, Any]] = []
    clusters: list[np.ndarray] = []
    distance_indices: list[int] = []
    moderation_indices: list[int] = []
    rows: list[dict[str, Any]] = []
    offset = 0
    for outcome in retained:
        part = work.loc[work["outcome"].eq(outcome)].copy()
        z_distance = _standardize(part[geography])
        z_area = _standardize(part[area])
        columns = [np.ones(len(part), dtype=float)]
        names = [f"{outcome}:intercept"]
        for predictor in controls:
            columns.append(_standardize(part[predictor]))
            names.append(f"{outcome}:z_{predictor}")
        columns.extend([z_area, z_distance, z_distance * z_area])
        area_name = f"{outcome}:z_{area}"
        distance_name = f"{outcome}:z_{geography}"
        moderation_name = f"{outcome}:z_{geography}:z_{area}"
        names.extend([area_name, distance_name, moderation_name])
        fit = _fit_single_beta_binomial(
            part["successes"].to_numpy(float),
            part["trials"].to_numpy(float),
            np.column_stack(columns),
            names,
            max_iter=int(probability_config.get("max_iter", 1000)),
        )
        fits.append(fit)
        clusters.append(part[cluster].to_numpy(str))
        distance_indices.append(offset + names.index(distance_name))
        moderation_indices.append(offset + names.index(moderation_name))
        rows.append(
            {
                "hypothesis": "H3_island_capacity_moderation",
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

    covariance, _, theta = _assemble_cluster_covariance(fits, clusters)
    distance_vector = theta[distance_indices]
    moderation_vector = theta[moderation_indices]
    moderation_cov = covariance[np.ix_(moderation_indices, moderation_indices)]
    moderation_se = np.sqrt(np.clip(np.diag(moderation_cov), 0.0, None))
    distance_se = np.sqrt(np.clip(np.diag(covariance[np.ix_(distance_indices, distance_indices)]), 0.0, None))
    for row, d, dse, interaction, ise in zip(
        rows, distance_vector, distance_se, moderation_vector, moderation_se, strict=True
    ):
        z = interaction / ise if ise > 0 else float("nan")
        row.update(
            {
                "distance_slope_log_odds_at_mean_area": float(d),
                "distance_slope_se": float(dse),
                "distance_by_area_log_odds": float(interaction),
                "distance_by_area_se": float(ise),
                "distance_by_area_p_value": _normal_two_sided_p(float(z)),
            }
        )

    statistic, df, p_value = _joint_test(moderation_vector, moderation_cov)
    dot = float(distance_vector @ moderation_vector)
    norm_product = float(np.linalg.norm(distance_vector) * np.linalg.norm(moderation_vector))
    cosine = dot / norm_product if norm_product > 0 else float("nan")
    supported = bool(math.isfinite(p_value) and p_value < float(ladder_config["alpha"]))
    if not supported or not math.isfinite(dot) or abs(dot) < 1e-12:
        classification = "unresolved"
    elif dot < 0:
        classification = "small_island_amplification"
    else:
        classification = "large_island_amplification"
    omnibus = {
        "hypothesis": "H3_island_capacity_moderation",
        "stratum": stratum,
        "context": context_value,
        "status": "fit",
        "n_retained_outcomes": len(retained),
        "retained_outcomes": "|".join(retained),
        "n_unique_islands": int(work.loc[work["outcome"].isin(retained), "island_id"].nunique()),
        "n_clusters": int(work.loc[work["outcome"].isin(retained), cluster].nunique()),
        "joint_wald_chisq": statistic,
        "joint_df": df,
        "p_value": p_value,
        "moderation_supported": supported,
        "distance_moderation_dot_product": dot,
        "distance_moderation_cosine": cosine,
        "directional_classification": classification,
        "all_optimizers_converged": all(bool(fit["success"]) for fit in fits),
    }
    return pd.DataFrame(rows), omnibus


def _fit_h3_between_contexts(
    data: pd.DataFrame,
    probability_config: dict[str, Any],
    ladder_config: dict[str, Any],
    *,
    stratum: str,
    context_a: str,
    context_b: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    geography = str(ladder_config["geography_column"])
    area = str(ladder_config["area_column"])
    context = str(ladder_config["context_column"])
    cluster = str(ladder_config["cluster_column"])
    controls = [str(x) for x in ladder_config["control_columns"]]
    threshold = int(ladder_config["minimum_islands_per_outcome"])
    outcomes = [str(x) for x in probability_config["model_outcomes"]]
    work = data.loc[
        data["stratum"].eq(stratum) & data[context].isin([context_a, context_b])
    ].copy()
    support = work.groupby(["outcome", context])["island_id"].nunique().unstack(fill_value=0)
    for value in (context_a, context_b):
        if value not in support.columns:
            support[value] = 0
    retained = [
        outcome
        for outcome in outcomes
        if outcome in support.index
        and int(support.loc[outcome, context_a]) >= threshold
        and int(support.loc[outcome, context_b]) >= threshold
    ]
    if len(retained) < int(ladder_config["minimum_outcomes_per_vector"]):
        return pd.DataFrame(), {
            "hypothesis": "H3_island_capacity_moderation",
            "stratum": stratum,
            "context_a": context_a,
            "context_b": context_b,
            "status": "not_testable",
            "n_retained_outcomes": len(retained),
            "retained_outcomes": "|".join(retained),
        }

    fits: list[dict[str, Any]] = []
    clusters: list[np.ndarray] = []
    target_indices: list[int] = []
    rows: list[dict[str, Any]] = []
    offset = 0
    for outcome in retained:
        part = work.loc[work["outcome"].eq(outcome)].copy()
        b = part[context].eq(context_b).to_numpy(float)
        z_distance = _standardize(part[geography])
        z_area = _standardize(part[area])
        columns = [np.ones(len(part), dtype=float), b]
        names = [f"{outcome}:intercept", f"{outcome}:context[{context_b}]"]
        for predictor in controls:
            z = _standardize(part[predictor])
            columns.extend([z, z * b])
            names.extend(
                [f"{outcome}:z_{predictor}", f"{outcome}:z_{predictor}:context[{context_b}]"]
            )
        columns.extend([z_area, z_area * b, z_distance, z_distance * b])
        names.extend(
            [
                f"{outcome}:z_{area}",
                f"{outcome}:z_{area}:context[{context_b}]",
                f"{outcome}:z_{geography}",
                f"{outcome}:z_{geography}:context[{context_b}]",
            ]
        )
        da = z_distance * z_area
        triple_name = f"{outcome}:z_{geography}:z_{area}:context[{context_b}]"
        columns.extend([da, da * b])
        names.extend([f"{outcome}:z_{geography}:z_{area}", triple_name])
        fit = _fit_single_beta_binomial(
            part["successes"].to_numpy(float),
            part["trials"].to_numpy(float),
            np.column_stack(columns),
            names,
            max_iter=int(probability_config.get("max_iter", 1000)),
        )
        fits.append(fit)
        clusters.append(part[cluster].to_numpy(str))
        target_indices.append(offset + names.index(triple_name))
        rows.append(
            {
                "hypothesis": "H3_island_capacity_moderation",
                "stratum": stratum,
                "context_a": context_a,
                "context_b": context_b,
                "outcome": outcome,
                "n_islands_context_a": int(support.loc[outcome, context_a]),
                "n_islands_context_b": int(support.loc[outcome, context_b]),
                "optimizer_success": bool(fit["success"]),
            }
        )
        offset += len(fit["names"])

    covariance, _, theta = _assemble_cluster_covariance(fits, clusters)
    vector = theta[target_indices]
    cov = covariance[np.ix_(target_indices, target_indices)]
    se = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    for row, estimate, stderr in zip(rows, vector, se, strict=True):
        z = estimate / stderr if stderr > 0 else float("nan")
        row.update(
            {
                "distance_by_area_difference_b_minus_a": float(estimate),
                "cluster_robust_se": float(stderr),
                "p_value": _normal_two_sided_p(float(z)),
            }
        )
    statistic, df, p_value = _joint_test(vector, cov)
    omnibus = {
        "hypothesis": "H3_island_capacity_moderation",
        "stratum": stratum,
        "context_a": context_a,
        "context_b": context_b,
        "status": "fit",
        "n_retained_outcomes": len(retained),
        "retained_outcomes": "|".join(retained),
        "n_unique_islands": int(work.loc[work["outcome"].isin(retained), "island_id"].nunique()),
        "n_clusters": int(work.loc[work["outcome"].isin(retained), cluster].nunique()),
        "joint_wald_chisq": statistic,
        "joint_df": df,
        "p_value": p_value,
        "moderation_difference_supported": bool(
            math.isfinite(p_value) and p_value < float(ladder_config["alpha"])
        ),
        "all_optimizers_converged": all(bool(fit["success"]) for fit in fits),
    }
    return pd.DataFrame(rows), omnibus


def run_ladder(
    counts: pd.DataFrame,
    covariates: pd.DataFrame,
    probability_config: dict[str, Any],
    ladder_config: dict[str, Any],
) -> dict[str, pd.DataFrame]:
    data = _prepare(counts, covariates, probability_config, ladder_config)
    strata = [str(x) for x in probability_config["strata"]]

    h1_slopes: list[pd.DataFrame] = []
    h1_rows: list[dict[str, Any]] = []
    for stratum in strata:
        slopes, omnibus = fit_h1_global_context_heterogeneity(
            data, probability_config, ladder_config, stratum=stratum
        )
        if not slopes.empty:
            h1_slopes.append(slopes)
        h1_rows.append(omnibus)

    h2_within_slopes, h2_between_slopes, h2_within, h2_between = run_probability_analysis(
        counts, covariates, probability_config
    )

    h3_within_slopes: list[pd.DataFrame] = []
    h3_within_rows: list[dict[str, Any]] = []
    h3_between_slopes: list[pd.DataFrame] = []
    h3_between_rows: list[dict[str, Any]] = []
    for stratum in strata:
        for context_value in [str(x) for x in ladder_config["contexts"]]:
            slopes, omnibus = _fit_h3_within_context(
                data,
                probability_config,
                ladder_config,
                stratum=stratum,
                context_value=context_value,
            )
            if not slopes.empty:
                h3_within_slopes.append(slopes)
            h3_within_rows.append(omnibus)
        for context_a, context_b in ladder_config["primary_between_contexts"]:
            slopes, omnibus = _fit_h3_between_contexts(
                data,
                probability_config,
                ladder_config,
                stratum=stratum,
                context_a=str(context_a),
                context_b=str(context_b),
            )
            if not slopes.empty:
                h3_between_slopes.append(slopes)
            h3_between_rows.append(omnibus)

    return {
        "h1_global_slopes": pd.concat(h1_slopes, ignore_index=True) if h1_slopes else pd.DataFrame(),
        "h1_global_omnibus": pd.DataFrame(h1_rows),
        "h2_within_slopes": h2_within_slopes,
        "h2_between_slopes": h2_between_slopes,
        "h2_within_omnibus": h2_within,
        "h2_between_omnibus": h2_between,
        "h3_within_slopes": pd.concat(h3_within_slopes, ignore_index=True) if h3_within_slopes else pd.DataFrame(),
        "h3_within_omnibus": pd.DataFrame(h3_within_rows),
        "h3_between_slopes": pd.concat(h3_between_slopes, ignore_index=True) if h3_between_slopes else pd.DataFrame(),
        "h3_between_omnibus": pd.DataFrame(h3_between_rows),
    }


@app.command("run")
def run(
    status_flora_csv: Path = typer.Option(..., exists=True),
    state_audit_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    probability_config_path: Path = typer.Option(..., exists=True),
    ladder_config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
) -> None:
    probability_config = yaml.safe_load(probability_config_path.read_text(encoding="utf-8"))
    ladder_config = yaml.safe_load(ladder_config_path.read_text(encoding="utf-8"))
    status_flora = pd.read_csv(status_flora_csv)
    state_audit = pd.read_csv(state_audit_csv)
    covariates = pd.read_csv(covariates_csv)
    counts = build_broad_counts(status_flora, state_audit, probability_config)
    results = run_ladder(counts, covariates, probability_config, ladder_config)
    output_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(output_dir / "all_data_ladder_counts.csv.gz", index=False)
    for name, frame in results.items():
        if not frame.empty:
            frame.insert(0, "evidence_scope", evidence_scope)
        frame.to_csv(output_dir / f"{name}.csv", index=False)
    manifest = {
        "contract": ladder_config["contract"],
        "evidence_scope": evidence_scope,
        "model_family": "beta_binomial_logit",
        "H1": "global context heterogeneity in distance-response vector",
        "H2": "planned direct biogeographic branching contrasts",
        "H3": "continuous distance-by-area moderation",
        "H4": "conditional native hierarchical assembly; not fitted by this command",
        "H5": "conditional independent mechanism gate; not fitted by this command",
    }
    (output_dir / "ladder_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
