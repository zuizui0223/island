"""Primary PR138 test: observed biogeographic structure of floral island syndromes.

The primary response is the observed broad trait composition of the status-resolved
island flora. Genus-fixed residuals are deliberately secondary: they decompose an
observed assemblage pattern into lineage-composition versus residual components, but
do not decide whether the assemblage-level pattern exists.

No pollinator predictor enters the primary model. Bombus, butterfly, bird and other
pollination syndromes are interpreted only after the biogeographic trait vector has
been estimated.
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

from island_v2.chapter1_context_analysis import (
    _chi_square_sf_integer_df,
    _fit_grouped_binomial_design,
)
from island_v2.flora_status_support import stratum_mask

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
    adjusted = np.clip(adjusted, 0.0, 1.0)
    restored = np.empty(n, dtype=float)
    restored[order] = adjusted
    out.loc[ok] = restored
    return out


def _one_sided_positive_p(z: float) -> float:
    if not math.isfinite(z):
        return float("nan")
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def _classify_signature(signature: object, positive: set[str], negative: set[str]) -> float:
    tokens = {token.strip() for token in str(signature or "").split("|") if token.strip()}
    if not tokens:
        return float("nan")
    if tokens <= positive:
        return 1.0
    if tokens <= negative:
        return 0.0
    return float("nan")


def build_observed_broad_counts(
    status_flora: pd.DataFrame,
    state_audit: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build island x stratum x broad-outcome direct-evidence binomial counts.

    Species whose direct multistate signature crosses the positive and negative broad
    sets are left unresolved for that binary outcome rather than forced to either side.
    """
    required_flora = {
        "island_id",
        "accepted_species",
        "origin_status",
        "endemic_status",
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
        trait_name = str(spec["trait_name"])
        positive = {str(x) for x in spec["positive_states"]}
        negative = {str(x) for x in spec["negative_states"]}
        states = audit.loc[
            audit["trait_name"].eq(trait_name),
            ["accepted_species", "canonical_signature"],
        ].drop_duplicates("accepted_species")
        states = states.copy()
        states["state"] = [
            _classify_signature(value, positive, negative)
            for value in states["canonical_signature"]
        ]
        states = states.dropna(subset=["state"])[["accepted_species", "state"]]
        joined = flora.merge(states, on="accepted_species", how="inner", validate="many_to_one")
        for stratum in [str(x) for x in config["strata"]]:
            subset = joined.loc[stratum_mask(joined, stratum), ["island_id", "accepted_species", "state"]]
            subset = subset.drop_duplicates(["island_id", "accepted_species"])
            if subset.empty:
                continue
            summary = (
                subset.groupby("island_id", as_index=False)
                .agg(successes=("state", "sum"), trials=("state", "size"))
            )
            summary["share"] = summary["successes"] / summary["trials"]
            summary["outcome"] = str(outcome)
            summary["stratum"] = stratum
            rows.append(summary)
    if not rows:
        return pd.DataFrame(
            columns=["island_id", "successes", "trials", "share", "outcome", "stratum"]
        )
    return pd.concat(rows, ignore_index=True)


def _standardized_masked(
    work: pd.DataFrame,
    mask: np.ndarray,
    column: str,
) -> np.ndarray:
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


def _directional_projection(
    slopes: np.ndarray,
    covariance: np.ndarray,
    directions: np.ndarray,
) -> tuple[float, float, float, float]:
    weights = directions / float(len(directions))
    estimate = float(weights @ slopes)
    variance = float(weights @ covariance @ weights)
    se = math.sqrt(max(variance, 0.0))
    z = estimate / se if se > 0 else float("nan")
    return estimate, se, z, _one_sided_positive_p(z)


def _prepare_counts(
    counts: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    geography = str(config["geography_column"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    required_counts = {"island_id", "outcome", "stratum", "successes", "trials"}
    required_cov = {"island_id", geography, context, cluster, *baseline}
    missing = required_counts - set(counts.columns)
    if missing:
        raise typer.BadParameter(f"broad counts missing columns: {sorted(missing)}")
    missing = required_cov - set(covariates.columns)
    if missing:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")
    data = counts.merge(
        covariates[["island_id", geography, context, cluster, *baseline]].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    numeric = ["successes", "trials", geography, *baseline]
    for column in numeric:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data[context] = data[context].fillna("").astype(str)
    data[cluster] = data[cluster].fillna("").astype(str)
    data = data.dropna(subset=numeric)
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
        data["stratum"].astype(str).eq(stratum)
        & data[context].eq(context_value)
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
    coefficients, fit, covariance = _fit_grouped_binomial_design(
        work["successes"].to_numpy(float),
        work["trials"].to_numpy(float),
        np.column_stack(columns),
        names,
        work[cluster].to_numpy(str),
    )
    if coefficients.empty:
        return coefficients, {
            "stratum": stratum,
            "support_tier": support_tier,
            "threshold": threshold,
            "context": context_value,
            "status": str(fit.get("status", "fit_failed")),
            "n_retained_outcomes": len(retained),
            "retained_outcomes": "|".join(retained),
        }
    beta = coefficients.set_index("predictor")["estimate_log_odds"]
    indices = [names.index(name) for name in geography_names]
    slopes = np.array([float(beta[name]) for name in geography_names])
    cov = covariance[np.ix_(indices, indices)]
    joint_stat, joint_df, joint_p = _joint_wald(slopes, cov)
    directions = np.array([
        float(config["broad_outcomes"][outcome]["northern_classic_direction"])
        for outcome in retained
    ])
    proj, proj_se, proj_z, proj_p = _directional_projection(slopes, cov, directions)
    slope_rows = []
    coef_indexed = coefficients.set_index("predictor")
    for outcome, name in zip(retained, geography_names, strict=True):
        row = coef_indexed.loc[name]
        slope_rows.append({
            "stratum": stratum,
            "support_tier": support_tier,
            "context": context_value,
            "outcome": outcome,
            "geography_slope_log_odds_per_response_sd": float(row["estimate_log_odds"]),
            "cluster_robust_se": float(row["cluster_robust_se"]),
            "p_value": float(row["p_value"]),
            "n_islands": int(counts.loc[outcome]),
        })
    return pd.DataFrame(slope_rows), {
        "stratum": stratum,
        "support_tier": support_tier,
        "threshold": threshold,
        "context": context_value,
        "status": "fit",
        "n_retained_outcomes": len(retained),
        "retained_outcomes": "|".join(retained),
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "joint_wald_chisq": joint_stat,
        "joint_df": joint_df,
        "p_value": joint_p,
        "classic_projection_estimate": proj,
        "classic_projection_se": proj_se,
        "classic_projection_z": proj_z,
        "classic_projection_one_sided_p": proj_p,
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
        data["stratum"].astype(str).eq(stratum)
        & data[context].isin([context_a, context_b])
    ].copy()
    support = work.groupby(["outcome", context])["island_id"].nunique().unstack(fill_value=0)
    for value in (context_a, context_b):
        if value not in support.columns:
            support[value] = 0
    retained = sorted(
        support.index[
            support[context_a].ge(threshold) & support[context_b].ge(threshold)
        ].astype(str)
    )
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
    coefficients, fit, covariance = _fit_grouped_binomial_design(
        work["successes"].to_numpy(float),
        work["trials"].to_numpy(float),
        np.column_stack(columns),
        names,
        work[cluster].to_numpy(str),
    )
    if coefficients.empty:
        return coefficients, {
            "stratum": stratum,
            "support_tier": support_tier,
            "threshold": threshold,
            "context_a": context_a,
            "context_b": context_b,
            "status": str(fit.get("status", "fit_failed")),
            "n_retained_outcomes": len(retained),
            "retained_outcomes": "|".join(retained),
        }
    beta = coefficients.set_index("predictor")["estimate_log_odds"]
    indices = [names.index(name) for name in interaction_names]
    vector = np.array([float(beta[name]) for name in interaction_names])
    cov = covariance[np.ix_(indices, indices)]
    joint_stat, joint_df, joint_p = _joint_wald(vector, cov)
    rows = []
    coef_indexed = coefficients.set_index("predictor")
    for outcome, name in zip(retained, interaction_names, strict=True):
        row = coef_indexed.loc[name]
        rows.append({
            "stratum": stratum,
            "support_tier": support_tier,
            "context_a": context_a,
            "context_b": context_b,
            "outcome": outcome,
            "slope_difference_b_minus_a": float(row["estimate_log_odds"]),
            "cluster_robust_se": float(row["cluster_robust_se"]),
            "p_value": float(row["p_value"]),
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
        "joint_wald_chisq": joint_stat,
        "joint_df": joint_df,
        "p_value": joint_p,
    }


def run_observed_pattern(
    broad_counts: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    data = _prepare_counts(broad_counts, covariates, config)
    contexts = [str(x) for x in config["contexts"]]
    strata = [str(x) for x in config["strata"]]
    tiers = {str(k): int(v) for k, v in config["support_tiers"].items()}
    within_slope_parts: list[pd.DataFrame] = []
    between_slope_parts: list[pd.DataFrame] = []
    within_rows: list[dict[str, Any]] = []
    between_rows: list[dict[str, Any]] = []
    for stratum in strata:
        for tier, threshold in tiers.items():
            for context_value in contexts:
                slopes, result = _fit_within(
                    data,
                    stratum=stratum,
                    context_value=context_value,
                    support_tier=tier,
                    threshold=threshold,
                    config=config,
                )
                if not slopes.empty:
                    within_slope_parts.append(slopes)
                within_rows.append(result)
            for context_a, context_b in itertools.combinations(contexts, 2):
                slopes, result = _fit_between(
                    data,
                    stratum=stratum,
                    context_a=context_a,
                    context_b=context_b,
                    support_tier=tier,
                    threshold=threshold,
                    config=config,
                )
                if not slopes.empty:
                    between_slope_parts.append(slopes)
                between_rows.append(result)
    within = pd.DataFrame(within_rows)
    between = pd.DataFrame(between_rows)
    if "p_value" in within.columns:
        within["q_joint_within_stratum_tier"] = within.groupby(
            ["stratum", "support_tier"], group_keys=False
        )["p_value"].transform(_bh)
        within["observed_vector_supported"] = within["q_joint_within_stratum_tier"].le(
            float(config["alpha"])
        )
    if "classic_projection_one_sided_p" in within.columns:
        within["q_classic_projection_within_stratum_tier"] = within.groupby(
            ["stratum", "support_tier"], group_keys=False
        )["classic_projection_one_sided_p"].transform(_bh)
        within["classic_direction_supported"] = within[
            "q_classic_projection_within_stratum_tier"
        ].le(float(config["alpha"]))
    if "p_value" in between.columns:
        between["q_joint_within_stratum_tier"] = between.groupby(
            ["stratum", "support_tier"], group_keys=False
        )["p_value"].transform(_bh)
        between["biogeographic_difference_supported"] = between[
            "q_joint_within_stratum_tier"
        ].le(float(config["alpha"]))
    return (
        pd.concat(within_slope_parts, ignore_index=True) if within_slope_parts else pd.DataFrame(),
        pd.concat(between_slope_parts, ignore_index=True) if between_slope_parts else pd.DataFrame(),
        within,
        between,
    )


@app.command("run")
def run(
    status_flora_csv: Path = typer.Option(..., exists=True),
    state_audit_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(
        Path("config/chapter1_pr136_biogeographic_pattern.yml"), exists=True
    ),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    counts = build_observed_broad_counts(
        pd.read_csv(status_flora_csv),
        pd.read_csv(state_audit_csv),
        config,
    )
    within_slopes, between_slopes, within, between = run_observed_pattern(
        counts,
        pd.read_csv(covariates_csv),
        config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(output_dir / "observed_broad_syndrome_counts.csv.gz", index=False)
    within_slopes.to_csv(output_dir / "observed_within_outcome_slopes.csv", index=False)
    between_slopes.to_csv(output_dir / "observed_between_outcome_slope_differences.csv", index=False)
    within.to_csv(output_dir / "observed_within_region_omnibus.csv", index=False)
    between.to_csv(output_dir / "observed_between_region_omnibus.csv", index=False)
    manifest = {
        "contract": str(config["contract"]),
        "primary_response": "observed status-stratified broad trait composition",
        "broad_outcomes": list(config["broad_outcomes"]),
        "pollinator_predictors_in_primary_model": False,
        "secondary_layer": "genus-fixed residual decomposition",
        "discussion_layer": "pollination-syndrome concordance only",
    }
    (output_dir / "observed_biogeographic_pattern_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
