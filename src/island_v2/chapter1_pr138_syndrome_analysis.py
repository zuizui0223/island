"""PR138 multivariate pollination-syndrome concordance analysis.

This module treats pollination syndromes as predeclared multivariate *trait-concordance
responses*. It never converts a floral phenotype into an observed pollinator identity.
The same machinery also scores the selfing syndrome separately, allowing attraction /
pollination-channel shifts to be distinguished from reproductive-assurance shifts.
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
from island_v2.chapter1_pr136_biogeographic_residual import _fit_weighted_clustered_design
from island_v2.flora_status_support import stratum_mask

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _normal_two_sided_p(z: float) -> float:
    return math.erfc(abs(float(z)) / math.sqrt(2.0)) if math.isfinite(z) else float("nan")


def _one_sided_positive_p(z: float) -> float:
    return 0.5 * math.erfc(float(z) / math.sqrt(2.0)) if math.isfinite(z) else float("nan")


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


def _joint_wald(beta: np.ndarray, covariance: np.ndarray) -> tuple[float, int, float]:
    rank = int(np.linalg.matrix_rank(covariance))
    if rank <= 0:
        return float("nan"), 0, float("nan")
    stat = float(beta @ np.linalg.pinv(covariance) @ beta)
    return stat, rank, _chi_square_sf_integer_df(stat, rank)


def _standardize(series: pd.Series) -> np.ndarray:
    x = pd.to_numeric(series, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def _tokens(value: object) -> set[str]:
    text = str(value or "").strip()
    if not text:
        return set()
    if text.startswith("[") and text.endswith("]"):
        try:
            parsed = json.loads(text)
            if isinstance(parsed, list):
                return {str(x).strip() for x in parsed if str(x).strip()}
        except json.JSONDecodeError:
            pass
    return {x.strip() for x in text.split("|") if x.strip()}


def _trait_concordance(value: object, preferred: set[str], opposed: set[str]) -> float:
    states = _tokens(value)
    if not states:
        return float("nan")
    has_preferred = bool(states & preferred)
    has_opposed = bool(states & opposed)
    if has_preferred and has_opposed:
        return float("nan")
    if states <= preferred:
        return 1.0
    if states <= opposed:
        return -1.0
    return float("nan")


def score_species_syndromes(
    trait_ledger: pd.DataFrame,
    syndrome_config: dict[str, Any],
) -> pd.DataFrame:
    required = {"accepted_species", "trait_name", "normalized_value"}
    missing = required - set(trait_ledger.columns)
    if missing:
        raise typer.BadParameter(f"trait ledger missing columns: {sorted(missing)}")
    ledger = trait_ledger.copy()
    ledger["accepted_species"] = ledger["accepted_species"].fillna("").astype(str)
    ledger["trait_name"] = ledger["trait_name"].fillna("").astype(str)
    ledger["normalized_value"] = ledger["normalized_value"].fillna("").astype(str)
    ledger = ledger.drop_duplicates(["accepted_species", "trait_name"], keep="first")
    lookup = ledger.set_index(["accepted_species", "trait_name"])["normalized_value"]
    species = sorted(ledger["accepted_species"].loc[ledger["accepted_species"].ne("")].unique())
    minimum = int(syndrome_config["score_definition"]["minimum_informative_traits"])

    rows: list[dict[str, Any]] = []
    for sp in species:
        for syndrome, spec in syndrome_config["syndromes"].items():
            weighted = 0.0
            weight_total = 0.0
            n_informative = 0
            used: list[str] = []
            for trait_name, trait_spec in spec["traits"].items():
                key = (sp, str(trait_name))
                if key not in lookup.index:
                    continue
                score = _trait_concordance(
                    lookup.loc[key],
                    {str(x) for x in trait_spec.get("preferred", [])},
                    {str(x) for x in trait_spec.get("opposed", [])},
                )
                if not math.isfinite(score):
                    continue
                weight = float(trait_spec.get("weight", 1.0))
                weighted += weight * score
                weight_total += weight
                n_informative += 1
                used.append(str(trait_name))
            if n_informative < minimum or weight_total <= 0:
                continue
            concordance = weighted / weight_total
            rows.append(
                {
                    "accepted_species": sp,
                    "syndrome": str(syndrome),
                    "syndrome_concordance": concordance,
                    "soft_membership": (concordance + 1.0) / 2.0,
                    "n_informative_traits": n_informative,
                    "informative_traits": "|".join(sorted(used)),
                }
            )
    return pd.DataFrame(rows)


def build_island_syndrome_scores(
    status_flora: pd.DataFrame,
    species_scores: pd.DataFrame,
    strata: list[str],
) -> pd.DataFrame:
    required = {"island_id", "accepted_species", "origin_status", "endemic_status", "floristic_status"}
    missing = required - set(status_flora.columns)
    if missing:
        raise typer.BadParameter(f"status flora missing columns: {sorted(missing)}")
    flora = status_flora.copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    merged = flora.merge(species_scores, on="accepted_species", how="inner", validate="many_to_many")
    rows: list[pd.DataFrame] = []
    for stratum in strata:
        subset = merged.loc[stratum_mask(merged, stratum)].copy()
        subset = subset.drop_duplicates(["island_id", "accepted_species", "syndrome"])
        if subset.empty:
            continue
        out = (
            subset.groupby(["island_id", "syndrome"], as_index=False)
            .agg(
                syndrome_score=("syndrome_concordance", "mean"),
                n_species=("accepted_species", "nunique"),
                mean_informative_traits=("n_informative_traits", "mean"),
            )
        )
        out["stratum"] = stratum
        rows.append(out)
    if not rows:
        return pd.DataFrame(columns=["island_id", "syndrome", "syndrome_score", "n_species", "stratum"])
    return pd.concat(rows, ignore_index=True)


def _prepare(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    pattern_config: dict[str, Any],
) -> pd.DataFrame:
    geography = str(pattern_config["geography_column"])
    context = str(pattern_config["context_column"])
    cluster = str(pattern_config["cluster_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    required = {"island_id", geography, context, cluster, *baseline}
    missing = required - set(covariates.columns)
    if missing:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")
    data = island_scores.merge(
        covariates[list(required)].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    for column in ["syndrome_score", "n_species", geography, *baseline]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data[context] = data[context].fillna("").astype(str)
    data[cluster] = data[cluster].fillna("").astype(str)
    data = data.dropna(subset=["syndrome_score", geography, *baseline])
    return data.loc[data[context].ne("") & data[cluster].ne("")].copy()


def _within_context(
    data: pd.DataFrame,
    *,
    stratum: str,
    context_value: str,
    support_tier: str,
    threshold: int,
    pattern_config: dict[str, Any],
    syndrome_config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    geography = str(pattern_config["geography_column"])
    context = str(pattern_config["context_column"])
    cluster = str(pattern_config["cluster_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    work = data.loc[data["stratum"].eq(stratum) & data[context].eq(context_value)].copy()
    counts = work.groupby("syndrome")["island_id"].nunique()
    syndromes = sorted(counts.loc[counts.ge(threshold)].index.astype(str))
    base = {
        "stratum": stratum,
        "support_tier": support_tier,
        "threshold": threshold,
        "context": context_value,
        "n_retained_syndromes": len(syndromes),
        "retained_syndromes": "|".join(syndromes),
    }
    if len(syndromes) < 2:
        return pd.DataFrame(), {**base, "status": "not_testable"}
    work = work.loc[work["syndrome"].isin(syndromes)].copy()
    names: list[str] = []
    columns: list[np.ndarray] = []
    slope_names: list[str] = []
    for syndrome in syndromes:
        mask = work["syndrome"].eq(syndrome).to_numpy()
        indicator = mask.astype(float)
        names.append(f"syndrome[{syndrome}]")
        columns.append(indicator)
        for predictor in baseline:
            z = np.zeros(len(work))
            z[mask] = _standardize(work.loc[mask, predictor])
            names.append(f"syndrome[{syndrome}]:z_{predictor}")
            columns.append(z)
        z_geo = np.zeros(len(work))
        z_geo[mask] = _standardize(work.loc[mask, geography])
        slope_name = f"syndrome[{syndrome}]:z_{geography}"
        names.append(slope_name)
        columns.append(z_geo)
        slope_names.append(slope_name)
    coefficients, covariance, fit = _fit_weighted_clustered_design(
        work["syndrome_score"].to_numpy(float),
        np.ones(len(work), dtype=float),
        np.column_stack(columns),
        names,
        work[cluster].to_numpy(str),
    )
    if coefficients.empty:
        return coefficients, {**base, "status": str(fit.get("status", "fit_failed"))}
    beta = coefficients.set_index("predictor")["estimate"]
    indices = [names.index(x) for x in slope_names]
    slopes = np.array([float(beta[x]) for x in slope_names])
    slope_cov = covariance[np.ix_(indices, indices)]
    stat, df, p_value = _joint_wald(slopes, slope_cov)

    slope_rows: list[dict[str, Any]] = []
    coef_index = coefficients.set_index("predictor")
    for syndrome, slope_name in zip(syndromes, slope_names, strict=True):
        row = coef_index.loc[slope_name]
        slope_rows.append(
            {
                "stratum": stratum,
                "support_tier": support_tier,
                "context": context_value,
                "syndrome": syndrome,
                "distance_slope": float(row["estimate"]),
                "cluster_robust_se": float(row["cluster_robust_se"]),
                "p_value": float(row["p_value"]),
                "n_islands": int(counts.loc[syndrome]),
            }
        )

    classic_names = ["large_bee_like", "generalized_accessible", "selfing_syndrome"]
    classic_projection = float("nan")
    classic_se = float("nan")
    classic_z = float("nan")
    classic_p = float("nan")
    if all(name in syndromes for name in classic_names):
        positions = [syndromes.index(x) for x in classic_names]
        classic_slopes = slopes[positions]
        classic_cov = slope_cov[np.ix_(positions, positions)]
        direction = np.array([-1.0, 1.0, 1.0]) / 3.0
        classic_projection = float(direction @ classic_slopes)
        classic_var = float(direction @ classic_cov @ direction)
        classic_se = math.sqrt(max(classic_var, 0.0))
        classic_z = classic_projection / classic_se if classic_se > 0 else float("nan")
        classic_p = _one_sided_positive_p(classic_z)

    return pd.DataFrame(slope_rows), {
        **base,
        "status": "fit",
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "joint_syndrome_vector_chisq": stat,
        "joint_syndrome_vector_df": df,
        "p_value": p_value,
        "northern_classic_projection": classic_projection,
        "northern_classic_projection_se": classic_se,
        "northern_classic_projection_z": classic_z,
        "northern_classic_projection_one_sided_p": classic_p,
    }


def _between_contexts(
    data: pd.DataFrame,
    *,
    stratum: str,
    context_a: str,
    context_b: str,
    support_tier: str,
    threshold: int,
    pattern_config: dict[str, Any],
) -> dict[str, Any]:
    geography = str(pattern_config["geography_column"])
    context = str(pattern_config["context_column"])
    cluster = str(pattern_config["cluster_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    work = data.loc[data["stratum"].eq(stratum) & data[context].isin([context_a, context_b])].copy()
    counts = work.groupby(["syndrome", context])["island_id"].nunique().unstack(fill_value=0)
    for value in [context_a, context_b]:
        if value not in counts.columns:
            counts[value] = 0
    syndromes = sorted(counts.index[(counts[context_a] >= threshold) & (counts[context_b] >= threshold)].astype(str))
    base = {
        "stratum": stratum,
        "support_tier": support_tier,
        "threshold": threshold,
        "context_a": context_a,
        "context_b": context_b,
        "n_retained_syndromes": len(syndromes),
        "retained_syndromes": "|".join(syndromes),
    }
    if len(syndromes) < 2:
        return {**base, "status": "not_testable"}
    work = work.loc[work["syndrome"].isin(syndromes)].copy()
    b_indicator = work[context].eq(context_b).to_numpy(float)
    names: list[str] = []
    columns: list[np.ndarray] = []
    interaction_names: list[str] = []
    for syndrome in syndromes:
        mask = work["syndrome"].eq(syndrome).to_numpy()
        indicator = mask.astype(float)
        names.append(f"syndrome[{syndrome}]")
        columns.append(indicator)
        for predictor in baseline:
            z = np.zeros(len(work))
            z[mask] = _standardize(work.loc[mask, predictor])
            names.append(f"syndrome[{syndrome}]:z_{predictor}")
            columns.append(z)
        names.append(f"syndrome[{syndrome}]:context[{context_b}]")
        columns.append(indicator * b_indicator)
        z_geo = np.zeros(len(work))
        z_geo[mask] = _standardize(work.loc[mask, geography])
        names.append(f"syndrome[{syndrome}]:z_{geography}")
        columns.append(z_geo)
        interaction = f"syndrome[{syndrome}]:z_{geography}:context[{context_b}]"
        names.append(interaction)
        columns.append(z_geo * b_indicator)
        interaction_names.append(interaction)
    coefficients, covariance, fit = _fit_weighted_clustered_design(
        work["syndrome_score"].to_numpy(float),
        np.ones(len(work), dtype=float),
        np.column_stack(columns),
        names,
        work[cluster].to_numpy(str),
    )
    if coefficients.empty:
        return {**base, "status": str(fit.get("status", "fit_failed"))}
    beta = coefficients.set_index("predictor")["estimate"]
    indices = [names.index(x) for x in interaction_names]
    vector = np.array([float(beta[x]) for x in interaction_names])
    cov = covariance[np.ix_(indices, indices)]
    stat, df, p_value = _joint_wald(vector, cov)
    return {
        **base,
        "status": "fit",
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "joint_context_difference_chisq": stat,
        "joint_context_difference_df": df,
        "p_value": p_value,
    }


def run_syndrome_analysis(
    trait_ledger: pd.DataFrame,
    status_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    pattern_config: dict[str, Any],
    syndrome_config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    species_scores = score_species_syndromes(trait_ledger, syndrome_config)
    strata = [str(x) for x in pattern_config["strata"]]
    island_scores = build_island_syndrome_scores(status_flora, species_scores, strata)
    data = _prepare(island_scores, covariates, pattern_config)
    contexts = [str(x) for x in pattern_config["contexts"]]
    tiers = {str(k): int(v) for k, v in pattern_config["support_tiers"].items()}

    slopes_parts: list[pd.DataFrame] = []
    within_rows: list[dict[str, Any]] = []
    between_rows: list[dict[str, Any]] = []
    for stratum in strata:
        for tier, threshold in tiers.items():
            for context_value in contexts:
                slopes, result = _within_context(
                    data,
                    stratum=stratum,
                    context_value=context_value,
                    support_tier=tier,
                    threshold=threshold,
                    pattern_config=pattern_config,
                    syndrome_config=syndrome_config,
                )
                if not slopes.empty:
                    slopes_parts.append(slopes)
                within_rows.append(result)
            for context_a, context_b in itertools.combinations(contexts, 2):
                between_rows.append(
                    _between_contexts(
                        data,
                        stratum=stratum,
                        context_a=context_a,
                        context_b=context_b,
                        support_tier=tier,
                        threshold=threshold,
                        pattern_config=pattern_config,
                    )
                )
    slopes = pd.concat(slopes_parts, ignore_index=True) if slopes_parts else pd.DataFrame()
    within = pd.DataFrame(within_rows)
    between = pd.DataFrame(between_rows)
    if "p_value" in within.columns:
        within["q_within_stratum_tier"] = within.groupby(
            ["stratum", "support_tier"], group_keys=False
        )["p_value"].apply(_bh)
        within["syndrome_vector_supported"] = within["q_within_stratum_tier"].le(0.05).fillna(False)
    if "northern_classic_projection_one_sided_p" in within.columns:
        within["q_classic_projection"] = within.groupby(
            ["stratum", "support_tier"], group_keys=False
        )["northern_classic_projection_one_sided_p"].apply(_bh)
        within["classic_projection_supported"] = within["q_classic_projection"].le(0.05).fillna(False)
    if "p_value" in between.columns:
        between["q_within_stratum_tier"] = between.groupby(
            ["stratum", "support_tier"], group_keys=False
        )["p_value"].apply(_bh)
        between["regional_syndrome_vector_difference_supported"] = between[
            "q_within_stratum_tier"
        ].le(0.05).fillna(False)
    return species_scores, island_scores, slopes, within, between


@app.command("run")
def run(
    trait_ledger_csv: Path = typer.Option(..., exists=True),
    status_flora_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    pattern_config_path: Path = typer.Option(..., exists=True),
    syndrome_config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    pattern_config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    syndrome_config = yaml.safe_load(syndrome_config_path.read_text(encoding="utf-8"))
    species_scores, island_scores, slopes, within, between = run_syndrome_analysis(
        pd.read_csv(trait_ledger_csv),
        pd.read_csv(status_flora_csv),
        pd.read_csv(covariates_csv),
        pattern_config,
        syndrome_config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    species_scores.to_csv(output_dir / "species_syndrome_concordance.csv.gz", index=False)
    island_scores.to_csv(output_dir / "island_syndrome_scores.csv.gz", index=False)
    slopes.to_csv(output_dir / "syndrome_distance_slopes.csv", index=False)
    within.to_csv(output_dir / "syndrome_within_region_omnibus.csv", index=False)
    between.to_csv(output_dir / "syndrome_between_region_omnibus.csv", index=False)
    manifest = {
        "contract": "chapter1_pr138_syndrome_analysis_v1",
        "syndrome_scores_are_pollinator_observations": False,
        "selfing_syndrome_separate": True,
        "n_species_scores": int(len(species_scores)),
        "n_island_scores": int(len(island_scores)),
    }
    (output_dir / "syndrome_analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
