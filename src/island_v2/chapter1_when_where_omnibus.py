"""Omnibus Chapter 1 tests for when and where isolation-associated filtering is detectable.

The module deliberately shifts the inferential target from individual atomic categories
to response vectors.

1. WHERE: within each biogeographic context and floristic-status stratum, test
   whether the vector of category-specific isolation slopes differs jointly from zero.
2. BETWEEN-WHERE: for each pair of contexts, test whether the vector of isolation
   slopes differs jointly between contexts.
3. WHEN: retain the same tests separately for all-native, native-nonendemic and
   endemic strata so persistence across floristic-status conditions is visible.

All models use stacked grouped-binomial rows and a cluster-robust covariance at the
predeclared spatial-block level. Category dependence within islands/blocks therefore
enters the sandwich covariance rather than being ignored by combining independent
p-values.
"""

from __future__ import annotations

import itertools
import json
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


def _scale_from_unique_islands(data: pd.DataFrame, columns: list[str]) -> dict[str, tuple[float, float]]:
    island_cov = data[["island_id", *columns]].drop_duplicates("island_id")
    scaling: dict[str, tuple[float, float]] = {}
    for column in columns:
        x = pd.to_numeric(island_cov[column], errors="coerce")
        mean = float(x.mean())
        sd = float(x.std(ddof=0))
        if not np.isfinite(sd) or sd <= 0:
            raise typer.BadParameter(f"constant or invalid omnibus predictor: {column}")
        scaling[column] = (mean, sd)
    return scaling


def _joint_wald(beta: np.ndarray, covariance: np.ndarray) -> tuple[float, int, float]:
    rank = int(np.linalg.matrix_rank(covariance))
    if rank <= 0:
        return float("nan"), 0, float("nan")
    statistic = float(beta @ np.linalg.pinv(covariance) @ beta)
    return statistic, rank, _chi_square_sf_integer_df(statistic, rank)


def _prepare_rows(
    composition: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    stratum: str,
    contexts: tuple[str, ...],
    baseline: list[str],
    isolation: str,
    context_column: str,
    cluster_column: str,
) -> pd.DataFrame:
    required = {
        "island_id", "stratum", "trait_name", "trait_state", "successes", "trials"
    }
    missing = required - set(composition.columns)
    if missing:
        raise typer.BadParameter(f"composition missing columns: {sorted(missing)}")
    joined = composition.loc[composition["stratum"].eq(stratum)].merge(
        covariates,
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    numeric = ["successes", "trials", isolation, *baseline]
    for column in numeric:
        joined[column] = pd.to_numeric(joined[column], errors="coerce")
    joined = joined.dropna(subset=numeric + [context_column, cluster_column]).copy()
    joined = joined.loc[
        joined[context_column].astype(str).isin(contexts)
        & (joined["trials"] > 0)
        & (joined["successes"] >= 0)
        & (joined["successes"] <= joined["trials"])
    ].copy()
    joined["response_id"] = (
        joined["trait_name"].astype(str) + "::" + joined["trait_state"].astype(str)
    )
    joined[context_column] = joined[context_column].astype(str)
    joined[cluster_column] = joined[cluster_column].astype(str)
    return joined


def _within_context_test(
    rows: pd.DataFrame,
    *,
    context: str,
    baseline: list[str],
    isolation: str,
    context_column: str,
    cluster_column: str,
    min_islands_per_response: int,
) -> dict[str, Any] | None:
    work = rows.loc[rows[context_column].eq(context)].copy()
    counts = work.groupby("response_id")["island_id"].nunique()
    responses = sorted(counts.loc[counts >= min_islands_per_response].index)
    work = work.loc[work["response_id"].isin(responses)].copy()
    if len(responses) < 2:
        return None
    scaling = _scale_from_unique_islands(work, [*baseline, isolation])
    response_levels = responses
    names: list[str] = []
    columns: list[np.ndarray] = []
    for response in response_levels:
        mask = work["response_id"].eq(response).to_numpy(float)
        names.append(f"response[{response}]")
        columns.append(mask)
        for predictor in baseline:
            mean, sd = scaling[predictor]
            z = (work[predictor].to_numpy(float) - mean) / sd
            names.append(f"response[{response}]:z_{predictor}")
            columns.append(mask * z)
        mean, sd = scaling[isolation]
        z_iso = (work[isolation].to_numpy(float) - mean) / sd
        names.append(f"response[{response}]:z_{isolation}")
        columns.append(mask * z_iso)
    X = np.column_stack(columns)
    coef, fit, covariance = _fit_grouped_binomial_design(
        work["successes"].to_numpy(float),
        work["trials"].to_numpy(float),
        X,
        names,
        work[cluster_column].to_numpy(str),
    )
    beta = coef.set_index("predictor")["estimate_log_odds"]
    iso_names = [f"response[{response}]:z_{isolation}" for response in response_levels]
    indices = [names.index(name) for name in iso_names]
    vector = np.array([float(beta[name]) for name in iso_names])
    cov = covariance[np.ix_(indices, indices)]
    statistic, df, p_value = _joint_wald(vector, cov)
    return {
        "context": context,
        "n_responses": len(response_levels),
        "responses": "|".join(response_levels),
        "min_response_islands": int(min(counts.loc[response_levels])),
        "max_response_islands": int(max(counts.loc[response_levels])),
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "joint_wald_chisq": statistic,
        "joint_df": df,
        "p_value": p_value,
    }


def _between_context_test(
    rows: pd.DataFrame,
    *,
    context_a: str,
    context_b: str,
    baseline: list[str],
    isolation: str,
    context_column: str,
    cluster_column: str,
    min_islands_per_response: int,
) -> dict[str, Any] | None:
    work = rows.loc[rows[context_column].isin([context_a, context_b])].copy()
    counts = (
        work.groupby(["response_id", context_column])["island_id"].nunique().unstack(fill_value=0)
    )
    if context_a not in counts.columns or context_b not in counts.columns:
        return None
    eligible = counts.index[
        (counts[context_a] >= min_islands_per_response)
        & (counts[context_b] >= min_islands_per_response)
    ]
    responses = sorted(eligible.astype(str))
    work = work.loc[work["response_id"].isin(responses)].copy()
    if len(responses) < 2:
        return None
    scaling = _scale_from_unique_islands(work, [*baseline, isolation])
    context_b_indicator = work[context_column].eq(context_b).to_numpy(float)
    names: list[str] = []
    columns: list[np.ndarray] = []
    interaction_names: list[str] = []
    for response in responses:
        mask = work["response_id"].eq(response).to_numpy(float)
        names.append(f"response[{response}]")
        columns.append(mask)
        for predictor in baseline:
            mean, sd = scaling[predictor]
            z = (work[predictor].to_numpy(float) - mean) / sd
            names.append(f"response[{response}]:z_{predictor}")
            columns.append(mask * z)
        names.append(f"response[{response}]:context[{context_b}]")
        columns.append(mask * context_b_indicator)
        mean, sd = scaling[isolation]
        z_iso = (work[isolation].to_numpy(float) - mean) / sd
        names.append(f"response[{response}]:z_{isolation}")
        columns.append(mask * z_iso)
        interaction = f"response[{response}]:z_{isolation}:context[{context_b}]"
        interaction_names.append(interaction)
        names.append(interaction)
        columns.append(mask * z_iso * context_b_indicator)
    X = np.column_stack(columns)
    coef, fit, covariance = _fit_grouped_binomial_design(
        work["successes"].to_numpy(float),
        work["trials"].to_numpy(float),
        X,
        names,
        work[cluster_column].to_numpy(str),
    )
    beta = coef.set_index("predictor")["estimate_log_odds"]
    indices = [names.index(name) for name in interaction_names]
    vector = np.array([float(beta[name]) for name in interaction_names])
    cov = covariance[np.ix_(indices, indices)]
    statistic, df, p_value = _joint_wald(vector, cov)
    return {
        "context_a": context_a,
        "context_b": context_b,
        "n_responses": len(responses),
        "responses": "|".join(responses),
        "min_islands_context_a": int(counts.loc[responses, context_a].min()),
        "min_islands_context_b": int(counts.loc[responses, context_b].min()),
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "joint_wald_chisq": statistic,
        "joint_df": df,
        "p_value": p_value,
    }


def run_when_where_omnibus(
    composition: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    baseline = [str(x) for x in config["baseline_covariates"]]
    isolation = str(config["isolation_column"])
    context_column = str(config["context_column"])
    cluster_column = str(config["cluster_column"])
    contexts = tuple(str(x) for x in config["contexts"])
    strata = tuple(str(x) for x in config["strata"])
    tiers = config.get("support_tiers", {"confirmatory": 50, "pilot": 30})

    within_rows: list[dict[str, Any]] = []
    between_rows: list[dict[str, Any]] = []
    for stratum in strata:
        base_rows = _prepare_rows(
            composition,
            covariates,
            stratum=stratum,
            contexts=contexts,
            baseline=baseline,
            isolation=isolation,
            context_column=context_column,
            cluster_column=cluster_column,
        )
        for tier, minimum in tiers.items():
            minimum = int(minimum)
            for context in contexts:
                result = _within_context_test(
                    base_rows,
                    context=context,
                    baseline=baseline,
                    isolation=isolation,
                    context_column=context_column,
                    cluster_column=cluster_column,
                    min_islands_per_response=minimum,
                )
                if result is not None:
                    within_rows.append({
                        "stratum": stratum,
                        "support_tier": str(tier),
                        "threshold": minimum,
                        **result,
                    })
            for context_a, context_b in itertools.combinations(contexts, 2):
                result = _between_context_test(
                    base_rows,
                    context_a=context_a,
                    context_b=context_b,
                    baseline=baseline,
                    isolation=isolation,
                    context_column=context_column,
                    cluster_column=cluster_column,
                    min_islands_per_response=minimum,
                )
                if result is not None:
                    between_rows.append({
                        "stratum": stratum,
                        "support_tier": str(tier),
                        "threshold": minimum,
                        **result,
                    })

    within = pd.DataFrame(within_rows)
    between = pd.DataFrame(between_rows)
    if not within.empty:
        within["q_within_stratum_tier"] = within.groupby(
            ["stratum", "support_tier"], group_keys=False
        )["p_value"].apply(_bh)
        within["where_supported"] = within["q_within_stratum_tier"].le(0.05)
    if not between.empty:
        between["q_within_stratum_tier"] = between.groupby(
            ["stratum", "support_tier"], group_keys=False
        )["p_value"].apply(_bh)
        between["between_where_supported"] = between["q_within_stratum_tier"].le(0.05)

    classification_rows: list[dict[str, Any]] = []
    if not within.empty:
        confirmatory = within.loc[within["support_tier"].eq("confirmatory")]
        for context in contexts:
            context_rows = confirmatory.loc[confirmatory["context"].eq(context)]
            all_native = context_rows.loc[context_rows["stratum"].eq("all_native")]
            nonendemic = context_rows.loc[context_rows["stratum"].eq("native_nonendemic")]
            endemic = context_rows.loc[context_rows["stratum"].eq("endemic")]
            def supported(frame: pd.DataFrame) -> bool:
                return bool(not frame.empty and frame.iloc[0].get("where_supported", False))
            if supported(all_native) and supported(nonendemic):
                when_class = "persists_in_native_nonendemic"
            elif supported(endemic) and not supported(nonendemic):
                when_class = "endemic_concentrated"
            elif supported(all_native):
                when_class = "all_native_only_or_unresolved"
            else:
                when_class = "no_confirmatory_vector_signal"
            classification_rows.append(
                {
                    "context": context,
                    "confirmatory_all_native": supported(all_native),
                    "confirmatory_native_nonendemic": supported(nonendemic),
                    "confirmatory_endemic": supported(endemic),
                    "when_class": when_class,
                }
            )
    return within, between, pd.DataFrame(classification_rows)


@app.command("run")
def run(
    composition_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/chapter1_when_where_omnibus.yml"), exists=True),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    within, between, classification = run_when_where_omnibus(
        pd.read_csv(composition_csv), pd.read_csv(covariates_csv), config
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    within.to_csv(output_dir / "when_where_within_context_omnibus.csv", index=False)
    between.to_csv(output_dir / "when_where_between_context_omnibus.csv", index=False)
    classification.to_csv(output_dir / "when_where_status_persistence.csv", index=False)
    manifest = {
        "contract": "chapter1_when_where_omnibus_v1",
        "primary_where_test": "within-context joint Wald of all supported atomic isolation slopes",
        "between_where_test": "pairwise joint Wald of atomic isolation-slope differences",
        "when_operationalization": "persistence across floristic-status strata",
        "cluster_robust": True,
        "pollinator_predictors": False,
        "config": config,
        "n_within_tests": int(len(within)),
        "n_between_tests": int(len(between)),
    }
    (output_dir / "chapter1_when_where_omnibus_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
