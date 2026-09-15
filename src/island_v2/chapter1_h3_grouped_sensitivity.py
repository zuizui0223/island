"""Matched grouped-binomial sensitivity for redesigned H3 area moderation."""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_all_data_ladder import _prepare
from island_v2.chapter1_all_data_probability import _chi_square_sf_integer_df, build_broad_counts
from island_v2.chapter1_context_analysis import _fit_grouped_binomial_design

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _z_masked(frame: pd.DataFrame, mask: np.ndarray, column: str) -> np.ndarray:
    values = pd.to_numeric(frame.loc[mask, column], errors="coerce").to_numpy(float)
    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError(f"constant predictor {column}")
    out = np.zeros(len(frame), dtype=float)
    out[mask] = (values - mean) / sd
    return out


def _joint(coefficients: pd.DataFrame, covariance: np.ndarray, names: list[str], all_names: list[str]):
    indexed = coefficients.set_index("predictor")["estimate_log_odds"]
    vector = np.array([float(indexed[name]) for name in names])
    indices = [all_names.index(name) for name in names]
    cov = covariance[np.ix_(indices, indices)]
    rank = int(np.linalg.matrix_rank(cov))
    stat = float(vector @ np.linalg.pinv(cov) @ vector) if rank > 0 else float("nan")
    p = _chi_square_sf_integer_df(stat, rank) if rank > 0 else float("nan")
    return stat, rank, p


def _within(
    data: pd.DataFrame,
    probability_config: dict[str, Any],
    ladder_config: dict[str, Any],
    *,
    stratum: str,
    context_value: str,
) -> dict[str, Any]:
    context = str(ladder_config["context_column"])
    cluster = str(ladder_config["cluster_column"])
    area = str(ladder_config["area_column"])
    distance = str(ladder_config["geography_column"])
    controls = [str(x) for x in ladder_config["control_columns"]]
    threshold = int(ladder_config["minimum_islands_per_outcome"])
    outcomes = [str(x) for x in probability_config["model_outcomes"]]
    work = data.loc[data["stratum"].eq(stratum) & data[context].eq(context_value)].copy()
    support = work.groupby("outcome")["island_id"].nunique()
    retained = [outcome for outcome in outcomes if int(support.get(outcome, 0)) >= threshold]
    if len(retained) < int(ladder_config["minimum_outcomes_per_vector"]):
        return {"stratum": stratum, "context": context_value, "status": "not_testable"}
    work = work.loc[work["outcome"].isin(retained)].copy().reset_index(drop=True)
    columns: list[np.ndarray] = []
    names: list[str] = []
    targets: list[str] = []
    for outcome in retained:
        mask = work["outcome"].eq(outcome).to_numpy()
        indicator = mask.astype(float)
        columns.append(indicator)
        names.append(f"{outcome}:intercept")
        for predictor in controls:
            z = _z_masked(work, mask, predictor)
            columns.append(z)
            names.append(f"{outcome}:z_{predictor}")
        za = _z_masked(work, mask, area)
        zd = _z_masked(work, mask, distance)
        columns.extend([za, zd, za * zd])
        target = f"{outcome}:z_{distance}:z_{area}"
        names.extend([f"{outcome}:z_{area}", f"{outcome}:z_{distance}", target])
        targets.append(target)
    coefficients, fit, covariance = _fit_grouped_binomial_design(
        work["successes"].to_numpy(float),
        work["trials"].to_numpy(float),
        np.column_stack(columns),
        names,
        work[cluster].to_numpy(str),
    )
    stat, df, p = _joint(coefficients, covariance, targets, names)
    return {
        "stratum": stratum,
        "context": context_value,
        "status": fit["status"],
        "n_retained_outcomes": len(retained),
        "retained_outcomes": "|".join(retained),
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(work[cluster].nunique()),
        "joint_wald_chisq": stat,
        "joint_df": df,
        "p_value": p,
    }


def _between(
    data: pd.DataFrame,
    probability_config: dict[str, Any],
    ladder_config: dict[str, Any],
    *,
    stratum: str,
    context_a: str,
    context_b: str,
) -> dict[str, Any]:
    context = str(ladder_config["context_column"])
    cluster = str(ladder_config["cluster_column"])
    area = str(ladder_config["area_column"])
    distance = str(ladder_config["geography_column"])
    controls = [str(x) for x in ladder_config["control_columns"]]
    threshold = int(ladder_config["minimum_islands_per_outcome"])
    outcomes = [str(x) for x in probability_config["model_outcomes"]]
    work = data.loc[data["stratum"].eq(stratum) & data[context].isin([context_a, context_b])].copy()
    support = work.groupby(["outcome", context])["island_id"].nunique().unstack(fill_value=0)
    for value in (context_a, context_b):
        if value not in support.columns:
            support[value] = 0
    retained = [
        outcome for outcome in outcomes
        if outcome in support.index
        and int(support.loc[outcome, context_a]) >= threshold
        and int(support.loc[outcome, context_b]) >= threshold
    ]
    if len(retained) < int(ladder_config["minimum_outcomes_per_vector"]):
        return {"stratum": stratum, "context_a": context_a, "context_b": context_b, "status": "not_testable"}
    work = work.loc[work["outcome"].isin(retained)].copy().reset_index(drop=True)
    b = work[context].eq(context_b).to_numpy(float)
    columns: list[np.ndarray] = []
    names: list[str] = []
    targets: list[str] = []
    for outcome in retained:
        mask = work["outcome"].eq(outcome).to_numpy()
        indicator = mask.astype(float)
        columns.extend([indicator, indicator * b])
        names.extend([f"{outcome}:intercept", f"{outcome}:context[{context_b}]"])
        for predictor in controls:
            z = _z_masked(work, mask, predictor)
            columns.extend([z, z * b])
            names.extend([f"{outcome}:z_{predictor}", f"{outcome}:z_{predictor}:context[{context_b}]"])
        za = _z_masked(work, mask, area)
        zd = _z_masked(work, mask, distance)
        columns.extend([za, za * b, zd, zd * b, za * zd, za * zd * b])
        target = f"{outcome}:z_{distance}:z_{area}:context[{context_b}]"
        names.extend(
            [
                f"{outcome}:z_{area}",
                f"{outcome}:z_{area}:context[{context_b}]",
                f"{outcome}:z_{distance}",
                f"{outcome}:z_{distance}:context[{context_b}]",
                f"{outcome}:z_{distance}:z_{area}",
                target,
            ]
        )
        targets.append(target)
    coefficients, fit, covariance = _fit_grouped_binomial_design(
        work["successes"].to_numpy(float),
        work["trials"].to_numpy(float),
        np.column_stack(columns),
        names,
        work[cluster].to_numpy(str),
    )
    stat, df, p = _joint(coefficients, covariance, targets, names)
    return {
        "stratum": stratum,
        "context_a": context_a,
        "context_b": context_b,
        "status": fit["status"],
        "n_retained_outcomes": len(retained),
        "retained_outcomes": "|".join(retained),
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(work[cluster].nunique()),
        "joint_wald_chisq": stat,
        "joint_df": df,
        "p_value": p,
    }


def run_grouped_h3(counts, covariates, probability_config, ladder_config):
    data = _prepare(counts, covariates, probability_config, ladder_config)
    within = []
    between = []
    for stratum in [str(x) for x in probability_config["strata"]]:
        for context_value in [str(x) for x in ladder_config["contexts"]]:
            within.append(_within(data, probability_config, ladder_config, stratum=stratum, context_value=context_value))
        for context_a, context_b in ladder_config["primary_between_contexts"]:
            between.append(_between(data, probability_config, ladder_config, stratum=stratum, context_a=str(context_a), context_b=str(context_b)))
    return pd.DataFrame(within), pd.DataFrame(between)


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
    counts = build_broad_counts(pd.read_csv(status_flora_csv), pd.read_csv(state_audit_csv), probability_config)
    within, between = run_grouped_h3(counts, pd.read_csv(covariates_csv), probability_config, ladder_config)
    output_dir.mkdir(parents=True, exist_ok=True)
    for frame in (within, between):
        if not frame.empty:
            frame.insert(0, "evidence_scope", evidence_scope)
    within.to_csv(output_dir / "h3_grouped_within_omnibus.csv", index=False)
    between.to_csv(output_dir / "h3_grouped_between_omnibus.csv", index=False)


if __name__ == "__main__":
    app()
