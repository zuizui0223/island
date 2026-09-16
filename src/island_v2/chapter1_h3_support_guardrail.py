"""Trait-support guardrail for redesigned H3 area moderation.

The primary beta-binomial model legitimately uses outcome-specific denominators, but
trait support itself can vary with island area. This module therefore repeats the
predeclared North--Tropical distance x area x context test under three alternative
information contracts that do not use trait outcomes to choose islands or weights:

1. capped information: cap each island x outcome denominator at 50 while preserving
   the empirical trait fraction;
2. equal island: one unit of information per island x outcome cell;
3. common support-area range + equal island: split by median denominator within each
   context x outcome, restrict both support halves to their overlapping area 5--95%
   range, then use equal island information.
"""
from __future__ import annotations

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
    x = pd.to_numeric(frame.loc[mask, column], errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError(f"constant predictor: {column}")
    out = np.zeros(len(frame), dtype=float)
    out[mask] = (x - mean) / sd
    return out


def _common_support_area_filter(
    frame: pd.DataFrame,
    *,
    context_column: str,
    area_column: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    parts = []
    audit_rows = []
    for (context_value, outcome), group in frame.groupby([context_column, "outcome"], sort=True):
        work = group.copy()
        support = pd.to_numeric(work["trials"], errors="coerce")
        area = pd.to_numeric(work[area_column], errors="coerce")
        median_support = float(support.median())
        low = work.loc[support.le(median_support)].copy()
        high = work.loc[support.gt(median_support)].copy()
        if len(low) < 20 or len(high) < 20:
            audit_rows.append(
                {
                    "context": context_value,
                    "outcome": outcome,
                    "status": "insufficient_support_halves",
                    "n_before": len(work),
                    "n_after": 0,
                }
            )
            continue
        low_area = pd.to_numeric(low[area_column], errors="coerce")
        high_area = pd.to_numeric(high[area_column], errors="coerce")
        overlap_low = max(float(low_area.quantile(0.05)), float(high_area.quantile(0.05)))
        overlap_high = min(float(low_area.quantile(0.95)), float(high_area.quantile(0.95)))
        if not np.isfinite(overlap_low) or not np.isfinite(overlap_high) or overlap_low >= overlap_high:
            audit_rows.append(
                {
                    "context": context_value,
                    "outcome": outcome,
                    "status": "no_common_area_range",
                    "n_before": len(work),
                    "n_after": 0,
                }
            )
            continue
        kept = work.loc[area.between(overlap_low, overlap_high)].copy()
        parts.append(kept)
        audit_rows.append(
            {
                "context": context_value,
                "outcome": outcome,
                "status": "retained",
                "median_support": median_support,
                "area_overlap_low": overlap_low,
                "area_overlap_high": overlap_high,
                "n_before": len(work),
                "n_after": len(kept),
            }
        )
    return (
        pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=frame.columns),
        pd.DataFrame(audit_rows),
    )


def _support_area_audit(frame: pd.DataFrame, *, context_column: str, area_column: str) -> pd.DataFrame:
    rows = []
    for (stratum, context_value, outcome), group in frame.groupby(
        ["stratum", context_column, "outcome"], sort=True
    ):
        support = pd.to_numeric(group["trials"], errors="coerce")
        area = pd.to_numeric(group[area_column], errors="coerce")
        rho = float(support.rank().corr(area.rank())) if len(group) > 2 else float("nan")
        rows.append(
            {
                "stratum": stratum,
                "context": context_value,
                "outcome": outcome,
                "n_islands": int(group["island_id"].nunique()),
                "median_trials": float(support.median()),
                "spearman_trials_area": rho,
            }
        )
    return pd.DataFrame(rows)


def _fit_between(
    data: pd.DataFrame,
    probability_config: dict[str, Any],
    ladder_config: dict[str, Any],
    *,
    stratum: str,
    context_a: str,
    context_b: str,
    mode: str,
) -> tuple[dict[str, Any], pd.DataFrame]:
    context = str(ladder_config["context_column"])
    cluster = str(ladder_config["cluster_column"])
    area = str(ladder_config["area_column"])
    distance = str(ladder_config["geography_column"])
    controls = [str(x) for x in ladder_config["control_columns"]]
    threshold = int(ladder_config["minimum_islands_per_outcome"])
    outcomes = [str(x) for x in probability_config["model_outcomes"]]
    work = data.loc[data["stratum"].eq(stratum) & data[context].isin([context_a, context_b])].copy()
    filter_audit = pd.DataFrame()
    if mode == "common_support_equal":
        work, filter_audit = _common_support_area_filter(
            work, context_column=context, area_column=area
        )
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
        return (
            {
                "stratum": stratum,
                "context_a": context_a,
                "context_b": context_b,
                "mode": mode,
                "status": "not_testable",
                "n_retained_outcomes": len(retained),
                "retained_outcomes": "|".join(retained),
            },
            filter_audit,
        )
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

    fraction = work["successes"].to_numpy(float) / work["trials"].to_numpy(float)
    original_trials = work["trials"].to_numpy(float)
    if mode == "capped_50":
        effective_trials = np.minimum(original_trials, 50.0)
    elif mode in {"equal_island", "common_support_equal"}:
        effective_trials = np.ones(len(work), dtype=float)
    else:
        raise ValueError(f"unknown support mode: {mode}")
    effective_successes = fraction * effective_trials
    coefficients, fit, covariance = _fit_grouped_binomial_design(
        effective_successes,
        effective_trials,
        np.column_stack(columns),
        names,
        work[cluster].to_numpy(str),
    )
    indexed = coefficients.set_index("predictor")["estimate_log_odds"]
    vector = np.array([float(indexed[name]) for name in targets])
    indices = [names.index(name) for name in targets]
    cov = covariance[np.ix_(indices, indices)]
    rank = int(np.linalg.matrix_rank(cov))
    statistic = float(vector @ np.linalg.pinv(cov) @ vector) if rank > 0 else float("nan")
    p_value = _chi_square_sf_integer_df(statistic, rank) if rank > 0 else float("nan")
    return (
        {
            "stratum": stratum,
            "context_a": context_a,
            "context_b": context_b,
            "mode": mode,
            "status": fit["status"],
            "n_retained_outcomes": len(retained),
            "retained_outcomes": "|".join(retained),
            "n_unique_islands": int(work["island_id"].nunique()),
            "n_clusters": int(work[cluster].nunique()),
            "joint_wald_chisq": statistic,
            "joint_df": rank,
            "p_value": p_value,
        },
        filter_audit,
    )


def run_guardrail(counts, covariates, probability_config, ladder_config):
    data = _prepare(counts, covariates, probability_config, ladder_config)
    support_audit = _support_area_audit(
        data,
        context_column=str(ladder_config["context_column"]),
        area_column=str(ladder_config["area_column"]),
    )
    rows = []
    filter_parts = []
    for stratum in [str(x) for x in probability_config["strata"]]:
        for context_a, context_b in ladder_config["primary_between_contexts"]:
            for mode in ("capped_50", "equal_island", "common_support_equal"):
                result, filter_audit = _fit_between(
                    data,
                    probability_config,
                    ladder_config,
                    stratum=stratum,
                    context_a=str(context_a),
                    context_b=str(context_b),
                    mode=mode,
                )
                rows.append(result)
                if not filter_audit.empty:
                    filter_audit = filter_audit.copy()
                    filter_audit.insert(0, "stratum", stratum)
                    filter_parts.append(filter_audit)
    return (
        pd.DataFrame(rows),
        support_audit,
        pd.concat(filter_parts, ignore_index=True) if filter_parts else pd.DataFrame(),
    )


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
    results, support, filters = run_guardrail(
        counts, pd.read_csv(covariates_csv), probability_config, ladder_config
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    for frame in (results, support, filters):
        if not frame.empty:
            frame.insert(0, "evidence_scope", evidence_scope)
    results.to_csv(output_dir / "h3_support_guardrail_results.csv", index=False)
    support.to_csv(output_dir / "h3_support_area_correlations.csv", index=False)
    filters.to_csv(output_dir / "h3_common_support_filter_audit.csv", index=False)


if __name__ == "__main__":
    app()
