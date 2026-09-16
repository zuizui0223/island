"""Conditional North--Tropical H3 interpretation at fixed area values.

After a supported distance x area x context interaction, this module evaluates the
conditional difference in distance-response vectors at z(log area) = -1, 0, +1.
These are fixed symmetric one-SD interpretation points, not searched thresholds.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_all_data_ladder import _prepare
from island_v2.chapter1_all_data_probability import (
    _assemble_cluster_covariance,
    _bh,
    _chi_square_sf_integer_df,
    _fit_single_beta_binomial,
    _normal_two_sided_p,
    _standardize,
    build_broad_counts,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _fit_one_stratum(
    data: pd.DataFrame,
    probability_config: dict[str, Any],
    ladder_config: dict[str, Any],
    *,
    stratum: str,
    context_a: str,
    context_b: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    geography = str(ladder_config["geography_column"])
    area = str(ladder_config["area_column"])
    context = str(ladder_config["context_column"])
    cluster = str(ladder_config["cluster_column"])
    controls = [str(x) for x in ladder_config["control_columns"]]
    threshold = int(ladder_config["minimum_islands_per_outcome"])
    outcomes = [str(x) for x in probability_config["model_outcomes"]]
    area_values = [float(x) for x in ladder_config["hypotheses"]["H3_island_capacity_moderation"]["conditional_area_z_values"]]

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
        return pd.DataFrame(), pd.DataFrame(
            [
                {
                    "stratum": stratum,
                    "context_a": context_a,
                    "context_b": context_b,
                    "status": "not_testable",
                    "n_retained_outcomes": len(retained),
                    "retained_outcomes": "|".join(retained),
                }
            ]
        )

    fits: list[dict[str, Any]] = []
    clusters: list[np.ndarray] = []
    distance_context_indices: list[int] = []
    triple_indices: list[int] = []
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
        distance_context_name = f"{outcome}:z_{geography}:context[{context_b}]"
        names.extend(
            [
                f"{outcome}:z_{area}",
                f"{outcome}:z_{area}:context[{context_b}]",
                f"{outcome}:z_{geography}",
                distance_context_name,
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
        distance_context_indices.append(offset + names.index(distance_context_name))
        triple_indices.append(offset + names.index(triple_name))
        offset += len(fit["names"])

    covariance, _, theta = _assemble_cluster_covariance(fits, clusters)
    d = theta[distance_context_indices]
    t = theta[triple_indices]
    c_dd = covariance[np.ix_(distance_context_indices, distance_context_indices)]
    c_tt = covariance[np.ix_(triple_indices, triple_indices)]
    c_dt = covariance[np.ix_(distance_context_indices, triple_indices)]

    slope_rows: list[dict[str, Any]] = []
    omnibus_rows: list[dict[str, Any]] = []
    for area_z in area_values:
        vector = d + area_z * t
        cov = c_dd + (area_z**2) * c_tt + area_z * (c_dt + c_dt.T)
        se = np.sqrt(np.clip(np.diag(cov), 0.0, None))
        for outcome, estimate, stderr in zip(retained, vector, se, strict=True):
            z = estimate / stderr if stderr > 0 else float("nan")
            slope_rows.append(
                {
                    "stratum": stratum,
                    "context_a": context_a,
                    "context_b": context_b,
                    "area_z": area_z,
                    "outcome": outcome,
                    "conditional_distance_difference_b_minus_a": float(estimate),
                    "cluster_robust_se": float(stderr),
                    "p_value": _normal_two_sided_p(float(z)),
                }
            )
        rank = int(np.linalg.matrix_rank(cov))
        statistic = float(vector @ np.linalg.pinv(cov) @ vector) if rank > 0 else float("nan")
        p_value = _chi_square_sf_integer_df(statistic, rank) if rank > 0 else float("nan")
        omnibus_rows.append(
            {
                "stratum": stratum,
                "context_a": context_a,
                "context_b": context_b,
                "area_z": area_z,
                "status": "fit",
                "n_retained_outcomes": len(retained),
                "retained_outcomes": "|".join(retained),
                "n_unique_islands": int(work.loc[work["outcome"].isin(retained), "island_id"].nunique()),
                "n_clusters": int(work.loc[work["outcome"].isin(retained), cluster].nunique()),
                "joint_wald_chisq": statistic,
                "joint_df": rank,
                "p_value": p_value,
                "all_optimizers_converged": all(bool(fit["success"]) for fit in fits),
            }
        )
    omnibus = pd.DataFrame(omnibus_rows)
    if not omnibus.empty:
        omnibus["q_value_across_fixed_area_points"] = _bh(omnibus["p_value"])
        omnibus["conditional_difference_supported"] = omnibus[
            "q_value_across_fixed_area_points"
        ].le(float(ladder_config["alpha"])).fillna(False)
    return pd.DataFrame(slope_rows), omnibus


def run_conditional_area(
    counts: pd.DataFrame,
    covariates: pd.DataFrame,
    probability_config: dict[str, Any],
    ladder_config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    data = _prepare(counts, covariates, probability_config, ladder_config)
    slope_parts: list[pd.DataFrame] = []
    omnibus_parts: list[pd.DataFrame] = []
    for stratum in [str(x) for x in probability_config["strata"]]:
        for context_a, context_b in ladder_config["primary_between_contexts"]:
            slopes, omnibus = _fit_one_stratum(
                data,
                probability_config,
                ladder_config,
                stratum=stratum,
                context_a=str(context_a),
                context_b=str(context_b),
            )
            if not slopes.empty:
                slope_parts.append(slopes)
            omnibus_parts.append(omnibus)
    return (
        pd.concat(slope_parts, ignore_index=True) if slope_parts else pd.DataFrame(),
        pd.concat(omnibus_parts, ignore_index=True) if omnibus_parts else pd.DataFrame(),
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
    counts = build_broad_counts(
        pd.read_csv(status_flora_csv), pd.read_csv(state_audit_csv), probability_config
    )
    slopes, omnibus = run_conditional_area(
        counts, pd.read_csv(covariates_csv), probability_config, ladder_config
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    for frame in (slopes, omnibus):
        if not frame.empty:
            frame.insert(0, "evidence_scope", evidence_scope)
    slopes.to_csv(output_dir / "h3_conditional_area_slopes.csv", index=False)
    omnibus.to_csv(output_dir / "h3_conditional_area_omnibus.csv", index=False)
    manifest = {
        "contract": ladder_config["contract"],
        "evidence_scope": evidence_scope,
        "area_z_values": ladder_config["hypotheses"]["H3_island_capacity_moderation"]["conditional_area_z_values"],
        "interpretation": "fixed symmetric one-SD points; no threshold search",
    }
    (output_dir / "h3_conditional_area_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    app()
