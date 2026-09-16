"""Pooled five-channel pollination-disruption bridge for Chapter 1 H5.

This module deliberately does not call cross-sectional channel states temporal
"decline".  A disrupted state means that a pollination channel was source-available
but adequately non-detected on the island under the frozen effort policy.  The bridge
asks whether that conservative channel-disruption signal aligns with two plant-side
responses: reproductive assurance and pollinator-facing floral architecture.
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

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _normal_two_sided_p(z: float) -> float:
    return math.erfc(abs(float(z)) / math.sqrt(2.0)) if math.isfinite(z) else float("nan")


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != "chapter1_h5_multichannel_two_route_v1":
        raise typer.BadParameter("unexpected multichannel two-route contract")
    return config


def _to_numeric(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = frame.copy()
    for column in columns:
        out[column] = pd.to_numeric(out[column], errors="coerce")
    return out


def build_composite_channel_exposure(
    observations: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {
        "island_id",
        "channel_id",
        "observation_state",
        "background_record_count",
        "background_spatial_units",
        "background_temporal_units",
        "distinct_dataset_count",
        "latest_background_year",
    }
    if missing := required - set(observations.columns):
        raise typer.BadParameter(f"channel observation table missing columns: {sorted(missing)}")

    channels = [str(x) for x in config["channels"]]
    work = observations.loc[observations["channel_id"].astype(str).isin(channels)].copy()
    work["island_id"] = work["island_id"].fillna("").astype(str)
    work["channel_id"] = work["channel_id"].fillna("").astype(str)
    if work.duplicated(["island_id", "channel_id"]).any():
        raise typer.BadParameter("channel observations must be unique by island_id x channel_id")

    numeric = [
        "background_record_count",
        "background_spatial_units",
        "background_temporal_units",
        "distinct_dataset_count",
        "latest_background_year",
    ]
    work = _to_numeric(work, numeric)
    gate = config["symmetric_effort_gate"]
    min_year = int(gate["reference_year"]) - int(gate["max_years_since_latest_background"])
    work["symmetric_effort_qualified"] = (
        work["background_record_count"].ge(float(gate["min_background_records"]))
        & work["background_spatial_units"].ge(float(gate["min_background_spatial_units"]))
        & work["background_temporal_units"].ge(float(gate["min_background_temporal_units"]))
        & work["distinct_dataset_count"].ge(float(gate["min_distinct_datasets"]))
        & work["latest_background_year"].ge(float(min_year))
    )
    work["retained_symmetric"] = (
        work["observation_state"].astype(str).eq("detected")
        & work["symmetric_effort_qualified"]
    )
    work["disrupted_strict"] = work["observation_state"].astype(str).eq("adequate_non_detection")
    work["evaluable_symmetric"] = work["retained_symmetric"] | work["disrupted_strict"]

    rows: list[dict[str, Any]] = []
    for island_id, group in work.groupby("island_id", sort=True):
        row: dict[str, Any] = {
            "island_id": island_id,
            "n_source_available_channels": int(group["channel_id"].nunique()),
            "n_evaluable_channels": int(group["evaluable_symmetric"].sum()),
            "n_retained_channels": int(group["retained_symmetric"].sum()),
            "n_disrupted_channels": int(group["disrupted_strict"].sum()),
        }
        row["any_channel_disrupted"] = float(row["n_disrupted_channels"] > 0)
        for channel in channels:
            part = group.loc[group["channel_id"].eq(channel)]
            row[f"eval_{channel}"] = float(
                bool(part["evaluable_symmetric"].iloc[0]) if len(part) else False
            )
        rows.append(row)
    return pd.DataFrame(rows), work


def build_plant_scores(scores: pd.DataFrame, stratum: str) -> pd.DataFrame:
    required = {"island_id", "stratum", "syndrome", "syndrome_score"}
    if missing := required - set(scores.columns):
        raise typer.BadParameter(f"syndrome score table missing columns: {sorted(missing)}")
    part = scores.loc[scores["stratum"].astype(str).eq(stratum)].copy()
    wide = part.pivot_table(index="island_id", columns="syndrome", values="syndrome_score", aggfunc="first")
    needed = {
        "selfing_core",
        "selfing_syndrome",
        "generalized_accessible",
        "large_bee_like",
        "butterfly_like",
        "bird_like",
    }
    if missing := needed - set(wide.columns):
        raise typer.BadParameter(f"syndrome score table lacks required axes: {sorted(missing)}")
    wide["attraction_shift"] = (-wide["large_bee_like"] + wide["generalized_accessible"]) / 2.0
    wide["shared_named_architecture"] = wide[
        ["large_bee_like", "butterfly_like", "bird_like"]
    ].mean(axis=1, skipna=False)
    return wide.reset_index()


def _z(values: pd.Series) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def _clustered_ols(
    data: pd.DataFrame,
    *,
    response: str,
    exposure: str,
    controls: list[str],
    eval_indicators: list[str],
    cluster_column: str,
    condition_on_selfing: bool,
) -> dict[str, Any]:
    columns = [response, exposure, cluster_column, *controls, "n_evaluable_channels", *eval_indicators]
    if condition_on_selfing:
        columns.append("selfing_core")
    work = data[columns].copy()
    numeric = [response, exposure, *controls, "n_evaluable_channels", *eval_indicators]
    if condition_on_selfing:
        numeric.append("selfing_core")
    work = _to_numeric(work, numeric)
    work = work.dropna(subset=numeric)
    work = work.loc[work[cluster_column].fillna("").astype(str).ne("")].copy()
    if len(work) < 20 or work[exposure].nunique() < 2:
        return {
            "status": "not_testable",
            "response": response,
            "condition_on_selfing_core": condition_on_selfing,
            "n_islands": int(len(work)),
            "n_disrupted": int(work[exposure].sum()) if len(work) else 0,
        }

    design: list[np.ndarray] = [np.ones(len(work)), work[exposure].to_numpy(float)]
    names = ["intercept", exposure]
    for control in [*controls, "n_evaluable_channels"]:
        try:
            design.append(_z(work[control]))
            names.append(f"z_{control}")
        except ValueError:
            continue

    # n_evaluable_channels equals the sum of all channel indicators.  Keep the count
    # and all but the first channel indicator to control composition without exact
    # collinearity.
    for indicator in eval_indicators[1:]:
        if work[indicator].nunique() > 1:
            design.append(work[indicator].to_numpy(float))
            names.append(indicator)
    if condition_on_selfing:
        try:
            design.append(_z(work["selfing_core"]))
            names.append("z_selfing_core")
        except ValueError:
            pass

    X = np.column_stack(design)
    y = work[response].to_numpy(float)
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

    index = names.index(exposure)
    estimate = float(beta[index])
    stderr = float(math.sqrt(max(float(covariance[index, index]), 0.0)))
    z_value = estimate / stderr if stderr > 0 else float("nan")
    return {
        "status": "fit",
        "response": response,
        "condition_on_selfing_core": condition_on_selfing,
        "n_islands": int(n),
        "n_clusters": int(g),
        "n_disrupted": int(work[exposure].sum()),
        "n_no_documented_disruption": int((1.0 - work[exposure]).sum()),
        "disruption_estimate": estimate,
        "cluster_robust_se": stderr,
        "z_value": float(z_value),
        "p_value": _normal_two_sided_p(z_value),
    }


def _bh_qvalues(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce")
    valid = p.dropna().sort_values()
    if valid.empty:
        return pd.Series(np.nan, index=values.index)
    m = len(valid)
    raw = valid.to_numpy(float) * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(raw[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out = pd.Series(np.nan, index=values.index, dtype=float)
    out.loc[valid.index] = adjusted
    return out


def run_bridge(
    syndrome_scores: pd.DataFrame,
    channel_observations: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
    *,
    evidence_scope: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    composite, channel_long = build_composite_channel_exposure(channel_observations, config)
    model = config["models"]
    context_column = str(model["context_column"])
    cluster_column = str(model["cluster_column"])
    controls = [str(x) for x in model["controls"]]
    required_cov = {"island_id", context_column, cluster_column, *controls}
    if missing := required_cov - set(covariates.columns):
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")

    eval_indicators = [f"eval_{x}" for x in config["channels"]]
    strata = [str(model["primary_flora_scope"]), *[str(x) for x in model["native_sensitivity_scopes"]]]
    min_channels_values = [
        int(config["composite_exposure"]["primary_min_evaluable_channels"]),
        *[int(x) for x in config["composite_exposure"]["sensitivity_min_evaluable_channels"]],
    ]
    min_channels_values = list(dict.fromkeys(min_channels_values))

    rows: list[dict[str, Any]] = []
    support_rows: list[dict[str, Any]] = []
    for stratum in strata:
        plant = build_plant_scores(syndrome_scores, stratum)
        data = (
            composite.merge(covariates[list(required_cov)].drop_duplicates("island_id"), on="island_id", how="left")
            .merge(plant, on="island_id", how="inner")
        )
        for minimum in min_channels_values:
            eligible = data.loc[data["n_evaluable_channels"].ge(minimum)].copy()
            for context in [str(x) for x in model["descriptive_contexts"]]:
                part = eligible.loc[eligible[context_column].astype(str).eq(context)].copy()
                support_rows.append(
                    {
                        "evidence_scope": evidence_scope,
                        "stratum": stratum,
                        "min_evaluable_channels": minimum,
                        "context": context,
                        "n_islands": int(len(part)),
                        "n_disrupted": int(part["any_channel_disrupted"].sum()) if len(part) else 0,
                        "n_no_documented_disruption": int((1.0 - part["any_channel_disrupted"]).sum()) if len(part) else 0,
                    }
                )
                specifications = [
                    ("selfing_core", False, "route_A_reproductive_assurance"),
                    ("selfing_syndrome", False, "route_A_full_selfing_syndrome"),
                    ("generalized_accessible", False, "route_B_accessibility_unconditional"),
                    ("generalized_accessible", True, "route_B_accessibility_conditional_on_selfing"),
                    ("shared_named_architecture", False, "route_B_shared_architecture_unconditional"),
                    ("shared_named_architecture", True, "route_B_shared_architecture_conditional_on_selfing"),
                    ("attraction_shift", True, "legacy_attraction_shift_conditional_on_selfing"),
                ]
                for response, conditional, role in specifications:
                    result = _clustered_ols(
                        part,
                        response=response,
                        exposure="any_channel_disrupted",
                        controls=controls,
                        eval_indicators=eval_indicators,
                        cluster_column=cluster_column,
                        condition_on_selfing=conditional,
                    )
                    result.update(
                        {
                            "evidence_scope": evidence_scope,
                            "stratum": stratum,
                            "min_evaluable_channels": minimum,
                            "context": context,
                            "analysis_role": role,
                        }
                    )
                    rows.append(result)

    results = pd.DataFrame(rows)
    support = pd.DataFrame(support_rows)
    primary = config["primary_decision"]
    primary_contexts = [str(x) for x in model["primary_contexts"]]
    primary_rows = results.loc[
        results["evidence_scope"].eq(str(primary["evidence_scope"]))
        & results["stratum"].eq(str(model["primary_flora_scope"]))
        & results["min_evaluable_channels"].eq(int(primary["min_evaluable_channels"]))
        & results["context"].isin(primary_contexts)
        & results["analysis_role"].isin(
            ["route_A_reproductive_assurance", "route_B_accessibility_conditional_on_selfing"]
        )
    ].copy()
    primary_rows["q_value"] = _bh_qvalues(primary_rows["p_value"])
    if not primary_rows.empty:
        key = list(zip(primary_rows.index, primary_rows["q_value"], strict=False))
        qmap = {index: q for index, q in key}
        results["q_value_primary_family"] = [qmap.get(index, np.nan) for index in results.index]
    else:
        results["q_value_primary_family"] = np.nan

    support_primary = support.loc[
        support["evidence_scope"].eq(str(primary["evidence_scope"]))
        & support["stratum"].eq(str(model["primary_flora_scope"]))
        & support["min_evaluable_channels"].eq(int(primary["min_evaluable_channels"]))
        & support["context"].isin(primary_contexts)
    ].copy()
    overlap = bool(
        len(support_primary) == len(primary_contexts)
        and support_primary["n_disrupted"].ge(int(primary["require_each_primary_context_disrupted_n"])).all()
        and support_primary["n_no_documented_disruption"].ge(int(primary["require_each_primary_context_retained_n"])).all()
    )
    direction_a = primary_rows.loc[primary_rows["analysis_role"].eq("route_A_reproductive_assurance"), "disruption_estimate"]
    direction_b = primary_rows.loc[
        primary_rows["analysis_role"].eq("route_B_accessibility_conditional_on_selfing"), "disruption_estimate"
    ]
    directional = bool(
        len(direction_a) == len(primary_contexts)
        and len(direction_b) == len(primary_contexts)
        and direction_a.gt(0).all()
        and direction_b.gt(0).all()
    )
    fdr_supported = bool(
        len(primary_rows) == 2 * len(primary_contexts)
        and primary_rows["q_value"].le(float(primary["fdr_alpha"])).all()
    )
    decision = {
        "contract": config["contract"],
        "evidence_scope": evidence_scope,
        "primary_min_evaluable_channels": int(primary["min_evaluable_channels"]),
        "primary_contexts": primary_contexts,
        "symmetric_state_overlap_passed": overlap,
        "two_route_directional_requirement_passed": directional,
        "primary_four_tests_fdr_supported": fdr_supported,
        "simple_global_pollinator_disruption_two_route_promoted": bool(overlap and directional and fdr_supported),
        "classification": (
            "simple_global_two_route_supported"
            if overlap and directional and fdr_supported
            else "simple_global_two_route_not_supported"
        ),
        "claim_boundary": "Cross-sectional functional-channel disruption is not temporal pollinator decline or effective-service loss.",
    }
    return results, support, decision


@app.command("run")
def run(
    syndrome_scores_csv: Path = typer.Option(..., exists=True),
    channel_observations_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    results, support, decision = run_bridge(
        pd.read_csv(syndrome_scores_csv),
        pd.read_csv(channel_observations_csv),
        pd.read_csv(covariates_csv),
        config,
        evidence_scope=evidence_scope,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "multichannel_two_route_models.csv", index=False)
    support.to_csv(output_dir / "multichannel_two_route_support.csv", index=False)
    (output_dir / "multichannel_two_route_decision.json").write_text(
        json.dumps(decision, indent=2) + "\n", encoding="utf-8"
    )
    primary = results.loc[
        results["stratum"].eq(str(config["models"]["primary_flora_scope"]))
        & results["min_evaluable_channels"].eq(int(config["primary_decision"]["min_evaluable_channels"]))
        & results["context"].isin([str(x) for x in config["models"]["primary_contexts"]])
        & results["analysis_role"].isin(
            ["route_A_reproductive_assurance", "route_B_accessibility_conditional_on_selfing"]
        )
    ]
    lines = [
        f"# Multichannel two-route bridge — {evidence_scope}",
        "",
        f"- classification: **{decision['classification']}**",
        f"- symmetric state overlap passed: **{decision['symmetric_state_overlap_passed']}**",
        f"- directional requirement passed: **{decision['two_route_directional_requirement_passed']}**",
        f"- four primary tests FDR supported: **{decision['primary_four_tests_fdr_supported']}**",
        "",
        "| context | route | estimate | p | q | n | disrupted |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in primary.iterrows():
        lines.append(
            "| {context} | {role} | {estimate:.5f} | {p:.5g} | {q:.5g} | {n} | {disrupted} |".format(
                context=row["context"],
                role=row["analysis_role"],
                estimate=float(row.get("disruption_estimate", float("nan"))),
                p=float(row.get("p_value", float("nan"))),
                q=float(row.get("q_value_primary_family", float("nan"))),
                n=int(row.get("n_islands", 0)),
                disrupted=int(row.get("n_disrupted", 0)),
            )
        )
    lines.extend(
        [
            "",
            "> The exposure is conservative source-to-island functional-channel disruption, not temporal decline or measured pollination service.",
        ]
    )
    (output_dir / "RESULT_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    app()
