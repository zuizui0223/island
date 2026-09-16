"""Identity-aware pollinator-disruption diagnostic for Chapter 1 H5.

The unit is island x pollinator channel.  Only strict effort-qualified retained or
source-available disrupted states enter.  Islands are given equal total weight across
all channel rows.  Channel x biogeographic-context fixed effects prevent the common
disruption coefficient from being estimated by static differences among channel
identities or regions.

This module is post hoc.  It can narrow mechanism interpretations but cannot promote
pollinator causation.
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


def _bh_qvalues(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce")
    valid = p.dropna().sort_values()
    out = pd.Series(np.nan, index=values.index, dtype=float)
    if valid.empty:
        return out
    m = len(valid)
    raw = valid.to_numpy(float) * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(raw[::-1])[::-1]
    out.loc[valid.index] = np.clip(adjusted, 0.0, 1.0)
    return out


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != "chapter1_h5_identity_matched_two_route_v1":
        raise typer.BadParameter("unexpected identity-matched two-route contract")
    return config


def strict_channel_rows(observations: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
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
    if work.duplicated(["island_id", "channel_id"]).any():
        raise typer.BadParameter("channel observations must be unique by island_id x channel_id")
    numeric = [
        "background_record_count",
        "background_spatial_units",
        "background_temporal_units",
        "distinct_dataset_count",
        "latest_background_year",
    ]
    for column in numeric:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    gate = config["symmetric_effort_gate"]
    min_year = int(gate["reference_year"]) - int(gate["max_years_since_latest_background"])
    qualified = (
        work["background_record_count"].ge(float(gate["min_background_records"]))
        & work["background_spatial_units"].ge(float(gate["min_background_spatial_units"]))
        & work["background_temporal_units"].ge(float(gate["min_background_temporal_units"]))
        & work["distinct_dataset_count"].ge(float(gate["min_distinct_datasets"]))
        & work["latest_background_year"].ge(float(min_year))
    )
    retained = work["observation_state"].astype(str).eq("detected") & qualified
    disrupted = work["observation_state"].astype(str).eq("adequate_non_detection")
    work["strict_state"] = np.where(retained, "retained", np.where(disrupted, "disrupted", "missing"))
    work = work.loc[work["strict_state"].isin(["retained", "disrupted"])].copy()
    work["disrupted"] = work["strict_state"].eq("disrupted").astype(float)
    return work


def plant_scores(scores: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    required = {"island_id", "stratum", "syndrome", "syndrome_score"}
    if missing := required - set(scores.columns):
        raise typer.BadParameter(f"syndrome score table missing columns: {sorted(missing)}")
    scope = str(config["models"]["flora_scope"])
    part = scores.loc[scores["stratum"].astype(str).eq(scope)].copy()
    wide = part.pivot_table(index="island_id", columns="syndrome", values="syndrome_score", aggfunc="first")
    needed = {"selfing_core", "large_bee_like", "butterfly_like", "bird_like"}
    if missing := needed - set(wide.columns):
        raise typer.BadParameter(f"syndrome score table lacks required axes: {sorted(missing)}")
    return wide.reset_index()


def _z(values: pd.Series) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    sd = float(np.std(x, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - float(np.mean(x))) / sd


def _weighted_clustered_context_effects(
    data: pd.DataFrame,
    *,
    response: str,
    controls: list[str],
    context_column: str,
    cluster_column: str,
    primary_contexts: list[str],
    condition_on_selfing: bool,
) -> list[dict[str, Any]]:
    columns = [
        response,
        "disrupted",
        "island_id",
        "channel_id",
        context_column,
        cluster_column,
        *controls,
    ]
    if condition_on_selfing and response != "selfing_core":
        columns.append("selfing_core")
    work = data.loc[data[context_column].astype(str).isin(primary_contexts), columns].copy()
    numeric = [response, "disrupted", *controls]
    if condition_on_selfing and response != "selfing_core":
        numeric.append("selfing_core")
    for column in numeric:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.dropna(subset=numeric)
    work = work.loc[work[cluster_column].fillna("").astype(str).ne("")].copy()
    if len(work) < 20:
        return []

    island_n = work.groupby("island_id").size()
    weights = work["island_id"].map(lambda x: 1.0 / float(island_n.loc[x])).to_numpy(float)

    parts: list[np.ndarray] = [np.ones(len(work), dtype=float)]
    names = ["intercept"]
    cell = work[context_column].astype(str) + "::" + work["channel_id"].astype(str)
    levels = sorted(cell.unique())
    for level in levels[1:]:
        parts.append(cell.eq(level).to_numpy(float))
        names.append(f"cell[{level}]")

    for context in primary_contexts:
        parts.append(
            work["disrupted"].to_numpy(float)
            * work[context_column].astype(str).eq(context).to_numpy(float)
        )
        names.append(f"disrupted[{context}]")

    for control in controls:
        try:
            parts.append(_z(work[control]))
            names.append(f"z_{control}")
        except ValueError:
            continue
    if condition_on_selfing and response != "selfing_core":
        try:
            parts.append(_z(work["selfing_core"]))
            names.append("z_selfing_core")
        except ValueError:
            pass

    X = np.column_stack(parts)
    y = work[response].to_numpy(float)
    sqrt_w = np.sqrt(weights)
    Xw = X * sqrt_w[:, None]
    yw = y * sqrt_w
    bread = np.linalg.pinv(Xw.T @ Xw)
    beta = bread @ Xw.T @ yw
    residual = y - X @ beta

    labels = work[cluster_column].astype(str).to_numpy()
    unique = np.unique(labels)
    meat = np.zeros((X.shape[1], X.shape[1]), dtype=float)
    for label in unique:
        mask = labels == label
        score = X[mask].T @ (weights[mask] * residual[mask])
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    n = len(work)
    k = X.shape[1]
    g = len(unique)
    if g > 1 and n > k:
        covariance *= (g / (g - 1.0)) * ((n - 1.0) / (n - k))

    rows: list[dict[str, Any]] = []
    for context in primary_contexts:
        name = f"disrupted[{context}]"
        index = names.index(name)
        estimate = float(beta[index])
        stderr = float(math.sqrt(max(float(covariance[index, index]), 0.0)))
        z_value = estimate / stderr if stderr > 0 else float("nan")
        part = work.loc[work[context_column].astype(str).eq(context)]
        rows.append(
            {
                "context": context,
                "response": response,
                "condition_on_selfing_core": bool(condition_on_selfing),
                "n_channel_rows": int(len(part)),
                "n_islands": int(part["island_id"].nunique()),
                "n_clusters": int(part[cluster_column].nunique()),
                "n_disrupted_rows": int(part["disrupted"].sum()),
                "estimate": estimate,
                "cluster_robust_se": stderr,
                "z_value": float(z_value),
                "p_value": _normal_two_sided_p(z_value),
            }
        )
    return rows


def run_diagnostic(
    scores: pd.DataFrame,
    observations: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
    *,
    evidence_scope: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    strict = strict_channel_rows(observations, config)
    plant = plant_scores(scores, config)
    model = config["models"]
    context_column = str(model["context_column"])
    cluster_column = str(model["cluster_column"])
    controls = [str(x) for x in model["controls"]]
    contexts = [str(x) for x in model["primary_contexts"]]
    required_cov = {"island_id", context_column, cluster_column, *controls}
    if missing := required_cov - set(covariates.columns):
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")
    data = (
        strict.merge(covariates[list(required_cov)].drop_duplicates("island_id"), on="island_id", how="left")
        .merge(plant, on="island_id", how="left")
    )

    support = (
        data.loc[data[context_column].astype(str).isin(contexts)]
        .groupby([context_column, "channel_id", "strict_state"], dropna=False)
        .size()
        .rename("n")
        .reset_index()
    )

    rows: list[dict[str, Any]] = []
    route_a = _weighted_clustered_context_effects(
        data,
        response="selfing_core",
        controls=controls,
        context_column=context_column,
        cluster_column=cluster_column,
        primary_contexts=contexts,
        condition_on_selfing=False,
    )
    for row in route_a:
        row["route"] = "A_reproductive_assurance"
        row["evidence_scope"] = evidence_scope
        rows.append(row)

    mapping = {str(k): str(v) for k, v in config["identity_matched_architecture"].items() if isinstance(v, str)}
    matched = data.loc[data["channel_id"].astype(str).isin(mapping)].copy()
    matched["channel_matched_architecture"] = [
        row[mapping[str(row["channel_id"])]] for _, row in matched.iterrows()
    ]
    route_b = _weighted_clustered_context_effects(
        matched,
        response="channel_matched_architecture",
        controls=controls,
        context_column=context_column,
        cluster_column=cluster_column,
        primary_contexts=contexts,
        condition_on_selfing=True,
    )
    for row in route_b:
        row["route"] = "B_identity_matched_architecture"
        row["evidence_scope"] = evidence_scope
        rows.append(row)

    results = pd.DataFrame(rows)
    results["q_value_primary_family"] = _bh_qvalues(results["p_value"])
    decision = {
        "contract": config["contract"],
        "evidence_scope": evidence_scope,
        "n_primary_tests": int(len(results)),
        "any_primary_fdr_supported": bool(results["q_value_primary_family"].le(0.05).fillna(False).any()),
        "causal_promotion": False,
        "classification": "posthoc_identity_aware_diagnostic_only",
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
    results, support, decision = run_diagnostic(
        pd.read_csv(syndrome_scores_csv),
        pd.read_csv(channel_observations_csv),
        pd.read_csv(covariates_csv),
        config,
        evidence_scope=evidence_scope,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "identity_matched_two_route_models.csv", index=False)
    support.to_csv(output_dir / "identity_matched_channel_support.csv", index=False)
    (output_dir / "identity_matched_two_route_decision.json").write_text(
        json.dumps(decision, indent=2) + "\n", encoding="utf-8"
    )
    lines = [f"# Identity-matched two-route diagnostic — {evidence_scope}", ""]
    for _, row in results.iterrows():
        lines.append(
            f"- {row['route']} / {row['context']}: beta={row['estimate']:.6g}, "
            f"p={row['p_value']:.6g}, q={row['q_value_primary_family']:.6g}, "
            f"islands={int(row['n_islands'])}, disrupted channel rows={int(row['n_disrupted_rows'])}"
        )
    lines.extend([
        "",
        "Post-hoc diagnostic only; strict channel states are not temporal abundance or service decline.",
    ])
    (output_dir / "RESULT_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    app()
