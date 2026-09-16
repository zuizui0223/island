"""Independent Bombus occurrence bridge for Chapter 1 H5.

The canonical exact-island GBIF Search data are independent of focal plant outcomes.
Only the frozen retained/disrupted states are used. Insufficient-effort and unresolved
islands are missing, never absence. Because the Search may stop after a confirmatory
detection, raw Bombus record counts are not treated as abundance.

This is a post-hoc falsification bridge. It can narrow the upstream mechanism claim but
cannot promote a causal pollinator-loss mechanism.
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
    if not isinstance(config, dict) or config.get("contract") != "chapter1_h5_bombus_upstream_bridge_v1":
        raise typer.BadParameter("unexpected Bombus upstream bridge contract")
    return config


def project_bombus_state(observations: pd.DataFrame) -> pd.DataFrame:
    required = {"island_id", "channel_id", "observation_state"}
    if missing := required - set(observations.columns):
        raise typer.BadParameter(f"Bombus Search table missing columns: {sorted(missing)}")
    work = observations.loc[observations["channel_id"].astype(str).eq("bombus")].copy()
    mapping = {"detected": "retained", "adequate_non_detection": "disrupted"}
    work["bombus_state"] = work["observation_state"].astype(str).map(mapping)
    work["bombus_disrupted"] = work["bombus_state"].map({"retained": 0.0, "disrupted": 1.0})
    return work


def build_plant_scores(scores: pd.DataFrame, *, stratum: str) -> pd.DataFrame:
    required = {"island_id", "stratum", "syndrome", "syndrome_score"}
    if missing := required - set(scores.columns):
        raise typer.BadParameter(f"syndrome score table missing columns: {sorted(missing)}")
    part = scores.loc[scores["stratum"].astype(str).eq(stratum)].copy()
    wide = part.pivot_table(
        index="island_id", columns="syndrome", values="syndrome_score", aggfunc="first"
    )
    needed = {"selfing_core", "large_bee_like", "generalized_accessible"}
    if missing := needed - set(wide.columns):
        raise typer.BadParameter(f"syndrome score table lacks required axes: {sorted(missing)}")
    wide["attraction_shift"] = (
        -pd.to_numeric(wide["large_bee_like"], errors="coerce")
        + pd.to_numeric(wide["generalized_accessible"], errors="coerce")
    ) / 2.0
    return wide.reset_index()


def _z(values: pd.Series) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def _fit_clustered_ols(
    data: pd.DataFrame,
    *,
    response: str,
    controls: list[str],
    context_column: str,
    cluster_column: str,
    contexts: list[str],
    condition_on_selfing: bool,
) -> dict[str, Any]:
    columns = [response, "bombus_disrupted", context_column, cluster_column, *controls]
    if condition_on_selfing:
        columns.append("selfing_core")
    work = data.loc[data[context_column].astype(str).isin(contexts), columns].copy()
    numeric = [response, "bombus_disrupted", *controls]
    if condition_on_selfing:
        numeric.append("selfing_core")
    for column in numeric:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.dropna(subset=numeric)
    work = work.loc[
        work[context_column].fillna("").astype(str).ne("")
        & work[cluster_column].fillna("").astype(str).ne("")
    ].copy()
    if len(work) < 20 or work["bombus_disrupted"].nunique() < 2:
        return {
            "status": "not_testable",
            "response": response,
            "condition_on_selfing_core": condition_on_selfing,
            "n_islands": int(len(work)),
            "n_disrupted": int(work["bombus_disrupted"].sum()) if len(work) else 0,
            "n_retained": int((1.0 - work["bombus_disrupted"]).sum()) if len(work) else 0,
        }

    observed_contexts = sorted(work[context_column].astype(str).unique())
    design_parts: list[np.ndarray] = [
        np.ones(len(work), dtype=float),
        work["bombus_disrupted"].to_numpy(float),
    ]
    names = ["intercept", "bombus_disrupted"]
    for value in observed_contexts[1:]:
        design_parts.append(work[context_column].astype(str).eq(value).to_numpy(float))
        names.append(f"context[{value}]")
    for control in controls:
        try:
            design_parts.append(_z(work[control]))
            names.append(f"z_{control}")
        except ValueError:
            continue
    if condition_on_selfing:
        try:
            design_parts.append(_z(work["selfing_core"]))
            names.append("z_selfing_core")
        except ValueError:
            pass

    X = np.column_stack(design_parts)
    y = work[response].to_numpy(float)
    bread = np.linalg.pinv(X.T @ X)
    beta = bread @ (X.T @ y)
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

    index = names.index("bombus_disrupted")
    estimate = float(beta[index])
    stderr = float(math.sqrt(max(float(covariance[index, index]), 0.0)))
    z_value = estimate / stderr if stderr > 0 else float("nan")
    return {
        "status": "fit",
        "response": response,
        "condition_on_selfing_core": condition_on_selfing,
        "contexts": "|".join(contexts),
        "contexts_present": "|".join(observed_contexts),
        "n_islands": int(n),
        "n_clusters": int(g),
        "n_disrupted": int(work["bombus_disrupted"].sum()),
        "n_retained": int((1.0 - work["bombus_disrupted"]).sum()),
        "bombus_disrupted_estimate": estimate,
        "cluster_robust_se": stderr,
        "z_value": float(z_value),
        "p_value": _normal_two_sided_p(z_value),
    }


def state_support(data: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    context = str(config["models"]["context_column"])
    rows: list[dict[str, Any]] = []
    for value, part in data.dropna(subset=["bombus_state"]).groupby(context, sort=True):
        counts = part["bombus_state"].value_counts()
        rows.append(
            {
                "context": str(value),
                "n_retained": int(counts.get("retained", 0)),
                "n_disrupted": int(counts.get("disrupted", 0)),
                "n_evaluable": int(len(part)),
            }
        )
    return pd.DataFrame(rows)


def overlap_gate(support: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    spec = config["support_and_overlap"]
    primary = [str(x) for x in config["models"]["primary_contexts"]]
    minimum = int(spec["context_min_each_state_for_confirmatory_interpretation"])
    rows = support.loc[support["context"].isin(primary)].copy()
    rows["passes_context_overlap"] = rows["n_retained"].ge(minimum) & rows["n_disrupted"].ge(minimum)
    n_pass = int(rows["passes_context_overlap"].sum())
    required = 2 if bool(spec.get("require_two_contexts_meet_context_overlap_gate", False)) else 1
    return {
        "primary_contexts": primary,
        "minimum_each_state_per_context": minimum,
        "n_primary_contexts_passing": n_pass,
        "required_primary_contexts_passing": required,
        "overlap_gate_passed": bool(n_pass >= required),
    }


def fit_bridge(
    syndrome_scores: pd.DataFrame,
    observations: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
    *,
    evidence_scope: str,
    h3_scores: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    plant = build_plant_scores(syndrome_scores, stratum=str(config["plant_input"]["flora_scope"]))
    bombus = project_bombus_state(observations)
    model = config["models"]
    context = str(model["context_column"])
    cluster = str(model["cluster_column"])
    controls = [str(x) for x in model["controls"]]
    required_cov = {"island_id", context, cluster, *controls}
    if missing := required_cov - set(covariates.columns):
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")
    data = (
        bombus.merge(covariates[list(required_cov)].drop_duplicates("island_id"), on="island_id", how="left")
        .merge(plant, on="island_id", how="left")
    )
    support = state_support(data, config)
    gate = overlap_gate(support, config)

    rows: list[dict[str, Any]] = []
    context_sets = [
        [str(x) for x in model["global_contexts"]],
        [str(x) for x in model["primary_contexts"]],
    ]
    for contexts in context_sets:
        for response, conditional in [
            ("selfing_core", False),
            ("attraction_shift", False),
            ("attraction_shift", True),
            ("large_bee_like", False),
            ("generalized_accessible", False),
        ]:
            result = _fit_clustered_ols(
                data,
                response=response,
                controls=controls,
                context_column=context,
                cluster_column=cluster,
                contexts=contexts,
                condition_on_selfing=conditional,
            )
            result["evidence_scope"] = evidence_scope
            result["analysis_role"] = "global_adjusted" if len(contexts) == 4 else "north_tropical_adjusted"
            rows.append(result)

    secondary_rows: list[dict[str, Any]] = []
    if h3_scores is not None and not h3_scores.empty:
        h3 = h3_scores.loc[
            h3_scores["outcome"].astype(str).isin(
                [str(x) for x in config["secondary_taxonomic_check"]["outcomes"]]
            )
        ].copy()
        h3 = h3.rename(columns={str(config["secondary_taxonomic_check"]["stage"]): "h3_response"})
        joined = (
            bombus.merge(covariates[list(required_cov)].drop_duplicates("island_id"), on="island_id", how="left")
            .merge(h3[["island_id", "outcome", "h3_response"]], on="island_id", how="inner")
        )
        for outcome, part in joined.groupby("outcome", sort=True):
            for contexts in context_sets:
                result = _fit_clustered_ols(
                    part,
                    response="h3_response",
                    controls=controls,
                    context_column=context,
                    cluster_column=cluster,
                    contexts=contexts,
                    condition_on_selfing=False,
                )
                result["evidence_scope"] = evidence_scope
                result["h3_outcome"] = str(outcome)
                result["analysis_role"] = "post_genus_residual"
                secondary_rows.append(result)

    results = pd.DataFrame(rows)
    secondary = pd.DataFrame(secondary_rows)
    primary_rows = results.loc[
        results["analysis_role"].eq("north_tropical_adjusted")
        & (
            results["response"].eq("selfing_core")
            | (results["response"].eq("attraction_shift") & results["condition_on_selfing_core"].eq(True))
        )
    ].copy()
    directional = bool(
        len(primary_rows) == 2
        and primary_rows["status"].eq("fit").all()
        and primary_rows["bombus_disrupted_estimate"].gt(0).all()
    )
    significant = bool(
        directional
        and primary_rows["p_value"].le(0.05).all()
    )
    decision = {
        "contract": config["contract"],
        "evidence_scope": evidence_scope,
        **gate,
        "parallel_pathway_directional_consistency": directional,
        "both_primary_coefficients_p_le_0_05": significant,
        "bombus_upstream_mechanism_promoted": False,
        "classification": (
            "overlap_gate_failed"
            if not gate["overlap_gate_passed"]
            else "posthoc_bridge_not_promotable_by_contract"
        ),
        "claim": (
            "Canonical exact-island Bombus occurrence evidence does not close the upstream causal arrow."
        ),
    }
    return results, secondary, support, decision


@app.command("run")
def run(
    syndrome_scores_csv: Path = typer.Option(..., exists=True),
    bombus_observations_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
    h3_scores_csv: Path | None = typer.Option(None),
) -> None:
    config = load_config(config_path)
    results, secondary, support, decision = fit_bridge(
        pd.read_csv(syndrome_scores_csv),
        pd.read_csv(bombus_observations_csv),
        pd.read_csv(covariates_csv),
        config,
        evidence_scope=evidence_scope,
        h3_scores=pd.read_csv(h3_scores_csv) if h3_scores_csv else None,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "bombus_upstream_bridge_models.csv", index=False)
    secondary.to_csv(output_dir / "bombus_upstream_post_genus_models.csv", index=False)
    support.to_csv(output_dir / "bombus_upstream_state_support.csv", index=False)
    (output_dir / "bombus_upstream_bridge_decision.json").write_text(
        json.dumps(decision, indent=2) + "\n", encoding="utf-8"
    )
    primary = results.loc[results["analysis_role"].eq("north_tropical_adjusted")]
    lines = [
        f"# Bombus upstream bridge — {evidence_scope}",
        "",
        f"- overlap gate: **{decision['overlap_gate_passed']}**",
        f"- classification: **{decision['classification']}**",
        "- exact-island Search states: retained = detected; disrupted = source-available + adequate non-detection.",
        "- insufficient-effort and unresolved islands are missing, never absence.",
        "",
        "## North–Tropical adjusted diagnostic",
        "",
    ]
    for row in primary.to_dict("records"):
        if row.get("status") != "fit":
            continue
        suffix = " conditional on selfing_core" if row.get("condition_on_selfing_core") else ""
        lines.append(
            f"- {row['response']}{suffix}: beta_disrupted={row['bombus_disrupted_estimate']:.6g}, "
            f"SE={row['cluster_robust_se']:.6g}, p={row['p_value']:.6g}, "
            f"n={int(row['n_islands'])} ({int(row['n_disrupted'])} disrupted)."
        )
    lines.extend([
        "",
        "> This post-hoc bridge cannot prove historical pollinator decline or direct pollinator-mediated selection.",
        "",
    ])
    (output_dir / "RESULT_SUMMARY.md").write_text("\n".join(lines), encoding="utf-8")
    typer.echo(json.dumps(decision, indent=2))


if __name__ == "__main__":
    app()
