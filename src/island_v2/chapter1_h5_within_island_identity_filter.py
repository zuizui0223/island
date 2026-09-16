"""Within-island identity-preserving H5 association test.

The unit is island x identity-matched pollinator channel. Only islands that contain
both a strict retained and a strict disrupted matched channel are informative. Island
fixed effects remove all island-level confounding; channel x context fixed effects
remove static differences among syndrome templates and channel identities.

The contract is frozen before focal-outcome inspection. This remains an occurrence-
state association test and cannot identify effective pollination service or causation.
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

from island_v2.chapter1_all_data_probability import _chi_square_sf_integer_df
from island_v2.chapter1_h5_identity_matched_two_route import strict_channel_rows

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _normal_two_sided_p(z: float) -> float:
    return math.erfc(abs(float(z)) / math.sqrt(2.0)) if math.isfinite(z) else float("nan")


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != "chapter1_h5_within_island_identity_filter_v1":
        raise typer.BadParameter("unexpected within-island identity-filter contract")
    return config


def _strict_config(config: dict[str, Any]) -> dict[str, Any]:
    return {
        "channels": list(config["matched_channels"]),
        "symmetric_effort_gate": config["strict_channel_state"]["symmetric_effort_gate"],
    }


def _plant_wide(scores: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    required = {"island_id", "stratum", "syndrome", "syndrome_score"}
    if missing := required - set(scores.columns):
        raise typer.BadParameter(f"syndrome score table missing columns: {sorted(missing)}")
    part = scores.loc[scores["stratum"].astype(str).eq(str(config["flora_scope"]))].copy()
    wide = part.pivot_table(index="island_id", columns="syndrome", values="syndrome_score", aggfunc="first")
    axes = sorted(set(str(x) for x in config["matched_channels"].values()))
    if missing := set(axes) - set(wide.columns):
        raise typer.BadParameter(f"syndrome score table lacks required axes: {sorted(missing)}")
    out = wide[axes].apply(pd.to_numeric, errors="coerce").dropna().reset_index()
    return out


def prepare_rows(
    scores: pd.DataFrame,
    observations: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    strict = strict_channel_rows(observations, _strict_config(config))
    plant = _plant_wide(scores, config)
    elig = config["eligibility"]
    context_column = str(elig["context_column"])
    cluster_column = str(elig["cluster_column"])
    contexts = [str(x) for x in elig["primary_contexts"]]
    required_cov = {"island_id", context_column, cluster_column}
    if missing := required_cov - set(covariates.columns):
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")

    cov = covariates[["island_id", context_column, cluster_column]].copy()
    multiplicity = cov.groupby("island_id", dropna=False)[[context_column, cluster_column]].nunique(dropna=False)
    conflicts = multiplicity.gt(1).any(axis=1)
    if bool(conflicts.any()):
        examples = [str(x) for x in conflicts.index[conflicts][:5]]
        raise typer.BadParameter(f"conflicting covariate rows for island_id: {examples}")
    cov = cov.drop_duplicates("island_id")

    data = strict.merge(cov, on="island_id", how="left").merge(plant, on="island_id", how="inner")
    data = data.loc[data[context_column].astype(str).isin(contexts)].copy()
    data = data.loc[data[cluster_column].fillna("").astype(str).ne("")].copy()

    min_channels = int(elig["minimum_evaluable_matched_channels_per_island"])
    audit_rows: list[dict[str, Any]] = []
    keep: list[str] = []
    for island_id, group in data.groupby("island_id", sort=False):
        n_channels = int(group["channel_id"].nunique())
        states = set(group["strict_state"].astype(str))
        mixed = {"retained", "disrupted"}.issubset(states)
        context_values = sorted(set(group[context_column].astype(str)))
        block_values = sorted(set(group[cluster_column].astype(str)))
        eligible = n_channels >= min_channels and mixed and len(context_values) == 1 and len(block_values) == 1
        audit_rows.append(
            {
                "island_id": island_id,
                "context": context_values[0] if len(context_values) == 1 else "|".join(context_values),
                "n_matched_channels": n_channels,
                "n_retained": int(group["strict_state"].eq("retained").sum()),
                "n_disrupted": int(group["strict_state"].eq("disrupted").sum()),
                "eligible_mixed_island": bool(eligible),
            }
        )
        if eligible:
            keep.append(str(island_id))
    data = data.loc[data["island_id"].astype(str).isin(keep)].copy().reset_index(drop=True)
    return data, pd.DataFrame(audit_rows)


def _response_for_mapping(data: pd.DataFrame, mapping: dict[str, str]) -> np.ndarray:
    response = []
    for row in data.itertuples(index=False):
        channel = str(row.channel_id)
        axis = str(mapping[channel])
        response.append(float(getattr(row, axis)))
    return np.asarray(response, dtype=float)


def _fit_mapping(
    data: pd.DataFrame,
    config: dict[str, Any],
    *,
    mapping_name: str,
    mapping: dict[str, str],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    elig = config["eligibility"]
    contexts = [str(x) for x in elig["primary_contexts"]]
    context_column = str(elig["context_column"])
    cluster_column = str(elig["cluster_column"])
    minimum = int(elig["minimum_mixed_islands_per_context"])

    counts = data.groupby(context_column)["island_id"].nunique().to_dict()
    support_passed = all(int(counts.get(context, 0)) >= minimum for context in contexts)
    if not support_passed:
        empty = pd.DataFrame(
            [
                {
                    "mapping": mapping_name,
                    "context": context,
                    "n_mixed_islands": int(counts.get(context, 0)),
                    "estimate": np.nan,
                    "cluster_robust_se": np.nan,
                    "z_value": np.nan,
                    "p_value": np.nan,
                    "support_passed": False,
                }
                for context in contexts
            ]
        )
        return empty, {
            "mapping": mapping_name,
            "support_passed": False,
            "joint_df": 0,
            "joint_wald_chisq": np.nan,
            "joint_p_value": np.nan,
        }

    work = data.copy()
    y = _response_for_mapping(work, mapping)
    island_levels = sorted(work["island_id"].astype(str).unique())
    cell = work[context_column].astype(str) + "::" + work["channel_id"].astype(str)
    cell_levels = sorted(cell.unique())

    parts: list[np.ndarray] = []
    names: list[str] = []
    for island in island_levels:
        parts.append(work["island_id"].astype(str).eq(island).to_numpy(float))
        names.append(f"island[{island}]")
    for level in cell_levels[1:]:
        parts.append(cell.eq(level).to_numpy(float))
        names.append(f"cell[{level}]")
    targets: list[str] = []
    for context in contexts:
        name = f"disrupted[{context}]"
        parts.append(
            work["disrupted"].to_numpy(float)
            * work[context_column].astype(str).eq(context).to_numpy(float)
        )
        names.append(name)
        targets.append(name)

    X = np.column_stack(parts)
    island_n = work.groupby("island_id").size()
    weights = work["island_id"].map(lambda x: 1.0 / float(island_n.loc[x])).to_numpy(float)
    sqrt_w = np.sqrt(weights)
    Xw = X * sqrt_w[:, None]
    yw = y * sqrt_w
    bread = np.linalg.pinv(Xw.T @ Xw)
    beta = bread @ Xw.T @ yw
    residual = y - X @ beta

    labels = work[cluster_column].astype(str).to_numpy()
    unique_clusters = np.unique(labels)
    meat = np.zeros((X.shape[1], X.shape[1]), dtype=float)
    for label in unique_clusters:
        mask = labels == label
        score = X[mask].T @ (weights[mask] * residual[mask])
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    n = len(work)
    rank_x = int(np.linalg.matrix_rank(Xw))
    g = len(unique_clusters)
    if g > 1 and n > rank_x:
        covariance *= (g / (g - 1.0)) * ((n - 1.0) / (n - rank_x))

    rows: list[dict[str, Any]] = []
    indices = []
    estimates = []
    for context, target in zip(contexts, targets, strict=True):
        index = names.index(target)
        indices.append(index)
        estimate = float(beta[index])
        estimates.append(estimate)
        stderr = float(math.sqrt(max(float(covariance[index, index]), 0.0)))
        z_value = estimate / stderr if stderr > 0 else float("nan")
        part = work.loc[work[context_column].astype(str).eq(context)]
        rows.append(
            {
                "mapping": mapping_name,
                "context": context,
                "n_channel_rows": int(len(part)),
                "n_mixed_islands": int(part["island_id"].nunique()),
                "n_clusters": int(part[cluster_column].nunique()),
                "estimate": estimate,
                "cluster_robust_se": stderr,
                "z_value": z_value,
                "p_value": _normal_two_sided_p(z_value),
                "support_passed": True,
            }
        )

    target_cov = covariance[np.ix_(indices, indices)]
    joint_rank = int(np.linalg.matrix_rank(target_cov))
    vector = np.asarray(estimates, dtype=float)
    joint = float(vector @ np.linalg.pinv(target_cov) @ vector) if joint_rank > 0 else float("nan")
    joint_p = _chi_square_sf_integer_df(joint, joint_rank) if joint_rank > 0 else float("nan")
    return pd.DataFrame(rows), {
        "mapping": mapping_name,
        "support_passed": True,
        "joint_df": joint_rank,
        "joint_wald_chisq": joint,
        "joint_p_value": joint_p,
        "effect_norm": float(np.linalg.norm(vector)),
    }


def run_analysis(
    scores: pd.DataFrame,
    observations: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
    *,
    evidence_scope: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    data, audit = prepare_rows(scores, observations, covariates, config)
    matched = {str(k): str(v) for k, v in config["matched_channels"].items()}
    mappings: dict[str, dict[str, str]] = {"matched": matched}
    mappings.update(
        {
            str(name): {str(k): str(v) for k, v in mapping.items()}
            for name, mapping in config["falsification_mappings"].items()
            if isinstance(mapping, dict)
        }
    )
    model_parts = []
    omnibus_rows = []
    for name, mapping in mappings.items():
        rows, omnibus = _fit_mapping(data, config, mapping_name=name, mapping=mapping)
        rows.insert(0, "evidence_scope", evidence_scope)
        model_parts.append(rows)
        omnibus["evidence_scope"] = evidence_scope
        omnibus_rows.append(omnibus)
    models = pd.concat(model_parts, ignore_index=True)
    omnibus = pd.DataFrame(omnibus_rows)

    matched_models = models.loc[models["mapping"].eq("matched")]
    matched_omnibus = omnibus.loc[omnibus["mapping"].eq("matched")].iloc[0]
    decision = {
        "contract": config["contract"],
        "evidence_scope": evidence_scope,
        "n_eligible_mixed_islands": int(data["island_id"].nunique()),
        "support_gate_passed": bool(matched_omnibus["support_passed"]),
        "matched_joint_p_value": float(matched_omnibus["joint_p_value"])
        if pd.notna(matched_omnibus["joint_p_value"])
        else None,
        "matched_context_estimates": {
            str(row.context): float(row.estimate) if pd.notna(row.estimate) else None
            for row in matched_models.itertuples()
        },
        "causal_promotion": False,
    }
    return models, omnibus, audit, decision


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
    models, omnibus, audit, decision = run_analysis(
        pd.read_csv(syndrome_scores_csv),
        pd.read_csv(channel_observations_csv),
        pd.read_csv(covariates_csv),
        config,
        evidence_scope=evidence_scope,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    models.to_csv(output_dir / "within_island_identity_models.csv", index=False)
    omnibus.to_csv(output_dir / "within_island_identity_omnibus.csv", index=False)
    audit.to_csv(output_dir / "within_island_identity_support_audit.csv", index=False)
    (output_dir / "within_island_identity_decision.json").write_text(
        json.dumps(decision, indent=2) + "\n", encoding="utf-8"
    )
    lines = [f"# Within-island identity filter — {evidence_scope}", ""]
    for row in omnibus.itertuples():
        lines.append(
            f"- {row.mapping}: support={row.support_passed}, chi2={row.joint_wald_chisq:.6g}, "
            f"df={int(row.joint_df)}, p={row.joint_p_value:.6g}"
            if pd.notna(row.joint_p_value)
            else f"- {row.mapping}: support={row.support_passed}, not testable"
        )
    lines.extend([
        "",
        "Island fixed effects absorb all island-level covariates; occurrence state is not effective service.",
    ])
    (output_dir / "RESULT_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.echo("\n".join(lines))


if __name__ == "__main__":
    app()
