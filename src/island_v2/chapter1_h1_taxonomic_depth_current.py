"""Taxonomic representation depth of the current seven-response Chapter 1 H1.

This is a representation-depth audit, not a causal mechanism test. It uses the same
species support at observed, family-residual and genus-residual stages. Family/genus
expectations are leave-one-species-out means from the fixed scored species pool and
never impute missing trait states.

A beta-binomial H1 gate is first re-fit on the exact taxonomically eligible species
support. Taxonomic attenuation is interpreted only where that common-support H1
remains supported.
"""
from __future__ import annotations

import json
import math
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_all_data_probability import (
    _bh,
    _chi_square_sf_integer_df,
    build_broad_counts,
    run_probability_analysis,
)
from island_v2.chapter1_h3_observed_taxonomic_depth import (
    STAGES,
    _assemble_cluster_covariance,
    _fit_ols_component,
    _z,
    build_atomic_taxonomic_residuals,
    build_island_stage_scores,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def load_config(path: Path) -> dict[str, Any]:
    cfg = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(cfg, dict) or cfg.get("contract") != "chapter1_h1_taxonomic_depth_current_v1":
        raise typer.BadParameter("unexpected current H1 taxonomic-depth contract")
    return cfg


def _common_support_audit(
    state_audit: pd.DataFrame,
    species_residuals: pd.DataFrame,
    probability_config: dict[str, Any],
) -> pd.DataFrame:
    audit = state_audit.copy()
    audit["accepted_species"] = audit["accepted_species"].astype(str)
    audit["trait_name"] = audit["trait_name"].astype(str)
    pieces: list[pd.DataFrame] = []
    for outcome in probability_config["model_outcomes"]:
        spec = probability_config["broad_outcomes"][outcome]
        eligible = set(
            species_residuals.loc[
                species_residuals["outcome"].astype(str).eq(str(outcome)),
                "accepted_species",
            ].astype(str)
        )
        part = audit.loc[
            audit["trait_name"].eq(str(spec["trait_name"]))
            & audit["accepted_species"].isin(eligible)
        ].copy()
        if not part.empty:
            pieces.append(part)
    if not pieces:
        return audit.iloc[0:0].copy()
    return pd.concat(pieces, ignore_index=True)


def _prepare_model_data(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    cfg: dict[str, Any],
) -> pd.DataFrame:
    context = str(cfg["context_column"])
    cluster = str(cfg["cluster_column"])
    distance = str(cfg["distance_column"])
    controls = [str(x) for x in cfg["controls"]]
    needed = ["island_id", context, cluster, distance, *controls]
    missing = set(needed) - set(covariates.columns)
    if missing:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")
    data = island_scores.merge(
        covariates[needed].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    for column in [*STAGES, "n_species", distance, *controls]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data[context] = data[context].fillna("").astype(str)
    data[cluster] = data[cluster].fillna("").astype(str)
    data = data.dropna(subset=[*STAGES, distance, *controls])
    return data.loc[data[context].isin([str(x) for x in cfg["contexts"]]) & data[cluster].ne("")].copy()


def fit_within_stage(
    data: pd.DataFrame,
    *,
    context_value: str,
    stage: str,
    cfg: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any], np.ndarray]:
    context = str(cfg["context_column"])
    cluster = str(cfg["cluster_column"])
    distance = str(cfg["distance_column"])
    controls = [str(x) for x in cfg["controls"]]
    threshold = int(cfg["island_model"]["minimum_islands_per_context_outcome"])
    outcomes = [str(x) for x in cfg["atomic_outcomes"]]

    work = data.loc[data[context].eq(context_value)].copy()
    counts = work.groupby("outcome")["island_id"].nunique()
    retained = [x for x in outcomes if int(counts.get(x, 0)) >= threshold]
    if len(retained) < int(cfg["island_model"]["minimum_outcomes_per_vector"]):
        return pd.DataFrame(), {
            "context": context_value,
            "stage": stage,
            "status": "not_testable",
            "n_retained_outcomes": len(retained),
            "retained_outcomes": "|".join(retained),
        }, np.array([], dtype=float)

    fits: list[dict[str, Any]] = []
    indices: list[int] = []
    rows: list[dict[str, Any]] = []
    offset = 0
    for outcome in retained:
        part = work.loc[work["outcome"].eq(outcome)].copy()
        columns = [np.ones(len(part), dtype=float)]
        names = [f"{outcome}:intercept"]
        for predictor in controls:
            columns.append(_z(part[predictor]))
            names.append(f"{outcome}:z_{predictor}")
        columns.append(_z(part[distance]))
        distance_name = f"{outcome}:z_{distance}"
        names.append(distance_name)
        fit = _fit_ols_component(
            part[stage].to_numpy(float),
            np.column_stack(columns),
            names,
            part[cluster].to_numpy(str),
        )
        fits.append(fit)
        indices.append(offset + names.index(distance_name))
        rows.append(
            {
                "context": context_value,
                "stage": stage,
                "outcome": outcome,
                "n_islands": int(part["island_id"].nunique()),
                "n_clusters": int(part[cluster].nunique()),
                "median_species_support": float(part["n_species"].median()),
            }
        )
        offset += len(names)

    theta, covariance, _ = _assemble_cluster_covariance(fits)
    vector = theta[indices]
    vector_cov = covariance[np.ix_(indices, indices)]
    se = np.sqrt(np.clip(np.diag(vector_cov), 0.0, None))
    for row, estimate, stderr in zip(rows, vector, se, strict=True):
        z = float(estimate / stderr) if stderr > 0 else float("nan")
        row.update(
            isolation_slope=float(estimate),
            cluster_robust_se=float(stderr),
            p_value=math.erfc(abs(z) / math.sqrt(2.0)) if math.isfinite(z) else float("nan"),
        )

    rank = int(np.linalg.matrix_rank(vector_cov))
    statistic = float(vector @ np.linalg.pinv(vector_cov) @ vector) if rank > 0 else float("nan")
    omnibus = {
        "context": context_value,
        "stage": stage,
        "status": "fit",
        "n_retained_outcomes": len(retained),
        "retained_outcomes": "|".join(retained),
        "n_unique_islands": int(work.loc[work["outcome"].isin(retained), "island_id"].nunique()),
        "n_clusters": int(work.loc[work["outcome"].isin(retained), cluster].nunique()),
        "joint_wald_chisq": statistic,
        "joint_df": rank,
        "p_value": _chi_square_sf_integer_df(statistic, rank),
        "vector_norm": float(np.linalg.norm(vector)),
    }
    return pd.DataFrame(rows), omnibus, vector



def _fit_stage_vector_point(
    data: pd.DataFrame,
    *,
    context_value: str,
    stage: str,
    retained: list[str],
    cfg: dict[str, Any],
) -> np.ndarray:
    context = str(cfg["context_column"])
    distance = str(cfg["distance_column"])
    controls = [str(x) for x in cfg["controls"]]
    work = data.loc[data[context].eq(context_value)].copy()
    values: list[float] = []
    for outcome in retained:
        part = work.loc[work["outcome"].eq(outcome)].copy()
        if part.empty:
            return np.array([], dtype=float)
        columns = [np.ones(len(part), dtype=float)]
        try:
            for predictor in controls:
                columns.append(_z(part[predictor]))
            columns.append(_z(part[distance]))
        except ValueError:
            return np.array([], dtype=float)
        X = np.column_stack(columns)
        y = part[stage].to_numpy(float)
        beta = np.linalg.pinv(X.T @ X) @ (X.T @ y)
        values.append(float(beta[-1]))
    return np.asarray(values, dtype=float)


def _bootstrap_attenuation(
    data: pd.DataFrame,
    stage_omnibus: pd.DataFrame,
    cfg: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    spec = cfg.get("attenuation_bootstrap", {})
    if not bool(spec.get("enabled", False)):
        return pd.DataFrame(), pd.DataFrame()
    context_col = str(cfg["context_column"])
    cluster_col = str(cfg["cluster_column"])
    draws = int(spec["draws"])
    seed = int(spec["seed"])
    rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []

    for context_index, context_value in enumerate([str(x) for x in cfg["contexts"]]):
        observed_row = stage_omnibus.loc[
            stage_omnibus["context"].eq(context_value)
            & stage_omnibus["stage"].eq("observed_score")
            & stage_omnibus["status"].eq("fit")
        ]
        if observed_row.empty:
            continue
        retained = [
            x for x in str(observed_row.iloc[0]["retained_outcomes"]).split("|") if x
        ]
        part = data.loc[data[context_col].eq(context_value)].copy()
        labels = sorted(part[cluster_col].astype(str).unique())
        if len(labels) < 2:
            continue
        by_cluster = {
            label: part.loc[part[cluster_col].astype(str).eq(label)].copy()
            for label in labels
        }
        rng = np.random.default_rng(seed + context_index)
        for draw in range(draws):
            sampled = rng.choice(labels, size=len(labels), replace=True)
            pieces: list[pd.DataFrame] = []
            for replicate, label in enumerate(sampled):
                frame = by_cluster[str(label)].copy()
                frame[cluster_col] = f"{label}__boot{replicate}"
                pieces.append(frame)
            boot = pd.concat(pieces, ignore_index=True)
            vectors: dict[str, np.ndarray] = {}
            valid = True
            for stage in STAGES:
                vector = _fit_stage_vector_point(
                    boot,
                    context_value=context_value,
                    stage=stage,
                    retained=retained,
                    cfg=cfg,
                )
                if len(vector) != len(retained):
                    valid = False
                    break
                vectors[stage] = vector
            if not valid:
                continue
            observed_norm = float(np.linalg.norm(vectors["observed_score"]))
            family_norm = float(np.linalg.norm(vectors["after_family_residual"]))
            genus_norm = float(np.linalg.norm(vectors["after_genus_residual"]))
            if observed_norm <= 0 or family_norm <= 0:
                continue
            rows.append(
                {
                    "context": context_value,
                    "draw": draw,
                    "observed_norm": observed_norm,
                    "family_norm": family_norm,
                    "genus_norm": genus_norm,
                    "family_total_attenuation": 1.0 - family_norm / observed_norm,
                    "genus_total_attenuation": 1.0 - genus_norm / observed_norm,
                    "incremental_genus_attenuation": 1.0 - genus_norm / family_norm,
                }
            )

    draws_frame = pd.DataFrame(rows)
    if draws_frame.empty:
        return draws_frame, pd.DataFrame()
    for context_value, part in draws_frame.groupby("context", sort=False):
        item: dict[str, Any] = {
            "context": context_value,
            "valid_draws": int(len(part)),
        }
        for column in (
            "family_total_attenuation",
            "genus_total_attenuation",
            "incremental_genus_attenuation",
        ):
            item[f"{column}_median"] = float(part[column].median())
            item[f"{column}_ci_low"] = float(part[column].quantile(0.025))
            item[f"{column}_ci_high"] = float(part[column].quantile(0.975))
        summaries.append(item)
    return draws_frame, pd.DataFrame(summaries)


def _run_common_support_h1(
    flora: pd.DataFrame,
    common_audit: pd.DataFrame,
    covariates: pd.DataFrame,
    probability_config: dict[str, Any],
    *,
    flora_scope: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cfg = deepcopy(probability_config)
    cfg["strata"] = ["all_observed"]
    counts = build_broad_counts(flora, common_audit, cfg)
    counts["stratum"] = flora_scope
    cfg["strata"] = [flora_scope]
    slopes, _, omnibus, _ = run_probability_analysis(counts, covariates, cfg)
    return slopes, omnibus


def _classify_context(
    gate_row: pd.Series | None,
    stage_rows: pd.DataFrame,
    *,
    alpha: float,
) -> str:
    if gate_row is None or not bool(gate_row.get("vector_supported", False)):
        return "common_support_H1_not_supported"
    by_stage = stage_rows.set_index("stage") if not stage_rows.empty else pd.DataFrame()
    if any(stage not in by_stage.index for stage in STAGES):
        return "taxonomic_decomposition_not_testable"
    supported = {
        stage: bool(float(by_stage.loc[stage, "q_value"]) <= alpha)
        for stage in STAGES
        if math.isfinite(float(by_stage.loc[stage, "q_value"]))
    }
    if not supported.get("observed_score", False):
        return "linear_decomposition_observed_stage_not_supported"
    if not supported.get("after_family_residual", False):
        return "compatible_with_family_level_structuring"
    if not supported.get("after_genus_residual", False):
        return "compatible_with_genus_level_structuring"
    return "below_genus_or_non_taxonomic_residual_retained"


def run_analysis(
    state_audit: pd.DataFrame,
    taxonomy: pd.DataFrame,
    flora: pd.DataFrame,
    covariates: pd.DataFrame,
    probability_config: dict[str, Any],
    depth_config: dict[str, Any],
    *,
    evidence_scope: str,
    flora_scope: str,
) -> dict[str, Any]:
    residuals = build_atomic_taxonomic_residuals(
        state_audit, taxonomy, probability_config, depth_config
    )
    common_audit = _common_support_audit(state_audit, residuals, probability_config)
    gate_slopes, gate = _run_common_support_h1(
        flora,
        common_audit,
        covariates,
        probability_config,
        flora_scope=flora_scope,
    )

    island_scores = build_island_stage_scores(flora, residuals)
    data = _prepare_model_data(island_scores, covariates, depth_config)

    slope_parts: list[pd.DataFrame] = []
    omnibus_rows: list[dict[str, Any]] = []
    for stage in STAGES:
        for context_value in [str(x) for x in depth_config["contexts"]]:
            slopes, result, _ = fit_within_stage(
                data,
                context_value=context_value,
                stage=stage,
                cfg=depth_config,
            )
            if not slopes.empty:
                slope_parts.append(slopes)
            omnibus_rows.append(result)

    stage_slopes = pd.concat(slope_parts, ignore_index=True) if slope_parts else pd.DataFrame()
    stage_omnibus = pd.DataFrame(omnibus_rows)
    if "p_value" in stage_omnibus.columns:
        stage_omnibus["q_value"] = stage_omnibus.groupby("stage", group_keys=False)[
            "p_value"
        ].transform(_bh)
        stage_omnibus["vector_supported"] = (
            stage_omnibus["q_value"]
            .le(float(depth_config["multiplicity"]["alpha"]))
            .fillna(False)
        )

    bootstrap_draws, bootstrap_summary = _bootstrap_attenuation(
        data, stage_omnibus, depth_config
    )

    attenuation_rows: list[dict[str, Any]] = []
    classifications: list[dict[str, Any]] = []
    gate_index = gate.set_index("context") if not gate.empty and "context" in gate.columns else None
    for context_value in [str(x) for x in depth_config["contexts"]]:
        part = stage_omnibus.loc[stage_omnibus["context"].eq(context_value)].copy()
        norms = {
            str(row.stage): float(row.vector_norm)
            for row in part.itertuples()
            if getattr(row, "status", "") == "fit" and math.isfinite(float(row.vector_norm))
        }
        observed = norms.get("observed_score", float("nan"))
        attenuation_rows.append(
            {
                "context": context_value,
                "observed_norm": observed,
                "family_norm": norms.get("after_family_residual", float("nan")),
                "genus_norm": norms.get("after_genus_residual", float("nan")),
                "family_attenuation": (
                    1.0 - norms["after_family_residual"] / observed
                    if observed > 0 and "after_family_residual" in norms
                    else float("nan")
                ),
                "genus_attenuation": (
                    1.0 - norms["after_genus_residual"] / observed
                    if observed > 0 and "after_genus_residual" in norms
                    else float("nan")
                ),
            }
        )
        gate_row = gate_index.loc[context_value] if gate_index is not None and context_value in gate_index.index else None
        classifications.append(
            {
                "context": context_value,
                "classification": _classify_context(
                    gate_row,
                    part,
                    alpha=float(depth_config["multiplicity"]["alpha"]),
                ),
                "common_support_H1_q": (
                    float(gate_row["q_value"])
                    if gate_row is not None and pd.notna(gate_row.get("q_value"))
                    else float("nan")
                ),
            }
        )

    support = (
        residuals.groupby("outcome", as_index=False)
        .agg(
            n_species_common_support=("accepted_species", "nunique"),
            n_families=("family", "nunique"),
            n_genera=("genus", "nunique"),
            median_family_size=("family_n", "median"),
            median_genus_size=("genus_n", "median"),
        )
    )
    manifest = {
        "contract": depth_config["contract"],
        "evidence_scope": evidence_scope,
        "flora_scope": flora_scope,
        "n_species_outcome_rows": int(len(residuals)),
        "n_unique_islands": int(island_scores["island_id"].nunique()) if not island_scores.empty else 0,
        "common_support_gate_model": "current seven-response beta-binomial H1",
        "decomposition_model": "equal-island clustered linear representation-depth audit",
        "claim_boundary": (
            "Persistence below genus does not prove within-lineage evolution; attenuation at "
            "genus does not prove dispersal or colonization filtering. The audit localizes "
            "taxonomic representation depth only. Paired spatial-block bootstrap "
            "quantifies attenuation uncertainty but does not convert the audit into a "
            "causal assembly test."
        ),
    }
    return {
        "common_support_gate_slopes": gate_slopes,
        "common_support_gate": gate,
        "stage_slopes": stage_slopes,
        "stage_omnibus": stage_omnibus,
        "attenuation": pd.DataFrame(attenuation_rows),
        "attenuation_bootstrap": bootstrap_draws,
        "attenuation_bootstrap_summary": bootstrap_summary,
        "classification": pd.DataFrame(classifications),
        "support": support,
        "manifest": manifest,
    }


@app.command()
def main(
    state_audit_csv: Path = typer.Option(..., exists=True),
    taxonomy_csv: Path = typer.Option(..., exists=True),
    status_flora_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    probability_config_path: Path = typer.Option(..., exists=True),
    depth_config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    flora_scope: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
) -> None:
    probability_config = yaml.safe_load(probability_config_path.read_text(encoding="utf-8"))
    depth_config = load_config(depth_config_path)
    outputs = run_analysis(
        pd.read_csv(state_audit_csv),
        pd.read_csv(taxonomy_csv),
        pd.read_csv(status_flora_csv),
        pd.read_csv(covariates_csv),
        probability_config,
        depth_config,
        evidence_scope=evidence_scope,
        flora_scope=flora_scope,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, frame in outputs.items():
        if isinstance(frame, pd.DataFrame):
            frame.to_csv(output_dir / f"{name}.csv", index=False)
    (output_dir / "manifest.json").write_text(
        json.dumps(outputs["manifest"], indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(outputs["classification"].to_csv(index=False))


if __name__ == "__main__":
    app()
