"""Atomic H4 taxonomic-depth analysis for redesigned Chapter 1.

The module converts the six atomic trait audits to species-level 0/1 response scores,
builds the existing outcome-blind GIFT source-matched family/genus decomposition, and
then asks whether the H3 distance-by-area vector persists after family and genus
composition are removed. Equal island weight is used deliberately as an area-support
guardrail.
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

from island_v2.chapter1_all_data_probability import _bh, _chi_square_sf_integer_df
from island_v2.chapter1_taxonomic_depth_decomposition import (
    build_source_group_contract,
    build_taxonomic_decomposition,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _truthy(series: pd.Series) -> pd.Series:
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def _classify(value: object, positive: set[str], negative: set[str]) -> float:
    tokens = {x.strip() for x in str(value or "").split("|") if x.strip()}
    if not tokens:
        return float("nan")
    if tokens <= positive:
        return 1.0
    if tokens <= negative:
        return 0.0
    return float("nan")


def atomic_species_scores(state_audit: pd.DataFrame, probability_config: dict[str, Any]) -> pd.DataFrame:
    required = {"accepted_species", "trait_name", "resolved_for_primary", "canonical_signature"}
    missing = required - set(state_audit.columns)
    if missing:
        raise typer.BadParameter(f"state audit missing columns: {sorted(missing)}")
    audit = state_audit.loc[_truthy(state_audit["resolved_for_primary"])].copy()
    rows: list[pd.DataFrame] = []
    for outcome in probability_config["model_outcomes"]:
        spec = probability_config["broad_outcomes"][str(outcome)]
        part = audit.loc[
            audit["trait_name"].astype(str).eq(str(spec["trait_name"])),
            ["accepted_species", "canonical_signature"],
        ].drop_duplicates("accepted_species")
        positive = {str(x) for x in spec["positive_states"]}
        negative = {str(x) for x in spec["negative_states"]}
        part = part.copy()
        part["syndrome_concordance"] = [
            _classify(value, positive, negative) for value in part["canonical_signature"]
        ]
        part = part.dropna(subset=["syndrome_concordance"])
        part["syndrome"] = str(outcome)
        rows.append(part[["accepted_species", "syndrome", "syndrome_concordance"]])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def _z(values: pd.Series) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("constant predictor")
    return (x - mean) / sd


def _z_masked(frame: pd.DataFrame, mask: np.ndarray, column: str) -> np.ndarray:
    x = pd.to_numeric(frame.loc[mask, column], errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError(f"constant predictor: {column}")
    out = np.zeros(len(frame), dtype=float)
    out[mask] = (x - mean) / sd
    return out


def _fit_clustered_ols(y: np.ndarray, design: np.ndarray, clusters: np.ndarray):
    y = np.asarray(y, dtype=float)
    design = np.asarray(design, dtype=float)
    beta = np.linalg.pinv(design.T @ design) @ (design.T @ y)
    residual = y - design @ beta
    bread = np.linalg.pinv(design.T @ design)
    labels = np.asarray(clusters).astype(str)
    unique = np.unique(labels)
    meat = np.zeros((design.shape[1], design.shape[1]), dtype=float)
    for label in unique:
        mask = labels == label
        score = design[mask].T @ residual[mask]
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    n = len(y)
    k = design.shape[1]
    g = len(unique)
    if g > 1 and n > k:
        covariance *= (g / (g - 1.0)) * ((n - 1.0) / (n - k))
    return beta, covariance, g


def _joint(vector: np.ndarray, covariance: np.ndarray) -> tuple[float, int, float]:
    rank = int(np.linalg.matrix_rank(covariance))
    if rank <= 0:
        return float("nan"), 0, float("nan")
    stat = float(vector @ np.linalg.pinv(covariance) @ vector)
    return stat, rank, _chi_square_sf_integer_df(stat, rank)


def fit_stage_model(
    decomposition: pd.DataFrame,
    covariates: pd.DataFrame,
    probability_config: dict[str, Any],
    ladder_config: dict[str, Any],
    h4_config: dict[str, Any],
    *,
    source_mode: str,
    stratum: str,
    context_value: str,
    stage: str,
) -> dict[str, Any]:
    outcomes = [str(x) for x in h4_config["atomic_outcomes"]]
    context = str(ladder_config["context_column"])
    cluster = str(ladder_config["cluster_column"])
    distance = str(ladder_config["geography_column"])
    area = str(ladder_config["area_column"])
    controls = [str(x) for x in ladder_config["control_columns"]]
    threshold = int(h4_config["minimum_support"]["islands_per_context_outcome"])

    work = decomposition.loc[
        decomposition["source_mode"].eq(source_mode)
        & decomposition["stratum"].eq(stratum)
        & decomposition["syndrome"].isin(outcomes)
    ].copy()
    required_cov = ["island_id", context, cluster, distance, area, *controls]
    work = work.merge(
        covariates[required_cov].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    work = work.loc[work[context].astype(str).eq(context_value)].copy()
    for column in [stage, distance, area, *controls]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.dropna(subset=[stage, distance, area, *controls])
    support = work.groupby("syndrome")["island_id"].nunique()
    retained = [outcome for outcome in outcomes if int(support.get(outcome, 0)) >= threshold]
    if len(retained) < 2:
        return {
            "source_mode": source_mode,
            "stratum": stratum,
            "context": context_value,
            "stage": stage,
            "status": "not_testable",
            "n_retained_outcomes": len(retained),
            "retained_outcomes": "|".join(retained),
        }
    work = work.loc[work["syndrome"].isin(retained)].reset_index(drop=True)
    columns: list[np.ndarray] = []
    names: list[str] = []
    distance_indices: list[int] = []
    moderation_indices: list[int] = []
    for outcome in retained:
        mask = work["syndrome"].eq(outcome).to_numpy()
        indicator = mask.astype(float)
        columns.append(indicator)
        names.append(f"{outcome}:intercept")
        for predictor in controls:
            columns.append(_z_masked(work, mask, predictor))
            names.append(f"{outcome}:z_{predictor}")
        za = _z_masked(work, mask, area)
        zd = _z_masked(work, mask, distance)
        columns.extend([za, zd, za * zd])
        names.extend(
            [
                f"{outcome}:z_{area}",
                f"{outcome}:z_{distance}",
                f"{outcome}:z_{distance}:z_{area}",
            ]
        )
        distance_indices.append(len(names) - 2)
        moderation_indices.append(len(names) - 1)
    design = np.column_stack(columns)
    beta, covariance, n_clusters = _fit_clustered_ols(
        work[stage].to_numpy(float), design, work[cluster].to_numpy(str)
    )
    distance_vector = beta[distance_indices]
    moderation_vector = beta[moderation_indices]
    distance_cov = covariance[np.ix_(distance_indices, distance_indices)]
    moderation_cov = covariance[np.ix_(moderation_indices, moderation_indices)]
    d_stat, d_df, d_p = _joint(distance_vector, distance_cov)
    m_stat, m_df, m_p = _joint(moderation_vector, moderation_cov)
    return {
        "source_mode": source_mode,
        "stratum": stratum,
        "context": context_value,
        "stage": stage,
        "status": "fit",
        "n_retained_outcomes": len(retained),
        "retained_outcomes": "|".join(retained),
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(n_clusters),
        "distance_joint_wald": d_stat,
        "distance_joint_df": d_df,
        "distance_p_value": d_p,
        "moderation_joint_wald": m_stat,
        "moderation_joint_df": m_df,
        "moderation_p_value": m_p,
    }


def classify_depth(stage_results: pd.DataFrame, alpha: float) -> pd.DataFrame:
    work = stage_results.loc[stage_results["status"].eq("fit")].copy()
    if work.empty:
        return pd.DataFrame()
    work["moderation_q_across_source_modes"] = work.groupby(
        ["stratum", "context", "stage"], group_keys=False
    )["moderation_p_value"].transform(_bh)
    work["moderation_supported"] = work["moderation_q_across_source_modes"].le(alpha)
    rows: list[dict[str, Any]] = []
    for (stratum, context, source_mode), group in work.groupby(
        ["stratum", "context", "source_mode"], sort=True
    ):
        status = group.set_index("stage")["moderation_supported"].to_dict()
        observed = bool(status.get("observed_score", False))
        family = bool(status.get("after_family_residual", False))
        genus = bool(status.get("after_genus_residual", False))
        if not observed:
            classification = "observed_not_supported"
        elif not family:
            classification = "compatible_with_family_or_deeper_sorting"
        elif not genus:
            classification = "compatible_with_genus_level_assembly_beyond_family"
        else:
            classification = "residual_below_genus_retained"
        rows.append(
            {
                "stratum": stratum,
                "context": context,
                "source_mode": source_mode,
                "classification": classification,
            }
        )
    per_mode = pd.DataFrame(rows)
    aggregate: list[dict[str, Any]] = []
    for (stratum, context), group in per_mode.groupby(["stratum", "context"], sort=True):
        values = sorted(set(group["classification"]))
        aggregate.append(
            {
                "stratum": stratum,
                "context": context,
                "n_source_modes": int(group["source_mode"].nunique()),
                "source_mode_classifications": "|".join(
                    f"{r.source_mode}:{r.classification}" for r in group.itertuples()
                ),
                "robust_classification": values[0] if len(values) == 1 else "source_definition_sensitive",
            }
        )
    return per_mode.merge(pd.DataFrame(aggregate), on=["stratum", "context"], how="left")


def run_h4(
    state_audit: pd.DataFrame,
    taxonomy: pd.DataFrame,
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    assignments: pd.DataFrame,
    covariates: pd.DataFrame,
    probability_config: dict[str, Any],
    ladder_config: dict[str, Any],
    h4_config: dict[str, Any],
    source_config: dict[str, Any],
):
    scores = atomic_species_scores(state_audit, probability_config)
    taxonomy = taxonomy[["accepted_species", "genus", "family"]].copy().fillna("")
    source_modes = [str(x) for x in source_config["source_assignment"]["primary_modes"]]
    strata = [str(h4_config["primary_stratum"]), *[str(x) for x in h4_config["sensitivity_strata"]]]
    axes = [str(x) for x in h4_config["atomic_outcomes"]]
    positions, availability, _, source_audit = build_source_group_contract(
        scores,
        taxonomy,
        gift_flora,
        axes=axes,
        minimum_source_scored_species=int(h4_config["source_group_training"]["minimum_scored_species_per_group"]),
    )
    minimum = h4_config["minimum_support"]
    decomposition = build_taxonomic_decomposition(
        scores,
        taxonomy,
        status_flora,
        positions,
        availability,
        assignments,
        covariates,
        axes=axes,
        source_modes=source_modes,
        strata=strata,
        minimum_species=int(minimum["observed_species_per_island_response"]),
        minimum_families=int(minimum["represented_families"]),
        minimum_genera=int(minimum["represented_genera"]),
    )
    stage_rows = []
    for source_mode in source_modes:
        for stratum in strata:
            for context_value in [str(x) for x in h4_config["contexts"]]:
                for stage in [str(x) for x in h4_config["stages"]]:
                    stage_rows.append(
                        fit_stage_model(
                            decomposition,
                            covariates,
                            probability_config,
                            ladder_config,
                            h4_config,
                            source_mode=source_mode,
                            stratum=stratum,
                            context_value=context_value,
                            stage=stage,
                        )
                    )
    stage_results = pd.DataFrame(stage_rows)
    classification = classify_depth(stage_results, float(ladder_config["alpha"]))
    manifest = {
        "contract": h4_config["contract"],
        "n_atomic_species_scores": int(len(scores)),
        "n_decomposition_rows": int(len(decomposition)),
        "source_modes": source_modes,
        "source_audit": source_audit,
        "equal_island_weight": True,
    }
    return scores, decomposition, stage_results, classification, manifest


@app.command("run")
def run(
    state_audit_csv: Path = typer.Option(..., exists=True),
    taxonomy_csv: Path = typer.Option(..., exists=True),
    status_flora_csv: Path = typer.Option(..., exists=True),
    gift_flora_csv: Path = typer.Option(..., exists=True),
    source_assignments_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    probability_config_path: Path = typer.Option(..., exists=True),
    ladder_config_path: Path = typer.Option(..., exists=True),
    h4_config_path: Path = typer.Option(..., exists=True),
    source_config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
) -> None:
    probability_config = yaml.safe_load(probability_config_path.read_text(encoding="utf-8"))
    ladder_config = yaml.safe_load(ladder_config_path.read_text(encoding="utf-8"))
    h4_config = yaml.safe_load(h4_config_path.read_text(encoding="utf-8"))
    source_config = yaml.safe_load(source_config_path.read_text(encoding="utf-8"))
    scores, decomposition, stage_results, classification, manifest = run_h4(
        pd.read_csv(state_audit_csv),
        pd.read_csv(taxonomy_csv),
        pd.read_csv(status_flora_csv),
        pd.read_csv(gift_flora_csv),
        pd.read_csv(source_assignments_csv),
        pd.read_csv(covariates_csv),
        probability_config,
        ladder_config,
        h4_config,
        source_config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    for frame in (scores, decomposition, stage_results, classification):
        if not frame.empty and "evidence_scope" not in frame.columns:
            frame.insert(0, "evidence_scope", evidence_scope)
    scores.to_csv(output_dir / "atomic_species_scores.csv.gz", index=False)
    decomposition.to_csv(output_dir / "atomic_taxonomic_decomposition.csv.gz", index=False)
    stage_results.to_csv(output_dir / "h4_atomic_stage_results.csv", index=False)
    classification.to_csv(output_dir / "h4_atomic_depth_classification.csv", index=False)
    manifest["evidence_scope"] = evidence_scope
    (output_dir / "h4_atomic_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    app()
