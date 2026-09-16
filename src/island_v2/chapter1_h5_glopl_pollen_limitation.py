"""Frozen GloPL pollen-limitation H5 gate for Chapter 1.

The geographic admission gate was completed without reading pollen-limitation outcomes.
This module is the separately frozen next stage. It validates the pinned GloPL bytes
before reading the predeclared effect columns, reuses exact preflight row-to-island
assignments, aggregates repeated rows to publication x island cells, and tests whether
experimental pollen limitation covaries with mainland isolation differently in northern
mid-latitude versus tropical islands.

A supported result is an experimental pollen-limitation association, not evidence of
pollinator abundance decline, visitor identity, historical extinction, floral selection,
or mediation of the plant-trait H2 response.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT = "chapter1_h5_glopl_pollen_limitation_v1"
CONTEXTS = ("northern_midlatitude", "tropical")
RAW_COVARIATES = (
    "log_distance_to_continent_km",
    "log_island_area_km2",
    "climate_pc1",
    "climate_pc2",
    "climate_pc3",
    "climate_pc4",
)
Z_COVARIATES = tuple(f"z_{x}" for x in RAW_COVARIATES)
EFFECT_COLUMNS = (
    "PL_Effect_Size",
    "PL_Effect_Size_Type1",
    "PL_Effect_Size_Type2",
    "Constant_added",
    "Level_of_Supplementation",
)


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected GloPL pollen-limitation contract")
    return config


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_source_sha256(path: Path, config: dict[str, Any]) -> str:
    observed = _sha256(path)
    expected = str(config["source"]["sha256"])
    if observed != expected:
        raise ValueError(f"source SHA-256 mismatch: expected {expected}, observed {observed}")
    return observed


def read_effect_columns(path: Path, config: dict[str, Any]) -> pd.DataFrame:
    allowed = [str(x) for x in config["columns_allowed_after_freeze"]]
    if allowed != list(EFFECT_COLUMNS):
        raise typer.BadParameter("unexpected frozen GloPL effect-column contract")
    frame = pd.read_csv(path, usecols=allowed, dtype=str, encoding="latin-1").fillna("")
    frame = frame[allowed].reset_index(drop=True)
    frame.insert(0, "row_id", np.arange(len(frame), dtype=int))
    return frame


def _truthy(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().str.casefold().isin(
        {"true", "1", "yes", "y", "t"}
    )


def _unique_island_covariates(covariates: pd.DataFrame) -> pd.DataFrame:
    required = {"island_id", "spatial_block", *RAW_COVARIATES}
    if missing := required - set(covariates.columns):
        raise typer.BadParameter(f"GloPL island covariates missing columns: {sorted(missing)}")
    cols = ["island_id", "spatial_block", *RAW_COVARIATES]
    work = covariates[cols].copy()
    work["island_id"] = work["island_id"].astype(str)
    multiplicity = work.groupby("island_id", dropna=False)[cols[1:]].nunique(dropna=False)
    conflicts = multiplicity.gt(1).any(axis=1)
    if bool(conflicts.any()):
        examples = [str(x) for x in conflicts.index[conflicts][:5]]
        raise typer.BadParameter(f"conflicting GloPL island covariates: {examples}")
    return work.drop_duplicates("island_id")


def _standardize_in_place(work: pd.DataFrame) -> pd.DataFrame:
    out = work.copy()
    for raw, z_name in zip(RAW_COVARIATES, Z_COVARIATES, strict=True):
        values = pd.to_numeric(out[raw], errors="coerce").to_numpy(float)
        mean = float(np.mean(values))
        sd = float(np.std(values, ddof=0))
        if not np.isfinite(sd) or sd <= 0:
            raise typer.BadParameter(f"constant or invalid GloPL predictor: {raw}")
        out[z_name] = (values - mean) / sd
    return out


def prepare_analysis_rows(
    effects: pd.DataFrame,
    matched: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    required_effect = {"row_id", *EFFECT_COLUMNS}
    if missing := required_effect - set(effects.columns):
        raise typer.BadParameter(f"GloPL effect table missing columns: {sorted(missing)}")
    required_match = {
        "row_id",
        "island_id",
        "study_key",
        "analysis_regime",
        "boundary_only_match",
        "multi_polygon_match",
    }
    if missing := required_match - set(matched.columns):
        raise typer.BadParameter(f"GloPL preflight match missing columns: {sorted(missing)}")
    if effects["row_id"].duplicated().any() or matched["row_id"].duplicated().any():
        raise typer.BadParameter("GloPL row_id must be unique in effect and preflight tables")

    work = effects.merge(
        matched[list(required_match)], on="row_id", how="left", validate="one_to_one"
    )
    work["PL_Effect_Size"] = pd.to_numeric(work["PL_Effect_Size"], errors="coerce")
    work["study_key"] = work["study_key"].fillna("").astype(str).str.strip()
    work["analysis_regime"] = work["analysis_regime"].fillna("").astype(str)
    boundary = work["boundary_only_match"].fillna(False).astype(bool)
    multi = work["multi_polygon_match"].fillna(False).astype(bool)
    finite = np.isfinite(work["PL_Effect_Size"].to_numpy(float))
    eligible = (
        work["island_id"].notna()
        & ~boundary
        & ~multi
        & work["analysis_regime"].isin(CONTEXTS)
        & work["study_key"].ne("")
        & finite
    )
    work = work.loc[eligible].copy()
    work["island_id"] = work["island_id"].astype(str)

    island_cov = _unique_island_covariates(covariates)
    work = work.merge(island_cov, on="island_id", how="left", validate="many_to_one")
    for column in RAW_COVARIATES:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work["spatial_block"] = work["spatial_block"].fillna("").astype(str)
    work = work.dropna(subset=list(RAW_COVARIATES))
    work = work.loc[work["spatial_block"].ne("")].copy()
    work = _standardize_in_place(work)
    return work.sort_values("row_id").reset_index(drop=True)


def aggregate_publication_island(rows: pd.DataFrame) -> pd.DataFrame:
    required = {
        "study_key",
        "island_id",
        "analysis_regime",
        "spatial_block",
        "PL_Effect_Size",
        *Z_COVARIATES,
    }
    if missing := required - set(rows.columns):
        raise typer.BadParameter(f"publication-island aggregation missing columns: {sorted(missing)}")
    if rows.empty:
        return pd.DataFrame()
    fixed = ["analysis_regime", "spatial_block", *Z_COVARIATES]
    mult = rows.groupby(["study_key", "island_id"], dropna=False)[fixed].nunique(dropna=False)
    if bool(mult.gt(1).any(axis=1).any()):
        raise typer.BadParameter("publication x island cell has conflicting island covariates")

    aggregations: dict[str, Any] = {
        "analysis_regime": ("analysis_regime", "first"),
        "spatial_block": ("spatial_block", "first"),
        "PL_Effect_Size": ("PL_Effect_Size", "mean"),
        "n_effect_rows": ("PL_Effect_Size", "size"),
    }
    for z_name in Z_COVARIATES:
        aggregations[z_name] = (z_name, "first")
    out = (
        rows.groupby(["study_key", "island_id"], as_index=False)
        .agg(**aggregations)
        .reset_index(drop=True)
    )
    n_cells = out.groupby("study_key")["island_id"].transform("size").astype(float)
    out["analysis_weight"] = 1.0 / n_cells
    return out


def build_primary_design(frame: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    required = {"analysis_regime", *Z_COVARIATES}
    if missing := required - set(frame.columns):
        raise typer.BadParameter(f"GloPL model frame missing columns: {sorted(missing)}")
    tropical = frame["analysis_regime"].astype(str).eq("tropical").to_numpy(float)
    distance = pd.to_numeric(frame["z_log_distance_to_continent_km"], errors="coerce").to_numpy(float)
    names = [
        "intercept",
        "context_tropical",
        "z_log_distance_to_continent_km",
        "z_log_distance_to_continent_km:context_tropical",
        "z_log_island_area_km2",
        "z_climate_pc1",
        "z_climate_pc2",
        "z_climate_pc3",
        "z_climate_pc4",
    ]
    columns = [
        np.ones(len(frame), dtype=float),
        tropical,
        distance,
        distance * tropical,
        *[
            pd.to_numeric(frame[name], errors="coerce").to_numpy(float)
            for name in Z_COVARIATES[1:]
        ],
    ]
    return np.column_stack(columns), names


def _normal_two_sided_p(z_value: float) -> float:
    if not math.isfinite(z_value):
        return float("nan")
    return math.erfc(abs(float(z_value)) / math.sqrt(2.0))


def _normal_upper_p(z_value: float) -> float:
    if not math.isfinite(z_value):
        return float("nan")
    return 0.5 * math.erfc(float(z_value) / math.sqrt(2.0))


def _normal_lower_p(z_value: float) -> float:
    if not math.isfinite(z_value):
        return float("nan")
    return 0.5 * math.erfc(-float(z_value) / math.sqrt(2.0))


def _fit_clustered_wls(
    frame: pd.DataFrame,
    *,
    weights: np.ndarray,
    config: dict[str, Any],
) -> dict[str, Any]:
    if frame.empty:
        return {"evaluable": False, "reason": "no_rows"}
    counts = frame.groupby("analysis_regime")["island_id"].size().to_dict()
    min_cells = int(config["model"]["minimum_publication_island_cells_per_context"])
    if any(int(counts.get(context, 0)) < min_cells for context in CONTEXTS):
        return {
            "evaluable": False,
            "reason": "below_minimum_publication_island_cells_per_context",
            "n_cells_by_context": {context: int(counts.get(context, 0)) for context in CONTEXTS},
        }
    clusters = frame["spatial_block"].astype(str).to_numpy()
    n_clusters = int(pd.Series(clusters).nunique())
    min_clusters = int(config["model"]["minimum_clusters"])
    if n_clusters < min_clusters:
        return {
            "evaluable": False,
            "reason": "below_minimum_clusters",
            "n_clusters": n_clusters,
        }
    design, names = build_primary_design(frame)
    y = pd.to_numeric(frame["PL_Effect_Size"], errors="coerce").to_numpy(float)
    weights = np.asarray(weights, dtype=float)
    if (
        len(y) != len(weights)
        or not np.isfinite(y).all()
        or not np.isfinite(design).all()
        or not np.isfinite(weights).all()
        or np.any(weights <= 0)
    ):
        return {"evaluable": False, "reason": "nonfinite_model_input"}
    n, p = design.shape
    if n <= p + 2:
        return {"evaluable": False, "reason": "insufficient_rows_for_design"}
    xtwx = design.T @ (weights[:, None] * design)
    beta = np.linalg.pinv(xtwx) @ (design.T @ (weights * y))
    residual = y - design @ beta
    bread = np.linalg.pinv(xtwx)
    meat = np.zeros((p, p), dtype=float)
    for cluster in np.unique(clusters):
        mask = clusters == cluster
        score = design[mask].T @ (weights[mask] * residual[mask])
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    if n_clusters > 1 and n > p:
        covariance *= (n_clusters / (n_clusters - 1.0)) * ((n - 1.0) / (n - p))
    if not np.isfinite(covariance).all():
        return {"evaluable": False, "reason": "nonfinite_covariance"}
    se = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
    idx = {name: i for i, name in enumerate(names)}
    d = idx["z_log_distance_to_continent_km"]
    interaction_idx = idx["z_log_distance_to_continent_km:context_tropical"]
    north = float(beta[d])
    north_se = float(se[d])
    interaction = float(beta[interaction_idx])
    interaction_se = float(se[interaction_idx])
    north_z = north / north_se if north_se > 0 else float("nan")
    int_z = interaction / interaction_se if interaction_se > 0 else float("nan")
    tropical = north + interaction
    tropical_var = float(
        covariance[d, d]
        + covariance[interaction_idx, interaction_idx]
        + 2.0 * covariance[d, interaction_idx]
    )
    tropical_se = math.sqrt(max(tropical_var, 0.0))
    tropical_z = tropical / tropical_se if tropical_se > 0 else float("nan")
    return {
        "evaluable": True,
        "n_publication_island_cells": int(len(frame)),
        "n_cells_by_context": {context: int(counts.get(context, 0)) for context in CONTEXTS},
        "n_unique_islands": int(frame["island_id"].nunique()),
        "n_publications": int(frame["study_key"].nunique()),
        "n_clusters": n_clusters,
        "northern_distance_slope": north,
        "northern_distance_slope_se": north_se,
        "northern_distance_slope_two_sided_p": _normal_two_sided_p(north_z),
        "northern_one_sided_positive_p": _normal_upper_p(north_z),
        "interaction_tropical_minus_northern": interaction,
        "interaction_se": interaction_se,
        "interaction_two_sided_p": _normal_two_sided_p(int_z),
        "interaction_one_sided_negative_p": _normal_lower_p(int_z),
        "tropical_distance_slope": tropical,
        "tropical_distance_slope_se": tropical_se,
        "tropical_distance_slope_two_sided_p": _normal_two_sided_p(tropical_z),
        "coefficient_names": names,
        "coefficients": {name: float(beta[i]) for i, name in enumerate(names)},
    }


def fit_publication_weighted(rows: pd.DataFrame, config: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, Any]]:
    cells = aggregate_publication_island(rows)
    if cells.empty:
        return cells, {"evaluable": False, "reason": "no_publication_island_cells"}
    result = _fit_clustered_wls(
        cells,
        weights=cells["analysis_weight"].to_numpy(float),
        config=config,
    )
    return cells, result


def _fit_equal_island(rows: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    cells = aggregate_publication_island(rows)
    if cells.empty:
        return {"evaluable": False, "reason": "no_publication_island_cells"}
    fixed = ["analysis_regime", "spatial_block", *Z_COVARIATES]
    mult = cells.groupby("island_id")[fixed].nunique(dropna=False)
    if bool(mult.gt(1).any(axis=1).any()):
        return {"evaluable": False, "reason": "conflicting_island_covariates"}
    aggregations: dict[str, Any] = {
        "analysis_regime": ("analysis_regime", "first"),
        "spatial_block": ("spatial_block", "first"),
        "PL_Effect_Size": ("PL_Effect_Size", "mean"),
        "study_key": ("study_key", lambda x: "|".join(sorted(set(map(str, x))))),
    }
    for z_name in Z_COVARIATES:
        aggregations[z_name] = (z_name, "first")
    island = cells.groupby("island_id", as_index=False).agg(**aggregations)
    # The frozen equal-island guardrail requested HC3. With one row per island,
    # island-level spatial-block clustering is not needed; compute HC3 directly.
    counts = island.groupby("analysis_regime")["island_id"].size().to_dict()
    min_cells = int(config["model"]["minimum_publication_island_cells_per_context"])
    if any(int(counts.get(context, 0)) < min_cells for context in CONTEXTS):
        return {
            "evaluable": False,
            "reason": "below_minimum_islands_per_context_equal_island",
            "n_cells_by_context": {context: int(counts.get(context, 0)) for context in CONTEXTS},
        }
    X, names = build_primary_design(island)
    y = island["PL_Effect_Size"].to_numpy(float)
    n, p = X.shape
    if n <= p + 2 or np.linalg.matrix_rank(X) < p:
        return {"evaluable": False, "reason": "equal_island_design_not_full_rank"}
    xtx_inv = np.linalg.inv(X.T @ X)
    beta = xtx_inv @ X.T @ y
    residual = y - X @ beta
    leverage = np.einsum("ij,jk,ik->i", X, xtx_inv, X)
    denom = np.clip(1.0 - leverage, 1e-8, None)
    scaled = residual / denom
    meat = X.T @ ((scaled**2)[:, None] * X)
    cov = xtx_inv @ meat @ xtx_inv
    se = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    idx = {name: i for i, name in enumerate(names)}
    d = idx["z_log_distance_to_continent_km"]
    j = idx["z_log_distance_to_continent_km:context_tropical"]
    north = float(beta[d])
    interaction = float(beta[j])
    north_se = float(se[d])
    int_se = float(se[j])
    north_z = north / north_se if north_se > 0 else float("nan")
    int_z = interaction / int_se if int_se > 0 else float("nan")
    return {
        "evaluable": True,
        "n_islands": int(len(island)),
        "n_islands_by_context": {context: int(counts.get(context, 0)) for context in CONTEXTS},
        "northern_distance_slope": north,
        "northern_one_sided_positive_p": _normal_upper_p(north_z),
        "interaction_tropical_minus_northern": interaction,
        "interaction_one_sided_negative_p": _normal_lower_p(int_z),
        "interaction_two_sided_p": _normal_two_sided_p(int_z),
    }


def classify_h5_result(
    primary: dict[str, Any], sensitivities: dict[str, dict[str, Any]], config: dict[str, Any]
) -> dict[str, Any]:
    alpha = float(config["alpha"])
    north_supported = bool(
        primary.get("evaluable")
        and float(primary.get("northern_distance_slope", float("nan"))) > 0
        and float(primary.get("northern_one_sided_positive_p", 1.0)) <= alpha
    )
    context_supported = bool(
        north_supported
        and float(primary.get("interaction_tropical_minus_northern", float("nan"))) < 0
        and float(primary.get("interaction_one_sided_negative_p", 1.0)) <= alpha
    )
    required_sens = ["supplemental_only", "no_zero_constant"]
    sensitivity_direction_consistent = True
    for name in required_sens:
        result = sensitivities.get(name, {})
        sensitivity_direction_consistent = sensitivity_direction_consistent and bool(
            result.get("evaluable")
            and float(result.get("northern_distance_slope", float("nan"))) > 0
            and float(result.get("interaction_tropical_minus_northern", float("nan"))) < 0
        )
    full = bool(context_supported and sensitivity_direction_consistent)
    if full:
        classification = "context_specific_pollen_limitation_gradient_supported"
    elif context_supported:
        classification = "primary_context_specific_gradient_sensitivity_fragile_not_promoted"
    elif north_supported:
        classification = "north_pollen_limitation_gradient_supported_context_specificity_not_established"
    elif primary.get("evaluable"):
        classification = "pollen_limitation_gradient_not_supported"
    else:
        classification = "pollen_limitation_gradient_not_evaluable"
    return {
        "classification": classification,
        "north_service_gradient_supported": north_supported,
        "context_specific_service_gradient_supported": context_supported,
        "sensitivity_direction_consistent": sensitivity_direction_consistent,
        "full_mechanism_promotion": full,
        "global_pollinator_causal_mechanism_identified": False,
    }


def _sensitivity_rows(rows: pd.DataFrame, name: str) -> pd.DataFrame:
    if name == "supplemental_only":
        return rows.loc[rows["PL_Effect_Size_Type2"].fillna("").astype(str).str.strip().eq("Sup")].copy()
    if name == "no_zero_constant":
        return rows.loc[~_truthy(rows["Constant_added"])].copy()
    raise ValueError(f"unknown sensitivity: {name}")


@app.command("run")
def run(
    glopl_csv: Path = typer.Option(..., exists=True),
    preflight_match_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    source_sha = validate_source_sha256(glopl_csv, config)
    # Outcome columns are read only after the pinned source bytes are verified.
    effects = read_effect_columns(glopl_csv, config)
    matched = pd.read_csv(preflight_match_csv, dtype=str).fillna("")
    for column in ("boundary_only_match", "multi_polygon_match"):
        if column in matched.columns:
            matched[column] = matched[column].astype(str).str.casefold().isin({"true", "1", "yes"})
    covariates = pd.read_csv(covariates_csv, dtype=str).fillna("")
    rows = prepare_analysis_rows(effects, matched, covariates, config)

    primary_cells, primary = fit_publication_weighted(rows, config)
    sensitivities: dict[str, dict[str, Any]] = {}
    sensitivity_cells: dict[str, pd.DataFrame] = {}
    for name in ("supplemental_only", "no_zero_constant"):
        cells, result = fit_publication_weighted(_sensitivity_rows(rows, name), config)
        sensitivity_cells[name] = cells
        sensitivities[name] = result
    sensitivities["equal_island"] = _fit_equal_island(rows, config)
    decision = classify_h5_result(primary, sensitivities, config)

    result = {
        "contract": CONTRACT,
        "status": "completed",
        "source_sha256": source_sha,
        "n_exact_island_finite_effect_rows": int(len(rows)),
        "primary": primary,
        "sensitivities": sensitivities,
        "decision": decision,
        "claim_ceiling": config["claim_ceiling"],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    primary_cells.to_csv(output_dir / "primary_publication_island_cells.csv", index=False)
    for name, cells in sensitivity_cells.items():
        cells.to_csv(output_dir / f"{name}_publication_island_cells.csv", index=False)
    (output_dir / "RESULT.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    lines = [
        "# H5 GloPL experimental pollen-limitation gate",
        "",
        f"- exact-island finite-effect rows: {len(rows)}",
        f"- primary evaluable: {primary.get('evaluable')}",
    ]
    if primary.get("evaluable"):
        lines.extend(
            [
                f"- primary cells / islands / publications / blocks: {primary['n_publication_island_cells']} / {primary['n_unique_islands']} / {primary['n_publications']} / {primary['n_clusters']}",
                f"- North distance slope: {primary['northern_distance_slope']:.8f}; one-sided positive p={primary['northern_one_sided_positive_p']:.8g}",
                f"- Tropical-minus-North interaction: {primary['interaction_tropical_minus_northern']:.8f}; one-sided negative p={primary['interaction_one_sided_negative_p']:.8g}; two-sided p={primary['interaction_two_sided_p']:.8g}",
                f"- Tropical distance slope: {primary['tropical_distance_slope']:.8f}; two-sided p={primary['tropical_distance_slope_two_sided_p']:.8g}",
            ]
        )
    for name in ("supplemental_only", "no_zero_constant", "equal_island"):
        item = sensitivities[name]
        if item.get("evaluable"):
            lines.append(
                f"- {name}: North slope={item['northern_distance_slope']:.8f}; interaction={item['interaction_tropical_minus_northern']:.8f}"
            )
        else:
            lines.append(f"- {name}: not evaluable ({item.get('reason', 'unknown')})")
    lines.extend(
        [
            "",
            f"**Classification:** {decision['classification']}",
            "",
            "This analysis measures experimental pollen limitation, not pollinator abundance decline, visitor identity, historical extinction, or mediation of the floral-trait H2 response.",
        ]
    )
    (output_dir / "RESULT_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.echo(json.dumps(result, indent=2))


if __name__ == "__main__":
    app()
