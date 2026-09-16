"""Frozen GloPL pollen-limitation H5 gate for Chapter 1.

The geographic admission gate was completed without reading pollen-limitation outcomes.
This stage validates the pinned GloPL bytes before reading the frozen effect columns,
reuses exact preflight row-to-island assignments, collapses repeated observations to
publication x island cells, and tests isolation-associated experimental pollen limitation.
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


def _normalise_row_id(frame: pd.DataFrame, *, label: str) -> pd.DataFrame:
    out = frame.copy()
    numeric = pd.to_numeric(out["row_id"], errors="coerce")
    if numeric.isna().any() or not np.equal(numeric.to_numpy(float), np.floor(numeric.to_numpy(float))).all():
        raise typer.BadParameter(f"{label} row_id contains non-integer values")
    out["row_id"] = numeric.astype(int)
    if out["row_id"].duplicated().any():
        raise typer.BadParameter(f"{label} row_id must be unique")
    return out


def _unique_island_covariates(covariates: pd.DataFrame) -> pd.DataFrame:
    cols = ["island_id", "spatial_block", *RAW_COVARIATES]
    if missing := set(cols) - set(covariates.columns):
        raise typer.BadParameter(f"GloPL island covariates missing columns: {sorted(missing)}")
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
    required_match = {
        "row_id",
        "island_id",
        "study_key",
        "analysis_regime",
        "boundary_only_match",
        "multi_polygon_match",
    }
    if missing := required_effect - set(effects.columns):
        raise typer.BadParameter(f"GloPL effect table missing columns: {sorted(missing)}")
    if missing := required_match - set(matched.columns):
        raise typer.BadParameter(f"GloPL preflight match missing columns: {sorted(missing)}")

    effects = _normalise_row_id(effects, label="effect")
    matched = _normalise_row_id(matched, label="preflight")
    work = effects.merge(
        matched[list(required_match)], on="row_id", how="left", validate="one_to_one"
    )
    work["PL_Effect_Size"] = pd.to_numeric(work["PL_Effect_Size"], errors="coerce")
    work["study_key"] = work["study_key"].fillna("").astype(str).str.strip()
    work["analysis_regime"] = work["analysis_regime"].fillna("").astype(str)
    island_text = work["island_id"].fillna("").astype(str).str.strip()
    finite = np.isfinite(work["PL_Effect_Size"].to_numpy(float))
    eligible = (
        island_text.ne("")
        & ~work["boundary_only_match"].fillna(False).astype(bool)
        & ~work["multi_polygon_match"].fillna(False).astype(bool)
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
    multiplicity = rows.groupby(["study_key", "island_id"], dropna=False)[fixed].nunique(dropna=False)
    if bool(multiplicity.gt(1).any(axis=1).any()):
        raise typer.BadParameter("publication x island cell has conflicting island covariates")
    aggregations: dict[str, Any] = {
        "analysis_regime": ("analysis_regime", "first"),
        "spatial_block": ("spatial_block", "first"),
        "PL_Effect_Size": ("PL_Effect_Size", "mean"),
        "n_effect_rows": ("PL_Effect_Size", "size"),
    }
    for name in Z_COVARIATES:
        aggregations[name] = (name, "first")
    out = rows.groupby(["study_key", "island_id"], as_index=False).agg(**aggregations)
    n_cells = out.groupby("study_key")["island_id"].transform("size").astype(float)
    out["analysis_weight"] = 1.0 / n_cells
    return out.reset_index(drop=True)


def build_primary_design(frame: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    if missing := {"analysis_regime", *Z_COVARIATES} - set(frame.columns):
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
        np.ones(len(frame)),
        tropical,
        distance,
        distance * tropical,
        *[pd.to_numeric(frame[x], errors="coerce").to_numpy(float) for x in Z_COVARIATES[1:]],
    ]
    return np.column_stack(columns), names


def _p2(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0)) if math.isfinite(z) else float("nan")


def _upper(z: float) -> float:
    return 0.5 * math.erfc(z / math.sqrt(2.0)) if math.isfinite(z) else float("nan")


def _lower(z: float) -> float:
    return 0.5 * math.erfc(-z / math.sqrt(2.0)) if math.isfinite(z) else float("nan")


def _fit_clustered_wls(frame: pd.DataFrame, weights: np.ndarray, config: dict[str, Any]) -> dict[str, Any]:
    if frame.empty:
        return {"evaluable": False, "reason": "no_rows"}
    counts = frame.groupby("analysis_regime")["island_id"].size().to_dict()
    minimum = int(config["model"]["minimum_publication_island_cells_per_context"])
    if any(int(counts.get(c, 0)) < minimum for c in CONTEXTS):
        return {"evaluable": False, "reason": "below_minimum_publication_island_cells_per_context", "n_cells_by_context": {c: int(counts.get(c, 0)) for c in CONTEXTS}}
    clusters = frame["spatial_block"].astype(str).to_numpy()
    n_clusters = int(pd.Series(clusters).nunique())
    if n_clusters < int(config["model"]["minimum_clusters"]):
        return {"evaluable": False, "reason": "below_minimum_clusters", "n_clusters": n_clusters}
    X, names = build_primary_design(frame)
    y = frame["PL_Effect_Size"].to_numpy(float)
    w = np.asarray(weights, dtype=float)
    n, p = X.shape
    if not np.isfinite(X).all() or not np.isfinite(y).all() or not np.isfinite(w).all() or np.any(w <= 0):
        return {"evaluable": False, "reason": "nonfinite_model_input"}
    if n <= p + 2 or np.linalg.matrix_rank(X) < p:
        return {"evaluable": False, "reason": "design_not_full_rank"}
    bread = np.linalg.inv(X.T @ (w[:, None] * X))
    beta = bread @ (X.T @ (w * y))
    residual = y - X @ beta
    meat = np.zeros((p, p))
    for cluster in np.unique(clusters):
        mask = clusters == cluster
        score = X[mask].T @ (w[mask] * residual[mask])
        meat += np.outer(score, score)
    cov = bread @ meat @ bread
    cov *= (n_clusters / (n_clusters - 1.0)) * ((n - 1.0) / (n - p))
    se = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    idx = {name: i for i, name in enumerate(names)}
    d = idx["z_log_distance_to_continent_km"]
    j = idx["z_log_distance_to_continent_km:context_tropical"]
    north, interaction = float(beta[d]), float(beta[j])
    north_se, interaction_se = float(se[d]), float(se[j])
    nz = north / north_se if north_se > 0 else float("nan")
    iz = interaction / interaction_se if interaction_se > 0 else float("nan")
    tropical = north + interaction
    tropical_var = float(cov[d, d] + cov[j, j] + 2.0 * cov[d, j])
    tropical_se = math.sqrt(max(tropical_var, 0.0))
    tz = tropical / tropical_se if tropical_se > 0 else float("nan")
    return {
        "evaluable": True,
        "n_publication_island_cells": int(len(frame)),
        "n_cells_by_context": {c: int(counts.get(c, 0)) for c in CONTEXTS},
        "n_unique_islands": int(frame["island_id"].nunique()),
        "n_publications": int(frame["study_key"].nunique()),
        "n_clusters": n_clusters,
        "northern_distance_slope": north,
        "northern_distance_slope_se": north_se,
        "northern_distance_slope_two_sided_p": _p2(nz),
        "northern_one_sided_positive_p": _upper(nz),
        "interaction_tropical_minus_northern": interaction,
        "interaction_se": interaction_se,
        "interaction_two_sided_p": _p2(iz),
        "interaction_one_sided_negative_p": _lower(iz),
        "tropical_distance_slope": tropical,
        "tropical_distance_slope_se": tropical_se,
        "tropical_distance_slope_two_sided_p": _p2(tz),
        "coefficients": {name: float(beta[i]) for i, name in enumerate(names)},
    }


def fit_publication_weighted(rows: pd.DataFrame, config: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, Any]]:
    cells = aggregate_publication_island(rows)
    if cells.empty:
        return cells, {"evaluable": False, "reason": "no_publication_island_cells"}
    return cells, _fit_clustered_wls(cells, cells["analysis_weight"].to_numpy(float), config)


def _fit_equal_island(rows: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    cells = aggregate_publication_island(rows)
    if cells.empty:
        return {"evaluable": False, "reason": "no_publication_island_cells"}
    fixed = ["analysis_regime", "spatial_block", *Z_COVARIATES]
    if bool(cells.groupby("island_id")[fixed].nunique(dropna=False).gt(1).any(axis=1).any()):
        return {"evaluable": False, "reason": "conflicting_island_covariates"}
    aggs: dict[str, Any] = {
        "analysis_regime": ("analysis_regime", "first"),
        "spatial_block": ("spatial_block", "first"),
        "PL_Effect_Size": ("PL_Effect_Size", "mean"),
        "study_key": ("study_key", lambda x: "|".join(sorted(set(map(str, x))))),
    }
    for name in Z_COVARIATES:
        aggs[name] = (name, "first")
    island = cells.groupby("island_id", as_index=False).agg(**aggs)
    counts = island.groupby("analysis_regime")["island_id"].size().to_dict()
    minimum = int(config["model"]["minimum_publication_island_cells_per_context"])
    if any(int(counts.get(c, 0)) < minimum for c in CONTEXTS):
        return {"evaluable": False, "reason": "below_minimum_islands_per_context_equal_island", "n_islands_by_context": {c: int(counts.get(c, 0)) for c in CONTEXTS}}
    X, names = build_primary_design(island)
    y = island["PL_Effect_Size"].to_numpy(float)
    n, p = X.shape
    if n <= p + 2 or np.linalg.matrix_rank(X) < p:
        return {"evaluable": False, "reason": "equal_island_design_not_full_rank"}
    inv = np.linalg.inv(X.T @ X)
    beta = inv @ X.T @ y
    residual = y - X @ beta
    leverage = np.einsum("ij,jk,ik->i", X, inv, X)
    scaled = residual / np.clip(1.0 - leverage, 1e-8, None)
    cov = inv @ (X.T @ ((scaled**2)[:, None] * X)) @ inv
    se = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    idx = {name: i for i, name in enumerate(names)}
    d, j = idx["z_log_distance_to_continent_km"], idx["z_log_distance_to_continent_km:context_tropical"]
    north, interaction = float(beta[d]), float(beta[j])
    nz = north / float(se[d]) if se[d] > 0 else float("nan")
    iz = interaction / float(se[j]) if se[j] > 0 else float("nan")
    return {
        "evaluable": True,
        "n_islands": int(len(island)),
        "n_islands_by_context": {c: int(counts.get(c, 0)) for c in CONTEXTS},
        "northern_distance_slope": north,
        "northern_one_sided_positive_p": _upper(nz),
        "interaction_tropical_minus_northern": interaction,
        "interaction_one_sided_negative_p": _lower(iz),
        "interaction_two_sided_p": _p2(iz),
    }


def classify_h5_result(primary: dict[str, Any], sensitivities: dict[str, dict[str, Any]], config: dict[str, Any]) -> dict[str, Any]:
    alpha = float(config["alpha"])
    north = bool(primary.get("evaluable") and float(primary.get("northern_distance_slope", float("nan"))) > 0 and float(primary.get("northern_one_sided_positive_p", 1.0)) <= alpha)
    context = bool(north and float(primary.get("interaction_tropical_minus_northern", float("nan"))) < 0 and float(primary.get("interaction_one_sided_negative_p", 1.0)) <= alpha)
    consistent = all(
        bool(
            sensitivities.get(name, {}).get("evaluable")
            and float(sensitivities[name].get("northern_distance_slope", float("nan"))) > 0
            and float(sensitivities[name].get("interaction_tropical_minus_northern", float("nan"))) < 0
        )
        for name in ("supplemental_only", "no_zero_constant")
    )
    full = bool(context and consistent)
    if full:
        classification = "context_specific_pollen_limitation_gradient_supported"
    elif context:
        classification = "primary_context_specific_gradient_sensitivity_fragile_not_promoted"
    elif north:
        classification = "north_pollen_limitation_gradient_supported_context_specificity_not_established"
    elif primary.get("evaluable"):
        classification = "pollen_limitation_gradient_not_supported"
    else:
        classification = "pollen_limitation_gradient_not_evaluable"
    return {
        "classification": classification,
        "north_service_gradient_supported": north,
        "context_specific_service_gradient_supported": context,
        "sensitivity_direction_consistent": consistent,
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
    effects = read_effect_columns(glopl_csv, config)
    matched = pd.read_csv(preflight_match_csv, dtype=str).fillna("")
    for column in ("boundary_only_match", "multi_polygon_match"):
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
        lines += [
            f"- primary cells / islands / publications / blocks: {primary['n_publication_island_cells']} / {primary['n_unique_islands']} / {primary['n_publications']} / {primary['n_clusters']}",
            f"- North distance slope: {primary['northern_distance_slope']:.8f}; one-sided positive p={primary['northern_one_sided_positive_p']:.8g}",
            f"- Tropical-minus-North interaction: {primary['interaction_tropical_minus_northern']:.8f}; one-sided negative p={primary['interaction_one_sided_negative_p']:.8g}; two-sided p={primary['interaction_two_sided_p']:.8g}",
            f"- Tropical distance slope: {primary['tropical_distance_slope']:.8f}; two-sided p={primary['tropical_distance_slope_two_sided_p']:.8g}",
        ]
    for name in ("supplemental_only", "no_zero_constant", "equal_island"):
        item = sensitivities[name]
        if item.get("evaluable"):
            lines.append(f"- {name}: North slope={item['northern_distance_slope']:.8f}; interaction={item['interaction_tropical_minus_northern']:.8f}")
        else:
            lines.append(f"- {name}: not evaluable ({item.get('reason', 'unknown')})")
    lines += [
        "",
        f"**Classification:** {decision['classification']}",
        "",
        "This analysis measures experimental pollen limitation, not pollinator abundance decline, visitor identity, historical extinction, or mediation of the floral-trait H2 response.",
    ]
    (output_dir / "RESULT_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.echo(json.dumps(result, indent=2))


if __name__ == "__main__":
    app()
