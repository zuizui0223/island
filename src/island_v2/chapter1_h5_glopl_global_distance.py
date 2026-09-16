"""All-site GloPL global mainland-to-island distance extension for Chapter 1 H5.

The preflight reads only coordinates and publication identity, then places every valid
GloPL site on the same continuous distance-to-seeded-major-continent axis already used
by Chapter 1. Outcome columns are read only after the frozen geography/support gates are
recorded. The resulting associations do not identify a causal island effect.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import typer
import yaml
from scipy.stats import chi2
from shapely.geometry import Point
from shapely.ops import unary_union

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT = "chapter1_h5_glopl_global_distance_v1"
CONTEXTS = (
    "northern_midlatitude",
    "northern_high_latitude",
    "tropical",
    "southern_extratropical",
)
MEASUREMENT_COLUMNS = (
    "PL_Effect_Size_Type1",
    "PL_Effect_Size_Type2",
    "Constant_added",
    "Level_of_Supplementation",
)
EFFECT_COLUMNS = ("PL_Effect_Size", *MEASUREMENT_COLUMNS)


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected all-site GloPL contract")
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


def read_preflight_metadata(path: Path, config: dict[str, Any]) -> pd.DataFrame:
    allowed = [str(x) for x in config["preflight"]["allowed_columns"]]
    frame = pd.read_csv(path, usecols=allowed, dtype=str, encoding="latin-1").fillna("")
    return frame[allowed].reset_index(drop=True)


def build_study_key(frame: pd.DataFrame) -> pd.Series:
    doi = frame["DOI"].fillna("").astype(str).str.strip().str.casefold()
    author = frame["Author"].fillna("").astype(str).str.strip()
    year = frame["Year"].fillna("").astype(str).str.strip()
    fallback = pd.Series("", index=frame.index, dtype="object")
    valid_fallback = author.ne("") & year.ne("")
    fallback.loc[valid_fallback] = "author_year:" + author.loc[valid_fallback] + "|" + year.loc[valid_fallback]
    return pd.Series(np.where(doi.ne(""), "doi:" + doi, fallback), index=frame.index, dtype="object")


def assign_context(latitude: pd.Series, config: dict[str, Any]) -> pd.Series:
    lat = pd.to_numeric(latitude, errors="coerce")
    tropic = float(config["contexts"]["tropical_abs_lat_lt"])
    arctic = float(config["contexts"]["northern_high_lat_ge"])
    out = pd.Series("unresolved", index=lat.index, dtype="object")
    out.loc[lat.abs() < tropic] = "tropical"
    out.loc[(lat >= tropic) & (lat < arctic)] = "northern_midlatitude"
    out.loc[lat >= arctic] = "northern_high_latitude"
    out.loc[lat <= -tropic] = "southern_extratropical"
    out.loc[lat.isna()] = "unresolved"
    return out


def select_seeded_continent_union(
    land: gpd.GeoDataFrame,
    seeds: dict[str, Any],
    *,
    projected_crs: str,
):
    if land.crs is None:
        raise typer.BadParameter("Natural Earth land layer has no CRS")
    land_wgs = land.to_crs(4326)
    selected: set[Any] = set()
    for name, coords in seeds.items():
        lon, lat = float(coords[0]), float(coords[1])
        covered = land_wgs.geometry.covers(Point(lon, lat))
        if not bool(covered.any()):
            raise typer.BadParameter(f"Natural Earth land did not contain continent seed {name}")
        selected.add(land_wgs.index[covered][0])
    selected_land = land_wgs.loc[sorted(selected)].to_crs(projected_crs)
    return unary_union(selected_land.geometry.tolist())


def compute_distance_to_continent_union(
    points: gpd.GeoSeries,
    continent_union,
    *,
    projected_crs: str,
) -> np.ndarray:
    if points.crs is None:
        raise typer.BadParameter("GloPL point series has no CRS")
    projected = points.to_crs(projected_crs)
    return projected.distance(continent_union).to_numpy(float) / 1000.0


def _site_key(latitude: pd.Series, longitude: pd.Series) -> pd.Series:
    lat = pd.to_numeric(latitude, errors="coerce")
    lon = pd.to_numeric(longitude, errors="coerce")
    values: list[str] = []
    for la, lo in zip(lat, lon, strict=True):
        if not np.isfinite(la) or not np.isfinite(lo):
            values.append("")
        else:
            values.append(f"{la:.15g}|{lo:.15g}")
    return pd.Series(values, index=latitude.index, dtype="object")


def build_preflight_table(
    metadata: pd.DataFrame,
    land: gpd.GeoDataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    work = metadata.copy().reset_index(drop=True)
    work.insert(0, "row_id", np.arange(len(work), dtype=int))
    work["latitude"] = pd.to_numeric(work["Latitude"], errors="coerce")
    work["longitude"] = pd.to_numeric(work["Longitude"], errors="coerce")
    valid = (
        work["latitude"].between(-90.0, 90.0, inclusive="both")
        & work["longitude"].between(-180.0, 180.0, inclusive="both")
    )
    work["site_key"] = _site_key(work["latitude"], work["longitude"])
    work["study_key"] = build_study_key(work)
    work["analysis_regime"] = assign_context(work["latitude"], config)
    work["distance_to_major_continent_km"] = np.nan
    if bool(valid.any()):
        points = gpd.GeoSeries(
            gpd.points_from_xy(work.loc[valid, "longitude"], work.loc[valid, "latitude"]),
            crs="EPSG:4326",
            index=work.index[valid],
        )
        union = select_seeded_continent_union(
            land,
            config["geography"]["continent_seeds"],
            projected_crs=str(config["geography"]["crs_for_distance"]),
        )
        work.loc[valid, "distance_to_major_continent_km"] = compute_distance_to_continent_union(
            points,
            union,
            projected_crs=str(config["geography"]["crs_for_distance"]),
        )
    work["log1p_distance_to_major_continent_km"] = np.log1p(
        pd.to_numeric(work["distance_to_major_continent_km"], errors="coerce")
    )
    work["valid_coordinate"] = valid
    work["distance_zero"] = pd.to_numeric(
        work["distance_to_major_continent_km"], errors="coerce"
    ).fillna(np.inf).le(1e-9)
    return work


def summarize_preflight(table: pd.DataFrame) -> dict[str, Any]:
    eligible = table.loc[
        table["valid_coordinate"].astype(bool)
        & table["site_key"].astype(str).ne("")
        & table["study_key"].astype(str).ne("")
        & table["analysis_regime"].isin(CONTEXTS)
    ].copy()
    by_context: dict[str, dict[str, int]] = {}
    for context in CONTEXTS:
        part = eligible.loc[eligible["analysis_regime"].eq(context)]
        by_context[context] = {
            "n_rows": int(len(part)),
            "n_sites": int(part["site_key"].nunique()),
            "n_studies": int(part["study_key"].nunique()),
            "n_zero_distance_sites": int(part.loc[part["distance_zero"], "site_key"].nunique()),
            "n_positive_distance_sites": int(part.loc[~part["distance_zero"], "site_key"].nunique()),
        }
    return {
        "n_rows_total": int(len(table)),
        "n_rows_valid_coordinates": int(table["valid_coordinate"].sum()),
        "n_unique_sites": int(eligible["site_key"].nunique()),
        "n_unique_studies": int(eligible["study_key"].nunique()),
        "n_zero_distance_sites": int(eligible.loc[eligible["distance_zero"], "site_key"].nunique()),
        "n_positive_distance_sites": int(eligible.loc[~eligible["distance_zero"], "site_key"].nunique()),
        "by_context": by_context,
    }


def evaluate_preflight_gate(summary: dict[str, Any], config: dict[str, Any]) -> dict[str, Any]:
    gate = config["preflight"]["support_gates"]
    global_pass = (
        int(summary["n_unique_sites"]) >= int(gate["global_min_unique_sites"])
        and int(summary["n_unique_studies"]) >= int(gate["global_min_unique_studies"])
    )
    pair_pass = all(
        int(summary["by_context"][context]["n_sites"])
        >= int(gate["north_tropical_min_unique_sites_per_context"])
        and int(summary["by_context"][context]["n_studies"])
        >= int(gate["north_tropical_min_unique_studies_per_context"])
        for context in ("northern_midlatitude", "tropical")
    )
    four_pass = all(
        int(summary["by_context"][context]["n_sites"])
        >= int(gate["four_context_min_unique_sites_per_context"])
        and int(summary["by_context"][context]["n_studies"])
        >= int(gate["four_context_min_unique_studies_per_context"])
        for context in CONTEXTS
    )
    return {
        "global_gradient_admitted": bool(global_pass),
        "north_tropical_pair_admitted": bool(pair_pass),
        "four_context_heterogeneity_admitted": bool(four_pass),
        "threshold_relaxation_after_preflight": bool(
            config["preflight"]["threshold_relaxation_after_preflight"]
        ),
    }


def read_effect_columns(path: Path, config: dict[str, Any]) -> pd.DataFrame:
    columns = [str(x) for x in config["analysis"]["columns_read_after_preflight_pass"]]
    frame = pd.read_csv(path, usecols=columns, dtype=str, encoding="latin-1").fillna("")
    frame = frame[columns].reset_index(drop=True)
    frame.insert(0, "row_id", np.arange(len(frame), dtype=int))
    return frame


def _truthy(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().str.casefold().isin(
        {"true", "1", "yes", "y", "t"}
    )


def prepare_analysis_rows(
    effects: pd.DataFrame,
    preflight: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    needed_preflight = {
        "row_id",
        "site_key",
        "study_key",
        "analysis_regime",
        "log1p_distance_to_major_continent_km",
        "valid_coordinate",
    }
    if missing := needed_preflight - set(preflight.columns):
        raise typer.BadParameter(f"global GloPL preflight missing columns: {sorted(missing)}")
    if len(effects) != len(preflight):
        raise typer.BadParameter("GloPL effect/preflight row counts differ")
    left = effects.copy()
    right = preflight[list(needed_preflight)].copy()
    left["row_id"] = pd.to_numeric(left["row_id"], errors="raise").astype(int)
    right["row_id"] = pd.to_numeric(right["row_id"], errors="raise").astype(int)
    if left["row_id"].duplicated().any() or right["row_id"].duplicated().any():
        raise typer.BadParameter("global GloPL row_id must be unique")
    work = left.merge(right, on="row_id", how="left", validate="one_to_one")
    work["PL_Effect_Size"] = pd.to_numeric(work["PL_Effect_Size"], errors="coerce")
    finite = np.isfinite(work["PL_Effect_Size"].to_numpy(float))
    valid_coordinate = work["valid_coordinate"].astype(str).str.casefold().isin(
        {"true", "1", "yes"}
    )
    eligible = (
        finite
        & valid_coordinate
        & work["study_key"].fillna("").astype(str).ne("")
        & work["site_key"].fillna("").astype(str).ne("")
        & work["analysis_regime"].isin(CONTEXTS)
    )
    work = work.loc[eligible].copy()
    for column in MEASUREMENT_COLUMNS:
        work[column] = work[column].fillna("").astype(str).str.strip()
    work["Constant_added_bool"] = _truthy(work["Constant_added"])
    return work.reset_index(drop=True)


def aggregate_measurement_cells(rows: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    if rows.empty:
        return pd.DataFrame()
    group_cols = [
        "study_key",
        "site_key",
        "analysis_regime",
        "log1p_distance_to_major_continent_km",
        *MEASUREMENT_COLUMNS,
    ]
    if missing := set(group_cols + ["PL_Effect_Size"]) - set(rows.columns):
        raise typer.BadParameter(f"global GloPL aggregation missing: {sorted(missing)}")
    out = (
        rows.groupby(group_cols, as_index=False, dropna=False)
        .agg(PL_Effect_Size=("PL_Effect_Size", "mean"), n_effect_rows=("PL_Effect_Size", "size"))
        .reset_index(drop=True)
    )
    n_cells = out.groupby("study_key")["site_key"].transform("size").astype(float)
    out["analysis_weight"] = float(config["analysis"]["publication_total_weight"]) / n_cells
    return out


def _measurement_dummies(frame: pd.DataFrame) -> tuple[list[np.ndarray], list[str]]:
    columns: list[np.ndarray] = []
    names: list[str] = []
    for variable in MEASUREMENT_COLUMNS:
        values = frame[variable].fillna("").astype(str)
        levels = sorted(values.unique().tolist())
        for level in levels[1:]:
            columns.append(values.eq(level).to_numpy(float))
            names.append(f"measure:{variable}={level}")
    return columns, names


def _context_dummies(frame: pd.DataFrame) -> tuple[list[np.ndarray], list[str], list[str]]:
    values = frame["analysis_regime"].astype(str)
    reference = "northern_midlatitude"
    levels = [context for context in CONTEXTS if context != reference and bool(values.eq(context).any())]
    columns = [values.eq(level).to_numpy(float) for level in levels]
    names = [f"context_{level}" for level in levels]
    return columns, names, levels


def build_context_design(frame: pd.DataFrame, config: dict[str, Any]) -> tuple[np.ndarray, list[str]]:
    distance = pd.to_numeric(
        frame["z_log1p_distance_to_major_continent_km"], errors="coerce"
    ).to_numpy(float)
    context_cols, context_names, context_levels = _context_dummies(frame)
    measure_cols, measure_names = _measurement_dummies(frame)
    columns: list[np.ndarray] = [np.ones(len(frame)), *context_cols, distance]
    names = ["intercept", *context_names, "z_distance"]
    for col, level in zip(context_cols, context_levels, strict=True):
        columns.append(distance * col)
        names.append(f"z_distance:context_{level}")
    columns.extend(measure_cols)
    names.extend(measure_names)
    return np.column_stack(columns), names


def build_global_design(frame: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    distance = pd.to_numeric(
        frame["z_log1p_distance_to_major_continent_km"], errors="coerce"
    ).to_numpy(float)
    context_cols, context_names, _ = _context_dummies(frame)
    measure_cols, measure_names = _measurement_dummies(frame)
    columns = [np.ones(len(frame)), *context_cols, distance, *measure_cols]
    names = ["intercept", *context_names, "z_distance", *measure_names]
    return np.column_stack(columns), names


def _p2(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0)) if math.isfinite(z) else float("nan")


def _upper(z: float) -> float:
    return 0.5 * math.erfc(z / math.sqrt(2.0)) if math.isfinite(z) else float("nan")


def _lower(z: float) -> float:
    return 0.5 * math.erfc(-z / math.sqrt(2.0)) if math.isfinite(z) else float("nan")


def _clustered_wls(
    frame: pd.DataFrame,
    X: np.ndarray,
    names: list[str],
) -> dict[str, Any]:
    y = pd.to_numeric(frame["PL_Effect_Size"], errors="coerce").to_numpy(float)
    w = pd.to_numeric(frame["analysis_weight"], errors="coerce").to_numpy(float)
    clusters = frame["study_key"].astype(str).to_numpy()
    n, p = X.shape
    n_clusters = int(pd.Series(clusters).nunique())
    if n <= p + 2 or n_clusters <= p + 1:
        return {"evaluable": False, "reason": "insufficient_rows_or_clusters", "n": n, "p": p, "n_clusters": n_clusters}
    if np.linalg.matrix_rank(X) < p:
        return {"evaluable": False, "reason": "design_not_full_rank", "n": n, "p": p, "n_clusters": n_clusters}
    if not np.isfinite(X).all() or not np.isfinite(y).all() or not np.isfinite(w).all() or np.any(w <= 0):
        return {"evaluable": False, "reason": "nonfinite_input"}
    xtwx = X.T @ (w[:, None] * X)
    bread = np.linalg.inv(xtwx)
    beta = bread @ (X.T @ (w * y))
    residual = y - X @ beta
    meat = np.zeros((p, p), dtype=float)
    for cluster in np.unique(clusters):
        mask = clusters == cluster
        score = X[mask].T @ (w[mask] * residual[mask])
        meat += np.outer(score, score)
    cov = bread @ meat @ bread
    cov *= (n_clusters / (n_clusters - 1.0)) * ((n - 1.0) / (n - p))
    if not np.isfinite(cov).all():
        return {"evaluable": False, "reason": "nonfinite_covariance"}
    se = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    return {
        "evaluable": True,
        "beta": beta,
        "covariance": cov,
        "se": se,
        "names": names,
        "n_cells": int(n),
        "n_publications": n_clusters,
        "n_sites": int(frame["site_key"].nunique()),
    }


def _standardize_distance(cells: pd.DataFrame) -> tuple[pd.DataFrame, float, float]:
    out = cells.copy()
    values = pd.to_numeric(out["log1p_distance_to_major_continent_km"], errors="coerce").to_numpy(float)
    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise typer.BadParameter("global GloPL distance has zero/invalid variance")
    out["z_log1p_distance_to_major_continent_km"] = (values - mean) / sd
    return out, mean, sd


def _reuse_distance_scaling(cells: pd.DataFrame, mean: float, sd: float) -> pd.DataFrame:
    out = cells.copy()
    values = pd.to_numeric(out["log1p_distance_to_major_continent_km"], errors="coerce").to_numpy(float)
    out["z_log1p_distance_to_major_continent_km"] = (values - mean) / sd
    return out


def fit_global_gradient(cells: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    support = config["support_after_effect_availability"]
    if len(cells) < int(support["minimum_cells_for_global_fit"]) or cells["study_key"].nunique() < int(
        support["minimum_publications_for_global_fit"]
    ):
        return {"evaluable": False, "reason": "global_effect_support_failed"}
    X, names = build_global_design(cells)
    fit = _clustered_wls(cells, X, names)
    if not fit["evaluable"]:
        return fit
    idx = fit["names"].index("z_distance")
    slope = float(fit["beta"][idx])
    se = float(fit["se"][idx])
    z = slope / se if se > 0 else float("nan")
    return {
        "evaluable": True,
        "distance_slope": slope,
        "distance_slope_se": se,
        "two_sided_p": _p2(z),
        "one_sided_positive_p": _upper(z),
        "n_cells": fit["n_cells"],
        "n_publications": fit["n_publications"],
        "n_sites": fit["n_sites"],
    }


def fit_context_gradient(cells: pd.DataFrame, config: dict[str, Any], *, preflight_four_pass: bool) -> dict[str, Any]:
    support = config["support_after_effect_availability"]
    counts = cells.groupby("analysis_regime").agg(n_cells=("site_key", "size"), n_publications=("study_key", "nunique"))
    effect_support = all(
        context in counts.index
        and int(counts.loc[context, "n_cells"]) >= int(support["minimum_cells_per_four_context"])
        and int(counts.loc[context, "n_publications"]) >= int(support["minimum_publications_per_four_context"])
        for context in CONTEXTS
    )
    X, names = build_context_design(cells, config)
    fit = _clustered_wls(cells, X, names)
    if not preflight_four_pass or not effect_support or not fit["evaluable"]:
        return {
            "evaluable": False,
            "reason": "four_context_support_failed" if (not preflight_four_pass or not effect_support) else fit.get("reason"),
            "effect_support": bool(effect_support),
            "preflight_support": bool(preflight_four_pass),
        }
    name_to_idx = {name: i for i, name in enumerate(fit["names"])}
    base_idx = name_to_idx["z_distance"]
    base = float(fit["beta"][base_idx])
    slopes: dict[str, float] = {"northern_midlatitude": base}
    interaction_indices: list[int] = []
    for context in CONTEXTS:
        if context == "northern_midlatitude":
            continue
        name = f"z_distance:context_{context}"
        idx = name_to_idx[name]
        interaction_indices.append(idx)
        slopes[context] = base + float(fit["beta"][idx])
    b = fit["beta"][interaction_indices]
    cov = fit["covariance"][np.ix_(interaction_indices, interaction_indices)]
    rank = int(np.linalg.matrix_rank(cov))
    if rank <= 0:
        joint_p = float("nan")
        stat = float("nan")
    else:
        stat = float(b.T @ np.linalg.pinv(cov) @ b)
        joint_p = float(chi2.sf(stat, rank))
    return {
        "evaluable": True,
        "context_slopes": slopes,
        "joint_wald": stat,
        "joint_df": rank,
        "joint_p": joint_p,
        "n_cells": fit["n_cells"],
        "n_publications": fit["n_publications"],
        "n_sites": fit["n_sites"],
    }


def fit_north_tropical_pair(cells: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    pair = cells.loc[cells["analysis_regime"].isin(["northern_midlatitude", "tropical"])].copy()
    support = config["support_after_effect_availability"]
    counts = pair.groupby("analysis_regime").agg(n_cells=("site_key", "size"), n_publications=("study_key", "nunique"))
    if any(
        context not in counts.index
        or int(counts.loc[context, "n_cells"]) < int(support["minimum_cells_per_primary_pair_context"])
        or int(counts.loc[context, "n_publications"]) < int(support["minimum_publications_per_primary_pair_context"])
        for context in ("northern_midlatitude", "tropical")
    ):
        return {"evaluable": False, "reason": "pair_effect_support_failed"}
    X, names = build_context_design(pair, config)
    fit = _clustered_wls(pair, X, names)
    if not fit["evaluable"]:
        return fit
    index = {name: i for i, name in enumerate(fit["names"])}
    d = index["z_distance"]
    interaction = index["z_distance:context_tropical"]
    north = float(fit["beta"][d])
    north_se = float(fit["se"][d])
    inter = float(fit["beta"][interaction])
    inter_se = float(fit["se"][interaction])
    north_z = north / north_se if north_se > 0 else float("nan")
    inter_z = inter / inter_se if inter_se > 0 else float("nan")
    tropical = north + inter
    return {
        "evaluable": True,
        "northern_distance_slope": north,
        "northern_distance_slope_se": north_se,
        "northern_one_sided_positive_p": _upper(north_z),
        "northern_two_sided_p": _p2(north_z),
        "interaction_tropical_minus_northern": inter,
        "interaction_se": inter_se,
        "interaction_one_sided_negative_p": _lower(inter_z),
        "interaction_two_sided_p": _p2(inter_z),
        "tropical_distance_slope": tropical,
        "n_cells": fit["n_cells"],
        "n_publications": fit["n_publications"],
        "n_sites": fit["n_sites"],
    }


def classify_global_result(
    global_result: dict[str, Any],
    pair_result: dict[str, Any],
    heterogeneity_result: dict[str, Any],
    sensitivities: dict[str, dict[str, float]],
    config: dict[str, Any],
) -> dict[str, Any]:
    alpha = 0.05
    sensitivity_values = list(sensitivities.values())
    global_direction = bool(sensitivity_values) and all(
        float(result.get("global_distance_slope", float("nan"))) > 0
        for result in sensitivity_values
    )
    pair_direction = bool(sensitivity_values) and all(
        float(result.get("north_slope", float("nan"))) > 0
        and float(result.get("interaction", float("nan"))) < 0
        for result in sensitivity_values
    )
    global_supported = bool(
        global_result.get("evaluable")
        and float(global_result.get("distance_slope", float("nan"))) > 0
        and float(global_result.get("one_sided_positive_p", 1.0)) <= alpha
        and global_direction
    )
    pair_supported = bool(
        pair_result.get("evaluable")
        and float(pair_result.get("northern_distance_slope", float("nan"))) > 0
        and float(pair_result.get("northern_one_sided_positive_p", 1.0)) <= alpha
        and float(pair_result.get("interaction_tropical_minus_northern", float("nan"))) < 0
        and float(pair_result.get("interaction_one_sided_negative_p", 1.0)) <= alpha
        and pair_direction
    )
    heterogeneity_supported = bool(
        heterogeneity_result.get("evaluable")
        and float(heterogeneity_result.get("joint_p", 1.0)) <= alpha
    )
    return {
        "global_gradient_supported": global_supported,
        "north_tropical_specificity_supported": pair_supported,
        "four_context_heterogeneity_supported": heterogeneity_supported,
        "causal_pollination_mechanism_identified": False,
        "classification": (
            "global_and_context_specific_pollen_limitation_gradients_supported"
            if global_supported and pair_supported
            else "global_pollen_limitation_gradient_supported_context_specificity_not_established"
            if global_supported
            else "context_specific_pattern_supported_without_global_gradient"
            if pair_supported
            else "global_pollen_limitation_distance_gradient_not_supported"
        ),
    }


def _analysis_for_cells(
    cells: pd.DataFrame,
    config: dict[str, Any],
    *,
    preflight_four_pass: bool,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    return (
        fit_global_gradient(cells, config),
        fit_north_tropical_pair(cells, config),
        fit_context_gradient(cells, config, preflight_four_pass=preflight_four_pass),
    )


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, allow_nan=True) + "\n", encoding="utf-8")


@app.command("preflight")
def preflight(
    glopl_csv: Path = typer.Option(..., exists=True),
    land_geojson: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    source_sha = validate_source_sha256(glopl_csv, config)
    metadata = read_preflight_metadata(glopl_csv, config)
    land = gpd.read_file(land_geojson)
    table = build_preflight_table(metadata, land, config)
    summary = summarize_preflight(table)
    decision = evaluate_preflight_gate(summary, config)
    result = {
        "contract": CONTRACT,
        "status": "completed_outcome_blind_global_geography_preflight",
        "source_sha256": source_sha,
        "geography_sha256": _sha256(land_geojson),
        "summary": summary,
        "decision": decision,
        "forbidden_outcome_columns_materialized": False,
        "claim_ceiling": config["claim_ceiling"],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_dir / "GLOBAL_SITE_PREFLIGHT.csv.gz", index=False, compression="gzip")
    _write_json(output_dir / "PREFLIGHT.json", result)
    lines = [
        "# GloPL all-site global-distance preflight",
        "",
        f"- rows: {summary['n_rows_total']}",
        f"- valid-coordinate rows: {summary['n_rows_valid_coordinates']}",
        f"- unique sites: {summary['n_unique_sites']}",
        f"- unique studies: {summary['n_unique_studies']}",
        f"- zero-distance sites: {summary['n_zero_distance_sites']}",
        f"- positive-distance sites: {summary['n_positive_distance_sites']}",
        f"- global gate: {decision['global_gradient_admitted']}",
        f"- North-Tropical gate: {decision['north_tropical_pair_admitted']}",
        f"- four-context gate: {decision['four_context_heterogeneity_admitted']}",
    ]
    (output_dir / "PREFLIGHT_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.echo(json.dumps(result, indent=2))


@app.command("analyse")
def analyse(
    glopl_csv: Path = typer.Option(..., exists=True),
    preflight_csv: Path = typer.Option(..., exists=True),
    preflight_json: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    source_sha = validate_source_sha256(glopl_csv, config)
    preflight_result = json.loads(preflight_json.read_text(encoding="utf-8"))
    decision = preflight_result["decision"]
    if not decision.get("global_gradient_admitted") or not decision.get("north_tropical_pair_admitted"):
        raise typer.BadParameter("frozen GloPL global preflight gate did not admit primary analyses")
    preflight_table = pd.read_csv(preflight_csv, dtype=str).fillna("")
    effects = read_effect_columns(glopl_csv, config)
    rows = prepare_analysis_rows(effects, preflight_table, config)
    cells = aggregate_measurement_cells(rows, config)
    cells, distance_mean, distance_sd = _standardize_distance(cells)

    global_result, pair_result, context_result = _analysis_for_cells(
        cells,
        config,
        preflight_four_pass=bool(decision.get("four_context_heterogeneity_admitted")),
    )

    sensitivities: dict[str, dict[str, Any]] = {}
    for name, mask in {
        "supplemental_only": rows["PL_Effect_Size_Type2"].astype(str).eq("Sup"),
        "no_zero_constant": ~rows["Constant_added_bool"].astype(bool),
    }.items():
        sens_rows = rows.loc[mask].copy()
        sens_cells = aggregate_measurement_cells(sens_rows, config)
        sens_cells = _reuse_distance_scaling(sens_cells, distance_mean, distance_sd)
        g, p, h = _analysis_for_cells(
            sens_cells,
            config,
            preflight_four_pass=bool(decision.get("four_context_heterogeneity_admitted")),
        )
        sensitivities[name] = {
            "global": g,
            "pair": p,
            "heterogeneity": h,
            "global_distance_slope": g.get("distance_slope"),
            "north_slope": p.get("northern_distance_slope"),
            "interaction": p.get("interaction_tropical_minus_northern"),
        }

    decision_result = classify_global_result(
        global_result,
        pair_result,
        context_result,
        sensitivities,
        config,
    )
    result = {
        "contract": CONTRACT,
        "status": "completed",
        "source_sha256": source_sha,
        "preflight": preflight_result,
        "n_finite_effect_rows": int(len(rows)),
        "n_measurement_cells": int(len(cells)),
        "n_unique_sites": int(cells["site_key"].nunique()),
        "n_publications": int(cells["study_key"].nunique()),
        "distance_standardization": {"mean": distance_mean, "sd": distance_sd},
        "global_gradient": global_result,
        "north_tropical_pair": pair_result,
        "four_context_heterogeneity": context_result,
        "sensitivities": sensitivities,
        "decision": decision_result,
        "claim_ceiling": config["claim_ceiling"],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    cells.to_csv(output_dir / "PRIMARY_MEASUREMENT_CELLS.csv.gz", index=False, compression="gzip")
    rows.to_csv(output_dir / "ELIGIBLE_EFFECT_ROWS.csv.gz", index=False, compression="gzip")
    _write_json(output_dir / "RESULT.json", result)
    lines = [
        "# GloPL all-site global-distance result",
        "",
        f"- finite effect rows: {len(rows)}",
        f"- measurement cells: {len(cells)}",
        f"- unique sites: {cells['site_key'].nunique()}",
        f"- publications: {cells['study_key'].nunique()}",
        f"- global slope: {global_result.get('distance_slope')}",
        f"- global one-sided p: {global_result.get('one_sided_positive_p')}",
        f"- North slope: {pair_result.get('northern_distance_slope')}",
        f"- North one-sided p: {pair_result.get('northern_one_sided_positive_p')}",
        f"- Tropical-minus-North interaction: {pair_result.get('interaction_tropical_minus_northern')}",
        f"- interaction one-sided p: {pair_result.get('interaction_one_sided_negative_p')}",
        f"- four-context joint p: {context_result.get('joint_p')}",
        f"- classification: {decision_result['classification']}",
    ]
    (output_dir / "RESULT_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.echo(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    app()
