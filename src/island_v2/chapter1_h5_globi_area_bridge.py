"""Area-conditioned GloBI source-breadth bridge for redesigned Chapter 1 H5.

This module reuses the outcome-blind GloBI V2 island-enrichment table and asks a new,
predeclared question matched to the redesigned H3/H4 ladder: does independently
sampled source-genus interaction breadth show a distance-by-area signal, and does that
moderation differ between northern-midlatitude and tropical contexts?

The analysis uses equal island weight. It does not treat GloBI breadth as true ecological
specialization and cannot identify historical pollinator loss or effective service.
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
    return math.erfc(abs(float(z)) / math.sqrt(2.0))


def _bh(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    out = pd.Series(np.nan, index=series.index, dtype=float)
    valid = values.dropna()
    if valid.empty:
        return out
    ordered = valid.sort_values()
    n = len(ordered)
    adjusted = ordered.to_numpy(float) * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out.loc[ordered.index] = adjusted
    return out


def _z(series: pd.Series) -> np.ndarray:
    x = pd.to_numeric(series, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def _fit_clustered_ols(
    y: np.ndarray,
    design: np.ndarray,
    clusters: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, int]:
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


def _prepare(
    enrichment: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    spec = config["independent_globi_area_bridge"]
    required_enrichment = {
        "island_id",
        "metric",
        "min_independent_references",
        "stratum",
        "source_mode",
        "source_matching",
        "entry_enrichment",
    }
    if missing := required_enrichment - set(enrichment.columns):
        raise typer.BadParameter(f"GloBI enrichment missing columns: {sorted(missing)}")
    required_cov = {
        "island_id",
        "analysis_regime",
        "spatial_block",
        "log_distance_to_continent_km",
        "log_island_area_km2",
        "climate_pc1",
        "climate_pc2",
        "climate_pc3",
        "climate_pc4",
    }
    if missing := required_cov - set(covariates.columns):
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")
    data = enrichment.merge(
        covariates[list(required_cov)].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    data = data.loc[
        data["metric"].astype(str).eq(str(spec["response_metric"]))
        & data["stratum"].astype(str).isin([str(x) for x in spec["strata"]])
        & data["source_mode"].astype(str).isin([str(x) for x in spec["source_modes"]])
        & data["analysis_regime"].astype(str).isin([str(x) for x in spec["contexts"]])
    ].copy()
    numeric = [
        "entry_enrichment",
        "min_independent_references",
        "log_distance_to_continent_km",
        "log_island_area_km2",
        "climate_pc1",
        "climate_pc2",
        "climate_pc3",
        "climate_pc4",
    ]
    for column in numeric:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data = data.dropna(subset=numeric)
    data = data.loc[data["spatial_block"].fillna("").astype(str).ne("")].copy()
    return data


def fit_within_cells(data: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    spec = config["independent_globi_area_bridge"]
    minimum = int(spec["support"]["minimum_islands"])
    rows: list[dict[str, Any]] = []
    group_cols = [
        "stratum",
        "source_mode",
        "source_matching",
        "min_independent_references",
        "analysis_regime",
    ]
    controls = ["climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"]
    for key, part in data.groupby(group_cols, sort=True):
        part = part.drop_duplicates("island_id").copy()
        n_islands = int(part["island_id"].nunique())
        if n_islands < minimum:
            rows.append(
                {
                    **dict(zip(group_cols, key, strict=True)),
                    "status": "not_testable",
                    "n_islands": n_islands,
                }
            )
            continue
        zd = _z(part["log_distance_to_continent_km"])
        za = _z(part["log_island_area_km2"])
        columns = [np.ones(len(part), dtype=float)]
        names = ["intercept"]
        for control in controls:
            columns.append(_z(part[control]))
            names.append(f"z_{control}")
        columns.extend([zd, za, zd * za])
        names.extend(["z_distance", "z_area", "z_distance:z_area"])
        beta, covariance, n_clusters = _fit_clustered_ols(
            part["entry_enrichment"].to_numpy(float),
            np.column_stack(columns),
            part["spatial_block"].to_numpy(str),
        )
        index = names.index("z_distance:z_area")
        estimate = float(beta[index])
        stderr = float(math.sqrt(max(float(covariance[index, index]), 0.0)))
        z_value = estimate / stderr if stderr > 0 else float("nan")
        rows.append(
            {
                **dict(zip(group_cols, key, strict=True)),
                "status": "fit",
                "n_islands": n_islands,
                "n_clusters": int(n_clusters),
                "distance_by_area": estimate,
                "cluster_robust_se": stderr,
                "z_value": z_value,
                "p_value": _normal_two_sided_p(z_value) if math.isfinite(z_value) else float("nan"),
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty and "p_value" in result:
        fit = result["status"].eq("fit")
        result.loc[fit, "q_source_mode_family"] = (
            result.loc[fit]
            .groupby(
                ["stratum", "source_matching", "min_independent_references", "analysis_regime"],
                group_keys=False,
            )["p_value"]
            .transform(_bh)
        )
    return result


def fit_between_cells(data: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    spec = config["independent_globi_area_bridge"]
    contexts = [str(x) for x in spec["contexts"]]
    if len(contexts) != 2:
        raise ValueError("between-context bridge currently requires exactly two contexts")
    context_a, context_b = contexts
    minimum = int(spec["support"]["minimum_islands"])
    rows: list[dict[str, Any]] = []
    group_cols = ["stratum", "source_mode", "source_matching", "min_independent_references"]
    controls = ["climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"]
    for key, part in data.groupby(group_cols, sort=True):
        part = part.drop_duplicates(["island_id", "analysis_regime"]).copy()
        support = part.groupby("analysis_regime")["island_id"].nunique()
        if int(support.get(context_a, 0)) < minimum or int(support.get(context_b, 0)) < minimum:
            rows.append(
                {
                    **dict(zip(group_cols, key, strict=True)),
                    "context_a": context_a,
                    "context_b": context_b,
                    "status": "not_testable",
                    "n_islands_a": int(support.get(context_a, 0)),
                    "n_islands_b": int(support.get(context_b, 0)),
                }
            )
            continue
        indicator = part["analysis_regime"].astype(str).eq(context_b).to_numpy(float)
        zd = _z(part["log_distance_to_continent_km"])
        za = _z(part["log_island_area_km2"])
        columns = [np.ones(len(part), dtype=float), indicator]
        names = ["intercept", f"context[{context_b}]"]
        for control in controls:
            zc = _z(part[control])
            columns.extend([zc, zc * indicator])
            names.extend([f"z_{control}", f"z_{control}:context"])
        columns.extend(
            [
                zd,
                zd * indicator,
                za,
                za * indicator,
                zd * za,
                zd * za * indicator,
            ]
        )
        names.extend(
            [
                "z_distance",
                "z_distance:context",
                "z_area",
                "z_area:context",
                "z_distance:z_area",
                "z_distance:z_area:context",
            ]
        )
        beta, covariance, n_clusters = _fit_clustered_ols(
            part["entry_enrichment"].to_numpy(float),
            np.column_stack(columns),
            part["spatial_block"].to_numpy(str),
        )
        index = names.index("z_distance:z_area:context")
        estimate = float(beta[index])
        stderr = float(math.sqrt(max(float(covariance[index, index]), 0.0)))
        z_value = estimate / stderr if stderr > 0 else float("nan")
        rows.append(
            {
                **dict(zip(group_cols, key, strict=True)),
                "context_a": context_a,
                "context_b": context_b,
                "status": "fit",
                "n_islands_a": int(support.get(context_a, 0)),
                "n_islands_b": int(support.get(context_b, 0)),
                "n_clusters": int(n_clusters),
                "distance_by_area_by_context": estimate,
                "cluster_robust_se": stderr,
                "z_value": z_value,
                "p_value": _normal_two_sided_p(z_value) if math.isfinite(z_value) else float("nan"),
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty and "p_value" in result:
        fit = result["status"].eq("fit")
        result.loc[fit, "q_source_mode_family"] = (
            result.loc[fit]
            .groupby(["stratum", "source_matching", "min_independent_references"], group_keys=False)[
                "p_value"
            ]
            .transform(_bh)
        )
    return result


def classify_primary(within: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    spec = config["independent_globi_area_bridge"]
    primary_matching = str(spec["primary_source_matching"])
    threshold = int(spec["primary_min_independent_references"])
    source_modes = [str(x) for x in spec["source_modes"]]
    alpha = float(spec["support"]["alpha"])
    north = "northern_midlatitude"
    rows: list[dict[str, Any]] = []
    for stratum in [str(x) for x in spec["strata"]]:
        target = within.loc[
            within["stratum"].astype(str).eq(stratum)
            & within["source_matching"].astype(str).eq(primary_matching)
            & within["min_independent_references"].eq(threshold)
            & within["analysis_regime"].astype(str).eq(north)
            & within["status"].eq("fit")
        ].copy()
        all_modes = set(target["source_mode"].astype(str)) == set(source_modes)
        primary_negative_supported = bool(
            all_modes
            and target["distance_by_area"].lt(0).all()
            and target["q_source_mode_family"].le(alpha).all()
        )
        effort_sign_stable = True
        for effort in [
            int(spec["primary_min_independent_references"]),
            *[int(x) for x in spec["sensitivity_min_independent_references"]],
        ]:
            part = within.loc[
                within["stratum"].astype(str).eq(stratum)
                & within["source_matching"].astype(str).eq(primary_matching)
                & within["min_independent_references"].eq(effort)
                & within["analysis_regime"].astype(str).eq(north)
                & within["status"].eq("fit")
            ]
            if set(part["source_mode"].astype(str)) != set(source_modes) or not part[
                "distance_by_area"
            ].lt(0).all():
                effort_sign_stable = False
        promoted = bool(primary_negative_supported and effort_sign_stable)
        rows.append(
            {
                "stratum": stratum,
                "n_primary_source_modes": int(target["source_mode"].nunique()),
                "primary_all_negative": bool(len(target) and target["distance_by_area"].lt(0).all()),
                "primary_all_fdr_supported": bool(
                    len(target) and target["q_source_mode_family"].le(alpha).all()
                ),
                "effort_threshold_sign_stable": bool(effort_sign_stable),
                "promoted": promoted,
                "classification": (
                    "simple_pollination_breadth_filter_supported"
                    if promoted
                    else "simple_pollination_breadth_filter_not_promoted"
                ),
            }
        )
    return pd.DataFrame(rows)


def run_bridge(
    enrichment: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    data = _prepare(enrichment, covariates, config)
    within = fit_within_cells(data, config)
    between = fit_between_cells(data, config)
    classification = classify_primary(within, config)
    return within, between, classification


@app.command("run")
def run(
    enrichment_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if config.get("contract") != "chapter1_h5_pollination_triangulation_v2":
        raise typer.BadParameter("unexpected H5 triangulation contract")
    enrichment = pd.read_csv(enrichment_csv)
    covariates = pd.read_csv(covariates_csv)
    within, between, classification = run_bridge(enrichment, covariates, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    within.to_csv(output_dir / "h5_globi_area_within.csv", index=False)
    between.to_csv(output_dir / "h5_globi_area_between.csv", index=False)
    classification.to_csv(output_dir / "h5_globi_area_primary_classification.csv", index=False)
    manifest = {
        "contract": config["contract"],
        "input_contract": config["independent_globi_area_bridge"]["input_contract"],
        "equal_island_weight": True,
        "n_within_rows": int(len(within)),
        "n_between_rows": int(len(between)),
        "primary_promoted_strata": int(classification.get("promoted", pd.Series(dtype=bool)).sum()),
        "claim_ceiling": config["claim_ceiling"],
    }
    (output_dir / "h5_globi_area_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(classification.to_csv(index=False))


if __name__ == "__main__":
    app()
