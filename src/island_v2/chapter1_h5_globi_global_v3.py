"""Global four-context GloBI source-breadth analysis for redesigned Chapter 1 H5.

The frozen GloBI genus predictor is unchanged. This module adds an all-observed
contemporary-flora layer across the four H1-H3 analysis regimes, while reusing the
frozen source-backed native enrichment as a biological sensitivity. All-observed
results are explicitly observational and are never promoted as native assembly.
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
from scipy import sparse

from island_v2.chapter1_all_data_probability import _bh, _chi_square_sf_integer_df
from island_v2.chapter1_globi_source_breadth_v2 import (
    compute_effort_matched_enrichment,
    reference_effort_bins,
)
from island_v2.chapter1_pr138_lineage_representation_bridge import (
    _availability_matrices,
    _genus,
    _source_assignment_matrix,
    broad_source_availability,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _normal_two_sided_p(z: float) -> float:
    return math.erfc(abs(float(z)) / math.sqrt(2.0))


def _z(values: pd.Series) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
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
    X = np.asarray(design, dtype=float)
    beta = np.linalg.pinv(X.T @ X) @ (X.T @ y)
    residual = y - X @ beta
    bread = np.linalg.pinv(X.T @ X)
    labels = np.asarray(clusters).astype(str)
    unique = np.unique(labels)
    meat = np.zeros((X.shape[1], X.shape[1]), dtype=float)
    for label in unique:
        mask = labels == label
        score = X[mask].T @ residual[mask]
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    n = len(y)
    k = X.shape[1]
    g = len(unique)
    if g > 1 and n > k:
        covariance *= (g / (g - 1.0)) * ((n - 1.0) / (n - k))
    return beta, covariance, g


def _joint(vector: np.ndarray, covariance: np.ndarray) -> tuple[float, int, float]:
    rank = int(np.linalg.matrix_rank(covariance))
    if rank <= 0:
        return float("nan"), 0, float("nan")
    statistic = float(vector @ np.linalg.pinv(covariance) @ vector)
    return statistic, rank, _chi_square_sf_integer_df(statistic, rank)


def _all_observed_count_matrix(
    status_flora: pd.DataFrame,
    island_index: dict[str, int],
    genus_index: dict[str, int],
) -> sparse.csr_matrix:
    required = {"island_id", "accepted_species"}
    if missing := required - set(status_flora.columns):
        raise typer.BadParameter(f"status flora missing columns: {sorted(missing)}")
    work = status_flora[["island_id", "accepted_species"]].drop_duplicates().copy()
    work["island_id"] = work["island_id"].astype(str)
    work["genus"] = work["accepted_species"].astype(str).map(_genus)
    work = work.loc[
        work["island_id"].isin(island_index) & work["genus"].isin(genus_index)
    ].copy()
    grouped = work.groupby(["island_id", "genus"], as_index=False).agg(
        n_island_species=("accepted_species", "nunique")
    )
    rows = grouped["island_id"].map(island_index).to_numpy(int)
    cols = grouped["genus"].map(genus_index).to_numpy(int)
    return sparse.csr_matrix(
        (grouped["n_island_species"].to_numpy(float), (rows, cols)),
        shape=(len(island_index), len(genus_index)),
    )


def build_all_observed_enrichment(
    genus_breadth: pd.DataFrame,
    gift_flora: pd.DataFrame,
    assignments: pd.DataFrame,
    status_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    metric = str(config["response_metric"])
    thresholds = [int(x) for x in config["reference_thresholds"]]
    source_modes = [str(x) for x in config["source_modes"]]
    matchings = [
        str(config["primary_source_matching"]),
        str(config["sensitivity_source_matching"]),
    ]
    minimum = int(config["minimum_represented_genera"])
    islands = sorted(covariates["island_id"].astype(str).unique())
    island_index = {island: index for index, island in enumerate(islands)}
    parts: list[pd.DataFrame] = []

    for threshold in thresholds:
        eligible = genus_breadth.copy()
        eligible["n_independent_references"] = pd.to_numeric(
            eligible["n_independent_references"], errors="coerce"
        )
        eligible[metric] = pd.to_numeric(eligible[metric], errors="coerce")
        eligible = eligible.loc[
            eligible["n_independent_references"].ge(threshold)
        ].dropna(subset=[metric]).copy()
        genera = sorted(set(eligible["genus"].astype(str)) - {""})
        if not genera:
            continue
        genus_index = {genus: index for index, genus in enumerate(genera)}
        ordered = eligible.drop_duplicates("genus").set_index("genus").loc[genera]
        positions = ordered[metric].to_numpy(float)
        effort_bins = reference_effort_bins(
            ordered["n_independent_references"].to_numpy(float)
        )
        availability = broad_source_availability(gift_flora, set(genera))
        entities = sorted(
            set(pd.to_numeric(assignments["entity_ID"], errors="coerce").dropna().astype(int))
            | set(pd.to_numeric(availability["entity_ID"], errors="coerce").dropna().astype(int))
        )
        entity_index = {entity: index for index, entity in enumerate(entities)}
        presence, richness = _availability_matrices(availability, entity_index, genus_index)
        counts = _all_observed_count_matrix(status_flora, island_index, genus_index).toarray()

        for source_mode in source_modes:
            assignment = _source_assignment_matrix(
                assignments, island_index, entity_index, source_mode=source_mode
            )
            prevalence = (assignment @ presence).toarray().astype(np.int16)
            source_richness = (assignment @ richness).toarray().astype(np.float32)
            for matching in matchings:
                # The non-effort-matched sensitivity is only needed at the primary
                # reference threshold; this keeps the global extension bounded.
                if matching != str(config["primary_source_matching"]) and threshold != int(
                    config["primary_reference_threshold"]
                ):
                    continue
                rows: list[dict[str, Any]] = []
                for island_position, island_id in enumerate(islands):
                    result = compute_effort_matched_enrichment(
                        prevalence[island_position],
                        source_richness[island_position],
                        counts[island_position],
                        positions,
                        effort_bins,
                        matching=matching,
                        minimum_represented_genera=minimum,
                    )
                    if result is None:
                        continue
                    rows.append(
                        {
                            "island_id": island_id,
                            "metric": metric,
                            "min_independent_references": threshold,
                            "stratum": "all_observed",
                            "source_mode": source_mode,
                            "source_matching": matching,
                            **result,
                        }
                    )
                if rows:
                    parts.append(pd.DataFrame(rows))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def _prepare_model_data(
    enrichment: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
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
        data["metric"].astype(str).eq(str(config["response_metric"]))
        & data["analysis_regime"].astype(str).isin([str(x) for x in config["contexts"]])
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


def fit_within_context(data: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    minimum = int(config["model"]["pilot_min_islands"])
    confirmatory = int(config["model"]["confirmatory_min_islands"])
    controls = ["climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"]
    groups = [
        "stratum",
        "source_mode",
        "source_matching",
        "min_independent_references",
        "analysis_regime",
    ]
    rows: list[dict[str, Any]] = []
    for key, part in data.groupby(groups, sort=True):
        part = part.drop_duplicates("island_id").copy()
        n = int(part["island_id"].nunique())
        base = {**dict(zip(groups, key, strict=True)), "n_islands": n}
        if n < minimum:
            rows.append({**base, "status": "not_testable", "support_class": "below_pilot"})
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
        beta, covariance, clusters = _fit_clustered_ols(
            part["entry_enrichment"].to_numpy(float),
            np.column_stack(columns),
            part["spatial_block"].to_numpy(str),
        )
        d = names.index("z_distance")
        da = names.index("z_distance:z_area")
        d_se = float(math.sqrt(max(float(covariance[d, d]), 0.0)))
        da_se = float(math.sqrt(max(float(covariance[da, da]), 0.0)))
        d_z = float(beta[d] / d_se) if d_se > 0 else float("nan")
        da_z = float(beta[da] / da_se) if da_se > 0 else float("nan")
        rows.append(
            {
                **base,
                "status": "fit",
                "support_class": "confirmatory" if n >= confirmatory else "pilot",
                "n_clusters": int(clusters),
                "distance_slope": float(beta[d]),
                "distance_se": d_se,
                "distance_p": _normal_two_sided_p(d_z),
                "distance_by_area": float(beta[da]),
                "distance_by_area_se": da_se,
                "distance_by_area_p": _normal_two_sided_p(da_z),
            }
        )
    out = pd.DataFrame(rows)
    fit = out["status"].eq("fit") if not out.empty else pd.Series(dtype=bool)
    if not out.empty and fit.any():
        family = [
            "stratum",
            "source_matching",
            "min_independent_references",
            "analysis_regime",
        ]
        out.loc[fit, "distance_q"] = (
            out.loc[fit].groupby(family, group_keys=False)["distance_p"].transform(_bh)
        )
        out.loc[fit, "distance_by_area_q"] = (
            out.loc[fit]
            .groupby(family, group_keys=False)["distance_by_area_p"]
            .transform(_bh)
        )
    return out


def fit_global_heterogeneity(data: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    contexts = [str(x) for x in config["contexts"]]
    minimum = int(config["model"]["confirmatory_min_islands"])
    controls = ["climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"]
    groups = ["stratum", "source_mode", "source_matching", "min_independent_references"]
    rows: list[dict[str, Any]] = []
    for key, part in data.groupby(groups, sort=True):
        part = part.drop_duplicates(["island_id", "analysis_regime"]).copy()
        support = part.groupby("analysis_regime")["island_id"].nunique()
        if any(int(support.get(context, 0)) < minimum for context in contexts):
            rows.append(
                {
                    **dict(zip(groups, key, strict=True)),
                    "status": "not_testable",
                    **{f"n_{context}": int(support.get(context, 0)) for context in contexts},
                }
            )
            continue
        columns: list[np.ndarray] = []
        names: list[str] = []
        distance_indices: list[int] = []
        moderation_indices: list[int] = []
        for context in contexts:
            mask = part["analysis_regime"].astype(str).eq(context).to_numpy()
            indicator = mask.astype(float)
            columns.append(indicator)
            names.append(f"{context}:intercept")
            local = part.loc[mask]
            for control in controls:
                z = np.zeros(len(part), dtype=float)
                z[mask] = _z(local[control])
                columns.append(z)
                names.append(f"{context}:z_{control}")
            zd = np.zeros(len(part), dtype=float)
            za = np.zeros(len(part), dtype=float)
            zd[mask] = _z(local["log_distance_to_continent_km"])
            za[mask] = _z(local["log_island_area_km2"])
            columns.extend([zd, za, zd * za])
            names.extend(
                [f"{context}:z_distance", f"{context}:z_area", f"{context}:z_distance:z_area"]
            )
            distance_indices.append(len(names) - 3)
            moderation_indices.append(len(names) - 1)
        beta, covariance, clusters = _fit_clustered_ols(
            part["entry_enrichment"].to_numpy(float),
            np.column_stack(columns),
            part["spatial_block"].to_numpy(str),
        )
        distance = beta[distance_indices]
        moderation = beta[moderation_indices]
        distance_cov = covariance[np.ix_(distance_indices, distance_indices)]
        moderation_cov = covariance[np.ix_(moderation_indices, moderation_indices)]
        # Equality contrasts against the first context.
        contrast = np.zeros((len(contexts) - 1, len(contexts)), dtype=float)
        for row_index in range(len(contexts) - 1):
            contrast[row_index, 0] = -1.0
            contrast[row_index, row_index + 1] = 1.0
        d_diff = contrast @ distance
        d_diff_cov = contrast @ distance_cov @ contrast.T
        da_diff = contrast @ moderation
        da_diff_cov = contrast @ moderation_cov @ contrast.T
        d_stat, d_df, d_p = _joint(d_diff, d_diff_cov)
        da_stat, da_df, da_p = _joint(da_diff, da_diff_cov)
        d0_stat, d0_df, d0_p = _joint(distance, distance_cov)
        da0_stat, da0_df, da0_p = _joint(moderation, moderation_cov)
        rows.append(
            {
                **dict(zip(groups, key, strict=True)),
                "status": "fit",
                **{f"n_{context}": int(support.get(context, 0)) for context in contexts},
                "n_clusters": int(clusters),
                "distance_vector": "|".join(f"{x:.8g}" for x in distance),
                "moderation_vector": "|".join(f"{x:.8g}" for x in moderation),
                "distance_nonzero_wald": d0_stat,
                "distance_nonzero_df": d0_df,
                "distance_nonzero_p": d0_p,
                "distance_heterogeneity_wald": d_stat,
                "distance_heterogeneity_df": d_df,
                "distance_heterogeneity_p": d_p,
                "moderation_nonzero_wald": da0_stat,
                "moderation_nonzero_df": da0_df,
                "moderation_nonzero_p": da0_p,
                "moderation_heterogeneity_wald": da_stat,
                "moderation_heterogeneity_df": da_df,
                "moderation_heterogeneity_p": da_p,
            }
        )
    out = pd.DataFrame(rows)
    fit = out["status"].eq("fit") if not out.empty else pd.Series(dtype=bool)
    if not out.empty and fit.any():
        family = ["stratum", "source_matching", "min_independent_references"]
        for pcol, qcol in [
            ("distance_heterogeneity_p", "distance_heterogeneity_q"),
            ("moderation_heterogeneity_p", "moderation_heterogeneity_q"),
            ("distance_nonzero_p", "distance_nonzero_q"),
            ("moderation_nonzero_p", "moderation_nonzero_q"),
        ]:
            out.loc[fit, qcol] = (
                out.loc[fit].groupby(family, group_keys=False)[pcol].transform(_bh)
            )
    return out


def classify(within: pd.DataFrame, global_tests: pd.DataFrame, config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    alpha = float(config["multiplicity"]["alpha"])
    threshold = int(config["primary_reference_threshold"])
    matching = str(config["primary_source_matching"])
    modes = [str(x) for x in config["source_modes"]]
    contexts = [str(x) for x in config["contexts"]]
    scope = str(config["flora_scopes"]["primary"])
    target = within.loc[
        within["stratum"].eq(scope)
        & within["source_matching"].eq(matching)
        & within["min_independent_references"].eq(threshold)
        & within["status"].eq("fit")
    ].copy()
    context_rows: list[dict[str, Any]] = []
    for context in contexts:
        part = target.loc[target["analysis_regime"].eq(context)]
        modes_complete = set(part["source_mode"].astype(str)) == set(modes)
        distance_supported = bool(
            modes_complete
            and part["distance_slope"].gt(0).all()
            and part["distance_q"].le(alpha).all()
        )
        moderation_supported = bool(
            modes_complete
            and part["distance_by_area"].lt(0).all()
            and part["distance_by_area_q"].le(alpha).all()
        )
        context_rows.append(
            {
                "context": context,
                "n_source_modes": int(part["source_mode"].nunique()),
                "all_positive_distance": bool(len(part) and part["distance_slope"].gt(0).all()),
                "distance_filter_promoted": distance_supported,
                "all_negative_distance_by_area": bool(
                    len(part) and part["distance_by_area"].lt(0).all()
                ),
                "small_island_amplification_promoted": moderation_supported,
            }
        )
    context_class = pd.DataFrame(context_rows)
    universal = bool(len(context_class) == len(contexts) and context_class["distance_filter_promoted"].all())
    small_island = bool(
        len(context_class) == len(contexts)
        and context_class["small_island_amplification_promoted"].all()
    )
    g = global_tests.loc[
        global_tests["stratum"].eq(scope)
        & global_tests["source_matching"].eq(matching)
        & global_tests["min_independent_references"].eq(threshold)
        & global_tests["status"].eq("fit")
    ]
    modes_complete = set(g["source_mode"].astype(str)) == set(modes)
    distance_heterogeneity = bool(
        modes_complete and g["distance_heterogeneity_q"].le(alpha).all()
    )
    moderation_heterogeneity = bool(
        modes_complete and g["moderation_heterogeneity_q"].le(alpha).all()
    )
    overall = pd.DataFrame(
        [
            {
                "flora_scope": scope,
                "universal_positive_breadth_filter": universal,
                "universal_small_island_amplification": small_island,
                "distance_context_heterogeneity_robust": distance_heterogeneity,
                "moderation_context_heterogeneity_robust": moderation_heterogeneity,
                "mechanism_promoted": bool(universal),
                "classification": (
                    "global_source_breadth_filter_supported"
                    if universal
                    else "global_source_breadth_filter_not_promoted"
                ),
            }
        ]
    )
    return context_class, overall


def run_global(
    genus_breadth: pd.DataFrame,
    native_enrichment: pd.DataFrame,
    gift_flora: pd.DataFrame,
    assignments: pd.DataFrame,
    status_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> dict[str, pd.DataFrame]:
    observed = build_all_observed_enrichment(
        genus_breadth, gift_flora, assignments, status_flora, covariates, config
    )
    native = native_enrichment.loc[
        native_enrichment["metric"].astype(str).eq(str(config["response_metric"]))
    ].copy()
    enrichment = pd.concat([observed, native], ignore_index=True)
    model_data = _prepare_model_data(enrichment, covariates, config)
    within = fit_within_context(model_data, config)
    global_tests = fit_global_heterogeneity(model_data, config)
    context_class, overall = classify(within, global_tests, config)
    return {
        "all_observed_enrichment": observed,
        "within_context": within,
        "global_tests": global_tests,
        "context_classification": context_class,
        "overall_classification": overall,
    }


@app.command("run")
def run(
    genus_breadth_csv: Path = typer.Option(..., exists=True),
    native_enrichment_csv: Path = typer.Option(..., exists=True),
    gift_flora_csv: Path = typer.Option(..., exists=True),
    source_assignments_csv: Path = typer.Option(..., exists=True),
    status_flora_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if config.get("contract") != "chapter1_h5_globi_global_v3":
        raise typer.BadParameter("unexpected global GloBI contract")
    outputs = run_global(
        pd.read_csv(genus_breadth_csv),
        pd.read_csv(native_enrichment_csv),
        pd.read_csv(gift_flora_csv),
        pd.read_csv(source_assignments_csv),
        pd.read_csv(status_flora_csv),
        pd.read_csv(covariates_csv),
        config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, frame in outputs.items():
        suffix = ".csv.gz" if len(frame) > 100_000 else ".csv"
        frame.to_csv(output_dir / f"{name}{suffix}", index=False)
    manifest = {
        "contract": config["contract"],
        "predictor_contract": config["predictor_contract"],
        "contexts": config["contexts"],
        "primary_flora_scope": config["flora_scopes"]["primary"],
        "native_sensitivities": config["flora_scopes"]["sensitivities"],
        "equal_island_weight": True,
        "all_observed_is_native_assembly": False,
        "n_all_observed_enrichment_rows": int(len(outputs["all_observed_enrichment"])),
    }
    (output_dir / "chapter1_h5_globi_global_v3_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
