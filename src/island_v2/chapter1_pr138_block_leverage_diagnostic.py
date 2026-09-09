"""Diagnose the one influential PR138 spatial block without excluding it.

The leave-one-block analysis identified a single spatial block whose deletion removes
FDR support for the scalar northern classic projection. This module asks *why* that
block is influential by separating sample-size/evidence leverage from ecological
contrast. The focus block remains in the primary analysis.

Outputs are descriptive/influence diagnostics. In particular, within-focus-block
syndrome slopes are deliberately reported without p-values because the focus block is
one spatial cluster and cannot support the primary cluster-robust inferential design.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _quantiles(series: pd.Series) -> dict[str, float]:
    x = _numeric(series).dropna()
    if x.empty:
        return {"min": np.nan, "q25": np.nan, "median": np.nan, "q75": np.nan, "max": np.nan, "sd": np.nan}
    q = x.quantile([0.0, 0.25, 0.5, 0.75, 1.0])
    return {
        "min": float(q.loc[0.0]),
        "q25": float(q.loc[0.25]),
        "median": float(q.loc[0.5]),
        "q75": float(q.loc[0.75]),
        "max": float(q.loc[1.0]),
        "sd": float(x.std(ddof=0)),
    }


def _descriptive_slope(x: pd.Series, y: pd.Series) -> tuple[int, float, float]:
    frame = pd.DataFrame({"x": _numeric(x), "y": _numeric(y)}).dropna()
    if len(frame) < 3:
        return len(frame), np.nan, np.nan
    x_values = frame["x"].to_numpy(float)
    y_values = frame["y"].to_numpy(float)
    x_sd = float(np.std(x_values, ddof=0))
    y_sd = float(np.std(y_values, ddof=0))
    if not np.isfinite(x_sd) or x_sd <= 0:
        return len(frame), np.nan, np.nan
    zx = (x_values - float(np.mean(x_values))) / x_sd
    design = np.column_stack([np.ones(len(frame)), zx])
    beta = np.linalg.pinv(design.T @ design) @ (design.T @ y_values)
    if not np.isfinite(y_sd) or y_sd <= 0:
        correlation = np.nan
    else:
        correlation = float(np.corrcoef(x_values, y_values)[0, 1])
    return len(frame), float(beta[1]), correlation


def _floristic_by_island(status_flora: pd.DataFrame) -> pd.DataFrame:
    required = {"island_id", "accepted_species", "origin_status", "endemic_status"}
    missing = required - set(status_flora.columns)
    if missing:
        raise typer.BadParameter(f"status flora missing columns: {sorted(missing)}")
    work = status_flora.copy()
    work["island_id"] = work["island_id"].astype(str)
    work["accepted_species"] = work["accepted_species"].astype(str)
    native = work.loc[work["origin_status"].fillna("").astype(str).eq("native")].copy()
    if native.empty:
        return pd.DataFrame(columns=["island_id", "native_species_count", "endemic_native_count", "endemic_native_fraction"])
    native["is_endemic"] = native["endemic_status"].fillna("").astype(str).eq("endemic").astype(int)
    summary = (
        native.drop_duplicates(["island_id", "accepted_species"])
        .groupby("island_id", as_index=False)
        .agg(
            native_species_count=("accepted_species", "size"),
            endemic_native_count=("is_endemic", "sum"),
        )
    )
    summary["endemic_native_fraction"] = summary["endemic_native_count"] / summary["native_species_count"]
    return summary


def _evidence_by_island(island_scores: pd.DataFrame, stratum: str) -> pd.DataFrame:
    required = {"island_id", "stratum", "syndrome", "syndrome_score", "n_species"}
    missing = required - set(island_scores.columns)
    if missing:
        raise typer.BadParameter(f"island scores missing columns: {sorted(missing)}")
    work = island_scores.loc[island_scores["stratum"].astype(str).eq(stratum)].copy()
    if work.empty:
        return pd.DataFrame(columns=["island_id", "syndrome_axes_available", "mean_species_support", "median_species_support", "min_species_support"])
    return (
        work.groupby("island_id", as_index=False)
        .agg(
            syndrome_axes_available=("syndrome", "nunique"),
            mean_species_support=("n_species", "mean"),
            median_species_support=("n_species", "median"),
            min_species_support=("n_species", "min"),
        )
    )


def _group_summary(
    islands: pd.DataFrame,
    *,
    group_name: str,
    distance_column: str,
) -> dict[str, Any]:
    row: dict[str, Any] = {"group": group_name, "n_islands": int(islands["island_id"].nunique())}
    for prefix, column in (
        ("latitude", "island_latitude"),
        ("longitude", "island_longitude"),
        ("distance", distance_column),
        ("log_area", "log_island_area_km2"),
        ("native_species", "native_species_count"),
        ("endemic_fraction", "endemic_native_fraction"),
        ("syndrome_axes", "syndrome_axes_available"),
        ("species_support", "mean_species_support"),
    ):
        if column not in islands.columns:
            continue
        q = _quantiles(islands[column])
        for stat, value in q.items():
            row[f"{prefix}_{stat}"] = value
    return row


def build_block_leverage_diagnostic(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    status_flora: pd.DataFrame,
    full_within: pd.DataFrame,
    block_deletion: pd.DataFrame,
    pattern_config: dict[str, Any],
    diagnostic_config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    focus_block = str(diagnostic_config["focus_block"])
    cluster = str(pattern_config["cluster_column"])
    context_column = str(pattern_config["context_column"])
    north = str(diagnostic_config["contexts"]["focus_primary_context"])
    strata = [str(x) for x in diagnostic_config["strata"]]
    distance_column = "distance_to_continent_km" if "distance_to_continent_km" in covariates.columns else str(pattern_config["geography_column"])

    required_cov = {"island_id", cluster, context_column, distance_column}
    missing = required_cov - set(covariates.columns)
    if missing:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")

    cov = covariates.copy()
    cov["island_id"] = cov["island_id"].astype(str)
    cov[cluster] = cov[cluster].fillna("").astype(str)
    cov[context_column] = cov[context_column].fillna("").astype(str)
    cov = cov.drop_duplicates("island_id")
    flora = _floristic_by_island(status_flora)

    group_rows: list[dict[str, Any]] = []
    syndrome_rows: list[dict[str, Any]] = []
    influence_parts: list[pd.DataFrame] = []

    for stratum in strata:
        evidence = _evidence_by_island(island_scores, stratum)
        northern = cov.loc[cov[context_column].eq(north)].copy()
        northern = northern.merge(flora, on="island_id", how="left", validate="one_to_one")
        northern = northern.merge(evidence, on="island_id", how="left", validate="one_to_one")
        focus = northern.loc[northern[cluster].eq(focus_block)].copy()
        rest = northern.loc[northern[cluster].ne(focus_block)].copy()
        scored_ids = set(evidence["island_id"].astype(str))
        focus_scored = focus.loc[focus["island_id"].isin(scored_ids)].copy()
        rest_scored = rest.loc[rest["island_id"].isin(scored_ids)].copy()

        for name, frame in (
            ("focus_all_geography", focus),
            ("focus_syndrome_scored", focus_scored),
            ("other_northern_all_geography", rest),
            ("other_northern_syndrome_scored", rest_scored),
        ):
            row = _group_summary(frame, group_name=name, distance_column=distance_column)
            row["stratum"] = stratum
            row["focus_block"] = focus_block
            group_rows.append(row)

        scores = island_scores.loc[
            island_scores["stratum"].astype(str).eq(stratum)
        ].copy()
        scores = scores.merge(
            cov[["island_id", cluster, context_column, distance_column]],
            on="island_id",
            how="left",
            validate="many_to_one",
        )
        scores = scores.loc[scores[context_column].eq(north)].copy()
        for syndrome, subset in scores.groupby("syndrome", sort=True):
            focus_scores = subset.loc[subset[cluster].eq(focus_block)]
            rest_scores = subset.loc[subset[cluster].ne(focus_block)]
            nf, slope_f, corr_f = _descriptive_slope(
                focus_scores[distance_column], focus_scores["syndrome_score"]
            )
            nr, slope_r, corr_r = _descriptive_slope(
                rest_scores[distance_column], rest_scores["syndrome_score"]
            )
            syndrome_rows.append(
                {
                    "stratum": stratum,
                    "syndrome": str(syndrome),
                    "focus_block": focus_block,
                    "focus_n_islands": nf,
                    "focus_mean_score": float(_numeric(focus_scores["syndrome_score"]).mean()) if nf else np.nan,
                    "focus_distance_slope_per_sd_descriptive": slope_f,
                    "focus_distance_correlation_descriptive": corr_f,
                    "rest_n_islands": nr,
                    "rest_mean_score": float(_numeric(rest_scores["syndrome_score"]).mean()) if nr else np.nan,
                    "rest_distance_slope_per_sd_descriptive": slope_r,
                    "rest_distance_correlation_descriptive": corr_r,
                }
            )

        deletion = block_deletion.loc[block_deletion["stratum"].astype(str).eq(stratum)].copy()
        full_row = full_within.loc[
            full_within["stratum"].astype(str).eq(stratum)
            & full_within["support_tier"].astype(str).eq("confirmatory")
            & full_within["context"].astype(str).eq(north)
        ]
        full_projection = (
            float(full_row.iloc[0]["northern_classic_projection"])
            if not full_row.empty
            else np.nan
        )
        deletion["full_classic_projection"] = full_projection
        deletion["deletion_delta_classic_projection"] = (
            full_projection - _numeric(deletion["north_classic_projection"])
        )
        deletion["absolute_projection_influence"] = deletion["deletion_delta_classic_projection"].abs()
        deletion["influence_rank_desc"] = deletion["absolute_projection_influence"].rank(method="min", ascending=False)
        deletion["block_size_rank_desc"] = _numeric(deletion["deleted_block_islands"]).rank(method="min", ascending=False)
        deletion["influence_percentile"] = deletion["absolute_projection_influence"].rank(method="average", pct=True)
        deletion["size_percentile"] = _numeric(deletion["deleted_block_islands"]).rank(method="average", pct=True)
        deletion["is_focus_block"] = deletion["deleted_block"].astype(str).eq(focus_block)
        influence_parts.append(deletion)

    groups = pd.DataFrame(group_rows)
    syndrome_descriptive = pd.DataFrame(syndrome_rows)
    influence = pd.concat(influence_parts, ignore_index=True) if influence_parts else pd.DataFrame()

    focus_influence = influence.loc[influence["is_focus_block"]].copy() if not influence.empty else pd.DataFrame()
    if not focus_influence.empty:
        focus_influence = focus_influence[
            [
                "stratum",
                "deleted_block",
                "deleted_block_islands",
                "full_classic_projection",
                "north_classic_projection",
                "north_classic_q",
                "deletion_delta_classic_projection",
                "absolute_projection_influence",
                "influence_rank_desc",
                "block_size_rank_desc",
                "influence_percentile",
                "size_percentile",
                "north_large_bee_slope",
                "north_generalized_slope",
                "north_tropical_vector_q",
            ]
        ]
    return groups, syndrome_descriptive, influence, focus_influence


@app.command("run")
def run(
    island_scores_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    status_flora_csv: Path = typer.Option(..., exists=True),
    full_within_csv: Path = typer.Option(..., exists=True),
    block_deletion_csv: Path = typer.Option(..., exists=True),
    pattern_config_path: Path = typer.Option(..., exists=True),
    diagnostic_config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    pattern_config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    diagnostic_config = yaml.safe_load(diagnostic_config_path.read_text(encoding="utf-8"))
    groups, syndrome, influence, focus = build_block_leverage_diagnostic(
        pd.read_csv(island_scores_csv),
        pd.read_csv(covariates_csv),
        pd.read_csv(status_flora_csv),
        pd.read_csv(full_within_csv),
        pd.read_csv(block_deletion_csv),
        pattern_config,
        diagnostic_config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    groups.to_csv(output_dir / "focus_vs_other_northern_summary.csv", index=False)
    syndrome.to_csv(output_dir / "focus_block_syndrome_descriptive.csv", index=False)
    influence.to_csv(output_dir / "all_block_influence_rank.csv", index=False)
    focus.to_csv(output_dir / "focus_block_influence.csv", index=False)
    manifest = {
        "contract": str(diagnostic_config["contract"]),
        "focus_block": str(diagnostic_config["focus_block"]),
        "focus_block_removed_from_primary": False,
        "within_focus_block_inference_claimed": False,
        "n_influence_rows": int(len(influence)),
    }
    (output_dir / "block_leverage_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
