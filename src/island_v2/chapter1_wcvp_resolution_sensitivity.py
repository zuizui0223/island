"""Audit WCVP regional-native H1 against TDWG-L3 geographic resolution.

The concern is that remote oceanic islands may map to geographically small WGSRPD
Level-3 units while near-source islands may inherit much broader Level-3 units.
This can make regional-native compatibility more precise at high isolation.

We quantify Level-3 geographic area from the pinned official TDWG WGSRPD Level-3
GeoJSON and re-fit the regional-native-compatible H1 with log Level-3 area as an
additional continuous covariate on top of the existing island-area and climate controls.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import typer
import yaml
from scipy.stats import spearmanr

from island_v2.chapter1_all_data_probability import (
    build_broad_counts,
    run_probability_analysis,
)
from island_v2.chapter1_wcvp_native_compatibility import (
    build_regional_native_compatible_flora,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def tdwg_l3_areas(level3_geojson: Path) -> pd.DataFrame:
    frame = gpd.read_file(level3_geojson)
    code_col = next(
        (
            col
            for col in ("LEVEL3_COD", "LEVEL3_CODE", "code", "CODE")
            if col in frame.columns
        ),
        None,
    )
    name_col = next(
        (
            col
            for col in ("LEVEL3_NAM", "LEVEL3_NAME", "name", "NAME")
            if col in frame.columns
        ),
        None,
    )
    if code_col is None:
        raise typer.BadParameter(
            f"TDWG L3 code column not found; columns={list(frame.columns)}"
        )
    work = frame[[code_col, *( [name_col] if name_col else [] ), "geometry"]].copy()
    work[code_col] = work[code_col].fillna("").astype(str).str.strip()
    work = work.loc[work[code_col].ne("")].copy()
    work = work.to_crs(6933)
    work["tdwg_l3_area_km2"] = work.geometry.area / 1_000_000.0
    agg = (
        work.groupby(code_col, as_index=False)
        .agg(tdwg_l3_area_km2=("tdwg_l3_area_km2", "sum"))
        .rename(columns={code_col: "tdwg_l3_code"})
    )
    if name_col:
        names = (
            frame[[code_col, name_col]]
            .dropna()
            .drop_duplicates(code_col)
            .rename(
                columns={
                    code_col: "tdwg_l3_code",
                    name_col: "tdwg_l3_name_geometry",
                }
            )
        )
        names["tdwg_l3_code"] = names["tdwg_l3_code"].astype(str).str.strip()
        agg = agg.merge(names, on="tdwg_l3_code", how="left", validate="one_to_one")
    return agg


def build_resolution_covariates(
    covariates: pd.DataFrame,
    island_tdwg: pd.DataFrame,
    areas: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    mapping = island_tdwg[
        ["island_id", "tdwg_l3_code", "tdwg_l3_name", "tdwg_match_status"]
    ].drop_duplicates("island_id").copy()
    mapping["island_id"] = mapping["island_id"].astype(str)
    mapping["tdwg_l3_code"] = mapping["tdwg_l3_code"].fillna("").astype(str)
    mapped = mapping.merge(
        areas,
        on="tdwg_l3_code",
        how="left",
        validate="many_to_one",
    )
    cov = covariates.copy()
    cov["island_id"] = cov["island_id"].astype(str)
    cov = cov.merge(
        mapped[
            [
                "island_id",
                "tdwg_l3_code",
                "tdwg_l3_name",
                "tdwg_match_status",
                "tdwg_l3_area_km2",
            ]
        ],
        on="island_id",
        how="left",
        validate="one_to_one",
    )
    cov["log_tdwg_l3_area_km2"] = np.log1p(
        pd.to_numeric(cov["tdwg_l3_area_km2"], errors="coerce")
    )
    if "log_island_area_km2" in cov.columns:
        cov["log_l3_to_island_area_ratio"] = (
            cov["log_tdwg_l3_area_km2"]
            - pd.to_numeric(cov["log_island_area_km2"], errors="coerce")
        )

    rows: list[dict[str, Any]] = []
    context_col = "analysis_regime"
    for context, part in cov.loc[
        cov["tdwg_match_status"].astype(str).eq("accepted")
    ].groupby(context_col, dropna=False):
        complete = part.dropna(
            subset=["log_distance_to_continent_km", "log_tdwg_l3_area_km2"]
        )
        if len(complete) < 3:
            continue
        rho, p = spearmanr(
            pd.to_numeric(
                complete["log_distance_to_continent_km"], errors="coerce"
            ),
            pd.to_numeric(
                complete["log_tdwg_l3_area_km2"], errors="coerce"
            ),
        )
        ratio_rho = ratio_p = float("nan")
        if "log_l3_to_island_area_ratio" in complete.columns:
            rcomplete = complete.dropna(subset=["log_l3_to_island_area_ratio"])
            if len(rcomplete) >= 3:
                ratio_rho, ratio_p = spearmanr(
                    pd.to_numeric(
                        rcomplete["log_distance_to_continent_km"],
                        errors="coerce",
                    ),
                    pd.to_numeric(
                        rcomplete["log_l3_to_island_area_ratio"],
                        errors="coerce",
                    ),
                )
        rows.append(
            {
                "context": str(context),
                "n_islands": int(len(complete)),
                "spearman_distance_vs_log_l3_area": float(rho),
                "p_distance_vs_log_l3_area": float(p),
                "spearman_distance_vs_log_l3_to_island_ratio": float(ratio_rho),
                "p_distance_vs_log_l3_to_island_ratio": float(ratio_p),
                "median_l3_area_km2": float(
                    complete["tdwg_l3_area_km2"].median()
                ),
            }
        )
    return cov, pd.DataFrame(rows)


def run_resolution_adjusted_h1(
    status_flora: pd.DataFrame,
    state_audit: pd.DataFrame,
    covariates: pd.DataFrame,
    wcvp_ranges: pd.DataFrame,
    island_tdwg: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    dict[str, Any],
]:
    upgraded, audit = build_regional_native_compatible_flora(
        status_flora,
        wcvp_ranges,
        island_tdwg,
    )
    native = upgraded.loc[upgraded["origin_status"].astype(str).eq("native")].copy()

    complete_islands = set(
        covariates.loc[
            pd.to_numeric(
                covariates["log_tdwg_l3_area_km2"], errors="coerce"
            ).notna(),
            "island_id",
        ].astype(str)
    )
    native_complete = native.loc[
        native["island_id"].astype(str).isin(complete_islands)
    ].copy()

    base_cfg = dict(config)
    base_cfg["strata"] = ["all_observed"]
    counts = build_broad_counts(native_complete, state_audit, base_cfg)
    counts["stratum"] = "regional_native_compatible"
    base_cfg["strata"] = ["regional_native_compatible"]

    matched_slopes, _, matched_within, _ = run_probability_analysis(
        counts,
        covariates,
        base_cfg,
    )

    adjusted_cfg = dict(base_cfg)
    baseline = [str(x) for x in adjusted_cfg["baseline_covariates"]]
    if "log_tdwg_l3_area_km2" not in baseline:
        baseline.append("log_tdwg_l3_area_km2")
    adjusted_cfg["baseline_covariates"] = baseline
    adjusted_slopes, _, adjusted_within, _ = run_probability_analysis(
        counts,
        covariates,
        adjusted_cfg,
    )
    return (
        matched_slopes,
        matched_within,
        adjusted_slopes,
        adjusted_within,
        audit,
    )


@app.command()
def main(
    status_flora_csv: Path = typer.Option(..., exists=True),
    state_audit_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    wcvp_ranges_csv: Path = typer.Option(..., exists=True),
    island_tdwg_csv: Path = typer.Option(..., exists=True),
    tdwg_level3_geojson: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
) -> None:
    areas = tdwg_l3_areas(tdwg_level3_geojson)
    covariates, correlation = build_resolution_covariates(
        pd.read_csv(covariates_csv),
        pd.read_csv(island_tdwg_csv),
        areas,
    )
    (
        matched_slopes,
        matched_omnibus,
        adjusted_slopes,
        adjusted_omnibus,
        audit,
    ) = run_resolution_adjusted_h1(
        pd.read_csv(status_flora_csv),
        pd.read_csv(state_audit_csv),
        covariates,
        pd.read_csv(wcvp_ranges_csv),
        pd.read_csv(island_tdwg_csv),
        yaml.safe_load(config_path.read_text(encoding="utf-8")),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    for frame in (
        matched_slopes,
        matched_omnibus,
        adjusted_slopes,
        adjusted_omnibus,
        correlation,
    ):
        if not frame.empty:
            frame.insert(0, "evidence_scope", evidence_scope)
    areas.to_csv(output_dir / "tdwg_l3_area.csv", index=False)
    correlation.to_csv(output_dir / "resolution_distance_audit.csv", index=False)
    matched_slopes.to_csv(
        output_dir / "resolution_matched_baseline_within_slopes.csv",
        index=False,
    )
    matched_omnibus.to_csv(
        output_dir / "resolution_matched_baseline_within_omnibus.csv",
        index=False,
    )
    adjusted_slopes.to_csv(
        output_dir / "resolution_adjusted_within_slopes.csv",
        index=False,
    )
    adjusted_omnibus.to_csv(
        output_dir / "resolution_adjusted_within_omnibus.csv",
        index=False,
    )
    pd.DataFrame([{"evidence_scope": evidence_scope, **audit}]).to_csv(
        output_dir / "regional_native_coverage.csv",
        index=False,
    )
    typer.echo(correlation.to_csv(index=False))
    typer.echo(matched_omnibus.to_csv(index=False))
    typer.echo(adjusted_omnibus.to_csv(index=False))


if __name__ == "__main__":
    app()
