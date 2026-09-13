"""Reporting-only transport diagnostics for the frozen Chapter 1 island Search.

This audit quantifies whether Search API failures or fixed-budget truncation are
associated with island geometry, area, or isolation. It is deliberately downstream of
canonical observation-state assignment: it never rewrites rows, changes support tiers,
or affects the N1 gate.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import geopandas as gpd
import pandas as pd
import typer
import yaml

from island_v2.chapter1_nee_island_search_observation import exact_query_wkts

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_OBSERVATION = {
    "island_id",
    "channel_id",
    "observation_state",
    "search_complete",
    "search_truncated",
    "search_error",
}
REQUIRED_COVARIATES = {
    "island_id",
    "area_km2",
    "distance_to_continent_km",
    "log_island_area_km2",
    "log_distance_to_continent_km",
}


def load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("transport audit config must be a mapping")
    if payload.get("contract") != "chapter1_nee_island_search_transport_audit_v1":
        raise typer.BadParameter("unexpected transport audit contract")
    return payload


def _truth(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def _error_class(value: object) -> str:
    text = "" if value is None else str(value).strip()
    lower = text.lower()
    if not text:
        return "none"
    if "query is too long" in lower or "request-uri too large" in lower or "uri too long" in lower:
        return "query_too_long"
    if "429" in lower or "too many requests" in lower:
        return "http_429"
    if "400" in lower or "bad request" in lower:
        return "http_400"
    if "server disconnected" in lower or "remoteprotocolerror" in lower:
        return "server_disconnect"
    if "timeout" in lower or "timed out" in lower:
        return "timeout"
    return "other"


def _decile(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.isna().any():
        raise ValueError("transport audit decile variable contains missing/nonnumeric values")
    rank = numeric.rank(method="average", pct=True)
    return rank.map(lambda value: min(10, max(1, int(math.ceil(float(value) * 10.0)))))


def prepare_audit_frame(
    observation: pd.DataFrame,
    covariates: pd.DataFrame,
    islands: gpd.GeoDataFrame,
) -> pd.DataFrame:
    missing_obs = REQUIRED_OBSERVATION.difference(observation.columns)
    if missing_obs:
        raise ValueError(f"observation table missing columns: {sorted(missing_obs)}")
    missing_cov = REQUIRED_COVARIATES.difference(covariates.columns)
    if missing_cov:
        raise ValueError(f"covariate table missing columns: {sorted(missing_cov)}")
    if "island_id" not in islands.columns:
        raise ValueError("island geometry requires island_id")

    obs = observation.copy()
    obs["island_id"] = obs["island_id"].fillna("").astype(str).str.strip()
    obs["channel_id"] = obs["channel_id"].fillna("").astype(str).str.strip()
    if obs[["island_id", "channel_id"]].duplicated().any():
        raise ValueError("observation table must be unique by island_id x channel_id")

    cov = covariates[list(REQUIRED_COVARIATES)].copy()
    cov["island_id"] = cov["island_id"].fillna("").astype(str).str.strip()
    if cov["island_id"].duplicated().any():
        raise ValueError("covariate table must be unique by island_id")

    geo = islands[["island_id", "geometry"]].copy()
    geo["island_id"] = geo["island_id"].astype(str)
    if geo["island_id"].duplicated().any():
        raise ValueError("island geometry must be unique by island_id")
    geo["exact_polygon_wkt_length"] = [
        int(sum(len(wkt) for wkt in exact_query_wkts(geometry))) for geometry in geo.geometry
    ]
    geo = pd.DataFrame(geo.drop(columns="geometry"))

    frame = obs.merge(cov, on="island_id", how="left", validate="many_to_one")
    frame = frame.merge(geo, on="island_id", how="left", validate="many_to_one")
    required_joined = [
        "area_km2",
        "distance_to_continent_km",
        "log_island_area_km2",
        "log_distance_to_continent_km",
        "exact_polygon_wkt_length",
    ]
    if frame[required_joined].isna().any().any():
        missing = frame.loc[frame[required_joined].isna().any(axis=1), "island_id"].unique().tolist()
        raise ValueError(f"audit geography missing for observation islands: {missing[:10]}")

    frame["search_error_flag"] = frame["search_error"].fillna("").astype(str).str.len().gt(0)
    frame["search_error_class"] = frame["search_error"].map(_error_class)
    frame["search_complete_flag"] = _truth(frame["search_complete"])
    frame["search_truncated_flag"] = _truth(frame["search_truncated"])
    frame["log1p_area_km2"] = pd.to_numeric(frame["log_island_area_km2"], errors="raise")
    frame["log1p_distance_to_continent_km"] = pd.to_numeric(
        frame["log_distance_to_continent_km"], errors="raise"
    )
    frame["exact_polygon_wkt_length"] = pd.to_numeric(
        frame["exact_polygon_wkt_length"], errors="raise"
    )
    return frame


def audit_transport(
    observation: pd.DataFrame,
    covariates: pd.DataFrame,
    islands: gpd.GeoDataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    frame = prepare_audit_frame(observation, covariates, islands)
    channels = [str(value) for value in config["inputs"]["channels"]]
    unexpected = sorted(set(frame["channel_id"]).difference(channels))
    if unexpected:
        raise ValueError(f"unexpected channels in observation table: {unexpected}")

    summary_rows: list[dict[str, Any]] = []
    decile_rows: list[dict[str, Any]] = []
    decile_variables = {
        "exact_polygon_wkt_length": "exact_polygon_wkt_length",
        "log1p_area_km2": "log1p_area_km2",
        "log1p_distance_to_continent_km": "log1p_distance_to_continent_km",
    }

    for channel in channels:
        work = frame.loc[frame["channel_id"].eq(channel)].copy()
        states = work["observation_state"].value_counts().to_dict()
        errors = work["search_error_flag"]
        truncated = work["search_truncated_flag"]
        complete = work["search_complete_flag"]

        def mean_where(column: str, mask: pd.Series) -> float | None:
            if not mask.any():
                return None
            return float(pd.to_numeric(work.loc[mask, column], errors="raise").mean())

        error_classes = work.loc[errors, "search_error_class"].value_counts().to_dict()
        summary_rows.append(
            {
                "channel_id": channel,
                "n_source_available": int(len(work)),
                "n_detected": int(states.get("detected", 0)),
                "n_adequate_non_detection": int(states.get("adequate_non_detection", 0)),
                "n_insufficient_effort": int(states.get("insufficient_effort", 0)),
                "n_unresolved": int(states.get("unresolved", 0)),
                "n_search_errors": int(errors.sum()),
                "n_truncated_searches": int(truncated.sum()),
                "n_complete_searches": int(complete.sum()),
                "search_error_fraction": float(errors.mean()) if len(work) else 0.0,
                "truncation_fraction": float(truncated.mean()) if len(work) else 0.0,
                "n_error_http_400": int(error_classes.get("http_400", 0)),
                "n_error_http_429": int(error_classes.get("http_429", 0)),
                "n_error_query_too_long": int(error_classes.get("query_too_long", 0)),
                "n_error_server_disconnect": int(error_classes.get("server_disconnect", 0)),
                "n_error_timeout": int(error_classes.get("timeout", 0)),
                "n_error_other": int(error_classes.get("other", 0)),
                "mean_log1p_area_error": mean_where("log1p_area_km2", errors),
                "mean_log1p_area_nonerror": mean_where("log1p_area_km2", ~errors),
                "difference_log1p_area_error_minus_nonerror": (
                    None
                    if not errors.any() or not (~errors).any()
                    else float(work.loc[errors, "log1p_area_km2"].mean() - work.loc[~errors, "log1p_area_km2"].mean())
                ),
                "mean_log1p_distance_error": mean_where("log1p_distance_to_continent_km", errors),
                "mean_log1p_distance_nonerror": mean_where("log1p_distance_to_continent_km", ~errors),
                "difference_log1p_distance_error_minus_nonerror": (
                    None
                    if not errors.any() or not (~errors).any()
                    else float(
                        work.loc[errors, "log1p_distance_to_continent_km"].mean()
                        - work.loc[~errors, "log1p_distance_to_continent_km"].mean()
                    )
                ),
                "mean_wkt_length_error": mean_where("exact_polygon_wkt_length", errors),
                "mean_wkt_length_nonerror": mean_where("exact_polygon_wkt_length", ~errors),
            }
        )

        if work.empty:
            continue
        for diagnostic, column in decile_variables.items():
            work["_decile"] = _decile(work[column])
            for decile, group in work.groupby("_decile", sort=True):
                group_errors = group["search_error_flag"]
                decile_rows.append(
                    {
                        "channel_id": channel,
                        "diagnostic": diagnostic,
                        "decile": int(decile),
                        "n_rows": int(len(group)),
                        "n_search_errors": int(group_errors.sum()),
                        "search_error_fraction": float(group_errors.mean()),
                        "min_value": float(pd.to_numeric(group[column], errors="raise").min()),
                        "max_value": float(pd.to_numeric(group[column], errors="raise").max()),
                    }
                )

    summary = pd.DataFrame(summary_rows)
    deciles = pd.DataFrame(decile_rows)
    receipt = {
        "contract": config["contract"],
        "status": config["status"],
        "reporting_only": True,
        "n_observation_rows": int(len(frame)),
        "channels": channels,
        "transport_audit_cannot_promote_or_demote_a_channel": True,
        "transport_audit_cannot_change_N1_pass": True,
        "rerun_or_recode_permitted": False,
    }
    return summary, deciles, receipt


@app.command("run")
def run_command(
    observation_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    islands_gpkg: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(
        Path("config/chapter1_nee_island_search_transport_audit.yml"), exists=True
    ),
) -> None:
    config = load_config(config_path)
    observation = pd.read_csv(observation_csv, dtype=str).fillna("")
    covariates = pd.read_csv(covariates_csv)
    islands = gpd.read_file(islands_gpkg, layer="islands")
    summary, deciles, receipt = audit_transport(observation, covariates, islands, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_dir / "island_search_transport_channel_summary.csv", index=False)
    deciles.to_csv(output_dir / "island_search_transport_decile_diagnostics.csv", index=False)
    (output_dir / "island_search_transport_audit_receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    app()
