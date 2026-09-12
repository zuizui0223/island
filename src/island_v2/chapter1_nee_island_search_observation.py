"""Credential-independent exact-island GBIF Search acquisition for N1 channels.

This module is a conservative operational alternative to a complete GBIF bulk-download
campaign. It can establish a target detection from a quality-filtered exact-island
record. It can establish adequate non-detection only when every frozen background
subquery is exhaustively retrieved before the shared record budget is exhausted.
Truncation can therefore remove information but can never manufacture disruption.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import geopandas as gpd
import httpx
import pandas as pd
import typer
import yaml

from island_v2.chapter1_nee_channel_observation import (
    OUTPUT_COLUMNS as OBSERVATION_COLUMNS,
    classify_channel_observations,
    load_policy,
    target_species,
)
from island_v2.gbif_blocks import canonical_gbif_polygon, resolve_query_taxon, stable_wkt

app = typer.Typer(add_completion=False, no_args_is_help=True)

SEARCH_AUDIT_COLUMNS = OBSERVATION_COLUMNS + [
    "search_complete",
    "search_truncated",
    "search_error",
    "raw_records_examined",
    "background_subqueries_completed",
    "background_subqueries_total",
    "geometry_components",
]

INTRODUCED_TOKENS = ("INTRODUCED", "INVASIVE", "ALIEN", "NATURALISED", "NATURALIZED")


@dataclass(frozen=True)
class SearchMeta:
    complete: bool
    truncated: bool
    error: str
    raw_records_examined: int
    subqueries_completed: int
    subqueries_total: int
    geometry_components: int


def load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("island Search acquisition config must be a mapping")
    if payload.get("contract") != "chapter1_nee_island_search_acquisition_v1":
        raise typer.BadParameter("unexpected island Search acquisition contract")
    return payload


def _text(value: object) -> str:
    return "" if value is None else str(value).strip()


def eligible_islands(source_availability: pd.DataFrame, channel_id: str) -> list[str]:
    required = {"island_id", "channel_id", "source_state"}
    missing = required.difference(source_availability.columns)
    if missing:
        raise ValueError(f"source availability missing columns: {sorted(missing)}")
    work = source_availability.copy()
    for column in required:
        work[column] = work[column].fillna("").astype(str).str.strip()
    channel_rows = work.loc[work["channel_id"].eq(channel_id)]
    if channel_rows["island_id"].duplicated().any():
        raise ValueError("source availability must be unique by island_id x channel_id")
    return sorted(
        channel_rows.loc[channel_rows["source_state"].eq("available"), "island_id"].unique()
    )


def exact_query_wkts(geometry: object) -> list[str]:
    canonical = canonical_gbif_polygon(geometry)
    if canonical.geom_type == "Polygon":
        return [stable_wkt(canonical)]
    if canonical.geom_type != "MultiPolygon":
        raise ValueError(f"unsupported exact island geometry: {canonical.geom_type}")
    parts = list(canonical.geoms)
    parts.sort(key=lambda part: (-float(part.area), stable_wkt(part)))
    return [stable_wkt(part) for part in parts]


def _known_introduced(value: object) -> bool:
    text = _text(value).upper()
    return any(token in text for token in INTRODUCED_TOKENS)


def _usable_record(record: dict[str, Any]) -> bool:
    if _text(record.get("occurrenceStatus")).upper() == "ABSENT":
        return False
    if _text(record.get("basisOfRecord")).upper() == "LIVING_SPECIMEN":
        return False
    establishment = _text(record.get("establishmentMeans")).upper()
    if "CAPTIVE" in establishment or "CULTIVATED" in establishment:
        return False
    return True


def _record_row(record: dict[str, Any], island_id: str) -> dict[str, Any] | None:
    if not _usable_record(record):
        return None
    return {
        "island_id": island_id,
        "species": _text(record.get("species")),
        "dataset_key": _text(record.get("datasetKey") or record.get("dataset_key")),
        "year": record.get("year", ""),
        "decimal_latitude": record.get("decimalLatitude", ""),
        "decimal_longitude": record.get("decimalLongitude", ""),
        "establishment_means": _text(record.get("establishmentMeans")),
        "gbif_id": _text(record.get("key") or record.get("gbifID")),
    }


def _quality_target_detected(records: list[dict[str, Any]], targets: set[str]) -> bool:
    for record in records:
        if not _usable_record(record):
            continue
        species = _text(record.get("species"))
        if species not in targets:
            continue
        if _known_introduced(record.get("establishmentMeans")):
            continue
        try:
            lat = float(record.get("decimalLatitude"))
            lon = float(record.get("decimalLongitude"))
            year = int(record.get("year"))
        except (TypeError, ValueError):
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180 and 1500 <= year <= 2100):
            continue
        return True
    return False


def _get_json_with_retry(
    client: httpx.Client,
    endpoint: str,
    *,
    params: dict[str, Any],
    retries: int,
    backoff_seconds: float,
) -> dict[str, Any]:
    error: Exception | None = None
    for attempt in range(retries):
        try:
            response = client.get(endpoint, params=params)
            response.raise_for_status()
            payload = response.json()
            if not isinstance(payload, dict):
                raise ValueError("GBIF Search response is not an object")
            return payload
        except Exception as exc:  # noqa: BLE001 - network boundary is fail-closed upstream
            error = exc
            if attempt + 1 < retries:
                time.sleep(backoff_seconds * (attempt + 1))
    raise RuntimeError(f"GBIF island Search failed after {retries} attempts: {error}")


def resolve_background_queries(
    client: httpx.Client,
    channel_id: str,
    config: dict[str, Any],
) -> list[tuple[str, str]]:
    if channel_id not in config["background_queries"]:
        raise ValueError(f"unregistered Search background channel: {channel_id}")
    resolved: list[tuple[str, str]] = []
    for spec in config["background_queries"][channel_id]:
        taxon = str(spec["taxon"])
        rank = str(spec["rank"])
        match = resolve_query_taxon(client, taxon, rank, None, None)
        resolved.append((taxon, str(match["usage_key"])))
    return resolved


def scan_island(
    client: httpx.Client,
    *,
    island_id: str,
    geometry_wkts: list[str],
    background_queries: list[tuple[str, str]],
    targets: set[str],
    config: dict[str, Any],
) -> tuple[pd.DataFrame, SearchMeta]:
    """Retrieve a bounded exact-island background with explicit completeness semantics."""
    if not geometry_wkts:
        raise ValueError("exact island query requires at least one geometry component")
    search = config["search"]
    endpoint = str(search["endpoint"])
    page_size = int(search["page_size"])
    max_records = int(search["max_records_examined_per_island_channel"])
    retries = int(search["retries"])
    backoff = float(search["retry_backoff_seconds"])

    raw_examined = 0
    subqueries_completed = 0
    total_subqueries = len(background_queries) * len(geometry_wkts)
    rows: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    truncated = False
    error_text = ""

    try:
        stop_after_detection = bool(search.get("stop_after_first_confirmatory_target_detection", True))
        for _taxon_name, taxon_key in background_queries:
            for geometry_wkt in geometry_wkts:
                offset = 0
                component_complete = False
                while raw_examined < max_records:
                    remaining = max_records - raw_examined
                    limit = min(page_size, remaining)
                    payload = _get_json_with_retry(
                        client,
                        endpoint,
                        params={
                            "taxon_key": taxon_key,
                            "geometry": geometry_wkt,
                            "has_coordinate": "true",
                            "limit": limit,
                            "offset": offset,
                        },
                        retries=retries,
                        backoff_seconds=backoff,
                    )
                    results = payload.get("results", [])
                    if not isinstance(results, list):
                        raise ValueError("GBIF Search results is not a list")
                    raw_examined += len(results)

                    for record in results:
                        row = _record_row(record, island_id)
                        if row is None:
                            continue
                        gbif_id = str(row["gbif_id"])
                        dedupe_key = gbif_id if gbif_id else json.dumps(row, sort_keys=True, default=str)
                        if dedupe_key in seen_ids:
                            continue
                        seen_ids.add(dedupe_key)
                        rows.append(row)

                    if stop_after_detection and _quality_target_detected(results, targets):
                        table = pd.DataFrame(rows)
                        return table, SearchMeta(
                            complete=False,
                            truncated=False,
                            error="",
                            raw_records_examined=raw_examined,
                            subqueries_completed=subqueries_completed,
                            subqueries_total=total_subqueries,
                            geometry_components=len(geometry_wkts),
                        )

                    end = bool(payload.get("endOfRecords")) or len(results) < limit
                    if end:
                        component_complete = True
                        subqueries_completed += 1
                        break
                    if not results:
                        raise ValueError("GBIF Search returned zero rows without endOfRecords")
                    offset += len(results)

                if not component_complete:
                    truncated = True
                    break
            if truncated:
                break
    except Exception as exc:  # noqa: BLE001
        error_text = f"{type(exc).__name__}: {exc}"

    complete = not truncated and not error_text and subqueries_completed == total_subqueries
    return pd.DataFrame(rows), SearchMeta(
        complete=complete,
        truncated=truncated,
        error=error_text,
        raw_records_examined=raw_examined,
        subqueries_completed=subqueries_completed,
        subqueries_total=total_subqueries,
        geometry_components=len(geometry_wkts),
    )


def _empty_observation_row(island_id: str, channel_id: str) -> dict[str, Any]:
    return {
        "island_id": island_id,
        "channel_id": channel_id,
        "observation_state": "insufficient_effort",
        "channel_record_count": 0,
        "background_record_count": 0,
        "background_spatial_units": 0,
        "background_temporal_units": 0,
        "distinct_dataset_count": 0,
        "latest_background_year": pd.NA,
        "evidence_source": "GBIF_exact_island_Search_API",
        "quality_flags": "below_background_records|low_spatial_dispersion|low_temporal_dispersion|low_dataset_dispersion|stale_or_missing_background_recency",
    }


def finalize_observation(
    records: pd.DataFrame,
    meta: SearchMeta,
    *,
    island_id: str,
    channel_id: str,
    catalog: pd.DataFrame,
    policy: dict[str, Any],
) -> dict[str, Any]:
    if records.empty:
        row = _empty_observation_row(island_id, channel_id)
    else:
        classified = classify_channel_observations(
            records,
            catalog,
            channel_id,
            policy,
            effort_tier="primary",
            confirmatory_catalog_only=True,
            evidence_source="GBIF_exact_island_Search_API",
        )
        match = classified.loc[classified["island_id"].astype(str).eq(str(island_id))]
        row = match.iloc[0].to_dict() if not match.empty else _empty_observation_row(island_id, channel_id)

    detected = str(row["observation_state"]) == "detected"
    flags = [value for value in str(row.get("quality_flags", "")).split("|") if value]
    if meta.error and not detected:
        row["observation_state"] = "unresolved"
        flags.append("search_error")
    elif meta.truncated and not detected:
        row["observation_state"] = "insufficient_effort"
        flags.append("search_truncated_fixed_budget")
    elif not meta.complete and not detected and not meta.error and not meta.truncated:
        # The only expected non-complete, non-error path is optional early stop after a
        # target detection. Without a detection, fail closed rather than infer absence.
        row["observation_state"] = "unresolved"
        flags.append("search_incomplete_without_target")
    elif str(row["observation_state"]) == "adequate_non_detection" and not meta.complete:
        row["observation_state"] = "insufficient_effort"
        flags.append("adequate_non_detection_blocked_without_exhaustive_search")

    row["quality_flags"] = "|".join(sorted(set(flags)))
    row["search_complete"] = bool(meta.complete)
    row["search_truncated"] = bool(meta.truncated)
    row["search_error"] = meta.error
    row["raw_records_examined"] = int(meta.raw_records_examined)
    row["background_subqueries_completed"] = int(meta.subqueries_completed)
    row["background_subqueries_total"] = int(meta.subqueries_total)
    row["geometry_components"] = int(meta.geometry_components)
    return row


def scan_channel(
    *,
    islands_gpkg: Path,
    source_availability: pd.DataFrame,
    catalog: pd.DataFrame,
    channel_id: str,
    config: dict[str, Any],
    policy: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    island_ids = eligible_islands(source_availability, channel_id)
    if not island_ids:
        return pd.DataFrame(columns=SEARCH_AUDIT_COLUMNS), {
            "contract": config["contract"],
            "channel_id": channel_id,
            "n_source_available_islands": 0,
            "n_detected": 0,
            "n_adequate_non_detection": 0,
            "n_insufficient_effort": 0,
            "n_unresolved": 0,
            "n_complete_searches": 0,
            "n_truncated_searches": 0,
            "n_search_errors": 0,
            "N1_fitted": False,
        }

    islands = gpd.read_file(islands_gpkg, layer="islands")
    if "island_id" not in islands.columns:
        raise ValueError("frozen island GeoPackage requires island_id")
    islands["island_id"] = islands["island_id"].astype(str)
    indexed = islands.set_index("island_id", drop=False)
    missing = sorted(set(island_ids).difference(indexed.index))
    if missing:
        raise ValueError(f"source-available islands missing from frozen geometry: {missing[:10]}")

    targets = target_species(catalog, channel_id, confirmatory_only=True)
    if not targets:
        raise ValueError(f"no confirmatory target taxa for channel {channel_id}")

    timeout = float(config["search"]["request_timeout_seconds"])
    output_rows: list[dict[str, Any]] = []
    with httpx.Client(
        timeout=timeout,
        follow_redirects=True,
        headers={"User-Agent": "island-v2/nee-island-search-observation"},
    ) as client:
        background_queries = resolve_background_queries(client, channel_id, config)
        for island_id in island_ids:
            try:
                geometry_wkts = exact_query_wkts(indexed.loc[island_id].geometry)
            except Exception as exc:  # noqa: BLE001
                output_rows.append(
                    finalize_observation(
                        pd.DataFrame(),
                        SearchMeta(False, False, f"geometry_error={type(exc).__name__}: {exc}", 0, 0, len(background_queries), 0),
                        island_id=island_id,
                        channel_id=channel_id,
                        catalog=catalog,
                        policy=policy,
                    )
                )
                continue
            records, meta = scan_island(
                client,
                island_id=island_id,
                geometry_wkts=geometry_wkts,
                background_queries=background_queries,
                targets=targets,
                config=config,
            )
            output_rows.append(
                finalize_observation(
                    records,
                    meta,
                    island_id=island_id,
                    channel_id=channel_id,
                    catalog=catalog,
                    policy=policy,
                )
            )

    result = pd.DataFrame(output_rows, columns=SEARCH_AUDIT_COLUMNS)
    counts = result["observation_state"].value_counts().to_dict()
    receipt = {
        "contract": config["contract"],
        "channel_id": channel_id,
        "n_source_available_islands": int(len(result)),
        "n_detected": int(counts.get("detected", 0)),
        "n_adequate_non_detection": int(counts.get("adequate_non_detection", 0)),
        "n_insufficient_effort": int(counts.get("insufficient_effort", 0)),
        "n_unresolved": int(counts.get("unresolved", 0)),
        "n_complete_searches": int(result["search_complete"].astype(bool).sum()),
        "n_truncated_searches": int(result["search_truncated"].astype(bool).sum()),
        "n_search_errors": int(result["search_error"].astype(str).str.len().gt(0).sum()),
        "max_records_per_island_channel": int(config["search"]["max_records_examined_per_island_channel"]),
        "adequate_non_detection_requires_exhaustive_search": True,
        "uses_focal_plant_traits": False,
        "uses_N1_effects": False,
        "N1_fitted": False,
    }
    return result, receipt


@app.command("scan")
def scan_command(
    islands_gpkg: Path = typer.Option(..., exists=True),
    source_availability_csv: Path = typer.Option(..., exists=True),
    catalog_csv: Path = typer.Option(..., exists=True),
    channel_id: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/chapter1_nee_island_search_acquisition.yml")),
    policy_path: Path = typer.Option(Path("config/chapter1_nee_channel_observation_policy.yml")),
) -> None:
    config = load_config(config_path)
    policy = load_policy(policy_path)
    source = pd.read_csv(source_availability_csv, dtype=str).fillna("")
    catalog = pd.read_csv(catalog_csv, dtype=str).fillna("")
    result, receipt = scan_channel(
        islands_gpkg=islands_gpkg,
        source_availability=source,
        catalog=catalog,
        channel_id=channel_id,
        config=config,
        policy=policy,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_dir / f"{channel_id}_island_search_observation.csv", index=False)
    (output_dir / f"{channel_id}_island_search_observation_receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(receipt))


if __name__ == "__main__":
    app()
