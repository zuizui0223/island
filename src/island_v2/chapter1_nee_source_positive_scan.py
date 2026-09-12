"""Positive-only pollinator-channel evidence for frozen GIFT mainland source entities.

This scanner can establish source ``available`` from a confirmatory functional-channel
species occurrence inside a frozen GIFT 3.2 polygon.  It can never establish absence:
a capped search with no hit is always ``unresolved``.  The cap therefore controls API
cost without creating false biological negatives.
"""

from __future__ import annotations

import json
import re
import time
from pathlib import Path
from typing import Any

import httpx
import pandas as pd
import typer
import yaml
from shapely.geometry import shape

from island_v2.gbif_blocks import canonical_gbif_polygon, resolve_query_taxon, stable_wkt

app = typer.Typer(add_completion=False, no_args_is_help=True)

GBIF_OCCURRENCE_SEARCH = "https://api.gbif.org/v1/occurrence/search"

OUTPUT_COLUMNS = [
    "entity_ID",
    "channel_id",
    "source_state",
    "evidence_id",
    "evidence_type",
    "source_citation",
    "source_url",
    "review_status",
    "quality_flags",
    "n_records_examined",
]


def load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("source positive-scan config must be a mapping")
    if payload.get("contract") != "chapter1_nee_source_positive_scan_v1":
        raise typer.BadParameter("unexpected source positive-scan contract")
    return payload


def _text(value: object) -> str:
    return "" if value is None else str(value).strip()


def normalized_binomial_key(value: object) -> str:
    """Return a conservative case-folded genus + species key or blank."""
    text = re.sub(r"\s+", " ", _text(value))
    parts = text.split(" ")
    if len(parts) < 2:
        return ""
    genus, epithet = parts[0], parts[1]
    if not genus or not epithet or epithet.lower() in {"sp.", "sp", "spp.", "spp"}:
        return ""
    return f"{genus.casefold()} {epithet.casefold()}"


def confirmatory_target_keys(catalog: pd.DataFrame, channel_id: str) -> set[str]:
    required = {"channel_id", "pollinator_species", "catalog_tier"}
    missing = required.difference(catalog.columns)
    if missing:
        raise ValueError(f"functional catalog missing columns: {sorted(missing)}")
    rows = catalog.loc[
        catalog["channel_id"].astype(str).eq(channel_id)
        & catalog["catalog_tier"].astype(str).eq("confirmatory")
    ]
    keys = {normalized_binomial_key(value) for value in rows["pollinator_species"]}
    keys.discard("")
    if not keys:
        raise ValueError(f"no confirmatory target species for channel {channel_id}")
    return keys


def _introduced(value: object) -> bool:
    text = _text(value).upper()
    return any(token in text for token in ("INTRODUCED", "INVASIVE", "ALIEN", "NATURALISED", "NATURALIZED"))


def _captive(record: dict[str, Any]) -> bool:
    establishment = _text(record.get("establishmentMeans")).upper()
    if any(token in establishment for token in ("CAPTIVE", "CULTIVATED")):
        return True
    return _text(record.get("basisOfRecord")).upper() == "LIVING_SPECIMEN"


def _explicit_absent(record: dict[str, Any]) -> bool:
    return _text(record.get("occurrenceStatus")).upper() == "ABSENT"


def find_confirmatory_positive(
    records: list[dict[str, Any]],
    target_keys: set[str],
    *,
    allowed_families: set[str] | None = None,
) -> tuple[dict[str, Any] | None, int, list[str]]:
    """Return the first quality-filtered confirmatory target hit and audit flags."""
    flags: list[str] = []
    examined = 0
    for record in records:
        examined += 1
        species_key = normalized_binomial_key(record.get("species"))
        if not species_key or species_key not in target_keys:
            continue
        if allowed_families is not None:
            family = _text(record.get("family"))
            if family not in allowed_families:
                continue
        if _explicit_absent(record) or _captive(record) or _introduced(record.get("establishmentMeans")):
            continue
        if not _text(record.get("establishmentMeans")):
            flags.append("matched_target_unknown_establishment")
        return record, examined, flags
    return None, examined, flags


def gift_entity_wkt(client: httpx.Client, entity_id: str, config: dict[str, Any]) -> str:
    url = str(config["source_geometry"]["geojson_template"]).format(entity_id=entity_id)
    response = client.get(url)
    response.raise_for_status()
    payload = response.json()
    features = payload.get("features", []) if isinstance(payload, dict) else []
    if not features:
        raise ValueError(f"GIFT entity {entity_id} has no geometry feature")
    geometry = canonical_gbif_polygon(shape(features[0]["geometry"]))
    return stable_wkt(geometry)


def _get_json_with_retry(
    client: httpx.Client,
    url: str,
    *,
    params: dict[str, Any],
    retries: int,
    backoff_seconds: float,
) -> dict[str, Any]:
    error: Exception | None = None
    for attempt in range(retries):
        try:
            response = client.get(url, params=params)
            response.raise_for_status()
            payload = response.json()
            if not isinstance(payload, dict):
                raise ValueError("GBIF search response is not an object")
            return payload
        except Exception as exc:  # network boundary; fail closed to unresolved upstream
            error = exc
            if attempt + 1 < retries:
                time.sleep(backoff_seconds * (attempt + 1))
    raise RuntimeError(f"GBIF source search failed after {retries} attempts: {error}")


def scan_entity_channel(
    client: httpx.Client,
    *,
    entity_id: str,
    geometry_wkt: str,
    channel_id: str,
    taxon_key: str,
    target_keys: set[str],
    config: dict[str, Any],
) -> dict[str, Any]:
    """Scan a capped occurrence window; no-hit and errors remain unresolved."""
    search = config["GBIF_search"]
    channel = config["channels"][channel_id]
    allowed_families = (
        {str(value) for value in channel.get("analytical_family_filter", [])}
        if channel.get("analytical_family_filter")
        else None
    )
    page_size = int(search["page_size"])
    max_pages = int(search["max_pages_per_entity_channel"])
    examined_total = 0
    flags: list[str] = []
    try:
        for page in range(max_pages):
            payload = _get_json_with_retry(
                client,
                str(search["endpoint"]),
                params={
                    "taxon_key": taxon_key,
                    "geometry": geometry_wkt,
                    "has_coordinate": "true",
                    "limit": page_size,
                    "offset": page * page_size,
                },
                retries=int(search["retries"]),
                backoff_seconds=float(search["retry_backoff_seconds"]),
            )
            results = payload.get("results", [])
            if not isinstance(results, list):
                raise ValueError("GBIF results is not a list")
            hit, examined, hit_flags = find_confirmatory_positive(
                results, target_keys, allowed_families=allowed_families
            )
            examined_total += examined
            flags.extend(hit_flags)
            if hit is not None:
                gbif_id = _text(hit.get("key") or hit.get("gbifID"))
                species = _text(hit.get("species"))
                return {
                    "entity_ID": entity_id,
                    "channel_id": channel_id,
                    "source_state": "available",
                    "evidence_id": f"GBIF:{gbif_id}" if gbif_id else f"GBIF_species:{species}",
                    "evidence_type": config["positive_evidence"]["evidence_type"],
                    "source_citation": f"GBIF occurrence {gbif_id or species} within frozen GIFT 3.2 entity {entity_id}",
                    "source_url": f"https://www.gbif.org/occurrence/{gbif_id}" if gbif_id else "https://www.gbif.org/",
                    "review_status": "accepted",
                    "quality_flags": "|".join(sorted(set(flags))),
                    "n_records_examined": examined_total,
                }
            if bool(payload.get("endOfRecords")) or len(results) < page_size:
                break
        return {
            "entity_ID": entity_id,
            "channel_id": channel_id,
            "source_state": "unresolved",
            "evidence_id": "",
            "evidence_type": "",
            "source_citation": "",
            "source_url": "",
            "review_status": "pending",
            "quality_flags": "capped_positive_scan_no_confirmatory_hit",
            "n_records_examined": examined_total,
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "entity_ID": entity_id,
            "channel_id": channel_id,
            "source_state": "unresolved",
            "evidence_id": "",
            "evidence_type": "",
            "source_citation": "",
            "source_url": "",
            "review_status": "pending",
            "quality_flags": f"source_scan_error={type(exc).__name__}",
            "n_records_examined": examined_total,
        }


def scan_channel(
    entity_ids: list[str],
    catalog: pd.DataFrame,
    channel_id: str,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if channel_id not in config["channels"]:
        raise ValueError(f"unregistered source-scan channel: {channel_id}")
    targets = confirmatory_target_keys(catalog, channel_id)
    channel = config["channels"][channel_id]
    timeout = float(config["GBIF_search"]["request_timeout_seconds"])
    rows: list[dict[str, Any]] = []
    geometry_failures = 0
    with httpx.Client(timeout=timeout, follow_redirects=True, headers={"User-Agent": "island-v2/nee-source-positive-scan"}) as client:
        resolved = resolve_query_taxon(
            client,
            str(channel["gbif_acquisition_taxon"]),
            str(channel["gbif_rank"]),
            None,
            None,
        )
        taxon_key = str(resolved["usage_key"])
        for entity_id in sorted({str(value) for value in entity_ids if str(value).strip()}):
            try:
                geometry_wkt = gift_entity_wkt(client, entity_id, config)
            except Exception as exc:  # noqa: BLE001
                geometry_failures += 1
                rows.append(
                    {
                        "entity_ID": entity_id,
                        "channel_id": channel_id,
                        "source_state": "unresolved",
                        "evidence_id": "",
                        "evidence_type": "",
                        "source_citation": "",
                        "source_url": "",
                        "review_status": "pending",
                        "quality_flags": f"gift_geometry_error={type(exc).__name__}",
                        "n_records_examined": 0,
                    }
                )
                continue
            rows.append(
                scan_entity_channel(
                    client,
                    entity_id=entity_id,
                    geometry_wkt=geometry_wkt,
                    channel_id=channel_id,
                    taxon_key=taxon_key,
                    target_keys=targets,
                    config=config,
                )
            )
    result = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    counts = result["source_state"].value_counts().to_dict() if not result.empty else {}
    receipt = {
        "contract": config["contract"],
        "channel_id": channel_id,
        "n_source_entities": int(len(result)),
        "n_available": int(counts.get("available", 0)),
        "n_unresolved": int(counts.get("unresolved", 0)),
        "n_structurally_absent": 0,
        "n_geometry_failures": int(geometry_failures),
        "max_records_per_entity_channel": int(config["GBIF_search"]["max_records_examined_per_entity_channel"]),
        "no_hit_means": "unresolved",
        "uses_focal_plant_traits": False,
        "uses_island_channel_outcomes": False,
    }
    return result, receipt


@app.command("scan")
def scan_command(
    source_assignments_csv: Path = typer.Option(..., exists=True),
    catalog_csv: Path = typer.Option(..., exists=True),
    channel_id: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/chapter1_nee_source_positive_scan.yml")),
) -> None:
    config = load_config(config_path)
    assignments = pd.read_csv(source_assignments_csv, dtype=str).fillna("")
    if "entity_ID" not in assignments.columns:
        raise typer.BadParameter("source assignments require entity_ID")
    catalog = pd.read_csv(catalog_csv, dtype=str).fillna("")
    result, receipt = scan_channel(assignments["entity_ID"].tolist(), catalog, channel_id, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_dir / f"{channel_id}_source_entity_states.csv", index=False)
    (output_dir / f"{channel_id}_source_positive_scan_receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(receipt))


if __name__ == "__main__":
    app()
