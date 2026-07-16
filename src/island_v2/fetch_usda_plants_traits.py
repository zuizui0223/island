"""Acquire direct flower-colour records from the public USDA PLANTS API.

The API is not versioned, so every response is snapshotted with retrieval time
and checksums.  Only direct ``Flower Color`` filter values for source records
explicitly classified as dicots or monocots are mapped.  Non-angiosperms,
infraspecific names, and unresolved identifiers remain in an exclusion audit.
"""

from __future__ import annotations

import hashlib
import json
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

API_ROOT = "https://plantsservices.sc.egov.usda.gov/api"
RESULTS_URL = f"{API_ROOT}/characteristicSearchResults"
FILTERS_URL = f"{API_ROOT}/characteristicSearchFilterResults"
PROFILE_ROOT = f"{API_ROOT}/PlantProfile"
DATABASE_URL = "https://plants.usda.gov"
HELP_URL = "https://plants.sc.egov.usda.gov/help"
USER_AGENT = "island-floral-traits/0.1 (source-backed USDA PLANTS acquisition)"
JsonGetter = Callable[[str], Any]

COLOR_MAP = {
    "Blue": "blue_purple",
    "Brown": "green_brown_inconspicuous",
    "Green": "green_brown_inconspicuous",
    "Orange": "yellow_orange",
    "Purple": "blue_purple",
    "Red": "red_pink",
    "White": "white",
    "Yellow": "yellow_orange",
}
ANGIOSPERM_GROUPS = {"Dicot", "Monocot"}


def _fetch_json(url: str, *, timeout_seconds: float = 180, retries: int = 4) -> Any:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
                return json.loads(response.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            if attempt >= retries:
                raise
            time.sleep(min(30.0, 2.0**attempt))
    raise RuntimeError("unreachable")


def _rows(value: Any, label: str) -> list[dict[str, Any]]:
    if isinstance(value, list) and all(isinstance(row, dict) for row in value):
        return value
    raise ValueError(f"USDA PLANTS {label} response is not a list of records")


def _text(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return " ".join(str(value).split())


def _integer(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def map_usda_flower_colors(raw_values: Iterable[str]) -> str | None:
    """Map the complete direct USDA value set without guessing unknown labels."""
    values = {_text(value) for value in raw_values if _text(value)}
    if not values or not values <= set(COLOR_MAP):
        return None
    mapped = {COLOR_MAP[value] for value in values}
    return next(iter(mapped)) if len(mapped) == 1 else "multicolored_variable"


def _filter_by_name(filters: list[dict[str, Any]], name: str) -> dict[str, Any]:
    matches = [row for row in filters if _text(row.get("name")) == name]
    if len(matches) != 1:
        raise ValueError(f"USDA PLANTS expected one {name!r} filter; found {len(matches)}")
    return matches[0]


def _value_members(filter_row: dict[str, Any]) -> dict[int, set[str]]:
    members: dict[int, set[str]] = defaultdict(set)
    for option in _rows(filter_row.get("filters"), f"{filter_row.get('name')} options"):
        value = _text(option.get("display"))
        for raw_id in option.get("plantIds") or []:
            plant_id = _integer(raw_id)
            if plant_id is not None:
                members[plant_id].add(value)
    return members


def _citation(accessed_on: str) -> str:
    return (
        "Natural Resources Conservation Service. PLANTS Database. "
        "United States Department of Agriculture. "
        f"Accessed {accessed_on} from {DATABASE_URL}."
    )


def _audit_row(
    *,
    source_plant_ids: Iterable[int],
    resolved_plant_id: int | None,
    scientific_name: str,
    raw_values: Iterable[str],
    groups: Iterable[str],
    reason: str,
) -> dict[str, str]:
    return {
        "source_plant_ids": "|".join(str(value) for value in sorted(source_plant_ids)),
        "resolved_plant_id": "" if resolved_plant_id is None else str(resolved_plant_id),
        "scientific_name": scientific_name,
        "raw_flower_colors": "|".join(sorted(raw_values)),
        "source_groups": "|".join(sorted(groups)),
        "reason": reason,
    }


def prepare_usda_rows(
    plant_rows: Iterable[dict[str, Any]],
    filter_rows: Iterable[dict[str, Any]],
    *,
    profile_resolutions: Iterable[dict[str, Any]] = (),
    accessed_on: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Validate USDA bulk responses and create canonical rows plus exclusions."""
    plants = list(plant_rows)
    filters = list(filter_rows)
    plants_by_id = {
        plant_id: row for row in plants if (plant_id := _integer(row.get("id"))) is not None
    }
    profiles_by_requested_id = {
        requested_id: row
        for row in profile_resolutions
        if (requested_id := _integer(row.get("requested_id"))) is not None
    }
    flower_by_id = _value_members(_filter_by_name(filters, "Flower Color"))
    group_by_id = _value_members(_filter_by_name(filters, "Group"))

    resolved: dict[int, dict[str, Any]] = defaultdict(
        lambda: {"source_ids": set(), "raw_values": set(), "groups": set(), "profile_urls": set()}
    )
    audit: list[dict[str, str]] = []
    for source_id, raw_values in sorted(flower_by_id.items()):
        groups = group_by_id.get(source_id, set())
        if len(groups) != 1:
            audit.append(
                _audit_row(
                    source_plant_ids=[source_id],
                    resolved_plant_id=None,
                    scientific_name="",
                    raw_values=raw_values,
                    groups=groups,
                    reason="missing_or_ambiguous_source_group",
                )
            )
            continue
        if not groups <= ANGIOSPERM_GROUPS:
            audit.append(
                _audit_row(
                    source_plant_ids=[source_id],
                    resolved_plant_id=None,
                    scientific_name="",
                    raw_values=raw_values,
                    groups=groups,
                    reason="non_angiosperm_source_group",
                )
            )
            continue

        resolved_id = source_id if source_id in plants_by_id else None
        profile_url = ""
        if resolved_id is None:
            resolution = profiles_by_requested_id.get(source_id, {})
            profile_url = _text(resolution.get("url"))
            response = resolution.get("response")
            if isinstance(response, dict):
                resolved_id = _integer(response.get("Id"))
        if resolved_id is None or resolved_id not in plants_by_id:
            audit.append(
                _audit_row(
                    source_plant_ids=[source_id],
                    resolved_plant_id=resolved_id,
                    scientific_name="",
                    raw_values=raw_values,
                    groups=groups,
                    reason="missing_taxon_resolution",
                )
            )
            continue

        item = resolved[resolved_id]
        item["source_ids"].add(source_id)
        item["raw_values"].update(raw_values)
        item["groups"].update(groups)
        if profile_url:
            item["profile_urls"].add(profile_url)

    prepared: list[dict[str, str]] = []
    for resolved_id, source in sorted(resolved.items()):
        plant = plants_by_id[resolved_id]
        scientific_name = _text(plant.get("scientificNameWithoutAuthor"))
        if len(scientific_name.split()) != 2 or "×" in scientific_name:
            audit.append(
                _audit_row(
                    source_plant_ids=source["source_ids"],
                    resolved_plant_id=resolved_id,
                    scientific_name=scientific_name,
                    raw_values=source["raw_values"],
                    groups=source["groups"],
                    reason="non_binomial_or_hybrid_taxon",
                )
            )
            continue
        value = map_usda_flower_colors(source["raw_values"])
        if value is None:
            audit.append(
                _audit_row(
                    source_plant_ids=source["source_ids"],
                    resolved_plant_id=resolved_id,
                    scientific_name=scientific_name,
                    raw_values=source["raw_values"],
                    groups=source["groups"],
                    reason="unmapped_source_trait_value",
                )
            )
            continue

        symbol = _text(plant.get("symbol"))
        evidence = {
            "characteristic_name": "Flower Color",
            "raw_values": sorted(source["raw_values"]),
            "source_plant_ids": sorted(source["source_ids"]),
            "resolved_plant_id": resolved_id,
            "symbol": symbol,
            "scientificNameWithoutAuthor": scientific_name,
            "source_groups": sorted(source["groups"]),
            "filter_url": FILTERS_URL,
            "profile_resolution_urls": sorted(source["profile_urls"]),
        }
        prepared.append(
            {
                "scientific_name": scientific_name,
                "trait_name": "flower_primary_color",
                "trait_value": value,
                "source_record_id": f"usda_plants:{accessed_on}:{resolved_id}:flower_color",
                "source_citation": _citation(accessed_on),
                "source_url": FILTERS_URL,
                "evidence_excerpt": json.dumps(evidence, ensure_ascii=False, separators=(",", ":")),
                "usda_resolved_plant_id": str(resolved_id),
                "usda_symbol": symbol,
                "usda_raw_values": "|".join(sorted(source["raw_values"])),
                "usda_profile_url": (
                    f"https://plants.sc.egov.usda.gov/plant-profile/{symbol}/characteristics"
                    if symbol
                    else ""
                ),
            }
        )

    prepared_frame = pd.DataFrame(
        prepared,
        columns=[
            "scientific_name",
            "trait_name",
            "trait_value",
            "source_record_id",
            "source_citation",
            "source_url",
            "evidence_excerpt",
            "usda_resolved_plant_id",
            "usda_symbol",
            "usda_raw_values",
            "usda_profile_url",
        ],
    )
    if not prepared_frame.empty:
        prepared_frame = prepared_frame.sort_values(
            ["scientific_name", "source_record_id"]
        ).reset_index(drop=True)
    audit_frame = pd.DataFrame(
        audit,
        columns=[
            "source_plant_ids",
            "resolved_plant_id",
            "scientific_name",
            "raw_flower_colors",
            "source_groups",
            "reason",
        ],
    )
    return prepared_frame, audit_frame


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def fetch_usda_plants_traits(
    *,
    output_dir: Path,
    max_workers: int = 4,
    getter: JsonGetter = _fetch_json,
    accessed_on: str | None = None,
) -> dict[str, Any]:
    """Snapshot the public USDA endpoints and prepare direct flower-colour rows."""
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir = output_dir / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    accessed_on = accessed_on or datetime.now(timezone.utc).date().isoformat()

    plants = _rows(getter(RESULTS_URL), "characteristic search results")
    filters = _rows(getter(FILTERS_URL), "characteristic search filters")
    flower_ids = set(_value_members(_filter_by_name(filters, "Flower Color")))
    known_ids = {plant_id for row in plants if (plant_id := _integer(row.get("id"))) is not None}
    missing_ids = sorted(flower_ids - known_ids)
    profile_urls = [f"{PROFILE_ROOT}/{plant_id}" for plant_id in missing_ids]
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        responses = list(executor.map(getter, profile_urls))
    profile_resolutions = [
        {"requested_id": plant_id, "url": url, "response": response}
        for plant_id, url, response in zip(missing_ids, profile_urls, responses, strict=True)
    ]

    prepared, audit = prepare_usda_rows(
        plants,
        filters,
        profile_resolutions=profile_resolutions,
        accessed_on=accessed_on,
    )
    _write_json(metadata_dir / "characteristicSearchResults.json", plants)
    _write_json(metadata_dir / "characteristicSearchFilterResults.json", filters)
    _write_json(metadata_dir / "missingPlantProfileResolutions.json", profile_resolutions)
    prepared_path = output_dir / "usda_plants_prepared_canonical_long.csv.gz"
    audit_path = output_dir / "usda_plants_unmapped_or_excluded_records.csv.gz"
    prepared.to_csv(prepared_path, index=False, compression="gzip")
    audit.to_csv(audit_path, index=False, compression="gzip")

    audit_reasons = Counter(audit["reason"]) if not audit.empty else Counter()
    manifest: dict[str, Any] = {
        "version": "1.0",
        "source_snapshot": "USDA PLANTS public current API (unversioned)",
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "accessed_on": accessed_on,
        "results_url": RESULTS_URL,
        "filters_url": FILTERS_URL,
        "profile_root": PROFILE_ROOT,
        "database_url": DATABASE_URL,
        "help_and_citation_url": HELP_URL,
        "n_characteristic_result_rows": len(plants),
        "n_filter_definitions": len(filters),
        "n_flower_color_source_ids": len(flower_ids),
        "n_missing_ids_resolved_by_profile": len(profile_resolutions),
        "n_prepared_canonical_rows": len(prepared),
        "n_prepared_species": int(prepared["scientific_name"].nunique())
        if not prepared.empty
        else 0,
        "n_unmapped_or_excluded_rows": len(audit),
        "audit_reason_counts": dict(sorted(audit_reasons.items())),
        "direct_value_map": COLOR_MAP,
        "rules": {
            "versioning": (
                "the public API exposes no release version; exact JSON snapshots, retrieval "
                "timestamp, query URLs, and file checksums are archived"
            ),
            "source_scope": "Flower Color records classified by USDA as Dicot or Monocot",
            "taxonomy": (
                "use characteristicSearchResults scientificNameWithoutAuthor; resolve only "
                "missing IDs through the exact PlantProfile/{id} endpoint; no infraspecific "
                "collapse or fuzzy match"
            ),
            "mapping": "only the eight enumerated USDA Flower Color values are mapped",
            "evidence": "evidence_excerpt serializes the exact source IDs, values, name, and group",
        },
        "files": {},
    }
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name != "usda_plants_fetch_manifest.json":
            manifest["files"][str(path.relative_to(output_dir))] = {
                "sha256": _sha256(path),
                "bytes": path.stat().st_size,
            }
    _write_json(output_dir / "usda_plants_fetch_manifest.json", manifest)
    return manifest


@app.callback()
def main() -> None:
    """Acquire public USDA PLANTS flower-colour records with complete audits."""


@app.command()
def run(
    output_dir: Path = typer.Option(...),
    max_workers: int = typer.Option(4, min=1, max=8),
) -> None:
    """Fetch, snapshot, validate, and prepare direct USDA PLANTS records."""
    typer.echo(
        json.dumps(
            fetch_usda_plants_traits(output_dir=output_dir, max_workers=max_workers),
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    app()
