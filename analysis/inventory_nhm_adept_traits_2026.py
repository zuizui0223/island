"""Inventory the 2026 NHM ADEPT angiosperm-trait resources.

The Natural History Museum exposes the data through a public CKAN datastore.
This inventory client keeps each immutable 1,000-row response as a gzip cache,
so interrupted runs resume without requesting completed pages again.  It does
not promote evidence: it only measures exact-name overlap and writes the raw
source records needed for a separate quote/value audit.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests

API_URL = "https://data.nhm.ac.uk/api/3/action/datastore_search"
PAGE_SIZE = 1_000
RESOURCES = {
    "bhl": "cc636680-abfb-4143-9b99-62dcb562f8da",
    "ecoflora": "378435a5-f6d4-42b6-be40-21d3bda73158",
    "flora_of_china": "6537b171-c40e-4936-90bd-f1cd4f9697e7",
    "flora_of_north_america": "dd2dccaf-8cd1-4091-940c-a505ae6cb03f",
    "flora_of_pakistan": "265a930d-a80f-4e46-b82d-0d36141d17ce",
}
IDENTITY_FIELDS = ["_id", "taxon", "source", "source_id", "wfo", "bhl"]
TRAIT_FIELDS = [
    "inflorescence arrangement",
    "flower architecture",
    "flower symmetry",
    "flower shape",
    "flower colour",
    "petal colour",
    "corolla colour",
    "corolla tube min_ length [cm]",
    "corolla tube max_ length [cm]",
    "labellum colour",
    "reproduction system",
    "pollination",
    "reward",
]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_cached(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def _fetch_page(
    session: requests.Session,
    *,
    resource_id: str,
    offset: int,
    after: list[Any] | None,
    cache_path: Path,
) -> dict[str, Any]:
    if cache_path.exists():
        return _read_cached(cache_path)
    params = {
        "resource_id": resource_id,
        "limit": PAGE_SIZE,
        "fields": ",".join(IDENTITY_FIELDS + TRAIT_FIELDS),
    }
    # NHM's CKAN deployment rejects deep numeric offsets.  Its stable cursor
    # includes the datastore version, so it also protects a resumed inventory
    # from silently crossing resource revisions.
    if after is None:
        params["offset"] = offset
    else:
        params["after"] = json.dumps(after)
    last_error: Exception | None = None
    for attempt in range(5):
        try:
            response = session.get(API_URL, params=params, timeout=90)
            response.raise_for_status()
            payload = response.json()
            if not payload.get("success"):
                raise RuntimeError(f"CKAN response was not successful: {payload}")
            result = payload["result"]
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(cache_path, "wt", encoding="utf-8", newline="") as handle:
                json.dump(result, handle, ensure_ascii=False, sort_keys=True)
            return result
        except (requests.RequestException, ValueError, RuntimeError) as exc:
            last_error = exc
            if attempt == 4:
                break
            time.sleep(2**attempt)
    raise RuntimeError(f"failed to fetch resource page at offset {offset}") from last_error


def inventory_resource(
    *,
    resource_name: str,
    master_csv: Path,
    output_dir: Path,
) -> dict[str, Any]:
    resource_id = RESOURCES[resource_name]
    cache_dir = output_dir / "page_cache" / resource_name
    session = requests.Session()
    session.headers.update(
        {"User-Agent": "island-trait-research/1.0 (reproducible CKAN client)"}
    )

    records: list[dict[str, Any]] = []
    offset = 0
    total: int | None = None
    after: list[Any] | None = None
    while total is None or offset < total:
        cache_path = cache_dir / f"offset-{offset:06d}.json.gz"
        result = _fetch_page(
            session,
            resource_id=resource_id,
            offset=offset,
            after=after,
            cache_path=cache_path,
        )
        total = int(result["total"])
        page = result.get("records", [])
        if not page and offset < total:
            raise RuntimeError(f"empty page before total at offset {offset}/{total}")
        records.extend(page)
        after = result.get("after")
        offset += len(page)

    frame = pd.DataFrame(records).fillna("")
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    accepted = set(master["accepted_species"])
    frame["accepted_name_exact"] = frame["taxon"].isin(accepted)
    frame["has_target_trait"] = frame.reindex(columns=TRAIT_FIELDS).ne("").any(axis=1)
    selected = frame.loc[frame["accepted_name_exact"] & frame["has_target_trait"]].copy()
    selected = selected.sort_values(["taxon", "source_id", "_id"], kind="stable")

    output_dir.mkdir(parents=True, exist_ok=True)
    records_path = output_dir / f"{resource_name}_exact_master_records.csv.gz"
    selected.to_csv(records_path, index=False, compression="gzip")

    counts = {
        field: int(selected[field].ne("").sum())
        for field in TRAIT_FIELDS
        if field in selected
    }
    report = {
        "contract": "nhm_adept_2026_inventory_v1",
        "resource_name": resource_name,
        "resource_id": resource_id,
        "api_url": API_URL,
        "page_size": PAGE_SIZE,
        "source_rows": len(frame),
        "exact_master_rows_with_target_trait": len(selected),
        "exact_master_species_with_target_trait": int(selected["taxon"].nunique()),
        "target_field_nonblank_rows": counts,
        "requests_required": len(list(cache_dir.glob("offset-*.json.gz"))),
        "resumed_without_refetch": True,
        "promotion_performed": False,
        "output": {
            "path": str(records_path),
            "sha256": _sha256(records_path),
            "bytes": records_path.stat().st_size,
        },
    }
    report_path = output_dir / f"{resource_name}_inventory_summary.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resource", choices=sorted(RESOURCES), default="bhl")
    parser.add_argument(
        "--master-csv",
        type=Path,
        default=Path("data/v2/staging/gbif/collected/island_taxa.csv"),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    report = inventory_resource(
        resource_name=args.resource,
        master_csv=args.master_csv,
        output_dir=args.output_dir,
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
