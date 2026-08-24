"""Acquire exact self-fertility records from Ken Fern's public plant databases.

Useful Tropical Plants and Useful Temperate Plants expose a browsable species
index and a structured ``Self-fertile`` property.  This adapter consumes an
already-prioritised strict-unresolved queue, caches every source page, and
retains only exact accepted-name/family matches whose cultivation status
explicitly includes wild plants.  Positive self-fertility is SC evidence and
is never interpreted as autonomous selfing.

Both sites and Plants For A Future are assigned one Ken Fern species lineage
downstream, so overlapping republications cannot count as independent sources.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import html
import json
import re
import time
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd
import requests

USER_AGENT = (
    "island-trait-research/0.1 (scientific evidence audit; https://github.com/zuizui0223/island)"
)
H1_RE = re.compile(r'<div class="latin_name"><h1>(.*?)</h1>', re.IGNORECASE | re.DOTALL)
FAMILY_RE = re.compile(r'<div class="family"><h4>(.*?)</h4>', re.IGNORECASE | re.DOTALL)
ROW_RE = re.compile(
    r"<tr[^>]*>\s*<td[^>]*>(.*?)</td>\s*<td[^>]*>(.*?)</td>\s*</tr>", re.IGNORECASE | re.DOTALL
)
PROPERTY_TABLE_RE = re.compile(
    r'<table[^>]+class="[^"]*tableProperties[^"]*"[^>]*>(.*?)</table>',
    re.IGNORECASE | re.DOTALL,
)


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _visible(fragment: str) -> str:
    return " ".join(html.unescape(re.sub(r"<[^>]+>", " ", fragment)).split())


def parse_page(species: str, family: str, url: str, payload: bytes) -> dict[str, object]:
    raw = payload.decode("utf-8", errors="replace")
    h1_match = H1_RE.search(raw)
    family_match = FAMILY_RE.search(raw)
    matched_name = _visible(h1_match.group(1)) if h1_match else ""
    matched_family = _visible(family_match.group(1)) if family_match else ""
    table_match = PROPERTY_TABLE_RE.search(raw)
    table = table_match.group(1) if table_match else ""
    properties = {_visible(key): _visible(value) for key, value in ROW_RE.findall(table)}
    self_fertile = properties.get("Self-fertile", "")
    cultivation = properties.get("Cultivation Status", "")
    identity_exact = matched_name.casefold() == species.casefold()
    family_exact = matched_family.casefold() == family.casefold()
    wild_explicit = bool(re.search(r"\bwild\b", cultivation, re.IGNORECASE))
    if not identity_exact:
        status = "identity_mismatch"
    elif not family_exact:
        status = "family_mismatch"
    elif self_fertile.casefold() != "yes":
        status = "not_positive_self_fertile"
    elif not wild_explicit:
        status = "cultivated_or_ornamental_only"
    else:
        status = "accepted_wild_self_fertile"
    excerpt = (
        f"Self-fertile: {self_fertile}; Cultivation Status: {cultivation}" if self_fertile else ""
    )
    return {
        "accepted_species": species,
        "source_url": url,
        "matched_page_name": matched_name,
        "matched_family": matched_family,
        "identity_exact": identity_exact,
        "family_exact": family_exact,
        "self_fertile_raw": self_fertile,
        "cultivation_status": cultivation,
        "wild_explicit": wild_explicit,
        "supporting_excerpt": excerpt,
        "normalized_value": "SC" if status == "accepted_wild_self_fertile" else "",
        "trait_name": "self_incompatibility" if status == "accepted_wild_self_fertile" else "",
        "axis": "reproductive_assurance" if status == "accepted_wild_self_fertile" else "",
        "status": status,
        "page_sha256": _sha256(payload),
        "content_fingerprint": _sha256(excerpt.casefold().encode("utf-8")),
    }


def fetch_one(
    species: str, family: str, url: str, cache_dir: Path, timeout: float
) -> dict[str, object]:
    cache_path = cache_dir / f"{hashlib.sha256(url.encode('utf-8')).hexdigest()}.html"
    from_cache = cache_path.exists()
    try:
        if from_cache:
            payload = cache_path.read_bytes()
            http_status = 200
        else:
            response = requests.get(url, headers={"User-Agent": USER_AGENT}, timeout=timeout)
            http_status = response.status_code
            response.raise_for_status()
            payload = response.content
            cache_path.write_bytes(payload)
            time.sleep(0.15)
        row = parse_page(species, family, url, payload)
        row.update(
            {
                "http_status": http_status,
                "retrieved_at_utc": datetime.now(UTC).isoformat(),
                "from_cache": from_cache,
                "cache_file": cache_path.name,
                "error": "",
            }
        )
        return row
    except Exception as exc:  # noqa: BLE001 - every queued page remains accounted for
        return {
            "accepted_species": species,
            "source_url": url,
            "matched_page_name": "",
            "matched_family": "",
            "identity_exact": False,
            "family_exact": False,
            "self_fertile_raw": "",
            "cultivation_status": "",
            "wild_explicit": False,
            "supporting_excerpt": "",
            "normalized_value": "",
            "trait_name": "",
            "axis": "",
            "status": "fetch_error",
            "page_sha256": "",
            "content_fingerprint": "",
            "http_status": "",
            "retrieved_at_utc": datetime.now(UTC).isoformat(),
            "from_cache": from_cache,
            "cache_file": cache_path.name,
            "error": f"{type(exc).__name__}: {exc}",
        }


def build(
    *,
    queue_csv: Path,
    cache_dir: Path,
    output_dir: Path,
    max_workers: int,
    timeout: float,
) -> dict[str, object]:
    queue = pd.read_csv(queue_csv, dtype=str).fillna("")
    required = {"accepted_species", "family", "n_islands", "source_root", "source_url"}
    missing = required - set(queue.columns)
    if missing:
        raise ValueError(f"queue missing columns: {sorted(missing)}")
    queue["n_islands"] = pd.to_numeric(queue["n_islands"], errors="raise").astype(int)
    queue = queue.drop_duplicates(["accepted_species", "source_root"]).sort_values(
        ["n_islands", "accepted_species", "source_root"], ascending=[False, True, True]
    )
    cache_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    tasks = [(row.accepted_species, row.family, row.source_url) for row in queue.itertuples()]
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        rows = list(
            executor.map(
                lambda task: fetch_one(task[0], task[1], task[2], cache_dir, timeout), tasks
            )
        )
    # Preserve task order to avoid cross-pairing species present on both sites.
    inventory = queue.reset_index(drop=True).join(
        pd.DataFrame(rows).drop(columns=["accepted_species", "source_url"]).reset_index(drop=True)
    )
    accepted = inventory.loc[inventory["status"].eq("accepted_wild_self_fertile")].copy()
    # A Yes/No disagreement across the two Ken Fern databases is an internal
    # lineage conflict.  Do not select the positive row merely by site order.
    raw_values = inventory.groupby("accepted_species")["self_fertile_raw"].agg(
        lambda values: {value.casefold() for value in values if value}
    )
    conflicted = {species for species, values in raw_values.items() if len(values) > 1}
    inventory["lineage_internal_conflict"] = inventory["accepted_species"].isin(conflicted)
    accepted = accepted.loc[~accepted["accepted_species"].isin(conflicted)].copy()
    accepted = accepted.drop_duplicates("accepted_species")
    inventory.to_csv(output_dir / "ken_fern_fetch_inventory.csv.gz", index=False)
    accepted.to_csv(output_dir / "ken_fern_self_fertile_candidates.csv.gz", index=False)
    summary: dict[str, object] = {
        "contract": "ken_fern_exact_wild_species_self_fertility_acquisition_v1",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "queued_pages": len(queue),
        "queued_species": int(queue["accepted_species"].nunique()),
        "fetched_or_cached": int(inventory["http_status"].astype(str).eq("200").sum()),
        "fetch_errors": int(inventory["status"].eq("fetch_error").sum()),
        "identity_and_family_exact": int(
            (
                inventory["identity_exact"].astype(bool) & inventory["family_exact"].astype(bool)
            ).sum()
        ),
        "accepted_species": len(accepted),
        "accepted_n_islands_sum": int(accepted["n_islands"].sum()),
        "lineage_internal_conflict_species": len(conflicted),
        "status_counts": inventory["status"].value_counts().sort_index().to_dict(),
        "mapping": "wild species Self-fertile=Yes -> self_incompatibility=SC; never autonomous selfing",
        "source_lineage_policy": "PFAF, Useful Tropical Plants and Useful Temperate Plants collapse to one Ken Fern species lineage",
        "source_access": {
            "tropical": "https://tropical.theferns.info/",
            "temperate": "https://temperate.theferns.info/",
            "robots": "both robots.txt permit all paths",
            "license": "CC BY-NC-SA 3.0",
        },
    }
    (output_dir / "ken_fern_acquisition_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue-csv", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-workers", type=int, default=4)
    parser.add_argument("--timeout", type=float, default=45.0)
    args = parser.parse_args()
    print(json.dumps(build(**vars(args)), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
