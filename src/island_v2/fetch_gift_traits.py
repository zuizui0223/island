"""Fetch direct, public floral trait records from the versioned GIFT API.

Only raw (non-derived) records are mapped.  Every API response is archived and
records that fail taxonomic, reference, or ontology checks remain in an audit
table rather than being silently discarded or guessed.
"""

from __future__ import annotations

import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

VERSIONS_URL = "https://gift.uni-goettingen.de/api/index.php?query=versions"
API_ROOT = "https://gift.uni-goettingen.de/api/extended"
DEFAULT_VERSION = "3.2"
TARGET_TRAITS = {
    "3.21.1": "Flower_colour",
    "3.4.1": "Reproduction_sexual_1",
    "3.6.1": "Pollination_syndrome_1",
    "3.6.2": "Pollination_syndrome_2",
    "3.6.3": "Pollination_syndrome_3",
}
USER_AGENT = "island-floral-traits/0.1 (source-backed GIFT acquisition)"
JsonGetter = Callable[[str], Any]

COLOR_MAP = {
    "white": "white",
    "yellow": "yellow_orange",
    "orange": "yellow_orange",
    "blue": "blue_purple",
    "violet": "blue_purple",
    "pink": "red_pink",
    "red": "red_pink",
    "green": "green_brown_inconspicuous",
    "brown": "green_brown_inconspicuous",
    "black": "other_described",
    "grey": "other_described",
}
ANIMAL_TOKENS = {
    "ant",
    "bat",
    "bee",
    "beetle",
    "bird",
    "butterfly",
    "fly",
    "insect",
    "moth",
}
SPECIFIC_GUILDS = {
    "ant": "beetles_wasps_ants",
    "bat": "bats",
    "beetle": "beetles_wasps_ants",
    "bird": "birds",
    "butterfly": "butterflies",
    "fly": "flies",
    "moth": "moths",
}


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
    if value in (None, ""):
        return []
    if isinstance(value, list) and all(isinstance(row, dict) for row in value):
        return value
    raise ValueError(f"GIFT {label} response is not a list of records")


def _text(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return " ".join(str(value).split())


def _flag(value: Any) -> bool:
    return _text(value).casefold() in {"1", "true", "yes"}


def _tokens(value: str) -> set[str]:
    return {token.strip().casefold() for token in value.split("/") if token.strip()}


def _map_pollination_mode(tokens: set[str]) -> str | None:
    if not tokens or "other" in tokens:
        return None
    animals = tokens & ANIMAL_TOKENS
    abiotic = tokens & {"wind", "water"}
    unknown = tokens - ANIMAL_TOKENS - {"wind", "water"}
    if unknown:
        return None
    if animals:
        return "mixed" if abiotic else "biotic"
    if tokens == {"wind"}:
        return "abiotic_wind"
    return None


def _map_coarse_guild(tokens: set[str]) -> str | None:
    if tokens and tokens <= {"wind", "water"}:
        return "wind_water"
    if tokens == {"bird"}:
        return "birds"
    if tokens & {"bird", "insect"}:
        return "mixed_or_generalist"
    return None


def _map_detailed_guild(tokens: set[str]) -> str | None:
    if tokens and tokens <= {"wind", "water"}:
        return "wind_water"
    animals = tokens & ANIMAL_TOKENS
    if not animals:
        return None
    if "insect" in animals or "bee" in animals or "other" in tokens:
        return "mixed_or_generalist"
    guilds = {SPECIFIC_GUILDS[token] for token in animals if token in SPECIFIC_GUILDS}
    if len(guilds) == 1 and not (tokens & {"wind", "water"}):
        return next(iter(guilds))
    return "mixed_or_generalist"


def map_gift_trait_value(trait_id: str, raw_value: str) -> list[tuple[str, str]]:
    """Map only direct semantic correspondences into the repository ontology."""
    raw_value = _text(raw_value)
    tokens = _tokens(raw_value)
    if trait_id == "3.21.1":
        if not tokens or not tokens <= set(COLOR_MAP):
            return []
        categories = {COLOR_MAP[token] for token in tokens}
        value = next(iter(categories)) if len(categories) == 1 else "multicolored_variable"
        return [("flower_primary_color", value)]
    if trait_id == "3.4.1":
        value = {
            "bisexual": "hermaphroditic",
            "monoecious": "monoecious",
            "dioecious": "dioecious",
        }.get(raw_value.casefold())
        return [("sex_system", value)] if value else []
    if trait_id not in {"3.6.2", "3.6.3"}:
        return []

    mapped: list[tuple[str, str]] = []
    mode = _map_pollination_mode(tokens)
    if mode:
        mapped.append(("pollen_vector_mode", mode))
    guild = _map_coarse_guild(tokens) if trait_id == "3.6.2" else _map_detailed_guild(tokens)
    if guild:
        mapped.append(("pollination_functional_guild", guild))
    return mapped


def _record_audit(row: dict[str, Any], reason: str) -> dict[str, str]:
    return {
        "trait_derived_ID": _text(row.get("trait_derived_ID")),
        "expected_ref_ID": _text(row.get("_expected_ref_ID")),
        "returned_ref_ID": _text(row.get("ref_ID")),
        "expected_trait_ID": _text(row.get("_expected_trait_ID")),
        "returned_trait_ID": _text(row.get("trait_ID")),
        "trait_value": _text(row.get("trait_value")),
        "original_name": f"{_text(row.get('genus'))} {_text(row.get('species_epithet'))}".strip(),
        "work_ID": _text(row.get("work_ID")),
        "query_url": _text(row.get("_query_url")),
        "reason": reason,
    }


def prepare_gift_rows(
    raw_rows: Iterable[dict[str, Any]],
    *,
    references_by_id: dict[str, dict[str, Any]],
    species_by_id: dict[str, dict[str, Any]],
    version: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Validate raw records and produce canonical rows plus an exclusion audit."""
    prepared: list[dict[str, str]] = []
    audit: list[dict[str, str]] = []
    for row in raw_rows:
        expected_ref = _text(row.get("_expected_ref_ID"))
        expected_trait = _text(row.get("_expected_trait_ID"))
        returned_ref = _text(row.get("ref_ID"))
        returned_trait = _text(row.get("trait_ID"))
        if returned_ref != expected_ref or returned_trait != expected_trait:
            audit.append(_record_audit(row, "api_query_scope_mismatch"))
            continue
        if _flag(row.get("derived")) or _flag(row.get("bias_deriv")):
            audit.append(_record_audit(row, "derived_or_bias_derived"))
            continue
        reference = references_by_id.get(returned_ref)
        if not reference or _text(reference.get("restricted")) != "0":
            audit.append(_record_audit(row, "missing_or_restricted_reference"))
            continue
        citation = _text(reference.get("ref_long"))
        if not citation:
            audit.append(_record_audit(row, "missing_reference_citation"))
            continue
        if not (_flag(row.get("matched")) and _flag(row.get("resolved"))):
            audit.append(_record_audit(row, "unmatched_or_unresolved_taxon"))
            continue
        if any(_flag(row.get(column)) for column in ("cf_genus", "cf_species", "aff_species")):
            audit.append(_record_audit(row, "uncertain_source_name"))
            continue
        if _text(row.get("subtaxon")):
            audit.append(_record_audit(row, "infraspecific_source_scope"))
            continue
        work_id = _text(row.get("work_ID"))
        species = species_by_id.get(work_id)
        scientific_name = _text((species or {}).get("work_species"))
        if len(scientific_name.split()) != 2:
            audit.append(_record_audit(row, "missing_or_non_binomial_resolved_name"))
            continue
        trait_id = returned_trait
        raw_value = _text(row.get("trait_value"))
        mappings = map_gift_trait_value(trait_id, raw_value)
        if not mappings:
            audit.append(_record_audit(row, "unmapped_source_trait_value"))
            continue
        record_id = _text(row.get("trait_derived_ID"))
        if not record_id:
            audit.append(_record_audit(row, "missing_trait_record_id"))
            continue
        original_name = f"{_text(row.get('genus'))} {_text(row.get('species_epithet'))}".strip()
        evidence_fields = {
            "trait_derived_ID": record_id,
            "ref_ID": returned_ref,
            "trait_ID": trait_id,
            "trait_value": raw_value,
            "genus": _text(row.get("genus")),
            "species_epithet": _text(row.get("species_epithet")),
            "author": _text(row.get("author")),
            "matched": _text(row.get("matched")),
            "resolved": _text(row.get("resolved")),
            "work_ID": work_id,
            "work_species": scientific_name,
        }
        excerpt = json.dumps(evidence_fields, ensure_ascii=False, separators=(",", ":"))
        for trait_name, trait_value in mappings:
            prepared.append(
                {
                    "scientific_name": scientific_name,
                    "trait_name": trait_name,
                    "trait_value": trait_value,
                    "source_record_id": f"gift:{version}:{record_id}:{trait_name}",
                    "source_citation": citation,
                    "source_url": _text(row.get("_query_url")),
                    "evidence_excerpt": excerpt,
                    "gift_trait_ID": trait_id,
                    "gift_raw_value": raw_value,
                    "gift_original_name": original_name,
                    "gift_work_ID": work_id,
                    "gift_ref_ID": returned_ref,
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
            "gift_trait_ID",
            "gift_raw_value",
            "gift_original_name",
            "gift_work_ID",
            "gift_ref_ID",
        ],
    )
    if not prepared_frame.empty:
        prepared_frame = prepared_frame.drop_duplicates(
            ["scientific_name", "trait_name", "trait_value", "source_record_id"]
        ).sort_values(["scientific_name", "trait_name", "source_record_id"])
    audit_frame = pd.DataFrame(
        audit,
        columns=[
            "trait_derived_ID",
            "expected_ref_ID",
            "returned_ref_ID",
            "expected_trait_ID",
            "returned_trait_ID",
            "trait_value",
            "original_name",
            "work_ID",
            "query_url",
            "reason",
        ],
    )
    return prepared_frame.reset_index(drop=True), audit_frame.reset_index(drop=True)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _query_url(endpoint: str, **parameters: Any) -> str:
    return f"{endpoint}?{urllib.parse.urlencode(parameters)}"


def _target_pairs(reference_traits: list[dict[str, Any]]) -> list[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    for row in reference_traits:
        bias = _text(row.get("bias"))
        if bias not in {"", "0"}:
            continue
        ref_id = _text(row.get("ref_ID"))
        for column, value in row.items():
            trait_id = _text(value)
            if str(column).startswith("trait") and trait_id in TARGET_TRAITS and ref_id:
                pairs.add((ref_id, trait_id))
    return sorted(pairs, key=lambda pair: (int(pair[0]), pair[1]))


def fetch_gift_direct_traits(
    *,
    output_dir: Path,
    version: str = DEFAULT_VERSION,
    max_workers: int = 4,
    getter: JsonGetter = _fetch_json,
) -> dict[str, Any]:
    """Fetch a complete, version-pinned set of selected direct GIFT records."""
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir = output_dir / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)

    versions = _rows(getter(VERSIONS_URL), "versions")
    if version not in {_text(row.get("version")) for row in versions}:
        raise ValueError(f"GIFT version {version!r} is not a published stable version")
    endpoint = f"{API_ROOT}/index{version}.php"
    metadata_urls = {
        "traits_meta": _query_url(endpoint, query="traits_meta"),
        "reference_traits": _query_url(endpoint, query="reference_traits"),
        "references": _query_url(endpoint, query="references"),
    }
    metadata = {name: _rows(getter(url), name) for name, url in metadata_urls.items()}
    references_by_id = {
        _text(row.get("ref_ID")): row for row in metadata["references"] if row.get("ref_ID")
    }
    all_pairs = _target_pairs(metadata["reference_traits"])
    pairs = [
        pair
        for pair in all_pairs
        if _text((references_by_id.get(pair[0]) or {}).get("restricted")) == "0"
    ]
    if not pairs:
        raise ValueError("GIFT metadata produced no public target trait/reference pairs")

    species_urls = [
        _query_url(endpoint, query="species", startat=start, limit=100_000)
        for start in range(0, 600_000, 100_000)
    ]
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        species_responses = list(executor.map(getter, species_urls))
    species_rows = [
        row
        for index, response in enumerate(species_responses)
        for row in _rows(response, f"species chunk {index}")
    ]
    species_by_id = {_text(row.get("work_ID")): row for row in species_rows if row.get("work_ID")}

    pair_urls = [
        (
            pair,
            _query_url(
                endpoint,
                query="traits_raw",
                traitid=pair[1],
                deriv=0,
                biasderiv=0,
                refid=pair[0],
            ),
        )
        for pair in pairs
    ]
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        raw_responses = list(executor.map(getter, [url for _, url in pair_urls]))
    raw_rows: list[dict[str, Any]] = []
    response_counts: dict[str, int] = {}
    for ((ref_id, trait_id), url), response in zip(pair_urls, raw_responses, strict=True):
        rows = _rows(response, f"traits_raw {trait_id}/{ref_id}")
        response_counts[f"{trait_id}:{ref_id}"] = len(rows)
        for row in rows:
            raw_rows.append(
                {
                    **row,
                    "_expected_ref_ID": ref_id,
                    "_expected_trait_ID": trait_id,
                    "_query_url": url,
                }
            )

    prepared, audit = prepare_gift_rows(
        raw_rows,
        references_by_id=references_by_id,
        species_by_id=species_by_id,
        version=version,
    )

    _write_json(metadata_dir / "versions.json", versions)
    for name, rows in metadata.items():
        _write_json(metadata_dir / f"{name}.json", rows)
    pd.DataFrame(species_rows).to_csv(
        metadata_dir / "species.csv.gz", index=False, compression="gzip"
    )
    raw_path = output_dir / "gift_raw_api_records.csv.gz"
    prepared_path = output_dir / "gift_prepared_canonical_long.csv.gz"
    audit_path = output_dir / "gift_unmapped_or_excluded_records.csv.gz"
    pd.DataFrame(raw_rows).to_csv(raw_path, index=False, compression="gzip")
    prepared.to_csv(prepared_path, index=False, compression="gzip")
    audit.to_csv(audit_path, index=False, compression="gzip")

    audit_reasons = Counter(audit["reason"]) if not audit.empty else Counter()
    manifest = {
        "version": "1.0",
        "gift_version": version,
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "api_endpoint": endpoint,
        "versions_url": VERSIONS_URL,
        "metadata_urls": metadata_urls,
        "target_traits": TARGET_TRAITS,
        "n_reference_trait_pairs_metadata": len(all_pairs),
        "n_public_unbiased_reference_trait_pairs_fetched": len(pairs),
        "n_restricted_or_missing_reference_pairs_excluded": len(all_pairs) - len(pairs),
        "n_species_metadata_rows": len(species_rows),
        "n_raw_api_rows": len(raw_rows),
        "n_prepared_canonical_rows": len(prepared),
        "n_prepared_species": int(prepared["scientific_name"].nunique())
        if not prepared.empty
        else 0,
        "n_prepared_species_trait_pairs": int(
            len(prepared.drop_duplicates(["scientific_name", "trait_name"]))
        ),
        "prepared_species_by_trait": {
            str(trait): int(count)
            for trait, count in prepared.groupby("trait_name")["scientific_name"]
            .nunique()
            .sort_values(ascending=False)
            .items()
        }
        if not prepared.empty
        else {},
        "n_unmapped_or_excluded_rows": len(audit),
        "audit_reason_counts": dict(sorted(audit_reasons.items())),
        "raw_response_counts": response_counts,
        "rules": {
            "version_pin": "published GIFT stable version endpoint only",
            "source_scope": "public references with reference_traits bias=0 or missing",
            "derivation": "traits_raw deriv=0 and biasderiv=0; returned fields rechecked",
            "taxonomy": (
                "matched=1, resolved=1, no cf/aff flags or infraspecific scope; use GIFT "
                "version-pinned work_species"
            ),
            "mapping": (
                "only enumerated categorical correspondences; unsupported and ambiguous values "
                "remain in the audit"
            ),
            "evidence": "evidence_excerpt is a compact serialization of exact raw GIFT fields",
        },
        "files": {},
    }
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name != "gift_fetch_manifest.json":
            manifest["files"][str(path.relative_to(output_dir))] = {
                "sha256": _sha256(path),
                "bytes": path.stat().st_size,
            }
    manifest_path = output_dir / "gift_fetch_manifest.json"
    _write_json(manifest_path, manifest)
    return manifest


@app.callback()
def main() -> None:
    """Acquire direct public GIFT floral traits with complete exclusion audits."""


@app.command()
def run(
    output_dir: Path = typer.Option(...),
    version: str = typer.Option(DEFAULT_VERSION),
    max_workers: int = typer.Option(4, min=1, max=8),
) -> None:
    """Fetch, validate, map, and archive the selected GIFT records."""
    manifest = fetch_gift_direct_traits(
        output_dir=output_dir,
        version=version,
        max_workers=max_workers,
    )
    typer.echo(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    app()
