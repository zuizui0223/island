"""Resolve high-value coverage gaps with local, species-direct public-web retrieval."""

from __future__ import annotations

import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2.europe_pmc_discovery import (
    FULL_TEXT_ERROR_COLUMNS,
    FULL_TEXT_SOURCE_COLUMNS,
    SOURCE_COLUMNS as DISCOVERY_SOURCE_COLUMNS,
    discover_europe_pmc_full_text_sources,
    discover_europe_pmc_sources,
    extract_species_direct_text_passages,
)
from island_v2.v1_category_search import (
    SOURCE_COLUMNS,
    extract_evidence_from_sources,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT_VERSION = "direct_web_reacquisition_europe_pmc_full_text_v1"
TARGET_AXIS = "reproductive_assurance"
WFO_RELEASE = "2026-06"
WFO_RELEASE_DOI = "10.5281/zenodo.20782718"
WFO_MATCH_URL = "https://list.worldfloraonline.org/matching_rest.php"
WFO_GQL_URL = "https://list.worldfloraonline.org/gql.php"
USER_AGENT = "island-floral-traits/0.1 (local direct-web reacquisition)"
QUALITY_RANK = {"": 0, "low": 1, "medium": 2, "high": 3}
DIRECT_EVIDENCE_COLUMNS = [
    "accepted_species",
    "axis",
    "trait_name",
    "normalized_value",
    "evidence_quality",
    "evidence_scope",
    "source_group",
    "source_provider",
    "source_url",
    "source_citation",
    "source_excerpt",
    "source_record_id",
    "source_lineage",
    "name_resolution_lineage",
    "name_match_method",
    "inference_rule",
]


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _valid_species_name(value: object) -> bool:
    name = _text(value)
    return bool(
        re.fullmatch(
            r"[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+(?: [a-z][A-Za-z.-]+)?",
            name,
        )
    )


def _wfo_json_request(
    url: str,
    *,
    payload: dict[str, str] | None = None,
    timeout_seconds: float = 45,
    max_attempts: int = 3,
) -> dict[str, Any]:
    body = json.dumps(payload).encode("utf-8") if payload is not None else None
    headers = {"Accept": "application/json", "User-Agent": USER_AGENT}
    if body is not None:
        headers["Content-Type"] = "application/json"
    for attempt in range(1, max_attempts + 1):
        request = urllib.request.Request(url, data=body, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
                value = json.loads(response.read().decode("utf-8"))
                if not isinstance(value, dict):
                    raise ValueError("WFO response was not a JSON object")
                return value
        except urllib.error.HTTPError as exc:
            if exc.code not in {408, 425, 429, 500, 502, 503, 504} or attempt == max_attempts:
                raise
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            if attempt == max_attempts:
                raise
        time.sleep(min(30.0, 2.0 ** (attempt - 1)))
    raise RuntimeError("unreachable")


def _canonical_species_name(full_name: object) -> str:
    text = _text(full_name)
    match = re.match(
        r"^([A-Z][A-Za-z.-]+)\s+([a-z][A-Za-z.-]+)(?:\s+(.*))?$",
        text,
    )
    if not match:
        return ""
    remainder = _text(match.group(3))
    if remainder and re.match(
        r"^(?:subsp\.|ssp\.|var\.|f\.)(?:\s|$)",
        remainder,
        re.IGNORECASE,
    ):
        return ""
    return f"{match.group(1)} {match.group(2)}"


def _wfo_aliases_from_payloads(
    accepted_species: str,
    match_payload: dict[str, Any],
    graph_payload: dict[str, Any],
    *,
    max_aliases: int,
) -> list[dict[str, str]]:
    match = match_payload.get("match")
    parsed = match_payload.get("parsedName")
    if (
        not isinstance(match, dict)
        or not isinstance(parsed, dict)
        or _text(parsed.get("canonical_form")).casefold()
        != accepted_species.casefold()
    ):
        return []
    observed_wfo_id = _text(match.get("wfo_id"))
    taxon = (graph_payload.get("data") or {}).get("taxonNameById")
    if not isinstance(taxon, dict):
        return []
    usage = taxon.get("currentPreferredUsage")
    if not isinstance(usage, dict):
        return []
    accepted_wfo_id = _text((usage.get("hasName") or {}).get("id"))
    usage_id = _text(usage.get("id"))
    rows: list[dict[str, str]] = []
    for synonym in usage.get("hasSynonym") or []:
        if not isinstance(synonym, dict):
            continue
        search_name = _canonical_species_name(synonym.get("fullNameStringPlain"))
        synonym_id = _text(synonym.get("id"))
        if (
            not search_name
            or search_name.casefold() == accepted_species.casefold()
            or not synonym_id
        ):
            continue
        rows.append(
            {
                "accepted_species": accepted_species,
                "search_name": search_name,
                "name_role": "synonym",
                "name_resolution_lineage": (
                    f"wfo:{WFO_RELEASE_DOI}:{synonym_id}"
                ),
                "observed_wfo_id": observed_wfo_id,
                "accepted_wfo_id": accepted_wfo_id,
                "accepted_usage_id": usage_id,
                "synonym_wfo_id": synonym_id,
                "wfo_release": WFO_RELEASE,
                "wfo_release_doi": WFO_RELEASE_DOI,
            }
        )
    return list(
        {
            row["search_name"].casefold(): row
            for row in rows
        }.values()
    )[:max_aliases]


def resolve_wfo_exact_synonyms(
    species_names: list[str],
    *,
    min_interval_seconds: float = 1.0,
    max_aliases_per_species: int = 20,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Resolve exact WFO accepted usages and their species-rank synonyms."""

    aliases: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    for species in species_names:
        started = time.monotonic()
        try:
            parameters = urllib.parse.urlencode(
                {
                    "input_string": species,
                    "fuzzy_names": "0",
                    "fuzzy_authors": "0",
                    "check_homonyms": "true",
                    "check_rank": "true",
                }
            )
            match_payload = _wfo_json_request(f"{WFO_MATCH_URL}?{parameters}")
            match = match_payload.get("match")
            wfo_id = _text(match.get("wfo_id")) if isinstance(match, dict) else ""
            if not wfo_id:
                raise ValueError("WFO exact match returned no wfo_id")
            query = (
                "query { taxonNameById(nameId: "
                f'"{wfo_id}") {{ id stableUri fullNameStringPlain '
                "currentPreferredUsage { id hasName { id stableUri "
                "fullNameStringPlain } hasSynonym { id stableUri "
                "fullNameStringPlain } } } }"
            )
            elapsed = time.monotonic() - started
            if elapsed < min_interval_seconds:
                time.sleep(min_interval_seconds - elapsed)
            graph_payload = _wfo_json_request(
                WFO_GQL_URL,
                payload={"query": query},
            )
            aliases.extend(
                _wfo_aliases_from_payloads(
                    species,
                    match_payload,
                    graph_payload,
                    max_aliases=max_aliases_per_species,
                )
            )
        except Exception as exc:  # noqa: BLE001
            errors.append(
                {
                    "accepted_species": species,
                    "source": "wfo_plant_list",
                    "error": _text(exc),
                }
            )
        elapsed = time.monotonic() - started
        if elapsed < min_interval_seconds:
            time.sleep(min_interval_seconds - elapsed)
    return pd.DataFrame(aliases), pd.DataFrame(
        errors,
        columns=["accepted_species", "source", "error"],
    )


def select_reproductive_targets(
    coverage: pd.DataFrame,
    *,
    max_species: int,
) -> pd.DataFrame:
    """Select two-axis-complete species whose reproductive axis is unresolved."""

    required = {"accepted_species", "axis", "after_quality"}
    missing = required.difference(coverage.columns)
    if missing:
        raise ValueError(f"coverage is missing columns: {sorted(missing)}")
    if max_species < 1:
        raise ValueError("max_species must be positive")
    table = coverage[list(required)].copy().fillna("")
    table["accepted_species"] = table["accepted_species"].map(_text)
    table["axis"] = table["axis"].map(_text)
    table["after_quality"] = table["after_quality"].map(lambda value: _text(value).casefold())
    duplicate = table.duplicated(["accepted_species", "axis"], keep=False)
    if duplicate.any():
        raise ValueError("coverage contains duplicate accepted_species x axis rows")

    quality = table.pivot(
        index="accepted_species",
        columns="axis",
        values="after_quality",
    ).fillna("")
    for axis in (
        "flower_colour",
        "floral_structural_complexity",
        TARGET_AXIS,
    ):
        if axis not in quality:
            quality[axis] = ""
    selected = quality.loc[
        quality[TARGET_AXIS].eq("")
        & quality["flower_colour"].ne("")
        & quality["floral_structural_complexity"].ne("")
    ].copy()
    selected = selected.loc[selected.index.map(_valid_species_name)].copy()
    selected["existing_quality_score"] = (
        selected["flower_colour"].map(QUALITY_RANK)
        + selected["floral_structural_complexity"].map(QUALITY_RANK)
    )
    selected = (
        selected.reset_index()
        .sort_values(
            ["existing_quality_score", "accepted_species"],
            ascending=[False, True],
        )
        .head(max_species)
        .reset_index(drop=True)
    )
    selected.insert(1, "target_axis", TARGET_AXIS)
    return selected[
        [
            "accepted_species",
            "target_axis",
            "flower_colour",
            "floral_structural_complexity",
            "existing_quality_score",
        ]
    ]


def _shared_full_text_sources(full_text_sources: pd.DataFrame) -> pd.DataFrame:
    if full_text_sources.empty:
        return pd.DataFrame(columns=SOURCE_COLUMNS)
    rows: list[dict[str, str]] = []
    group_columns = [
        "accepted_species",
        "source_url",
        "source_citation",
        "source_type",
        "evidence_scope",
    ]
    for keys, group in full_text_sources.groupby(group_columns, sort=False):
        row = dict(zip(group_columns, keys, strict=True))
        row["source_text"] = " ".join(
            dict.fromkeys(group["source_text"].map(_text))
        )
        rows.append(row)
    return pd.DataFrame(rows, columns=SOURCE_COLUMNS)


def _search_name_table(
    targets: pd.DataFrame,
    aliases: pd.DataFrame | None,
) -> pd.DataFrame:
    rows = [
        {
            "accepted_species": species,
            "search_name": species,
            "name_role": "accepted",
            "name_resolution_lineage": "",
        }
        for species in targets["accepted_species"].map(_text)
    ]
    if aliases is not None and not aliases.empty:
        required = {"accepted_species", "search_name"}
        missing = required.difference(aliases.columns)
        if missing:
            raise ValueError(f"aliases are missing columns: {sorted(missing)}")
        targets_set = set(targets["accepted_species"].map(_text))
        for record in aliases.fillna("").to_dict("records"):
            accepted = _text(record.get("accepted_species"))
            search_name = _text(record.get("search_name"))
            if accepted not in targets_set or not _valid_species_name(search_name):
                continue
            rows.append(
                {
                    "accepted_species": accepted,
                    "search_name": search_name,
                    "name_role": _text(record.get("name_role")) or "synonym",
                    "name_resolution_lineage": _text(
                        record.get("name_resolution_lineage")
                    ),
                }
            )
    table = pd.DataFrame(rows).drop_duplicates(
        ["accepted_species", "search_name"],
        keep="first",
    )
    ambiguous = (
        table.groupby("search_name")["accepted_species"].nunique().loc[lambda value: value > 1]
    )
    if len(ambiguous):
        raise ValueError(
            "one search name maps to multiple accepted species: "
            f"{ambiguous.index.tolist()[:10]}"
        )
    return table.reset_index(drop=True)


def _remap_discovery_sources(
    records: tuple[dict[str, str], ...],
    search_names: pd.DataFrame,
) -> tuple[dict[str, str], ...]:
    metadata = search_names.set_index("search_name")
    remapped: list[dict[str, str]] = []
    for record in records:
        searched = _text(record.get("accepted_species"))
        if searched not in metadata.index:
            raise ValueError(f"discovery source has an unplanned search name: {searched}")
        mapping = metadata.loc[searched]
        row = dict(record)
        row["accepted_species"] = _text(mapping["accepted_species"])
        row["searched_name"] = searched
        row["name_resolution_lineage"] = _text(
            mapping["name_resolution_lineage"]
        )
        if searched.casefold() != row["accepted_species"].casefold():
            row["source_citation"] = " | ".join(
                part
                for part in (
                    _text(row.get("source_citation")),
                    f"searched_exact_synonym={searched}",
                    row["name_resolution_lineage"],
                )
                if part
            )
        remapped.append(row)
    return tuple(remapped)


def _title_abstract_direct_sources(
    records: tuple[dict[str, str], ...],
    aliases_by_species: dict[str, list[str]],
) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for record in records:
        species = _text(record.get("accepted_species"))
        doi = _text(record.get("doi")).casefold()
        pmcid = _text(record.get("pmcid")).upper()
        source_code = _text(record.get("europe_pmc_source"))
        source_id = _text(record.get("europe_pmc_id"))
        lineage = (
            f"doi:{doi}"
            if doi
            else (
                f"pmcid:{pmcid.casefold()}"
                if pmcid
                else f"europe-pmc:{source_code.casefold()}:{source_id.casefold()}"
            )
        )
        passages = extract_species_direct_text_passages(
            _text(record.get("source_text")),
            species,
            aliases=aliases_by_species.get(species, []),
        )
        for passage in passages:
            rows.append(
                {
                    "accepted_species": species,
                    "source_text": passage["text"],
                    "source_url": _text(record.get("source_url")),
                    "source_citation": _text(record.get("source_citation")),
                    "source_type": "europe_pmc_title_abstract_direct",
                    "evidence_scope": "species_direct",
                    "title": _text(record.get("title")),
                    "pmcid": pmcid,
                    "doi": doi,
                    "license": _text(record.get("license")),
                    "is_open_access": _text(record.get("is_open_access")),
                    "full_text_xml_url": "",
                    "full_text_xml_sha256": "",
                    "passage_sha256": passage["passage_sha256"],
                    "source_record_id": (
                        f"EPMC:{source_code}:{source_id}:"
                        f"{passage['passage_sha256']}"
                    ),
                    "source_lineage": lineage,
                    "searched_name": _text(record.get("searched_name")) or species,
                    "matched_name": passage["matched_name"],
                    "name_match_method": (
                        "accepted_exact"
                        if passage["matched_name"].casefold() == species.casefold()
                        else "synonym_exact"
                    ),
                    "name_resolution_lineage": _text(
                        record.get("name_resolution_lineage")
                    ),
                    "retrieved_at_utc": _text(record.get("retrieved_at_utc")),
                }
            )
    return pd.DataFrame(rows, columns=FULL_TEXT_SOURCE_COLUMNS)


def _direct_evidence(
    extracted: pd.DataFrame,
    full_text_sources: pd.DataFrame,
) -> pd.DataFrame:
    if extracted.empty:
        return pd.DataFrame(columns=DIRECT_EVIDENCE_COLUMNS)
    metadata = (
        full_text_sources.sort_values("source_record_id")
        .drop_duplicates(["accepted_species", "source_url"])
        .set_index(["accepted_species", "source_url"])
    )
    rows: list[dict[str, str]] = []
    for record in extracted.fillna("").to_dict("records"):
        if record["field"] not in {"mating_system", "self_incompatibility"}:
            continue
        key = (_text(record["species"]), _text(record["source_url"]))
        if key not in metadata.index:
            raise ValueError(f"full-text evidence has no source metadata: {key}")
        source = metadata.loc[key]
        rows.append(
            {
                "accepted_species": key[0],
                "axis": TARGET_AXIS,
                "trait_name": _text(record["field"]),
                "normalized_value": _text(record["value"]),
                "evidence_quality": "medium",
                "evidence_scope": "species_direct",
                "source_group": "direct_web_fulltext",
                "source_provider": "Europe PMC",
                "source_url": key[1],
                "source_citation": _text(record["source_citation"]),
                "source_excerpt": _text(record["source_excerpt"]),
                "source_record_id": (
                    f"{_text(source['source_record_id'])}:"
                    f"{_text(record['rule_id'])}"
                ),
                "source_lineage": _text(source["source_lineage"]),
                "name_resolution_lineage": _text(
                    source.get("name_resolution_lineage", "")
                ),
                "name_match_method": _text(source.get("name_match_method"))
                or "accepted_exact",
                "inference_rule": "",
            }
        )
    if not rows:
        return pd.DataFrame(columns=DIRECT_EVIDENCE_COLUMNS)
    return (
        pd.DataFrame(rows, columns=DIRECT_EVIDENCE_COLUMNS)
        .drop_duplicates(
            [
                "accepted_species",
                "axis",
                "trait_name",
                "normalized_value",
                "source_lineage",
            ]
        )
        .sort_values(
            [
                "accepted_species",
                "trait_name",
                "normalized_value",
                "source_lineage",
            ]
        )
        .reset_index(drop=True)
    )


def run_europe_pmc_full_text_pilot(
    *,
    coverage_csv: Path,
    output_dir: Path,
    max_species: int = 100,
    min_interval_seconds: float = 0.25,
    aliases_csv: Path | None = None,
    target_species_csv: Path | None = None,
    resolve_wfo_synonyms: bool = False,
    max_aliases_per_species: int = 20,
) -> dict[str, Any]:
    """Run one bounded local pilot over two-axis-complete unresolved species."""

    coverage = pd.read_csv(coverage_csv, dtype=str).fillna("")
    targets = select_reproductive_targets(
        coverage,
        max_species=len(coverage),
    )
    if target_species_csv is not None:
        requested = pd.read_csv(target_species_csv, dtype=str).fillna("")
        if "accepted_species" not in requested.columns:
            raise ValueError("target_species_csv must contain accepted_species")
        requested_species = set(requested["accepted_species"].map(_text))
        targets = targets.loc[
            targets["accepted_species"].isin(requested_species)
        ].copy()
    targets = targets.head(max_species).reset_index(drop=True)
    if targets.empty:
        raise ValueError("no requested species satisfy the target selection contract")
    supplied_aliases = (
        pd.read_csv(aliases_csv, dtype=str).fillna("")
        if aliases_csv is not None
        else None
    )
    wfo_aliases, wfo_errors = (
        resolve_wfo_exact_synonyms(
            targets["accepted_species"].tolist(),
            min_interval_seconds=max(1.0, min_interval_seconds),
            max_aliases_per_species=max_aliases_per_species,
        )
        if resolve_wfo_synonyms
        else (
            pd.DataFrame(),
            pd.DataFrame(columns=["accepted_species", "source", "error"]),
        )
    )
    alias_parts = [
        table
        for table in (supplied_aliases, wfo_aliases)
        if table is not None and not table.empty
    ]
    aliases = pd.concat(alias_parts, ignore_index=True) if alias_parts else None
    search_names = _search_name_table(targets, aliases)
    discovery = discover_europe_pmc_sources(
        search_names["search_name"].tolist(),
        min_interval_seconds=min_interval_seconds,
    )
    remapped_sources = _remap_discovery_sources(discovery.sources, search_names)
    aliases_by_species = {
        accepted: group.loc[
            group["name_role"].ne("accepted"),
            "search_name",
        ].tolist()
        for accepted, group in search_names.groupby("accepted_species")
    }
    full_text = discover_europe_pmc_full_text_sources(
        remapped_sources,
        aliases_by_species=aliases_by_species,
        min_interval_seconds=min_interval_seconds,
    )
    title_abstract_direct = _title_abstract_direct_sources(
        remapped_sources,
        aliases_by_species,
    )

    discovery_sources = pd.DataFrame(
        remapped_sources,
        columns=[
            *DISCOVERY_SOURCE_COLUMNS,
            "searched_name",
            "name_resolution_lineage",
        ],
    )
    discovery_errors = pd.DataFrame(
        discovery.errors,
        columns=[
            "accepted_species",
            "query",
            "search_api_url",
            "error_class",
            "error_code",
            "http_status",
            "message",
            "attempts",
            "retry_after_seconds",
            "retrieved_at_utc",
        ],
    )
    full_text_sources = pd.DataFrame(
        full_text.sources,
        columns=FULL_TEXT_SOURCE_COLUMNS,
    )
    full_text_errors = pd.DataFrame(
        full_text.errors,
        columns=FULL_TEXT_ERROR_COLUMNS,
    )
    all_direct_sources = pd.concat(
        [title_abstract_direct, full_text_sources],
        ignore_index=True,
    ).drop_duplicates(
        ["accepted_species", "source_record_id"],
    )
    shared_sources = _shared_full_text_sources(all_direct_sources)
    extracted = extract_evidence_from_sources(shared_sources)
    evidence = _direct_evidence(extracted, all_direct_sources)

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "targets": output_dir / "target_species.csv",
        "search_names": output_dir / "search_names.csv",
        "wfo_aliases": output_dir / "wfo_exact_aliases.csv.gz",
        "wfo_errors": output_dir / "wfo_alias_errors.csv.gz",
        "discovery_sources": output_dir / "europe_pmc_discovery_sources.csv.gz",
        "discovery_errors": output_dir / "europe_pmc_discovery_errors.csv.gz",
        "title_abstract_direct": (
            output_dir / "europe_pmc_title_abstract_direct_sources.csv.gz"
        ),
        "full_text_sources": output_dir / "europe_pmc_full_text_sources.csv.gz",
        "full_text_errors": output_dir / "europe_pmc_full_text_errors.csv.gz",
        "evidence": output_dir / "direct_web_fulltext_evidence.csv.gz",
    }
    targets.to_csv(paths["targets"], index=False)
    search_names.to_csv(paths["search_names"], index=False)
    wfo_aliases.to_csv(paths["wfo_aliases"], index=False, compression="gzip")
    wfo_errors.to_csv(paths["wfo_errors"], index=False, compression="gzip")
    discovery_sources.to_csv(paths["discovery_sources"], index=False, compression="gzip")
    discovery_errors.to_csv(paths["discovery_errors"], index=False, compression="gzip")
    title_abstract_direct.to_csv(
        paths["title_abstract_direct"],
        index=False,
        compression="gzip",
    )
    full_text_sources.to_csv(paths["full_text_sources"], index=False, compression="gzip")
    full_text_errors.to_csv(paths["full_text_errors"], index=False, compression="gzip")
    evidence.to_csv(paths["evidence"], index=False, compression="gzip")

    pure_gain_by_axis = {TARGET_AXIS: int(evidence["accepted_species"].nunique())}
    summary = {
        "contract_version": CONTRACT_VERSION,
        "input": {
            "coverage_csv": str(coverage_csv),
            "coverage_sha256": _file_sha256(coverage_csv),
            "denominator_species": int(coverage["accepted_species"].nunique()),
            "denominator_species_axis": int(len(coverage)),
        },
        "selection": {
            "rule": (
                "reproductive assurance unresolved; flower colour and floral "
                "structural complexity already filled; exact species names only"
            ),
            "n_target_species": len(targets),
            "n_search_names": len(search_names),
            "n_exact_synonyms": int(search_names["name_role"].ne("accepted").sum()),
            "n_wfo_aliases": len(wfo_aliases),
            "n_wfo_errors": len(wfo_errors),
        },
        "retrieval": {
            "n_title_abstract_sources": len(discovery_sources),
            "n_species_with_title_abstract_sources": int(
                discovery_sources["accepted_species"].nunique()
            ),
            "n_discovery_errors": len(discovery_errors),
            "n_title_abstract_direct_passages": len(title_abstract_direct),
            "n_species_with_title_abstract_direct_passages": int(
                title_abstract_direct["accepted_species"].nunique()
            ),
            "n_full_text_passages": len(full_text_sources),
            "n_species_with_full_text_passages": int(
                full_text_sources["accepted_species"].nunique()
            ),
            "n_full_text_errors": len(full_text_errors),
        },
        "accepted": {
            "quality": "medium",
            "n_evidence_rows": len(evidence),
            "n_species_axis_pure_gain": int(evidence["accepted_species"].nunique()),
            "pure_gain_by_axis": pure_gain_by_axis,
            "family_inference_used": False,
            "global_fallback_used": False,
        },
        "files": {
            key: {
                "path": str(path),
                "sha256": _file_sha256(path),
                "rows": int(
                    {
                        "targets": len(targets),
                        "search_names": len(search_names),
                        "wfo_aliases": len(wfo_aliases),
                        "wfo_errors": len(wfo_errors),
                        "discovery_sources": len(discovery_sources),
                        "discovery_errors": len(discovery_errors),
                        "title_abstract_direct": len(title_abstract_direct),
                        "full_text_sources": len(full_text_sources),
                        "full_text_errors": len(full_text_errors),
                        "evidence": len(evidence),
                    }[key]
                ),
            }
            for key, path in paths.items()
        },
    }
    summary_path = output_dir / "direct_web_fulltext_summary.json"
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


def _coverage_snapshot(coverage: pd.DataFrame, quality_column: str) -> dict[str, Any]:
    filled = coverage[quality_column].ne("")
    by_axis: dict[str, Any] = {}
    for axis, group in coverage.groupby("axis", sort=True):
        counts = group[quality_column].value_counts()
        by_axis[_text(axis)] = {
            "denominator": len(group),
            "filled_species": int(group[quality_column].ne("").sum()),
            "fill_rate": float(group[quality_column].ne("").mean()),
            "quality_counts": {
                quality: int(counts.get(quality, 0))
                for quality in ("high", "medium", "low")
            },
        }
    species_counts = (
        coverage.assign(filled=filled.astype(int))
        .groupby("accepted_species")["filled"]
        .sum()
        .value_counts()
    )
    counts = coverage.loc[filled, quality_column].value_counts()
    return {
        "quality_counts": {
            quality: int(counts.get(quality, 0))
            for quality in ("high", "medium", "low")
        },
        "filled_species_axis": int(filled.sum()),
        "unresolved_species_axis": int((~filled).sum()),
        "by_axis": by_axis,
        "species_by_filled_axis_count": {
            str(value): int(species_counts.get(value, 0))
            for value in range(4)
        },
    }


def combine_direct_web_evidence(
    *,
    coverage_csv: Path,
    evidence_csvs: list[Path],
    output_dir: Path,
) -> dict[str, Any]:
    """Apply direct-web evidence to the fixed ledger with quality precedence."""

    if not evidence_csvs:
        raise ValueError("at least one evidence CSV is required")
    coverage = pd.read_csv(coverage_csv, dtype=str).fillna("")
    required_coverage = {"accepted_species", "axis", "after_quality"}
    if missing := required_coverage.difference(coverage.columns):
        raise ValueError(f"coverage is missing columns: {sorted(missing)}")
    if coverage.duplicated(["accepted_species", "axis"]).any():
        raise ValueError("coverage contains duplicate accepted_species x axis rows")
    coverage["after_quality"] = coverage["after_quality"].map(
        lambda value: _text(value).casefold()
    )
    parts: list[pd.DataFrame] = []
    source_inputs: list[dict[str, Any]] = []
    for path in evidence_csvs:
        part = pd.read_csv(path, dtype=str).fillna("")
        required = {
            "accepted_species",
            "axis",
            "evidence_quality",
            "source_lineage",
        }
        if missing := required.difference(part.columns):
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        part["input_file"] = str(path)
        parts.append(part)
        source_inputs.append(
            {
                "path": str(path),
                "sha256": _file_sha256(path),
                "rows": len(part),
            }
        )
    evidence = pd.concat(parts, ignore_index=True).fillna("")
    evidence["evidence_quality"] = evidence["evidence_quality"].map(
        lambda value: _text(value).casefold()
    )
    invalid_quality = sorted(
        set(evidence["evidence_quality"]).difference({"high", "medium", "low"})
    )
    if invalid_quality:
        raise ValueError(f"invalid evidence quality values: {invalid_quality}")
    universe = pd.MultiIndex.from_frame(
        coverage[["accepted_species", "axis"]]
    )
    evidence_keys = pd.MultiIndex.from_frame(
        evidence[["accepted_species", "axis"]]
    )
    if not evidence_keys.isin(universe).all():
        raise ValueError("direct-web evidence contains species-axis outside the ledger")
    evidence = evidence.drop_duplicates(
        ["accepted_species", "axis", "source_lineage", "evidence_quality"]
    )
    ranked = evidence.assign(
        quality_rank=evidence["evidence_quality"].map(QUALITY_RANK)
    ).sort_values(
        [
            "accepted_species",
            "axis",
            "quality_rank",
            "source_lineage",
        ],
        ascending=[True, True, False, True],
    )
    best = ranked.drop_duplicates(["accepted_species", "axis"], keep="first")
    candidate = best[
        ["accepted_species", "axis", "evidence_quality", "quality_rank"]
    ]
    merged = coverage.merge(
        candidate,
        on=["accepted_species", "axis"],
        how="left",
    ).fillna({"evidence_quality": "", "quality_rank": 0})
    merged["baseline_quality_rank"] = merged["after_quality"].map(QUALITY_RANK)
    promote = merged["quality_rank"].astype(int) > merged["baseline_quality_rank"]
    merged["direct_web_quality"] = ""
    merged.loc[promote, "direct_web_quality"] = merged.loc[
        promote,
        "evidence_quality",
    ]
    merged["combined_quality"] = merged["after_quality"]
    merged.loc[promote, "combined_quality"] = merged.loc[
        promote,
        "evidence_quality",
    ]

    before = _coverage_snapshot(merged, "after_quality")
    after = _coverage_snapshot(merged, "combined_quality")
    pure_gain = merged.loc[
        merged["after_quality"].eq("")
        & merged["direct_web_quality"].ne("")
    ]
    upgrades = merged.loc[
        merged["after_quality"].ne("")
        & merged["direct_web_quality"].ne("")
    ]
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "direct_web_combined_evidence.csv.gz"
    coverage_path = output_dir / "direct_web_combined_coverage.csv.gz"
    evidence.to_csv(evidence_path, index=False, compression="gzip")
    merged.to_csv(coverage_path, index=False, compression="gzip")
    summary = {
        "contract_version": "direct_web_combined_coverage_v1",
        "denominator": {
            "species": int(coverage["accepted_species"].nunique()),
            "axes": int(coverage["axis"].nunique()),
            "species_axis": len(coverage),
        },
        "quality_precedence": ["high", "medium", "low"],
        "before": before,
        "after": after,
        "net_change": {
            "filled_species_axis": (
                after["filled_species_axis"] - before["filled_species_axis"]
            ),
            "unresolved_species_axis": (
                after["unresolved_species_axis"]
                - before["unresolved_species_axis"]
            ),
            "quality_counts": {
                quality: (
                    after["quality_counts"][quality]
                    - before["quality_counts"][quality]
                )
                for quality in ("high", "medium", "low")
            },
            "pure_gain_species_axis": len(pure_gain),
            "quality_upgrades": len(upgrades),
            "pure_gain_by_axis": {
                axis: int(pure_gain["axis"].eq(axis).sum())
                for axis in sorted(coverage["axis"].unique())
            },
        },
        "policy": {
            "family_inference_allowed": False,
            "global_fallback_allowed": False,
            "source_lineage_deduplicated": True,
        },
        "source_inputs": source_inputs,
        "files": {
            evidence_path.name: {
                "sha256": _file_sha256(evidence_path),
                "rows": len(evidence),
            },
            coverage_path.name: {
                "sha256": _file_sha256(coverage_path),
                "rows": len(merged),
            },
        },
    }
    summary_path = output_dir / "direct_web_combined_summary.json"
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


@app.command("europe-pmc-fulltext")
def europe_pmc_fulltext_command(
    coverage_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    output_dir: Path = typer.Option(..., file_okay=False),
    max_species: int = typer.Option(100, min=1, max=10_000),
    min_interval_seconds: float = typer.Option(0.25, min=0.0, max=10.0),
    aliases_csv: Path | None = typer.Option(None, exists=True, dir_okay=False),
    target_species_csv: Path | None = typer.Option(
        None,
        exists=True,
        dir_okay=False,
    ),
    resolve_wfo_synonyms: bool = typer.Option(False),
    max_aliases_per_species: int = typer.Option(20, min=1, max=100),
) -> None:
    """Run the bounded Europe PMC OA full-text reacquisition lane locally."""

    summary = run_europe_pmc_full_text_pilot(
        coverage_csv=coverage_csv,
        output_dir=output_dir,
        max_species=max_species,
        min_interval_seconds=min_interval_seconds,
        aliases_csv=aliases_csv,
        target_species_csv=target_species_csv,
        resolve_wfo_synonyms=resolve_wfo_synonyms,
        max_aliases_per_species=max_aliases_per_species,
    )
    typer.echo(json.dumps(summary, ensure_ascii=False, indent=2))


@app.command("combine")
def combine_command(
    coverage_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    evidence_csv: list[Path] = typer.Option(
        ...,
        "--evidence-csv",
        exists=True,
        dir_okay=False,
    ),
    output_dir: Path = typer.Option(..., file_okay=False),
) -> None:
    """Combine one or more local direct-web evidence artifacts."""

    summary = combine_direct_web_evidence(
        coverage_csv=coverage_csv,
        evidence_csvs=evidence_csv,
        output_dir=output_dir,
    )
    typer.echo(json.dumps(summary, ensure_ascii=False, indent=2))
