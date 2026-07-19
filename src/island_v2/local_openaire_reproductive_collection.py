"""Resumable local OpenAIRE Graph discovery for reproductive assurance.

OpenAIRE title/description metadata is discovery evidence only. Deterministic
matches are emitted as unreviewed leads and never become biological values
without the strict local review/materialization step.
"""

from __future__ import annotations

import fcntl
import hashlib
import html
import json
import math
import os
import re
import ssl
import time
from contextlib import contextmanager
from datetime import UTC, datetime, timedelta
from email.utils import parsedate_to_datetime
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, Iterator
from urllib.parse import quote, urlsplit, urlunsplit

import httpx
import pandas as pd
import truststore
import typer

from island_v2.local_reproductive_assurance_collection import (
    ALLOWED_DIRECT_VALUES,
    CHECKPOINT_COLUMNS,
    CONTRACT_VERSION,
    RA_TRAITS,
    TERMINAL_STATUSES,
    _atomic_write_csv,
    _atomic_write_text,
    _bool,
    _explicit_mapping_supported,
    _file_sha256,
    _nonnegative_int,
    _sha256_text,
    _text,
    _validate_queue,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

OPENAIRE_API_URL = "https://api.openaire.eu/graph/v3/research-products"
PROVIDER = "openaire"
PROVIDER_VERSION = "openaire_graph_v3_three_group_paginated_exact_species_ra_v1"
SOURCE_TYPE = "openaire_title_description"
COLLECTOR_CONTRACT = "local_openaire_reproductive_collection_v1"
GROUP_LEDGER_CONTRACT = "local_openaire_query_group_checkpoint_v1"
PAGINATION_CONTRACT = "all_pages_sequential_until_num_found"

# One AND plus three ORs is exactly four logical operators, the documented
# OpenAIRE maximum for one search expression.
QUERY_GROUPS: dict[str, tuple[str, ...]] = {
    "compatibility": (
        "self-compatible",
        "self-incompatible",
        "self-pollination",
        "self-fertilization",
    ),
    "mating": ("selfing", "autogamy", "outcrossing", "mating system"),
    "assurance": (
        "cleistogamy",
        "mixed mating",
        "autonomous selfing",
        "delayed selfing",
    ),
}
QUERY_GROUP_ORDER = tuple(QUERY_GROUPS)
GROUP_TERMINAL_STATUSES = frozenset({"group_pass_complete", "group_no_result"})
REPRODUCTIVE_TERMS = re.compile(
    r"self[- ]?incompatib|self[- ]?compatib|autonomous self|delayed self|"
    r"mixed[- ]mating|mating system|breeding system|cleistogam|autogam|outcross|"
    r"self[- ]fertili|self[- ]pollinat|obligate self",
    flags=re.IGNORECASE,
)
RETRYABLE_HTTP = frozenset({408, 425, 429})
GROUP_STATUSES = frozenset({"pending", "running", "retry", *GROUP_TERMINAL_STATUSES})
SNAPSHOT_INCONSISTENCY_ERRORS = frozenset(
    {
        "invalid_response:numFound_changed_across_pages",
        "invalid_response:empty_page_before_numFound",
        "invalid_response:short_page_before_numFound",
        "invalid_response:paginated_result_count_incomplete",
    }
)

GROUP_LEDGER_COLUMNS = [
    "species",
    "query_group",
    "query_terms",
    "query",
    "query_sha256",
    "provider",
    "provider_version",
    "status",
    "terminal",
    "attempts",
    "source_count",
    "candidate_count",
    "page_size",
    "num_found",
    "pages_completed",
    "result_count",
    "last_error",
    "requested_at_utc",
    "updated_at",
    "last_batch_id",
    "source_run_id",
    "raw_response_path",
    "raw_response_hash",
    "contract_version",
]

SOURCE_COLUMNS = [
    "accepted_species",
    "source_text",
    "source_url",
    "source_citation",
    "source_type",
    "evidence_scope",
    "openaire_record_id",
    "openaire_doi",
    "openaire_title",
    "openaire_description",
    "openaire_publication_date",
    "openaire_year",
    "openaire_publication",
    "openaire_query_group",
    "raw_response_page",
    "raw_response_row",
    "local_file_path",
    "local_file_hash",
    "retrieval_date",
    "retrieved_at_utc",
    "source_run_id",
    "source_text_hash",
    "source_dedup_key",
]

LEAD_COLUMNS = [
    "source_row_index",
    "accepted_species",
    "trait",
    "provisional_value",
    "source_url",
    "source_citation",
    "source_text",
    "local_file_path",
    "local_file_hash",
    "retrieval_date",
    "review_status",
    "review_note",
]

ERROR_COLUMNS = [
    "species",
    "query_group",
    "provider",
    "batch_id",
    "attempt",
    "error_class",
    "http_status",
    "message",
    "updated_at",
]


class OpenAireTransientError(RuntimeError):
    """A provider, response, or cache failure eligible for a later retry."""

    def __init__(
        self,
        message: str,
        *,
        http_status: int | None = None,
        stop_batch: bool = False,
    ) -> None:
        super().__init__(message)
        self.http_status = http_status
        self.stop_batch = stop_batch


@app.callback()
def main() -> None:
    """Collect exact-species OpenAIRE metadata without assigning traits."""


def _utc_now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _atomic_write_bytes(value: bytes, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        handle.write(value)
        handle.flush()
        os.fsync(handle.fileno())
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


@contextmanager
def _collection_lock(output_dir: Path) -> Iterator[None]:
    """Hold a nonblocking process lock for one OpenAIRE output directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    lock_path = output_dir / ".openaire_collection.lock"
    with lock_path.open("a+", encoding="utf-8") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"another OpenAIRE collector holds {lock_path.resolve()}") from exc
        handle.seek(0)
        handle.truncate()
        handle.write(
            json.dumps(
                {"pid": os.getpid(), "acquired_at_utc": _utc_now()},
                sort_keys=True,
            )
            + "\n"
        )
        handle.flush()
        os.fsync(handle.fileno())
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _species_is_valid(species: str) -> bool:
    value = _text(species)
    return len(value.split()) >= 2 and not any(token in value for token in ('"', "\\", "\r", "\n"))


def _strict_nonnegative_int(value: object, field: str) -> int:
    text = _text(value)
    if not re.fullmatch(r"\d+", text):
        raise ValueError(f"OpenAIRE group ledger {field} is not a nonnegative integer")
    return int(text)


def _strict_bool(value: object, field: str) -> bool:
    if isinstance(value, bool):
        return value
    text = _text(value).casefold()
    if text not in {"true", "false"}:
        raise ValueError(f"OpenAIRE group ledger {field} is not boolean")
    return text == "true"


def _parse_utc_timestamp(value: object, field: str) -> datetime:
    text = _text(value)
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ValueError(f"OpenAIRE {field} timestamp is invalid") from exc
    if not text.endswith("Z") or parsed.tzinfo is None or parsed.utcoffset() != timedelta(0):
        raise ValueError(f"OpenAIRE {field} timestamp is not UTC")
    return parsed


def _strip_markup(value: object) -> str:
    text = html.unescape(str(value or "")).replace("\xa0", " ")
    text = re.sub(r"<[^>]+>", " ", text)
    return _text(text)


def _doi(value: object) -> str:
    doi = _text(value).casefold()
    doi = re.sub(r"^(?:doi:\s*|https?://(?:dx\.)?doi\.org/)", "", doi)
    return doi.rstrip(".,")


def _species_pattern(species: str) -> re.Pattern[str]:
    tokens = _text(species).split()
    expression = r"\s+".join(re.escape(token) for token in tokens)
    return re.compile(rf"(?<!\w){expression}(?!\w)", flags=re.IGNORECASE)


def _explicit_mappings(source_text: str) -> list[tuple[str, str]]:
    return [
        (trait, value)
        for trait in RA_TRAITS
        for value in sorted(ALLOWED_DIRECT_VALUES[trait])
        if _explicit_mapping_supported(trait, value, source_text)
    ]


def _request_query(species: str, query_group: str) -> str:
    if query_group not in QUERY_GROUPS:
        raise ValueError(f"unknown OpenAIRE query group: {query_group}")
    terms = " OR ".join(f'"{term}"' for term in QUERY_GROUPS[query_group])
    query = f'"{_text(species)}" AND ({terms})'
    if len(re.findall(r"\b(?:AND|OR)\b", query)) > 4:
        raise AssertionError("OpenAIRE query exceeds four logical operators")
    return query


def _source_citation(title: str, publication: str, year: str, doi: str) -> str:
    venue = publication or "OpenAIRE indexed work"
    citation = f"{title}. {venue} ({year})."
    if doi:
        citation += f" DOI:{doi}"
    return _text(citation)


def ensure_openaire_checkpoint_rows(
    checkpoint: pd.DataFrame,
    current_queue: pd.DataFrame,
) -> pd.DataFrame:
    """Append external species x trait x provider keys and recover running rows."""
    missing = set(CHECKPOINT_COLUMNS).difference(checkpoint.columns)
    if missing:
        raise ValueError(f"provider checkpoint missing columns: {sorted(missing)}")
    frame = checkpoint[CHECKPOINT_COLUMNS].copy().fillna("")
    if frame.duplicated(["species", "trait", "provider"]).any():
        raise ValueError("provider checkpoint contains duplicate keys")
    queue_species = set(_validate_queue(current_queue))
    known_species = sorted({_text(value) for value in frame["species"] if _text(value)})
    if not queue_species.issubset(set(known_species)):
        raise ValueError("current queue contains species outside the provider checkpoint")

    existing = set(frame[["species", "trait", "provider"]].itertuples(index=False, name=None))
    additions: list[dict[str, object]] = []
    for species in known_species:
        for trait in RA_TRAITS:
            key = (species, trait, PROVIDER)
            if key in existing:
                continue
            if not _species_is_valid(species):
                status = "invalid_species_name"
                basis = "species_name_has_fewer_than_two_tokens"
            elif species not in queue_species:
                status = "provider_seed_covered"
                basis = "reproductive_assurance_already_completed_before_openaire_pass"
            else:
                status = "pending"
                basis = "openaire_literature_provider_unattempted"
            additions.append(
                {
                    "species": species,
                    "trait": trait,
                    "provider": PROVIDER,
                    "provider_version": PROVIDER_VERSION,
                    "status": status,
                    "terminal": status in TERMINAL_STATUSES,
                    "attempts": 0,
                    "candidate_count": 0,
                    "accepted_record_count": 0,
                    "legacy_status": "",
                    "legacy_enabled": "",
                    "last_error": (
                        "permanent:invalid_species_name" if status == "invalid_species_name" else ""
                    ),
                    "updated_at": "",
                    "last_wave_id": "",
                    "source_run_id": "local_openaire",
                    "migration_basis": basis,
                    "contract_version": CONTRACT_VERSION,
                }
            )
    if additions:
        frame = pd.concat([frame, pd.DataFrame(additions)], ignore_index=True)

    frame["terminal"] = frame["terminal"].map(_bool).astype(bool)
    provider_mask = frame["provider"].eq(PROVIDER)
    expected_terminal = frame["status"].isin(TERMINAL_STATUSES)
    if not frame.loc[provider_mask, "terminal"].eq(expected_terminal.loc[provider_mask]).all():
        raise ValueError("OpenAIRE external checkpoint terminal/status invariant failed")
    interrupted = provider_mask & frame["status"].eq("running") & ~frame["terminal"]
    if interrupted.any():
        frame.loc[interrupted, "status"] = "retry"
        frame.loc[interrupted, "terminal"] = False
        frame.loc[interrupted, "last_error"] = "transient:interrupted_previous_local_run"
        frame.loc[interrupted, "migration_basis"] = frame.loc[interrupted, "migration_basis"].map(
            lambda value: f"{_text(value)}|recovered_interrupted_run".strip("|")
        )
    for column in ("attempts", "candidate_count", "accepted_record_count"):
        frame[column] = frame[column].map(_nonnegative_int).astype("int64")
    if frame.duplicated(["species", "trait", "provider"]).any():
        raise AssertionError("extended provider checkpoint contains duplicate keys")
    return frame.sort_values(["species", "trait", "provider"]).reset_index(drop=True)


def ensure_query_group_ledger(
    ledger: pd.DataFrame,
    current_queue: pd.DataFrame,
    page_size: int = 10,
) -> pd.DataFrame:
    """Create and recover the durable species x OpenAIRE query-group ledger."""
    if ledger.empty and not set(GROUP_LEDGER_COLUMNS).issubset(ledger.columns):
        frame = pd.DataFrame(columns=GROUP_LEDGER_COLUMNS)
    else:
        missing = set(GROUP_LEDGER_COLUMNS).difference(ledger.columns)
        if missing:
            raise ValueError(f"OpenAIRE group ledger missing columns: {sorted(missing)}")
        frame = ledger[GROUP_LEDGER_COLUMNS].copy().fillna("")
    if frame.duplicated(["species", "query_group"]).any():
        raise ValueError("OpenAIRE group ledger contains duplicate keys")
    queue_species = _validate_queue(current_queue)
    existing = set(frame[["species", "query_group"]].itertuples(index=False, name=None))
    additions: list[dict[str, object]] = []
    for species in queue_species:
        for query_group in QUERY_GROUP_ORDER:
            if (species, query_group) in existing:
                continue
            additions.append(
                {
                    "species": species,
                    "query_group": query_group,
                    "query_terms": "|".join(QUERY_GROUPS[query_group]),
                    "query": _request_query(species, query_group),
                    "query_sha256": _sha256_text(_request_query(species, query_group)),
                    "provider": PROVIDER,
                    "provider_version": PROVIDER_VERSION,
                    "status": "pending",
                    "terminal": False,
                    "attempts": 0,
                    "source_count": 0,
                    "candidate_count": 0,
                    "page_size": int(page_size),
                    "num_found": 0,
                    "pages_completed": 0,
                    "result_count": 0,
                    "last_error": "",
                    "requested_at_utc": "",
                    "updated_at": "",
                    "last_batch_id": "",
                    "source_run_id": "local_openaire",
                    "raw_response_path": "",
                    "raw_response_hash": "",
                    "contract_version": GROUP_LEDGER_CONTRACT,
                }
            )
    if additions:
        frame = pd.concat([frame, pd.DataFrame(additions)], ignore_index=True)
    if not frame.empty:
        invalid_provider = ~frame["provider"].eq(PROVIDER)
        invalid_group = ~frame["query_group"].isin(QUERY_GROUP_ORDER)
        invalid_terms = frame.apply(
            lambda row: (
                _text(row["query_terms"])
                != "|".join(QUERY_GROUPS.get(_text(row["query_group"]), ()))
            ),
            axis=1,
        )
        invalid_query_hash = frame.apply(
            lambda row: (
                _text(row["query_sha256"])
                != _sha256_text(_request_query(_text(row["species"]), _text(row["query_group"])))
                if _text(row["query_group"]) in QUERY_GROUPS
                else True
            ),
            axis=1,
        )
        invalid_query = frame.apply(
            lambda row: (
                _text(row["query"])
                != _request_query(_text(row["species"]), _text(row["query_group"]))
                if _text(row["query_group"]) in QUERY_GROUPS
                else True
            ),
            axis=1,
        )
        invalid_version = ~frame["provider_version"].eq(PROVIDER_VERSION)
        invalid_contract = ~frame["contract_version"].eq(GROUP_LEDGER_CONTRACT)
        invalid_status = ~frame["status"].isin(GROUP_STATUSES)
        invalid_page_size = frame["page_size"].map(_nonnegative_int).ne(int(page_size))
        if (
            invalid_provider.any()
            or invalid_group.any()
            or invalid_terms.any()
            or invalid_query.any()
            or invalid_query_hash.any()
            or invalid_version.any()
            or invalid_contract.any()
            or invalid_status.any()
            or invalid_page_size.any()
        ):
            raise ValueError("OpenAIRE group ledger violates the query-group contract")
    frame["terminal"] = frame["terminal"].map(lambda value: _strict_bool(value, "terminal"))
    interrupted = frame["status"].eq("running") & ~frame["terminal"]
    if interrupted.any():
        frame.loc[interrupted, "status"] = "retry"
        frame.loc[interrupted, "last_error"] = "transient:interrupted_previous_local_run"
    for column in (
        "attempts",
        "source_count",
        "candidate_count",
        "page_size",
        "num_found",
        "pages_completed",
        "result_count",
    ):
        frame[column] = (
            frame[column]
            .map(lambda value, field=column: _strict_nonnegative_int(value, field))
            .astype("int64")
        )
    for requested_at in frame["requested_at_utc"].map(_text):
        if requested_at:
            _parse_utc_timestamp(requested_at, "group request")
    expected_terminal = frame["status"].isin(GROUP_TERMINAL_STATUSES)
    if not frame["terminal"].eq(expected_terminal).all():
        raise ValueError("OpenAIRE group ledger terminal/status invariant failed")
    terminal = frame["terminal"]
    zero_page = frame["pages_completed"].eq(0)
    page_derived_nonzero = (
        frame[
            [
                "source_count",
                "candidate_count",
                "num_found",
                "result_count",
            ]
        ]
        .gt(0)
        .any(axis=1)
    )
    has_raw_path = frame["raw_response_path"].map(_text).ne("")
    has_raw_hash = frame["raw_response_hash"].map(_text).ne("")
    has_any_raw_provenance = has_raw_path | has_raw_hash
    has_raw_provenance = has_raw_path & frame["raw_response_hash"].map(_text).str.fullmatch(
        r"[0-9a-f]{64}"
    )
    completed_page_span = frame["pages_completed"] * int(page_size)
    if (
        frame.loc[terminal, "pages_completed"].lt(1).any()
        or frame.loc[terminal, "result_count"].ne(frame.loc[terminal, "num_found"]).any()
        or (
            frame.loc[terminal & frame["status"].eq("group_pass_complete"), "source_count"]
            .lt(1)
            .any()
        )
        or (frame.loc[terminal & frame["status"].eq("group_no_result"), "source_count"].ne(0).any())
        or (
            frame.loc[terminal & frame["status"].eq("group_no_result"), "candidate_count"]
            .ne(0)
            .any()
        )
        or (~has_raw_provenance.loc[terminal]).any()
        or completed_page_span.loc[terminal].lt(frame.loc[terminal, "num_found"]).any()
    ):
        raise ValueError("OpenAIRE terminal group ledger progress is inconsistent")
    nonterminal = ~terminal
    if (
        (zero_page & page_derived_nonzero).any()
        or (zero_page & has_any_raw_provenance).any()
        or ((~zero_page) & ~has_raw_provenance).any()
        or (nonterminal & ~zero_page & completed_page_span.ge(frame["num_found"])).any()
        or (frame["status"].eq("pending") & ~zero_page).any()
    ):
        raise ValueError("OpenAIRE nonterminal group ledger progress is inconsistent")
    if frame.duplicated(["species", "query_group"]).any():
        raise AssertionError("extended OpenAIRE group ledger contains duplicate keys")
    return frame.sort_values(["species", "query_group"]).reset_index(drop=True)


def _checkpoint_species_to_process(
    checkpoint: pd.DataFrame,
    current_queue: pd.DataFrame,
    max_species: int,
) -> list[str]:
    queue_order = _validate_queue(current_queue)
    eligible = checkpoint.loc[
        checkpoint["provider"].eq(PROVIDER)
        & checkpoint["trait"].isin(RA_TRAITS)
        & checkpoint["status"].isin({"pending", "retry"})
        & ~checkpoint["terminal"].map(_bool)
    ]
    available = set(eligible["species"].map(_text))
    return [species for species in queue_order if species in available][:max_species]


def _set_external_state(
    checkpoint: pd.DataFrame,
    species: str,
    *,
    status_by_trait: dict[str, str],
    candidate_count_by_trait: dict[str, int],
    batch_id: str,
    updated_at: str,
    error: str = "",
    increment_attempts: bool = False,
    eligible_statuses: frozenset[str] | None = None,
    source_run_id: str | None = None,
) -> pd.DataFrame:
    frame = checkpoint.copy()
    for trait in RA_TRAITS:
        mask = (
            frame["species"].eq(species) & frame["trait"].eq(trait) & frame["provider"].eq(PROVIDER)
        )
        if int(mask.sum()) != 1:
            raise ValueError(f"missing OpenAIRE checkpoint key: {(species, trait)}")
        current = frame.loc[mask].iloc[0]
        if _bool(current["terminal"]):
            continue
        if eligible_statuses is not None and _text(current["status"]) not in eligible_statuses:
            continue
        if increment_attempts:
            frame.loc[mask, "attempts"] = (
                frame.loc[mask, "attempts"].map(_nonnegative_int) + 1
            ).astype("int64")
        status = status_by_trait[trait]
        frame.loc[mask, "status"] = status
        frame.loc[mask, "terminal"] = status in TERMINAL_STATUSES
        frame.loc[mask, "candidate_count"] = int(candidate_count_by_trait.get(trait, 0))
        frame.loc[mask, "last_error"] = error
        frame.loc[mask, "updated_at"] = updated_at
        frame.loc[mask, "last_wave_id"] = batch_id
        frame.loc[mask, "source_run_id"] = source_run_id or batch_id
        frame.loc[mask, "provider_version"] = PROVIDER_VERSION
    return frame


def _set_group_state(
    ledger: pd.DataFrame,
    species: str,
    query_group: str,
    *,
    status: str,
    batch_id: str,
    updated_at: str,
    source_count: int | None = None,
    candidate_count: int | None = None,
    num_found: int | None = None,
    pages_completed: int | None = None,
    result_count: int | None = None,
    error: str = "",
    increment_attempts: bool = False,
    requested_at: str | None = None,
    source_run_id: str | None = None,
    raw_response_path: str = "",
    raw_response_hash: str = "",
) -> pd.DataFrame:
    frame = ledger.copy()
    mask = frame["species"].eq(species) & frame["query_group"].eq(query_group)
    if int(mask.sum()) != 1:
        raise ValueError(f"missing OpenAIRE group key: {(species, query_group)}")
    if _bool(frame.loc[mask, "terminal"].iloc[0]):
        return frame
    if increment_attempts:
        frame.loc[mask, "attempts"] = (
            frame.loc[mask, "attempts"].map(_nonnegative_int) + 1
        ).astype("int64")
    frame.loc[mask, "status"] = status
    frame.loc[mask, "terminal"] = status in GROUP_TERMINAL_STATUSES
    if source_count is not None:
        frame.loc[mask, "source_count"] = int(source_count)
    if candidate_count is not None:
        frame.loc[mask, "candidate_count"] = int(candidate_count)
    if num_found is not None:
        frame.loc[mask, "num_found"] = int(num_found)
    if pages_completed is not None:
        frame.loc[mask, "pages_completed"] = int(pages_completed)
    if result_count is not None:
        frame.loc[mask, "result_count"] = int(result_count)
    frame.loc[mask, "last_error"] = error
    if requested_at is not None:
        frame.loc[mask, "requested_at_utc"] = requested_at
    frame.loc[mask, "updated_at"] = updated_at
    frame.loc[mask, "last_batch_id"] = batch_id
    if source_run_id is not None:
        frame.loc[mask, "source_run_id"] = source_run_id
    if raw_response_path:
        frame.loc[mask, "raw_response_path"] = raw_response_path
    if raw_response_hash:
        frame.loc[mask, "raw_response_hash"] = raw_response_hash
    frame.loc[mask, "provider_version"] = PROVIDER_VERSION
    frame.loc[mask, "contract_version"] = GROUP_LEDGER_CONTRACT
    return frame


def _all_groups_terminal(ledger: pd.DataFrame, species: str) -> bool:
    rows = ledger.loc[ledger["species"].eq(species)]
    return len(rows) == len(QUERY_GROUPS) and rows["terminal"].map(_bool).all()


def _request_identity(
    species: str,
    query_group: str,
    page: int,
    page_size: int,
) -> dict[str, object]:
    return {
        "endpoint": OPENAIRE_API_URL,
        "search": _request_query(species, query_group),
        "query_group": query_group,
        "query_terms": list(QUERY_GROUPS[query_group]),
        "type": "publication",
        "page": int(page),
        "page_size": int(page_size),
        "pagination_contract": PAGINATION_CONTRACT,
    }


def _raw_cache_path(
    raw_dir: Path,
    species: str,
    query_group: str,
    page: int,
    page_size: int,
) -> Path:
    identity = json.dumps(
        {
            "provider_version": PROVIDER_VERSION,
            "species": _text(species),
            **_request_identity(_text(species), query_group, page, page_size),
        },
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    return raw_dir / f"{_sha256_text(identity)}.response.json"


def _receipt_path(raw_path: Path) -> Path:
    return raw_path.with_name(raw_path.name.removesuffix(".response.json") + ".receipt.json")


def _write_raw_receipt(
    path: Path,
    *,
    species: str,
    query_group: str,
    page: int,
    retrieved_at: str,
    batch_id: str,
    page_size: int,
    authenticated: bool,
    raw_content: bytes,
) -> None:
    _atomic_write_bytes(raw_content, path)
    receipt = {
        "contract_version": COLLECTOR_CONTRACT,
        "species": species,
        "query_group": query_group,
        "provider": PROVIDER,
        "provider_version": PROVIDER_VERSION,
        "request": _request_identity(species, query_group, page, page_size),
        "authentication_mode": "bearer_token" if authenticated else "unauthenticated",
        "retrieved_at_utc": retrieved_at,
        "source_run_id": batch_id,
        "http_status": 200,
        "response_path": str(path.resolve()),
        "response_bytes": len(raw_content),
        "response_sha256": hashlib.sha256(raw_content).hexdigest(),
    }
    _atomic_write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        _receipt_path(path),
    )


def _load_raw_receipt(
    path: Path,
    species: str,
    query_group: str,
    page: int,
    page_size: int,
) -> dict[str, Any]:
    receipt = json.loads(_receipt_path(path).read_text(encoding="utf-8"))
    if not isinstance(receipt, dict):
        raise ValueError("OpenAIRE raw receipt must contain an object")
    expected = {
        "species": species,
        "query_group": query_group,
        "provider": PROVIDER,
        "provider_version": PROVIDER_VERSION,
        "contract_version": COLLECTOR_CONTRACT,
    }
    for field, value in expected.items():
        if _text(receipt.get(field)) != value:
            raise ValueError(f"OpenAIRE raw receipt {field} mismatch")
    if type(receipt.get("http_status")) is not int or receipt.get("http_status") != 200:
        raise ValueError("OpenAIRE raw receipt HTTP status mismatch")
    if _text(receipt.get("authentication_mode")) not in {
        "unauthenticated",
        "bearer_token",
    }:
        raise ValueError("OpenAIRE raw receipt authentication mode mismatch")
    if _text(receipt.get("response_path")) != str(path.resolve()):
        raise ValueError("OpenAIRE raw receipt response path mismatch")
    retrieved_at = _text(receipt.get("retrieved_at_utc"))
    _parse_utc_timestamp(retrieved_at, "raw receipt retrieval")
    if not _text(receipt.get("source_run_id")):
        raise ValueError("OpenAIRE raw receipt source run id is blank")
    request = receipt.get("request")
    if not isinstance(request, dict) or request != _request_identity(
        species, query_group, page, page_size
    ):
        raise ValueError("OpenAIRE raw receipt request fingerprint mismatch")
    content = path.read_bytes()
    if _text(receipt.get("response_sha256")) != hashlib.sha256(content).hexdigest():
        raise ValueError("OpenAIRE raw response hash mismatch")
    if type(receipt.get("response_bytes")) is not int or receipt.get("response_bytes") != len(
        content
    ):
        raise ValueError("OpenAIRE raw response byte count mismatch")
    payload = json.loads(content)
    if not isinstance(payload, dict):
        raise ValueError("OpenAIRE raw response payload is invalid")
    return {**receipt, "payload": payload}


def _retry_after_seconds(value: object) -> float | None:
    text = _text(value)
    if not text:
        return None
    try:
        return max(float(text), 0.0)
    except ValueError:
        try:
            target = parsedate_to_datetime(text)
            if target.tzinfo is None:
                target = target.replace(tzinfo=UTC)
            return max((target - datetime.now(UTC)).total_seconds(), 0.0)
        except (TypeError, ValueError, OverflowError):
            return None


def _effective_min_interval(access_token: str, requested: float | None) -> float:
    if access_token:
        interval = 0.51 if requested is None else float(requested)
        if not math.isfinite(interval) or interval < 0.5:
            raise ValueError("authenticated OpenAIRE interval must be at least 0.5 seconds")
        return interval
    interval = 60.1 if requested is None else float(requested)
    if not math.isfinite(interval) or interval < 60.0:
        raise ValueError(
            "OpenAIRE unauthenticated collection must remain at or below 60 requests/hour"
        )
    return interval


def _last_request_time(ledger: pd.DataFrame) -> datetime | None:
    observed: list[datetime] = []
    for value in ledger["requested_at_utc"].map(_text):
        if not value:
            continue
        observed.append(_parse_utc_timestamp(value, "group request"))
    return max(observed) if observed else None


def _query_openaire(
    client: httpx.Client,
    species: str,
    query_group: str,
    *,
    page: int,
    page_size: int,
    access_token: str,
    max_retries: int,
    backoff_seconds: float,
    before_request: Any,
    sleeper: Any,
) -> bytes:
    params = {
        "search": _request_query(species, query_group),
        "type": "publication",
        "page": int(page),
        "pageSize": int(page_size),
    }
    headers = {"Authorization": f"Bearer {access_token}"} if access_token else None
    attempt = 0
    while True:
        before_request()
        try:
            response = client.get(OPENAIRE_API_URL, params=params, headers=headers)
        except httpx.TransportError as exc:
            if attempt >= max_retries:
                raise OpenAireTransientError(
                    f"transport:{type(exc).__name__}:{exc}",
                    stop_batch=True,
                ) from exc
            sleeper(backoff_seconds * (2**attempt))
            attempt += 1
            continue

        status = response.status_code
        if status in RETRYABLE_HTTP or status >= 500:
            retry_after = _retry_after_seconds(response.headers.get("Retry-After"))
            if retry_after is not None and retry_after > 60:
                raise OpenAireTransientError(
                    f"http_status={status}:retry_after={retry_after}",
                    http_status=status,
                    stop_batch=True,
                )
            if attempt >= max_retries:
                raise OpenAireTransientError(
                    f"http_status={status}",
                    http_status=status,
                    stop_batch=True,
                )
            delay = retry_after if retry_after is not None else backoff_seconds * (2**attempt)
            sleeper(delay)
            attempt += 1
            continue
        if 400 <= status < 500:
            raise OpenAireTransientError(
                f"http_status={status}:provider_query_or_authentication_error",
                http_status=status,
                stop_batch=True,
            )
        try:
            response.raise_for_status()
        except httpx.HTTPError as exc:
            raise OpenAireTransientError(
                f"invalid_response:{type(exc).__name__}:{exc}",
                stop_batch=True,
            ) from exc
        return response.content


def _response_page(
    payload: object,
    page: int,
    page_size: int,
) -> tuple[list[object], int]:
    if not isinstance(payload, dict):
        raise OpenAireTransientError("invalid_response:top_level_not_object")
    header = payload.get("header")
    if not isinstance(header, dict):
        raise OpenAireTransientError("invalid_response:header_not_object")
    num_found = header.get("numFound")
    if isinstance(num_found, bool) or not isinstance(num_found, int) or num_found < 0:
        raise OpenAireTransientError("invalid_response:numFound_not_nonnegative_integer")
    response_page = header.get("page")
    if isinstance(response_page, bool) or not isinstance(response_page, int):
        raise OpenAireTransientError("invalid_response:page_not_integer")
    if response_page != int(page):
        raise OpenAireTransientError("invalid_response:unexpected_page")
    response_page_size = header.get("pageSize")
    if isinstance(response_page_size, bool) or not isinstance(response_page_size, int):
        raise OpenAireTransientError("invalid_response:page_size_not_integer")
    if response_page_size != int(page_size):
        raise OpenAireTransientError("invalid_response:unexpected_page_size")
    items = payload.get("results")
    if not isinstance(items, list):
        raise OpenAireTransientError("invalid_response:results_not_list")
    if not all(isinstance(item, dict) for item in items):
        raise OpenAireTransientError("invalid_response:result_member_not_object")
    if len(items) > num_found or len(items) > page_size:
        raise OpenAireTransientError("invalid_response:result_count_inconsistent")
    return items, num_found


def _publication_date(item: dict[str, Any]) -> str:
    candidates = [_text(item.get("publicationDate"))]
    instances = item.get("instances")
    if isinstance(instances, list):
        candidates.extend(
            _text(instance.get("publicationDate"))
            for instance in instances
            if isinstance(instance, dict)
        )
    for value in candidates:
        if re.fullmatch(r"(?:17|18|19|20)\d{2}", value):
            return value
        for date_format in ("%Y-%m", "%Y-%m-%d"):
            try:
                parsed = datetime.strptime(value, date_format)
            except ValueError:
                continue
            if 1700 <= parsed.year <= 2099:
                return value
    return ""


def _result_doi(item: dict[str, Any]) -> str:
    collections: list[object] = [item.get("pids")]
    instances = item.get("instances")
    if isinstance(instances, list):
        collections.extend(
            instance.get("pids") for instance in instances if isinstance(instance, dict)
        )
    for collection in collections:
        if not isinstance(collection, list):
            continue
        for identifier in collection:
            if not isinstance(identifier, dict):
                continue
            if _text(identifier.get("scheme")).casefold() == "doi":
                value = _doi(identifier.get("value"))
                if re.fullmatch(r"10\.\d{4,9}/\S+", value):
                    return value
    original_ids = item.get("originalIds")
    if isinstance(original_ids, list):
        for identifier in original_ids:
            value = _doi(identifier)
            if re.fullmatch(r"10\.\d{4,9}/\S+", value):
                return value
    return ""


def _http_url(value: object) -> str:
    candidate = _text(value)
    if not candidate or any(character.isspace() for character in candidate):
        return ""
    parts = urlsplit(candidate)
    return candidate if parts.scheme.casefold() in {"http", "https"} and bool(parts.netloc) else ""


def _result_url(item: dict[str, Any], doi: str) -> str:
    if doi:
        return f"https://doi.org/{doi}"
    instances = item.get("instances")
    if isinstance(instances, list):
        for instance in instances:
            if not isinstance(instance, dict):
                continue
            urls = instance.get("urls")
            if isinstance(urls, list):
                for value in urls:
                    if url := _http_url(value):
                        return url
    documentation = item.get("documentationUrls")
    if isinstance(documentation, list):
        for value in documentation:
            if url := _http_url(value):
                return url
    record_id = _text(item.get("id"))
    return (
        f"https://explore.openaire.eu/search/publication?pid={quote(record_id, safe='')}"
        if record_id
        else ""
    )


def _publication(item: dict[str, Any]) -> str:
    container = item.get("container")
    if isinstance(container, dict) and _strip_markup(container.get("name")):
        return _strip_markup(container.get("name"))
    if publisher := _strip_markup(item.get("publisher")):
        return publisher
    sources = item.get("sources")
    if isinstance(sources, list):
        for source in sources:
            if value := _strip_markup(source):
                return value
    return ""


def _source_key(doi: str, source_url: str, source_text_hash: str) -> str:
    if doi:
        return f"doi:{doi}"
    if source_url:
        parts = urlsplit(source_url)
        canonical_path = parts.path.rstrip("/") or "/"
        canonical = urlunsplit(
            (
                parts.scheme.casefold(),
                parts.netloc.casefold(),
                canonical_path,
                parts.query,
                "",
            )
        )
        return f"url:{canonical}"
    return f"text:{source_text_hash}"


def openaire_sources_from_receipt(
    receipt: dict[str, Any],
    raw_path: Path,
    page_size: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Extract dated exact-species sources and unreviewed deterministic leads."""
    species = _text(receipt.get("species"))
    query_group = _text(receipt.get("query_group"))
    request = receipt.get("request")
    if not isinstance(request, dict):
        raise ValueError("OpenAIRE receipt request is invalid")
    page = _nonnegative_int(request.get("page"))
    retrieved_at = _text(receipt.get("retrieved_at_utc"))
    retrieval_date = retrieved_at[:10]
    batch_id = _text(receipt.get("source_run_id"))
    items, _num_found = _response_page(receipt.get("payload"), page, page_size)
    local_hash = _file_sha256(raw_path)
    pattern = _species_pattern(species)
    records: dict[str, dict[str, str]] = {}
    mappings: dict[str, set[tuple[str, str]]] = {}

    for raw_row, item in enumerate(items):
        if not isinstance(item, dict):
            continue
        title = _strip_markup(item.get("mainTitle"))
        descriptions = item.get("descriptions")
        description = (
            " ".join(_strip_markup(value) for value in descriptions)
            if isinstance(descriptions, list)
            else _strip_markup(descriptions)
        )
        source_text = _text(f"{title}. {description}")
        if not source_text or pattern.search(source_text) is None:
            continue
        if REPRODUCTIVE_TERMS.search(source_text) is None:
            continue
        publication_date = _publication_date(item)
        if not publication_date:
            continue
        year = publication_date[:4]
        doi = _result_doi(item)
        source_url = _result_url(item, doi)
        if not source_url:
            continue
        publication = _publication(item)
        if not publication:
            continue
        text_hash = _sha256_text(source_text)
        dedup_key = _source_key(doi, source_url, text_hash)
        if dedup_key not in records:
            records[dedup_key] = {
                "accepted_species": species,
                "source_text": source_text,
                "source_url": source_url,
                "source_citation": _source_citation(title, publication, year, doi),
                "source_type": SOURCE_TYPE,
                "evidence_scope": "species_direct",
                "openaire_record_id": _text(item.get("id")),
                "openaire_doi": doi,
                "openaire_title": title,
                "openaire_description": description,
                "openaire_publication_date": publication_date,
                "openaire_year": year,
                "openaire_publication": publication,
                "openaire_query_group": query_group,
                "raw_response_page": str(page),
                "raw_response_row": str(raw_row),
                "local_file_path": str(raw_path.resolve()),
                "local_file_hash": local_hash,
                "retrieval_date": retrieval_date,
                "retrieved_at_utc": retrieved_at,
                "source_run_id": batch_id,
                "source_text_hash": text_hash,
                "source_dedup_key": dedup_key,
            }
            mappings[dedup_key] = set(_explicit_mappings(source_text))

    source_frame = pd.DataFrame(records.values(), columns=SOURCE_COLUMNS).fillna("")
    leads: list[dict[str, str | int]] = []
    for source_row_index, source in source_frame.iterrows():
        for trait, value in sorted(mappings[_text(source["source_dedup_key"])]):
            leads.append(
                {
                    "source_row_index": source_row_index,
                    "accepted_species": species,
                    "trait": trait,
                    "provisional_value": value,
                    "source_url": source["source_url"],
                    "source_citation": source["source_citation"],
                    "source_text": source["source_text"],
                    "local_file_path": source["local_file_path"],
                    "local_file_hash": source["local_file_hash"],
                    "retrieval_date": source["retrieval_date"],
                    "review_status": "unreviewed",
                    "review_note": (
                        "Rule match only; inspect the exact target-species claim before acceptance."
                    ),
                }
            )
    return source_frame, pd.DataFrame(leads, columns=LEAD_COLUMNS).fillna("")


def _archive_one_sided_cache(raw_path: Path) -> None:
    receipt_path = _receipt_path(raw_path)
    if raw_path.exists() == receipt_path.exists():
        return
    incomplete_dir = raw_path.parent / "incomplete"
    incomplete_dir.mkdir(parents=True, exist_ok=True)
    existing = raw_path if raw_path.exists() else receipt_path
    content = existing.read_bytes()
    digest = hashlib.sha256(content).hexdigest()
    archived = incomplete_dir / f"{existing.name}.{digest}.incomplete"
    if not archived.exists():
        _atomic_write_bytes(content, archived)


def _archive_invalid_cache(raw_path: Path) -> None:
    invalid_dir = raw_path.parent / "invalid"
    invalid_dir.mkdir(parents=True, exist_ok=True)
    for existing in (raw_path, _receipt_path(raw_path)):
        if not existing.exists():
            continue
        content = existing.read_bytes()
        digest = hashlib.sha256(content).hexdigest()
        archived = invalid_dir / f"{existing.name}.{digest}.invalid"
        if not archived.exists():
            _atomic_write_bytes(content, archived)


def _invalid_marker_path(raw_path: Path) -> Path:
    return raw_path.with_name(raw_path.name + ".invalidated")


def _invalidate_cache(raw_path: Path, reason: str) -> None:
    _archive_invalid_cache(raw_path)
    _atomic_write_text(_text(reason) + "\n", _invalid_marker_path(raw_path))


def _append_errors(current: pd.DataFrame, path: Path) -> pd.DataFrame:
    prior = (
        pd.read_csv(path, dtype=str).fillna("")
        if path.exists()
        else pd.DataFrame(columns=ERROR_COLUMNS)
    )
    combined = pd.concat([prior, current], ignore_index=True, sort=False).fillna("")
    for column in ERROR_COLUMNS:
        if column not in combined.columns:
            combined[column] = ""
    combined = combined[ERROR_COLUMNS].drop_duplicates(
        ["species", "query_group", "batch_id", "attempt", "message"], keep="first"
    )
    _atomic_write_csv(combined, path)
    return combined.reset_index(drop=True)


def _leads_from_sources(sources: pd.DataFrame) -> pd.DataFrame:
    leads: list[dict[str, str | int]] = []
    for index, source in sources.reset_index(drop=True).iterrows():
        for trait, value in _explicit_mappings(_text(source["source_text"])):
            leads.append(
                {
                    "source_row_index": index,
                    "accepted_species": source["accepted_species"],
                    "trait": trait,
                    "provisional_value": value,
                    "source_url": source["source_url"],
                    "source_citation": source["source_citation"],
                    "source_text": source["source_text"],
                    "local_file_path": source["local_file_path"],
                    "local_file_hash": source["local_file_hash"],
                    "retrieval_date": source["retrieval_date"],
                    "review_status": "unreviewed",
                    "review_note": (
                        "Rule match only; inspect the exact target-species claim before acceptance."
                    ),
                }
            )
    frame = pd.DataFrame(leads, columns=LEAD_COLUMNS).fillna("")
    if not frame.empty:
        frame = frame.drop_duplicates(
            ["source_row_index", "trait", "provisional_value"], keep="first"
        ).sort_values(["source_row_index", "trait", "provisional_value"])
    return frame.reset_index(drop=True)


def _deduplicate_sources(source_frames: list[pd.DataFrame]) -> pd.DataFrame:
    populated = [frame for frame in source_frames if not frame.empty]
    if not populated:
        return pd.DataFrame(columns=SOURCE_COLUMNS)
    combined = pd.concat(populated, ignore_index=True, sort=False).fillna("")
    for column in SOURCE_COLUMNS:
        if column not in combined.columns:
            combined[column] = ""
    return (
        combined[SOURCE_COLUMNS]
        .drop_duplicates(["accepted_species", "source_dedup_key"], keep="first")
        .sort_values(
            ["accepted_species", "source_dedup_key", "local_file_path"],
            kind="stable",
        )
        .reset_index(drop=True)
    )


def _group_evidence_paths(
    output_dir: Path,
    species: str,
    query_group: str,
) -> tuple[Path, Path]:
    directory = output_dir / "query_group_evidence" / _sha256_text(species) / query_group
    return directory / "source_evidence.csv", directory / "candidate_leads.csv"


def _write_group_evidence(
    output_dir: Path,
    species: str,
    query_group: str,
    sources: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    source_path, lead_path = _group_evidence_paths(output_dir, species, query_group)
    deduplicated = _deduplicate_sources([sources])
    leads = _leads_from_sources(deduplicated)
    _atomic_write_csv(deduplicated, source_path)
    _atomic_write_csv(leads, lead_path)
    return deduplicated, leads


def _append_group_evidence(
    output_dir: Path,
    species: str,
    query_group: str,
    sources: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    source_path, _lead_path = _group_evidence_paths(output_dir, species, query_group)
    prior = (
        pd.read_csv(source_path, dtype=str).fillna("")
        if source_path.exists()
        else pd.DataFrame(columns=SOURCE_COLUMNS)
    )
    return _write_group_evidence(
        output_dir,
        species,
        query_group,
        _deduplicate_sources([prior, sources]),
    )


def _rebuild_cumulative_evidence(
    output_dir: Path,
    group_ledger: pd.DataFrame,
    cumulative_source_path: Path,
    cumulative_lead_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    for row in group_ledger.sort_values(["species", "query_group"]).to_dict("records"):
        source_path, _lead_path = _group_evidence_paths(
            output_dir,
            _text(row["species"]),
            _text(row["query_group"]),
        )
        if not source_path.exists():
            continue
        frame = pd.read_csv(source_path, dtype=str).fillna("")
        missing = set(SOURCE_COLUMNS).difference(frame.columns)
        if missing:
            raise ValueError(f"OpenAIRE group evidence missing columns: {sorted(missing)}")
        frames.append(frame[SOURCE_COLUMNS])
    sources = _deduplicate_sources(frames)
    leads = _leads_from_sources(sources)
    _atomic_write_csv(sources, cumulative_source_path)
    _atomic_write_csv(leads, cumulative_lead_path)
    return sources, leads


def _sources_from_completed_pages(
    row: dict[str, Any],
    raw_dir: Path,
    page_size: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    species = _text(row["species"])
    query_group = _text(row["query_group"])
    pages_completed = _nonnegative_int(row["pages_completed"])
    source_frames: list[pd.DataFrame] = []
    observed_num_found: int | None = None
    result_count = 0
    final_path: Path | None = None
    final_receipt: dict[str, Any] | None = None

    for page in range(1, pages_completed + 1):
        raw_path = _raw_cache_path(raw_dir, species, query_group, page, page_size)
        if (
            not raw_path.exists()
            or not _receipt_path(raw_path).exists()
            or _invalid_marker_path(raw_path).exists()
        ):
            raise ValueError(
                "OpenAIRE completed group page is missing, one-sided, or invalidated: "
                f"{species} {query_group} page {page}"
            )
        receipt = _load_raw_receipt(raw_path, species, query_group, page, page_size)
        items, num_found = _response_page(receipt.get("payload"), page, page_size)
        if observed_num_found is None:
            observed_num_found = num_found
        elif observed_num_found != num_found:
            raise ValueError("OpenAIRE completed group numFound changed across cached pages")
        if page * page_size < num_found and len(items) != page_size:
            raise ValueError("OpenAIRE completed group has a short cached page before numFound")
        sources, _page_leads = openaire_sources_from_receipt(receipt, raw_path, page_size)
        source_frames.append(sources)
        result_count += len(items)
        final_path = raw_path
        final_receipt = receipt

    sources = _deduplicate_sources(source_frames)
    leads = _leads_from_sources(sources)
    expected_num_found = _nonnegative_int(row["num_found"])
    expected_result_count = _nonnegative_int(row["result_count"])
    if observed_num_found != expected_num_found or result_count != expected_result_count:
        raise ValueError("OpenAIRE completed group cached counts do not match its ledger")
    if _bool(row["terminal"]):
        if result_count != expected_num_found:
            raise ValueError("OpenAIRE terminal group cached result count is incomplete")
    elif result_count >= expected_num_found:
        raise ValueError("OpenAIRE nonterminal group cache is already complete")
    if len(sources) != _nonnegative_int(row["source_count"]) or len(leads) != _nonnegative_int(
        row["candidate_count"]
    ):
        raise ValueError("OpenAIRE completed group evidence counts do not match its ledger")
    if final_path is None or final_receipt is None:
        raise ValueError("OpenAIRE completed group has no final cached page")
    if _text(row["raw_response_path"]) != str(final_path.resolve()) or _text(
        row["raw_response_hash"]
    ) != _file_sha256(final_path):
        raise ValueError("OpenAIRE completed group final raw provenance mismatch")
    if _text(row["source_run_id"]) != _text(final_receipt.get("source_run_id")):
        raise ValueError("OpenAIRE completed group source run id mismatch")
    return sources, leads


def _validate_completed_group_caches(
    group_ledger: pd.DataFrame,
    raw_dir: Path,
    output_dir: Path,
    page_size: int,
) -> None:
    """Hard-fail if durable page progress cannot be reproduced from exact bytes."""
    validated: list[tuple[str, str, pd.DataFrame]] = []
    for row in group_ledger.to_dict("records"):
        species = _text(row["species"])
        query_group = _text(row["query_group"])
        pages_completed = _nonnegative_int(row["pages_completed"])
        if pages_completed == 0:
            if (
                any(
                    _nonnegative_int(row[column])
                    for column in ("source_count", "candidate_count", "num_found", "result_count")
                )
                or _text(row["raw_response_path"])
                or _text(row["raw_response_hash"])
            ):
                raise ValueError("OpenAIRE zero-page group ledger contains page-derived state")
            validated.append((species, query_group, pd.DataFrame(columns=SOURCE_COLUMNS)))
            continue
        sources, _leads = _sources_from_completed_pages(row, raw_dir, page_size)
        validated.append((species, query_group, sources))

    for species, query_group, sources in validated:
        _write_group_evidence(output_dir, species, query_group, sources)


def _reset_group_snapshot(
    group_ledger: pd.DataFrame,
    *,
    species: str,
    query_group: str,
    current_page: int,
    raw_dir: Path,
    output_dir: Path,
    group_checkpoint_path: Path,
    cumulative_source_path: Path,
    cumulative_lead_path: Path,
    page_size: int,
    batch_id: str,
    reason: str,
) -> pd.DataFrame:
    mask = group_ledger["species"].eq(species) & group_ledger["query_group"].eq(query_group)
    if int(mask.sum()) != 1:
        raise ValueError(f"missing OpenAIRE group key: {(species, query_group)}")
    prior_pages = _nonnegative_int(group_ledger.loc[mask, "pages_completed"].iloc[0])
    for page in range(1, max(prior_pages, current_page) + 1):
        raw_path = _raw_cache_path(raw_dir, species, query_group, page, page_size)
        if raw_path.exists() or _receipt_path(raw_path).exists():
            _invalidate_cache(raw_path, reason)
    _write_group_evidence(
        output_dir,
        species,
        query_group,
        pd.DataFrame(columns=SOURCE_COLUMNS),
    )
    frame = group_ledger.copy()
    frame.loc[mask, "status"] = "retry"
    frame.loc[mask, "terminal"] = False
    for column in (
        "source_count",
        "candidate_count",
        "num_found",
        "pages_completed",
        "result_count",
    ):
        frame.loc[mask, column] = 0
    frame.loc[mask, "last_error"] = f"transient:snapshot_reset:{reason}"
    frame.loc[mask, "updated_at"] = _utc_now()
    frame.loc[mask, "last_batch_id"] = batch_id
    frame.loc[mask, "source_run_id"] = batch_id
    frame.loc[mask, "raw_response_path"] = ""
    frame.loc[mask, "raw_response_hash"] = ""
    _atomic_write_csv(frame, group_checkpoint_path)
    _rebuild_cumulative_evidence(
        output_dir,
        frame,
        cumulative_source_path,
        cumulative_lead_path,
    )
    return frame


def _write_batch_outputs(
    source_frames: list[pd.DataFrame],
    errors: list[dict[str, object]],
    source_path: Path,
    lead_path: Path,
    error_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sources = _deduplicate_sources(source_frames)
    leads = _leads_from_sources(sources)
    error_frame = pd.DataFrame(errors, columns=ERROR_COLUMNS).fillna("")
    _atomic_write_csv(sources, source_path)
    _atomic_write_csv(leads, lead_path)
    _atomic_write_csv(error_frame, error_path)
    return sources, leads, error_frame


def _collect_openaire_batch_locked(
    *,
    queue: pd.DataFrame,
    provider_checkpoint: pd.DataFrame,
    provider_checkpoint_path: Path,
    output_dir: Path,
    max_species: int,
    group_checkpoint_path: Path | None = None,
    page_size: int = 10,
    min_interval_seconds: float | None = None,
    access_token: str | None = None,
    max_retries: int = 2,
    backoff_seconds: float = 2.0,
    client: httpx.Client | None = None,
    sleeper: Any = time.sleep,
) -> dict[str, Any]:
    """Run a bounded, resumable three-group OpenAIRE pass locally."""
    if max_species < 1:
        raise ValueError("OpenAIRE max_species must be positive")
    if page_size < 1 or page_size > 100:
        raise ValueError("OpenAIRE page_size must be between 1 and 100")
    if max_retries < 0:
        raise ValueError("OpenAIRE max_retries must be nonnegative")
    if not math.isfinite(float(backoff_seconds)) or backoff_seconds < 0:
        raise ValueError("OpenAIRE backoff_seconds must be finite and nonnegative")
    token = _text(access_token)
    effective_interval = _effective_min_interval(token, min_interval_seconds)
    checkpoint = ensure_openaire_checkpoint_rows(provider_checkpoint, queue)
    _atomic_write_csv(checkpoint, provider_checkpoint_path)
    group_checkpoint_path = group_checkpoint_path or output_dir / "query_group_checkpoint.csv"
    prior_groups = (
        pd.read_csv(group_checkpoint_path, dtype=str).fillna("")
        if group_checkpoint_path.exists()
        else pd.DataFrame(columns=GROUP_LEDGER_COLUMNS)
    )
    group_ledger = ensure_query_group_ledger(prior_groups, queue, page_size)
    _atomic_write_csv(group_ledger, group_checkpoint_path)
    raw_dir = output_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    cumulative_source_path = output_dir / "source_evidence.csv"
    cumulative_lead_path = output_dir / "candidate_leads.csv"
    cumulative_error_path = output_dir / "lookup_errors.csv"
    _validate_completed_group_caches(group_ledger, raw_dir, output_dir, page_size)
    _rebuild_cumulative_evidence(
        output_dir,
        group_ledger,
        cumulative_source_path,
        cumulative_lead_path,
    )
    species_to_process = _checkpoint_species_to_process(checkpoint, queue, max_species)

    started_at = _utc_now()
    batch_id = f"openaire_{started_at.replace('-', '').replace(':', '').replace('.', '')}"
    batch_dir = output_dir / "batches" / batch_id
    batch_dir.mkdir(parents=True, exist_ok=False)
    source_path = batch_dir / "source_evidence.csv"
    lead_path = batch_dir / "candidate_leads.csv"
    error_path = batch_dir / "lookup_errors.csv"

    source_frames: list[pd.DataFrame] = []
    errors: list[dict[str, object]] = []
    _write_batch_outputs(source_frames, errors, source_path, lead_path, error_path)
    processed_species: list[str] = []
    requested_groups = 0
    reused_raw = 0
    completed_groups = 0
    stop_batch = False
    external_attempted_species: set[str] = set()

    owns_client = client is None
    if client is None:
        context = truststore.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
        client = httpx.Client(
            timeout=45.0,
            follow_redirects=True,
            verify=context,
            headers={
                "User-Agent": (
                    "island-floral-v2/0.1 "
                    "(+https://github.com/zuizui0223/island; OpenAIRE metadata cache)"
                )
            },
        )

    try:
        for species in species_to_process:
            pending_groups = [
                group
                for group in QUERY_GROUP_ORDER
                if not group_ledger.loc[
                    group_ledger["species"].eq(species) & group_ledger["query_group"].eq(group),
                    "terminal",
                ]
                .map(_bool)
                .all()
            ]
            checkpoint = _set_external_state(
                checkpoint,
                species,
                status_by_trait={trait: "running" for trait in RA_TRAITS},
                candidate_count_by_trait={},
                batch_id=batch_id,
                updated_at=_utc_now(),
                increment_attempts=False,
                eligible_statuses=frozenset({"pending", "retry"}),
            )
            _atomic_write_csv(checkpoint, provider_checkpoint_path)
            processed_species.append(species)

            for query_group in pending_groups:
                group_failed = False
                while True:
                    group_mask = group_ledger["species"].eq(species) & group_ledger[
                        "query_group"
                    ].eq(query_group)
                    group_row = group_ledger.loc[group_mask].iloc[0]
                    page = _nonnegative_int(group_row["pages_completed"]) + 1
                    raw_path = _raw_cache_path(raw_dir, species, query_group, page, page_size)
                    invalid_marker = _invalid_marker_path(raw_path)
                    cache_ready = (
                        raw_path.exists()
                        and _receipt_path(raw_path).exists()
                        and not invalid_marker.exists()
                    )
                    cached_receipt: dict[str, Any] | None = None
                    if cache_ready:
                        try:
                            cached_receipt = _load_raw_receipt(
                                raw_path, species, query_group, page, page_size
                            )
                            _response_page(cached_receipt.get("payload"), page, page_size)
                        except (
                            OpenAireTransientError,
                            ValueError,
                            json.JSONDecodeError,
                            OSError,
                        ) as exc:
                            _invalidate_cache(
                                raw_path,
                                f"cached_validation:{type(exc).__name__}:{exc}",
                            )
                            cache_ready = False
                            cached_receipt = None
                    else:
                        _archive_one_sided_cache(raw_path)

                    if cache_ready:
                        group_ledger = _set_group_state(
                            group_ledger,
                            species,
                            query_group,
                            status="running",
                            batch_id=batch_id,
                            updated_at=_utc_now(),
                        )
                        _atomic_write_csv(group_ledger, group_checkpoint_path)

                    def before_request() -> None:
                        nonlocal checkpoint, group_ledger, requested_groups
                        last_request = _last_request_time(group_ledger)
                        if last_request is not None:
                            elapsed = (datetime.now(UTC) - last_request).total_seconds()
                            if elapsed < effective_interval:
                                sleeper(effective_interval - elapsed)
                        requested_at = _utc_now()
                        if species not in external_attempted_species:
                            checkpoint = _set_external_state(
                                checkpoint,
                                species,
                                status_by_trait={trait: "running" for trait in RA_TRAITS},
                                candidate_count_by_trait={},
                                batch_id=batch_id,
                                updated_at=_utc_now(),
                                increment_attempts=True,
                                eligible_statuses=frozenset({"running"}),
                            )
                            _atomic_write_csv(checkpoint, provider_checkpoint_path)
                            external_attempted_species.add(species)
                        group_ledger = _set_group_state(
                            group_ledger,
                            species,
                            query_group,
                            status="running",
                            batch_id=batch_id,
                            updated_at=_utc_now(),
                            increment_attempts=True,
                            requested_at=requested_at,
                        )
                        _atomic_write_csv(group_ledger, group_checkpoint_path)
                        requested_groups += 1

                    page_dir = (
                        batch_dir
                        / "query_groups"
                        / _sha256_text(species)[:16]
                        / query_group
                        / f"page_{page:05d}"
                    )
                    page_dir.mkdir(parents=True, exist_ok=False)
                    try:
                        if cache_ready:
                            if cached_receipt is None:
                                raise AssertionError("validated OpenAIRE receipt is missing")
                            receipt = cached_receipt
                            reused_raw += 1
                        else:
                            raw_content = _query_openaire(
                                client,
                                species,
                                query_group,
                                page=page,
                                page_size=page_size,
                                access_token=token,
                                max_retries=max_retries,
                                backoff_seconds=backoff_seconds,
                                before_request=before_request,
                                sleeper=sleeper,
                            )
                            retrieved_at = _utc_now()
                            _write_raw_receipt(
                                raw_path,
                                species=species,
                                query_group=query_group,
                                page=page,
                                retrieved_at=retrieved_at,
                                batch_id=batch_id,
                                page_size=page_size,
                                authenticated=bool(token),
                                raw_content=raw_content,
                            )
                            invalid_marker.unlink(missing_ok=True)
                            receipt = _load_raw_receipt(
                                raw_path, species, query_group, page, page_size
                            )
                        items, num_found = _response_page(receipt.get("payload"), page, page_size)
                        previous_num_found = _nonnegative_int(group_row["num_found"])
                        if page > 1 and previous_num_found != num_found:
                            raise OpenAireTransientError(
                                "invalid_response:numFound_changed_across_pages",
                                stop_batch=True,
                            )
                        result_count = _nonnegative_int(group_row["result_count"]) + len(items)
                        group_complete = page * page_size >= num_found
                        if not items and not group_complete:
                            raise OpenAireTransientError(
                                "invalid_response:empty_page_before_numFound",
                                stop_batch=True,
                            )
                        if not group_complete and len(items) != page_size:
                            raise OpenAireTransientError(
                                "invalid_response:short_page_before_numFound",
                                stop_batch=True,
                            )
                        if group_complete and result_count != num_found:
                            raise OpenAireTransientError(
                                "invalid_response:paginated_result_count_incomplete",
                                stop_batch=True,
                            )
                        sources, leads = openaire_sources_from_receipt(receipt, raw_path, page_size)
                        group_sources, group_leads = _append_group_evidence(
                            output_dir,
                            species,
                            query_group,
                            sources,
                        )
                        source_count = len(group_sources)
                        candidate_count = len(group_leads)
                        status = (
                            "group_pass_complete"
                            if group_complete and source_count
                            else ("group_no_result" if group_complete else "running")
                        )
                        if not sources.empty:
                            source_frames.append(sources)
                        _atomic_write_csv(sources, page_dir / "source_evidence.csv")
                        _atomic_write_csv(leads, page_dir / "candidate_leads.csv")
                        _atomic_write_csv(
                            pd.DataFrame(columns=ERROR_COLUMNS),
                            page_dir / "lookup_errors.csv",
                        )
                        _rebuild_cumulative_evidence(
                            output_dir,
                            group_ledger,
                            cumulative_source_path,
                            cumulative_lead_path,
                        )
                        _write_batch_outputs(
                            source_frames,
                            errors,
                            source_path,
                            lead_path,
                            error_path,
                        )
                        group_ledger = _set_group_state(
                            group_ledger,
                            species,
                            query_group,
                            status=status,
                            batch_id=batch_id,
                            updated_at=_utc_now(),
                            source_count=source_count,
                            candidate_count=candidate_count,
                            num_found=num_found,
                            pages_completed=page,
                            result_count=result_count,
                            source_run_id=_text(receipt.get("source_run_id")) or batch_id,
                            raw_response_path=str(raw_path.resolve()),
                            raw_response_hash=_file_sha256(raw_path),
                        )
                        _atomic_write_csv(group_ledger, group_checkpoint_path)
                        if group_complete:
                            completed_groups += 1
                            break
                    except (
                        OpenAireTransientError,
                        ValueError,
                        json.JSONDecodeError,
                        OSError,
                    ) as exc:
                        stop_batch = bool(getattr(exc, "stop_batch", False))
                        message = _text(exc)
                        attempt = int(group_ledger.loc[group_mask, "attempts"].iloc[0])
                        error = {
                            "species": species,
                            "query_group": query_group,
                            "provider": PROVIDER,
                            "batch_id": batch_id,
                            "attempt": attempt,
                            "error_class": "transient",
                            "http_status": getattr(exc, "http_status", None) or "",
                            "message": message,
                            "updated_at": _utc_now(),
                        }
                        errors.append(error)
                        error_frame = pd.DataFrame([error], columns=ERROR_COLUMNS)
                        if not (page_dir / "source_evidence.csv").exists():
                            _atomic_write_csv(
                                pd.DataFrame(columns=SOURCE_COLUMNS),
                                page_dir / "source_evidence.csv",
                            )
                        if not (page_dir / "candidate_leads.csv").exists():
                            _atomic_write_csv(
                                pd.DataFrame(columns=LEAD_COLUMNS),
                                page_dir / "candidate_leads.csv",
                            )
                        _atomic_write_csv(error_frame, page_dir / "lookup_errors.csv")
                        _append_errors(error_frame, cumulative_error_path)
                        if message in SNAPSHOT_INCONSISTENCY_ERRORS:
                            retained: list[pd.DataFrame] = []
                            for frame in source_frames:
                                keep = ~(
                                    frame["accepted_species"].eq(species)
                                    & frame["openaire_query_group"].eq(query_group)
                                )
                                if keep.any():
                                    retained.append(frame.loc[keep].copy())
                            source_frames[:] = retained
                        elif (
                            raw_path.exists()
                            and _receipt_path(raw_path).exists()
                            and (
                                isinstance(exc, json.JSONDecodeError)
                                or message.startswith("invalid_response:")
                                or message.startswith("OpenAIRE raw")
                                or "raw receipt" in message
                            )
                        ):
                            _invalidate_cache(
                                raw_path,
                                f"{type(exc).__name__}:{message}",
                            )
                        _write_batch_outputs(
                            source_frames,
                            errors,
                            source_path,
                            lead_path,
                            error_path,
                        )
                        if message in SNAPSHOT_INCONSISTENCY_ERRORS:
                            group_ledger = _reset_group_snapshot(
                                group_ledger,
                                species=species,
                                query_group=query_group,
                                current_page=page,
                                raw_dir=raw_dir,
                                output_dir=output_dir,
                                group_checkpoint_path=group_checkpoint_path,
                                cumulative_source_path=cumulative_source_path,
                                cumulative_lead_path=cumulative_lead_path,
                                page_size=page_size,
                                batch_id=batch_id,
                                reason=message,
                            )
                        else:
                            group_ledger = _set_group_state(
                                group_ledger,
                                species,
                                query_group,
                                status="retry",
                                batch_id=batch_id,
                                updated_at=_utc_now(),
                                error=f"transient:{type(exc).__name__}:{message}",
                            )
                            _atomic_write_csv(group_ledger, group_checkpoint_path)
                        group_failed = True
                        break
                if group_failed:
                    break
            if stop_batch:
                break
    finally:
        if owns_client and client is not None:
            client.close()

    batch_sources, batch_leads, errors_frame = _write_batch_outputs(
        source_frames,
        errors,
        source_path,
        lead_path,
        error_path,
    )
    cumulative_sources, cumulative_leads = _rebuild_cumulative_evidence(
        output_dir,
        group_ledger,
        cumulative_source_path,
        cumulative_lead_path,
    )
    cumulative_errors = _append_errors(pd.DataFrame(columns=ERROR_COLUMNS), cumulative_error_path)

    for species in processed_species:
        species_sources = cumulative_sources.loc[cumulative_sources["accepted_species"].eq(species)]
        counts = {
            trait: sum(
                trait == mapped_trait
                for source_text in species_sources["source_text"].map(_text)
                for mapped_trait, _value in _explicit_mappings(source_text)
            )
            for trait in RA_TRAITS
        }
        if _all_groups_terminal(group_ledger, species):
            states = {
                trait: (
                    "provider_pass_complete_with_candidates"
                    if counts[trait]
                    else (
                        "provider_pass_complete"
                        if not species_sources.empty
                        else "provider_no_result"
                    )
                )
                for trait in RA_TRAITS
            }
            error = ""
        else:
            states = {trait: "retry" for trait in RA_TRAITS}
            failed = group_ledger.loc[
                group_ledger["species"].eq(species) & ~group_ledger["terminal"].map(_bool),
                "last_error",
            ].map(_text)
            error = next((value for value in failed if value), "transient:query_groups_incomplete")
        checkpoint = _set_external_state(
            checkpoint,
            species,
            status_by_trait=states,
            candidate_count_by_trait=counts,
            batch_id=batch_id,
            updated_at=_utc_now(),
            error=error,
            eligible_statuses=frozenset({"running", "retry", "pending"}),
        )
    _atomic_write_csv(checkpoint, provider_checkpoint_path)

    report: dict[str, Any] = {
        "contract_version": COLLECTOR_CONTRACT,
        "provider": PROVIDER,
        "provider_version": PROVIDER_VERSION,
        "query_groups": {key: list(value) for key, value in QUERY_GROUPS.items()},
        "pagination_contract": PAGINATION_CONTRACT,
        "page_size": int(page_size),
        "authentication_mode": "bearer_token" if token else "unauthenticated",
        "effective_min_interval_seconds": effective_interval,
        "batch_id": batch_id,
        "started_at_utc": started_at,
        "finished_at_utc": _utc_now(),
        "scheduled_species": len(species_to_process),
        "processed_species": len(processed_species),
        "network_queries": requested_groups,
        "cached_raw_responses_reused": reused_raw,
        "query_groups_completed_this_batch": completed_groups,
        "source_records": len(batch_sources),
        "candidate_leads": len(batch_leads),
        "transient_errors": len(errors_frame),
        "cumulative_source_records": len(cumulative_sources),
        "cumulative_candidate_leads": len(cumulative_leads),
        "cumulative_errors": len(cumulative_errors),
        "remaining_openaire_species": len(
            _checkpoint_species_to_process(checkpoint, queue, len(queue))
        ),
        "biological_values_materialized": False,
        "external_acquisition_started": requested_groups > 0,
        "global_fallback_used": False,
        "genus_or_family_inference_used": False,
        "outputs": {},
    }
    for path in (
        source_path,
        lead_path,
        error_path,
        cumulative_source_path,
        cumulative_lead_path,
        cumulative_error_path,
        group_checkpoint_path,
        provider_checkpoint_path,
    ):
        report["outputs"][str(path.resolve())] = {
            "path": str(path.resolve()),
            "bytes": path.stat().st_size,
            "sha256": _file_sha256(path),
        }
    report_path = batch_dir / "campaign_status.json"
    _atomic_write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        report_path,
    )
    _atomic_write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        output_dir / "campaign_status.json",
    )
    return report


def collect_openaire_batch(
    *,
    queue: pd.DataFrame,
    provider_checkpoint: pd.DataFrame,
    provider_checkpoint_path: Path,
    output_dir: Path,
    max_species: int,
    group_checkpoint_path: Path | None = None,
    page_size: int = 10,
    min_interval_seconds: float | None = None,
    access_token: str | None = None,
    max_retries: int = 2,
    backoff_seconds: float = 2.0,
    client: httpx.Client | None = None,
    sleeper: Any = time.sleep,
) -> dict[str, Any]:
    """Serialize one resumable OpenAIRE pass for an output directory."""
    with _collection_lock(output_dir.resolve()):
        return _collect_openaire_batch_locked(
            queue=queue,
            provider_checkpoint=provider_checkpoint,
            provider_checkpoint_path=provider_checkpoint_path,
            output_dir=output_dir,
            max_species=max_species,
            group_checkpoint_path=group_checkpoint_path,
            page_size=page_size,
            min_interval_seconds=min_interval_seconds,
            access_token=access_token,
            max_retries=max_retries,
            backoff_seconds=backoff_seconds,
            client=client,
            sleeper=sleeper,
        )


@app.command("collect")
def collect_command(
    queue_csv: Path = typer.Option(..., exists=True),
    provider_checkpoint_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    max_species: int = typer.Option(1, min=1, max=100),
    page_size: int = typer.Option(10, min=1, max=100),
    min_interval_seconds: float | None = typer.Option(None, min=0.5, max=3600.0),
    max_retries: int = typer.Option(2, min=0, max=5),
    backoff_seconds: float = typer.Option(2.0, min=0.0, max=3600.0),
) -> None:
    """Run one bounded, checkpointed OpenAIRE metadata batch locally."""
    queue = pd.read_csv(queue_csv, dtype=str).fillna("")
    checkpoint = pd.read_csv(provider_checkpoint_csv, dtype=str).fillna("")
    report = collect_openaire_batch(
        queue=queue,
        provider_checkpoint=checkpoint,
        provider_checkpoint_path=provider_checkpoint_csv,
        output_dir=output_dir,
        max_species=max_species,
        page_size=page_size,
        min_interval_seconds=min_interval_seconds,
        access_token=os.environ.get("OPENAIRE_ACCESS_TOKEN"),
        max_retries=max_retries,
        backoff_seconds=backoff_seconds,
    )
    typer.echo(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    app()
