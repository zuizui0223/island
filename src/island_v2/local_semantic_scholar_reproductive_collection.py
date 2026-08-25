"""Resumable local Semantic Scholar discovery for reproductive assurance.

The Semantic Scholar Graph API is used only to discover exact-species,
dated title/abstract text.  Regex matches are written as unreviewed leads;
this collector never assigns or materializes a biological trait value.
"""

from __future__ import annotations

import hashlib
import html
import json
import os
import re
import ssl
import time
from datetime import UTC, datetime
from email.utils import parsedate_to_datetime
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any

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

SEMANTIC_SCHOLAR_API_URL = "https://api.semanticscholar.org/graph/v1/paper/search/bulk"
PROVIDER = "semantic_scholar"
PROVIDER_VERSION = "semantic_scholar_graph_bulk_exact_phrase_boolean_ra_v3"
SOURCE_TYPE = "semantic_scholar_title_abstract"
COLLECTOR_CONTRACT = "local_semantic_scholar_reproductive_collection_v1"
PAGINATION_CONTRACT = "first_bulk_response_only_no_token_follow"
PAPER_FIELDS = (
    "title",
    "abstract",
    "year",
    "venue",
    "url",
    "externalIds",
    "publicationDate",
    "openAccessPdf",
    "journal",
)
FIELDS_PARAM = ",".join(PAPER_FIELDS)
REPRODUCTIVE_TERMS = re.compile(
    r"self[- ]?incompatib|self[- ]?compatib|autonomous self|delayed self|"
    r"mixed[- ]mating|mating system|breeding system|cleistogam|autogam|outcross|"
    r"self[- ]fertili|self[- ]pollinat|obligate self",
    flags=re.IGNORECASE,
)
RETRYABLE_HTTP = frozenset({408, 425, 429})

SOURCE_COLUMNS = [
    "accepted_species",
    "source_text",
    "source_url",
    "source_citation",
    "source_type",
    "evidence_scope",
    "semantic_scholar_paper_id",
    "semantic_scholar_doi",
    "semantic_scholar_title",
    "semantic_scholar_abstract",
    "semantic_scholar_year",
    "semantic_scholar_publication_date",
    "semantic_scholar_publication",
    "semantic_scholar_open_access_url",
    "raw_response_row",
    "local_file_path",
    "local_file_hash",
    "retrieval_date",
    "retrieved_at_utc",
    "source_run_id",
    "source_text_hash",
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
    "provider",
    "batch_id",
    "attempt",
    "error_class",
    "http_status",
    "message",
    "updated_at",
]


class SemanticScholarTransientError(RuntimeError):
    """A provider or cache failure that remains eligible for a later retry."""

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
    """Collect exact-species Semantic Scholar metadata without assigning traits."""


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


def _species_is_valid(species: str) -> bool:
    value = _text(species)
    return len(value.split()) >= 2 and not any(
        token in value for token in ('"', "\\", "\r", "\n")
    )


def _strip_markup(value: object) -> str:
    text = html.unescape(str(value or "")).replace("\xa0", " ")
    text = re.sub(r"<[^>]+>", " ", text)
    return _text(text)


def _doi(value: object) -> str:
    doi = _text(value).casefold()
    doi = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", doi)
    return doi.rstrip(".,")


def _publication_year(item: dict[str, Any]) -> str:
    publication_date = _text(item.get("publicationDate"))
    match = re.fullmatch(r"((?:17|18|19|20)\d{2})(?:-\d{2}(?:-\d{2})?)?", publication_date)
    if match:
        return match.group(1)
    year = _text(item.get("year"))
    return year if re.fullmatch(r"(?:17|18|19|20)\d{2}", year) else ""


def _species_pattern(species: str) -> re.Pattern[str]:
    # Python's Unicode-aware ``\w`` supplies the boundary.  This rejects
    # Alpha beta2 while retaining names containing non-ASCII letters.
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


def _request_query(species: str) -> str:
    # Bulk search supports exact phrase, Boolean AND (``+``), and OR (``|``)
    # syntax.  Local exact-species and reproductive-text filters remain the
    # authority for emitting review leads.
    concepts = (
        "(self pollination) | (self fertilization) | (self compatibility) | "
        "(self incompatibility) | selfing | autogamy | outcrossing | "
        "(mating system) | (breeding system) | (mixed mating) | "
        "(delayed selfing) | cleistogamy"
    )
    return f'"{_text(species)}" + ({concepts})'


def _source_citation(title: str, publication: str, year: str, doi: str) -> str:
    venue = publication or "Semantic Scholar indexed work"
    citation = f"{title}. {venue} ({year})."
    if doi:
        citation += f" DOI:{doi}"
    return _text(citation)


def ensure_semantic_scholar_checkpoint_rows(
    checkpoint: pd.DataFrame,
    current_queue: pd.DataFrame,
) -> pd.DataFrame:
    """Append all provider keys once and recover interrupted rows as retry."""
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

    existing = set(
        frame[["species", "trait", "provider"]].itertuples(index=False, name=None)
    )
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
                basis = "reproductive_assurance_already_completed_before_semantic_scholar_pass"
            else:
                status = "pending"
                basis = "semantic_scholar_literature_provider_unattempted"
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
                        "permanent:invalid_species_name"
                        if status == "invalid_species_name"
                        else ""
                    ),
                    "updated_at": "",
                    "last_wave_id": "",
                    "source_run_id": "local_semantic_scholar",
                    "migration_basis": basis,
                    "contract_version": CONTRACT_VERSION,
                }
            )
    if additions:
        frame = pd.concat([frame, pd.DataFrame(additions)], ignore_index=True, sort=False)

    frame["terminal"] = frame["terminal"].map(_bool).astype(bool)
    interrupted = frame["provider"].eq(PROVIDER) & frame["status"].eq("running")
    if interrupted.any():
        frame.loc[interrupted, "status"] = "retry"
        frame.loc[interrupted, "terminal"] = False
        frame.loc[interrupted, "last_error"] = "transient:interrupted_previous_local_run"
        frame.loc[interrupted, "migration_basis"] = frame.loc[
            interrupted, "migration_basis"
        ].map(lambda value: f"{_text(value)}|recovered_interrupted_run".strip("|"))

    for column in ("attempts", "candidate_count", "accepted_record_count"):
        frame[column] = frame[column].map(_nonnegative_int).astype("int64")
    if frame.duplicated(["species", "trait", "provider"]).any():
        raise AssertionError("extended provider checkpoint contains duplicate keys")
    return frame.sort_values(["species", "trait", "provider"]).reset_index(drop=True)


def _checkpoint_species_to_query(
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


def _set_species_state(
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
            frame["species"].eq(species)
            & frame["trait"].eq(trait)
            & frame["provider"].eq(PROVIDER)
        )
        if int(mask.sum()) != 1:
            raise ValueError(f"missing Semantic Scholar checkpoint key: {(species, trait)}")
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


def _request_identity(species: str) -> dict[str, object]:
    return {
        "endpoint": SEMANTIC_SCHOLAR_API_URL,
        "query": _request_query(species),
        "fields": FIELDS_PARAM,
        "pagination_contract": PAGINATION_CONTRACT,
    }


def _raw_cache_path(raw_dir: Path, species: str) -> Path:
    identity = json.dumps(
        {
            "provider_version": PROVIDER_VERSION,
            "species": _text(species),
            **_request_identity(_text(species)),
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
    retrieved_at: str,
    batch_id: str,
    payload: object,
    raw_content: bytes | None = None,
) -> None:
    content = (
        raw_content
        if raw_content is not None
        else json.dumps(
            payload,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    )
    _atomic_write_bytes(content, path)
    receipt = {
        "contract_version": COLLECTOR_CONTRACT,
        "species": species,
        "provider": PROVIDER,
        "provider_version": PROVIDER_VERSION,
        "request": _request_identity(species),
        "retrieved_at_utc": retrieved_at,
        "source_run_id": batch_id,
        "http_status": 200,
        "response_path": str(path.resolve()),
        "response_bytes": len(content),
        "response_sha256": hashlib.sha256(content).hexdigest(),
    }
    _atomic_write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        _receipt_path(path),
    )


def _load_raw_receipt(path: Path, species: str) -> dict[str, Any]:
    receipt = json.loads(_receipt_path(path).read_text(encoding="utf-8"))
    if not isinstance(receipt, dict):
        raise ValueError("Semantic Scholar raw receipt must contain an object")
    expected_scalars = {
        "species": species,
        "provider": PROVIDER,
        "contract_version": COLLECTOR_CONTRACT,
        "provider_version": PROVIDER_VERSION,
    }
    for field, expected in expected_scalars.items():
        if _text(receipt.get(field)) != expected:
            raise ValueError(f"Semantic Scholar raw receipt {field} mismatch")
    request = receipt.get("request")
    if not isinstance(request, dict):
        raise ValueError("Semantic Scholar raw receipt request is invalid")
    expected_request = _request_identity(species)
    if _text(request.get("endpoint")) != expected_request["endpoint"]:
        raise ValueError("Semantic Scholar raw receipt endpoint mismatch")
    if _text(request.get("query")) != expected_request["query"]:
        raise ValueError("Semantic Scholar raw receipt query mismatch")
    if _text(request.get("fields")) != FIELDS_PARAM:
        raise ValueError("Semantic Scholar raw receipt fields mismatch")
    if _text(request.get("pagination_contract")) != PAGINATION_CONTRACT:
        raise ValueError("Semantic Scholar raw receipt pagination contract mismatch")
    content = path.read_bytes()
    if _text(receipt.get("response_sha256")) != hashlib.sha256(content).hexdigest():
        raise ValueError("Semantic Scholar raw response hash mismatch")
    if _nonnegative_int(receipt.get("response_bytes")) != len(content):
        raise ValueError("Semantic Scholar raw response byte count mismatch")
    payload = json.loads(content)
    if not isinstance(payload, dict):
        raise ValueError("Semantic Scholar raw response payload is invalid")
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


def _query_semantic_scholar(
    client: httpx.Client,
    species: str,
    *,
    max_retries: int,
    backoff_seconds: float,
    min_interval_seconds: float,
    sleeper: Any,
) -> tuple[object, bytes]:
    params = {
        "query": _request_query(species),
        "fields": FIELDS_PARAM,
    }
    attempt = 0
    while True:
        try:
            response = client.get(SEMANTIC_SCHOLAR_API_URL, params=params)
        except httpx.TransportError as exc:
            if attempt >= max_retries:
                raise SemanticScholarTransientError(
                    f"transport:{type(exc).__name__}:{exc}"
                ) from exc
            sleeper(max(min_interval_seconds, backoff_seconds * (2**attempt)))
            attempt += 1
            continue

        status = response.status_code
        if status in RETRYABLE_HTTP or status >= 500:
            retry_after = _retry_after_seconds(response.headers.get("Retry-After"))
            if retry_after is not None and retry_after > 60:
                raise SemanticScholarTransientError(
                    f"http_status={status}:retry_after={retry_after}",
                    http_status=status,
                    stop_batch=True,
                )
            if attempt >= max_retries:
                raise SemanticScholarTransientError(
                    f"http_status={status}",
                    http_status=status,
                )
            delay = retry_after if retry_after is not None else backoff_seconds * (2**attempt)
            sleeper(max(min_interval_seconds, delay))
            attempt += 1
            continue
        if 400 <= status < 500:
            # The query shape is collector-wide.  Never make one species
            # terminal because the provider rejects a request/configuration.
            raise SemanticScholarTransientError(
                f"http_status={status}",
                http_status=status,
                stop_batch=True,
            )
        try:
            response.raise_for_status()
            payload = json.loads(response.content)
        except (httpx.HTTPError, ValueError, json.JSONDecodeError) as exc:
            raise SemanticScholarTransientError(
                f"invalid_response:{type(exc).__name__}:{exc}"
            ) from exc
        # Schema validation happens after the exact HTTP bytes and their
        # receipt have been persisted.  Malformed provider responses therefore
        # remain auditable while their checkpoint keys stay retryable.
        return payload, response.content


def _response_items(payload: object) -> list[object]:
    if not isinstance(payload, dict):
        raise SemanticScholarTransientError("invalid_response:top_level_not_object")
    total = payload.get("total")
    if isinstance(total, bool) or not isinstance(total, int) or total < 0:
        raise SemanticScholarTransientError("invalid_response:total_not_nonnegative_integer")
    if "data" not in payload:
        if total == 0:
            return []
        raise SemanticScholarTransientError(
            "invalid_response:data_missing_for_positive_total"
        )
    items = payload["data"]
    if not isinstance(items, list):
        raise SemanticScholarTransientError("invalid_response:data_not_list")
    return items


def semantic_scholar_sources_from_receipt(
    receipt: dict[str, Any],
    raw_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Extract dated exact-species sources and unreviewed regex leads."""
    species = _text(receipt.get("species"))
    retrieved_at = _text(receipt.get("retrieved_at_utc"))
    retrieval_date = retrieved_at[:10]
    batch_id = _text(receipt.get("source_run_id"))
    payload = receipt.get("payload")
    items = _response_items(payload)
    local_hash = _file_sha256(raw_path)
    pattern = _species_pattern(species)
    sources: list[dict[str, str]] = []
    mappings_by_source: list[list[tuple[str, str]]] = []

    for raw_row, item in enumerate(items):
        if not isinstance(item, dict):
            continue
        title = _strip_markup(item.get("title"))
        abstract = _strip_markup(item.get("abstract"))
        source_text = _text(f"{title}. {abstract}")
        if not source_text or pattern.search(source_text) is None:
            continue
        if REPRODUCTIVE_TERMS.search(source_text) is None:
            continue
        year = _publication_year(item)
        if not year:
            continue
        external_ids = item.get("externalIds")
        doi = _doi(external_ids.get("DOI")) if isinstance(external_ids, dict) else ""
        journal = item.get("journal")
        publication = _strip_markup(journal.get("name")) if isinstance(journal, dict) else ""
        if not publication:
            publication = _strip_markup(item.get("venue"))
        open_access = item.get("openAccessPdf")
        open_access_url = (
            _text(open_access.get("url")) if isinstance(open_access, dict) else ""
        )
        source_url = (
            f"https://doi.org/{doi}"
            if doi
            else (_text(item.get("url")) or open_access_url)
        )
        if not source_url:
            continue
        citation = _source_citation(title, publication, year, doi)
        mappings = _explicit_mappings(source_text)
        sources.append(
            {
                "accepted_species": species,
                "source_text": source_text,
                "source_url": source_url,
                "source_citation": citation,
                "source_type": SOURCE_TYPE,
                "evidence_scope": "species_direct",
                "semantic_scholar_paper_id": _text(item.get("paperId")),
                "semantic_scholar_doi": doi,
                "semantic_scholar_title": title,
                "semantic_scholar_abstract": abstract,
                "semantic_scholar_year": year,
                "semantic_scholar_publication_date": _text(item.get("publicationDate")),
                "semantic_scholar_publication": publication,
                "semantic_scholar_open_access_url": open_access_url,
                "raw_response_row": str(raw_row),
                "local_file_path": str(raw_path.resolve()),
                "local_file_hash": local_hash,
                "retrieval_date": retrieval_date,
                "retrieved_at_utc": retrieved_at,
                "source_run_id": batch_id,
                "source_text_hash": _sha256_text(source_text),
            }
        )
        mappings_by_source.append(mappings)

    source_frame = pd.DataFrame(sources, columns=SOURCE_COLUMNS).fillna("")
    leads: list[dict[str, str | int]] = []
    for index, mappings in enumerate(mappings_by_source):
        source = source_frame.iloc[index]
        for trait, value in mappings:
            leads.append(
                {
                    "source_row_index": index,
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
    lead_frame = pd.DataFrame(leads, columns=LEAD_COLUMNS).fillna("")
    if not lead_frame.empty:
        lead_frame = lead_frame.drop_duplicates(
            [
                "accepted_species",
                "source_url",
                "source_text",
                "trait",
                "provisional_value",
            ],
            keep="first",
        ).reset_index(drop=True)
    return source_frame, lead_frame


def _append_deduplicated(
    current: pd.DataFrame,
    path: Path,
    *,
    columns: list[str],
    keys: list[str],
) -> pd.DataFrame:
    prior = (
        pd.read_csv(path, dtype=str).fillna("")
        if path.exists()
        else pd.DataFrame(columns=columns)
    )
    combined = pd.concat([prior, current], ignore_index=True, sort=False).fillna("")
    for column in columns:
        if column not in combined.columns:
            combined[column] = ""
    combined = combined[columns].drop_duplicates(keys, keep="first")
    _atomic_write_csv(combined, path)
    return combined


def collect_semantic_scholar_batch(
    *,
    queue: pd.DataFrame,
    provider_checkpoint: pd.DataFrame,
    provider_checkpoint_path: Path,
    output_dir: Path,
    max_species: int,
    min_interval_seconds: float = 1.05,
    max_retries: int = 2,
    backoff_seconds: float = 2.0,
    client: httpx.Client | None = None,
    sleeper: Any = time.sleep,
) -> dict[str, Any]:
    """Run one bounded pass and durably persist every state transition."""
    checkpoint = ensure_semantic_scholar_checkpoint_rows(provider_checkpoint, queue)
    _atomic_write_csv(checkpoint, provider_checkpoint_path)
    species_to_query = _checkpoint_species_to_query(checkpoint, queue, max_species)
    started_at = _utc_now()
    batch_id = f"semantic_scholar_{started_at.replace('-', '').replace(':', '').replace('.', '')}"
    batch_dir = output_dir / "batches" / batch_id
    raw_dir = output_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    batch_dir.mkdir(parents=True, exist_ok=False)

    source_frames: list[pd.DataFrame] = []
    lead_frames: list[pd.DataFrame] = []
    errors: list[dict[str, object]] = []
    outcomes: dict[str, dict[str, Any]] = {}
    queried = 0
    reused_raw = 0
    last_request_at: float | None = None

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
                    "(+https://github.com/zuizui0223/island; Semantic Scholar metadata cache)"
                )
            },
        )

    try:
        for species in species_to_query:
            raw_path = _raw_cache_path(raw_dir, species)
            checkpoint = _set_species_state(
                checkpoint,
                species,
                status_by_trait={trait: "running" for trait in RA_TRAITS},
                candidate_count_by_trait={},
                batch_id=batch_id,
                updated_at=_utc_now(),
                increment_attempts=not raw_path.exists(),
                eligible_statuses=frozenset({"pending", "retry"}),
            )
            _atomic_write_csv(checkpoint, provider_checkpoint_path)
            stop_batch = False
            try:
                if raw_path.exists():
                    receipt = _load_raw_receipt(raw_path, species)
                    reused_raw += 1
                else:
                    if last_request_at is not None:
                        elapsed = time.monotonic() - last_request_at
                        if elapsed < min_interval_seconds:
                            sleeper(min_interval_seconds - elapsed)
                    queried += 1
                    try:
                        payload, raw_content = _query_semantic_scholar(
                            client,
                            species,
                            max_retries=max_retries,
                            backoff_seconds=backoff_seconds,
                            min_interval_seconds=min_interval_seconds,
                            sleeper=sleeper,
                        )
                    finally:
                        last_request_at = time.monotonic()
                    retrieved_at = _utc_now()
                    _write_raw_receipt(
                        raw_path,
                        species=species,
                        retrieved_at=retrieved_at,
                        batch_id=batch_id,
                        payload=payload,
                        raw_content=raw_content,
                    )
                    receipt = _load_raw_receipt(raw_path, species)
                sources, leads = semantic_scholar_sources_from_receipt(receipt, raw_path)
                if not sources.empty:
                    source_frames.append(sources)
                if not leads.empty:
                    lead_frames.append(leads)
                counts = (
                    leads.groupby("trait").size().astype(int).to_dict()
                    if not leads.empty
                    else {}
                )
                states = {
                    trait: (
                        "provider_pass_complete_with_candidates"
                        if counts.get(trait, 0)
                        else ("provider_pass_complete" if not sources.empty else "provider_no_result")
                    )
                    for trait in RA_TRAITS
                }
                outcomes[species] = {
                    "states": states,
                    "counts": counts,
                    "error": "",
                    "source_run_id": _text(receipt.get("source_run_id")) or batch_id,
                    "updated_at": _utc_now(),
                }
            except (
                SemanticScholarTransientError,
                ValueError,
                json.JSONDecodeError,
                OSError,
            ) as exc:
                stop_batch = bool(getattr(exc, "stop_batch", False))
                outcomes[species] = {
                    "states": {trait: "retry" for trait in RA_TRAITS},
                    "counts": {},
                    "error": f"transient:{type(exc).__name__}:{exc}",
                    "source_run_id": batch_id,
                    "updated_at": _utc_now(),
                }
                attempt = int(
                    checkpoint.loc[
                        checkpoint["species"].eq(species)
                        & checkpoint["trait"].eq(RA_TRAITS[0])
                        & checkpoint["provider"].eq(PROVIDER),
                        "attempts",
                    ].iloc[0]
                )
                errors.append(
                    {
                        "species": species,
                        "provider": PROVIDER,
                        "batch_id": batch_id,
                        "attempt": attempt,
                        "error_class": "transient",
                        "http_status": getattr(exc, "http_status", None) or "",
                        "message": _text(exc),
                        "updated_at": _utc_now(),
                    }
                )
            if stop_batch:
                break
    finally:
        if owns_client and client is not None:
            client.close()

    sources = (
        pd.concat(source_frames, ignore_index=True, sort=False).fillna("")
        if source_frames
        else pd.DataFrame(columns=SOURCE_COLUMNS)
    )
    leads = (
        pd.concat(lead_frames, ignore_index=True, sort=False).fillna("")
        if lead_frames
        else pd.DataFrame(columns=LEAD_COLUMNS)
    )
    errors_frame = pd.DataFrame(errors, columns=ERROR_COLUMNS).fillna("")
    if not sources.empty:
        sources = sources.drop_duplicates(
            ["accepted_species", "source_url", "source_text_hash"], keep="first"
        ).reset_index(drop=True)
    if not leads.empty:
        key_to_index = {
            (row.accepted_species, row.source_url, row.source_text_hash): index
            for index, row in enumerate(sources.itertuples(index=False))
        }
        for index, row in leads.iterrows():
            matching = sources.loc[
                sources["accepted_species"].eq(row["accepted_species"])
                & sources["source_url"].eq(row["source_url"])
                & sources["source_text"].eq(row["source_text"])
            ]
            if len(matching) != 1:
                raise AssertionError(
                    "Semantic Scholar lead does not resolve to one batch source row"
                )
            source_row = matching.iloc[0]
            leads.loc[index, "source_row_index"] = key_to_index[
                (
                    source_row["accepted_species"],
                    source_row["source_url"],
                    source_row["source_text_hash"],
                )
            ]
        leads = leads.drop_duplicates(
            ["source_row_index", "trait", "provisional_value"], keep="first"
        ).sort_values(["source_row_index", "trait", "provisional_value"])

    source_path = batch_dir / "source_evidence.csv"
    lead_path = batch_dir / "candidate_leads.csv"
    error_path = batch_dir / "lookup_errors.csv"
    _atomic_write_csv(sources, source_path)
    _atomic_write_csv(leads, lead_path)
    _atomic_write_csv(errors_frame, error_path)
    cumulative_sources = _append_deduplicated(
        sources,
        output_dir / "source_evidence.csv",
        columns=SOURCE_COLUMNS,
        keys=["accepted_species", "source_url", "source_text_hash"],
    )
    cumulative_errors = _append_deduplicated(
        errors_frame,
        output_dir / "lookup_errors.csv",
        columns=ERROR_COLUMNS,
        keys=["species", "batch_id", "attempt", "error_class", "message"],
    )

    # Raw response, receipt, and batch source/lead files are all durable before
    # any provider key becomes terminal.
    for species, outcome in outcomes.items():
        checkpoint = _set_species_state(
            checkpoint,
            species,
            status_by_trait=outcome["states"],
            candidate_count_by_trait=outcome["counts"],
            batch_id=batch_id,
            updated_at=outcome["updated_at"],
            error=outcome["error"],
            eligible_statuses=frozenset({"running"}),
            source_run_id=outcome["source_run_id"],
        )
    _atomic_write_csv(checkpoint, provider_checkpoint_path)

    report: dict[str, Any] = {
        "contract_version": COLLECTOR_CONTRACT,
        "provider": PROVIDER,
        "provider_version": PROVIDER_VERSION,
        "pagination_contract": PAGINATION_CONTRACT,
        "batch_id": batch_id,
        "started_at_utc": started_at,
        "finished_at_utc": _utc_now(),
        "requested_species": len(outcomes),
        "scheduled_species": len(species_to_query),
        "network_queries": queried,
        "cached_raw_responses_reused": reused_raw,
        "source_records": int(len(sources)),
        "candidate_leads": int(len(leads)),
        "candidate_species": int(leads["accepted_species"].nunique()) if len(leads) else 0,
        "transient_errors": int(len(errors_frame)),
        "permanent_errors": 0,
        "cumulative_source_records": int(len(cumulative_sources)),
        "cumulative_errors": int(len(cumulative_errors)),
        "remaining_semantic_scholar_species": len(
            _checkpoint_species_to_query(checkpoint, queue, len(queue))
        ),
        "biological_values_materialized": False,
        "external_acquisition_started": bool(outcomes),
        "global_fallback_used": False,
        "genus_or_family_inference_used": False,
        "outputs": {},
    }
    for path in (source_path, lead_path, error_path, provider_checkpoint_path):
        key = path.name if path.parent == batch_dir else str(path.resolve())
        report["outputs"][key] = {
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


@app.command("collect")
def collect_command(
    queue_csv: Path = typer.Option(..., exists=True),
    provider_checkpoint_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    max_species: int = typer.Option(25, min=1, max=500),
    min_interval_seconds: float = typer.Option(1.05, min=1.0, max=60.0),
    max_retries: int = typer.Option(2, min=0, max=5),
    backoff_seconds: float = typer.Option(2.0, min=0.0, max=60.0),
) -> None:
    """Run one bounded, checkpointed Semantic Scholar metadata batch locally."""
    queue = pd.read_csv(queue_csv, dtype=str).fillna("")
    checkpoint = pd.read_csv(provider_checkpoint_csv, dtype=str).fillna("")
    report = collect_semantic_scholar_batch(
        queue=queue,
        provider_checkpoint=checkpoint,
        provider_checkpoint_path=provider_checkpoint_csv,
        output_dir=output_dir,
        max_species=max_species,
        min_interval_seconds=min_interval_seconds,
        max_retries=max_retries,
        backoff_seconds=backoff_seconds,
    )
    typer.echo(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    app()
