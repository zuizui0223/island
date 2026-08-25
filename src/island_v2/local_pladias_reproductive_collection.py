"""Checkpointed exact-taxon Pladias breeding-system discovery.

Only the ``Generative reproduction type`` field on an exact binomial taxon
page is collected.  Pladias values become unreviewed candidate leads; this
collector never materializes a biological record by itself.
"""

from __future__ import annotations

import hashlib
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
from urllib.parse import quote, unquote, urljoin, urlsplit

import httpx
import pandas as pd
import truststore
import typer
from bs4 import BeautifulSoup, Tag

from island_v2.local_reproductive_assurance_collection import (
    CHECKPOINT_COLUMNS,
    CONTRACT_VERSION,
    RA_TRAITS,
    TERMINAL_STATUSES,
    _atomic_write_csv,
    _atomic_write_text,
    _bool,
    _file_sha256,
    _nonnegative_int,
    _sha256_text,
    _text,
    _validate_queue,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

PLADIAS_BASE_URL = "https://pladias.cz/en/taxon/data"
PROVIDER = "pladias"
PROVIDER_VERSION = "pladias_exact_taxon_generative_reproduction_v1"
SOURCE_TYPE = "pladias_generative_reproduction_type"
COLLECTOR_CONTRACT = "local_pladias_reproductive_collection_v1"
EXACT_PAGE_CONTRACT = "url_encoded_exact_binomial_no_taxonomic_fallback"
FIELD_LABEL = "Generative reproduction type"
PUBLICATION = "Chrtek J. Jr. (2018) Generative reproduction type. – www.pladias.cz."
PUBLICATION_YEAR = "2018"

# The mapping is deliberately provider-field specific.  Pladias defines
# facultative allogamy as outcrossing prevailing while selfing remains possible,
# facultative autogamy as mainly selfing with rare outcrossing, and mixed mating
# as both modes being common.  Its unqualified allogamy/autogamy categories are
# the obligate endpoints.  The Island ontology has no separate obligate-
# outcrossing value, so the requested mapping uses predominantly_outcrossing.
CATEGORY_MAPPINGS: dict[str, tuple[str, str]] = {
    "allogamy": ("mating_system", "predominantly_outcrossing"),
    "facultative allogamy": ("mating_system", "predominantly_outcrossing"),
    "autogamy": ("mating_system", "obligate_selfing"),
    "facultative autogamy": ("mating_system", "predominantly_selfing"),
    "mixed mating": ("mating_system", "mixed_mating"),
}
CATEGORY_DEFINITIONS = {
    "allogamy": "Pladias unqualified allogamy category (obligate outcrossing endpoint).",
    "facultative allogamy": "Pladias: outcrossing prevails, but selfing is possible.",
    "autogamy": "Pladias unqualified autogamy category (obligate autogamy endpoint).",
    "facultative autogamy": "Pladias: mainly selfing, while outcrossing is rare.",
    "mixed mating": "Pladias: both outcrossing and selfing are common.",
}
RETRYABLE_HTTP = frozenset({408, 425, 429})
NO_RESULT_HTTP = frozenset({404, 410})
REDIRECT_HTTP = frozenset({301, 302, 303, 307, 308})

SOURCE_COLUMNS = [
    "accepted_species",
    "source_text",
    "source_url",
    "source_citation",
    "source_type",
    "evidence_scope",
    "pladias_field_label",
    "pladias_category",
    "pladias_category_definition",
    "pladias_publication",
    "pladias_year",
    "raw_response_status",
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
    "mapping_basis",
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


class PladiasTransientError(RuntimeError):
    """A provider/cache failure that remains eligible for a later retry."""

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
    """Collect exact-binomial Pladias breeding-system fields locally."""


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
    return len(value.split()) == 2 and not any(
        token in value for token in ('"', "\\", "/", "\r", "\n")
    )


def _request_url(species: str) -> str:
    return f"{PLADIAS_BASE_URL}/{quote(_text(species), safe='')}"


def _request_identity(species: str) -> dict[str, object]:
    return {
        "url": _request_url(species),
        "species": _text(species),
        "field": FIELD_LABEL,
        "exact_page_contract": EXACT_PAGE_CONTRACT,
    }


def _response_url_is_exact(value: object, species: str) -> bool:
    try:
        parsed = urlsplit(_text(value))
    except ValueError:
        return False
    expected_path = f"/en/taxon/data/{_text(species)}"
    return (
        parsed.scheme.casefold() == "https"
        and (parsed.hostname or "").casefold() in {"pladias.cz", "www.pladias.cz"}
        and unquote(parsed.path).rstrip("/") == expected_path
        and not parsed.query
    )


def _response_url_is_no_result_redirect(value: object) -> bool:
    """Recognize Pladias' exact-name miss redirect to its generic taxon index."""
    try:
        parsed = urlsplit(_text(value))
    except ValueError:
        return False
    return (
        parsed.scheme.casefold() == "https"
        and (parsed.hostname or "").casefold() in {"pladias.cz", "www.pladias.cz"}
        and parsed.path.rstrip("/") == "/en/taxon"
        and bool(re.fullmatch(r"_fid=[A-Za-z0-9]+", parsed.query))
        and not parsed.fragment
    )


def _redirect_location_is_no_result(value: object, species: str) -> bool:
    """Recognize only Pladias' observed exact-name-miss redirect target.

    Pladias currently emits an absolute ``http://`` Location even though the
    exact-page request is HTTPS; following it causes another redirect before
    the generic index is downloaded.  Resolving the Location locally lets the
    collector record the actual first response without requesting either
    generic page.  Relative Locations remain supported, but no other host,
    path, query, credentials, port, or fragment is accepted.
    """
    location = _text(value)
    if not location:
        return False
    try:
        parsed = urlsplit(urljoin(_request_url(species), location))
        port = parsed.port
    except ValueError:
        return False
    scheme = parsed.scheme.casefold()
    expected_port = 80 if scheme == "http" else 443 if scheme == "https" else None
    return (
        expected_port is not None
        and (parsed.hostname or "").casefold() in {"pladias.cz", "www.pladias.cz"}
        and parsed.username is None
        and parsed.password is None
        and port in {None, expected_port}
        and parsed.path.rstrip("/") == "/en/taxon"
        and bool(re.fullmatch(r"_fid=[A-Za-z0-9]+", parsed.query))
        and not parsed.fragment
    )


def _raw_cache_path(raw_dir: Path, species: str) -> Path:
    identity = json.dumps(
        {
            "provider_version": PROVIDER_VERSION,
            **_request_identity(species),
        },
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    return raw_dir / f"{_sha256_text(identity)}.response.html"


def _receipt_path(raw_path: Path) -> Path:
    return raw_path.with_name(raw_path.name.removesuffix(".response.html") + ".receipt.json")


def _write_raw_receipt(
    path: Path,
    *,
    species: str,
    retrieved_at: str,
    batch_id: str,
    http_status: int,
    response_url: str,
    raw_content: bytes,
    response_location: str = "",
) -> None:
    _atomic_write_bytes(raw_content, path)
    receipt = {
        "contract_version": COLLECTOR_CONTRACT,
        "species": species,
        "provider": PROVIDER,
        "provider_version": PROVIDER_VERSION,
        "request": _request_identity(species),
        "retrieved_at_utc": retrieved_at,
        "source_run_id": batch_id,
        "http_status": int(http_status),
        "response_url": response_url,
        "response_location": _text(response_location),
        "response_path": str(path.resolve()),
        "response_bytes": len(raw_content),
        "response_sha256": hashlib.sha256(raw_content).hexdigest(),
    }
    _atomic_write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        _receipt_path(path),
    )


def _load_raw_receipt(path: Path, species: str) -> dict[str, Any]:
    receipt = json.loads(_receipt_path(path).read_text(encoding="utf-8"))
    if not isinstance(receipt, dict):
        raise ValueError("Pladias raw receipt must contain an object")
    expected_scalars = {
        "species": species,
        "provider": PROVIDER,
        "provider_version": PROVIDER_VERSION,
        "contract_version": COLLECTOR_CONTRACT,
    }
    for field, expected in expected_scalars.items():
        if _text(receipt.get(field)) != expected:
            raise ValueError(f"Pladias raw receipt {field} mismatch")
    request = receipt.get("request")
    if not isinstance(request, dict) or request != _request_identity(species):
        raise ValueError("Pladias raw receipt request mismatch")
    retrieved_at = _text(receipt.get("retrieved_at_utc"))
    try:
        parsed_retrieval = datetime.fromisoformat(retrieved_at.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ValueError("Pladias raw receipt retrieval timestamp is invalid") from exc
    if (
        not retrieved_at.endswith("Z")
        or parsed_retrieval.utcoffset() != datetime.now(UTC).utcoffset()
    ):
        raise ValueError("Pladias raw receipt retrieval timestamp is not UTC")
    if not _text(receipt.get("source_run_id")):
        raise ValueError("Pladias raw receipt source run ID is blank")
    if _text(receipt.get("response_path")) != str(path.resolve()):
        raise ValueError("Pladias raw receipt response path mismatch")
    status = _nonnegative_int(receipt.get("http_status"))
    if status not in {200, *NO_RESULT_HTTP, *REDIRECT_HTTP}:
        raise ValueError("Pladias raw receipt has non-cacheable HTTP status")
    exact_response = _response_url_is_exact(receipt.get("response_url"), species)
    legacy_followed_miss = status == 200 and _response_url_is_no_result_redirect(
        receipt.get("response_url")
    )
    cached_first_hop_miss = (
        status in REDIRECT_HTTP
        and exact_response
        and _redirect_location_is_no_result(receipt.get("response_location"), species)
    )
    if not (
        legacy_followed_miss
        or cached_first_hop_miss
        or (exact_response and status in {200, *NO_RESULT_HTTP})
    ):
        raise ValueError("Pladias raw receipt response URL is not the exact taxon URL")
    content = path.read_bytes()
    if _text(receipt.get("response_sha256")) != hashlib.sha256(content).hexdigest():
        raise ValueError("Pladias raw response hash mismatch")
    if _nonnegative_int(receipt.get("response_bytes")) != len(content):
        raise ValueError("Pladias raw response byte count mismatch")
    return {**receipt, "payload": content}


def _normalized_category(value: object) -> str:
    normalized = re.sub(r"[\u2010-\u2015\u2212]", "-", _text(value).casefold())
    return normalized.strip(" .:;,|/")


def _category_mappings(category: str) -> list[tuple[str, str, str]]:
    """Map only explicit text from the one Pladias field."""
    folded = _normalized_category(category)
    has_si = bool(re.search(r"\bself[- ]incompatib\w*\b", folded))
    has_sc = bool(re.search(r"\bself[- ]compatib\w*\b", folded))
    base = re.sub(r"\bself[- ](?:incompatib|compatib)\w*\b", " ", folded)
    segments = [
        _text(segment).casefold() for segment in re.split(r"[,;|/()]+", base) if _text(segment)
    ]

    mappings: list[tuple[str, str, str]] = []
    recognized = [segment for segment in segments if segment in CATEGORY_MAPPINGS]
    recognized_values = {CATEGORY_MAPPINGS[segment] for segment in recognized}
    selected_category = ""
    if ("mating_system", "mixed_mating") in recognized_values:
        selected_category = "mixed mating"
    elif len(recognized_values) == 1 and recognized:
        selected_category = recognized[0]
    if selected_category:
        trait, value = CATEGORY_MAPPINGS[selected_category]
        definitions = " | ".join(
            dict.fromkeys(CATEGORY_DEFINITIONS[segment] for segment in recognized)
        )
        mappings.append((trait, value, definitions))
    if has_si and has_sc:
        mappings.append(
            (
                "self_incompatibility",
                "mixed_or_variable",
                "The exact Pladias field explicitly contains both compatibility states.",
            )
        )
    elif has_si:
        mappings.append(
            (
                "self_incompatibility",
                "SI",
                "The exact Pladias field explicitly states self-incompatibility.",
            )
        )
    elif has_sc:
        mappings.append(
            (
                "self_incompatibility",
                "SC",
                "The exact Pladias field explicitly states self-compatibility.",
            )
        )
    return mappings


def _extract_exact_field(content: bytes, species: str) -> str:
    soup = BeautifulSoup(content, "html.parser")
    headings = [_text(node.get_text(" ", strip=True)) for node in soup.select("h1")]
    if not any(name.casefold() == species.casefold() for name in headings):
        raise PladiasTransientError("invalid_response:exact_species_heading_missing")

    # Current Pladias markup stores each feature in ``li.featureDetail``.  The
    # modal below it repeats the label and definitions, so the value must come
    # from the leading summary ``div`` only, never from modal prose.
    for item in soup.select("li.featureDetail"):
        if not isinstance(item, Tag):
            continue
        summary = item.find("div", class_=lambda value: value and "d-flex" in value)
        if not isinstance(summary, Tag):
            continue
        labels = [_text(node.get_text(" ", strip=True)) for node in summary.select("span")]
        if FIELD_LABEL.casefold() not in {label.casefold() for label in labels}:
            continue
        value = summary.find("b")
        return _text(value.get_text(" ", strip=True)) if isinstance(value, Tag) else ""

    # Fail-closed compatibility for semantically equivalent table/list markup.
    for row in soup.select("tr"):
        cells = row.find_all(["th", "td"], recursive=False)
        if (
            len(cells) >= 2
            and _text(cells[0].get_text(" ", strip=True)).rstrip(":").casefold()
            == FIELD_LABEL.casefold()
        ):
            return _text(" ".join(cell.get_text(" ", strip=True) for cell in cells[1:]))
    for term in soup.select("dt"):
        if _text(term.get_text(" ", strip=True)).rstrip(":").casefold() != FIELD_LABEL.casefold():
            continue
        definition = term.find_next_sibling("dd")
        return _text(definition.get_text(" ", strip=True)) if isinstance(definition, Tag) else ""
    return ""


def pladias_sources_from_receipt(
    receipt: dict[str, Any],
    raw_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Extract the exact field and emit review-only candidate mappings."""
    species = _text(receipt.get("species"))
    status = _nonnegative_int(receipt.get("http_status"))
    if status in NO_RESULT_HTTP:
        return pd.DataFrame(columns=SOURCE_COLUMNS), pd.DataFrame(columns=LEAD_COLUMNS)
    if status in REDIRECT_HTTP:
        if _redirect_location_is_no_result(receipt.get("response_location"), species):
            return pd.DataFrame(columns=SOURCE_COLUMNS), pd.DataFrame(columns=LEAD_COLUMNS)
        raise PladiasTransientError("invalid_response:unrecognized_cached_redirect")
    if status != 200:
        raise PladiasTransientError(f"invalid_response:http_status={status}")
    if _response_url_is_no_result_redirect(receipt.get("response_url")):
        return pd.DataFrame(columns=SOURCE_COLUMNS), pd.DataFrame(columns=LEAD_COLUMNS)
    payload = receipt.get("payload")
    if not isinstance(payload, bytes):
        raise PladiasTransientError("invalid_response:payload_not_bytes")
    category = _extract_exact_field(payload, species)
    if not category:
        return pd.DataFrame(columns=SOURCE_COLUMNS), pd.DataFrame(columns=LEAD_COLUMNS)

    retrieved_at = _text(receipt.get("retrieved_at_utc"))
    source_url = _request_url(species)
    source_text = f"{species}. {FIELD_LABEL}: {category}."
    mappings = _category_mappings(category)
    definitions = " | ".join(dict.fromkeys(mapping[2] for mapping in mappings))
    source = pd.DataFrame(
        [
            {
                "accepted_species": species,
                "source_text": source_text,
                "source_url": source_url,
                "source_citation": PUBLICATION,
                "source_type": SOURCE_TYPE,
                "evidence_scope": "species_direct",
                "pladias_field_label": FIELD_LABEL,
                "pladias_category": category,
                "pladias_category_definition": definitions,
                "pladias_publication": PUBLICATION,
                "pladias_year": PUBLICATION_YEAR,
                "raw_response_status": str(status),
                "local_file_path": str(raw_path.resolve()),
                "local_file_hash": _file_sha256(raw_path),
                "retrieval_date": retrieved_at[:10],
                "retrieved_at_utc": retrieved_at,
                "source_run_id": _text(receipt.get("source_run_id")),
                "source_text_hash": _sha256_text(source_text),
            }
        ],
        columns=SOURCE_COLUMNS,
    ).fillna("")
    leads = pd.DataFrame(
        [
            {
                "source_row_index": 0,
                "accepted_species": species,
                "trait": trait,
                "provisional_value": value,
                "source_url": source_url,
                "source_citation": PUBLICATION,
                "source_text": source_text,
                "mapping_basis": basis,
                "local_file_path": str(raw_path.resolve()),
                "local_file_hash": _file_sha256(raw_path),
                "retrieval_date": retrieved_at[:10],
                "review_status": "unreviewed",
                "review_note": (
                    "Exact Pladias field mapping only; inspect the cached taxon page before acceptance."
                ),
            }
            for trait, value, basis in mappings
        ],
        columns=LEAD_COLUMNS,
    ).fillna("")
    return source, leads


def ensure_pladias_checkpoint_rows(
    checkpoint: pd.DataFrame,
    current_queue: pd.DataFrame,
) -> pd.DataFrame:
    """Append provider keys once and recover interrupted local rows."""
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
                basis = "species_name_is_not_an_exact_binomial"
            elif species not in queue_species:
                status = "provider_seed_covered"
                basis = "reproductive_assurance_already_completed_before_pladias_pass"
            else:
                status = "pending"
                basis = "pladias_exact_taxon_field_unattempted"
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
                    "source_run_id": "local_pladias",
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
        frame.loc[interrupted, "migration_basis"] = frame.loc[interrupted, "migration_basis"].map(
            lambda value: f"{_text(value)}|recovered_interrupted_run".strip("|")
        )
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
            frame["species"].eq(species) & frame["trait"].eq(trait) & frame["provider"].eq(PROVIDER)
        )
        if int(mask.sum()) != 1:
            raise ValueError(f"missing Pladias checkpoint key: {(species, trait)}")
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


def _query_pladias(
    client: httpx.Client,
    species: str,
    *,
    max_retries: int,
    backoff_seconds: float,
    min_interval_seconds: float,
    sleeper: Any,
) -> httpx.Response:
    attempt = 0
    while True:
        try:
            # Never fetch Pladias' generic taxon index.  An exact-name miss is
            # conclusively signalled by the first response's tightly validated
            # Location, while an exact taxon page is returned directly as 200.
            response = client.get(_request_url(species), follow_redirects=False)
        except httpx.TransportError as exc:
            if attempt >= max_retries:
                raise PladiasTransientError(f"transport:{type(exc).__name__}:{exc}") from exc
            sleeper(max(min_interval_seconds, backoff_seconds * (2**attempt)))
            attempt += 1
            continue
        status = response.status_code
        if status in RETRYABLE_HTTP or status >= 500:
            retry_after = _retry_after_seconds(response.headers.get("Retry-After"))
            if retry_after is not None and retry_after > 60:
                raise PladiasTransientError(
                    f"http_status={status}:retry_after={retry_after}",
                    http_status=status,
                    stop_batch=True,
                )
            if attempt >= max_retries:
                raise PladiasTransientError(
                    f"http_status={status}",
                    http_status=status,
                    stop_batch=status in RETRYABLE_HTTP,
                )
            delay = retry_after if retry_after is not None else backoff_seconds * (2**attempt)
            sleeper(max(min_interval_seconds, delay))
            attempt += 1
            continue
        if status in REDIRECT_HTTP:
            if _response_url_is_exact(response.url, species) and _redirect_location_is_no_result(
                response.headers.get("Location"), species
            ):
                return response
            raise PladiasTransientError(
                "invalid_response:redirected_outside_exact_taxon_page",
                http_status=status,
                stop_batch=True,
            )
        if status not in {200, *NO_RESULT_HTTP}:
            raise PladiasTransientError(
                f"http_status={status}",
                http_status=status,
                stop_batch=True,
            )
        if not (
            _response_url_is_exact(response.url, species)
            or _response_url_is_no_result_redirect(response.url)
        ):
            raise PladiasTransientError(
                "invalid_response:redirected_outside_exact_taxon_page",
                http_status=status,
                stop_batch=True,
            )
        return response


def _append_deduplicated(
    current: pd.DataFrame,
    path: Path,
    *,
    columns: list[str],
    keys: list[str],
) -> pd.DataFrame:
    prior = (
        pd.read_csv(path, dtype=str).fillna("") if path.exists() else pd.DataFrame(columns=columns)
    )
    combined = pd.concat([prior, current], ignore_index=True, sort=False).fillna("")
    for column in columns:
        if column not in combined.columns:
            combined[column] = ""
    combined = combined[columns].drop_duplicates(keys, keep="first")
    _atomic_write_csv(combined, path)
    return combined


def _archive_cache(raw_path: Path, kind: str) -> None:
    target_dir = raw_path.parent / kind
    target_dir.mkdir(parents=True, exist_ok=True)
    for existing in (raw_path, _receipt_path(raw_path)):
        if not existing.exists():
            continue
        content = existing.read_bytes()
        archived = target_dir / f"{existing.name}.{hashlib.sha256(content).hexdigest()}.{kind}"
        if not archived.exists():
            _atomic_write_bytes(content, archived)


def collect_pladias_batch(
    *,
    queue: pd.DataFrame,
    provider_checkpoint: pd.DataFrame,
    provider_checkpoint_path: Path,
    output_dir: Path,
    max_species: int,
    min_interval_seconds: float = 0.51,
    max_retries: int = 2,
    backoff_seconds: float = 2.0,
    client: httpx.Client | None = None,
    sleeper: Any = time.sleep,
) -> dict[str, Any]:
    """Run one bounded, resumable exact-page Pladias pass."""
    checkpoint = ensure_pladias_checkpoint_rows(provider_checkpoint, queue)
    _atomic_write_csv(checkpoint, provider_checkpoint_path)
    species_to_query = _checkpoint_species_to_query(checkpoint, queue, max_species)
    started_at = _utc_now()
    batch_id = f"pladias_{started_at.replace('-', '').replace(':', '').replace('.', '')}"
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
                    "(+https://github.com/zuizui0223/island; Pladias exact-taxon cache)"
                )
            },
        )

    try:
        for species in species_to_query:
            raw_path = _raw_cache_path(raw_dir, species)
            cache_ready = raw_path.exists() and _receipt_path(raw_path).exists()
            receipt: dict[str, Any] | None = None
            if cache_ready:
                try:
                    receipt = _load_raw_receipt(raw_path, species)
                    pladias_sources_from_receipt(receipt, raw_path)
                except (PladiasTransientError, ValueError, json.JSONDecodeError, OSError):
                    _archive_cache(raw_path, "invalid")
                    cache_ready = False
                    receipt = None
            elif raw_path.exists() or _receipt_path(raw_path).exists():
                _archive_cache(raw_path, "incomplete")
            checkpoint = _set_species_state(
                checkpoint,
                species,
                status_by_trait={trait: "running" for trait in RA_TRAITS},
                candidate_count_by_trait={},
                batch_id=batch_id,
                updated_at=_utc_now(),
                increment_attempts=not cache_ready,
                eligible_statuses=frozenset({"pending", "retry"}),
            )
            _atomic_write_csv(checkpoint, provider_checkpoint_path)
            stop_batch = False
            try:
                if cache_ready:
                    if receipt is None:
                        raise AssertionError("validated Pladias cache receipt is missing")
                    reused_raw += 1
                else:
                    if last_request_at is not None:
                        elapsed = time.monotonic() - last_request_at
                        if elapsed < min_interval_seconds:
                            sleeper(min_interval_seconds - elapsed)
                    queried += 1
                    try:
                        response = _query_pladias(
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
                        http_status=response.status_code,
                        response_url=str(response.url),
                        raw_content=response.content,
                        response_location=response.headers.get("Location", ""),
                    )
                    receipt = _load_raw_receipt(raw_path, species)
                sources, leads = pladias_sources_from_receipt(receipt, raw_path)
                if not sources.empty:
                    source_frames.append(sources)
                if not leads.empty:
                    lead_frames.append(leads)
                counts = (
                    leads.groupby("trait").size().astype(int).to_dict() if not leads.empty else {}
                )
                states = {
                    trait: (
                        "provider_pass_complete_with_candidates"
                        if counts.get(trait, 0)
                        else (
                            "provider_pass_complete" if not sources.empty else "provider_no_result"
                        )
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
            except (PladiasTransientError, ValueError, json.JSONDecodeError, OSError) as exc:
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
        for index, row in leads.iterrows():
            matching = sources.loc[
                sources["accepted_species"].eq(row["accepted_species"])
                & sources["source_url"].eq(row["source_url"])
                & sources["source_text"].eq(row["source_text"])
            ]
            if len(matching) != 1:
                raise AssertionError("Pladias lead does not resolve to one batch source row")
            leads.loc[index, "source_row_index"] = int(matching.index[0])
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

    # Batch evidence and receipts are durable before terminal checkpoint state.
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
        "exact_page_contract": EXACT_PAGE_CONTRACT,
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
        "remaining_pladias_species": len(
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
    min_interval_seconds: float = typer.Option(0.51, min=0.5, max=60.0),
    max_retries: int = typer.Option(2, min=0, max=5),
    backoff_seconds: float = typer.Option(2.0, min=0.0, max=60.0),
) -> None:
    """Run one bounded, checkpointed Pladias exact-taxon batch locally."""
    queue = pd.read_csv(queue_csv, dtype=str).fillna("")
    checkpoint = pd.read_csv(provider_checkpoint_csv, dtype=str).fillna("")
    report = collect_pladias_batch(
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
