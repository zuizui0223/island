"""Conservative Europe PMC title/abstract discovery for reproductive evidence.

Only the public Europe PMC REST API is queried.  The adapter freezes article
metadata when the full focal species name and a breeding-system or reproductive
term occur in the returned title/abstract. It never scrapes article pages and
never assigns or infers a plant trait.
"""

from __future__ import annotations

import csv
import hashlib
import html
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from datetime import UTC, datetime
from email.utils import parsedate_to_datetime
from html.parser import HTMLParser
from pathlib import Path
from typing import Any
from xml.etree import ElementTree

EUROPE_PMC_REST_ROOT = "https://www.ebi.ac.uk/europepmc/webservices/rest"
EUROPE_PMC_SEARCH_URL = f"{EUROPE_PMC_REST_ROOT}/search"
EUROPE_PMC_ARTICLE_ROOT = "https://europepmc.org/article"
USER_AGENT = "island-floral-traits/0.1 (Europe PMC REST discovery)"
CONTRACT_VERSION = "europe_pmc_title_abstract_reproduction_v2"
FULL_TEXT_CONTRACT_VERSION = "europe_pmc_species_direct_full_text_v1"
MAX_FULL_TEXT_BYTES = 25_000_000

REPRODUCTIVE_QUERY_TERMS = (
    '"breeding system"',
    '"mating system"',
    "selfing",
    "outcross*",
    "autogam*",
    '"self incompatibility"',
    '"self-incompatibility"',
    '"self compatibility"',
    '"self-compatibility"',
    "cleistogam*",
    '"reproductive biology"',
)

REPRODUCTIVE_TEXT_PATTERNS = tuple(
    re.compile(pattern, flags=re.IGNORECASE)
    for pattern in (
        r"\bbreeding systems?\b",
        r"\bmating systems?\b",
        r"\bselfing\b",
        r"\boutcross\w*\b",
        r"\bautogam\w*\b",
        r"\bself[- ]incompatib\w*\b",
        r"\bself[- ]compatib\w*\b",
        r"\bcleistogam\w*\b",
        r"\breproductive biology\b",
    )
)

FULL_TEXT_DIRECT_PATTERNS = tuple(
    re.compile(pattern, flags=re.IGNORECASE)
    for pattern in (
        r"\bbreeding systems?\b",
        r"\bmating systems?\b",
        r"\bselfing\b",
        r"\boutcross\w*\b",
        r"\bautogam\w*\b",
        r"\bself[- ]incompatib\w*\b",
        r"\bself[- ]compatib\w*\b",
        r"\bcleistogam\w*\b",
        r"\bself[- ]fertili[sz]\w*\b",
        r"\bcross[- ]fertili[sz]\w*\b",
        r"\bautonomous self[- ]pollination\b",
    )
)

RETRYABLE_HTTP_STATUSES = {408, 425, 429, 500, 502, 503, 504}

SOURCE_COLUMNS = [
    "accepted_species",
    "source_text",
    "title",
    "abstract",
    "source_url",
    "article_url",
    "article_api_url",
    "search_api_url",
    "source_citation",
    "source_type",
    "evidence_scope",
    "pmid",
    "pmcid",
    "doi",
    "europe_pmc_source",
    "europe_pmc_id",
    "author_string",
    "journal_title",
    "publication_year",
    "first_publication_date",
    "license",
    "is_open_access",
    "query",
    "search_result_rank",
    "retrieved_at_utc",
    "source_text_sha256",
    "metadata_sha256",
]

ERROR_COLUMNS = [
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
]

FULL_TEXT_SOURCE_COLUMNS = [
    "accepted_species",
    "source_text",
    "source_url",
    "source_citation",
    "source_type",
    "evidence_scope",
    "title",
    "pmcid",
    "doi",
    "license",
    "is_open_access",
    "full_text_xml_url",
    "full_text_xml_sha256",
    "passage_sha256",
    "source_record_id",
    "source_lineage",
    "searched_name",
    "matched_name",
    "name_match_method",
    "name_resolution_lineage",
    "retrieved_at_utc",
]

FULL_TEXT_ERROR_COLUMNS = [
    "accepted_species",
    "pmcid",
    "full_text_xml_url",
    "error_class",
    "error_code",
    "http_status",
    "message",
    "attempts",
    "retry_after_seconds",
    "retrieved_at_utc",
]


@dataclass(frozen=True)
class HttpResponse:
    """Small provider-neutral HTTP response used by the retry loop."""

    status: int
    headers: Mapping[str, str]
    body: bytes


Transport = Callable[[str, float], HttpResponse]
SleepHook = Callable[[float], None]
ClockHook = Callable[[], float]


@dataclass(frozen=True)
class DiscoveryBatch:
    """Returned source records and fully accounted discovery errors."""

    requested_species: tuple[str, ...]
    sources: tuple[dict[str, str], ...]
    errors: tuple[dict[str, str], ...]


@dataclass(frozen=True)
class FullTextBatch:
    """Species-direct passages and fully accounted full-text misses."""

    requested_articles: tuple[tuple[str, str], ...]
    sources: tuple[dict[str, str], ...]
    errors: tuple[dict[str, str], ...]


class EuropePmcRequestError(RuntimeError):
    """A classified Europe PMC request failure."""

    def __init__(
        self,
        message: str,
        *,
        retryable: bool,
        code: str,
        attempts: int,
        http_status: int | None = None,
        retry_after_seconds: float | None = None,
    ) -> None:
        super().__init__(message)
        self.retryable = retryable
        self.code = code
        self.attempts = attempts
        self.http_status = http_status
        self.retry_after_seconds = retry_after_seconds


@dataclass
class BoundedRateLimiter:
    """Serial request pacing with injectable sleep and monotonic-clock hooks."""

    min_interval_seconds: float = 0.25
    sleeper: SleepHook = time.sleep
    clock: ClockHook = time.monotonic
    _last_request_at: float | None = field(default=None, init=False)

    def __post_init__(self) -> None:
        if not 0 <= self.min_interval_seconds <= 10:
            raise ValueError("min_interval_seconds must be between 0 and 10")

    def wait(self) -> None:
        now = self.clock()
        if self._last_request_at is not None:
            remaining = self.min_interval_seconds - (now - self._last_request_at)
            if remaining > 0:
                self.sleeper(remaining)
        self._last_request_at = self.clock()


class _PlainTextParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.parts: list[str] = []

    def handle_data(self, data: str) -> None:
        self.parts.append(data)


def _text(value: object) -> str:
    if value is None:
        return ""
    return " ".join(str(value).replace("\xa0", " ").split())


def _plain_text(value: object) -> str:
    decoded = html.unescape(str(value or ""))
    parser = _PlainTextParser()
    try:
        parser.feed(decoded)
        parser.close()
    except Exception:  # pragma: no cover - HTMLParser is deliberately tolerant
        return _text(re.sub(r"<[^>]+>", " ", decoded))
    text = _text(" ".join(parser.parts))
    return re.sub(r"\s+([,.;:!?])", r"\1", text)


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256_text(value: str) -> str:
    return _sha256_bytes(value.encode("utf-8"))


def _canonical_json_sha256(value: object) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return _sha256_bytes(encoded)


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _normalise_species(species: object) -> str:
    value = _text(species)
    if len(value.split()) < 2:
        raise ValueError("species must contain at least a genus and specific epithet")
    if any(character in value for character in ('"', "\\", "\r", "\n")):
        raise ValueError("species contains an unsafe Europe PMC query character")
    return value


def build_query(species: str) -> str:
    """Build one exact, title/abstract-scoped reproductive literature query."""

    focal = _normalise_species(species)
    term_query = " OR ".join(f"TITLE_ABS:{term}" for term in REPRODUCTIVE_QUERY_TERMS)
    return f'TITLE_ABS:"{focal}" AND ({term_query})'


def build_search_api_url(species: str, *, page_size: int = 25) -> str:
    if not 1 <= page_size <= 100:
        raise ValueError("page_size must be between 1 and 100")
    parameters = [
        ("query", build_query(species)),
        ("resultType", "core"),
        ("pageSize", str(page_size)),
        ("format", "json"),
    ]
    return f"{EUROPE_PMC_SEARCH_URL}?{urllib.parse.urlencode(parameters)}"


def _stdlib_transport(url: str, timeout_seconds: float) -> HttpResponse:
    request = urllib.request.Request(
        url,
        headers={"Accept": "*/*", "User-Agent": USER_AGENT},
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
            return HttpResponse(
                status=int(response.status),
                headers={str(key): str(value) for key, value in response.headers.items()},
                body=response.read(),
            )
    except urllib.error.HTTPError as exc:
        return HttpResponse(
            status=int(exc.code),
            headers={str(key): str(value) for key, value in (exc.headers or {}).items()},
            body=exc.read(),
        )


def _header(headers: Mapping[str, str], name: str) -> str:
    wanted = name.casefold()
    for key, value in headers.items():
        if str(key).casefold() == wanted:
            return _text(value)
    return ""


def _retry_after_seconds(
    value: object,
    *,
    now: datetime | None = None,
) -> float | None:
    raw = _text(value)
    if not raw:
        return None
    try:
        return max(0.0, float(raw))
    except ValueError:
        pass
    try:
        target = parsedate_to_datetime(raw)
    except (TypeError, ValueError, OverflowError):
        return None
    if target.tzinfo is None:
        target = target.replace(tzinfo=UTC)
    reference = now or datetime.now(UTC)
    if reference.tzinfo is None:
        reference = reference.replace(tzinfo=UTC)
    return max(0.0, (target - reference).total_seconds())


def _backoff_delay(attempt: int, base_seconds: float, maximum_seconds: float) -> float:
    return min(maximum_seconds, base_seconds * (2 ** max(0, attempt - 1)))


def _validate_request_bounds(
    *,
    timeout_seconds: float,
    max_attempts: int,
    backoff_seconds: float,
    max_backoff_seconds: float,
) -> None:
    if not 1 <= timeout_seconds <= 120:
        raise ValueError("timeout_seconds must be between 1 and 120")
    if not 1 <= max_attempts <= 5:
        raise ValueError("max_attempts must be between 1 and 5")
    if not 0 <= backoff_seconds <= 30:
        raise ValueError("backoff_seconds must be between 0 and 30")
    if not backoff_seconds <= max_backoff_seconds <= 120:
        raise ValueError("max_backoff_seconds must be between backoff_seconds and 120")


def _request_json(
    url: str,
    *,
    transport: Transport,
    timeout_seconds: float,
    max_attempts: int,
    backoff_seconds: float,
    max_backoff_seconds: float,
    sleeper: SleepHook,
    before_request: Callable[[], None],
) -> tuple[dict[str, Any], int]:
    _validate_request_bounds(
        timeout_seconds=timeout_seconds,
        max_attempts=max_attempts,
        backoff_seconds=backoff_seconds,
        max_backoff_seconds=max_backoff_seconds,
    )
    for attempt in range(1, max_attempts + 1):
        before_request()
        try:
            response = transport(url, timeout_seconds)
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            if attempt >= max_attempts:
                raise EuropePmcRequestError(
                    f"Europe PMC transport failure: {exc}",
                    retryable=True,
                    code="transport_error",
                    attempts=attempt,
                ) from exc
            sleeper(_backoff_delay(attempt, backoff_seconds, max_backoff_seconds))
            continue

        if response.status != 200:
            retryable = response.status in RETRYABLE_HTTP_STATUSES
            retry_after = _retry_after_seconds(_header(response.headers, "Retry-After"))
            if not retryable:
                raise EuropePmcRequestError(
                    f"Europe PMC returned terminal HTTP {response.status}",
                    retryable=False,
                    code=f"http_{response.status}",
                    attempts=attempt,
                    http_status=response.status,
                    retry_after_seconds=retry_after,
                )
            if attempt >= max_attempts or (
                retry_after is not None and retry_after > max_backoff_seconds
            ):
                raise EuropePmcRequestError(
                    f"Europe PMC returned retryable HTTP {response.status}",
                    retryable=True,
                    code=f"http_{response.status}",
                    attempts=attempt,
                    http_status=response.status,
                    retry_after_seconds=retry_after,
                )
            delay = (
                retry_after
                if retry_after is not None
                else _backoff_delay(attempt, backoff_seconds, max_backoff_seconds)
            )
            sleeper(delay)
            continue

        try:
            payload = json.loads(response.body.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            if attempt >= max_attempts:
                raise EuropePmcRequestError(
                    f"Europe PMC returned invalid JSON: {exc}",
                    retryable=True,
                    code="invalid_json",
                    attempts=attempt,
                    http_status=response.status,
                ) from exc
            sleeper(_backoff_delay(attempt, backoff_seconds, max_backoff_seconds))
            continue
        if not isinstance(payload, dict):
            raise EuropePmcRequestError(
                "Europe PMC JSON response is not an object",
                retryable=False,
                code="invalid_response_schema",
                attempts=attempt,
                http_status=response.status,
            )
        return payload, attempt
    raise RuntimeError("unreachable")


def _utc_timestamp() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _error_row(
    species: str,
    *,
    query: str,
    search_api_url: str,
    error_class: str,
    error_code: str,
    message: str,
    attempts: int,
    retrieved_at_utc: str,
    http_status: int | None = None,
    retry_after_seconds: float | None = None,
) -> dict[str, str]:
    return {
        "accepted_species": species,
        "query": query,
        "search_api_url": search_api_url,
        "error_class": error_class,
        "error_code": error_code,
        "http_status": str(http_status) if http_status is not None else "",
        "message": _text(message),
        "attempts": str(attempts),
        "retry_after_seconds": (
            f"{retry_after_seconds:g}" if retry_after_seconds is not None else ""
        ),
        "retrieved_at_utc": retrieved_at_utc,
    }


def _contains_full_species(text: str, species: str) -> bool:
    pattern = re.compile(
        rf"(?<![\w-]){re.escape(species)}(?![\w-])",
        flags=re.IGNORECASE,
    )
    return bool(pattern.search(text))


def _contains_reproductive_term(text: str) -> bool:
    return any(pattern.search(text) for pattern in REPRODUCTIVE_TEXT_PATTERNS)


def _join_source_text(title: str, abstract: str) -> str:
    if title and abstract:
        separator = " " if title.endswith((".", "?", "!")) else ". "
        return f"{title}{separator}{abstract}"
    return title or abstract


def _clean_doi(value: object) -> str:
    doi = _text(value)
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if doi.casefold().startswith(prefix):
            doi = doi[len(prefix) :]
            break
    return doi.strip().casefold()


def _nested_journal_title(record: Mapping[str, Any]) -> str:
    journal_info = record.get("journalInfo")
    if not isinstance(journal_info, Mapping):
        return ""
    journal = journal_info.get("journal")
    if not isinstance(journal, Mapping):
        return ""
    return _plain_text(journal.get("title"))


def _license(record: Mapping[str, Any]) -> str:
    direct = _text(record.get("license"))
    if direct:
        return direct
    url_list = record.get("fullTextUrlList")
    if not isinstance(url_list, Mapping):
        return ""
    entries = url_list.get("fullTextUrl")
    if not isinstance(entries, list):
        return ""
    for entry in entries:
        if isinstance(entry, Mapping) and _text(entry.get("license")):
            return _text(entry.get("license"))
    return ""


def _sentence(value: str) -> str:
    value = value.strip()
    if not value:
        return ""
    return value if value.endswith((".", "?", "!")) else f"{value}."


def _citation(
    *,
    authors: str,
    year: str,
    title: str,
    journal: str,
    doi: str,
    pmid: str,
    pmcid: str,
) -> str:
    parts = [_sentence(authors)] if authors else []
    if year:
        parts.append(f"({year}).")
    if title:
        parts.append(_sentence(title))
    if journal:
        parts.append(_sentence(journal))
    identifiers = []
    if doi:
        identifiers.append(f"DOI:{doi}")
    if pmid:
        identifiers.append(f"PMID:{pmid}")
    if pmcid:
        identifiers.append(f"PMCID:{pmcid}")
    if identifiers:
        parts.append("; ".join(identifiers))
    return _text(" ".join(part for part in parts if part))


def _source_record(
    species: str,
    *,
    query: str,
    search_api_url: str,
    result: Mapping[str, Any],
    rank: int,
    retrieved_at_utc: str,
) -> dict[str, str] | None:
    title = _plain_text(result.get("title"))
    abstract = _plain_text(result.get("abstractText"))
    source_text = _join_source_text(title, abstract)
    if not source_text:
        return None
    if not _contains_full_species(source_text, species):
        return None
    if not _contains_reproductive_term(source_text):
        return None

    article_source = _text(result.get("source")).upper()
    article_id = _text(result.get("id"))
    if not article_source or not article_id:
        return None
    source_path = urllib.parse.quote(article_source, safe="")
    id_path = urllib.parse.quote(article_id, safe="")
    article_url = f"{EUROPE_PMC_ARTICLE_ROOT}/{source_path}/{id_path}"
    article_api_url = (
        f"{EUROPE_PMC_REST_ROOT}/article/{source_path}/{id_path}?resultType=core&format=json"
    )

    pmid = _text(result.get("pmid"))
    if not pmid and article_source == "MED" and article_id.isdigit():
        pmid = article_id
    pmcid = _text(result.get("pmcid")).upper()
    doi = _clean_doi(result.get("doi"))
    authors = _plain_text(result.get("authorString"))
    journal = _nested_journal_title(result)
    year = _text(result.get("pubYear"))
    first_publication_date = _text(result.get("firstPublicationDate"))
    license_value = _license(result)

    row = {
        "accepted_species": species,
        "source_text": source_text,
        "title": title,
        "abstract": abstract,
        "source_url": article_url,
        "article_url": article_url,
        "article_api_url": article_api_url,
        "search_api_url": search_api_url,
        "source_citation": _citation(
            authors=authors,
            year=year,
            title=title,
            journal=journal,
            doi=doi,
            pmid=pmid,
            pmcid=pmcid,
        ),
        "source_type": "europe_pmc_title_abstract",
        "evidence_scope": "species_direct",
        "pmid": pmid,
        "pmcid": pmcid,
        "doi": doi,
        "europe_pmc_source": article_source,
        "europe_pmc_id": article_id,
        "author_string": authors,
        "journal_title": journal,
        "publication_year": year,
        "first_publication_date": first_publication_date,
        "license": license_value,
        "is_open_access": _text(result.get("isOpenAccess")),
        "query": query,
        "search_result_rank": str(rank),
        "retrieved_at_utc": retrieved_at_utc,
        "source_text_sha256": _sha256_text(source_text),
        "metadata_sha256": "",
    }
    stable_metadata = {
        key: value
        for key, value in row.items()
        if key not in {"retrieved_at_utc", "metadata_sha256"}
    }
    row["metadata_sha256"] = _canonical_json_sha256(stable_metadata)
    return row


def discover_species_sources(
    species: str,
    *,
    page_size: int = 25,
    transport: Transport = _stdlib_transport,
    timeout_seconds: float = 30,
    max_attempts: int = 3,
    backoff_seconds: float = 1,
    max_backoff_seconds: float = 30,
    sleeper: SleepHook = time.sleep,
    before_request: Callable[[], None] = lambda: None,
    retrieved_at_utc: str | None = None,
) -> DiscoveryBatch:
    """Return species-direct Europe PMC sources or a classified accounted error."""

    retrieved = retrieved_at_utc or _utc_timestamp()
    raw_species = _text(species)
    try:
        focal = _normalise_species(raw_species)
        query = build_query(focal)
        search_api_url = build_search_api_url(focal, page_size=page_size)
    except ValueError as exc:
        error = _error_row(
            raw_species,
            query="",
            search_api_url="",
            error_class="terminal",
            error_code="invalid_species_name",
            message=str(exc),
            attempts=0,
            retrieved_at_utc=retrieved,
        )
        return DiscoveryBatch((raw_species,), (), (error,))

    try:
        payload, attempts = _request_json(
            search_api_url,
            transport=transport,
            timeout_seconds=timeout_seconds,
            max_attempts=max_attempts,
            backoff_seconds=backoff_seconds,
            max_backoff_seconds=max_backoff_seconds,
            sleeper=sleeper,
            before_request=before_request,
        )
    except EuropePmcRequestError as exc:
        error = _error_row(
            focal,
            query=query,
            search_api_url=search_api_url,
            error_class="retryable" if exc.retryable else "terminal",
            error_code=exc.code,
            message=str(exc),
            attempts=exc.attempts,
            retrieved_at_utc=retrieved,
            http_status=exc.http_status,
            retry_after_seconds=exc.retry_after_seconds,
        )
        return DiscoveryBatch((focal,), (), (error,))

    result_list = payload.get("resultList")
    if not isinstance(result_list, Mapping):
        error = _error_row(
            focal,
            query=query,
            search_api_url=search_api_url,
            error_class="terminal",
            error_code="invalid_response_schema",
            message="Europe PMC response is missing resultList",
            attempts=attempts,
            retrieved_at_utc=retrieved,
            http_status=200,
        )
        return DiscoveryBatch((focal,), (), (error,))
    results = result_list.get("result", [])
    if isinstance(results, Mapping):
        results = [results]
    if not isinstance(results, list) or not all(isinstance(item, Mapping) for item in results):
        error = _error_row(
            focal,
            query=query,
            search_api_url=search_api_url,
            error_class="terminal",
            error_code="invalid_response_schema",
            message="Europe PMC resultList.result is not a list of objects",
            attempts=attempts,
            retrieved_at_utc=retrieved,
            http_status=200,
        )
        return DiscoveryBatch((focal,), (), (error,))

    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for rank, result in enumerate(results, start=1):
        row = _source_record(
            focal,
            query=query,
            search_api_url=search_api_url,
            result=result,
            rank=rank,
            retrieved_at_utc=retrieved,
        )
        if row is None:
            continue
        identity = (row["europe_pmc_source"], row["europe_pmc_id"])
        if identity in seen:
            continue
        seen.add(identity)
        rows.append(row)

    if rows:
        return DiscoveryBatch((focal,), tuple(rows), ())
    error_code = "no_results" if not results else "no_species_direct_reproductive_source"
    message = (
        "Europe PMC returned no matching records"
        if not results
        else "No result contained the full species name and a reproductive term in title/abstract"
    )
    error = _error_row(
        focal,
        query=query,
        search_api_url=search_api_url,
        error_class="terminal",
        error_code=error_code,
        message=message,
        attempts=attempts,
        retrieved_at_utc=retrieved,
        http_status=200,
    )
    return DiscoveryBatch((focal,), (), (error,))


def discover_europe_pmc_sources(
    species_names: Iterable[str],
    *,
    page_size: int = 25,
    transport: Transport = _stdlib_transport,
    timeout_seconds: float = 30,
    max_attempts: int = 3,
    backoff_seconds: float = 1,
    max_backoff_seconds: float = 30,
    min_interval_seconds: float = 0.25,
    sleeper: SleepHook = time.sleep,
    clock: ClockHook = time.monotonic,
    retrieved_at_utc: str | None = None,
) -> DiscoveryBatch:
    """Discover sources serially with a bounded, shared request-rate limiter."""

    requested = tuple(dict.fromkeys(_text(species) for species in species_names if _text(species)))
    retrieved = retrieved_at_utc or _utc_timestamp()
    limiter = BoundedRateLimiter(
        min_interval_seconds=min_interval_seconds,
        sleeper=sleeper,
        clock=clock,
    )
    sources: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    for species in requested:
        batch = discover_species_sources(
            species,
            page_size=page_size,
            transport=transport,
            timeout_seconds=timeout_seconds,
            max_attempts=max_attempts,
            backoff_seconds=backoff_seconds,
            max_backoff_seconds=max_backoff_seconds,
            sleeper=sleeper,
            before_request=limiter.wait,
            retrieved_at_utc=retrieved,
        )
        sources.extend(batch.sources)
        errors.extend(batch.errors)
    return DiscoveryBatch(requested, tuple(sources), tuple(errors))


def _write_csv(
    path: Path,
    columns: Sequence[str],
    rows: Sequence[Mapping[str, str]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def freeze_europe_pmc_sources(
    species_names: Iterable[str],
    output_dir: Path,
    **discovery_kwargs: Any,
) -> dict[str, Any]:
    """Freeze returned metadata, errors, hashes, and a no-inference manifest."""

    batch = discover_europe_pmc_sources(species_names, **discovery_kwargs)
    output_dir.mkdir(parents=True, exist_ok=True)
    source_path = output_dir / "europe_pmc_source_texts.csv"
    error_path = output_dir / "europe_pmc_discovery_errors.csv"
    _write_csv(source_path, SOURCE_COLUMNS, batch.sources)
    _write_csv(error_path, ERROR_COLUMNS, batch.errors)

    retrieved = ""
    if batch.sources:
        retrieved = batch.sources[0]["retrieved_at_utc"]
    elif batch.errors:
        retrieved = batch.errors[0]["retrieved_at_utc"]
    manifest = {
        "contract_version": CONTRACT_VERSION,
        "source": "Europe PMC REST API",
        "api_root": EUROPE_PMC_REST_ROOT,
        "retrieved_at_utc": retrieved,
        "n_requested_species": len(batch.requested_species),
        "n_source_records": len(batch.sources),
        "n_species_with_sources": len({row["accepted_species"] for row in batch.sources}),
        "n_errors": len(batch.errors),
        "error_class_counts": {
            kind: sum(row["error_class"] == kind for row in batch.errors)
            for kind in ("retryable", "terminal")
        },
        "query_terms": list(REPRODUCTIVE_QUERY_TERMS),
        "files": {
            source_path.name: {
                "sha256": _file_sha256(source_path),
                "rows": len(batch.sources),
            },
            error_path.name: {
                "sha256": _file_sha256(error_path),
                "rows": len(batch.errors),
            },
        },
        "policy": (
            "Europe PMC REST core JSON only; no article-page scraping; retain only full "
            "species-name plus reproductive-term title/abstract records; discovery assigns "
            "and infers no plant trait."
        ),
    }
    (output_dir / "europe_pmc_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest


def build_full_text_xml_url(pmcid: object) -> str:
    """Return the official Europe PMC OA full-text endpoint for one PMCID."""

    value = _text(pmcid).upper()
    if not re.fullmatch(r"PMC\d+", value):
        raise ValueError("pmcid must match PMC followed by digits")
    return f"{EUROPE_PMC_REST_ROOT}/{value}/fullTextXML"


def _request_bytes(
    url: str,
    *,
    transport: Transport,
    timeout_seconds: float,
    max_attempts: int,
    backoff_seconds: float,
    max_backoff_seconds: float,
    sleeper: SleepHook,
    before_request: Callable[[], None],
) -> tuple[bytes, int]:
    _validate_request_bounds(
        timeout_seconds=timeout_seconds,
        max_attempts=max_attempts,
        backoff_seconds=backoff_seconds,
        max_backoff_seconds=max_backoff_seconds,
    )
    for attempt in range(1, max_attempts + 1):
        before_request()
        try:
            response = transport(url, timeout_seconds)
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            if attempt >= max_attempts:
                raise EuropePmcRequestError(
                    f"Europe PMC transport failure: {exc}",
                    retryable=True,
                    code="transport_error",
                    attempts=attempt,
                ) from exc
            sleeper(_backoff_delay(attempt, backoff_seconds, max_backoff_seconds))
            continue

        if response.status == 200:
            if len(response.body) > MAX_FULL_TEXT_BYTES:
                raise EuropePmcRequestError(
                    "Europe PMC full-text XML exceeded 25 MB",
                    retryable=False,
                    code="full_text_too_large",
                    attempts=attempt,
                    http_status=response.status,
                )
            return response.body, attempt

        retryable = response.status in RETRYABLE_HTTP_STATUSES
        retry_after = _retry_after_seconds(_header(response.headers, "Retry-After"))
        if not retryable:
            raise EuropePmcRequestError(
                f"Europe PMC returned terminal HTTP {response.status}",
                retryable=False,
                code=f"http_{response.status}",
                attempts=attempt,
                http_status=response.status,
                retry_after_seconds=retry_after,
            )
        if attempt >= max_attempts or (
            retry_after is not None and retry_after > max_backoff_seconds
        ):
            raise EuropePmcRequestError(
                f"Europe PMC returned retryable HTTP {response.status}",
                retryable=True,
                code=f"http_{response.status}",
                attempts=attempt,
                http_status=response.status,
                retry_after_seconds=retry_after,
            )
        sleeper(
            retry_after
            if retry_after is not None
            else _backoff_delay(attempt, backoff_seconds, max_backoff_seconds)
        )
    raise RuntimeError("unreachable")


def _element_name(element: ElementTree.Element) -> str:
    return str(element.tag).rsplit("}", 1)[-1].casefold()


def _sentence_candidates(text: str) -> list[str]:
    normalized = _text(text)
    if not normalized:
        return []
    return [
        _text(part)
        for part in re.split(r"(?<=[.!?])\s+(?=[A-Z0-9(])", normalized)
        if _text(part)
    ]


def extract_species_direct_text_passages(
    text: str,
    accepted_species: str,
    *,
    aliases: Iterable[str] = (),
    max_passages: int = 50,
) -> tuple[dict[str, str], ...]:
    """Keep text sentences joining an exact focal name to a direct trait term."""

    species = _normalise_species(accepted_species)
    names = tuple(
        dict.fromkeys(
            _normalise_species(name)
            for name in (species, *aliases)
            if _text(name)
        )
    )
    rows: list[dict[str, str]] = []
    for unit in _sentence_candidates(text):
        matched_name = next(
            (name for name in names if _contains_full_species(unit, name)),
            "",
        )
        if not matched_name or not any(
            pattern.search(unit) for pattern in FULL_TEXT_DIRECT_PATTERNS
        ):
            continue
        rows.append(
            {
                "text": unit,
                "matched_name": matched_name,
                "element": "sentence",
                "passage_sha256": _sha256_text(unit),
            }
        )
        if len(rows) >= max_passages:
            break
    return tuple(rows)


def extract_species_direct_full_text_passages(
    xml_bytes: bytes,
    accepted_species: str,
    *,
    aliases: Iterable[str] = (),
    max_passages: int = 50,
) -> tuple[dict[str, str], ...]:
    """Keep only XML sentences/table rows joining an exact name to a direct trait term."""

    species = _normalise_species(accepted_species)
    names = tuple(
        dict.fromkeys(
            _normalise_species(name)
            for name in (species, *aliases)
            if _text(name)
        )
    )
    if not 1 <= max_passages <= 500:
        raise ValueError("max_passages must be between 1 and 500")
    if len(xml_bytes) > MAX_FULL_TEXT_BYTES:
        raise ValueError("Europe PMC full-text XML exceeded 25 MB")
    try:
        root = ElementTree.fromstring(xml_bytes)
    except ElementTree.ParseError as exc:
        raise ValueError(f"invalid Europe PMC full-text XML: {exc}") from exc

    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    wanted_elements = {"article-title", "title", "p", "caption", "tr"}
    for element in root.iter():
        element_name = _element_name(element)
        if element_name not in wanted_elements:
            continue
        text = _text(" ".join(element.itertext()))
        if not text:
            continue
        units = [text] if element_name == "tr" else _sentence_candidates(text)
        for unit in units:
            matched_name = next(
                (name for name in names if _contains_full_species(unit, name)),
                "",
            )
            if not matched_name or not any(pattern.search(unit) for pattern in FULL_TEXT_DIRECT_PATTERNS):
                continue
            passage_sha = _sha256_text(unit)
            if passage_sha in seen:
                continue
            seen.add(passage_sha)
            rows.append(
                {
                    "text": unit,
                    "matched_name": matched_name,
                    "element": element_name,
                    "passage_sha256": passage_sha,
                }
            )
            if len(rows) >= max_passages:
                return tuple(rows)
    return tuple(rows)


def _full_text_error_row(
    *,
    species: str,
    pmcid: str,
    url: str,
    error_class: str,
    error_code: str,
    message: str,
    attempts: int,
    retrieved_at_utc: str,
    http_status: int | None = None,
    retry_after_seconds: float | None = None,
) -> dict[str, str]:
    return {
        "accepted_species": species,
        "pmcid": pmcid,
        "full_text_xml_url": url,
        "error_class": error_class,
        "error_code": error_code,
        "http_status": str(http_status) if http_status is not None else "",
        "message": _text(message),
        "attempts": str(attempts),
        "retry_after_seconds": (
            f"{retry_after_seconds:g}" if retry_after_seconds is not None else ""
        ),
        "retrieved_at_utc": retrieved_at_utc,
    }


def discover_europe_pmc_full_text_sources(
    discovery_sources: Iterable[Mapping[str, str]],
    *,
    aliases_by_species: Mapping[str, Iterable[str]] | None = None,
    transport: Transport = _stdlib_transport,
    timeout_seconds: float = 45,
    max_attempts: int = 3,
    backoff_seconds: float = 1,
    max_backoff_seconds: float = 30,
    min_interval_seconds: float = 0.25,
    sleeper: SleepHook = time.sleep,
    clock: ClockHook = time.monotonic,
    retrieved_at_utc: str | None = None,
) -> FullTextBatch:
    """Retrieve OA XML and retain only exact-name, direct reproductive passages."""

    retrieved = retrieved_at_utc or _utc_timestamp()
    aliases_map = aliases_by_species or {}
    limiter = BoundedRateLimiter(
        min_interval_seconds=min_interval_seconds,
        sleeper=sleeper,
        clock=clock,
    )
    records: list[tuple[str, str, Mapping[str, str]]] = []
    seen_articles: set[tuple[str, str]] = set()
    errors: list[dict[str, str]] = []
    sources: list[dict[str, str]] = []

    for record in discovery_sources:
        species = _text(record.get("accepted_species"))
        pmcid = _text(record.get("pmcid")).upper()
        if not species:
            continue
        if not pmcid:
            errors.append(
                _full_text_error_row(
                    species=species,
                    pmcid="",
                    url="",
                    error_class="terminal",
                    error_code="missing_pmcid",
                    message="Europe PMC discovery record has no PMCID",
                    attempts=0,
                    retrieved_at_utc=retrieved,
                )
            )
            continue
        key = (species, pmcid)
        if key in seen_articles:
            continue
        seen_articles.add(key)
        try:
            url = build_full_text_xml_url(pmcid)
        except ValueError as exc:
            errors.append(
                _full_text_error_row(
                    species=species,
                    pmcid=pmcid,
                    url="",
                    error_class="terminal",
                    error_code="invalid_pmcid",
                    message=str(exc),
                    attempts=0,
                    retrieved_at_utc=retrieved,
                )
            )
            continue
        records.append((species, url, record))

    for species, url, record in records:
        pmcid = _text(record.get("pmcid")).upper()
        if _text(record.get("is_open_access")).casefold() not in {"y", "yes", "true", "1"}:
            errors.append(
                _full_text_error_row(
                    species=species,
                    pmcid=pmcid,
                    url=url,
                    error_class="terminal",
                    error_code="not_open_access",
                    message="Europe PMC record is not marked open access",
                    attempts=0,
                    retrieved_at_utc=retrieved,
                )
            )
            continue
        try:
            xml_bytes, attempts = _request_bytes(
                url,
                transport=transport,
                timeout_seconds=timeout_seconds,
                max_attempts=max_attempts,
                backoff_seconds=backoff_seconds,
                max_backoff_seconds=max_backoff_seconds,
                sleeper=sleeper,
                before_request=limiter.wait,
            )
            passages = extract_species_direct_full_text_passages(
                xml_bytes,
                species,
                aliases=aliases_map.get(species, ()),
            )
        except EuropePmcRequestError as exc:
            errors.append(
                _full_text_error_row(
                    species=species,
                    pmcid=pmcid,
                    url=url,
                    error_class="retryable" if exc.retryable else "terminal",
                    error_code=exc.code,
                    message=str(exc),
                    attempts=exc.attempts,
                    retrieved_at_utc=retrieved,
                    http_status=exc.http_status,
                    retry_after_seconds=exc.retry_after_seconds,
                )
            )
            continue
        except ValueError as exc:
            errors.append(
                _full_text_error_row(
                    species=species,
                    pmcid=pmcid,
                    url=url,
                    error_class="terminal",
                    error_code="invalid_full_text_xml",
                    message=str(exc),
                    attempts=1,
                    retrieved_at_utc=retrieved,
                    http_status=200,
                )
            )
            continue
        if not passages:
            errors.append(
                _full_text_error_row(
                    species=species,
                    pmcid=pmcid,
                    url=url,
                    error_class="terminal",
                    error_code="no_species_direct_reproductive_passage",
                    message=(
                        "No XML sentence or table row contained both an exact focal "
                        "name and a direct reproductive trait term"
                    ),
                    attempts=attempts,
                    retrieved_at_utc=retrieved,
                    http_status=200,
                )
            )
            continue

        xml_sha = _sha256_bytes(xml_bytes)
        doi = _clean_doi(record.get("doi"))
        lineage = f"doi:{doi}" if doi else f"pmcid:{pmcid.casefold()}"
        base_citation = _text(record.get("source_citation"))
        for passage in passages:
            passage_sha = passage["passage_sha256"]
            citation = " | ".join(
                part
                for part in (
                    base_citation,
                    f"Europe PMC OA full text {pmcid}",
                    f"full_text_xml_sha256={xml_sha}",
                    f"passage_sha256={passage_sha}",
                )
                if part
            )
            sources.append(
                {
                    "accepted_species": species,
                    "source_text": passage["text"],
                    "source_url": url,
                    "source_citation": citation,
                    "source_type": "europe_pmc_full_text_xml",
                    "evidence_scope": "species_direct",
                    "title": _plain_text(record.get("title")),
                    "pmcid": pmcid,
                    "doi": doi,
                    "license": _text(record.get("license")),
                    "is_open_access": _text(record.get("is_open_access")),
                    "full_text_xml_url": url,
                    "full_text_xml_sha256": xml_sha,
                    "passage_sha256": passage_sha,
                    "source_record_id": f"{pmcid}:{passage_sha}",
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
                    "retrieved_at_utc": retrieved,
                }
            )
    requested = tuple((species, _text(record.get("pmcid")).upper()) for species, _, record in records)
    return FullTextBatch(requested, tuple(sources), tuple(errors))


def freeze_europe_pmc_full_text_sources(
    discovery_sources: Iterable[Mapping[str, str]],
    output_dir: Path,
    **discovery_kwargs: Any,
) -> dict[str, Any]:
    """Freeze OA XML passages, accounted misses, hashes, and source-lineage policy."""

    batch = discover_europe_pmc_full_text_sources(discovery_sources, **discovery_kwargs)
    output_dir.mkdir(parents=True, exist_ok=True)
    source_path = output_dir / "europe_pmc_full_text_sources.csv"
    error_path = output_dir / "europe_pmc_full_text_errors.csv"
    _write_csv(source_path, FULL_TEXT_SOURCE_COLUMNS, batch.sources)
    _write_csv(error_path, FULL_TEXT_ERROR_COLUMNS, batch.errors)
    manifest = {
        "contract_version": FULL_TEXT_CONTRACT_VERSION,
        "source": "Europe PMC REST OA fullTextXML",
        "api_root": EUROPE_PMC_REST_ROOT,
        "n_requested_articles": len(batch.requested_articles),
        "n_source_passages": len(batch.sources),
        "n_species_with_passages": len(
            {row["accepted_species"] for row in batch.sources}
        ),
        "n_errors": len(batch.errors),
        "files": {
            source_path.name: {
                "sha256": _file_sha256(source_path),
                "rows": len(batch.sources),
            },
            error_path.name: {
                "sha256": _file_sha256(error_path),
                "rows": len(batch.errors),
            },
        },
        "policy": (
            "Official OA fullTextXML only; exact accepted species or explicitly supplied "
            "synonym must co-occur with a direct reproductive trait term in the same XML "
            "sentence or table row; no family inference, global fallback, or abbreviated-"
            "name expansion."
        ),
    }
    (output_dir / "europe_pmc_full_text_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest
