"""Official web-search discovery without treating snippets as trait evidence.

The adapter supports Google Custom Search JSON API for callers that already
have credentials and Brave Search API as an optional fallback.  Every HTTP
attempt, including retries, consumes one unit from an explicit hard query
budget.  Returned rows are discovery leads only: this module never emits a
trait value and never promotes a search snippet to evidence.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import html
import json
import os
import re
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from datetime import UTC, datetime
from email.utils import parsedate_to_datetime
from pathlib import Path
from typing import IO
from typing import Any

import typer

GOOGLE_CUSTOM_SEARCH_URL = "https://customsearch.googleapis.com/customsearch/v1"
BRAVE_SEARCH_URL = "https://api.search.brave.com/res/v1/web/search"
USER_AGENT = "island-floral-traits/0.1 (official search discovery)"
CONTRACT_VERSION = "official-search-discovery-v1"

GOOGLE_API_KEY_ENV = "GOOGLE_CUSTOM_SEARCH_API_KEY"
GOOGLE_CX_ENV = "GOOGLE_CUSTOM_SEARCH_CX"
BRAVE_API_KEY_ENV = "BRAVE_SEARCH_API_KEY"

DISCOVERY_FILENAME = "discovery.csv"
ERROR_FILENAME = "errors.csv"
REPORT_FILENAME = "report.json"

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Discover public search-result URLs without extracting or inferring plant traits.",
)

DISCOVERY_COLUMNS = [
    "species",
    "provider",
    "query",
    "rank",
    "result_url",
    "title",
    "snippet",
    "retrieved_at_utc",
    "discovery_sha256",
]

ERROR_COLUMNS = [
    "species",
    "provider",
    "query",
    "error_class",
    "error_code",
    "message",
]

QUERY_TERMS = (
    "flower",
    "floral",
    "pollination",
    "pollinator",
    '"breeding system"',
    '"mating system"',
    "selfing",
    "outcrossing",
    '"self-incompatibility"',
)

RETRYABLE_HTTP_STATUSES = {408, 425, 429, 500, 502, 503, 504}


@dataclass(frozen=True, repr=False)
class SearchRequest:
    """An HTTP request whose repr deliberately cannot expose credentials."""

    url: str
    headers: Mapping[str, str]


@dataclass(frozen=True)
class SearchHttpResponse:
    status: int
    headers: Mapping[str, str]
    body: bytes


Transport = Callable[[SearchRequest, float], SearchHttpResponse]
SleepHook = Callable[[float], None]
ClockHook = Callable[[], float]


@dataclass(frozen=True)
class SearchDiscoveryBatch:
    requested_species: tuple[str, ...]
    records: tuple[dict[str, str], ...]
    errors: tuple[dict[str, str], ...]
    queries_used: int
    max_queries: int


class SearchRequestError(RuntimeError):
    def __init__(
        self,
        message: str,
        *,
        retryable: bool,
        code: str,
    ) -> None:
        super().__init__(message)
        self.retryable = retryable
        self.code = code


@dataclass
class QueryBudget:
    """A deterministic hard cap over all real HTTP attempts."""

    max_queries: int
    used: int = field(default=0, init=False)

    def __post_init__(self) -> None:
        if not 0 <= self.max_queries <= 1_000_000:
            raise ValueError("max_queries must be between 0 and 1,000,000")

    @property
    def remaining(self) -> int:
        return self.max_queries - self.used

    def reserve(self) -> None:
        if self.remaining <= 0:
            raise SearchRequestError(
                "hard search query budget exhausted",
                retryable=False,
                code="budget_exhausted",
            )
        self.used += 1


@dataclass
class BoundedRateLimiter:
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


def _text(value: object) -> str:
    if value is None:
        return ""
    decoded = html.unescape(str(value)).replace("\xa0", " ")
    decoded = re.sub(r"<[^>]+>", " ", decoded)
    text = " ".join(decoded.split())
    return re.sub(r"\s+([,.;:!?])", r"\1", text)


def _normalise_species(species: object) -> str:
    value = _text(species)
    if len(value.split()) < 2:
        raise ValueError("species must contain at least a genus and specific epithet")
    if any(character in value for character in ('"', "\\", "\r", "\n")):
        raise ValueError("species contains an unsafe search-query character")
    return value


def build_exact_species_query(species: str) -> str:
    """Build the same exact-phrase reproductive discovery query for both APIs."""

    focal = _normalise_species(species)
    terms = " OR ".join(QUERY_TERMS)
    return f'"{focal}" ({terms})'


def _secret(value: str | None, label: str) -> str:
    cleaned = str(value or "").strip()
    if not cleaned:
        return ""
    if any(character in cleaned for character in ("\r", "\n", "\x00")):
        raise ValueError(f"{label} contains an unsafe character")
    return cleaned


def _provider_credentials(
    *,
    google_api_key: str | None,
    google_cx: str | None,
    brave_api_key: str | None,
) -> tuple[tuple[str, str, str], ...]:
    google_key = _secret(google_api_key, "google_api_key")
    cx = _secret(google_cx, "google_cx")
    brave_key = _secret(brave_api_key, "brave_api_key")
    if bool(google_key) != bool(cx):
        raise ValueError("Google Custom Search requires both google_api_key and google_cx")
    providers: list[tuple[str, str, str]] = []
    if google_key:
        providers.append(("google_custom_search", google_key, cx))
    if brave_key:
        providers.append(("brave_search", brave_key, ""))
    if not providers:
        raise ValueError("an explicit Google or Brave API key is required")
    return tuple(providers)


def _google_request(query: str, api_key: str, cx: str, count: int) -> SearchRequest:
    parameters = [
        ("key", api_key),
        ("cx", cx),
        ("q", query),
        ("num", str(count)),
        ("safe", "active"),
    ]
    return SearchRequest(
        url=f"{GOOGLE_CUSTOM_SEARCH_URL}?{urllib.parse.urlencode(parameters)}",
        headers={"Accept": "application/json", "User-Agent": USER_AGENT},
    )


def _brave_request(query: str, api_key: str, count: int) -> SearchRequest:
    parameters = [
        ("q", query),
        ("count", str(count)),
        ("safesearch", "moderate"),
    ]
    return SearchRequest(
        url=f"{BRAVE_SEARCH_URL}?{urllib.parse.urlencode(parameters)}",
        headers={
            "Accept": "application/json",
            "User-Agent": USER_AGENT,
            "X-Subscription-Token": api_key,
        },
    )


def _stdlib_transport(request: SearchRequest, timeout_seconds: float) -> SearchHttpResponse:
    http_request = urllib.request.Request(request.url, headers=dict(request.headers))
    try:
        with urllib.request.urlopen(http_request, timeout=timeout_seconds) as response:
            return SearchHttpResponse(
                status=int(response.status),
                headers={str(key): str(value) for key, value in response.headers.items()},
                body=response.read(),
            )
    except urllib.error.HTTPError as exc:
        return SearchHttpResponse(
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


def _retry_after_seconds(value: object) -> float | None:
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
    return max(0.0, (target - datetime.now(UTC)).total_seconds())


def _validate_bounds(
    *,
    results_per_query: int,
    timeout_seconds: float,
    max_attempts: int,
    backoff_seconds: float,
    max_backoff_seconds: float,
) -> None:
    if not 1 <= results_per_query <= 10:
        raise ValueError("results_per_query must be between 1 and 10")
    if not 1 <= timeout_seconds <= 120:
        raise ValueError("timeout_seconds must be between 1 and 120")
    if not 1 <= max_attempts <= 5:
        raise ValueError("max_attempts must be between 1 and 5")
    if not 0 <= backoff_seconds <= 30:
        raise ValueError("backoff_seconds must be between 0 and 30")
    if not backoff_seconds <= max_backoff_seconds <= 120:
        raise ValueError("max_backoff_seconds must be between backoff_seconds and 120")


def _backoff(attempt: int, base_seconds: float, maximum_seconds: float) -> float:
    return min(maximum_seconds, base_seconds * (2 ** max(0, attempt - 1)))


def _request_json(
    provider: str,
    request: SearchRequest,
    *,
    budget: QueryBudget,
    limiter: BoundedRateLimiter,
    transport: Transport,
    timeout_seconds: float,
    max_attempts: int,
    backoff_seconds: float,
    max_backoff_seconds: float,
    sleeper: SleepHook,
) -> dict[str, Any]:
    for attempt in range(1, max_attempts + 1):
        budget.reserve()
        limiter.wait()
        try:
            response = transport(request, timeout_seconds)
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            if attempt >= max_attempts:
                raise SearchRequestError(
                    f"{provider} transport failure",
                    retryable=True,
                    code="transport_error",
                ) from exc
            if budget.remaining <= 0:
                raise SearchRequestError(
                    "hard search query budget exhausted",
                    retryable=False,
                    code="budget_exhausted",
                ) from exc
            sleeper(_backoff(attempt, backoff_seconds, max_backoff_seconds))
            continue

        if response.status != 200:
            retryable = response.status in RETRYABLE_HTTP_STATUSES
            if not retryable:
                raise SearchRequestError(
                    f"{provider} returned terminal HTTP {response.status}",
                    retryable=False,
                    code=f"http_{response.status}",
                )
            retry_after = _retry_after_seconds(_header(response.headers, "Retry-After"))
            if attempt >= max_attempts:
                raise SearchRequestError(
                    f"{provider} returned retryable HTTP {response.status}",
                    retryable=True,
                    code=f"http_{response.status}",
                )
            if budget.remaining <= 0:
                raise SearchRequestError(
                    "hard search query budget exhausted",
                    retryable=False,
                    code="budget_exhausted",
                )
            if retry_after is not None and retry_after > max_backoff_seconds:
                raise SearchRequestError(
                    f"{provider} Retry-After exceeds the bounded backoff",
                    retryable=True,
                    code=f"http_{response.status}",
                )
            sleeper(
                retry_after
                if retry_after is not None
                else _backoff(attempt, backoff_seconds, max_backoff_seconds)
            )
            continue

        try:
            payload = json.loads(response.body.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            if attempt >= max_attempts:
                raise SearchRequestError(
                    f"{provider} returned invalid JSON",
                    retryable=True,
                    code="invalid_json",
                ) from exc
            if budget.remaining <= 0:
                raise SearchRequestError(
                    "hard search query budget exhausted",
                    retryable=False,
                    code="budget_exhausted",
                ) from exc
            sleeper(_backoff(attempt, backoff_seconds, max_backoff_seconds))
            continue
        if not isinstance(payload, dict):
            raise SearchRequestError(
                f"{provider} JSON response is not an object",
                retryable=False,
                code="invalid_response_schema",
            )
        return payload
    raise RuntimeError("unreachable")


def _clean_result_url(value: object) -> str:
    raw = _text(value)
    if not raw:
        return ""
    parsed = urllib.parse.urlsplit(raw)
    if parsed.scheme.casefold() not in {"http", "https"} or not parsed.netloc:
        return ""
    return urllib.parse.urlunsplit(
        (parsed.scheme.casefold(), parsed.netloc, parsed.path, parsed.query, "")
    )


def _raw_results(provider: str, payload: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    if provider == "google_custom_search":
        items = payload.get("items", [])
    else:
        web = payload.get("web")
        items = web.get("results", []) if isinstance(web, Mapping) else []
    if items is None:
        return []
    if not isinstance(items, list):
        raise SearchRequestError(
            f"{provider} result list has an invalid schema",
            retryable=False,
            code="invalid_response_schema",
        )
    return [item for item in items if isinstance(item, Mapping)]


def _record(
    *,
    species: str,
    provider: str,
    query: str,
    rank: int,
    raw: Mapping[str, Any],
    retrieved_at_utc: str,
) -> dict[str, str] | None:
    url_key = "link" if provider == "google_custom_search" else "url"
    snippet_key = "snippet" if provider == "google_custom_search" else "description"
    result_url = _clean_result_url(raw.get(url_key))
    if not result_url:
        return None
    row = {
        "species": species,
        "provider": provider,
        "query": query,
        "rank": str(rank),
        "result_url": result_url,
        "title": _text(raw.get("title")),
        "snippet": _text(raw.get(snippet_key)),
        "retrieved_at_utc": retrieved_at_utc,
        "discovery_sha256": "",
    }
    stable = {
        key: value
        for key, value in row.items()
        if key not in {"retrieved_at_utc", "discovery_sha256"}
    }
    encoded = json.dumps(
        stable,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    row["discovery_sha256"] = hashlib.sha256(encoded).hexdigest()
    return row


def _error(
    species: str,
    provider: str,
    query: str,
    error_class: str,
    error_code: str,
    message: str,
) -> dict[str, str]:
    return {
        "species": species,
        "provider": provider,
        "query": query,
        "error_class": error_class,
        "error_code": error_code,
        "message": _text(message),
    }


def _utc_timestamp() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def discover_official_search_results(
    species_names: Iterable[str],
    *,
    max_queries: int,
    google_api_key: str | None = None,
    google_cx: str | None = None,
    brave_api_key: str | None = None,
    results_per_query: int = 10,
    timeout_seconds: float = 30,
    max_attempts: int = 3,
    backoff_seconds: float = 1,
    max_backoff_seconds: float = 30,
    min_interval_seconds: float = 0.25,
    transport: Transport = _stdlib_transport,
    sleeper: SleepHook = time.sleep,
    clock: ClockHook = time.monotonic,
    retrieved_at_utc: str | None = None,
) -> SearchDiscoveryBatch:
    """Return official-search URL leads under an exact hard request budget."""

    providers = _provider_credentials(
        google_api_key=google_api_key,
        google_cx=google_cx,
        brave_api_key=brave_api_key,
    )
    _validate_bounds(
        results_per_query=results_per_query,
        timeout_seconds=timeout_seconds,
        max_attempts=max_attempts,
        backoff_seconds=backoff_seconds,
        max_backoff_seconds=max_backoff_seconds,
    )
    requested = tuple(dict.fromkeys(_text(species) for species in species_names if _text(species)))
    retrieved = retrieved_at_utc or _utc_timestamp()
    budget = QueryBudget(max_queries)
    limiter = BoundedRateLimiter(
        min_interval_seconds=min_interval_seconds,
        sleeper=sleeper,
        clock=clock,
    )
    records: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []

    for raw_species in requested:
        try:
            species = _normalise_species(raw_species)
            query = build_exact_species_query(species)
        except ValueError as exc:
            errors.append(
                _error(
                    raw_species,
                    "",
                    "",
                    "terminal",
                    "invalid_species_name",
                    str(exc),
                )
            )
            continue
        if budget.remaining <= 0:
            errors.append(
                _error(
                    species,
                    "",
                    query,
                    "terminal",
                    "budget_exhausted",
                    "hard search query budget exhausted before species lookup",
                )
            )
            continue

        for provider, api_key, cx in providers:
            request = (
                _google_request(query, api_key, cx, results_per_query)
                if provider == "google_custom_search"
                else _brave_request(query, api_key, results_per_query)
            )
            try:
                payload = _request_json(
                    provider,
                    request,
                    budget=budget,
                    limiter=limiter,
                    transport=transport,
                    timeout_seconds=timeout_seconds,
                    max_attempts=max_attempts,
                    backoff_seconds=backoff_seconds,
                    max_backoff_seconds=max_backoff_seconds,
                    sleeper=sleeper,
                )
                raw_results = _raw_results(provider, payload)
            except SearchRequestError as exc:
                errors.append(
                    _error(
                        species,
                        provider,
                        query,
                        "retryable" if exc.retryable else "terminal",
                        exc.code,
                        str(exc),
                    )
                )
                if exc.code == "budget_exhausted":
                    break
                continue

            provider_records = [
                record
                for rank, raw in enumerate(raw_results, start=1)
                if (
                    record := _record(
                        species=species,
                        provider=provider,
                        query=query,
                        rank=rank,
                        raw=raw,
                        retrieved_at_utc=retrieved,
                    )
                )
                is not None
            ]
            if provider_records:
                records.extend(provider_records)
                break
            errors.append(
                _error(
                    species,
                    provider,
                    query,
                    "terminal",
                    "no_results" if not raw_results else "no_valid_result_urls",
                    "search returned no valid discovery URLs",
                )
            )
            if budget.remaining <= 0:
                break

    return SearchDiscoveryBatch(
        requested_species=requested,
        records=tuple(records),
        errors=tuple(errors),
        queries_used=budget.used,
        max_queries=budget.max_queries,
    )


def read_species_csv(path: Path) -> tuple[str, ...]:
    """Read nonblank values from an exact ``species`` column in CSV or CSV.gz."""

    source = Path(path)
    if source.suffix.casefold() == ".gz":
        handle = gzip.open(source, mode="rt", encoding="utf-8-sig", newline="")
    else:
        handle = source.open(mode="r", encoding="utf-8-sig", newline="")
    with handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "species" not in reader.fieldnames:
            raise ValueError("species CSV must contain a column named 'species'")
        species = tuple(value for row in reader if (value := _text(row.get("species"))))
    if not species:
        raise ValueError("species CSV contains no nonblank species values")
    return species


def _atomic_text_write(path: Path, writer: Callable[[IO[str]], None]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            prefix=f".{path.name}.",
            suffix=".tmp",
            dir=path.parent,
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            writer(handle)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def _atomic_csv(
    path: Path,
    columns: list[str],
    rows: Iterable[Mapping[str, str]],
) -> None:
    def write(handle: IO[str]) -> None:
        csv_writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            extrasaction="raise",
            lineterminator="\n",
        )
        csv_writer.writeheader()
        csv_writer.writerows(rows)

    _atomic_text_write(path, write)


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    def write(handle: IO[str]) -> None:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")

    _atomic_text_write(path, write)


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def freeze_official_search_results(
    species_names: Iterable[str],
    output_dir: Path,
    **discovery_kwargs: Any,
) -> dict[str, Any]:
    """Run discovery and atomically freeze discovery-only CSVs and a report."""

    batch = discover_official_search_results(species_names, **discovery_kwargs)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    discovery_path = destination / DISCOVERY_FILENAME
    error_path = destination / ERROR_FILENAME
    report_path = destination / REPORT_FILENAME

    _atomic_csv(discovery_path, DISCOVERY_COLUMNS, batch.records)
    _atomic_csv(error_path, ERROR_COLUMNS, batch.errors)
    report: dict[str, Any] = {
        "contract_version": CONTRACT_VERSION,
        "completed_at_utc": _utc_timestamp(),
        "n_requested_species": len(batch.requested_species),
        "n_species_with_results": len({row["species"] for row in batch.records}),
        "n_discovery_records": len(batch.records),
        "n_errors": len(batch.errors),
        "queries_used": batch.queries_used,
        "max_queries": batch.max_queries,
        "budget_remaining": batch.max_queries - batch.queries_used,
        "files": {
            discovery_path.name: {
                "rows": len(batch.records),
                "sha256": _file_sha256(discovery_path),
            },
            error_path.name: {
                "rows": len(batch.errors),
                "sha256": _file_sha256(error_path),
            },
        },
        "policy": {
            "artifact_scope": "search_result_discovery_only",
            "trait_values_emitted": False,
            "snippets_are_trait_evidence": False,
        },
    }
    _atomic_json(report_path, report)
    return report


@app.command()
def main(
    species_csv: Path = typer.Option(
        ...,
        "--species-csv",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Input CSV or CSV.gz containing an exact species column.",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        file_okay=False,
        dir_okay=True,
        help="Directory for discovery.csv, errors.csv, and report.json.",
    ),
    max_queries: int = typer.Option(
        ...,
        "--max-queries",
        min=1,
        max=1_000_000,
        help="Required hard cap over every HTTP attempt, including retries.",
    ),
    google_api_key: str | None = typer.Option(
        None,
        "--google-api-key",
        envvar=GOOGLE_API_KEY_ENV,
        hide_input=True,
        show_default=False,
        show_envvar=True,
        help="Google Custom Search API key; prefer the named environment variable.",
    ),
    google_cx: str | None = typer.Option(
        None,
        "--google-cx",
        envvar=GOOGLE_CX_ENV,
        hide_input=True,
        show_default=False,
        show_envvar=True,
        help="Google Programmable Search engine ID; prefer the named environment variable.",
    ),
    brave_api_key: str | None = typer.Option(
        None,
        "--brave-api-key",
        envvar=BRAVE_API_KEY_ENV,
        hide_input=True,
        show_default=False,
        show_envvar=True,
        help="Brave Search API key; prefer the named environment variable.",
    ),
    results_per_query: int = typer.Option(10, min=1, max=10),
    timeout_seconds: float = typer.Option(30, min=1, max=120),
    max_attempts: int = typer.Option(3, min=1, max=5),
    backoff_seconds: float = typer.Option(1, min=0, max=30),
    max_backoff_seconds: float = typer.Option(30, min=0, max=120),
    min_interval_seconds: float = typer.Option(0.25, min=0, max=10),
) -> None:
    """Run a bounded discovery campaign from a species CSV."""

    try:
        species_names = read_species_csv(species_csv)
        report = freeze_official_search_results(
            species_names,
            output_dir,
            max_queries=max_queries,
            google_api_key=google_api_key,
            google_cx=google_cx,
            brave_api_key=brave_api_key,
            results_per_query=results_per_query,
            timeout_seconds=timeout_seconds,
            max_attempts=max_attempts,
            backoff_seconds=backoff_seconds,
            max_backoff_seconds=max_backoff_seconds,
            min_interval_seconds=min_interval_seconds,
        )
    except (OSError, UnicodeError, csv.Error, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from None
    typer.echo(json.dumps(report, ensure_ascii=False, sort_keys=True))
