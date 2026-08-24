"""Pluggable, discovery-only open-Web search for unresolved plant traits.

Search is deliberately separated from evidence extraction.  Adapters return
URLs and snippets, but snippets can never cross the evidence interface.  The
same task manifest can be replayed through Google Custom Search, Brave Search,
or deterministic fixtures without changing query generation or downstream
identity checks.
"""

# Typer's documented command style evaluates Option metadata in defaults.
# ruff: noqa: B008

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import urllib.parse
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, fields
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Protocol

import httpx
import typer
import yaml

GOOGLE_URL = "https://customsearch.googleapis.com/customsearch/v1"
BRAVE_URL = "https://api.search.brave.com/res/v1/web/search"
USER_AGENT = "island-floral-traits/0.1 (open-web discovery; contact via repository)"
CONTRACT = "open_web_search_discovery_v1"

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _text(value: object) -> str:
    return " ".join(str(value or "").replace("\xa0", " ").split())


def _now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _sha(payload: Mapping[str, object]) -> str:
    raw = json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _safe_url(value: object) -> str:
    raw = _text(value)
    if not raw:
        return ""
    parsed = urllib.parse.urlsplit(raw)
    if parsed.scheme.casefold() not in {"http", "https"} or not parsed.netloc:
        return ""
    return urllib.parse.urlunsplit(
        (parsed.scheme.casefold(), parsed.netloc.casefold(), parsed.path, parsed.query, "")
    )


@dataclass(frozen=True)
class SearchTask:
    task_id: str
    accepted_species: str
    trait_name: str
    searched_name: str
    name_kind: str
    language: str
    query: str
    priority_score: float
    potential_species_trait_unlocked: int
    expected_information_gain: float
    query_cost_units: int = 1


@dataclass(frozen=True)
class SearchHit:
    task_id: str
    accepted_species: str
    trait_name: str
    searched_name: str
    name_kind: str
    language: str
    query: str
    provider: str
    rank: int
    result_url: str
    result_domain: str
    title: str
    snippet_discovery_only: str
    retrieved_at_utc: str
    discovery_sha256: str


class SearchBackend(Protocol):
    """One real seam: adapters vary, query/evidence policy does not."""

    name: str
    cost_per_query_usd: float | None

    def search(self, task: SearchTask, *, count: int) -> Sequence[Mapping[str, object]]:
        """Return raw URL/title/snippet mappings; never trait evidence."""


class GoogleCustomSearchBackend:
    name = "google_custom_search"
    cost_per_query_usd: float | None = None

    def __init__(
        self,
        api_key: str,
        cx: str,
        *,
        client: httpx.Client | None = None,
    ) -> None:
        if not _text(api_key) or not _text(cx):
            raise ValueError("Google Custom Search requires an API key and engine ID")
        self._api_key = api_key
        self._cx = cx
        self._client = client or httpx.Client(timeout=30, follow_redirects=True)

    def search(self, task: SearchTask, *, count: int) -> Sequence[Mapping[str, object]]:
        response = self._client.get(
            GOOGLE_URL,
            params={
                "key": self._api_key,
                "cx": self._cx,
                "q": task.query,
                "num": min(10, count),
                "safe": "active",
                "hl": task.language if task.language != "la" else "en",
            },
            headers={"Accept": "application/json", "User-Agent": USER_AGENT},
        )
        response.raise_for_status()
        payload = response.json()
        return [
            {
                "url": item.get("link", ""),
                "title": item.get("title", ""),
                "snippet": item.get("snippet", ""),
            }
            for item in payload.get("items", [])
            if isinstance(item, Mapping)
        ]


class BraveSearchBackend:
    name = "brave_search"
    cost_per_query_usd: float | None = None

    def __init__(self, api_key: str, *, client: httpx.Client | None = None) -> None:
        if not _text(api_key):
            raise ValueError("Brave Search requires an API key")
        self._api_key = api_key
        self._client = client or httpx.Client(timeout=30, follow_redirects=True)

    def search(self, task: SearchTask, *, count: int) -> Sequence[Mapping[str, object]]:
        language = "en" if task.language == "la" else task.language
        response = self._client.get(
            BRAVE_URL,
            params={
                "q": task.query,
                "count": min(20, count),
                "safesearch": "moderate",
                "search_lang": language,
            },
            headers={
                "Accept": "application/json",
                "User-Agent": USER_AGENT,
                "X-Subscription-Token": self._api_key,
            },
        )
        response.raise_for_status()
        payload = response.json()
        web = payload.get("web", {}) if isinstance(payload, Mapping) else {}
        items = web.get("results", []) if isinstance(web, Mapping) else []
        return [
            {
                "url": item.get("url", ""),
                "title": item.get("title", ""),
                "snippet": item.get("description", ""),
            }
            for item in items
            if isinstance(item, Mapping)
        ]


class FixtureSearchBackend:
    name = "fixture"
    cost_per_query_usd = 0.0

    def __init__(self, fixture: Mapping[str, Sequence[Mapping[str, object]]]) -> None:
        self._fixture = fixture

    @classmethod
    def from_json(cls, path: Path) -> FixtureSearchBackend:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(payload, Mapping):
            raise TypeError("search fixture must be a query-to-results object")
        return cls(payload)  # type: ignore[arg-type]

    def search(self, task: SearchTask, *, count: int) -> Sequence[Mapping[str, object]]:
        return list(self._fixture.get(task.query, []))[:count]


def load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise TypeError("open-Web acquisition config must be a mapping")
    return payload


def _query(
    name: str,
    labels: Sequence[str],
    *,
    pdf: bool,
) -> str:
    quoted_labels = " OR ".join(f'"{_text(label)}"' for label in labels if _text(label))
    suffix = " filetype:pdf" if pdf else ""
    return f'"{_text(name)}" ({quoted_labels}){suffix}'


def build_query_tasks(
    queue: Sequence[Mapping[str, object]],
    *,
    config: Mapping[str, Any],
    synonyms: Mapping[str, Sequence[object]] | None = None,
    vernacular_names: Mapping[str, Sequence[tuple[str, str]]] | None = None,
    max_names_per_species: int = 3,
    include_pdf_query: bool = True,
) -> list[SearchTask]:
    """Expand prioritized species-trait rows into deterministic multilingual tasks."""

    languages = config.get("languages", {})
    aliases = synonyms or {}
    local_names = vernacular_names or {}
    tasks: list[SearchTask] = []
    seen: set[tuple[str, str, str, str]] = set()
    for row in queue:
        species = _text(row.get("accepted_species"))
        trait = _text(row.get("trait_name"))
        if not species or trait not in {
            name
            for language in languages.values()
            for name in (language.get("query_labels", {}) if isinstance(language, dict) else {})
        }:
            continue
        name_variants: list[tuple[str, str, str]] = [(species, "accepted", "")]
        scientific_aliases: list[tuple[str, str, str]] = []
        for raw_alias in aliases.get(species, ()):
            if isinstance(raw_alias, (tuple, list)) and len(raw_alias) >= 2:
                alias_name = _text(raw_alias[0])
                alias_kind = _text(raw_alias[1]).casefold()
            else:
                alias_name = _text(raw_alias)
                alias_kind = "synonym"
            if alias_name and alias_kind in {"synonym", "basionym"}:
                scientific_aliases.append((alias_name, alias_kind, ""))
        name_variants.extend(list(dict.fromkeys(scientific_aliases))[:max_names_per_species])
        name_variants.extend(
            (name, "vernacular", language)
            for name, language in list(dict.fromkeys(local_names.get(species, ())))[
                :max_names_per_species
            ]
        )
        distribution_languages = [
            _text(value).casefold()
            for value in re.split(
                r"[|,; ]+",
                _text(row.get("distribution_languages") or row.get("search_languages")),
            )
            if _text(value)
        ]
        requested_languages = list(dict.fromkeys(["en", "la", *distribution_languages]))
        for name, kind, name_language in name_variants[: 1 + max_names_per_species * 2]:
            languages_for_name = [name_language] if kind == "vernacular" else requested_languages
            for language in languages_for_name:
                language_cfg = languages.get(language, {})
                labels = (
                    language_cfg.get("query_labels", {}).get(trait, [])
                    if isinstance(language_cfg, dict)
                    else []
                )
                if not labels:
                    continue
                effective_language = language
                for pdf in [False, True] if include_pdf_query else [False]:
                    query = _query(name, labels[:3], pdf=pdf)
                    identity = (species, trait, effective_language, query)
                    if identity in seen:
                        continue
                    seen.add(identity)
                    stable = {
                        "accepted_species": species,
                        "trait_name": trait,
                        "searched_name": name,
                        "name_kind": kind,
                        "language": effective_language,
                        "query": query,
                    }
                    tasks.append(
                        SearchTask(
                            task_id=_sha(stable)[:24],
                            accepted_species=species,
                            trait_name=trait,
                            searched_name=name,
                            name_kind=kind,
                            language=effective_language,
                            query=query,
                            priority_score=float(row.get("priority_score") or 0),
                            potential_species_trait_unlocked=int(
                                float(row.get("potential_species_trait_unlocked") or 0)
                            ),
                            expected_information_gain=float(
                                row.get("expected_information_gain") or 0
                            ),
                        )
                    )
    return sorted(
        tasks,
        key=lambda item: (
            -item.priority_score,
            -item.expected_information_gain,
            item.accepted_species,
            item.trait_name,
            item.task_id,
        ),
    )


def discover(
    tasks: Iterable[SearchTask],
    *,
    backends: Sequence[SearchBackend],
    max_queries: int,
    results_per_query: int = 10,
) -> tuple[list[SearchHit], list[dict[str, object]], dict[str, object]]:
    """Run tasks with provider fallback under a hard request budget."""

    if not backends:
        raise ValueError("at least one search backend is required")
    if max_queries < 0:
        raise ValueError("max_queries cannot be negative")
    used = 0
    hits: list[SearchHit] = []
    errors: list[dict[str, object]] = []
    attempted_task_ids: list[str] = []
    by_backend: dict[str, int] = {backend.name: 0 for backend in backends}
    for task in tasks:
        if used >= max_queries:
            break
        task_attempted = False
        for backend in backends:
            if used >= max_queries:
                break
            used += 1
            task_attempted = True
            by_backend[backend.name] += 1
            try:
                rows = backend.search(task, count=results_per_query)
            except (httpx.HTTPError, ValueError, KeyError, json.JSONDecodeError) as exc:
                errors.append(
                    {
                        "task_id": task.task_id,
                        "accepted_species": task.accepted_species,
                        "trait_name": task.trait_name,
                        "provider": backend.name,
                        "error": type(exc).__name__,
                        "message": _text(exc),
                    }
                )
                continue
            retrieved = _now()
            backend_hits: list[SearchHit] = []
            for rank, row in enumerate(rows, start=1):
                url = _safe_url(row.get("url"))
                if not url:
                    continue
                stable = {
                    "task_id": task.task_id,
                    "provider": backend.name,
                    "rank": rank,
                    "result_url": url,
                }
                backend_hits.append(
                    SearchHit(
                        task_id=task.task_id,
                        accepted_species=task.accepted_species,
                        trait_name=task.trait_name,
                        searched_name=task.searched_name,
                        name_kind=task.name_kind,
                        language=task.language,
                        query=task.query,
                        provider=backend.name,
                        rank=rank,
                        result_url=url,
                        result_domain=urllib.parse.urlsplit(url)
                        .netloc.casefold()
                        .removeprefix("www."),
                        title=_text(row.get("title")),
                        snippet_discovery_only=_text(row.get("snippet")),
                        retrieved_at_utc=retrieved,
                        discovery_sha256=_sha(stable),
                    )
                )
            hits.extend(backend_hits)
            if backend_hits:
                break
        if task_attempted:
            attempted_task_ids.append(task.task_id)
    costs = {
        backend.name: (
            None
            if backend.cost_per_query_usd is None
            else round(by_backend[backend.name] * backend.cost_per_query_usd, 6)
        )
        for backend in backends
    }
    report: dict[str, object] = {
        "contract": CONTRACT,
        "completed_at_utc": _now(),
        "queries_used": used,
        "max_queries": max_queries,
        "hits": len(hits),
        "unique_urls": len({hit.result_url for hit in hits}),
        "provider_query_counts": by_backend,
        "provider_cost_usd": costs,
        "attempted_task_count": len(attempted_task_ids),
        "attempted_task_ids": attempted_task_ids,
        "snippets_are_evidence": False,
        "source_pages_must_be_fetched": True,
    }
    return hits, errors, report


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def load_strict_synonyms(path: Path | None) -> dict[str, list[str]]:
    """Load only synonym mappings corroborated by two named backbones."""

    if path is None:
        return {}
    grouped: dict[tuple[str, str], set[str]] = {}
    for row in read_csv_rows(path):
        accepted = _text(row.get("accepted_species"))
        synonym = _text(row.get("synonym"))
        backbone = _text(row.get("backbone"))
        if accepted and synonym and backbone:
            grouped.setdefault((accepted, synonym), set()).add(backbone.casefold())
    result: dict[str, list[str]] = {}
    for (accepted, synonym), backbones in grouped.items():
        if len(backbones) >= 2:
            result.setdefault(accepted, []).append(synonym)
    return {accepted: sorted(set(names)) for accepted, names in result.items()}


def load_strict_name_variants(
    path: Path | None,
) -> dict[str, list[tuple[str, str]]]:
    """Load two-backbone-agreed synonyms and basionyms with their name kind."""

    if path is None:
        return {}
    grouped: dict[tuple[str, str, str], set[str]] = {}
    for row in read_csv_rows(path):
        accepted = _text(row.get("accepted_species"))
        explicit_basionym = _text(row.get("basionym"))
        name = _text(row.get("synonym") or explicit_basionym or row.get("name"))
        kind = _text(row.get("name_kind")).casefold()
        if not kind:
            kind = "basionym" if explicit_basionym else "synonym"
        backbone = _text(row.get("backbone"))
        if accepted and name and kind in {"synonym", "basionym"} and backbone:
            grouped.setdefault((accepted, name, kind), set()).add(backbone.casefold())
    result: dict[str, list[tuple[str, str]]] = {}
    for (accepted, name, kind), backbones in grouped.items():
        if len(backbones) >= 2:
            result.setdefault(accepted, []).append((name, kind))
    return {
        accepted: sorted(set(names), key=lambda item: (item[1], item[0].casefold()))
        for accepted, names in result.items()
    }


def load_vernacular_names(path: Path | None) -> dict[str, list[tuple[str, str]]]:
    if path is None:
        return {}
    result: dict[str, list[tuple[str, str]]] = {}
    for row in read_csv_rows(path):
        accepted = _text(row.get("accepted_species"))
        name = _text(row.get("vernacular_name") or row.get("local_name"))
        language = _text(row.get("language")).casefold()
        if accepted and name and language:
            result.setdefault(accepted, []).append((name, language))
    return {
        accepted: sorted(set(names), key=lambda item: (item[1], item[0].casefold()))
        for accepted, names in result.items()
    }


def write_csv(
    path: Path,
    rows: Sequence[Mapping[str, object]],
    *,
    fieldnames: Sequence[str] = (),
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(rows[0]) if rows else list(fieldnames)
    if not columns:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _backends(name: str, fixture_json: Path | None) -> list[SearchBackend]:
    if name == "fixture":
        if fixture_json is None:
            raise ValueError("fixture backend requires --fixture-json")
        return [FixtureSearchBackend.from_json(fixture_json)]
    google_key = _text(os.getenv("GOOGLE_CUSTOM_SEARCH_API_KEY"))
    google_cx = _text(os.getenv("GOOGLE_CUSTOM_SEARCH_CX"))
    brave_key = _text(os.getenv("BRAVE_SEARCH_API_KEY"))
    result: list[SearchBackend] = []
    if name in {"auto", "google"} and google_key and google_cx:
        result.append(GoogleCustomSearchBackend(google_key, google_cx))
    if name in {"auto", "brave"} and brave_key:
        result.append(BraveSearchBackend(brave_key))
    if not result:
        raise ValueError(
            "no requested official-search credential is configured; "
            "set Brave or Google CSE secrets, or use the fixture backend"
        )
    return result


@app.command("build-queries")
def build_queries_command(
    queue_csv: Path = typer.Option(..., exists=True),
    config_yaml: Path = typer.Option(Path("config/open_web_trait_acquisition.yml"), exists=True),
    output_csv: Path = typer.Option(...),
    synonym_csv: Path | None = typer.Option(
        None,
        help="Optional synonym,accepted_species,backbone snapshot; two backbones must agree",
        exists=True,
    ),
    vernacular_csv: Path | None = typer.Option(
        None,
        help="Optional accepted_species,vernacular_name,language table",
        exists=True,
    ),
    max_names_per_species: int = typer.Option(3, min=0, max=10),
) -> None:
    tasks = build_query_tasks(
        read_csv_rows(queue_csv),
        config=load_config(config_yaml),
        synonyms=load_strict_name_variants(synonym_csv),
        vernacular_names=load_vernacular_names(vernacular_csv),
        max_names_per_species=max_names_per_species,
    )
    write_csv(
        output_csv,
        [asdict(task) for task in tasks],
        fieldnames=[field.name for field in fields(SearchTask)],
    )
    typer.echo(json.dumps({"tasks": len(tasks), "output": str(output_csv)}, sort_keys=True))


@app.command("discover")
def discover_command(
    tasks_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    backend: str = typer.Option("auto", help="auto, google, brave, or fixture"),
    fixture_json: Path | None = typer.Option(None),
    max_queries: int = typer.Option(..., min=0),
    results_per_query: int = typer.Option(10, min=1, max=20),
) -> None:
    tasks = [
        SearchTask(
            task_id=_text(row.get("task_id")),
            accepted_species=_text(row.get("accepted_species")),
            trait_name=_text(row.get("trait_name")),
            searched_name=_text(row.get("searched_name")),
            name_kind=_text(row.get("name_kind")),
            language=_text(row.get("language")),
            query=_text(row.get("query")),
            priority_score=float(row.get("priority_score") or 0),
            potential_species_trait_unlocked=int(
                float(row.get("potential_species_trait_unlocked") or 0)
            ),
            expected_information_gain=float(row.get("expected_information_gain") or 0),
            query_cost_units=int(float(row.get("query_cost_units") or 1)),
        )
        for row in read_csv_rows(tasks_csv)
    ]
    output_dir.mkdir(parents=True, exist_ok=True)
    try:
        adapters = _backends(backend, fixture_json)
        hits, errors, report = discover(
            tasks,
            backends=adapters,
            max_queries=max_queries,
            results_per_query=results_per_query,
        )
        status = "completed"
    except ValueError as exc:
        hits, errors = [], []
        report = {
            "contract": CONTRACT,
            "completed_at_utc": _now(),
            "status": "blocked_missing_or_invalid_search_credentials",
            "message": _text(exc),
            "queries_used": 0,
            "max_queries": max_queries,
            "hits": 0,
            "snippets_are_evidence": False,
        }
        status = "blocked"
    write_csv(
        output_dir / "search_hits.csv",
        [asdict(hit) for hit in hits],
        fieldnames=[field.name for field in fields(SearchHit)],
    )
    write_csv(
        output_dir / "search_errors.csv",
        errors,
        fieldnames=(
            "task_id",
            "accepted_species",
            "trait_name",
            "provider",
            "error",
            "message",
        ),
    )
    domain_rows: list[dict[str, object]] = []
    for domain in sorted({hit.result_domain for hit in hits}):
        domain_hits = [hit for hit in hits if hit.result_domain == domain]
        domain_rows.append(
            {
                "domain": domain,
                "hit_count": len(domain_hits),
                "unique_tasks": len({hit.task_id for hit in domain_hits}),
                "source_name": "",
                "source_type": "unclassified",
                "region": "",
                "language": "",
                "organization_type": "",
                "robots_access_notes": "not_yet_audited",
                "page_pattern": "",
                "scientific_name_displayed": "",
                "cultivar_descriptions_present": "",
                "provenance_quality": "pending",
                "recommended_evidence_tier": "D_discovery_only_until_reviewed",
                "success_rate": "",
                "false_positive_audit": "pending",
                "review_status": "pending_registry_review",
            }
        )
    write_csv(
        output_dir / "discovered_domain_registry.csv",
        domain_rows,
        fieldnames=(
            "domain",
            "hit_count",
            "unique_tasks",
            "source_name",
            "source_type",
            "region",
            "language",
            "organization_type",
            "robots_access_notes",
            "page_pattern",
            "scientific_name_displayed",
            "cultivar_descriptions_present",
            "provenance_quality",
            "recommended_evidence_tier",
            "success_rate",
            "false_positive_audit",
            "review_status",
        ),
    )
    (output_dir / "search_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps({"status": status, **report}, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    app()
