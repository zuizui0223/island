from __future__ import annotations

import csv
import hashlib
import json
import urllib.error
import urllib.parse
from pathlib import Path

from island_v2.europe_pmc_discovery import (
    HttpResponse,
    build_query,
    build_search_api_url,
    discover_species_sources,
    freeze_europe_pmc_sources,
)


class FakeTransport:
    def __init__(self, responses):
        self.responses = list(responses)
        self.calls: list[tuple[str, float]] = []

    def __call__(self, url: str, timeout_seconds: float) -> HttpResponse:
        self.calls.append((url, timeout_seconds))
        response = self.responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return response


def _json_response(payload: object, status: int = 200, headers=None) -> HttpResponse:
    return HttpResponse(
        status=status,
        headers=headers or {},
        body=json.dumps(payload).encode("utf-8"),
    )


def _matching_payload() -> dict[str, object]:
    return {
        "resultList": {
            "result": [
                {
                    "id": "12345678",
                    "source": "MED",
                    "pmid": "12345678",
                    "pmcid": "PMC7654321",
                    "doi": "https://doi.org/10.1000/PLANT.1",
                    "title": "Pollination biology of &lt;i&gt;Plantus alba&lt;/i&gt;.",
                    "abstractText": (
                        "The breeding system of <i>Plantus alba</i> was studied in the field."
                    ),
                    "authorString": "Botanist A, Ecologist B.",
                    "journalInfo": {"journal": {"title": "Plant Reproduction"}},
                    "pubYear": "2024",
                    "firstPublicationDate": "2024-06-01",
                    "license": "cc by",
                    "isOpenAccess": "Y",
                },
                {
                    "id": "87654321",
                    "source": "MED",
                    "title": "Pollination in Plantus species",
                    "abstractText": "P. alba was visited by insects.",
                },
            ]
        }
    }


def test_query_is_exact_quoted_species_and_title_abstract_reproductive_terms() -> None:
    query = build_query("Plantus alba")
    assert query.startswith('TITLE_ABS:"Plantus alba" AND (')
    assert "TITLE_ABS:pollinat*" in query
    assert 'TITLE_ABS:"breeding system"' in query
    assert 'TITLE_ABS:"self-incompatibility"' in query

    url = build_search_api_url("Plantus alba", page_size=17)
    parameters = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)
    assert parameters == {
        "query": [query],
        "resultType": ["core"],
        "pageSize": ["17"],
        "format": ["json"],
    }


def test_discovery_keeps_only_full_species_title_abstract_metadata() -> None:
    transport = FakeTransport([_json_response(_matching_payload())])
    batch = discover_species_sources(
        "Plantus alba",
        transport=transport,
        sleeper=lambda _seconds: None,
        retrieved_at_utc="2026-07-15T00:00:00Z",
    )

    assert not batch.errors
    assert len(batch.sources) == 1
    row = batch.sources[0]
    assert row["accepted_species"] == "Plantus alba"
    assert row["title"] == "Pollination biology of Plantus alba."
    assert row["abstract"].startswith("The breeding system of Plantus alba")
    assert row["pmid"] == "12345678"
    assert row["pmcid"] == "PMC7654321"
    assert row["doi"] == "10.1000/plant.1"
    assert row["license"] == "cc by"
    assert row["evidence_scope"] == "species_direct"
    assert row["source_url"] == "https://europepmc.org/article/MED/12345678"
    assert row["article_api_url"].startswith(
        "https://www.ebi.ac.uk/europepmc/webservices/rest/article/MED/12345678"
    )
    assert "DOI:10.1000/plant.1" in row["source_citation"]
    assert (
        row["source_text_sha256"] == hashlib.sha256(row["source_text"].encode("utf-8")).hexdigest()
    )
    assert len(row["metadata_sha256"]) == 64
    assert "trait_name" not in row
    assert "value" not in row


def test_retry_after_is_honoured_then_zero_hit_is_terminal() -> None:
    sleeps: list[float] = []
    transport = FakeTransport(
        [
            HttpResponse(429, {"Retry-After": "2"}, b"{}"),
            _json_response({"resultList": {"result": []}}),
        ]
    )
    batch = discover_species_sources(
        "Plantus alba",
        transport=transport,
        max_attempts=2,
        sleeper=sleeps.append,
        retrieved_at_utc="2026-07-15T00:00:00Z",
    )

    assert len(transport.calls) == 2
    assert sleeps == [2]
    assert batch.errors[0]["error_class"] == "terminal"
    assert batch.errors[0]["error_code"] == "no_results"
    assert batch.errors[0]["attempts"] == "2"


def test_terminal_http_does_not_retry_and_transport_exhaustion_remains_retryable() -> None:
    terminal_transport = FakeTransport([HttpResponse(400, {}, b"bad query")])
    terminal = discover_species_sources(
        "Plantus alba",
        transport=terminal_transport,
        max_attempts=3,
        sleeper=lambda _seconds: None,
        retrieved_at_utc="2026-07-15T00:00:00Z",
    )
    assert len(terminal_transport.calls) == 1
    assert terminal.errors[0]["error_class"] == "terminal"
    assert terminal.errors[0]["error_code"] == "http_400"

    retry_sleeps: list[float] = []
    retry_transport = FakeTransport(
        [urllib.error.URLError("temporary"), urllib.error.URLError("temporary")]
    )
    retryable = discover_species_sources(
        "Plantus alba",
        transport=retry_transport,
        max_attempts=2,
        backoff_seconds=0.5,
        sleeper=retry_sleeps.append,
        retrieved_at_utc="2026-07-15T00:00:00Z",
    )
    assert retry_sleeps == [0.5]
    assert retryable.errors[0]["error_class"] == "retryable"
    assert retryable.errors[0]["error_code"] == "transport_error"


def test_freeze_writes_hashed_sources_errors_and_no_inference_manifest(tmp_path: Path) -> None:
    transport = FakeTransport([_json_response(_matching_payload())])
    output_dir = tmp_path / "europe_pmc"
    manifest = freeze_europe_pmc_sources(
        ["Plantus alba"],
        output_dir,
        transport=transport,
        min_interval_seconds=0,
        sleeper=lambda _seconds: None,
        retrieved_at_utc="2026-07-15T00:00:00Z",
    )

    source_path = output_dir / "europe_pmc_source_texts.csv"
    error_path = output_dir / "europe_pmc_discovery_errors.csv"
    with source_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1
    assert rows[0]["evidence_scope"] == "species_direct"
    assert list(csv.DictReader(error_path.open(encoding="utf-8", newline=""))) == []
    assert manifest["n_requested_species"] == 1
    assert manifest["n_source_records"] == 1
    assert manifest["n_errors"] == 0
    assert "infers no plant trait" in manifest["policy"]
    assert (
        manifest["files"][source_path.name]["sha256"]
        == hashlib.sha256(source_path.read_bytes()).hexdigest()
    )
