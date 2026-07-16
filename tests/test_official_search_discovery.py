from __future__ import annotations

import csv
import gzip
import hashlib
import json
import tomllib
import urllib.parse
from pathlib import Path

import pytest
from typer.testing import CliRunner

import island_v2.official_search_discovery as discovery_module
from island_v2.official_search_discovery import (
    BRAVE_API_KEY_ENV,
    DISCOVERY_COLUMNS,
    DISCOVERY_FILENAME,
    ERROR_COLUMNS,
    ERROR_FILENAME,
    GOOGLE_API_KEY_ENV,
    GOOGLE_CX_ENV,
    REPORT_FILENAME,
    SearchDiscoveryBatch,
    SearchHttpResponse,
    app,
    build_exact_species_query,
    discover_official_search_results,
    read_species_csv,
)


class FakeTransport:
    def __init__(self, responses):
        self.responses = list(responses)
        self.requests = []

    def __call__(self, request, timeout_seconds):
        self.requests.append((request, timeout_seconds))
        response = self.responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return response


def _response(payload: object, status: int = 200, headers=None) -> SearchHttpResponse:
    return SearchHttpResponse(
        status=status,
        headers=headers or {},
        body=json.dumps(payload).encode("utf-8"),
    )


def test_exact_species_query_and_explicit_credentials_are_required() -> None:
    query = build_exact_species_query("Plantus alba")
    assert query.startswith('"Plantus alba" (')
    assert '"breeding system"' in query
    assert '"self-incompatibility"' in query

    with pytest.raises(ValueError, match="explicit Google or Brave API key"):
        discover_official_search_results(["Plantus alba"], max_queries=1)
    with pytest.raises(ValueError, match="both google_api_key and google_cx"):
        discover_official_search_results(
            ["Plantus alba"],
            max_queries=1,
            google_api_key="google-key",
        )


def test_google_records_are_discovery_only_and_credentials_are_not_output() -> None:
    transport = FakeTransport(
        [
            _response(
                {
                    "items": [
                        {
                            "link": "https://flora.example/plantus#flowers",
                            "title": "&lt;i&gt;Plantus alba&lt;/i&gt; profile",
                            "snippet": "Flowers are described on the source page.",
                        },
                        {
                            "link": "javascript:alert(1)",
                            "title": "unsafe",
                            "snippet": "ignored",
                        },
                    ]
                }
            )
        ]
    )
    batch = discover_official_search_results(
        ["Plantus alba"],
        max_queries=1,
        google_api_key="google-secret",
        google_cx="engine-id",
        transport=transport,
        min_interval_seconds=0,
        sleeper=lambda _seconds: None,
        retrieved_at_utc="2026-07-15T00:00:00Z",
    )

    assert batch.queries_used == 1
    assert len(batch.records) == 1
    row = batch.records[0]
    assert list(row) == DISCOVERY_COLUMNS
    assert row["provider"] == "google_custom_search"
    assert row["result_url"] == "https://flora.example/plantus"
    assert row["title"] == "Plantus alba profile"
    assert "trait_name" not in row
    assert "trait_value" not in row
    assert "evidence" not in row
    assert "google-secret" not in json.dumps(row)
    stable = {
        key: value
        for key, value in row.items()
        if key not in {"retrieved_at_utc", "discovery_sha256"}
    }
    expected_hash = hashlib.sha256(
        json.dumps(
            stable,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest()
    assert row["discovery_sha256"] == expected_hash

    request = transport.requests[0][0]
    parameters = urllib.parse.parse_qs(urllib.parse.urlsplit(request.url).query)
    assert parameters["key"] == ["google-secret"]
    assert parameters["cx"] == ["engine-id"]
    assert parameters["q"] == [build_exact_species_query("Plantus alba")]


def test_brave_is_optional_fallback_after_google_zero_results() -> None:
    transport = FakeTransport(
        [
            _response({}),
            _response(
                {
                    "web": {
                        "results": [
                            {
                                "url": "https://garden.example/plantus",
                                "title": "Plantus alba",
                                "description": "A public horticultural profile.",
                            }
                        ]
                    }
                }
            ),
        ]
    )
    batch = discover_official_search_results(
        ["Plantus alba"],
        max_queries=2,
        google_api_key="google-key",
        google_cx="engine-id",
        brave_api_key="brave-key",
        transport=transport,
        min_interval_seconds=0,
        sleeper=lambda _seconds: None,
        retrieved_at_utc="2026-07-15T00:00:00Z",
    )

    assert batch.queries_used == 2
    assert [record["provider"] for record in batch.records] == ["brave_search"]
    assert batch.errors[0]["provider"] == "google_custom_search"
    assert batch.errors[0]["error_code"] == "no_results"
    brave_request = transport.requests[1][0]
    assert brave_request.headers["X-Subscription-Token"] == "brave-key"
    assert "brave-key" not in brave_request.url


def test_hard_budget_counts_retry_attempts_and_prevents_extra_http() -> None:
    transport = FakeTransport(
        [
            SearchHttpResponse(429, {"Retry-After": "0"}, b"{}"),
            _response({"items": []}),
        ]
    )
    batch = discover_official_search_results(
        ["Plantus alba", "Plantus rubra"],
        max_queries=1,
        google_api_key="google-key",
        google_cx="engine-id",
        brave_api_key="brave-key",
        max_attempts=3,
        transport=transport,
        min_interval_seconds=0,
        sleeper=lambda _seconds: None,
        retrieved_at_utc="2026-07-15T00:00:00Z",
    )

    assert len(transport.requests) == 1
    assert batch.queries_used == batch.max_queries == 1
    assert not batch.records
    assert [error["error_code"] for error in batch.errors] == [
        "budget_exhausted",
        "budget_exhausted",
    ]


def test_retry_backoff_is_bounded_and_success_uses_two_budget_units() -> None:
    sleeps: list[float] = []
    transport = FakeTransport(
        [
            SearchHttpResponse(503, {}, b"{}"),
            _response(
                {
                    "web": {
                        "results": [
                            {
                                "url": "https://plants.example/alba",
                                "title": "Plantus alba",
                                "description": "Search lead only.",
                            }
                        ]
                    }
                }
            ),
        ]
    )
    batch = discover_official_search_results(
        ["Plantus alba"],
        max_queries=2,
        brave_api_key="brave-key",
        max_attempts=2,
        backoff_seconds=0.5,
        max_backoff_seconds=1,
        transport=transport,
        min_interval_seconds=0,
        sleeper=sleeps.append,
        retrieved_at_utc="2026-07-15T00:00:00Z",
    )

    assert batch.queries_used == 2
    assert sleeps == [0.5]
    assert len(batch.records) == 1
    assert not batch.errors


def test_species_csv_reader_supports_plain_and_gzip_and_requires_species(
    tmp_path: Path,
) -> None:
    plain = tmp_path / "species.csv"
    plain.write_text(
        "species,note\nPlantus alba,one\n,blank\nPlantus rubra,two\n",
        encoding="utf-8",
    )
    compressed = tmp_path / "species.csv.gz"
    with gzip.open(compressed, mode="wt", encoding="utf-8", newline="") as handle:
        handle.write("species\nPlantus caerulea\n")

    assert read_species_csv(plain) == ("Plantus alba", "Plantus rubra")
    assert read_species_csv(compressed) == ("Plantus caerulea",)

    missing = tmp_path / "missing.csv"
    missing.write_text("accepted_species\nPlantus alba\n", encoding="utf-8")
    with pytest.raises(ValueError, match="column named 'species'"):
        read_species_csv(missing)


def test_cli_uses_environment_credentials_and_atomically_writes_discovery_only_outputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    species_path = tmp_path / "species.csv.gz"
    with gzip.open(species_path, mode="wt", encoding="utf-8", newline="") as handle:
        handle.write("species\nPlantus alba\nPlantus rubra\n")
    output_dir = tmp_path / "output"
    captured: dict[str, object] = {}
    record = {
        "species": "Plantus alba",
        "provider": "google_custom_search",
        "query": build_exact_species_query("Plantus alba"),
        "rank": "1",
        "result_url": "https://flora.example/plantus",
        "title": "Plantus alba flora",
        "snippet": "A search lead, not trait evidence.",
        "retrieved_at_utc": "2026-07-15T00:00:00Z",
        "discovery_sha256": "a" * 64,
    }
    error = {
        "species": "Plantus rubra",
        "provider": "google_custom_search",
        "query": build_exact_species_query("Plantus rubra"),
        "error_class": "terminal",
        "error_code": "no_results",
        "message": "search returned no valid discovery URLs",
    }

    def fake_discovery(species_names, **kwargs):
        captured["species_names"] = tuple(species_names)
        captured["kwargs"] = kwargs
        return SearchDiscoveryBatch(
            requested_species=("Plantus alba", "Plantus rubra"),
            records=(record,),
            errors=(error,),
            queries_used=2,
            max_queries=2,
        )

    monkeypatch.setattr(
        discovery_module,
        "discover_official_search_results",
        fake_discovery,
    )
    google_secret = "google-secret-must-not-leak"
    cx_secret = "cx-secret-must-not-leak"
    result = CliRunner().invoke(
        app,
        [
            "--species-csv",
            str(species_path),
            "--output-dir",
            str(output_dir),
            "--max-queries",
            "2",
        ],
        env={
            GOOGLE_API_KEY_ENV: google_secret,
            GOOGLE_CX_ENV: cx_secret,
            BRAVE_API_KEY_ENV: "",
        },
    )

    assert result.exit_code == 0, result.output
    assert captured["species_names"] == ("Plantus alba", "Plantus rubra")
    kwargs = captured["kwargs"]
    assert isinstance(kwargs, dict)
    assert kwargs["max_queries"] == 2
    assert kwargs["google_api_key"] == google_secret
    assert kwargs["google_cx"] == cx_secret
    assert set(path.name for path in output_dir.iterdir()) == {
        DISCOVERY_FILENAME,
        ERROR_FILENAME,
        REPORT_FILENAME,
    }
    with (output_dir / DISCOVERY_FILENAME).open(encoding="utf-8", newline="") as handle:
        discovery_rows = list(csv.DictReader(handle))
    with (output_dir / ERROR_FILENAME).open(encoding="utf-8", newline="") as handle:
        error_rows = list(csv.DictReader(handle))
    assert list(discovery_rows[0]) == DISCOVERY_COLUMNS
    assert list(error_rows[0]) == ERROR_COLUMNS
    assert not ({"trait_name", "trait_value", "evidence_type"} & discovery_rows[0].keys())
    report = json.loads((output_dir / REPORT_FILENAME).read_text(encoding="utf-8"))
    assert report["queries_used"] == 2
    assert report["max_queries"] == 2
    assert report["policy"] == {
        "artifact_scope": "search_result_discovery_only",
        "snippets_are_trait_evidence": False,
        "trait_values_emitted": False,
    }
    for filename in (DISCOVERY_FILENAME, ERROR_FILENAME):
        assert (
            report["files"][filename]["sha256"]
            == hashlib.sha256((output_dir / filename).read_bytes()).hexdigest()
        )
    rendered = result.output + "".join(
        path.read_text(encoding="utf-8") for path in output_dir.iterdir()
    )
    assert google_secret not in rendered
    assert cx_secret not in rendered
    assert not list(output_dir.glob("*.tmp"))


def test_cli_requires_hard_budget_and_help_never_renders_secret_values(
    tmp_path: Path,
) -> None:
    species_path = tmp_path / "species.csv"
    species_path.write_text("species\nPlantus alba\n", encoding="utf-8")
    runner = CliRunner()
    missing_budget = runner.invoke(
        app,
        [
            "--species-csv",
            str(species_path),
            "--output-dir",
            str(tmp_path / "output"),
        ],
        env={BRAVE_API_KEY_ENV: "brave-secret"},
    )
    assert missing_budget.exit_code == 2
    assert "--max-queries" in missing_budget.output

    zero_budget = runner.invoke(
        app,
        [
            "--species-csv",
            str(species_path),
            "--output-dir",
            str(tmp_path / "zero-output"),
            "--max-queries",
            "0",
        ],
        env={BRAVE_API_KEY_ENV: "brave-secret"},
    )
    assert zero_budget.exit_code == 2
    assert "--max-queries" in zero_budget.output

    explicit_secret = "explicit-google-secret"
    incomplete_google = runner.invoke(
        app,
        [
            "--species-csv",
            str(species_path),
            "--output-dir",
            str(tmp_path / "incomplete-output"),
            "--max-queries",
            "1",
            "--google-api-key",
            explicit_secret,
        ],
        env={
            GOOGLE_API_KEY_ENV: "",
            GOOGLE_CX_ENV: "",
            BRAVE_API_KEY_ENV: "",
        },
    )
    assert incomplete_google.exit_code == 2
    assert "Google Custom Search requires" in incomplete_google.output
    assert "google_cx" in incomplete_google.output
    assert explicit_secret not in incomplete_google.output

    google_secret = "google-help-secret"
    brave_secret = "brave-help-secret"
    help_result = runner.invoke(
        app,
        ["--help"],
        env={
            GOOGLE_API_KEY_ENV: google_secret,
            GOOGLE_CX_ENV: "cx-help-secret",
            BRAVE_API_KEY_ENV: brave_secret,
        },
    )
    assert help_result.exit_code == 0
    assert "--google-api-key" in help_result.output
    assert "--google-cx" in help_result.output
    assert "--brave-api-key" in help_result.output
    assert google_secret not in help_result.output
    assert brave_secret not in help_result.output
    assert google_secret not in repr(app)
    assert brave_secret not in repr(app)


def test_pyproject_registers_official_search_discovery_entry_point() -> None:
    project = tomllib.loads(
        (Path(__file__).parents[1] / "pyproject.toml").read_text(encoding="utf-8")
    )
    assert project["project"]["scripts"]["island-v2-official-search-discovery"] == (
        "island_v2.official_search_discovery:app"
    )
