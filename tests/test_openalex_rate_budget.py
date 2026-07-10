from datetime import UTC, datetime
from pathlib import Path

import httpx
import pandas as pd
import pytest

from island_v2 import global_trait_campaign as campaign


class FakeClient:
    def __init__(self, responses):
        self.responses = list(responses)
        self.calls = 0

    def get(self, url, params):
        self.calls += 1
        return self.responses.pop(0)


def _response(status: int, *, headers=None, payload=None):
    return httpx.Response(
        status,
        headers=headers or {},
        json=payload or {"results": []},
        request=httpx.Request("GET", "https://api.openalex.org/works"),
    )


def _config():
    return campaign.load_config(Path("config/global_trait_campaign.yml"))


def _ledger(config):
    master = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "genus": "Alpha",
                "family": "FamA",
                "n_islands": 2,
                "n_records": 10,
            }
        ]
    )
    return campaign.reconcile_ledger(master, None, config)


def test_retry_after_parses_seconds_and_http_date():
    now = datetime(2026, 7, 10, 12, 0, tzinfo=UTC)

    assert campaign._retry_after_seconds("12", now=now) == 12
    assert campaign._retry_after_seconds("Fri, 10 Jul 2026 12:00:30 GMT", now=now) == 30
    assert campaign._retry_after_seconds("garbage", now=now) is None


def test_openalex_request_honours_short_retry_after_then_succeeds():
    sleeps = []
    client = FakeClient(
        [
            _response(429, headers={"Retry-After": "2"}),
            _response(200),
        ]
    )

    response = campaign._openalex_request(
        client,
        {"search": "test"},
        max_retries=1,
        backoff_seconds=5,
        max_backoff_seconds=60,
        sleeper=sleeps.append,
    )

    assert response.status_code == 200
    assert client.calls == 2
    assert sleeps == [2]


def test_openalex_request_does_not_burn_calls_on_daily_quota_429():
    client = FakeClient([_response(429)])

    with pytest.raises(campaign.OpenAlexTransientError, match="status=429"):
        campaign._openalex_request(
            client,
            {"search": "test"},
            max_retries=3,
            backoff_seconds=5,
            max_backoff_seconds=60,
            sleeper=lambda _: None,
        )

    assert client.calls == 1


def test_missing_api_key_returns_transient_errors_without_http(monkeypatch):
    monkeypatch.delenv("OPENALEX_API_KEY", raising=False)
    config = _config()
    batch = pd.DataFrame(
        {
            "accepted_species": ["Alpha one", "Beta two"],
            "genus": ["Alpha", "Beta"],
            "family": ["FamA", "FamB"],
            "n_islands": [2, 1],
            "n_records": [10, 5],
        }
    )

    sources, errors = campaign.fetch_openalex_sources(
        batch,
        "reproductive_openalex",
        config,
    )

    assert sources.empty
    assert len(errors) == 2
    assert all(":transient:missing_OPENALEX_API_KEY" in error for error in errors)


def test_transient_openalex_failure_does_not_consume_attempt_budget():
    config = _config()
    ledger = _ledger(config)
    task = "reproductive_openalex"
    batch = ledger[["accepted_species", "genus", "family", "n_islands", "n_records"]]
    candidates = pd.DataFrame(columns=campaign.COMMON_CANDIDATE_COLUMNS)

    result = campaign.apply_task_result(
        ledger,
        batch,
        candidates,
        ["openalex:Alpha one:transient:status=429 retry_after=missing"],
        set(),
        task,
        "wave_test",
        config,
    )

    row = result.iloc[0]
    assert row[f"{task}_status"] == "retry"
    assert int(row[f"{task}_attempts"]) == 0
    assert "status=429" in row[f"{task}_last_error"]


def test_permanent_failure_still_consumes_attempt_budget():
    config = _config()
    ledger = _ledger(config)
    task = "reproductive_openalex"
    batch = ledger[["accepted_species", "genus", "family", "n_islands", "n_records"]]
    candidates = pd.DataFrame(columns=campaign.COMMON_CANDIDATE_COLUMNS)

    result = campaign.apply_task_result(
        ledger,
        batch,
        candidates,
        ["openalex:Alpha one:permanent:status=400"],
        set(),
        task,
        "wave_test",
        config,
    )

    row = result.iloc[0]
    assert row[f"{task}_status"] == "retry"
    assert int(row[f"{task}_attempts"]) == 1
