import json
from pathlib import Path

import httpx
import pandas as pd
import pytest

import island_v2.local_openaire_reproductive_collection as openaire_module
from island_v2.local_openaire_reproductive_collection import (
    GROUP_LEDGER_COLUMNS,
    OpenAireTransientError,
    PROVIDER,
    QUERY_GROUPS,
    _collection_lock,
    _effective_min_interval,
    _file_sha256,
    _invalid_marker_path,
    _load_raw_receipt,
    _raw_cache_path,
    _receipt_path,
    _request_query,
    _write_raw_receipt,
    collect_openaire_batch,
    ensure_openaire_checkpoint_rows,
    ensure_query_group_ledger,
    openaire_sources_from_receipt,
)
from island_v2.local_reproductive_assurance_collection import (
    CHECKPOINT_COLUMNS,
    RA_TRAITS,
    build_reviewed_search_batch,
    materialize_reviewed_cache,
)


PAGE_SIZE = 2


def _queue(*species: str) -> pd.DataFrame:
    return pd.DataFrame({"accepted_species": list(species)})


def _base_checkpoint(*species: str) -> pd.DataFrame:
    rows = []
    for species in species:
        rows.append(
            {
                "species": species,
                "trait": "self_incompatibility",
                "provider": "europe_pmc",
                "provider_version": "europe_pmc_v2",
                "status": "provider_no_result",
                "terminal": True,
                "attempts": 1,
                "candidate_count": 0,
                "accepted_record_count": 0,
                "legacy_status": "no_hit",
                "legacy_enabled": "true",
                "last_error": "",
                "updated_at": "2026-07-17T00:00:00Z",
                "last_wave_id": "",
                "source_run_id": "29637034874",
                "migration_basis": "test",
                "contract_version": ("local_reproductive_assurance_provider_checkpoint_v2"),
            }
        )
    return pd.DataFrame(rows, columns=CHECKPOINT_COLUMNS)


def _result(
    species: str = "Alpha beta",
    *,
    suffix: str = "1",
    doi: str | None = "10.1000/example",
    description: str | None = None,
) -> dict:
    description = description or (
        f"In {species}, 346 seeds were obtained through "
        "self-pollination in a controlled breeding experiment."
    )
    result = {
        "id": f"openaire::{suffix}",
        "type": "publication",
        "mainTitle": f"Fertile clones of {species}",
        "descriptions": [description],
        "publicationDate": "1996-01-01",
        "publisher": "Example Society",
        "container": {"name": "Example Journal"},
        "pids": ([{"scheme": "doi", "value": doi}] if doi else []),
        "instances": [
            {
                "urls": [f"https://example.test/{suffix}"],
                "publicationDate": "1996-01-01",
            }
        ],
    }
    return result


def _payload(
    results: list[dict] | None = None,
    *,
    page: int = 1,
    page_size: int = PAGE_SIZE,
    num_found: int | None = None,
) -> dict:
    results = [] if results is None else results
    return {
        "header": {
            "numFound": len(results) if num_found is None else num_found,
            "maxScore": 1.0,
            "queryTime": 5,
            "page": page,
            "pageSize": page_size,
        },
        "results": results,
    }


def _write_fixture(
    raw: Path,
    species: str,
    query_group: str,
    payload: dict,
    *,
    page: int = 1,
    page_size: int = PAGE_SIZE,
    batch: str = "openaire_test",
):
    content = json.dumps(payload, separators=(",", ":")).encode()
    _write_raw_receipt(
        raw,
        species=species,
        query_group=query_group,
        page=page,
        retrieved_at="2026-07-19T03:00:00Z",
        batch_id=batch,
        page_size=page_size,
        authenticated=False,
        raw_content=content,
    )
    return _load_raw_receipt(raw, species, query_group, page, page_size)


def _axis_checkpoint(species: str) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "accepted_species": species,
                "genus": species.split()[0],
                "family": "Exampleaceae",
                "colour_vividness_evidence_traits": "flower_primary_color",
                "colour_vividness_direct": True,
                "floral_complexity_evidence_traits": "floral_form",
                "floral_complexity_direct": True,
                "reproductive_assurance_evidence_traits": "",
                "reproductive_assurance_direct": False,
                "all_three_axes_direct": False,
                "missing_direct_axes": "reproductive_assurance",
                "reproductive_assurance_priority": True,
                "provider_state_available": True,
                "checkpoint_scope": "test",
            }
        ]
    )


def _query_group(request: httpx.Request) -> str:
    query = request.url.params["search"]
    return next(
        group
        for group, terms in QUERY_GROUPS.items()
        if all(f'"{term}"' in query for term in terms)
    )


def test_three_query_groups_are_exact_and_have_at_most_four_operators():
    expected = {
        "compatibility": (
            '"Alpha beta" AND ("self-compatible" OR "self-incompatible" OR '
            '"self-pollination" OR "self-fertilization")'
        ),
        "mating": (
            '"Alpha beta" AND ("selfing" OR "autogamy" OR "outcrossing" OR "mating system")'
        ),
        "assurance": (
            '"Alpha beta" AND ("cleistogamy" OR "mixed mating" OR '
            '"autonomous selfing" OR "delayed selfing")'
        ),
    }
    assert {group: _request_query("Alpha beta", group) for group in QUERY_GROUPS} == expected
    for query in expected.values():
        assert query.count(" AND ") + query.count(" OR ") == 4


def test_external_checkpoint_and_internal_group_ledger_are_idempotent_and_recover():
    base = _base_checkpoint("Alpha beta", "Gamma delta", "Agave")
    external = ensure_openaire_checkpoint_rows(base, _queue("Alpha beta"))
    provider_rows = external.loc[external["provider"].eq(PROVIDER)]
    assert len(provider_rows) == 3 * len(RA_TRAITS)
    assert set(provider_rows.loc[provider_rows["species"].eq("Alpha beta"), "status"]) == {
        "pending"
    }
    assert set(provider_rows.loc[provider_rows["species"].eq("Gamma delta"), "status"]) == {
        "provider_seed_covered"
    }
    assert set(provider_rows.loc[provider_rows["species"].eq("Agave"), "status"]) == {
        "invalid_species_name"
    }

    ledger = ensure_query_group_ledger(
        pd.DataFrame(columns=GROUP_LEDGER_COLUMNS), _queue("Alpha beta"), PAGE_SIZE
    )
    assert len(ledger) == 3
    assert not ledger.duplicated(["species", "query_group"]).any()
    assert set(ledger["page_size"]) == {PAGE_SIZE}
    assert ledger["query_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    running = ledger["query_group"].eq("compatibility")
    ledger.loc[running, "status"] = "running"
    recovered = ensure_query_group_ledger(ledger, _queue("Alpha beta"), PAGE_SIZE)
    assert recovered.loc[recovered["query_group"].eq("compatibility"), "status"].item() == "retry"
    assert len(ensure_query_group_ledger(recovered, _queue("Alpha beta"), PAGE_SIZE)) == 3

    terminal = recovered["query_group"].eq("assurance")
    recovered.loc[terminal, "status"] = "group_no_result"
    recovered.loc[terminal, "terminal"] = True
    recovered.loc[terminal, "attempts"] = 7
    recovered.loc[terminal, "pages_completed"] = 1
    recovered.loc[terminal, "raw_response_path"] = "/tmp/openaire-terminal.response.json"
    recovered.loc[terminal, "raw_response_hash"] = "a" * 64
    immutable = ensure_query_group_ledger(recovered, _queue("Alpha beta"), PAGE_SIZE)
    assurance = immutable.loc[immutable["query_group"].eq("assurance")].iloc[0]
    assert assurance["status"] == "group_no_result"
    assert assurance["attempts"] == 7

    invalid = recovered.copy()
    invalid.loc[0, "query_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="query-group contract"):
        ensure_query_group_ledger(invalid, _queue("Alpha beta"), PAGE_SIZE)

    invalid_external = external.copy()
    invalid_key = invalid_external["provider"].eq(PROVIDER) & invalid_external["status"].eq(
        "pending"
    )
    invalid_external.loc[invalid_key, "terminal"] = True
    with pytest.raises(ValueError, match="external checkpoint terminal/status invariant"):
        ensure_openaire_checkpoint_rows(invalid_external, _queue("Alpha beta"))


def test_receipt_extracts_exact_dated_source_and_deduplicates_doi(tmp_path: Path):
    duplicate = _result(suffix="duplicate")
    payload = _payload([_result(), duplicate], num_found=2)
    raw = tmp_path / "raw.response.json"
    receipt = _write_fixture(raw, "Alpha beta", "compatibility", payload)
    sources, leads = openaire_sources_from_receipt(receipt, raw, PAGE_SIZE)
    assert len(sources) == 1
    assert sources.loc[0, "source_url"] == "https://doi.org/10.1000/example"
    assert sources.loc[0, "openaire_year"] == "1996"
    assert sources.loc[0, "openaire_publication"] == "Example Journal"
    assert sources.loc[0, "openaire_query_group"] == "compatibility"
    assert sources.loc[0, "raw_response_page"] == "1"
    assert sources.loc[0, "local_file_hash"] == _file_sha256(raw)
    assert set(leads["review_status"]) == {"unreviewed"}
    observed = set(leads[["trait", "provisional_value"]].itertuples(index=False, name=None))
    assert ("self_incompatibility", "SC") in observed


def test_exact_boundary_date_term_and_instance_url_filters_fail_closed(tmp_path: Path):
    suffix = _result(
        "Alpha beta2",
        suffix="suffix",
        doi="10.1000/suffix",
    )
    undated = _result(suffix="undated", doi="10.1000/undated")
    undated["publicationDate"] = None
    undated["instances"][0]["publicationDate"] = None
    no_term = _result(
        suffix="no-term",
        doi="10.1000/no-term",
        description="Alpha beta has blue flowers.",
    )
    no_term["mainTitle"] = "Flower colour of Alpha beta"
    invalid_date = _result(suffix="invalid-date", doi="10.1000/invalid-date")
    invalid_date["publicationDate"] = "1996-99-99"
    invalid_date["instances"][0]["publicationDate"] = "1996-99-99"
    no_publication = _result(suffix="no-publication", doi="10.1000/no-publication")
    no_publication["container"] = {}
    no_publication["publisher"] = None
    no_publication["sources"] = []
    payload = _payload(
        [suffix, undated, no_term, invalid_date, no_publication],
        page_size=10,
        num_found=5,
    )
    raw = tmp_path / "raw.response.json"
    receipt = _write_fixture(raw, "Alpha beta", "compatibility", payload, page_size=10)
    sources, leads = openaire_sources_from_receipt(receipt, raw, 10)
    assert sources.empty
    assert leads.empty

    fallback = _result(suffix="fallback", doi=None)
    payload = _payload([fallback], page_size=10)
    raw = tmp_path / "fallback.response.json"
    receipt = _write_fixture(raw, "Alpha beta", "compatibility", payload, page_size=10)
    sources, _leads = openaire_sources_from_receipt(receipt, raw, 10)
    assert sources.loc[0, "source_url"] == "https://example.test/fallback"
    assert sources.loc[0, "openaire_doi"] == ""

    invalid_doi = _result(suffix="invalid-doi", doi="not-a-doi")
    payload = _payload([invalid_doi], page_size=10)
    raw = tmp_path / "invalid-doi.response.json"
    receipt = _write_fixture(raw, "Alpha beta", "compatibility", payload, page_size=10)
    sources, _leads = openaire_sources_from_receipt(receipt, raw, 10)
    assert sources.loc[0, "source_url"] == "https://example.test/invalid-doi"
    assert sources.loc[0, "openaire_doi"] == ""


def test_non_object_result_member_fails_closed(tmp_path: Path):
    raw = tmp_path / "raw.response.json"
    receipt = _write_fixture(
        raw,
        "Alpha beta",
        "compatibility",
        _payload(["not-an-object"], num_found=1),
    )
    with pytest.raises(OpenAireTransientError, match="result_member_not_object"):
        openaire_sources_from_receipt(receipt, raw, PAGE_SIZE)


def test_rate_profile_requires_token_for_faster_than_sixty_seconds():
    assert _effective_min_interval("", None) == 60.1
    assert _effective_min_interval("real-token", None) == 0.51
    with pytest.raises(ValueError, match="60 requests/hour"):
        _effective_min_interval("", 59.9)
    with pytest.raises(ValueError, match="at least 0.5"):
        _effective_min_interval("real-token", 0.49)
    with pytest.raises(ValueError, match="60 requests/hour"):
        _effective_min_interval("", float("nan"))
    with pytest.raises(ValueError, match="at least 0.5"):
        _effective_min_interval("real-token", float("inf"))


def test_every_wire_retry_is_rate_limited_counted_and_checkpointed(tmp_path: Path):
    requests: list[tuple[str, int]] = []
    sleeps: list[float] = []
    compatibility_responses = 0

    def handler(request: httpx.Request) -> httpx.Response:
        nonlocal compatibility_responses
        group = _query_group(request)
        page = int(request.url.params["page"])
        requests.append((group, page))
        if group == "compatibility" and compatibility_responses == 0:
            compatibility_responses += 1
            return httpx.Response(503, request=request)
        return httpx.Response(
            200,
            json=_payload([], page=page, num_found=0),
            request=request,
        )

    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    client = httpx.Client(transport=httpx.MockTransport(handler))
    report = collect_openaire_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "openaire",
        max_species=1,
        page_size=PAGE_SIZE,
        min_interval_seconds=0.5,
        access_token="real-token",
        max_retries=1,
        backoff_seconds=0,
        client=client,
        sleeper=sleeps.append,
    )
    client.close()

    assert requests == [
        ("compatibility", 1),
        ("compatibility", 1),
        ("mating", 1),
        ("assurance", 1),
    ]
    assert report["network_queries"] == 4
    assert sum(delay > 0 for delay in sleeps) >= 3
    assert max(sleeps) >= 0.49
    ledger = pd.read_csv(tmp_path / "openaire" / "query_group_checkpoint.csv", dtype=str).fillna("")
    assert ledger.set_index("query_group")["attempts"].to_dict() == {
        "assurance": "1",
        "compatibility": "2",
        "mating": "1",
    }
    assert ledger["requested_at_utc"].str.endswith("Z").all()
    checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    assert set(checkpoint.loc[checkpoint["provider"].eq(PROVIDER), "attempts"]) == {"1"}


def test_paginated_collection_deduplicates_groups_and_terminalizes_after_all_three(
    tmp_path: Path,
):
    requests: list[tuple[str, int]] = []
    sleeps: list[float] = []

    def handler(request: httpx.Request) -> httpx.Response:
        group = _query_group(request)
        page = int(request.url.params["page"])
        requests.append((group, page))
        assert request.url.params["type"] == "publication"
        assert request.url.params["pageSize"] == str(PAGE_SIZE)
        assert "authorization" not in request.headers
        return httpx.Response(
            200,
            json=_payload([_result()], page=page, num_found=1),
            request=request,
        )

    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    _base_checkpoint("Alpha beta").to_csv(checkpoint_path, index=False)
    client = httpx.Client(transport=httpx.MockTransport(handler))
    report = collect_openaire_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=_base_checkpoint("Alpha beta"),
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "openaire",
        max_species=1,
        page_size=PAGE_SIZE,
        min_interval_seconds=None,
        max_retries=0,
        client=client,
        sleeper=sleeps.append,
    )
    client.close()
    assert requests == [(group, 1) for group in QUERY_GROUPS]
    assert len(sleeps) == 2
    assert all(delay > 59 for delay in sleeps)
    assert report["network_queries"] == 3
    assert report["source_records"] == 1
    assert report["authentication_mode"] == "unauthenticated"
    checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    rows = checkpoint.loc[checkpoint["provider"].eq(PROVIDER)]
    assert rows["terminal"].str.casefold().eq("true").all()
    ledger = pd.read_csv(tmp_path / "openaire" / "query_group_checkpoint.csv", dtype=str).fillna("")
    assert ledger["terminal"].str.casefold().eq("true").all()
    assert set(ledger["pages_completed"]) == {"1"}
    cumulative = pd.read_csv(tmp_path / "openaire" / "source_evidence.csv", dtype=str).fillna("")
    assert len(cumulative) == 1


def test_group_counts_are_post_dedup_across_pages(tmp_path: Path):
    def handler(request: httpx.Request) -> httpx.Response:
        group = _query_group(request)
        page = int(request.url.params["page"])
        if group != "compatibility":
            return httpx.Response(
                200,
                json=_payload([], page=page, num_found=0),
                request=request,
            )
        results = (
            [_result(suffix="one"), _result(suffix="two")]
            if page == 1
            else [_result(suffix="three")]
        )
        return httpx.Response(
            200,
            json=_payload(results, page=page, num_found=3),
            request=request,
        )

    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    client = httpx.Client(transport=httpx.MockTransport(handler))
    report = collect_openaire_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "openaire",
        max_species=1,
        page_size=PAGE_SIZE,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    client.close()
    assert report["source_records"] == 1
    ledger = pd.read_csv(tmp_path / "openaire" / "query_group_checkpoint.csv", dtype=str).fillna("")
    compatibility = ledger.loc[ledger["query_group"].eq("compatibility")].iloc[0]
    assert compatibility["result_count"] == "3"
    assert compatibility["source_count"] == "1"
    assert compatibility["candidate_count"] == "1"


def test_batch_root_outputs_exist_before_any_group_becomes_terminal(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    output_dir = tmp_path / "openaire"
    group_path = output_dir / "query_group_checkpoint.csv"
    real_write = openaire_module._atomic_write_csv
    observed_terminal_writes = 0

    def observing_write(frame: pd.DataFrame, path: Path) -> None:
        nonlocal observed_terminal_writes
        if (
            path == group_path
            and "terminal" in frame.columns
            and frame["terminal"].map(lambda value: str(value).casefold() == "true").any()
        ):
            observed_terminal_writes += 1
            batch_dirs = list((output_dir / "batches").iterdir())
            assert len(batch_dirs) == 1
            for name in (
                "source_evidence.csv",
                "candidate_leads.csv",
                "lookup_errors.csv",
            ):
                assert (batch_dirs[0] / name).exists()
        real_write(frame, path)

    monkeypatch.setattr(openaire_module, "_atomic_write_csv", observing_write)

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(
            200,
            json=_payload([], page=1, num_found=0),
            request=request,
        )

    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    client = httpx.Client(transport=httpx.MockTransport(handler))
    collect_openaire_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=output_dir,
        max_species=1,
        page_size=PAGE_SIZE,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    client.close()
    assert observed_terminal_writes >= 3


def test_partial_page_and_group_resume_never_requeries_completed_work(tmp_path: Path):
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    output_dir = tmp_path / "openaire"
    first_requests: list[tuple[str, int]] = []

    def first_handler(request: httpx.Request) -> httpx.Response:
        group = _query_group(request)
        page = int(request.url.params["page"])
        first_requests.append((group, page))
        if page == 2:
            return httpx.Response(503, request=request)
        return httpx.Response(
            200,
            json=_payload(
                [
                    _result(suffix="1", doi="10.1000/one"),
                    _result(suffix="2", doi="10.1000/two"),
                ],
                page=1,
                num_found=3,
            ),
            request=request,
        )

    first_client = httpx.Client(transport=httpx.MockTransport(first_handler))
    first = collect_openaire_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=output_dir,
        max_species=1,
        page_size=PAGE_SIZE,
        max_retries=0,
        client=first_client,
        sleeper=lambda _seconds: None,
    )
    first_client.close()
    assert first_requests == [("compatibility", 1), ("compatibility", 2)]
    assert first["transient_errors"] == 1
    checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    external = checkpoint.loc[checkpoint["provider"].eq(PROVIDER)]
    assert external["terminal"].str.casefold().eq("false").all()
    ledger = pd.read_csv(output_dir / "query_group_checkpoint.csv", dtype=str).fillna("")
    compatibility = ledger.loc[ledger["query_group"].eq("compatibility")].iloc[0]
    assert compatibility["status"] == "retry"
    assert compatibility["pages_completed"] == "1"
    assert set(ledger.loc[~ledger["query_group"].eq("compatibility"), "status"]) == {"pending"}

    second_requests: list[tuple[str, int]] = []

    def second_handler(request: httpx.Request) -> httpx.Response:
        group = _query_group(request)
        page = int(request.url.params["page"])
        second_requests.append((group, page))
        if group == "compatibility":
            return httpx.Response(
                200,
                json=_payload(
                    [_result(suffix="3", doi="10.1000/three")],
                    page=2,
                    num_found=3,
                ),
                request=request,
            )
        return httpx.Response(
            200,
            json=_payload([], page=1, num_found=0),
            request=request,
        )

    second_client = httpx.Client(transport=httpx.MockTransport(second_handler))
    second = collect_openaire_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=checkpoint,
        provider_checkpoint_path=checkpoint_path,
        output_dir=output_dir,
        max_species=1,
        page_size=PAGE_SIZE,
        max_retries=0,
        client=second_client,
        sleeper=lambda _seconds: None,
    )
    second_client.close()
    assert second_requests == [
        ("compatibility", 2),
        ("mating", 1),
        ("assurance", 1),
    ]
    assert second["remaining_openaire_species"] == 0
    result = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    assert result.loc[result["provider"].eq(PROVIDER), "terminal"].str.casefold().eq("true").all()
    ledger = pd.read_csv(output_dir / "query_group_checkpoint.csv", dtype=str).fillna("")
    assert ledger["terminal"].str.casefold().eq("true").all()
    assert ledger.loc[ledger["query_group"].eq("compatibility"), "pages_completed"].item() == "2"


def test_changed_num_found_resets_whole_group_and_restarts_at_page_one(tmp_path: Path):
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    output_dir = tmp_path / "openaire"
    first_requests: list[tuple[str, int]] = []

    def first_handler(request: httpx.Request) -> httpx.Response:
        group = _query_group(request)
        page = int(request.url.params["page"])
        first_requests.append((group, page))
        if page == 1:
            return httpx.Response(
                200,
                json=_payload(
                    [
                        _result(suffix="one", doi="10.1000/one"),
                        _result(suffix="two", doi="10.1000/two"),
                    ],
                    page=1,
                    num_found=3,
                ),
                request=request,
            )
        return httpx.Response(
            200,
            json=_payload(
                [
                    _result(suffix="three", doi="10.1000/three"),
                    _result(suffix="four", doi="10.1000/four"),
                ],
                page=2,
                num_found=4,
            ),
            request=request,
        )

    first_client = httpx.Client(transport=httpx.MockTransport(first_handler))
    first = collect_openaire_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=output_dir,
        max_species=1,
        page_size=PAGE_SIZE,
        max_retries=0,
        client=first_client,
        sleeper=lambda _seconds: None,
    )
    first_client.close()
    assert first_requests == [("compatibility", 1), ("compatibility", 2)]
    assert first["source_records"] == 0
    ledger = pd.read_csv(output_dir / "query_group_checkpoint.csv", dtype=str).fillna("")
    compatibility = ledger.loc[ledger["query_group"].eq("compatibility")].iloc[0]
    assert compatibility["status"] == "retry"
    assert compatibility["pages_completed"] == "0"
    assert compatibility["source_count"] == "0"
    assert compatibility["raw_response_path"] == ""
    assert pd.read_csv(output_dir / "source_evidence.csv").empty
    for page in (1, 2):
        raw = _raw_cache_path(output_dir / "raw", "Alpha beta", "compatibility", page, PAGE_SIZE)
        assert _invalid_marker_path(raw).exists()

    checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    second_requests: list[tuple[str, int]] = []

    def second_handler(request: httpx.Request) -> httpx.Response:
        group = _query_group(request)
        page = int(request.url.params["page"])
        second_requests.append((group, page))
        results = [_result()] if group == "compatibility" else []
        return httpx.Response(
            200,
            json=_payload(results, page=page, num_found=len(results)),
            request=request,
        )

    second_client = httpx.Client(transport=httpx.MockTransport(second_handler))
    second = collect_openaire_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=checkpoint,
        provider_checkpoint_path=checkpoint_path,
        output_dir=output_dir,
        max_species=1,
        page_size=PAGE_SIZE,
        max_retries=0,
        client=second_client,
        sleeper=lambda _seconds: None,
    )
    second_client.close()
    assert second_requests == [
        ("compatibility", 1),
        ("mating", 1),
        ("assurance", 1),
    ]
    assert second["remaining_openaire_species"] == 0


def test_interrupted_complete_raw_pages_are_reused_without_network_or_attempt(
    tmp_path: Path,
):
    output_dir = tmp_path / "openaire"
    raw_dir = output_dir / "raw"
    raw_dir.mkdir(parents=True)
    group_path = output_dir / "query_group_checkpoint.csv"
    group_ledger = ensure_query_group_ledger(
        pd.DataFrame(columns=GROUP_LEDGER_COLUMNS), _queue("Alpha beta"), PAGE_SIZE
    )
    group_ledger["status"] = "running"
    group_ledger["attempts"] = 1
    group_ledger["requested_at_utc"] = "2026-07-19T02:00:00Z"
    group_ledger.to_csv(group_path, index=False)
    for query_group in QUERY_GROUPS:
        raw = _raw_cache_path(raw_dir, "Alpha beta", query_group, 1, PAGE_SIZE)
        _write_fixture(
            raw,
            "Alpha beta",
            query_group,
            _payload([], page=1, num_found=0),
            batch=f"original_{query_group}",
        )

    external = ensure_openaire_checkpoint_rows(_base_checkpoint("Alpha beta"), _queue("Alpha beta"))
    mask = external["provider"].eq(PROVIDER)
    external.loc[mask, "status"] = "running"
    external.loc[mask, "terminal"] = False
    external.loc[mask, "attempts"] = 1
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    external.to_csv(checkpoint_path, index=False)

    def handler(_request: httpx.Request) -> httpx.Response:
        raise AssertionError("complete interrupted raw pages must be reused")

    client = httpx.Client(transport=httpx.MockTransport(handler))
    report = collect_openaire_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=external,
        provider_checkpoint_path=checkpoint_path,
        output_dir=output_dir,
        max_species=1,
        page_size=PAGE_SIZE,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    client.close()
    assert report["network_queries"] == 0
    assert report["cached_raw_responses_reused"] == 3
    result = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    external_rows = result.loc[result["provider"].eq(PROVIDER)]
    assert set(external_rows["attempts"]) == {"1"}
    assert external_rows["terminal"].str.casefold().eq("true").all()
    recovered_groups = pd.read_csv(group_path, dtype=str).fillna("")
    assert set(recovered_groups["attempts"]) == {"1"}
    assert recovered_groups["terminal"].str.casefold().eq("true").all()


@pytest.mark.parametrize("status", [400, 401, 403, 429, 503])
def test_provider_errors_stop_batch_and_leave_current_retry_later_pending(
    tmp_path: Path,
    status: int,
):
    requests = 0

    def handler(request: httpx.Request) -> httpx.Response:
        nonlocal requests
        requests += 1
        return httpx.Response(status, request=request)

    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta", "Gamma delta")
    base.to_csv(checkpoint_path, index=False)
    client = httpx.Client(transport=httpx.MockTransport(handler))
    collect_openaire_batch(
        queue=_queue("Alpha beta", "Gamma delta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "openaire",
        max_species=2,
        page_size=PAGE_SIZE,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    client.close()
    assert requests == 1
    checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    rows = checkpoint.loc[checkpoint["provider"].eq(PROVIDER)]
    assert set(rows.loc[rows["species"].eq("Alpha beta"), "status"]) == {"retry"}
    assert set(rows.loc[rows["species"].eq("Gamma delta"), "status"]) == {"pending"}


def test_authenticated_requests_use_bearer_without_persisting_token(tmp_path: Path):
    token = "real-secret-token-value"

    def handler(request: httpx.Request) -> httpx.Response:
        assert request.headers["authorization"] == f"Bearer {token}"
        return httpx.Response(
            200,
            json=_payload([], page=1, num_found=0),
            request=request,
        )

    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    client = httpx.Client(transport=httpx.MockTransport(handler))
    report = collect_openaire_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "openaire",
        max_species=1,
        page_size=PAGE_SIZE,
        min_interval_seconds=0.5,
        access_token=token,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    client.close()
    assert report["authentication_mode"] == "bearer_token"
    assert report["effective_min_interval_seconds"] == 0.5
    for path in (tmp_path / "openaire").rglob("*"):
        if path.is_file():
            assert token.encode() not in path.read_bytes()


def test_terminal_external_key_is_immutable(tmp_path: Path):
    base = ensure_openaire_checkpoint_rows(_base_checkpoint("Alpha beta"), _queue("Alpha beta"))
    accepted = base["provider"].eq(PROVIDER) & base["trait"].eq("self_incompatibility")
    base.loc[accepted, "status"] = "accepted_evidence"
    base.loc[accepted, "terminal"] = True
    base.loc[accepted, "accepted_record_count"] = 1
    base.loc[accepted, "attempts"] = 7
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base.to_csv(checkpoint_path, index=False)

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(
            200,
            json=_payload([], page=1, num_found=0),
            request=request,
        )

    client = httpx.Client(transport=httpx.MockTransport(handler))
    collect_openaire_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "openaire",
        max_species=1,
        page_size=PAGE_SIZE,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    client.close()
    result = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    row = result.loc[
        result["provider"].eq(PROVIDER) & result["trait"].eq("self_incompatibility")
    ].iloc[0]
    assert row["status"] == "accepted_evidence"
    assert row["accepted_record_count"] == "1"
    assert row["attempts"] == "7"


def test_malformed_200_keeps_group_and_external_retryable_and_caches_bytes(
    tmp_path: Path,
):
    raw_content = b'{"header":{"numFound":1},"results":['

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, content=raw_content, request=request)

    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    client = httpx.Client(transport=httpx.MockTransport(handler))
    collect_openaire_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "openaire",
        max_species=1,
        page_size=PAGE_SIZE,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    client.close()
    raw = next((tmp_path / "openaire" / "raw").glob("*.response.json"))
    assert raw.read_bytes() == raw_content
    assert next((tmp_path / "openaire" / "raw").glob("*.receipt.json"))
    checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    assert set(checkpoint.loc[checkpoint["provider"].eq(PROVIDER), "status"]) == {"retry"}


def test_strict_review_materialization_and_raw_tamper_revalidation(tmp_path: Path):
    def handler(request: httpx.Request) -> httpx.Response:
        group = _query_group(request)
        results = [_result()] if group == "compatibility" else []
        return httpx.Response(
            200,
            json=_payload(results, page=1, page_size=10),
            request=request,
        )

    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    output_dir = tmp_path / "openaire"
    client = httpx.Client(transport=httpx.MockTransport(handler))
    collect_openaire_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=output_dir,
        max_species=1,
        page_size=10,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    client.close()
    source_path = output_dir / "source_evidence.csv"
    source = pd.read_csv(source_path, dtype=str).fillna("")
    providers = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    excerpt = (
        "In Alpha beta, 346 seeds were obtained through "
        "self-pollination in a controlled breeding experiment."
    )
    selection = {
        "accepted": [
            {
                "source_row_index": 0,
                "accepted_species": "Alpha beta",
                "source_url": source.loc[0, "source_url"],
                "accepted_trait": "self_incompatibility",
                "accepted_value": "SC",
                "evidence_excerpt": excerpt,
                "review_reason": "Direct target-species self-pollination result.",
            }
        ]
    }
    candidates, review = build_reviewed_search_batch(
        queue=_queue("Alpha beta"),
        search_evidence=source,
        provider_checkpoint=providers,
        selection_payload=selection,
        local_file_path=str(source_path),
        local_file_hash=_file_sha256(source_path),
    )
    raw_path = Path(source.loc[0, "local_file_path"])
    result = materialize_reviewed_cache(
        queue=_queue("Alpha beta"),
        baseline_checkpoint=_axis_checkpoint("Alpha beta"),
        provider_checkpoint=providers,
        candidates=candidates,
        review_payload=review,
    )
    record = result["records"].iloc[0]
    assert record["provider"] == PROVIDER
    assert record["DOI"] == "10.1000/example"
    assert record["year"] == "1996"
    assert record["local_file_hash"] == _file_sha256(raw_path)

    raw_path.write_bytes(raw_path.read_bytes() + b" ")
    with pytest.raises(ValueError, match="raw-response provenance mismatch"):
        materialize_reviewed_cache(
            queue=_queue("Alpha beta"),
            baseline_checkpoint=_axis_checkpoint("Alpha beta"),
            provider_checkpoint=providers,
            candidates=candidates,
            review_payload=review,
        )


def test_receipt_hash_tamper_is_rejected(tmp_path: Path):
    raw = _raw_cache_path(tmp_path, "Alpha beta", "compatibility", 1, PAGE_SIZE)
    _write_fixture(
        raw,
        "Alpha beta",
        "compatibility",
        _payload([_result()], num_found=1),
    )
    raw.write_bytes(raw.read_bytes() + b" ")
    with pytest.raises(ValueError, match="hash mismatch"):
        _load_raw_receipt(raw, "Alpha beta", "compatibility", 1, PAGE_SIZE)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("authentication_mode", "cookie", "authentication mode mismatch"),
        ("response_path", "/tmp/wrong", "response path mismatch"),
        ("retrieved_at_utc", "2026-07-19T03:00:00+09:00", "timestamp is not UTC"),
        ("source_run_id", "", "source run id is blank"),
        ("http_status", "200", "HTTP status mismatch"),
        ("response_bytes", "2", "byte count mismatch"),
    ],
)
def test_receipt_contract_metadata_tamper_is_rejected(
    tmp_path: Path,
    field: str,
    value: str,
    message: str,
):
    raw = _raw_cache_path(tmp_path, "Alpha beta", "compatibility", 1, PAGE_SIZE)
    _write_fixture(
        raw,
        "Alpha beta",
        "compatibility",
        _payload([], num_found=0),
    )
    receipt_path = _receipt_path(raw)
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    receipt[field] = value
    receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        _load_raw_receipt(raw, "Alpha beta", "compatibility", 1, PAGE_SIZE)


def test_terminal_group_cache_tamper_hard_fails_without_recollection(tmp_path: Path):
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    output_dir = tmp_path / "openaire"

    def initial_handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(
            200,
            json=_payload([], page=1, num_found=0),
            request=request,
        )

    client = httpx.Client(transport=httpx.MockTransport(initial_handler))
    collect_openaire_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=output_dir,
        max_species=1,
        page_size=PAGE_SIZE,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    client.close()
    raw = _raw_cache_path(output_dir / "raw", "Alpha beta", "compatibility", 1, PAGE_SIZE)
    raw.write_bytes(raw.read_bytes() + b" ")
    checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")

    def forbidden_handler(_request: httpx.Request) -> httpx.Response:
        raise AssertionError("terminal OpenAIRE groups must never be recollected")

    forbidden = httpx.Client(transport=httpx.MockTransport(forbidden_handler))
    with pytest.raises(ValueError, match="hash mismatch"):
        collect_openaire_batch(
            queue=_queue("Alpha beta"),
            provider_checkpoint=checkpoint,
            provider_checkpoint_path=checkpoint_path,
            output_dir=output_dir,
            max_species=1,
            page_size=PAGE_SIZE,
            max_retries=0,
            client=forbidden,
            sleeper=lambda _seconds: None,
        )
    forbidden.close()
    assert len(list((output_dir / "batches").iterdir())) == 1


def test_output_directory_collection_lock_is_nonblocking(tmp_path: Path):
    output_dir = tmp_path / "openaire"
    with _collection_lock(output_dir):
        with pytest.raises(RuntimeError, match="another OpenAIRE collector"):
            with _collection_lock(output_dir):
                raise AssertionError("nested lock must not be acquired")
