import copy
import json
from pathlib import Path
from urllib.parse import unquote

import httpx
import pandas as pd
import pytest

from island_v2.local_doaj_reproductive_collection import (
    PAGINATION_CONTRACT,
    PROVIDER,
    PROVIDER_VERSION,
    _load_raw_receipt,
    _raw_cache_path,
    _request_query,
    _write_raw_receipt,
    collect_doaj_batch,
    doaj_sources_from_receipt,
    ensure_doaj_checkpoint_rows,
)
from island_v2.local_reproductive_assurance_collection import (
    CHECKPOINT_COLUMNS,
    RA_TRAITS,
    _file_sha256,
    build_reviewed_search_batch,
    materialize_reviewed_cache,
)


PAGE_SIZE = 10


def _queue(*species: str) -> pd.DataFrame:
    return pd.DataFrame({"accepted_species": list(species)})


def _base_checkpoint(*species: str) -> pd.DataFrame:
    rows = []
    for name in species:
        rows.append(
            {
                "species": name,
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


def _payload(species: str = "Alpha beta") -> dict:
    return {
        "total": 1,
        "page": 1,
        "pageSize": PAGE_SIZE,
        "results": [
            {
                "id": "doaj-record-1",
                "bibjson": {
                    "identifier": [
                        {"id": "10.1000/EXAMPLE", "type": "doi"},
                        {"id": "1234-5678", "type": "eissn"},
                    ],
                    "journal": {
                        "title": "Example Journal",
                        "publisher": "Example Society",
                    },
                    "year": "2024",
                    "link": [
                        {
                            "type": "fulltext",
                            "content_type": "pdf",
                            "url": "https://example.test/paper.pdf",
                        }
                    ],
                    "title": f"Breeding system of {species}",
                    "abstract": (f"{species} is self-incompatible and predominantly outcrossing."),
                },
            }
        ],
    }


def _fruit_payload(species: str = "Solanum melongena") -> dict:
    payload = _payload(species)
    record = payload["results"][0]
    record["id"] = "f96e38f102974b0485b1c4f1612e97ba"
    bibjson = record["bibjson"]
    bibjson["identifier"] = [
        {"id": "0373-5680", "type": "pissn"},
        {"id": "1851-7471", "type": "eissn"},
    ]
    bibjson["year"] = "2011"
    bibjson["journal"] = {
        "title": "Revista de la Sociedad Entomológica Argentina",
        "publisher": "Sociedad Entomológica Argentina",
    }
    bibjson["title"] = "Insects associated with eggplant flowers"
    bibjson["abstract"] = (
        f"This experiment on {species} compared spontaneous self pollination "
        "with pollination by biological agents. Of the Solanum melongena "
        "flowers submitted to the spontaneous self pollination test, 39% set "
        "fruits."
    )
    bibjson["link"] = [
        {
            "type": "fulltext",
            "url": "https://example.test/eggplant-fulltext",
        }
    ]
    return payload


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


def _write_fixture(raw: Path, species: str, payload: dict, batch: str = "doaj_test"):
    _write_raw_receipt(
        raw,
        species=species,
        retrieved_at="2026-07-19T01:00:00Z",
        batch_id=batch,
        payload=payload,
        page_size=PAGE_SIZE,
    )
    return _load_raw_receipt(raw, species, PAGE_SIZE)


def test_checkpoint_extension_is_idempotent_and_recovers_running_rows():
    base = _base_checkpoint("Alpha beta", "Gamma delta", "Agave")
    extended = ensure_doaj_checkpoint_rows(base, _queue("Alpha beta"))
    doaj = extended.loc[extended["provider"].eq(PROVIDER)]
    assert len(doaj) == 3 * len(RA_TRAITS)
    assert not extended.duplicated(["species", "trait", "provider"]).any()
    assert set(doaj.loc[doaj["species"].eq("Alpha beta"), "status"]) == {"pending"}
    assert set(doaj.loc[doaj["species"].eq("Gamma delta"), "status"]) == {"provider_seed_covered"}
    assert set(doaj.loc[doaj["species"].eq("Agave"), "status"]) == {"invalid_species_name"}
    running = extended["provider"].eq(PROVIDER) & extended["species"].eq("Alpha beta")
    extended.loc[running, "status"] = "running"
    recovered = ensure_doaj_checkpoint_rows(extended, _queue("Alpha beta"))
    rows = recovered.loc[running]
    assert set(rows["status"]) == {"retry"}
    assert rows["terminal"].eq(False).all()  # noqa: E712
    assert len(recovered) == len(extended)


def test_receipt_yields_dated_exact_species_unreviewed_leads(tmp_path: Path):
    species = "Crassula multicava"
    raw = tmp_path / "raw.response.json"
    receipt = _write_fixture(raw, species, _payload(species))
    sources, leads = doaj_sources_from_receipt(receipt, raw)
    assert len(sources) == 1
    assert sources.loc[0, "accepted_species"] == species
    assert sources.loc[0, "doaj_year"] == "2024"
    assert sources.loc[0, "doaj_record_id"] == "doaj-record-1"
    assert sources.loc[0, "source_url"] == "https://doi.org/10.1000/example"
    assert "DOI:10.1000/example" in sources.loc[0, "source_citation"]
    assert sources.loc[0, "evidence_scope"] == "species_direct"
    assert sources.loc[0, "local_file_hash"] == _file_sha256(raw)
    observed = set(leads[["trait", "provisional_value"]].itertuples(index=False, name=None))
    assert ("self_incompatibility", "SI") in observed
    assert ("mating_system", "predominantly_outcrossing") in observed
    assert set(leads["review_status"]) == {"unreviewed"}


def test_explicit_spontaneous_self_pollination_fruit_result_maps_two_traits(
    tmp_path: Path,
):
    species = "Solanum melongena"
    raw = tmp_path / "raw.response.json"
    receipt = _write_fixture(raw, species, _fruit_payload(species))
    sources, leads = doaj_sources_from_receipt(receipt, raw)
    assert len(sources) == 1
    observed = set(leads[["trait", "provisional_value"]].itertuples(index=False, name=None))
    assert ("self_incompatibility", "SC") in observed
    assert ("autonomous_selfing_capacity", "autonomous") in observed


def test_exact_species_date_term_and_url_filters_fail_closed(tmp_path: Path):
    payload = _payload()
    original = payload["results"][0]
    suffix = copy.deepcopy(original)
    suffix["id"] = "suffix"
    suffix["bibjson"]["title"] = "Breeding system of Alpha beta2"
    suffix["bibjson"]["abstract"] = "Alpha beta2 is self-incompatible."
    suffix["bibjson"]["identifier"] = [{"id": "10.1000/suffix", "type": "doi"}]
    undated = copy.deepcopy(original)
    undated["id"] = "undated"
    undated["bibjson"]["year"] = ""
    undated["created_date"] = "2024-01-01T00:00:00Z"
    no_term = copy.deepcopy(original)
    no_term["id"] = "no-term"
    no_term["bibjson"]["title"] = "Flower colour of Alpha beta"
    no_term["bibjson"]["abstract"] = "Alpha beta has blue flowers."
    no_url = copy.deepcopy(original)
    no_url["id"] = ""
    no_url["bibjson"]["identifier"] = []
    no_url["bibjson"]["link"] = []
    payload["results"] = [suffix, undated, no_term, no_url]
    payload["total"] = len(payload["results"])
    raw = tmp_path / "raw.response.json"
    receipt = _write_fixture(raw, "Alpha beta", payload)
    sources, leads = doaj_sources_from_receipt(receipt, raw)
    assert sources.empty
    assert leads.empty


def test_zero_results_are_valid_and_terminal_no_result(tmp_path: Path):
    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json={"total": 0, "results": []}, request=request)

    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    client = httpx.Client(transport=httpx.MockTransport(handler))
    report = collect_doaj_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "doaj",
        max_species=1,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    result = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    rows = result.loc[result["provider"].eq(PROVIDER)]
    assert set(rows["status"]) == {"provider_no_result"}
    assert rows["terminal"].str.casefold().eq("true").all()
    assert report["source_records"] == 0
    client.close()


def test_receipt_validates_request_fingerprint_exact_bytes_and_tamper(tmp_path: Path):
    raw = tmp_path / "raw.response.json"
    content = b'{"total":1, "results":[]}'
    kwargs = {
        "species": "Alpha beta",
        "retrieved_at": "2026-07-19T01:00:00Z",
        "batch_id": "doaj_test",
        "payload": {"total": 1, "results": []},
        "page_size": PAGE_SIZE,
        "raw_content": content,
    }
    _write_raw_receipt(raw, **kwargs)
    assert raw.read_bytes() == content
    receipt_path = raw.with_name(raw.name.removesuffix(".response.json") + ".receipt.json")
    metadata = json.loads(receipt_path.read_text(encoding="utf-8"))
    metadata["request"]["page_size"] = 20
    receipt_path.write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="page size mismatch"):
        _load_raw_receipt(raw, "Alpha beta", PAGE_SIZE)
    _write_raw_receipt(raw, **kwargs)
    metadata = json.loads(receipt_path.read_text(encoding="utf-8"))
    metadata["request"]["query"] = "different"
    receipt_path.write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="query mismatch"):
        _load_raw_receipt(raw, "Alpha beta", PAGE_SIZE)
    _write_raw_receipt(raw, **kwargs)
    metadata = json.loads(receipt_path.read_text(encoding="utf-8"))
    metadata["provider_version"] = "doaj_old"
    receipt_path.write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="provider_version mismatch"):
        _load_raw_receipt(raw, "Alpha beta", PAGE_SIZE)
    assert PROVIDER_VERSION.endswith("_v1")
    _write_raw_receipt(raw, **kwargs)
    raw.write_bytes(raw.read_bytes() + b" ")
    with pytest.raises(ValueError, match="response hash mismatch"):
        _load_raw_receipt(raw, "Alpha beta", PAGE_SIZE)


def test_collection_is_queue_bounded_rate_contract_resumable_and_writes_before_terminal(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    import island_v2.local_doaj_reproductive_collection as module

    requests = 0
    expected_query = _request_query("Alpha beta")

    def handler(request: httpx.Request) -> httpx.Response:
        nonlocal requests
        requests += 1
        encoded = request.url.path.split("/api/v4/search/articles/", 1)[1]
        assert unquote(encoded) == expected_query
        assert request.url.params["page"] == "1"
        assert request.url.params["pageSize"] == str(PAGE_SIZE)
        assert "authorization" not in request.headers
        return httpx.Response(200, json=_payload(), request=request)

    original_write = module._atomic_write_csv
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    output_dir = tmp_path / "doaj"

    def observing_write(frame: pd.DataFrame, path: Path, **kwargs):
        if path == checkpoint_path and "provider" in frame.columns:
            doaj = frame.loc[frame["provider"].eq(PROVIDER) & frame["species"].eq("Alpha beta")]
            if not doaj.empty and doaj["terminal"].map(module._bool).any():
                batches = list((output_dir / "batches").glob("*/source_evidence.csv"))
                assert batches, "terminal transition occurred before batch source output"
        return original_write(frame, path, **kwargs)

    monkeypatch.setattr(module, "_atomic_write_csv", observing_write)
    base = _base_checkpoint("Alpha beta", "Gamma delta")
    base.to_csv(checkpoint_path, index=False)
    client = httpx.Client(transport=httpx.MockTransport(handler))
    report = collect_doaj_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=output_dir,
        max_species=10,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    assert report["requested_species"] == 1
    assert report["scheduled_species"] == 1
    assert report["pagination_contract"] == PAGINATION_CONTRACT
    assert report["page_size"] == PAGE_SIZE
    assert report["biological_values_materialized"] is False
    assert requests == 1
    checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    alpha = checkpoint.loc[
        checkpoint["provider"].eq(PROVIDER) & checkpoint["species"].eq("Alpha beta")
    ]
    gamma = checkpoint.loc[
        checkpoint["provider"].eq(PROVIDER) & checkpoint["species"].eq("Gamma delta")
    ]
    assert alpha["terminal"].str.casefold().eq("true").all()
    assert set(alpha["attempts"]) == {"1"}
    assert set(gamma["status"]) == {"provider_seed_covered"}
    second = collect_doaj_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=checkpoint,
        provider_checkpoint_path=checkpoint_path,
        output_dir=output_dir,
        max_species=10,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    assert second["requested_species"] == 0
    assert requests == 1
    client.close()


@pytest.mark.parametrize(
    ("status", "payload"),
    [
        (408, None),
        (425, None),
        (429, None),
        (500, None),
        (503, None),
        (599, None),
        (403, None),
        (200, {"total": 1, "results": "not-a-list"}),
        (200, {"total": 1}),
        (200, {"total": "0", "results": []}),
    ],
)
def test_http_and_schema_failures_remain_retryable_and_cache_malformed_200(
    tmp_path: Path, status: int, payload: dict | None
):
    def handler(request: httpx.Request) -> httpx.Response:
        if payload is None:
            return httpx.Response(status, request=request)
        return httpx.Response(status, json=payload, request=request)

    client = httpx.Client(transport=httpx.MockTransport(handler))
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    report = collect_doaj_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "doaj",
        max_species=1,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    rows = checkpoint.loc[checkpoint["provider"].eq(PROVIDER)]
    assert set(rows["status"]) == {"retry"}
    assert rows["terminal"].str.casefold().eq("false").all()
    assert set(rows["attempts"]) == {"1"}
    assert report["transient_errors"] == 1
    if status == 200:
        assert len(list((tmp_path / "doaj" / "raw").glob("*.response.json"))) == 1
        assert len(list((tmp_path / "doaj" / "raw").glob("*.receipt.json"))) == 1
    client.close()


def test_invalid_json_200_is_cached_as_exact_bytes_before_retry(tmp_path: Path):
    raw_content = b'{"total":1,"results":['

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, content=raw_content, request=request)

    client = httpx.Client(transport=httpx.MockTransport(handler))
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    report = collect_doaj_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "doaj",
        max_species=1,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    raw = next((tmp_path / "doaj" / "raw").glob("*.response.json"))
    receipt = next((tmp_path / "doaj" / "raw").glob("*.receipt.json"))
    metadata = json.loads(receipt.read_text(encoding="utf-8"))
    assert raw.read_bytes() == raw_content
    assert metadata["response_sha256"] == _file_sha256(raw)
    assert metadata["response_bytes"] == len(raw_content)
    assert report["transient_errors"] == 1
    checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    assert set(checkpoint.loc[checkpoint["provider"].eq(PROVIDER), "status"]) == {"retry"}
    client.close()


def test_malformed_complete_cache_is_preserved_then_refreshed(tmp_path: Path):
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    output_dir = tmp_path / "doaj"

    def malformed(request: httpx.Request) -> httpx.Response:
        return httpx.Response(
            200,
            json={"total": 1, "results": "not-a-list"},
            request=request,
        )

    first_client = httpx.Client(transport=httpx.MockTransport(malformed))
    collect_doaj_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=output_dir,
        max_species=1,
        min_interval_seconds=0,
        max_retries=0,
        client=first_client,
        sleeper=lambda _seconds: None,
    )
    first_client.close()
    first_checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    requests = 0

    def valid(request: httpx.Request) -> httpx.Response:
        nonlocal requests
        requests += 1
        return httpx.Response(200, json=_payload(), request=request)

    second_client = httpx.Client(transport=httpx.MockTransport(valid))
    second = collect_doaj_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=first_checkpoint,
        provider_checkpoint_path=checkpoint_path,
        output_dir=output_dir,
        max_species=1,
        min_interval_seconds=0,
        max_retries=0,
        client=second_client,
        sleeper=lambda _seconds: None,
    )
    second_client.close()
    result = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    rows = result.loc[result["provider"].eq(PROVIDER)]
    assert requests == 1
    assert second["cached_raw_responses_reused"] == 0
    assert list((output_dir / "raw" / "invalid").glob("*.invalid"))
    assert set(rows["attempts"]) == {"2"}
    assert rows["terminal"].str.casefold().eq("true").all()


@pytest.mark.parametrize("status", [403, 429])
def test_provider_wide_4xx_stops_batch_and_leaves_later_species_pending(
    tmp_path: Path,
    status: int,
):
    requests = 0

    def handler(request: httpx.Request) -> httpx.Response:
        nonlocal requests
        requests += 1
        return httpx.Response(status, request=request)

    client = httpx.Client(transport=httpx.MockTransport(handler))
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta", "Gamma delta")
    base.to_csv(checkpoint_path, index=False)
    report = collect_doaj_batch(
        queue=_queue("Alpha beta", "Gamma delta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "doaj",
        max_species=2,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    rows = checkpoint.loc[checkpoint["provider"].eq(PROVIDER)]
    assert requests == 1
    assert report["requested_species"] == 1
    assert set(rows.loc[rows["species"].eq("Alpha beta"), "status"]) == {"retry"}
    assert set(rows.loc[rows["species"].eq("Gamma delta"), "status"]) == {"pending"}
    client.close()


def test_partial_terminal_key_is_never_overwritten(tmp_path: Path):
    checkpoint = ensure_doaj_checkpoint_rows(_base_checkpoint("Alpha beta"), _queue("Alpha beta"))
    accepted = (
        checkpoint["provider"].eq(PROVIDER)
        & checkpoint["species"].eq("Alpha beta")
        & checkpoint["trait"].eq("self_incompatibility")
    )
    checkpoint.loc[accepted, "status"] = "accepted_evidence"
    checkpoint.loc[accepted, "terminal"] = True
    checkpoint.loc[accepted, "accepted_record_count"] = 1
    checkpoint.loc[accepted, "attempts"] = 7
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    checkpoint.to_csv(checkpoint_path, index=False)

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json=_payload(), request=request)

    client = httpx.Client(transport=httpx.MockTransport(handler))
    collect_doaj_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=checkpoint,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "doaj",
        max_species=1,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    result = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    row = result.loc[
        result["provider"].eq(PROVIDER)
        & result["species"].eq("Alpha beta")
        & result["trait"].eq("self_incompatibility")
    ].iloc[0]
    assert row["status"] == "accepted_evidence"
    assert row["terminal"].casefold() == "true"
    assert row["accepted_record_count"] == "1"
    assert row["attempts"] == "7"
    client.close()


def test_interrupted_raw_is_reused_with_original_provenance_and_no_attempt(
    tmp_path: Path,
):
    checkpoint = ensure_doaj_checkpoint_rows(_base_checkpoint("Alpha beta"), _queue("Alpha beta"))
    mask = checkpoint["provider"].eq(PROVIDER)
    checkpoint.loc[mask, "status"] = "running"
    checkpoint.loc[mask, "terminal"] = False
    checkpoint.loc[mask, "attempts"] = 1
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    checkpoint.to_csv(checkpoint_path, index=False)
    raw_dir = tmp_path / "doaj" / "raw"
    raw_dir.mkdir(parents=True)
    raw = _raw_cache_path(raw_dir, "Alpha beta", PAGE_SIZE)
    _write_fixture(raw, "Alpha beta", _payload(), batch="doaj_interrupted")

    def handler(_request: httpx.Request) -> httpx.Response:
        raise AssertionError("valid interrupted raw cache must be reused")

    client = httpx.Client(transport=httpx.MockTransport(handler))
    report = collect_doaj_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=checkpoint,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "doaj",
        max_species=1,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    assert report["cached_raw_responses_reused"] == 1
    result = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    recovered = result.loc[result["provider"].eq(PROVIDER)]
    assert set(recovered["attempts"]) == {"1"}
    assert set(recovered["source_run_id"]) == {"doaj_interrupted"}
    assert recovered["terminal"].str.casefold().eq("true").all()
    client.close()


def test_incomplete_raw_is_archived_before_fresh_auditable_attempt(tmp_path: Path):
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    raw_dir = tmp_path / "doaj" / "raw"
    raw_dir.mkdir(parents=True)
    raw = _raw_cache_path(raw_dir, "Alpha beta", PAGE_SIZE)
    raw.write_bytes(b"partial-crash-body")
    requests = 0

    def handler(request: httpx.Request) -> httpx.Response:
        nonlocal requests
        requests += 1
        return httpx.Response(200, json=_payload(), request=request)

    client = httpx.Client(transport=httpx.MockTransport(handler))
    collect_doaj_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "doaj",
        max_species=1,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    assert requests == 1
    assert list((raw_dir / "incomplete").glob("*.incomplete"))
    assert _load_raw_receipt(raw, "Alpha beta", PAGE_SIZE)["payload"] == _payload()
    result = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    assert set(result.loc[result["provider"].eq(PROVIDER), "attempts"]) == {"1"}
    client.close()


def test_source_without_mapping_completes_pass_but_creates_no_lead(tmp_path: Path):
    payload = _payload()
    payload["results"][0]["bibjson"]["abstract"] = (
        "Alpha beta was included in a breeding system survey; no categorical result was reported."
    )

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json=payload, request=request)

    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    client = httpx.Client(transport=httpx.MockTransport(handler))
    report = collect_doaj_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "doaj",
        max_species=1,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    result = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    rows = result.loc[result["provider"].eq(PROVIDER)]
    assert report["source_records"] == 1
    assert report["candidate_leads"] == 0
    assert set(rows["status"]) == {"provider_pass_complete"}
    client.close()


def test_strict_review_materialization_preserves_raw_provenance_and_revalidates(
    tmp_path: Path,
):
    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json=_payload(), request=request)

    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    output_dir = tmp_path / "doaj"
    client = httpx.Client(transport=httpx.MockTransport(handler))
    report = collect_doaj_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=output_dir,
        max_species=1,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    client.close()
    source_path = output_dir / "batches" / report["batch_id"] / "source_evidence.csv"
    source = pd.read_csv(source_path, dtype=str).fillna("")
    providers = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    excerpt = "Alpha beta is self-incompatible and predominantly outcrossing."
    selection = {
        "accepted": [
            {
                "source_row_index": 0,
                "accepted_species": "Alpha beta",
                "source_url": source.loc[0, "source_url"],
                "accepted_trait": "self_incompatibility",
                "accepted_value": "SI",
                "evidence_excerpt": excerpt,
                "review_reason": "Exact target-species statement explicitly reports SI.",
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
    assert candidates.loc[0, "provider"] == PROVIDER
    assert candidates.loc[0, "local_file_path"] == str(raw_path.resolve())
    assert candidates.loc[0, "local_file_hash"] == _file_sha256(raw_path)
    result = materialize_reviewed_cache(
        queue=_queue("Alpha beta"),
        baseline_checkpoint=_axis_checkpoint("Alpha beta"),
        provider_checkpoint=providers,
        candidates=candidates,
        review_payload=review,
    )
    record = result["records"].iloc[0]
    assert record["provider"] == PROVIDER
    assert record["evidence_excerpt"] == excerpt
    assert record["local_file_path"] == str(raw_path.resolve())
    assert record["local_file_hash"] == _file_sha256(raw_path)
    assert record["DOI"] == "10.1000/example"
    assert record["year"] == "2024"

    ambiguous = {
        "accepted": [{**selection["accepted"][0], "evidence_excerpt": "self-incompatible"}]
    }
    ambiguous_candidates, ambiguous_review = build_reviewed_search_batch(
        queue=_queue("Alpha beta"),
        search_evidence=source,
        provider_checkpoint=providers,
        selection_payload=ambiguous,
        local_file_path=str(source_path),
        local_file_hash=_file_sha256(source_path),
    )
    with pytest.raises(ValueError, match="must name the exact species"):
        materialize_reviewed_cache(
            queue=_queue("Alpha beta"),
            baseline_checkpoint=_axis_checkpoint("Alpha beta"),
            provider_checkpoint=providers,
            candidates=ambiguous_candidates,
            review_payload=ambiguous_review,
        )

    raw_path.write_bytes(raw_path.read_bytes() + b" ")
    with pytest.raises(ValueError, match="raw-response provenance mismatch"):
        materialize_reviewed_cache(
            queue=_queue("Alpha beta"),
            baseline_checkpoint=_axis_checkpoint("Alpha beta"),
            provider_checkpoint=providers,
            candidates=candidates,
            review_payload=review,
        )
