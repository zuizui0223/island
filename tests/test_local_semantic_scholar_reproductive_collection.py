import json
from pathlib import Path

import httpx
import pandas as pd
import pytest

from island_v2.local_reproductive_assurance_collection import (
    CHECKPOINT_COLUMNS,
    RA_TRAITS,
    _file_sha256,
    build_reviewed_search_batch,
    materialize_reviewed_cache,
)
from island_v2.local_semantic_scholar_reproductive_collection import (
    FIELDS_PARAM,
    PROVIDER,
    PROVIDER_VERSION,
    _load_raw_receipt,
    _raw_cache_path,
    _request_query,
    _write_raw_receipt,
    collect_semantic_scholar_batch,
    ensure_semantic_scholar_checkpoint_rows,
    semantic_scholar_sources_from_receipt,
)


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
                "contract_version": "local_reproductive_assurance_provider_checkpoint_v2",
            }
        )
    return pd.DataFrame(rows, columns=CHECKPOINT_COLUMNS)


def _payload(species: str = "Alpha beta") -> dict:
    return {
        "total": 1,
        "offset": 0,
        "data": [
            {
                "paperId": "paper-1",
                "title": f"Breeding system of {species}",
                "abstract": (
                    f"{species} is self-incompatible and predominantly outcrossing."
                ),
                "year": 2024,
                "venue": "Example Journal",
                "url": "https://www.semanticscholar.org/paper/paper-1",
                "externalIds": {"DOI": "10.1000/example"},
                "publicationDate": "2024-04-02",
                "openAccessPdf": {"url": "https://example.test/paper.pdf"},
                "journal": {"name": "Example Journal", "volume": "1", "pages": "1-2"},
            }
        ],
    }


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


def test_checkpoint_extension_adds_all_original_species_four_keys_once():
    base = _base_checkpoint("Alpha beta", "Gamma delta", "Agave")
    extended = ensure_semantic_scholar_checkpoint_rows(base, _queue("Alpha beta"))
    semantic = extended.loc[extended["provider"].eq(PROVIDER)]
    assert len(semantic) == 3 * len(RA_TRAITS)
    assert not extended.duplicated(["species", "trait", "provider"]).any()
    assert set(semantic.loc[semantic["species"].eq("Alpha beta"), "status"]) == {"pending"}
    assert set(semantic.loc[semantic["species"].eq("Gamma delta"), "status"]) == {
        "provider_seed_covered"
    }
    assert set(semantic.loc[semantic["species"].eq("Agave"), "status"]) == {
        "invalid_species_name"
    }
    repeated = ensure_semantic_scholar_checkpoint_rows(extended, _queue("Alpha beta"))
    assert len(repeated) == len(extended)


def test_receipt_yields_dated_exact_unicode_species_review_leads(tmp_path: Path):
    species = "Crassula multicava"
    raw = tmp_path / "raw.response.json"
    _write_raw_receipt(
        raw,
        species=species,
        retrieved_at="2026-07-19T01:00:00Z",
        batch_id="semantic_test",
        payload=_payload(species),
    )
    receipt = _load_raw_receipt(raw, species)
    sources, leads = semantic_scholar_sources_from_receipt(receipt, raw)
    assert len(sources) == 1
    assert sources.loc[0, "accepted_species"] == species
    assert sources.loc[0, "semantic_scholar_year"] == "2024"
    assert sources.loc[0, "semantic_scholar_publication_date"] == "2024-04-02"
    assert sources.loc[0, "source_url"] == "https://doi.org/10.1000/example"
    assert "DOI:10.1000/example" in sources.loc[0, "source_citation"]
    assert sources.loc[0, "evidence_scope"] == "species_direct"
    assert sources.loc[0, "local_file_hash"] == _file_sha256(raw)
    observed = set(leads[["trait", "provisional_value"]].itertuples(index=False, name=None))
    assert ("self_incompatibility", "SI") in observed
    assert ("mating_system", "predominantly_outcrossing") in observed
    assert set(leads["review_status"]) == {"unreviewed"}


def test_exact_species_postfilter_and_publication_date_are_fail_closed(tmp_path: Path):
    payload = _payload()
    suffix = {**payload["data"][0]}
    suffix["title"] = "Breeding system of Alpha beta2"
    suffix["abstract"] = "Alpha beta2 is self-incompatible."
    suffix["paperId"] = "paper-suffix"
    suffix["externalIds"] = {"DOI": "10.1000/suffix"}
    undated = {**payload["data"][0]}
    undated["paperId"] = "paper-undated"
    undated["externalIds"] = {"DOI": "10.1000/undated"}
    undated["year"] = None
    undated["publicationDate"] = None
    payload["data"] = [suffix, undated]
    raw = tmp_path / "raw.response.json"
    _write_raw_receipt(
        raw,
        species="Alpha beta",
        retrieved_at="2026-07-19T01:00:00Z",
        batch_id="semantic_test",
        payload=payload,
    )
    receipt = _load_raw_receipt(raw, "Alpha beta")
    sources, leads = semantic_scholar_sources_from_receipt(receipt, raw)
    assert sources.empty
    assert leads.empty


def test_total_zero_without_data_is_a_valid_empty_first_page(tmp_path: Path):
    raw = tmp_path / "raw.response.json"
    _write_raw_receipt(
        raw,
        species="Alpha beta",
        retrieved_at="2026-07-19T01:00:00Z",
        batch_id="semantic_test",
        payload={"total": 0},
    )
    receipt = _load_raw_receipt(raw, "Alpha beta")
    sources, leads = semantic_scholar_sources_from_receipt(receipt, raw)
    assert sources.empty
    assert leads.empty


def test_raw_receipt_validates_request_fingerprint_and_exact_body_hash(tmp_path: Path):
    raw = tmp_path / "raw.response.json"
    content = b'{"offset":0, "total":1, "data":[]}'
    _write_raw_receipt(
        raw,
        species="Alpha beta",
        retrieved_at="2026-07-19T01:00:00Z",
        batch_id="semantic_test",
        payload={"offset": 0, "total": 1, "data": []},
        raw_content=content,
    )
    assert raw.read_bytes() == content
    receipt_path = raw.with_name(raw.name.removesuffix(".response.json") + ".receipt.json")
    metadata = json.loads(receipt_path.read_text(encoding="utf-8"))
    metadata["request"]["pagination_contract"] = "follow_all_pages"
    receipt_path.write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="pagination contract mismatch"):
        _load_raw_receipt(raw, "Alpha beta")
    _write_raw_receipt(
        raw,
        species="Alpha beta",
        retrieved_at="2026-07-19T01:00:00Z",
        batch_id="semantic_test",
        payload={"offset": 0, "total": 1, "data": []},
        raw_content=content,
    )
    metadata = json.loads(receipt_path.read_text(encoding="utf-8"))
    metadata["provider_version"] = "semantic_scholar_graph_bulk_exact_phrase_v2"
    receipt_path.write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="provider_version mismatch"):
        _load_raw_receipt(raw, "Alpha beta")
    assert PROVIDER_VERSION.endswith("_v3")
    _write_raw_receipt(
        raw,
        species="Alpha beta",
        retrieved_at="2026-07-19T01:00:00Z",
        batch_id="semantic_test",
        payload={"offset": 0, "total": 1, "data": []},
        raw_content=content,
    )
    raw.write_bytes(raw.read_bytes() + b" ")
    with pytest.raises(ValueError, match="response hash mismatch"):
        _load_raw_receipt(raw, "Alpha beta")


def test_collection_is_bounded_current_queue_only_and_resumable(tmp_path: Path):
    requests = 0
    expected_query = (
        '"Alpha beta" + ((self pollination) | (self fertilization) | '
        "(self compatibility) | (self incompatibility) | selfing | autogamy | "
        "outcrossing | (mating system) | (breeding system) | (mixed mating) | "
        "(delayed selfing) | cleistogamy)"
    )
    assert _request_query("Alpha beta") == expected_query

    def handler(request: httpx.Request) -> httpx.Response:
        nonlocal requests
        requests += 1
        assert request.url.path.endswith("/graph/v1/paper/search/bulk")
        assert "limit" not in request.url.params
        assert "token" not in request.url.params
        assert request.url.params["fields"] == FIELDS_PARAM
        assert request.url.params["query"] == expected_query
        assert "x-api-key" not in request.headers
        payload = _payload()
        payload["token"] = "ignored-next-page-token"
        return httpx.Response(200, json=payload, request=request)

    client = httpx.Client(transport=httpx.MockTransport(handler))
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta", "Gamma delta")
    base.to_csv(checkpoint_path, index=False)
    output_dir = tmp_path / "semantic"
    report = collect_semantic_scholar_batch(
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
    assert report["candidate_species"] == 1
    assert report["pagination_contract"] == "first_bulk_response_only_no_token_follow"
    assert report["biological_values_materialized"] is False
    assert requests == 1
    assert (output_dir / "source_evidence.csv").exists()
    batch = output_dir / "batches" / report["batch_id"]
    assert (batch / "source_evidence.csv").exists()
    assert (batch / "candidate_leads.csv").exists()
    semantic = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    semantic = semantic.loc[semantic["provider"].eq(PROVIDER)]
    alpha = semantic.loc[semantic["species"].eq("Alpha beta")]
    gamma = semantic.loc[semantic["species"].eq("Gamma delta")]
    assert alpha["terminal"].str.casefold().eq("true").all()
    assert set(alpha["attempts"]) == {"1"}
    assert set(gamma["status"]) == {"provider_seed_covered"}

    second = collect_semantic_scholar_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=pd.read_csv(checkpoint_path, dtype=str).fillna(""),
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
    ("status", "headers", "payload"),
    [
        (429, {}, None),
        (500, {}, None),
        (403, {}, None),
        (200, {}, {"total": 1, "data": "not-a-list"}),
        (200, {}, {"total": 1}),
        (200, {}, {"total": "0"}),
    ],
)
def test_http_and_schema_failures_remain_retryable(
    tmp_path: Path,
    status: int,
    headers: dict[str, str],
    payload: dict | None,
):
    def handler(request: httpx.Request) -> httpx.Response:
        if payload is None:
            return httpx.Response(status, headers=headers, request=request)
        return httpx.Response(status, json=payload, headers=headers, request=request)

    client = httpx.Client(transport=httpx.MockTransport(handler))
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    report = collect_semantic_scholar_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "semantic",
        max_species=1,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    semantic = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    semantic = semantic.loc[semantic["provider"].eq(PROVIDER)]
    assert set(semantic["status"]) == {"retry"}
    assert semantic["terminal"].str.casefold().eq("false").all()
    assert set(semantic["attempts"]) == {"1"}
    assert report["transient_errors"] == 1
    if status == 200:
        assert len(list((tmp_path / "semantic" / "raw").glob("*.response.json"))) == 1
        assert len(list((tmp_path / "semantic" / "raw").glob("*.receipt.json"))) == 1
    client.close()


def test_long_retry_after_stops_batch_without_terminalizing_other_species(tmp_path: Path):
    requests = 0

    def handler(request: httpx.Request) -> httpx.Response:
        nonlocal requests
        requests += 1
        return httpx.Response(429, headers={"Retry-After": "120"}, request=request)

    client = httpx.Client(transport=httpx.MockTransport(handler))
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta", "Gamma delta")
    base.to_csv(checkpoint_path, index=False)
    report = collect_semantic_scholar_batch(
        queue=_queue("Alpha beta", "Gamma delta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "semantic",
        max_species=2,
        min_interval_seconds=0,
        max_retries=2,
        client=client,
        sleeper=lambda _seconds: None,
    )
    semantic = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    semantic = semantic.loc[semantic["provider"].eq(PROVIDER)]
    assert requests == 1
    assert report["requested_species"] == 1
    assert set(semantic.loc[semantic["species"].eq("Alpha beta"), "status"]) == {"retry"}
    assert set(semantic.loc[semantic["species"].eq("Gamma delta"), "status"]) == {"pending"}
    client.close()


def test_partial_terminal_key_is_never_overwritten(tmp_path: Path):
    checkpoint = ensure_semantic_scholar_checkpoint_rows(
        _base_checkpoint("Alpha beta"), _queue("Alpha beta")
    )
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
    collect_semantic_scholar_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=checkpoint,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "semantic",
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


def test_interrupted_raw_is_reused_without_request_or_attempt_increment(tmp_path: Path):
    checkpoint = ensure_semantic_scholar_checkpoint_rows(
        _base_checkpoint("Alpha beta"), _queue("Alpha beta")
    )
    semantic = checkpoint["provider"].eq(PROVIDER)
    checkpoint.loc[semantic, "status"] = "running"
    checkpoint.loc[semantic, "terminal"] = False
    checkpoint.loc[semantic, "attempts"] = 1
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    checkpoint.to_csv(checkpoint_path, index=False)
    raw_dir = tmp_path / "semantic" / "raw"
    raw_dir.mkdir(parents=True)
    _write_raw_receipt(
        _raw_cache_path(raw_dir, "Alpha beta"),
        species="Alpha beta",
        retrieved_at="2026-07-19T01:00:00Z",
        batch_id="semantic_interrupted",
        payload=_payload(),
    )

    def handler(_request: httpx.Request) -> httpx.Response:
        raise AssertionError("valid interrupted raw cache must be reused")

    client = httpx.Client(transport=httpx.MockTransport(handler))
    report = collect_semantic_scholar_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=checkpoint,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "semantic",
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
    assert recovered["terminal"].str.casefold().eq("true").all()
    client.close()


def test_strict_review_materialization_preserves_raw_response_provenance(tmp_path: Path):
    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json=_payload(), request=request)

    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    output_dir = tmp_path / "semantic"
    client = httpx.Client(transport=httpx.MockTransport(handler))
    report = collect_semantic_scholar_batch(
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
