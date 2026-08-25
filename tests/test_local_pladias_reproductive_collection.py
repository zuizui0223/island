import json
from pathlib import Path
from urllib.parse import unquote

import httpx
import pandas as pd
import pytest

import island_v2.local_pladias_reproductive_collection as pladias_module
from island_v2.local_pladias_reproductive_collection import (
    EXACT_PAGE_CONTRACT,
    FIELD_LABEL,
    PROVIDER,
    PROVIDER_VERSION,
    PUBLICATION_YEAR,
    _category_mappings,
    _extract_exact_field,
    _load_raw_receipt,
    _request_url,
    _write_raw_receipt,
    collect_pladias_batch,
    ensure_pladias_checkpoint_rows,
    pladias_sources_from_receipt,
)
from island_v2.local_reproductive_assurance_collection import (
    CHECKPOINT_COLUMNS,
    RA_TRAITS,
    _explicit_mapping_supported,
    _file_sha256,
    _pladias_field_mapping_supported,
    build_reviewed_search_batch,
    materialize_reviewed_cache,
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


def _page(species: str = "Alpha beta", category: str = "facultative allogamy") -> bytes:
    return f"""<!doctype html>
    <html><head><title>{species} • Pladias</title></head><body>
      <h1><i>{species}</i></h1>
      <li class="list-group-item featureDetail">
        <div class="d-flex justify-content-between align-items-center">
          <div><span class="font-italic">{FIELD_LABEL}</span><span>:</span>
          <b>{category}</b></div>
        </div>
        <div class="modal"><h3>{FIELD_LABEL}</h3>
          <p>Definitions mention allogamy, autogamy, self-incompatibility,
          and mixed mating, but this prose is not the taxon value.</p>
        </div>
      </li>
    </body></html>""".encode()


def _write_fixture(
    raw: Path,
    species: str,
    category: str = "facultative allogamy",
    *,
    status: int = 200,
    batch: str = "pladias_test",
) -> dict:
    _write_raw_receipt(
        raw,
        species=species,
        retrieved_at="2026-07-19T01:00:00Z",
        batch_id=batch,
        http_status=status,
        response_url=_request_url(species),
        raw_content=_page(species, category) if status == 200 else b"not found",
    )
    return _load_raw_receipt(raw, species)


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


@pytest.mark.parametrize(
    ("category", "expected"),
    [
        ("allogamy", {("mating_system", "predominantly_outcrossing")}),
        ("facultative allogamy", {("mating_system", "predominantly_outcrossing")}),
        ("autogamy", {("mating_system", "obligate_selfing")}),
        ("facultative autogamy", {("mating_system", "predominantly_selfing")}),
        ("mixed mating", {("mating_system", "mixed_mating")}),
        (
            "allogamy self-incompatibility",
            {
                ("mating_system", "predominantly_outcrossing"),
                ("self_incompatibility", "SI"),
            },
        ),
        (
            "allogamy self-incompatibility, facultative allogamy",
            {
                ("mating_system", "predominantly_outcrossing"),
                ("self_incompatibility", "SI"),
            },
        ),
        (
            "allogamy self-incompatibility, mixed mating",
            {
                ("mating_system", "mixed_mating"),
                ("self_incompatibility", "SI"),
            },
        ),
        ("autogamy, facultative allogamy", set()),
        ("apomixis", set()),
    ],
)
def test_only_exact_pladias_categories_map(category: str, expected: set[tuple[str, str]]):
    observed = {(trait, value) for trait, value, _basis in _category_mappings(category)}
    assert observed == expected
    if "autogamy" in category:
        assert not any(trait == "autonomous_selfing_capacity" for trait, _value in observed)


def test_self_compatibility_is_never_inferred_from_autogamy():
    autogamy = _category_mappings("autogamy")
    assert not any(trait == "self_incompatibility" for trait, _value, _basis in autogamy)
    explicit = _category_mappings("autogamy; self-compatible")
    assert ("self_incompatibility", "SC") in {(trait, value) for trait, value, _basis in explicit}


def test_parser_takes_only_taxon_summary_field_and_requires_exact_heading():
    assert _extract_exact_field(_page(category="mixed mating"), "Alpha beta") == "mixed mating"
    wrong = _page(species="Gamma delta", category="allogamy")
    with pytest.raises(pladias_module.PladiasTransientError, match="heading_missing"):
        _extract_exact_field(wrong, "Alpha beta")


def test_receipt_yields_hashed_species_direct_unreviewed_leads(tmp_path: Path):
    species = "Allium sphaerocephalon"
    raw = tmp_path / "raw.response.html"
    receipt = _write_fixture(raw, species, "mixed mating")
    sources, leads = pladias_sources_from_receipt(receipt, raw)
    assert len(sources) == 1
    assert sources.loc[0, "accepted_species"] == species
    assert sources.loc[0, "pladias_category"] == "mixed mating"
    assert sources.loc[0, "pladias_year"] == PUBLICATION_YEAR
    assert sources.loc[0, "evidence_scope"] == "species_direct"
    assert sources.loc[0, "local_file_hash"] == _file_sha256(raw)
    assert sources.loc[0, "retrieval_date"] == "2026-07-19"
    assert sources.loc[0, "source_url"] == (
        "https://pladias.cz/en/taxon/data/Allium%20sphaerocephalon"
    )
    assert set(leads[["trait", "provisional_value"]].itertuples(index=False, name=None)) == {
        ("mating_system", "mixed_mating")
    }
    assert set(leads["review_status"]) == {"unreviewed"}


def test_allogamy_si_field_emits_two_leads_but_no_unreviewed_materialization(tmp_path: Path):
    species = "Alpha beta"
    raw = tmp_path / "raw.response.html"
    receipt = _write_fixture(raw, species, "allogamy self-incompatibility")
    sources, leads = pladias_sources_from_receipt(receipt, raw)
    assert len(sources) == 1
    assert set(leads[["trait", "provisional_value"]].itertuples(index=False, name=None)) == {
        ("mating_system", "predominantly_outcrossing"),
        ("self_incompatibility", "SI"),
    }
    assert "biological_values_materialized" not in sources.columns


def test_checkpoint_is_idempotent_recovers_running_and_rejects_non_binomial():
    base = _base_checkpoint("Alpha beta", "Gamma delta", "Agave", "Alpha beta subsp. minor")
    extended = ensure_pladias_checkpoint_rows(base, _queue("Alpha beta"))
    rows = extended.loc[extended["provider"].eq(PROVIDER)]
    assert len(rows) == 4 * len(RA_TRAITS)
    assert set(rows.loc[rows["species"].eq("Alpha beta"), "status"]) == {"pending"}
    assert set(rows.loc[rows["species"].eq("Gamma delta"), "status"]) == {"provider_seed_covered"}
    assert set(rows.loc[rows["species"].eq("Agave"), "status"]) == {"invalid_species_name"}
    assert set(rows.loc[rows["species"].eq("Alpha beta subsp. minor"), "status"]) == {
        "invalid_species_name"
    }
    mask = extended["provider"].eq(PROVIDER) & extended["species"].eq("Alpha beta")
    extended.loc[mask, "status"] = "running"
    recovered = ensure_pladias_checkpoint_rows(extended, _queue("Alpha beta"))
    assert set(recovered.loc[mask, "status"]) == {"retry"}
    assert len(recovered) == len(extended)


def test_receipt_rejects_wrong_url_hash_and_provider_version(tmp_path: Path):
    raw = tmp_path / "raw.response.html"
    _write_fixture(raw, "Alpha beta")
    receipt_path = raw.with_name("raw.receipt.json")
    metadata = json.loads(receipt_path.read_text())
    metadata["response_url"] = "https://pladias.cz/en/taxon/data/Alpha"
    receipt_path.write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="exact taxon URL"):
        _load_raw_receipt(raw, "Alpha beta")
    _write_fixture(raw, "Alpha beta")
    metadata = json.loads(receipt_path.read_text())
    metadata["provider_version"] = "old"
    receipt_path.write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="provider_version mismatch"):
        _load_raw_receipt(raw, "Alpha beta")
    assert PROVIDER_VERSION.endswith("_v1")
    _write_fixture(raw, "Alpha beta")
    metadata = json.loads(receipt_path.read_text())
    metadata["retrieved_at_utc"] = "not-a-date"
    receipt_path.write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="retrieval timestamp is invalid"):
        _load_raw_receipt(raw, "Alpha beta")
    _write_fixture(raw, "Alpha beta")
    raw.write_bytes(raw.read_bytes() + b" ")
    with pytest.raises(ValueError, match="response hash mismatch"):
        _load_raw_receipt(raw, "Alpha beta")


def test_collection_uses_only_encoded_exact_page_is_resumable_and_review_only(tmp_path: Path):
    requests = 0

    def handler(request: httpx.Request) -> httpx.Response:
        nonlocal requests
        requests += 1
        assert unquote(request.url.path) == "/en/taxon/data/Alpha beta"
        assert request.url.host == "pladias.cz"
        assert not request.url.query
        return httpx.Response(200, content=_page(), request=request)

    client = httpx.Client(transport=httpx.MockTransport(handler), follow_redirects=True)
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta", "Gamma delta")
    base.to_csv(checkpoint_path, index=False)
    report = collect_pladias_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "pladias",
        max_species=2,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    assert report["requested_species"] == 1
    assert report["exact_page_contract"] == EXACT_PAGE_CONTRACT
    assert report["candidate_leads"] == 1
    assert report["biological_values_materialized"] is False
    assert report["global_fallback_used"] is False
    assert report["genus_or_family_inference_used"] is False
    result = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    alpha = result.loc[result["provider"].eq(PROVIDER) & result["species"].eq("Alpha beta")]
    assert alpha["terminal"].str.casefold().eq("true").all()
    assert set(alpha.loc[alpha["trait"].eq("mating_system"), "status"]) == {
        "provider_pass_complete_with_candidates"
    }
    assert set(alpha.loc[~alpha["trait"].eq("mating_system"), "status"]) == {
        "provider_pass_complete"
    }
    second = collect_pladias_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=result,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "pladias",
        max_species=2,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    assert second["requested_species"] == 0
    assert requests == 1
    client.close()


def test_404_is_cached_terminal_no_result(tmp_path: Path):
    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(404, content=b"missing exact taxon", request=request)

    client = httpx.Client(transport=httpx.MockTransport(handler))
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    report = collect_pladias_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "pladias",
        max_species=1,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    result = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    rows = result.loc[result["provider"].eq(PROVIDER)]
    assert report["source_records"] == 0
    assert set(rows["status"]) == {"provider_no_result"}
    assert rows["terminal"].str.casefold().eq("true").all()
    raw = next((tmp_path / "pladias" / "raw").glob("*.response.html"))
    assert _load_raw_receipt(raw, "Alpha beta")["http_status"] == 404
    client.close()


def test_generic_taxon_index_redirect_is_cached_terminal_no_result(tmp_path: Path):
    requests = 0

    def handler(request: httpx.Request) -> httpx.Response:
        nonlocal requests
        requests += 1
        assert request.method == "GET"
        assert request.url.path.startswith("/en/taxon/data/")
        return httpx.Response(
            302,
            headers={"Location": "http://pladias.cz/en/taxon/?_fid=365k"},
            content=b"first-hop miss redirect",
            request=request,
        )

    client = httpx.Client(transport=httpx.MockTransport(handler), follow_redirects=True)
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta", "Gamma delta")
    base.to_csv(checkpoint_path, index=False)
    report = collect_pladias_batch(
        queue=_queue("Alpha beta", "Gamma delta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "pladias",
        max_species=2,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    result = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    rows = result.loc[result["provider"].eq(PROVIDER)]
    assert report["requested_species"] == 2
    assert report["source_records"] == 0
    assert report["transient_errors"] == 0
    assert set(rows["status"]) == {"provider_no_result"}
    assert rows["terminal"].str.casefold().eq("true").all()
    # The injected client follows redirects by default, but the collector's
    # per-request override must still stop after the exact-page response.
    assert requests == 2
    for raw in (tmp_path / "pladias" / "raw").glob("*.response.html"):
        receipt_path = raw.with_name(raw.name.removesuffix(".response.html") + ".receipt.json")
        species = json.loads(receipt_path.read_text())["species"]
        receipt = _load_raw_receipt(raw, species)
        assert receipt["http_status"] == 302
        assert receipt["response_url"] == _request_url(species)
        assert receipt["response_location"] == ("http://pladias.cz/en/taxon/?_fid=365k")
        assert receipt["payload"] == b"first-hop miss redirect"

    # Interrupted/nonterminal state must reuse the durable first-hop receipt,
    # not follow the redirect or recollect the exact endpoint.
    retry = result.copy()
    alpha = retry["provider"].eq(PROVIDER) & retry["species"].eq("Alpha beta")
    retry.loc[alpha, "status"] = "retry"
    retry.loc[alpha, "terminal"] = "false"
    retry.to_csv(checkpoint_path, index=False)
    resumed = collect_pladias_batch(
        queue=_queue("Alpha beta", "Gamma delta"),
        provider_checkpoint=retry,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "pladias",
        max_species=2,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    assert resumed["network_queries"] == 0
    assert resumed["cached_raw_responses_reused"] == 1
    assert requests == 2
    resumed_checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    resumed_alpha = resumed_checkpoint.loc[alpha]
    assert set(resumed_alpha["status"]) == {"provider_no_result"}
    assert resumed_alpha["terminal"].str.casefold().eq("true").all()
    client.close()


def test_legacy_followed_generic_index_receipt_remains_reusable(tmp_path: Path):
    raw = tmp_path / "raw.response.html"
    species = "Alpha beta"
    _write_raw_receipt(
        raw,
        species=species,
        retrieved_at="2026-07-19T01:00:00Z",
        batch_id="pladias_legacy",
        http_status=200,
        response_url="https://pladias.cz/en/taxon/?_fid=365k",
        raw_content=b"legacy generic taxon index",
    )
    receipt = _load_raw_receipt(raw, species)
    sources, leads = pladias_sources_from_receipt(receipt, raw)
    assert sources.empty
    assert leads.empty


@pytest.mark.parametrize(
    "location",
    [
        "https://example.org/en/taxon/?_fid=365k",
        "https://pladias.cz/en/taxon/data/Alpha%20beta",
        "https://pladias.cz/en/taxon/?_fid=365k&other=1",
        "https://user@pladias.cz/en/taxon/?_fid=365k",
        "",
    ],
)
def test_unrecognized_redirect_fails_closed_without_caching(
    tmp_path: Path,
    location: str,
):
    requests = 0

    def handler(request: httpx.Request) -> httpx.Response:
        nonlocal requests
        requests += 1
        headers = {"Location": location} if location else {}
        return httpx.Response(302, headers=headers, request=request)

    client = httpx.Client(transport=httpx.MockTransport(handler), follow_redirects=True)
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta", "Gamma delta")
    base.to_csv(checkpoint_path, index=False)
    report = collect_pladias_batch(
        queue=_queue("Alpha beta", "Gamma delta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "pladias",
        max_species=2,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    result = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    pladias = result.loc[result["provider"].eq(PROVIDER)]
    assert requests == 1
    assert report["transient_errors"] == 1
    assert set(pladias.loc[pladias["species"].eq("Alpha beta"), "status"]) == {"retry"}
    assert set(pladias.loc[pladias["species"].eq("Gamma delta"), "status"]) == {"pending"}
    assert not list((tmp_path / "pladias" / "raw").glob("*.response.html"))
    client.close()


@pytest.mark.parametrize("status", [403, 429, 500])
def test_provider_failures_remain_retryable_and_stop_when_provider_wide(
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
    report = collect_pladias_batch(
        queue=_queue("Alpha beta", "Gamma delta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "pladias",
        max_species=2,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    result = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    rows = result.loc[result["provider"].eq(PROVIDER)]
    assert report["transient_errors"] == (2 if status == 500 else 1)
    assert set(rows.loc[rows["species"].eq("Alpha beta"), "status"]) == {"retry"}
    if status in {403, 429}:
        assert requests == 1
        assert set(rows.loc[rows["species"].eq("Gamma delta"), "status"]) == {"pending"}
    else:
        assert requests == 2
        assert set(rows.loc[rows["species"].eq("Gamma delta"), "status"]) == {"retry"}
    client.close()


def test_shared_materialization_accepts_only_strict_reviewed_pladias_field(tmp_path: Path):
    species = "Alpha beta"
    raw = tmp_path / "raw.response.html"
    receipt = _write_fixture(raw, species, "facultative autogamy")
    source, _leads = pladias_sources_from_receipt(receipt, raw)
    provider = ensure_pladias_checkpoint_rows(_base_checkpoint(species), _queue(species))
    mask = provider["provider"].eq(PROVIDER)
    provider.loc[mask, "status"] = "provider_pass_complete_with_candidates"
    provider.loc[mask, "terminal"] = True
    provider.loc[mask, "attempts"] = 1
    provider.loc[mask, "updated_at"] = "2026-07-19T01:00:00Z"
    provider.loc[mask, "source_run_id"] = "pladias_test"
    source_path = tmp_path / "source_evidence.csv"
    source.to_csv(source_path, index=False)
    excerpt = source.loc[0, "source_text"]
    selection = {
        "accepted": [
            {
                "source_row_index": 0,
                "accepted_species": species,
                "source_url": source.loc[0, "source_url"],
                "accepted_trait": "mating_system",
                "accepted_value": "predominantly_selfing",
                "evidence_excerpt": excerpt,
                "review_reason": "Exact Pladias taxon field states facultative autogamy.",
            }
        ]
    }
    candidates, review = build_reviewed_search_batch(
        queue=_queue(species),
        search_evidence=source,
        provider_checkpoint=provider,
        selection_payload=selection,
        local_file_path=str(source_path),
        local_file_hash=_file_sha256(source_path),
    )
    assert candidates.loc[0, "provider"] == PROVIDER
    result = materialize_reviewed_cache(
        queue=_queue(species),
        baseline_checkpoint=_axis_checkpoint(species),
        provider_checkpoint=provider,
        candidates=candidates,
        review_payload=review,
    )
    record = result["records"].iloc[0]
    assert record["provider"] == PROVIDER
    assert record["trait"] == "mating_system"
    assert record["value"] == "predominantly_selfing"
    assert record["year"] == "2018"
    assert record["local_file_hash"] == _file_sha256(raw)

    disallowed = {
        "accepted": [
            {
                **selection["accepted"][0],
                "accepted_trait": "autonomous_selfing_capacity",
                "accepted_value": "autonomous",
                "review_reason": "Must fail: Pladias autogamy maps only to mating system.",
            }
        ]
    }
    with pytest.raises(ValueError, match="not an allowed exact-field mapping"):
        build_reviewed_search_batch(
            queue=_queue(species),
            search_evidence=source,
            provider_checkpoint=provider,
            selection_payload=disallowed,
            local_file_path=str(source_path),
            local_file_hash=_file_sha256(source_path),
        )


@pytest.mark.parametrize(
    ("value", "text", "expected"),
    [
        (
            "predominantly_outcrossing",
            "Alpha beta. Generative reproduction type: allogamy.",
            True,
        ),
        (
            "predominantly_outcrossing",
            "Alpha beta. Generative reproduction type: facultative allogamy.",
            True,
        ),
        (
            "predominantly_selfing",
            "Alpha beta. Generative reproduction type: facultative autogamy.",
            True,
        ),
        (
            "obligate_selfing",
            "Alpha beta. Generative reproduction type: autogamy.",
            True,
        ),
        (
            "obligate_selfing",
            "Alpha beta. Generative reproduction type: facultative autogamy.",
            False,
        ),
    ],
)
def test_shared_explicit_mapping_is_narrow_to_exact_pladias_field(
    value: str,
    text: str,
    expected: bool,
):
    assert _explicit_mapping_supported("mating_system", value, text) is expected


@pytest.mark.parametrize(
    ("text", "value", "expected"),
    [
        (
            "Alpha beta. Generative reproduction type: allogamy "
            "self-incompatibility, facultative allogamy.",
            "predominantly_outcrossing",
            True,
        ),
        (
            "Alpha beta. Generative reproduction type: allogamy "
            "self-incompatibility, mixed mating.",
            "mixed_mating",
            True,
        ),
        (
            "Alpha beta. Generative reproduction type: autogamy, facultative allogamy.",
            "obligate_selfing",
            False,
        ),
    ],
)
def test_shared_pladias_compound_category_mapping(text: str, value: str, expected: bool):
    assert _pladias_field_mapping_supported("mating_system", value, text) is expected


@pytest.mark.parametrize("spelling", ["selfpollinated", "self-pollinated", "self pollinated"])
def test_preanthesis_self_pollination_wording_is_explicit_autonomous_selfing(spelling: str):
    text = f"Alpha beta flowers could be {spelling} before anthesis."
    assert _explicit_mapping_supported("autonomous_selfing_capacity", "autonomous", text)
