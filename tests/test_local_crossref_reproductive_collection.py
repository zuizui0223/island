from pathlib import Path

import httpx
import pandas as pd
import pytest

from island_v2.local_crossref_reproductive_collection import (
    PROVIDER,
    _load_raw_receipt,
    _raw_cache_path,
    _write_raw_receipt,
    collect_crossref_batch,
    crossref_sources_from_receipt,
    ensure_crossref_checkpoint_rows,
)
from island_v2.local_reproductive_assurance_collection import CHECKPOINT_COLUMNS, RA_TRAITS
from island_v2.local_reproductive_assurance_collection import (
    _file_sha256,
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


def _payload(species: str = "Alpha beta") -> dict:
    return {
        "status": "ok",
        "message": {
            "items": [
                {
                    "title": [f"Breeding system of {species}"],
                    "abstract": (
                        f"<jats:p>{species} is self-incompatible and predominantly "
                        "outcrossing.</jats:p>"
                    ),
                    "published": {"date-parts": [[2024, 4, 2]]},
                    "DOI": "10.1000/example",
                    "container-title": ["Example Journal"],
                    "URL": "https://api.crossref.org/works/10.1000/example",
                }
            ]
        },
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


def test_checkpoint_extension_adds_four_crossref_keys_once():
    base = _base_checkpoint("Alpha beta", "Gamma delta", "Agave")
    extended = ensure_crossref_checkpoint_rows(base, _queue("Alpha beta"))
    crossref = extended.loc[extended["provider"].eq(PROVIDER)]
    assert len(crossref) == 3 * len(RA_TRAITS)
    assert not extended.duplicated(["species", "trait", "provider"]).any()
    assert set(crossref.loc[crossref["species"].eq("Alpha beta"), "status"]) == {"pending"}
    assert set(crossref.loc[crossref["species"].eq("Gamma delta"), "status"]) == {
        "provider_seed_covered"
    }
    assert set(crossref.loc[crossref["species"].eq("Agave"), "status"]) == {
        "invalid_species_name"
    }

    repeated = ensure_crossref_checkpoint_rows(extended, _queue("Alpha beta"))
    assert len(repeated) == len(extended)


def test_crossref_receipt_yields_dated_exact_species_review_leads(tmp_path: Path):
    raw = tmp_path / "raw.json"
    _write_raw_receipt(
        raw,
        species="Alpha beta",
        retrieved_at="2026-07-18T12:00:00Z",
        batch_id="crossref_test",
        payload=_payload(),
        rows=5,
    )
    receipt = _load_raw_receipt(raw, "Alpha beta", 5)
    sources, leads = crossref_sources_from_receipt(receipt, raw)
    assert len(sources) == 1
    assert sources.loc[0, "accepted_species"] == "Alpha beta"
    assert sources.loc[0, "evidence_scope"] == "species_direct"
    assert sources.loc[0, "crossref_year"] == "2024"
    assert "DOI:10.1000/example" in sources.loc[0, "source_citation"]
    assert len(sources.loc[0, "local_file_hash"]) == 64
    assert sources.loc[0, "retrieval_date"] == "2026-07-18"
    observed = set(leads[["trait", "provisional_value"]].itertuples(index=False, name=None))
    assert ("self_incompatibility", "SI") in observed
    assert ("mating_system", "predominantly_outcrossing") in observed
    assert set(leads["review_status"]) == {"unreviewed"}


def test_crossref_receipt_rejects_non_target_species_text(tmp_path: Path):
    raw = tmp_path / "raw.json"
    payload = _payload("Gamma delta")
    _write_raw_receipt(
        raw,
        species="Alpha beta",
        retrieved_at="2026-07-18T12:00:00Z",
        batch_id="crossref_test",
        payload=payload,
        rows=5,
    )
    receipt = _load_raw_receipt(raw, "Alpha beta", 5)
    sources, leads = crossref_sources_from_receipt(receipt, raw)
    assert sources.empty
    assert leads.empty


def test_created_only_date_and_species_suffix_are_not_direct_sources(tmp_path: Path):
    payload = _payload()
    created_only = payload["message"]["items"][0]
    created_only.pop("published")
    created_only["created"] = {"date-parts": [[2024, 4, 2]]}
    suffix = _payload()["message"]["items"][0]
    suffix["title"] = ["Breeding system of Alpha beta2"]
    suffix["abstract"] = "Alpha beta2 is self-incompatible."
    suffix["DOI"] = "10.1000/suffix"
    payload["message"]["items"].append(suffix)
    raw = tmp_path / "raw.response.json"
    _write_raw_receipt(
        raw,
        species="Alpha beta",
        retrieved_at="2026-07-18T12:00:00Z",
        batch_id="crossref_test",
        payload=payload,
        rows=5,
    )
    receipt = _load_raw_receipt(raw, "Alpha beta", 5)
    sources, leads = crossref_sources_from_receipt(receipt, raw)
    assert sources.empty
    assert leads.empty


def test_raw_receipt_validates_request_fingerprint_and_body_hash(tmp_path: Path):
    raw = tmp_path / "raw.response.json"
    _write_raw_receipt(
        raw,
        species="Alpha beta",
        retrieved_at="2026-07-18T12:00:00Z",
        batch_id="crossref_test",
        payload=_payload(),
        rows=5,
    )
    with pytest.raises(ValueError, match="row limit mismatch"):
        _load_raw_receipt(raw, "Alpha beta", 10)
    raw.write_bytes(raw.read_bytes() + b" ")
    with pytest.raises(ValueError, match="response hash mismatch"):
        _load_raw_receipt(raw, "Alpha beta", 5)


def test_collection_is_resumable_and_never_requeries_terminal_species(tmp_path: Path):
    requests = 0

    def handler(request: httpx.Request) -> httpx.Response:
        nonlocal requests
        requests += 1
        return httpx.Response(200, json=_payload(), request=request)

    client = httpx.Client(transport=httpx.MockTransport(handler))
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    output_dir = tmp_path / "crossref"
    report = collect_crossref_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=output_dir,
        max_species=1,
        rows=5,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    assert report["requested_species"] == 1
    assert report["candidate_species"] == 1
    assert requests == 1
    checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    crossref = checkpoint.loc[checkpoint["provider"].eq(PROVIDER)]
    assert crossref["terminal"].str.casefold().eq("true").all()
    assert set(crossref["attempts"]) == {"1"}
    assert (output_dir / "source_evidence.csv").exists()
    assert len(list((output_dir / "raw").glob("*.response.json"))) == 1
    assert len(list((output_dir / "raw").glob("*.receipt.json"))) == 1

    second = collect_crossref_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=checkpoint,
        provider_checkpoint_path=checkpoint_path,
        output_dir=output_dir,
        max_species=1,
        rows=5,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    assert second["requested_species"] == 0
    assert requests == 1
    client.close()


def test_crossref_review_preserves_raw_response_provenance_and_exact_excerpt(
    tmp_path: Path,
):
    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json=_payload(), request=request)

    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    output_dir = tmp_path / "crossref"
    client = httpx.Client(transport=httpx.MockTransport(handler))
    report = collect_crossref_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=output_dir,
        max_species=1,
        rows=5,
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
    assert candidates.loc[0, "local_file_path"].endswith(".response.json")
    assert candidates.loc[0, "local_file_hash"] == source.loc[0, "local_file_hash"]
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
    assert record["local_file_path"].endswith(".response.json")
    assert record["local_file_hash"] == _file_sha256(Path(record["local_file_path"]))
    crossref = result["provider_checkpoint"].loc[
        result["provider_checkpoint"]["provider"].eq(PROVIDER)
    ]
    accepted = crossref.loc[crossref["trait"].eq("self_incompatibility")].iloc[0]
    assert accepted["status"] == "accepted_evidence"
    assert set(crossref.loc[~crossref["trait"].eq("self_incompatibility"), "terminal"]) == {
        True
    }

    ambiguous_selection = {
        "accepted": [
            {
                **selection["accepted"][0],
                "evidence_excerpt": "self-incompatible",
            }
        ]
    }
    ambiguous_candidates, ambiguous_review = build_reviewed_search_batch(
        queue=_queue("Alpha beta"),
        search_evidence=source,
        provider_checkpoint=providers,
        selection_payload=ambiguous_selection,
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


def test_transient_failure_remains_retryable(tmp_path: Path):
    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(429, request=request)

    client = httpx.Client(transport=httpx.MockTransport(handler))
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    base = _base_checkpoint("Alpha beta")
    base.to_csv(checkpoint_path, index=False)
    report = collect_crossref_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=base,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "crossref",
        max_species=1,
        rows=5,
        min_interval_seconds=0,
        max_retries=0,
        client=client,
        sleeper=lambda _seconds: None,
    )
    checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    crossref = checkpoint.loc[checkpoint["provider"].eq(PROVIDER)]
    assert set(crossref["status"]) == {"retry"}
    assert crossref["terminal"].str.casefold().eq("false").all()
    assert set(crossref["attempts"]) == {"1"}
    assert report["transient_errors"] == 1
    client.close()


def test_partial_terminal_key_is_never_rewritten(tmp_path: Path):
    requests = 0

    def handler(request: httpx.Request) -> httpx.Response:
        nonlocal requests
        requests += 1
        return httpx.Response(200, json=_payload(), request=request)

    checkpoint = ensure_crossref_checkpoint_rows(
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
    client = httpx.Client(transport=httpx.MockTransport(handler))
    collect_crossref_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=checkpoint,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "crossref",
        max_species=1,
        rows=5,
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
    assert requests == 1
    client.close()


def test_interrupted_raw_receipt_is_reused_without_new_request_or_attempt(tmp_path: Path):
    checkpoint = ensure_crossref_checkpoint_rows(
        _base_checkpoint("Alpha beta"), _queue("Alpha beta")
    )
    crossref = checkpoint["provider"].eq(PROVIDER) & checkpoint["species"].eq("Alpha beta")
    checkpoint.loc[crossref, "status"] = "running"
    checkpoint.loc[crossref, "terminal"] = False
    checkpoint.loc[crossref, "attempts"] = 1
    checkpoint_path = tmp_path / "provider_checkpoint.csv"
    checkpoint.to_csv(checkpoint_path, index=False)
    raw_dir = tmp_path / "crossref" / "raw"
    raw_dir.mkdir(parents=True)
    _write_raw_receipt(
        _raw_cache_path(raw_dir, "Alpha beta", 5),
        species="Alpha beta",
        retrieved_at="2026-07-18T12:00:00Z",
        batch_id="crossref_interrupted",
        payload=_payload(),
        rows=5,
    )

    def handler(_request: httpx.Request) -> httpx.Response:
        raise AssertionError("raw receipt recovery must not make a network request")

    client = httpx.Client(transport=httpx.MockTransport(handler))
    report = collect_crossref_batch(
        queue=_queue("Alpha beta"),
        provider_checkpoint=checkpoint,
        provider_checkpoint_path=checkpoint_path,
        output_dir=tmp_path / "crossref",
        max_species=1,
        rows=5,
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
