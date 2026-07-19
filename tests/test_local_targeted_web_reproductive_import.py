from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pytest

import island_v2.local_targeted_web_reproductive_import as targeted
from island_v2.local_reproductive_assurance_collection import CHECKPOINT_COLUMNS, RA_TRAITS
from island_v2.local_targeted_web_reproductive_import import (
    APPROVAL_CONTRACT,
    INSPECTED_SPECIES,
    PROVIDER,
    build_targeted_web_review_batch,
    materialize_targeted_web_lane,
    run_targeted_web_import,
    run_targeted_web_round_import,
)


@dataclass
class Lane:
    root: Path
    queue: Path
    candidates: Path
    rejected: Path
    approval: Path
    baseline: Path
    providers: Path
    raw: Path


def _axis_rows(species_names: tuple[str, ...] = INSPECTED_SPECIES) -> pd.DataFrame:
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
                "provider_state_available": False,
                "checkpoint_scope": "materialized_axis_baseline_only",
            }
            for species in species_names
        ]
    )


def _make_lane(
    tmp_path: Path,
    *,
    species: str = "Penstemon digitalis",
    excerpt: str | None = None,
    accepted_trait: str = "self_incompatibility",
    accepted_value: str = "SC",
    approve: bool = True,
    inspected_species: tuple[str, ...] = INSPECTED_SPECIES,
    discovery_queue_hash: str = targeted.DISCOVERY_QUEUE_SHA256,
) -> Lane:
    root = tmp_path / "targeted_web_discovery"
    raw_dir = root / "raw"
    raw_dir.mkdir(parents=True)
    evidence = excerpt or f"{species} is self-compatible."
    raw = raw_dir / "source.txt"
    raw.write_text(evidence, encoding="utf-8")
    local_hash = targeted._file_sha256(raw)

    candidates = root / "candidate_records.csv"
    pd.DataFrame(
        [
            {
                "status": "candidate",
                "species": species,
                "trait": "reproductive_assurance",
                "category": "self-compatible",
                "raw_text": evidence,
                "url": "https://doi.org/10.1000/targeted-example",
                "retrieval_url": "https://example.org/source.txt",
                "doi": "10.1000/targeted-example",
                "publication": f"Breeding system of {species}",
                "venue": "Example Journal 1:1-2",
                "year": "2024",
                "source_type": "journal_article_primary_pollination_experiment",
                "evidence_basis": "Direct experiment.",
                "local_file": "raw/source.txt",
                "local_file_sha256": local_hash,
                "retrieval_date": "2026-07-19",
                "review_note": "Synthetic test fixture.",
            }
        ]
    ).to_csv(candidates, index=False)
    rejected = root / "rejected_or_held_hits.csv"
    pd.DataFrame(
        [
            {
                "status": "rejected_no_hit",
                "species": other,
                "possible_category": "",
                "source_or_hit": "",
                "doi": "",
                "reason": "No approved exact evidence.",
                "local_raw": "",
                "retrieval_date": "2026-07-19",
            }
            for other in inspected_species
            if other != species
        ]
    ).to_csv(rejected, index=False)

    queue = tmp_path / "queue.csv"
    pd.DataFrame({"accepted_species": inspected_species}).to_csv(queue, index=False)
    baseline = tmp_path / "three_axis_checkpoint.csv"
    _axis_rows(inspected_species).to_csv(baseline, index=False)
    providers = tmp_path / "provider_checkpoint.csv"
    pd.DataFrame(columns=CHECKPOINT_COLUMNS).to_csv(providers, index=False)

    accepted = []
    if approve:
        accepted.append(
            {
                "discovery_row_index": 0,
                "species": species,
                "source_url": "https://doi.org/10.1000/targeted-example",
                "local_file_sha256": local_hash,
                "accepted_trait": accepted_trait,
                "accepted_value": accepted_value,
                "evidence_excerpt": evidence,
                "review_reason": "Exact full species and mapped condition are explicit.",
            }
        )
    approval = root / "approval.json"
    approval.write_text(
        json.dumps(
            {
                "contract_version": APPROVAL_CONTRACT,
                "approval_id": "targeted-web-test-approval",
                "discovery_queue_sha256": discovery_queue_hash,
                "current_queue_sha256": targeted._file_sha256(queue),
                "candidate_file_sha256": targeted._file_sha256(candidates),
                "rejected_file_sha256": targeted._file_sha256(rejected),
                "inspected_species": list(inspected_species),
                "completed_at_utc": "2026-07-19T02:00:00+09:00",
                "accepted": accepted,
            },
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    return Lane(root, queue, candidates, rejected, approval, baseline, providers, raw)


def _run(lane: Lane, output: Path, *, providers: Path | None = None) -> dict[str, object]:
    return run_targeted_web_import(
        queue_csv=lane.queue,
        candidate_csv=lane.candidates,
        rejected_csv=lane.rejected,
        approval_json=lane.approval,
        baseline_checkpoint_csv=lane.baseline,
        provider_checkpoint_csv=providers or lane.providers,
        output_dir=output,
        expected_current_queue_hash=targeted._file_sha256(lane.queue),
    )


def _run_round(
    lane: Lane,
    output: Path,
    *,
    providers: Path,
    prior_records: Path,
    discovery_queue_hash: str,
) -> dict[str, object]:
    return run_targeted_web_round_import(
        queue_csv=lane.queue,
        candidate_csv=lane.candidates,
        rejected_csv=lane.rejected,
        approval_json=lane.approval,
        baseline_checkpoint_csv=lane.baseline,
        provider_checkpoint_csv=providers,
        output_dir=output,
        expected_current_queue_hash=targeted._file_sha256(lane.queue),
        expected_discovery_queue_hash=discovery_queue_hash,
        prior_records_csv=prior_records,
    )


def _update_approval(lane: Lane, **updates: object) -> None:
    payload = json.loads(lane.approval.read_text(encoding="utf-8"))
    payload.update(updates)
    lane.approval.write_text(json.dumps(payload), encoding="utf-8")


def test_exact_excerpt_can_be_verified_from_hashed_pdf_text(monkeypatch):
    excerpt = (
        "Silene dioica populations were studied. The obligately outbreeding "
        "nature of S. dioica was established."
    )

    class FakePage:
        def extract_text(self) -> str:
            return excerpt

    class FakeReader:
        pages = [FakePage()]

    monkeypatch.setattr(targeted, "PdfReader", lambda _stream: FakeReader())

    assert targeted._excerpt_is_exact(excerpt, "", b"%PDF-1.7\nfixture")


def test_invalid_pdf_fails_closed_during_exact_excerpt_check(monkeypatch):
    def raise_on_read(_stream):
        raise ValueError("invalid PDF")

    monkeypatch.setattr(targeted, "PdfReader", raise_on_read)

    assert not targeted._excerpt_is_exact(
        "Silene dioica is obligately outbreeding.",
        "",
        b"%PDF-1.7\ninvalid",
    )


def test_input_hashes_are_snapshotted_before_in_place_outputs(tmp_path: Path):
    lane = _make_lane(tmp_path)
    queue_hash = targeted._file_sha256(lane.queue)
    baseline_hash = targeted._file_sha256(lane.baseline)
    provider_hash = targeted._file_sha256(lane.providers)

    report = _run(lane, tmp_path)

    assert report["inputs"][str(lane.queue.resolve())]["sha256"] == queue_hash
    assert report["inputs"][str(lane.baseline.resolve())]["sha256"] == baseline_hash
    assert report["inputs"][str(lane.providers.resolve())]["sha256"] == provider_hash


def test_local_raw_hash_tamper_is_rejected(tmp_path: Path):
    lane = _make_lane(tmp_path)
    lane.raw.write_text("Penstemon digitalis is self-incompatible.", encoding="utf-8")
    with pytest.raises(ValueError, match="local raw hash mismatch"):
        _run(lane, tmp_path / "out")


def test_queue_hash_and_current_membership_are_both_required(tmp_path: Path):
    lane = _make_lane(tmp_path)
    with pytest.raises(ValueError, match="current queue hash"):
        run_targeted_web_import(
            queue_csv=lane.queue,
            candidate_csv=lane.candidates,
            rejected_csv=lane.rejected,
            approval_json=lane.approval,
            baseline_checkpoint_csv=lane.baseline,
            provider_checkpoint_csv=lane.providers,
            output_dir=tmp_path / "bad-hash",
            expected_current_queue_hash="0" * 64,
        )

    queue = pd.read_csv(lane.queue, dtype=str)
    queue.loc[queue["accepted_species"].eq("Penstemon digitalis"), "accepted_species"] = (
        "Outside species"
    )
    queue.to_csv(lane.queue, index=False)
    new_hash = targeted._file_sha256(lane.queue)
    _update_approval(lane, current_queue_sha256=new_hash)
    with pytest.raises(ValueError, match="outside the current queue"):
        _run(lane, tmp_path / "outside")


def test_abbreviated_only_species_is_rejected(tmp_path: Path):
    lane = _make_lane(
        tmp_path,
        species="Ipomoea hederifolia",
        excerpt="I. hederifolia is self-compatible.",
    )
    with pytest.raises(ValueError, match="exact full accepted species"):
        _run(lane, tmp_path / "out")


def test_held_or_rejected_species_cannot_enter_candidate_lane(tmp_path: Path):
    lane = _make_lane(tmp_path)
    rejected = pd.read_csv(lane.rejected, dtype=str).fillna("")
    rejected.loc[len(rejected)] = {
        **rejected.iloc[0].to_dict(),
        "species": "Penstemon digitalis",
        "reason": "Held for review.",
        "status": "held_citation_dependent",
    }
    rejected.to_csv(lane.rejected, index=False)
    _update_approval(lane, rejected_file_sha256=targeted._file_sha256(lane.rejected))
    with pytest.raises(ValueError, match="held/rejected species"):
        _run(lane, tmp_path / "out")


def test_approval_excerpt_must_be_exact_source_text(tmp_path: Path):
    lane = _make_lane(tmp_path)
    payload = json.loads(lane.approval.read_text(encoding="utf-8"))
    payload["accepted"][0]["evidence_excerpt"] = (
        "Penstemon digitalis has a self-compatible breeding system."
    )
    lane.approval.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="not exact source text"):
        _run(lane, tmp_path / "out")


def test_final_materialization_revalidates_raw_hash(tmp_path: Path):
    lane = _make_lane(tmp_path)
    queue_hash = targeted._file_sha256(lane.queue)
    source, candidates, review, _metadata = build_targeted_web_review_batch(
        queue=pd.read_csv(lane.queue, dtype=str).fillna(""),
        candidate_records=pd.read_csv(lane.candidates, dtype=str).fillna(""),
        rejected_hits=pd.read_csv(lane.rejected, dtype=str).fillna(""),
        approval_payload=json.loads(lane.approval.read_text(encoding="utf-8")),
        current_queue_hash=queue_hash,
        candidate_hash=targeted._file_sha256(lane.candidates),
        rejected_hash=targeted._file_sha256(lane.rejected),
        candidate_csv=lane.candidates,
        expected_current_queue_hash=queue_hash,
    )
    lane.raw.write_text("tampered after review batch construction", encoding="utf-8")
    with pytest.raises(ValueError, match="final raw-response provenance mismatch"):
        materialize_targeted_web_lane(
            queue=pd.read_csv(lane.queue, dtype=str).fillna(""),
            baseline_checkpoint=pd.read_csv(lane.baseline, dtype=str).fillna(""),
            provider_checkpoint=pd.read_csv(lane.providers, dtype=str).fillna(""),
            source_evidence=source,
            candidates=candidates,
            review_payload=review,
        )


def test_terminal_rows_are_idempotent_and_immutable(tmp_path: Path):
    lane = _make_lane(tmp_path)
    first = tmp_path / "first"
    _run(lane, first)
    second = tmp_path / "second"
    _run(lane, second, providers=first / "provider_checkpoint.csv")
    first_checkpoint = pd.read_csv(first / "provider_checkpoint.csv", dtype=str).fillna("")
    second_checkpoint = pd.read_csv(second / "provider_checkpoint.csv", dtype=str).fillna("")
    pd.testing.assert_frame_equal(first_checkpoint, second_checkpoint)

    _update_approval(lane, accepted=[])
    with pytest.raises(ValueError, match="terminal targeted-web checkpoint rows are immutable"):
        _run(lane, tmp_path / "changed", providers=first / "provider_checkpoint.csv")


def test_dynamic_round_appends_disjoint_blocks_and_preserves_existing_rows(
    tmp_path: Path,
):
    initial_lane = _make_lane(tmp_path / "initial")
    initial_output = tmp_path / "initial-output"
    _run(initial_lane, initial_output)

    round_species = ("Achillea millefolium", "Acorus calamus")
    discovery_hash = "a" * 64
    round_lane = _make_lane(
        tmp_path / "round",
        species="Achillea millefolium",
        inspected_species=round_species,
        discovery_queue_hash=discovery_hash,
    )
    cumulative_axes = pd.concat(
        [
            pd.read_csv(initial_output / "three_axis_checkpoint.csv.gz", dtype=str).fillna(""),
            _axis_rows(round_species),
        ],
        ignore_index=True,
    )
    cumulative_axes.to_csv(round_lane.baseline, index=False)

    round_output = tmp_path / "round-output"
    report = _run_round(
        round_lane,
        round_output,
        providers=initial_output / "provider_checkpoint.csv",
        prior_records=initial_output / "new_reproductive_assurance_records.csv",
        discovery_queue_hash=discovery_hash,
    )

    before = pd.read_csv(initial_output / "provider_checkpoint.csv", dtype=str).fillna("")
    after = pd.read_csv(round_output / "provider_checkpoint.csv", dtype=str).fillna("")
    before_targeted = before.loc[before["provider"].eq(PROVIDER)]
    preserved = after.loc[
        after["provider"].eq(PROVIDER) & after["species"].isin(set(INSPECTED_SPECIES))
    ]
    pd.testing.assert_frame_equal(
        targeted._canonical_checkpoint(before_targeted),
        targeted._canonical_checkpoint(preserved),
    )
    round_rows = after.loc[after["provider"].eq(PROVIDER) & after["species"].isin(round_species)]
    assert len(round_rows) == len(round_species) * len(RA_TRAITS)
    assert set(round_rows.groupby("species")["trait"].apply(frozenset)) == {frozenset(RA_TRAITS)}
    assert round_rows["terminal"].str.casefold().eq("true").all()

    assert report["targeted_web_provider_rows"] == 22 * len(RA_TRAITS)
    assert report["targeted_web_round_provider_rows"] == 2 * len(RA_TRAITS)
    assert report["targeted_web_existing_provider_rows"] == 20 * len(RA_TRAITS)
    assert report["targeted_web_appended_provider_rows"] == 2 * len(RA_TRAITS)
    assert report["targeted_web_total_provider_rows"] == 22 * len(RA_TRAITS)
    assert report["targeted_web_inspected_species"] == 2
    assert report["targeted_web_round_inspected_species"] == 2
    assert report["targeted_web_appended_species"] == 2
    assert report["targeted_web_total_inspected_species"] == 22
    assert report["new_record_count"] == 1
    assert report["total_record_count"] == 2
    assert report["targeted_web_approved_sources"] == 1
    assert report["targeted_web_import_performed_external_acquisition"] is False


def test_dynamic_round_species_are_exactly_approval_and_discovery_union(tmp_path: Path):
    round_species = ("Achillea millefolium", "Acorus calamus")
    discovery_hash = "b" * 64
    lane = _make_lane(
        tmp_path,
        species="Achillea millefolium",
        inspected_species=round_species,
        discovery_queue_hash=discovery_hash,
    )
    _update_approval(lane, inspected_species=["Achillea millefolium"])

    with pytest.raises(ValueError, match="exact inspected species set"):
        run_targeted_web_round_import(
            queue_csv=lane.queue,
            candidate_csv=lane.candidates,
            rejected_csv=lane.rejected,
            approval_json=lane.approval,
            baseline_checkpoint_csv=lane.baseline,
            provider_checkpoint_csv=lane.providers,
            output_dir=tmp_path / "out",
            expected_current_queue_hash=targeted._file_sha256(lane.queue),
            expected_discovery_queue_hash=discovery_hash,
            prior_records_csv=lane.providers,
        )


def test_append_overlap_requires_identical_terminal_state(tmp_path: Path):
    lane = _make_lane(tmp_path)
    output = tmp_path / "out"
    _run(lane, output)
    existing = pd.read_csv(output / "provider_checkpoint.csv", dtype=str).fillna("")
    candidates = pd.read_csv(output / "reviewed_candidates.csv", dtype=str).fillna("")
    review = json.loads((output / "targeted_web_review.json").read_text(encoding="utf-8"))
    planned = targeted._targeted_checkpoint_rows(
        candidates,
        review,
        inspected_species=INSPECTED_SPECIES,
        already_covered_species=frozenset(),
    )

    seeded, preserved, appended = targeted._seed_targeted_checkpoint(
        existing,
        planned,
        allow_append=True,
    )
    assert appended == 0
    assert len(seeded) == len(existing)
    assert len(preserved) == len(INSPECTED_SPECIES) * len(RA_TRAITS)

    changed = planned.copy()
    changed.loc[0, "source_run_id"] = "different-round"
    with pytest.raises(ValueError, match="overlapping targeted-web checkpoint rows"):
        targeted._seed_targeted_checkpoint(existing, changed, allow_append=True)


def test_append_rejects_nonterminal_existing_targeted_rows(tmp_path: Path):
    lane = _make_lane(tmp_path)
    output = tmp_path / "out"
    _run(lane, output)
    existing = pd.read_csv(output / "provider_checkpoint.csv", dtype=str).fillna("")
    existing.loc[0, ["status", "terminal"]] = ["pending", "false"]
    candidates = pd.read_csv(output / "reviewed_candidates.csv", dtype=str).fillna("")
    review = json.loads((output / "targeted_web_review.json").read_text(encoding="utf-8"))
    planned = targeted._targeted_checkpoint_rows(
        candidates,
        review,
        inspected_species=INSPECTED_SPECIES,
        already_covered_species=frozenset(),
    )

    with pytest.raises(ValueError, match="contains nonterminal state"):
        targeted._seed_targeted_checkpoint(existing, planned, allow_append=True)


def test_discovery_candidates_are_not_auto_accepted(tmp_path: Path):
    lane = _make_lane(tmp_path, approve=False)
    output = tmp_path / "out"
    report = _run(lane, output)
    assert pd.read_csv(output / "source_evidence.csv").empty
    assert pd.read_csv(output / "reviewed_candidates.csv").empty
    assert pd.read_csv(output / "new_reproductive_assurance_records.csv").empty
    checkpoint = pd.read_csv(output / "provider_checkpoint.csv", dtype=str).fillna("")
    assert len(checkpoint) == len(INSPECTED_SPECIES) * len(RA_TRAITS)
    assert set(checkpoint["status"]) == {"provider_pass_complete"}
    assert report["external_acquisition_started"] is True
    assert report["targeted_web_import_performed_external_acquisition"] is False


def test_species_completed_after_discovery_is_terminal_seed_covered(tmp_path: Path):
    lane = _make_lane(tmp_path)
    queue = pd.read_csv(lane.queue, dtype=str)
    queue = queue.loc[~queue["accepted_species"].eq("Solanum melongena")]
    queue.to_csv(lane.queue, index=False)
    axes = pd.read_csv(lane.baseline, dtype=str).fillna("")
    mask = axes["accepted_species"].eq("Solanum melongena")
    axes.loc[mask, "reproductive_assurance_direct"] = "True"
    axes.loc[mask, "reproductive_assurance_evidence_traits"] = "self_incompatibility"
    axes.to_csv(lane.baseline, index=False)
    _update_approval(
        lane,
        current_queue_sha256=targeted._file_sha256(lane.queue),
    )

    output = tmp_path / "out"
    report = _run(lane, output)
    checkpoint = pd.read_csv(output / "provider_checkpoint.csv", dtype=str).fillna("")
    solanum = checkpoint.loc[
        checkpoint["provider"].eq(PROVIDER) & checkpoint["species"].eq("Solanum melongena")
    ]
    assert len(solanum) == len(RA_TRAITS)
    assert set(solanum["status"]) == {"provider_seed_covered"}
    assert report["targeted_web_already_covered_species"] == ["Solanum melongena"]


def test_end_to_end_materialization_writes_checkpoint_last(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    lane = _make_lane(tmp_path)
    order: list[str] = []
    original_csv = targeted._atomic_write_csv
    original_text = targeted._atomic_write_text

    def tracked_csv(frame: pd.DataFrame, path: Path, *, gzip: bool = False) -> None:
        order.append(path.name)
        original_csv(frame, path, gzip=gzip)

    def tracked_text(value: str, path: Path) -> None:
        order.append(path.name)
        original_text(value, path)

    monkeypatch.setattr(targeted, "_atomic_write_csv", tracked_csv)
    monkeypatch.setattr(targeted, "_atomic_write_text", tracked_text)
    output = tmp_path / "out"
    report = _run(lane, output)

    assert order[-1] == "provider_checkpoint.csv"
    assert report["checkpoint_written_last"] is True
    sources = pd.read_csv(output / "source_evidence.csv", dtype=str).fillna("")
    candidates = pd.read_csv(output / "reviewed_candidates.csv", dtype=str).fillna("")
    records = pd.read_csv(output / "new_reproductive_assurance_records.csv", dtype=str).fillna("")
    providers = pd.read_csv(output / "provider_checkpoint.csv", dtype=str).fillna("")
    axes = pd.read_csv(output / "three_axis_checkpoint.csv.gz", dtype=str).fillna("")
    assert len(sources) == len(candidates) == len(records) == 1
    assert records.loc[0, "species"] == "Penstemon digitalis"
    assert records.loc[0, "trait"] == "self_incompatibility"
    assert records.loc[0, "value"] == "SC"
    assert records.loc[0, "provider"] == PROVIDER
    assert records.loc[0, "DOI"] == "10.1000/targeted-example"
    targeted_rows = providers.loc[providers["provider"].eq(PROVIDER)]
    assert len(targeted_rows) == len(INSPECTED_SPECIES) * len(RA_TRAITS)
    assert targeted_rows["terminal"].str.casefold().eq("true").all()
    assert set(targeted_rows["species"]) == set(INSPECTED_SPECIES)
    accepted = targeted_rows.loc[
        targeted_rows["species"].eq("Penstemon digitalis")
        & targeted_rows["trait"].eq("self_incompatibility")
    ].iloc[0]
    assert accepted["status"] == "accepted_evidence"
    assert accepted["accepted_record_count"] == "1"
    penstemon = axes.loc[axes["accepted_species"].eq("Penstemon digitalis")].iloc[0]
    assert penstemon["reproductive_assurance_direct"].casefold() == "true"
    assert penstemon["all_three_axes_direct"].casefold() == "true"
    assert report["targeted_web_global_fallback_used"] is False
    assert report["targeted_web_discovery_auto_accepted"] is False
    assert report["external_acquisition_started"] is True
    assert report["targeted_web_import_performed_external_acquisition"] is False
