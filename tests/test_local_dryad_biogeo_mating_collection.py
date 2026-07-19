import json
from pathlib import Path

import pandas as pd
import pytest

import island_v2.local_reproductive_assurance_collection as reproductive_module
from island_v2.local_dryad_biogeo_mating_collection import (
    ARTICLE_DOI,
    CHECKPOINT_COLUMNS,
    DATASET_DOI,
    PROVIDER,
    SOURCE_COLUMNS,
    TRAIT,
    _file_sha256,
    collect_dryad_biogeo_candidates,
    ensure_dryad_checkpoint_rows,
)
from island_v2.local_reproductive_assurance_collection import (
    build_reviewed_search_batch,
    materialize_reviewed_cache,
)


def _queue(*rows: tuple[str, str]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["accepted_species", "family"])


def _source_row(
    genus: str,
    species: str,
    family: str,
    *,
    si: str,
    mean_tm: str = "0.5",
    year: str = "2001",
    reference: str = "Primary reference.",
    references: str = "Full original literature reference.",
) -> dict[str, str]:
    return {
        "year": year,
        "reference": reference,
        "family": family,
        "genus": genus,
        "species": species,
        "m.d.g": "D",
        "latitude": "12.3",
        "hemisphere": "N",
        "biome": "6",
        "growth": "1",
        "life.history": "4",
        "si": si,
        "mean.tm": mean_tm,
        "References": references,
    }


def _write_source(path: Path, rows: list[dict[str, str]]) -> str:
    pd.DataFrame(rows).to_csv(path, index=False, encoding="latin-1")
    return _file_sha256(path)


def _axis(queue: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "accepted_species": row.accepted_species,
                "genus": row.accepted_species.split()[0],
                "family": row.family,
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
            for row in queue.itertuples(index=False)
        ]
    )


def _run(
    tmp_path: Path,
    *,
    original: pd.DataFrame,
    current: pd.DataFrame,
    rows: list[dict[str, str]],
    checkpoint: pd.DataFrame | None = None,
) -> tuple[dict, Path, Path]:
    source = tmp_path / "Database_BioGeo_mating.csv"
    expected_hash = _write_source(source, rows)
    output = tmp_path / "provider"
    checkpoint_path = output / "provider_checkpoint.csv"
    report = collect_dryad_biogeo_candidates(
        original_queue=original,
        current_queue=current,
        provider_checkpoint=(
            checkpoint if checkpoint is not None else pd.DataFrame(columns=CHECKPOINT_COLUMNS)
        ),
        provider_checkpoint_path=checkpoint_path,
        source_csv=source,
        output_dir=output,
        retrieval_date="2026-07-19",
        expected_csv_sha256=expected_hash,
    )
    return report, output, checkpoint_path


def _review_fixture(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    si: str = "0",
    mean_tm: str = "0.84",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Path]:
    queue = _queue(("Alpha beta", "Fabaceae"))
    _report, output, checkpoint_path = _run(
        tmp_path,
        original=queue,
        current=queue,
        rows=[
            _source_row(
                "Alpha",
                "beta",
                "Fabaceae",
                si=si,
                mean_tm=mean_tm,
                year="1956",
            )
        ],
    )
    search = pd.read_csv(output / "source_evidence.csv", dtype=str).fillna("")
    providers = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    monkeypatch.setattr(
        reproductive_module,
        "DRYAD_BIOGEO_CSV_SHA256",
        search.loc[0, "local_file_hash"],
    )
    return queue, search, providers, output / "source_evidence.csv"


def _selection(search: pd.DataFrame, *, value: str = "SI") -> dict[str, object]:
    return {
        "accepted": [
            {
                "source_row_index": 0,
                "accepted_species": "Alpha beta",
                "source_url": search.loc[0, "source_url"],
                "accepted_trait": "self_incompatibility",
                "accepted_value": value,
                "evidence_excerpt": search.loc[0, "source_text"],
                "review_reason": "Exact pinned Dryad si category reviewed.",
            }
        ]
    }


def _build_reviewed(
    queue: pd.DataFrame,
    search: pd.DataFrame,
    providers: pd.DataFrame,
    search_path: Path,
    *,
    value: str = "SI",
) -> tuple[pd.DataFrame, dict[str, object]]:
    return build_reviewed_search_batch(
        queue=queue,
        search_evidence=search,
        provider_checkpoint=providers,
        selection_payload=_selection(search, value=value),
        local_file_path=str(search_path),
        local_file_hash=_file_sha256(search_path),
    )


def test_static_scan_uses_only_explicit_si_and_preserves_provenance(tmp_path: Path):
    original = _queue(
        ("Alpha beta", "Fabaceae"),
        ("Gamma delta", "Lamiaceae"),
        ("Covered species", "Rosaceae"),
    )
    current = original.iloc[:2].copy()
    rows = [
        _source_row(
            "Alpha",
            "beta",
            "Fabaceae",
            si="0",
            mean_tm="0.999",
            year="1956",
            reference="Short primary citation.",
            references="Complete primary citation.",
        ),
        _source_row("Gamma", "delta", "Lamiaceae", si="", mean_tm="0.000"),
    ]
    report, output, checkpoint_path = _run(
        tmp_path,
        original=original,
        current=current,
        rows=rows,
    )

    leads = pd.read_csv(output / "candidate_leads.csv", dtype=str).fillna("")
    assert len(leads) == 1
    lead = leads.iloc[0]
    assert lead["accepted_species"] == "Alpha beta"
    assert lead["trait_name"] == TRAIT
    assert lead["provisional_candidate_value"] == "SI"
    assert lead["review_status"] == "unreviewed"
    assert lead["dataset_doi"] == DATASET_DOI
    assert lead["article_doi"] == ARTICLE_DOI
    assert lead["year"] == "1956"
    assert lead["original_literature_reference"] == "Complete primary citation."
    assert lead["source_mean_tm_raw_unmapped"] == "0.999"
    assert lead["local_file_hash"] == report["source_csv_sha256"]
    assert lead["retrieval_date"] == "2026-07-19"
    raw = json.loads(lead["source_text"])
    assert raw["si"] == "0"
    assert raw["mean.tm"] == "0.999"
    assert set(raw) == {
        "year",
        "reference",
        "family",
        "genus",
        "species",
        "m.d.g",
        "latitude",
        "hemisphere",
        "biome",
        "growth",
        "life.history",
        "si",
        "mean.tm",
        "References",
    }

    sources = pd.read_csv(output / "source_evidence.csv", dtype=str).fillna("")
    assert list(sources.columns) == SOURCE_COLUMNS
    gamma = sources.loc[sources["accepted_species"].eq("Gamma delta")].iloc[0]
    assert gamma["source_si_value"] == ""
    assert gamma["source_mean_tm_raw_unmapped"] == "0.000"
    assert gamma["mapping_status"] == "no_explicit_si_category_mean_tm_not_mapped"

    checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    dryad = checkpoint.loc[checkpoint["provider"].eq(PROVIDER)].set_index("species")
    assert len(dryad) == len(original)
    assert dryad.loc["Alpha beta", "status"] == "provider_pass_complete_with_candidates"
    assert dryad.loc["Alpha beta", "candidate_count"] == "1"
    assert dryad.loc["Gamma delta", "status"] == "provider_no_result"
    assert dryad.loc["Covered species", "status"] == "provider_seed_covered"
    assert set(dryad["terminal"].str.lower()) == {"true"}
    assert set(dryad["accepted_record_count"]) == {"0"}
    assert report["mean_tm_mapped"] is False
    assert report["biological_values_materialized"] is False
    assert report["all_original_provider_keys_terminal"] is True


def test_reviewed_dryad_candidate_builds_and_materializes_strictly(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    queue, search, providers, search_path = _review_fixture(tmp_path, monkeypatch)
    candidates, review = _build_reviewed(queue, search, providers, search_path)
    candidate = candidates.iloc[0]
    assert candidate["provider"] == PROVIDER
    assert candidate["source_type"] == "dryad_biogeo_mating_si_category"
    assert candidate["queue_family"] == "Fabaceae"
    assert candidate["source_family"] == "Fabaceae"
    assert candidate["provisional_candidate_value"] == "SI"
    assert candidate["source_mean_tm_raw_unmapped"] == "0.84"

    result = materialize_reviewed_cache(
        queue=queue,
        baseline_checkpoint=_axis(queue),
        provider_checkpoint=providers,
        candidates=candidates,
        review_payload=review,
    )
    record = result["records"].iloc[0]
    assert record["species"] == "Alpha beta"
    assert record["trait"] == "self_incompatibility"
    assert record["value"] == "SI"
    assert record["provider"] == PROVIDER
    assert record["DOI"] == DATASET_DOI
    assert record["year"] == "2018"
    assert json.loads(record["raw_text"])["si"] == "0"
    dryad = (
        result["provider_checkpoint"]
        .loc[result["provider_checkpoint"]["provider"].eq(PROVIDER)]
        .iloc[0]
    )
    assert dryad["status"] == "accepted_evidence"
    assert int(dryad["accepted_record_count"]) == 1
    assert result["three_axis_checkpoint"].loc[0, "reproductive_assurance_direct"]


@pytest.mark.parametrize(
    "mutation",
    ["si", "value", "species", "family", "hash", "source_type"],
)
def test_review_builder_rejects_altered_dryad_contract(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
):
    queue, search, providers, search_path = _review_fixture(tmp_path, monkeypatch)
    value = "SI"
    raw = json.loads(search.loc[0, "source_text"])
    if mutation == "si":
        raw["si"] = "1"
        search.loc[0, "source_text"] = json.dumps(raw, ensure_ascii=False, separators=(",", ":"))
    elif mutation == "value":
        value = "SC"
    elif mutation == "species":
        raw["genus"] = "Gamma"
        search.loc[0, "source_text"] = json.dumps(raw, ensure_ascii=False, separators=(",", ":"))
    elif mutation == "family":
        raw["family"] = "Wrongaceae"
        search.loc[0, "source_text"] = json.dumps(raw, ensure_ascii=False, separators=(",", ":"))
        search.loc[0, "source_family"] = "Wrongaceae"
    elif mutation == "hash":
        search.loc[0, "local_file_hash"] = "0" * 64
    elif mutation == "source_type":
        search.loc[0, "source_type"] = "curated_breeding_system_database"

    with pytest.raises(ValueError):
        _build_reviewed(queue, search, providers, search_path, value=value)


@pytest.mark.parametrize(
    "mutation",
    ["si", "value", "species", "family", "hash", "source_type"],
)
def test_materializer_revalidates_altered_dryad_contract(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
):
    queue, search, providers, search_path = _review_fixture(tmp_path, monkeypatch)
    candidates, review = _build_reviewed(queue, search, providers, search_path)
    candidates = candidates.copy()
    review = json.loads(json.dumps(review))
    raw = json.loads(candidates.loc[0, "source_text"])
    if mutation == "si":
        raw["si"] = "1"
        candidates.loc[0, "source_text"] = json.dumps(
            raw, ensure_ascii=False, separators=(",", ":")
        )
    elif mutation == "value":
        review["accepted"][0]["accepted_value"] = "SC"
    elif mutation == "species":
        raw["genus"] = "Gamma"
        candidates.loc[0, "source_text"] = json.dumps(
            raw, ensure_ascii=False, separators=(",", ":")
        )
    elif mutation == "family":
        candidates.loc[0, "source_family"] = "Wrongaceae"
    elif mutation == "hash":
        candidates.loc[0, "local_file_hash"] = "0" * 64
    elif mutation == "source_type":
        candidates.loc[0, "source_type"] = "europe_pmc_title_abstract"

    with pytest.raises(ValueError, match="exact pinned"):
        materialize_reviewed_cache(
            queue=queue,
            baseline_checkpoint=_axis(queue),
            provider_checkpoint=providers,
            candidates=candidates,
            review_payload=review,
        )


def test_mean_tm_without_explicit_si_cannot_build_reviewed_candidate(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    queue, search, providers, search_path = _review_fixture(
        tmp_path,
        monkeypatch,
        si="",
        mean_tm="0.999",
    )
    assert search.loc[0, "source_si_raw"] == ""
    assert search.loc[0, "source_mean_tm_raw_unmapped"] == "0.999"
    with pytest.raises(ValueError, match="exact pinned"):
        _build_reviewed(queue, search, providers, search_path)


def test_resume_skips_terminal_keys_without_duplicate_candidates(tmp_path: Path):
    original = _queue(("Alpha beta", "Fabaceae"))
    rows = [_source_row("Alpha", "beta", "Fabaceae", si="1")]
    first, output, checkpoint_path = _run(
        tmp_path,
        original=original,
        current=original,
        rows=rows,
    )
    checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    source = tmp_path / "Database_BioGeo_mating.csv"
    second = collect_dryad_biogeo_candidates(
        original_queue=original,
        current_queue=original,
        provider_checkpoint=checkpoint,
        provider_checkpoint_path=checkpoint_path,
        source_csv=source,
        output_dir=output,
        retrieval_date="2026-07-19",
        expected_csv_sha256=_file_sha256(source),
    )
    leads = pd.read_csv(output / "candidate_leads.csv", dtype=str).fillna("")
    assert first["scheduled_species"] == 1
    assert second["scheduled_species"] == 0
    assert len(leads) == 1
    resumed = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    row = resumed.loc[resumed["provider"].eq(PROVIDER)].iloc[0]
    assert row["attempts"] == "1"


@pytest.mark.parametrize(
    ("rows", "expected_error"),
    [
        (
            [
                _source_row("Alpha", "beta", "Fabaceae", si="0"),
                _source_row("Alpha", "beta", "Fabaceae", si="1"),
            ],
            "conflicting_si_values",
        ),
        ([_source_row("Alpha", "beta", "Wrongaceae", si="0")], "source_family_mismatch"),
        ([_source_row("Alpha", "beta", "Fabaceae", si="0.0")], "unsupported_si_value"),
    ],
)
def test_ambiguous_or_unsupported_rows_fail_closed(
    tmp_path: Path,
    rows: list[dict[str, str]],
    expected_error: str,
):
    queue = _queue(("Alpha beta", "Fabaceae"))
    report, output, checkpoint_path = _run(
        tmp_path,
        original=queue,
        current=queue,
        rows=rows,
    )
    leads = pd.read_csv(output / "candidate_leads.csv", dtype=str).fillna("")
    assert leads.empty
    errors = pd.read_csv(output / "lookup_errors.csv", dtype=str).fillna("")
    assert len(errors) == 1
    assert expected_error in errors.loc[0, "message"]
    checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    dryad = checkpoint.loc[checkpoint["provider"].eq(PROVIDER)].iloc[0]
    assert dryad["status"] == "permanent_failure"
    assert dryad["terminal"].lower() == "true"
    assert dryad["candidate_count"] == "0"
    assert report["new_permanent_failures"] == 1


def test_infraspecific_source_name_is_not_transferred_to_binomial(tmp_path: Path):
    queue = _queue(("Alpha beta", "Fabaceae"))
    rows = [_source_row("Alpha", "beta_subsp_minor", "Fabaceae", si="1")]
    _report, output, checkpoint_path = _run(
        tmp_path,
        original=queue,
        current=queue,
        rows=rows,
    )
    assert pd.read_csv(output / "candidate_leads.csv").empty
    assert pd.read_csv(output / "source_evidence.csv").empty
    checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    row = checkpoint.loc[checkpoint["provider"].eq(PROVIDER)].iloc[0]
    assert row["status"] == "provider_no_result"


def test_hash_mismatch_stops_before_checkpoint_creation(tmp_path: Path):
    queue = _queue(("Alpha beta", "Fabaceae"))
    source = tmp_path / "source.csv"
    _write_source(source, [_source_row("Alpha", "beta", "Fabaceae", si="1")])
    checkpoint_path = tmp_path / "out" / "provider_checkpoint.csv"
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        collect_dryad_biogeo_candidates(
            original_queue=queue,
            current_queue=queue,
            provider_checkpoint=pd.DataFrame(columns=CHECKPOINT_COLUMNS),
            provider_checkpoint_path=checkpoint_path,
            source_csv=source,
            output_dir=tmp_path / "out",
            retrieval_date="2026-07-19",
            expected_csv_sha256="0" * 64,
        )
    assert not checkpoint_path.exists()


def test_checkpoint_requires_current_queue_to_be_original_subset():
    original = _queue(("Alpha beta", "Fabaceae"))
    current = _queue(("Gamma delta", "Lamiaceae"))
    with pytest.raises(ValueError, match="outside original queue"):
        ensure_dryad_checkpoint_rows(
            pd.DataFrame(columns=CHECKPOINT_COLUMNS),
            original,
            current,
        )
