from __future__ import annotations

from pathlib import Path

import pandas as pd

from island_v2 import integrated_trait_coverage as coverage

DIRECT_COLUMNS = [
    "accepted_species",
    "matched_source_name",
    "trait_name",
    "reported_value",
    "raw_candidate_value",
    "source_url",
    "provider",
    "source_title",
    "excerpt",
    "retrieval_date",
    "stable_source_id",
    "evidence_scope",
    "evidence_quality",
]
CANDIDATE_COLUMNS = [
    "evidence_id",
    "source_name",
    "source_record_id",
    "accepted_species",
    "name_match_method",
    "trait_name",
    "raw_value",
    "standardized_value",
    "candidate_kind",
    "evidence_scope",
    "source_citation",
    "source_url",
    "evidence_excerpt",
    "inference_rule",
    "confidence",
    "review_status",
]


def _direct_row(
    species: str,
    trait: str,
    value: str,
    quality: str,
    *,
    url: str = "https://example.org/record",
) -> dict[str, str]:
    return {
        "accepted_species": species,
        "matched_source_name": species,
        "trait_name": trait,
        "reported_value": "",
        "raw_candidate_value": value,
        "source_url": url,
        "provider": "flora_or_monograph",
        "source_title": "Example flora",
        "excerpt": f"{species}: {trait} = {value}",
        "retrieval_date": "2026-07-24",
        "stable_source_id": f"{species}:{trait}",
        "evidence_scope": "species_direct",
        "evidence_quality": quality,
    }


def _bulk_row(
    source_name: str,
    species: str,
    trait: str,
    value: str,
    *,
    citation: str,
    scope: str = "species_direct",
    inference_rule: str = "",
) -> dict[str, str]:
    return {
        "evidence_id": f"{source_name}:{species}:{trait}",
        "source_name": source_name,
        "source_record_id": f"{source_name}:{species}:{trait}",
        "accepted_species": species,
        "name_match_method": "accepted_name_exact",
        "trait_name": trait,
        "raw_value": value,
        "standardized_value": value,
        "candidate_kind": "source_backed",
        "evidence_scope": scope,
        "source_citation": citation,
        "source_url": "https://plants.usda.gov/",
        "evidence_excerpt": f"{species}: {trait} = {value}",
        "inference_rule": inference_rule,
        "confidence": "medium",
        "review_status": "pending",
    }


def _write_three_pass(
    root: Path, shard_count: int = 1
) -> tuple[Path, Path, Path]:
    high = root / "high"
    medium = root / "medium"
    unresolved = root / "unresolved"
    high.mkdir()
    medium.mkdir()
    unresolved.mkdir()
    high_rows = [
        _direct_row("Alpha one", "flower_primary_color", "white", "high")
    ]
    medium_rows = [
        _direct_row("Alpha one", "floral_form", "tubular", "medium")
    ]
    unresolved_rows = [
        {"accepted_species": "Alpha one", "axis": "reproductive_assurance"},
        {"accepted_species": "Beta two", "axis": "flower_colour"},
        {"accepted_species": "Beta two", "axis": "floral_structural_complexity"},
        {"accepted_species": "Beta two", "axis": "reproductive_assurance"},
    ]
    for shard in range(shard_count):
        pd.DataFrame(
            high_rows if shard == 0 else [], columns=DIRECT_COLUMNS
        ).to_csv(
            high / f"accepted_high_evidence-{shard}.csv.gz",
            index=False,
        )
        pd.DataFrame(
            medium_rows if shard == 0 else [], columns=DIRECT_COLUMNS
        ).to_csv(
            medium / f"accepted_medium_evidence-{shard}.csv.gz",
            index=False,
        )
        pd.DataFrame(
            unresolved_rows if shard == 0 else [],
            columns=["accepted_species", "axis"],
        ).to_csv(
            unresolved / f"unresolved_after_pass3-{shard}.csv.gz",
            index=False,
        )
    # Production artifacts use a fixed basename inside distinct shard folders.
    for directory, prefix, basename in (
        (high, "accepted_high_evidence-", "accepted_high_evidence.csv.gz"),
        (medium, "accepted_medium_evidence-", "accepted_medium_evidence.csv.gz"),
        (unresolved, "unresolved_after_pass3-", "unresolved_after_pass3.csv.gz"),
    ):
        for path in list(directory.glob(f"{prefix}*.csv.gz")):
            shard_dir = directory / path.stem.split("-")[-1].split(".")[0]
            shard_dir.mkdir()
            path.rename(shard_dir / basename)
    return high, medium, unresolved


def _write_validated_low(root: Path) -> Path:
    low = root / "low"
    columns = [
        "accepted_species",
        "genus",
        "axis",
        "trait_name",
        "normalized_value",
        "evidence_quality",
        "evidence_scope",
        "inference_method",
        "n_direct_species",
        "masked_accuracy",
        "family_inference_used",
        "global_fallback_used",
    ]
    rows = {
        "flower_colour": [
            {
                "accepted_species": "Beta two",
                "genus": "Beta",
                "axis": "flower_colour",
                "trait_name": "flower_primary_color",
                "normalized_value": "white",
                "evidence_quality": "low",
                "evidence_scope": "genus_consensus",
                "inference_method": "validated_genus_consensus",
                "n_direct_species": "4",
                "masked_accuracy": "1.0",
                "family_inference_used": "False",
                "global_fallback_used": "False",
            }
        ],
        "floral_structural_complexity": [],
        "reproductive_assurance": [
            {
                "accepted_species": "Alpha one",
                "genus": "Alpha",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "normalized_value": "SI",
                "evidence_quality": "low",
                "evidence_scope": "genus_consensus",
                "inference_method": "validated_genus_consensus",
                "n_direct_species": "4",
                "masked_accuracy": "1.0",
                "family_inference_used": "False",
                "global_fallback_used": "False",
            }
        ],
    }
    for axis, axis_rows in rows.items():
        path = low / axis
        path.mkdir(parents=True)
        pd.DataFrame(axis_rows, columns=columns).to_csv(
            path / "accepted_validated_low_evidence.csv.gz", index=False
        )
    return low


def _write_all_master(root: Path) -> Path:
    all_master = root / "all-master"
    eol = all_master / "eol_traitbank"
    usda = all_master / "usda_plants_current"
    meyer = all_master / "meyer_2026_reproductive_traits"
    for path in (eol, usda, meyer):
        path.mkdir(parents=True)
    usda_citation = (
        "The PLANTS Database, United States Department of Agriculture, "
        "Natural Resources Conservation Service."
    )
    eol_row = _bulk_row(
        "eol_traitbank_all_zenodo_13305577",
        "Beta two",
        "flower_primary_color",
        "white",
        citation=usda_citation,
    )
    usda_row = _bulk_row(
        "usda_plants_current",
        "Beta two",
        "flower_primary_color",
        "white",
        citation=usda_citation,
    )
    rejected = _bulk_row(
        "usda_plants_current",
        "Beta two",
        "flower_primary_color",
        "red_pink",
        citation=usda_citation,
        scope="family_inference",
        inference_rule="family majority",
    )
    meyer_row = _bulk_row(
        "meyer_galloway_eckert_2026_dryad",
        "Beta two",
        "self_incompatibility",
        "SC",
        citation="Meyer et al. 2026. https://doi.org/10.5061/dryad.cc2fqz6hr",
    )
    pd.DataFrame([eol_row], columns=CANDIDATE_COLUMNS).to_csv(
        eol / "trait_candidates.csv", index=False
    )
    pd.DataFrame([usda_row, rejected], columns=CANDIDATE_COLUMNS).to_csv(
        usda / "trait_candidates.csv", index=False
    )
    pd.DataFrame([meyer_row], columns=CANDIDATE_COLUMNS).to_csv(
        meyer / "trait_candidates.csv.gz", index=False
    )
    pd.DataFrame(
        columns=[
            "species",
            "field",
            "value",
            "source_kind",
            "source_url",
            "source_record_id",
            "source_excerpt",
            "evidence_type",
            "confidence",
            "source_backed",
        ]
    ).to_csv(all_master / "all_species_trait_evidence.csv.gz", index=False)
    return all_master


def test_bulk_pending_rows_need_strict_contract_and_lineage_dedup(tmp_path: Path) -> None:
    all_master = _write_all_master(tmp_path)
    evidence, audit = coverage.load_bulk_artifact(all_master, {"sources": []})
    assert audit["accepted"].value_counts().to_dict() == {"True": 3, "False": 1}
    assert audit.loc[audit["accepted"].eq("False"), "reason"].item() == (
        "not_species_or_synonym_direct"
    )
    colour = evidence.loc[evidence["axis"].eq("flower_colour")]
    assert set(colour["source_group"]) == {"eol_traitbank", "usda_plants"}
    assert colour["source_lineage"].nunique() == 1
    deduped, duplicates = coverage.dedupe_source_lineages(evidence)
    assert len(deduped) == 2
    assert len(duplicates) == 2


def test_lineage_dedup_preserves_distinct_traits_on_the_same_axis() -> None:
    rows = [
        coverage._make_evidence(
            species="Alpha one",
            trait_name=trait,
            value=value,
            quality="medium",
            source_group="efloras",
            provider="efloras:2",
            url="https://example.org/same-treatment",
            record_id=f"alpha:{trait}",
            citation="Example flora",
            excerpt=f"{trait}={value}",
            scope="species_direct",
            name_match="accepted_name_exact",
            source_file=Path("artifact.csv"),
            contract="test",
            manifest={"sources": []},
        )
        for trait, value in (
            ("floral_form", "tubular"),
            ("floral_symmetry", "actinomorphic"),
        )
    ]
    deduped, duplicates = coverage.dedupe_source_lineages(
        pd.DataFrame(rows, columns=coverage.EVIDENCE_COLUMNS)
    )
    assert len(deduped) == 2
    assert duplicates.empty


def test_generic_bulk_candidate_preserves_explicit_original_source_lineage(
    tmp_path: Path,
) -> None:
    path = tmp_path / "trait_candidates.csv.gz"
    record = _bulk_row(
        "ferrer_2024_si_database",
        "Alpha one",
        "self_incompatibility",
        "SI",
        citation="Original study, 2001",
    )
    record["source_lineage"] = "citation-set:abc123"
    pd.DataFrame([record]).to_csv(path, index=False)

    evidence, audit = coverage.load_generic_bulk_candidates(path, {"sources": []})
    assert audit["accepted"].tolist() == ["True"]
    assert evidence["source_group"].tolist() == ["ferrer_2024_si"]
    assert evidence["source_lineage"].tolist() == ["citation-set:abc123"]
    assert evidence["lineage_method"].tolist() == [
        "candidate_explicit_source_lineage"
    ]


def test_aggregate_reports_before_after_low_upgrade_and_exact_denominator(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setattr(coverage, "THREE_PASS_SHARD_COUNT", 1)
    high, medium, unresolved = _write_three_pass(tmp_path)
    low = _write_validated_low(tmp_path)
    all_master = _write_all_master(tmp_path)
    targeted = tmp_path / "targeted"
    targeted.mkdir()
    pd.DataFrame(
        [_direct_row("Beta two", "floral_form", "rotate", "medium")],
        columns=DIRECT_COLUMNS,
    ).to_csv(targeted / "targeted_new_evidence.csv.gz", index=False)

    summary, frames = coverage.aggregate(
        three_pass_high_dir=high,
        three_pass_medium_dir=medium,
        three_pass_unresolved_dir=unresolved,
        expanded_artifact_dirs=[],
        targeted_dir=targeted,
        validated_low_dir=low,
        all_master_dir=all_master,
        source_manifest={"sources": []},
        expected_species=2,
    )

    assert summary["denominator"]["species_axis"] == 6
    assert summary["before_bulk_integration"]["quality_counts"] == {
        "high": 1,
        "medium": 2,
        "low": 2,
    }
    assert summary["after_bulk_integration"]["quality_counts"] == {
        "high": 1,
        "medium": 4,
        "low": 1,
    }
    assert summary["net_change"]["filled_species_axis"] == 1
    assert summary["net_change"]["validated_low_upgraded_to_direct"] == 1
    assert summary["after_bulk_integration"]["species_by_filled_axis_count"] == {
        "0": 0,
        "1": 0,
        "2": 0,
        "3": 2,
    }
    assert summary["after_bulk_integration"]["unresolved_species_axis"] == 0
    assert summary["checks"]["family_inference_absent"]
    assert summary["checks"]["global_fallback_absent"]
    assert len(frames["validated_low_upgrades"]) == 1


def test_validated_low_rejects_family_or_global_fallback(tmp_path: Path) -> None:
    low = _write_validated_low(tmp_path)
    path = low / "flower_colour" / "accepted_validated_low_evidence.csv.gz"
    frame = pd.read_csv(path, dtype=str).fillna("")
    frame.loc[0, "family_inference_used"] = "True"
    frame.to_csv(path, index=False)
    try:
        coverage.load_validated_low(low, {"sources": []})
    except ValueError as exc:
        assert "family inference" in str(exc)
    else:
        raise AssertionError("family inference should be rejected")


def test_workflow_pins_scoped_pilot_artifact_instead_of_raw_range_size() -> None:
    workflow = Path(
        ".github/workflows/build-integrated-trait-coverage.yml"
    ).read_text(encoding="utf-8")
    verification = workflow.split(
        "- name: Verify artifact-confirmed expanded-wave ranges",
        1,
    )[1].split("- name: Build integrated species-axis ledger", 1)[0]

    assert "expected_pilot_count = 715" in verification
    assert "9be8660b108a5d0926c157ba063135e0e" in verification
    assert "pilot targeted-master and shard plan disagree" in verification
    assert "does not contain 2,000 target species" not in verification
