from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from island_v2.materialize_trait_artifacts import materialize


def _candidate(source: str, *, excerpt_column: str = "evidence_excerpt") -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "accepted_species": "Aaa bbb",
                "trait_name": "flower_primary_color",
                "standardized_value": "white",
                "source_url": f"https://example.test/{source}",
                "source_record_id": f"{source}:1",
                excerpt_column: "Flowers white.",
            }
        ]
    )


def _artifact_dirs(tmp_path: Path) -> tuple[Path, Path]:
    eol = tmp_path / "eol"
    eol_path = eol / "ingested" / "bulk_trait_candidates_release.csv"
    eol_path.parent.mkdir(parents=True)
    _candidate("eol").to_csv(eol_path, index=False)

    austraits = tmp_path / "austraits"
    austraits_path = (
        austraits
        / "home"
        / "runner"
        / "work"
        / "island"
        / "island"
        / "data"
        / "v2"
        / "staging"
        / "traits"
        / "bulk"
        / "austraits_all"
        / "trait_candidates.csv.gz"
    )
    austraits_path.parent.mkdir(parents=True)
    _candidate("austraits").to_csv(austraits_path, index=False, compression="gzip")
    return eol, austraits


def test_materialize_copies_two_bulk_sources_and_records_evidence(tmp_path: Path) -> None:
    eol, austraits = _artifact_dirs(tmp_path)
    staging = tmp_path / "staging"
    output = tmp_path / "output"

    manifest = materialize(
        eol_artifact_dir=eol,
        austraits_artifact_dir=austraits,
        staging_root=staging,
        output_dir=output,
    )

    assert set(manifest["sources"]) == {"eol_traitbank", "austraits_all"}
    assert (staging / "eol_traitbank" / "trait_candidates.csv").exists()
    assert (staging / "bulk" / "austraits_all" / "trait_candidates.csv.gz").exists()
    assert manifest["sources"]["eol_traitbank"]["n_species"] == 1
    assert set(manifest["sources"]["eol_traitbank"]["empty_evidence_fields"].values()) == {0}
    persisted = json.loads((output / "materialization_manifest.json").read_text(encoding="utf-8"))
    assert persisted["source_count"] == 2


def test_materialize_optionally_adds_combined_efloras_candidates(tmp_path: Path) -> None:
    eol, austraits = _artifact_dirs(tmp_path)
    efloras = tmp_path / "efloras"
    efloras.mkdir()
    _candidate("efloras", excerpt_column="source_excerpt").rename(
        columns={"standardized_value": "candidate_value"}
    ).to_csv(efloras / "all_trait_candidates.csv", index=False)

    manifest = materialize(
        eol_artifact_dir=eol,
        austraits_artifact_dir=austraits,
        efloras_artifact_dir=efloras,
        staging_root=tmp_path / "staging",
        output_dir=tmp_path / "output",
    )

    assert set(manifest["sources"]) == {"eol_traitbank", "austraits_all", "efloras"}
    assert manifest["sources"]["efloras"]["excerpt_column"] == "source_excerpt"


def test_materialize_optionally_adds_gift_candidates(tmp_path: Path) -> None:
    eol, austraits = _artifact_dirs(tmp_path)
    gift = tmp_path / "gift"
    gift_path = (
        gift
        / "home"
        / "runner"
        / "work"
        / "island"
        / "island"
        / "data"
        / "v2"
        / "staging"
        / "traits"
        / "bulk"
        / "gift_v3_2"
        / "bulk_trait_candidates_gift_v3_2_direct.csv"
    )
    gift_path.parent.mkdir(parents=True)
    _candidate("gift").to_csv(gift_path, index=False)

    staging = tmp_path / "staging"
    manifest = materialize(
        eol_artifact_dir=eol,
        austraits_artifact_dir=austraits,
        gift_artifact_dir=gift,
        staging_root=staging,
        output_dir=tmp_path / "output",
    )

    assert set(manifest["sources"]) == {
        "eol_traitbank",
        "austraits_all",
        "gift_v3_2_direct",
    }
    assert manifest["sources"]["gift_v3_2_direct"]["n_species"] == 1
    assert (staging / "bulk" / "gift_v3_2" / "bulk_trait_candidates_gift_v3_2_direct.csv").exists()


def test_materialize_rejects_empty_source_evidence(tmp_path: Path) -> None:
    eol, austraits = _artifact_dirs(tmp_path)
    path = eol / "ingested" / "bulk_trait_candidates_release.csv"
    broken = _candidate("eol")
    broken.loc[0, "source_url"] = ""
    broken.to_csv(path, index=False)

    with pytest.raises(ValueError, match="empty evidence fields"):
        materialize(
            eol_artifact_dir=eol,
            austraits_artifact_dir=austraits,
            staging_root=tmp_path / "staging",
            output_dir=tmp_path / "output",
        )
