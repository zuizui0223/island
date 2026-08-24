from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd
import pytest

from analysis.prepare_palmpedia_archive_source_package_20260813 import build
from island_v2.open_web_common import reviewed_source_package_evidence

ROOT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "palmpedia_archive_checkpoint_20260813"
)


def test_palmpedia_package_passes_only_reviewed_trait_gates() -> None:
    evidence = pd.read_csv(
        ROOT / "palmpedia_archive_source_package_evidence_20260813.csv.gz",
        dtype=str,
    ).fillna("")
    audit = pd.read_csv(
        ROOT / "palmpedia_archive_source_package_audit_20260813.csv",
        dtype=str,
    ).fillna("")
    pages = pd.read_csv(
        ROOT / "palmpedia_archive_manual_page_audit_200_20260813.csv",
        dtype=str,
    ).fillna("")

    selected, scopes, gate = reviewed_source_package_evidence(evidence, audit)

    assert len(evidence) == 243
    assert evidence["source_record_id"].is_unique
    assert len(audit) == 243
    assert audit["candidate_id"].is_unique
    assert gate["package_gate_passed"] is True
    assert gate["precision"] == pytest.approx(239 / 243)
    assert gate["cultivar_contamination_rate"] == 0.0
    assert gate["approved_traits"] == [
        "flower_primary_color",
        "flower_size_class",
        "inflorescence_display",
    ]
    assert len(selected) == 231
    assert selected["accepted_species"].nunique() == 151
    assert selected[["accepted_species", "trait_name"]].drop_duplicates().shape[0] == 177
    assert selected[["accepted_species", "axis"]].drop_duplicates().shape[0] == 168
    assert selected["quality"].eq("medium").all()
    assert selected["evidence_scope"].eq("species_direct").all()
    assert selected["source_url"].str.startswith("https://web.archive.org/").all()
    assert not set(selected["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )

    form_scope = scopes.loc[scopes["trait_name"].eq("floral_form")].iloc[0]
    assert form_scope["reviewed"] == 10
    assert form_scope["accepted_correct"] == 8
    assert form_scope["precision"] == pytest.approx(0.8)
    assert bool(form_scope["production_approved"]) is False

    assert len(pages) == 200
    assert pages["accepted_species"].is_unique
    assert pages["identity_correct"].eq("true").all()
    assert pages["audit_stratum"].value_counts().to_dict() == {
        "trait_candidate_page": 154,
        "no_supported_candidate_page": 46,
    }


def test_reviewed_multistate_values_and_known_rejections_are_preserved() -> None:
    evidence = pd.read_csv(
        ROOT / "palmpedia_archive_source_package_evidence_20260813.csv.gz",
        dtype=str,
    ).fillna("")
    audit = pd.read_csv(
        ROOT / "palmpedia_archive_source_package_audit_20260813.csv",
        dtype=str,
    ).fillna("")

    def values(species: str, trait: str) -> set[str]:
        return set(
            evidence.loc[
                evidence["accepted_species"].eq(species)
                & evidence["trait_name"].eq(trait),
                "normalized_value",
            ]
        )

    assert values("Dypsis onilahensis", "flower_primary_color") == {
        "green_brown_inconspicuous|white|yellow_orange"
    }
    assert values("Dypsis nossibensis", "flower_primary_color") == {
        "blue_purple|white"
    }
    assert values("Saribus papuanus", "inflorescence_display") == {
        "few_flowered|solitary"
    }

    rejected = audit.loc[audit["accepted_correct"].eq("false")]
    assert len(rejected) == 4
    assert set(rejected["accepted_species"]) == {
        "Acanthophoenix crinita",
        "Astrocaryum vulgare",
        "Roystonea princeps",
        "Sabal antillensis",
    }
    assert rejected["audit_reason"].str.contains(
        "calyx|anther|fruit-stage", regex=True
    ).all()


def test_palmpedia_package_rebuild_is_byte_reproducible(tmp_path: Path) -> None:
    manifest = json.loads(
        (ROOT / "palmpedia_archive_checkpoint_manifest_20260813.json").read_text(
            encoding="utf-8"
        )
    )
    for name in manifest["input_hashes"]:
        shutil.copy2(ROOT / name, tmp_path / name)

    rebuilt = build(tmp_path)

    assert rebuilt["input_hashes"] == manifest["input_hashes"]
    assert rebuilt["output_hashes"] == manifest["output_hashes"]
    assert rebuilt["review"]["source_package_gate"] == manifest["review"][
        "source_package_gate"
    ]
