from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from analysis.acquire_wfo_global_flora_traits_20260810 import (
    AUDIT_SEED,
    SOURCE_META,
    deterministic_audit_sample,
)
from analysis.acquire_wfo_kew_africa_traits_20260810 import BACKBONE_SHA256
from island_v2.open_web_common import reviewed_source_package_evidence

PACKAGE = Path(
    "data/v2/staging/traits/direct_llm_pilot/"
    "20260810_wfo_global_six_source_acquisition"
)


def test_committed_holdout_is_reproducible_and_independent() -> None:
    evidence = pd.read_csv(PACKAGE / "wfo_global_six_evidence.csv.gz", dtype=str).fillna("")
    reviewed = pd.read_csv(
        PACKAGE / "wfo_global_six_manual_audit_200.csv", dtype=str
    ).fillna("")

    regenerated = deterministic_audit_sample(evidence)

    assert AUDIT_SEED == "wfo-global-six-floras-holdout-audit-20260810-v1"
    assert regenerated["candidate_id"].tolist() == reviewed["candidate_id"].tolist()
    assert reviewed["candidate_id"].is_unique
    assert reviewed["reviewer"].ne("").all()
    assert reviewed["reviewed_at_utc"].ne("").all()
    assert reviewed["audit_reason"].ne("").all()


def test_committed_package_scales_only_passing_traits() -> None:
    evidence = pd.read_csv(PACKAGE / "wfo_global_six_evidence.csv.gz", dtype=str).fillna("")
    audit = pd.read_csv(
        PACKAGE / "wfo_global_six_manual_audit_200.csv", dtype=str
    ).fillna("")

    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)

    assert summary["package_gate_passed"] is True
    assert summary["approval_unit"] == "trait_name"
    assert summary["overall_precision_gate_applied"] is False
    assert summary["reviewed"] == 200
    assert summary["accepted_correct"] == 187
    assert summary["precision"] == 0.935
    assert summary["cultivar_contamination_rate"] == 0.0
    assert set(scopes.loc[scopes["production_approved"], "trait_name"]) == {
        "floral_form",
        "floral_symmetry",
        "flower_size_class",
        "inflorescence_display",
    }
    assert not scopes.loc[
        scopes["trait_name"].isin({"flower_primary_color", "tube_depth_class"}),
        "production_approved",
    ].any()
    assert len(selected) == 5_434
    assert selected.drop_duplicates(["accepted_species", "trait_name"]).shape[0] == 4_742
    assert selected.drop_duplicates(["accepted_species", "axis"]).shape[0] == 3_788
    assert selected["accepted_species"].nunique() == 3_788


def test_committed_summary_pins_six_archives_and_fixed_universe() -> None:
    summary = json.loads(
        (PACKAGE / "wfo_global_six_source_package_summary.json").read_text(
            encoding="utf-8"
        )
    )
    review = json.loads(
        (PACKAGE / "wfo_global_six_review_summary.json").read_text(encoding="utf-8")
    )

    assert summary["fixed_universe_species"] == 106_295
    assert summary["backbone"]["sha256"] == BACKBONE_SHA256
    assert summary["source_archives"] == SOURCE_META
    assert summary["matched_description_rows"] == 22_981
    assert summary["incremental_evidence_rows"] == 6_939
    assert summary["incremental_species_trait"] == 6_099
    assert review["approved_traits"] == [
        "floral_form",
        "floral_symmetry",
        "flower_size_class",
        "inflorescence_display",
    ]
    assert review["selected_species_trait"] == 4_742
    assert review["selected_species_axis"] == 3_788


def test_committed_source_package_manifest_is_exact() -> None:
    entries = {}
    for line in (PACKAGE / "file_manifest.sha256").read_text(
        encoding="ascii"
    ).splitlines():
        digest, filename = line.split("  ", maxsplit=1)
        entries[filename] = digest

    files = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in PACKAGE.iterdir()
        if path.is_file() and path.name != "file_manifest.sha256"
    }
    assert entries == files
