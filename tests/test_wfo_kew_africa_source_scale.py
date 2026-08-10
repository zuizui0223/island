from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from analysis.acquire_wfo_kew_africa_traits_20260810 import (
    AUDIT_SEED,
    BACKBONE_SHA256,
    SOURCE_META,
    _lineage_columns,
    _unique_identity_map,
    deterministic_audit_sample,
)
from island_v2.open_web_common import reviewed_source_package_evidence

PACKAGE = Path(
    "data/v2/staging/traits/direct_llm_pilot/"
    "20260810_wfo_kew_africa_source_acquisition"
)


def test_identity_map_is_fail_closed_on_cross_species_wfo_ids() -> None:
    identity = _unique_identity_map(
        {
            "wfo-1": [
                {
                    "accepted_species": "Alpha beta",
                    "matched_wfo_name": "Alpha beta",
                    "name_match_method": "accepted",
                },
                {
                    "accepted_species": "Alpha beta",
                    "matched_wfo_name": "Alpha oldname",
                    "name_match_method": "synonym",
                },
            ],
            "wfo-2": [
                {
                    "accepted_species": "Alpha beta",
                    "matched_wfo_name": "Alpha beta",
                    "name_match_method": "accepted",
                },
                {
                    "accepted_species": "Gamma delta",
                    "matched_wfo_name": "Gamma delta",
                    "name_match_method": "accepted",
                },
            ],
        }
    )

    assert identity["wfo-1"]["accepted_species"] == "Alpha beta"
    assert "wfo-2" not in identity


def test_lineage_collapses_identical_cross_provider_content() -> None:
    inventory = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha beta",
                "provider": "flora_a",
                "description_sha256": "same",
            },
            {
                "accepted_species": "Alpha beta",
                "provider": "flora_b",
                "description_sha256": "same",
            },
            {
                "accepted_species": "Alpha beta",
                "provider": "flora_b",
                "description_sha256": "different",
            },
        ]
    )

    resolved = _lineage_columns(inventory)

    assert resolved.loc[0, "source_lineage"] == resolved.loc[1, "source_lineage"]
    assert resolved.loc[2, "source_lineage"] != resolved.loc[0, "source_lineage"]
    assert resolved.loc[0, "lineage_method"].startswith(
        "accepted_species_normalized_content_fingerprint"
    )


def test_committed_holdout_is_reproducible_and_independent() -> None:
    evidence = pd.read_csv(PACKAGE / "wfo_kew_africa_evidence.csv.gz", dtype=str).fillna("")
    reviewed = pd.read_csv(
        PACKAGE / "wfo_kew_africa_manual_audit_200.csv", dtype=str
    ).fillna("")

    regenerated = deterministic_audit_sample(evidence)

    assert AUDIT_SEED == "wfo-kew-africa-holdout-audit-20260810-v1"
    assert regenerated["candidate_id"].tolist() == reviewed["candidate_id"].tolist()
    assert reviewed["candidate_id"].is_unique
    assert reviewed["reviewer"].ne("").all()
    assert reviewed["reviewed_at_utc"].ne("").all()


def test_committed_package_passes_manual_trait_gate() -> None:
    evidence = pd.read_csv(PACKAGE / "wfo_kew_africa_evidence.csv.gz", dtype=str).fillna("")
    audit = pd.read_csv(
        PACKAGE / "wfo_kew_africa_manual_audit_200.csv", dtype=str
    ).fillna("")

    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)

    assert summary["package_gate_passed"] is True
    assert summary["reviewed"] == 200
    assert summary["accepted_correct"] == 197
    assert summary["precision"] == 0.985
    assert summary["cultivar_contamination_rate"] == 0.0
    assert summary["selected_evidence_rows"] == 1_485
    assert set(scopes.loc[scopes["production_approved"], "trait_name"]) == {
        "floral_form",
        "flower_primary_color",
        "flower_size_class",
        "inflorescence_display",
    }
    assert selected.drop_duplicates(["accepted_species", "trait_name"]).shape[0] == 1_266
    assert selected.drop_duplicates(["accepted_species", "axis"]).shape[0] == 1_248
    assert selected["accepted_species"].nunique() == 1_078


def test_committed_summary_pins_source_archives_and_fixed_universe() -> None:
    summary = json.loads(
        (PACKAGE / "wfo_kew_africa_source_package_summary.json").read_text(
            encoding="utf-8"
        )
    )
    review = json.loads(
        (PACKAGE / "wfo_kew_africa_review_summary.json").read_text(encoding="utf-8")
    )

    assert summary["fixed_universe_species"] == 106_295
    assert summary["backbone"]["sha256"] == BACKBONE_SHA256
    assert summary["source_archives"] == SOURCE_META
    assert summary["identity"] == {
        "ambiguous_wfo_ids_excluded": 4_683,
        "exact_fixed_backbone_rows": 106_277,
        "resolved_fixed_species": 95_207,
        "resolved_wfo_ids": 333_723,
        "strict_synonym_cluster_rows": 232_857,
    }
    assert summary["matched_description_rows"] == 37_779
    assert summary["incremental_species_trait"] == 1_271
    assert review["selected_species_trait"] == 1_266
    assert review["approved_traits"] == [
        "floral_form",
        "flower_primary_color",
        "flower_size_class",
        "inflorescence_display",
    ]


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
