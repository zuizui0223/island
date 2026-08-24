from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS
from island_v2.official_flora_colour_recovery_checkpoint import (
    apply_reviewed_audit,
    build_checkpoint,
    high_precision_candidate,
    observed_colour_states,
)
from island_v2.open_web_common import reviewed_source_package_evidence


def _row(**patch: str) -> dict[str, str]:
    row = {column: "" for column in EVIDENCE_COLUMNS}
    row.update(
        {
            "accepted_species": "Examplea alba",
            "axis": "flower_colour",
            "trait_name": "flower_primary_color",
            "normalized_value": "red_pink|white",
            "quality": "high",
            "source_group": "latest_public_web",
            "source_provider": "flora_of_panama",
            "source_url": "https://legacy.tropicos.org/Name/00000001",
            "source_record_id": "upstream-1",
            "source_citation": "Flora of Panama, NYBG",
            "source_excerpt": "Petals white with pink margins.",
            "evidence_scope": "species_direct",
            "name_match_method": (
                "wfo_june_2026_exact_accepted_species_family_consistent"
            ),
            "source_lineage": "upstream-lineage-1",
            "lineage_method": "wfo_taxon_id",
            "source_run_id": "wfo-global-six-floras-20260810",
            "source_artifact": "wfo-global-six-official-dwca-20260810",
            "source_file": "source.csv.gz",
            "acceptance_contract": (
                "official_wfo_provider_dwca_exact_species_description_quote_v1"
            ),
        }
    )
    row.update(patch)
    return row


def test_colour_gate_keeps_complete_petaloid_state_sets() -> None:
    assert observed_colour_states(
        "Flowers cream; corolla lavender with pinkish and green markings."
    ) == {
        "blue_purple",
        "green_brown_inconspicuous",
        "red_pink",
        "white",
    }
    assert high_precision_candidate(_row())[0]
    assert not high_precision_candidate(
        _row(
            normalized_value="white",
            source_excerpt="Petals white with pink margins.",
        )
    )[0]


def test_colour_gate_rejects_non_target_organs_buds_and_dry_colours() -> None:
    bad = [
        ("green_brown_inconspicuous", "Fruiting perianth enclosing a brown carpel."),
        ("white", "Petals hyaline and densely covered with white hairs."),
        ("yellow_orange", "Flowers punctate with orange glands."),
        ("white|yellow_orange", "Flowers yellow, buds whitish."),
        ("green_brown_inconspicuous|yellow_orange", "Petals yellow, drying brown."),
        (
            "blue_purple|red_pink|white",
            "Corolla white, tube purplish; filaments as long as the corolla, pinkish.",
        ),
    ]
    for normalized_value, excerpt in bad:
        assert not high_precision_candidate(
            _row(normalized_value=normalized_value, source_excerpt=excerpt)
        )[0]


def test_colour_gate_repairs_known_pdf_line_breaks() -> None:
    assert observed_colour_states(
        "Flowers white; petals with a yel- low crest and laven- der markings."
    ) == {"blue_purple", "white", "yellow_orange"}


def test_checkpoint_excludes_completed_pairs_and_deduplicates_mirrors() -> None:
    rows = []
    for index in range(202):
        rows.append(
            _row(
                accepted_species=f"Examplea species{index:03d}",
                source_url=f"https://legacy.tropicos.org/Name/{index:08d}",
                source_record_id=f"upstream-{index}",
                source_excerpt="Petals white with pink margins.",
            )
        )
    rows.append(
        {
            **rows[-1],
            "source_provider": "nybg_floraneotropica",
            "source_record_id": "mirror-record",
        }
    )
    baseline = pd.DataFrame(
        [
            {
                "accepted_species": "Examplea species000",
                "trait_name": "flower_primary_color",
            }
        ]
    )
    evidence, audit, selection = build_checkpoint(
        pd.DataFrame(rows), baseline, pd.DataFrame(columns=baseline.columns)
    )

    assert len(evidence) == 201
    assert len(audit) == 200
    assert audit["source_url"].nunique() == 200
    assert set(evidence["quality"]) == {"high"}
    assert set(evidence["trait_name"]) == {"flower_primary_color"}
    assert (
        selection["selection_reason"].eq("already_acquired_direct_species_trait").sum()
        == 1
    )


def test_reviewed_audit_cannot_change_deterministic_sample() -> None:
    candidates = pd.DataFrame(
        [
            _row(
                accepted_species=f"Examplea species{index:03d}",
                source_url=f"https://legacy.tropicos.org/Name/{index:08d}",
                source_record_id=f"upstream-{index}",
            )
            for index in range(200)
        ]
    )
    _, template, _ = build_checkpoint(
        candidates,
        pd.DataFrame(columns=["accepted_species", "trait_name"]),
        pd.DataFrame(columns=["accepted_species", "trait_name"]),
    )
    reviewed = template.copy()
    reviewed["accepted_correct"] = "true"
    reviewed["cultivar_status"] = "wild_or_unspecified_official_flora_record"
    reviewed["reviewer"] = "reviewer"
    reviewed["reviewed_at_utc"] = "2026-08-11T00:00:00Z"
    reviewed["audit_reason"] = "accepted after source-text review"
    assert len(apply_reviewed_audit(template, reviewed)) == 200


def test_committed_official_flora_colour_package_is_reviewed_and_hashed() -> None:
    root = Path(
        "data/v2/staging/traits/open_web_pilot/"
        "official_flora_colour_recovery_checkpoint_20260811"
    )
    manifest = json.loads(
        (root / "official_flora_colour_manifest_20260811.json").read_text(
            encoding="utf-8"
        )
    )
    evidence = pd.read_csv(
        root / "official_flora_colour_evidence_20260811.csv.gz", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        root / "official_flora_colour_audit_200_20260811.csv", dtype=str
    ).fillna("")
    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)

    assert manifest["prior_formal_run_id"] == "31445281207"
    assert manifest["audit_review_complete"] is True
    assert manifest["audit_rows"] == manifest["audit_unique_urls"] == 200
    assert len(evidence) == len(selected) == 590
    assert evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0] == 539
    assert evidence[["accepted_species", "axis"]].drop_duplicates().shape[0] == 539
    assert summary["package_gate_passed"] is True
    assert summary["precision"] == 1.0
    assert summary["cultivar_contamination_rate"] == 0.0
    assert set(scopes.loc[scopes["production_approved"], "trait_name"]) == {
        "flower_primary_color"
    }
    for record in manifest["files"]:
        path = root / record["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"]
