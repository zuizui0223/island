from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS
from island_v2.open_web_common import reviewed_source_package_evidence
from island_v2.wfo_kew_africa_recovery_checkpoint import (
    apply_reviewed_audit,
    build_checkpoint,
    high_precision_candidate,
    observed_display_states,
    observed_form_states,
)


def _row(**patch: str) -> dict[str, str]:
    row = {column: "" for column in EVIDENCE_COLUMNS}
    row.update(
        {
            "accepted_species": "Examplea alba",
            "axis": "floral_structural_complexity",
            "trait_name": "floral_form",
            "normalized_value": "bell_campanulate",
            "quality": "high",
            "source_group": "latest_public_web",
            "source_provider": "kew_flora_tropical_east_africa",
            "source_url": "https://list.worldfloraonline.org/wfo-0000000001",
            "source_record_id": "upstream-1",
            "source_citation": "Flora of Tropical East Africa, Kew",
            "source_excerpt": "Corolla 5 mm long, campanulate and white.",
            "evidence_scope": "species_direct",
            "name_match_method": (
                "wfo_plant_list_exact_accepted_or_synonym_species_family_consistent"
            ),
            "source_lineage": "upstream-lineage-1",
            "lineage_method": "wfo_taxon_id",
            "source_run_id": "wfo-kew-africa-three-floras-20260810",
            "source_artifact": "wfo-kew-africa-official-dwca-20260810",
            "source_file": "source.csv.gz",
            "acceptance_contract": (
                "official_wfo_provider_dwca_exact_species_description_quote_v1"
            ),
        }
    )
    row.update(patch)
    return row


def test_form_gate_uses_nearest_governing_organ_and_complete_state_set() -> None:
    assert observed_form_states("Corolla broadly campanulate.") == {
        "bell_campanulate"
    }
    assert observed_form_states(
        "Male flowers clustered; receptacle-tube campanulate; petals yellow."
    ) == set()
    assert observed_form_states(
        "Corolla broadly campanulate, slightly 2-labiate."
    ) == {"bell_campanulate", "bilabiate"}
    assert not high_precision_candidate(
        _row(source_excerpt="Corolla broadly campanulate, slightly 2-labiate.")
    )[0]
    assert not high_precision_candidate(
        _row(
            normalized_value="funnel_trumpet",
            source_excerpt=(
                "Disk-florets with funnel-shaped corollas; ray-florets present."
            ),
        )
    )[0]
    assert high_precision_candidate(_row())[0]


def test_display_gate_preserves_explicit_secondary_states() -> None:
    assert observed_display_states("Inflorescence a 1-3-flowered cyme.") == {
        "few_flowered",
        "umbel_corymb",
    }
    assert observed_display_states(
        "Flowers in terminal heads which elongate into spikes."
    ) == {"composite_display", "raceme_spike_panicle"}
    assert not high_precision_candidate(
        _row(
            trait_name="inflorescence_display",
            normalized_value="solitary",
            source_excerpt="Flowers solitary or fasciculate in the leaf axils.",
        )
    )[0]
    assert not high_precision_candidate(
        _row(
            trait_name="inflorescence_display",
            normalized_value="umbel_corymb",
            source_excerpt="Flowers in few-many-flowered terminal cymes.",
        )
    )[0]
    assert observed_display_states("Inflorescence a (1-)3(-5)-flowered cyme.") == {
        "few_flowered",
        "umbel_corymb",
    }


def test_checkpoint_excludes_completed_pairs_and_deduplicates_mirrors() -> None:
    rows = []
    for index in range(210):
        species = f"Examplea species{index:03d}"
        trait = "floral_form" if index < 55 else "inflorescence_display"
        value = "bell_campanulate" if trait == "floral_form" else "umbel_corymb"
        excerpt = (
            "Corolla campanulate and white."
            if trait == "floral_form"
            else "Flowers in terminal cymes."
        )
        rows.append(
            _row(
                accepted_species=species,
                trait_name=trait,
                normalized_value=value,
                source_excerpt=excerpt,
                source_url=f"https://list.worldfloraonline.org/wfo-{index:010d}",
                source_record_id=f"upstream-{index}",
            )
        )
    rows.append(
        {
            **rows[-1],
            "source_provider": "kew_flora_zambesiaca_current",
            "source_record_id": "mirror-record",
        }
    )
    baseline = pd.DataFrame(
        [{"accepted_species": "Examplea species000", "trait_name": "floral_form"}]
    )
    evidence, audit, selection = build_checkpoint(
        pd.DataFrame(rows), baseline, pd.DataFrame(columns=baseline.columns)
    )

    assert len(evidence) == 209
    assert len(audit) == 200
    assert audit["source_url"].nunique() == 200
    assert set(evidence["quality"]) == {"high"}
    assert set(evidence["trait_name"]) == {
        "floral_form",
        "inflorescence_display",
    }
    assert (
        selection["selection_reason"].eq("already_acquired_direct_species_trait").sum()
        == 1
    )


def test_reviewed_audit_cannot_change_deterministic_sample() -> None:
    candidates = pd.DataFrame(
        [
            _row(
                accepted_species=f"Examplea species{index:03d}",
                trait_name=("floral_form" if index < 50 else "inflorescence_display"),
                normalized_value=(
                    "bell_campanulate" if index < 50 else "umbel_corymb"
                ),
                source_excerpt=(
                    "Corolla campanulate."
                    if index < 50
                    else "Flowers in terminal cymes."
                ),
                source_url=f"https://list.worldfloraonline.org/wfo-{index:010d}",
                source_record_id=f"upstream-{index}",
            )
            for index in range(200)
        ]
    )
    evidence, template, _ = build_checkpoint(
        candidates,
        pd.DataFrame(columns=["accepted_species", "trait_name"]),
        pd.DataFrame(columns=["accepted_species", "trait_name"]),
    )
    assert len(evidence) == 200
    reviewed = template.copy()
    reviewed["accepted_correct"] = "true"
    reviewed["cultivar_status"] = "wild_or_unspecified_official_flora_record"
    reviewed["reviewer"] = "reviewer"
    reviewed["reviewed_at_utc"] = "2026-08-11T00:00:00Z"
    reviewed["audit_reason"] = "accepted after source-text review"
    assert len(apply_reviewed_audit(template, reviewed)) == 200


def test_committed_wfo_africa_recovery_package_is_reviewed_and_hashed() -> None:
    root = Path(
        "data/v2/staging/traits/open_web_pilot/"
        "wfo_kew_africa_recovery_checkpoint_20260811"
    )
    manifest = json.loads(
        (root / "wfo_kew_africa_recovery_manifest_20260811.json").read_text(
            encoding="utf-8"
        )
    )
    evidence = pd.read_csv(
        root / "wfo_kew_africa_recovery_evidence_20260811.csv.gz", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        root / "wfo_kew_africa_recovery_audit_200_20260811.csv", dtype=str
    ).fillna("")
    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)

    assert manifest["prior_formal_run_id"] == "31440559943"
    assert manifest["audit_review_complete"] is True
    assert manifest["audit_rows"] == manifest["audit_unique_urls"] == 200
    assert len(evidence) == len(selected) == 450
    assert evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0] == 371
    assert evidence[["accepted_species", "axis"]].drop_duplicates().shape[0] == 361
    assert summary["package_gate_passed"] is True
    assert summary["precision"] == 1.0
    assert summary["cultivar_contamination_rate"] == 0.0
    assert set(scopes.loc[scopes["production_approved"], "trait_name"]) == {
        "floral_form",
        "inflorescence_display",
    }
    for record in manifest["files"]:
        path = root / record["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"]
