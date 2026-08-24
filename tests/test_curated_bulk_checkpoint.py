from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import (
    validate_individually_reviewed_evidence,
)

ROOT = Path("data/v2/staging/traits/open_web_pilot")


def test_curated_bulk_checkpoint_is_reviewed_and_conflicts_are_retained() -> None:
    evidence = pd.read_csv(
        ROOT / "curated_bulk_direct_evidence_20260811.csv", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        ROOT / "curated_bulk_direct_manual_audit_20260811.csv", dtype=str
    ).fillna("")
    conflicts = pd.read_csv(
        ROOT / "curated_bulk_source_conflicts_20260811.csv", dtype=str
    ).fillna("")
    manifest = json.loads(
        (ROOT / "curated_bulk_checkpoint_manifest_20260811.json").read_text(
            encoding="utf-8"
        )
    )
    master = pd.read_csv(
        "data/v2/staging/gbif/collected/island_taxa.csv", dtype=str
    ).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence,
        audit,
        master,
        Path("config/trait_ontology.yml"),
    )
    assert len(accepted) == len(evidence) == len(audit) == 295
    assert evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0] == 253
    assert set(evidence["evidence_quality"]) == {"medium"}
    assert set(evidence["evidence_scope"]) == {"species_direct"}
    assert set(evidence["name_match_method"]) == {"accepted_name_exact"}
    assert set(evidence["inference_rule"]) == {""}
    assert set(conflicts["classification"].value_counts().to_dict()) == {
        "independent_source_agreement",
        "unresolved_direct_conflict",
    }
    assert int((conflicts["classification"] == "independent_source_agreement").sum()) == 1
    assert int((conflicts["classification"] == "unresolved_direct_conflict").sum()) == 41
    assert manifest["evidence_rows"] == 295
    assert manifest["unique_species_trait"] == 253
    assert manifest["selected_rows_by_source"] == {
        "eol_traitbank_flower_colour": 76,
        "meyer_2026": 117,
        "razanajatovo_2016": 102,
    }


def test_curated_bulk_checkpoint_excludes_completed_and_non_exact_records() -> None:
    selection = pd.read_csv(
        ROOT / "curated_bulk_selection_audit_20260811.csv.gz", dtype=str
    ).fillna("")
    selected = selection.loc[selection["selected"].eq("true")]
    assert len(selected) == 295
    assert set(selected["selection_reason"]) == {"selected_exact_direct_unresolved"}
    assert set(selected["name_match_method"]) == {"accepted_name_exact"}
    assert "already_resolved_in_baseline" in set(selection["selection_reason"])
    assert "not_exact_accepted_name" in set(selection["selection_reason"])
    assert "outside_checkpoint_trait_scope" in set(selection["selection_reason"])
