import hashlib
import json
from pathlib import Path

import pandas as pd

from island_v2.high_leverage_direct_checkpoint import (
    pfaf_rows,
    resolved_direct_pairs,
)
from island_v2.open_web_common import reviewed_source_package_evidence
from island_v2.open_web_finalize import validate_individually_reviewed_evidence


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def test_only_resolved_direct_cells_suppress_reacquisition() -> None:
    ledger = pd.DataFrame(
        [
            {
                "accepted_species": "Testus resolved",
                "trait_name": "self_incompatibility",
                "resolution_status": "resolved",
            },
            {
                "accepted_species": "Testus conflict",
                "trait_name": "self_incompatibility",
                "resolution_status": "unresolved",
            },
        ]
    )
    assert resolved_direct_pairs(ledger) == {
        ("Testus resolved", "self_incompatibility")
    }


def test_pfaf_maps_self_fertile_only_to_sc() -> None:
    candidates = pd.DataFrame(
        [
            {
                "accepted_species": "Testus exactus",
                "source_url": "https://pfaf.org/user/Plant.aspx?LatinName=Testus+exactus",
                "matched_page_name": "Testus exactus",
                "identity_exact": "True",
                "supporting_excerpt": "The plant is self-fertile",
                "normalized_value": "SC",
                "trait_name": "self_incompatibility",
                "status": "accepted_self_fertile_statement",
                "cultivar_qualified": "False",
                "nonsexual_qualified": "False",
                "dioecious_conflict": "False",
                "page_sha256": "a" * 64,
                "retrieved_at_utc": "2026-08-11T04:00:00+00:00",
            }
        ]
    )
    rows = pfaf_rows(
        candidates, master_names={"Testus exactus"}, completed_pairs=set()
    )
    assert len(rows) == 1
    assert rows[0]["trait_name"] == "self_incompatibility"
    assert rows[0]["normalized_value"] == "SC"
    assert rows[0]["evidence_quality"] == "medium"
    assert not rows[0]["inference_rule"]


def test_primary_rows_keep_three_reproductive_concepts_separate(tmp_path: Path) -> None:
    # The production hashes are tested by the committed checkpoint.  This unit
    # test exercises the filtering contract without fabricating source bytes.
    source = Path("src/island_v2/high_leverage_direct_checkpoint.py").read_text(
        encoding="utf-8"
    )
    assert 'trait="mating_system"' in source
    assert 'trait="self_incompatibility"' in source
    assert 'trait="autonomous_selfing_capacity"' in source
    assert "outcrossing_rate_to_mating_system_only" in source
    assert "self_fertile_to_self_incompatibility_only" in source


def test_committed_checkpoint_passes_individual_review_gate() -> None:
    root = Path(
        "data/v2/staging/traits/open_web_pilot/"
        "high_leverage_direct_checkpoint_20260811"
    )
    evidence = pd.read_csv(
        root / "high_leverage_direct_evidence_20260811.csv", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        root / "high_leverage_direct_manual_audit_20260811.csv", dtype=str
    ).fillna("")
    master = pd.DataFrame({"accepted_species": evidence["accepted_species"].unique()})

    accepted = validate_individually_reviewed_evidence(
        evidence,
        audit,
        master,
        Path("config/trait_ontology.yml"),
    )

    assert len(evidence) == 81
    assert len(accepted) == 81
    assert evidence["accepted_species"].nunique() == 81
    assert evidence["candidate_id"].is_unique
    assert evidence["source_excerpt"].ne("").all()
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert set(evidence["evidence_quality"]) == {"high", "medium"}
    assert not set(evidence["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )
    assert audit["decision"].eq("accept").all()
    assert audit["cultivar_contamination"].eq("false").all()


def test_combined_source_package_is_complete_and_hash_pinned() -> None:
    root = Path(
        "data/v2/staging/traits/open_web_pilot/"
        "combined_next_source_package_20260811"
    )
    evidence_path = root / "combined_next_evidence_20260811.csv.gz"
    audit_path = root / "combined_next_reviewed_audit_600_20260811.csv"
    manifest = json.loads(
        (root / "combined_next_source_package_manifest_20260811.json").read_text(
            encoding="utf-8"
        )
    )
    evidence = pd.read_csv(evidence_path, dtype=str).fillna("")
    audit = pd.read_csv(audit_path, dtype=str).fillna("")

    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)

    assert len(evidence) == len(selected) == 1062
    assert len(audit) == 600
    assert summary["package_gate_passed"]
    assert summary["precision"] == 1.0
    assert summary["cultivar_contamination_rate"] == 0.0
    assert scopes["production_approved"].all()
    assert evidence["source_record_id"].is_unique
    assert audit["candidate_id"].is_unique
    assert _sha256(evidence_path) == manifest["files"][evidence_path.name]["sha256"]
    assert _sha256(audit_path) == manifest["files"][audit_path.name]["sha256"]
