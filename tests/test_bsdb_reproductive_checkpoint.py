from __future__ import annotations

import pandas as pd
import pytest

from island_v2.bsdb_reproductive_checkpoint import (
    apply_reviewed_audit,
    build_checkpoint,
)


def epithet(index: int) -> str:
    return f"species{chr(97 + index // 26)}{chr(97 + index % 26)}"


def source_rows() -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for index in range(205):
        rows.append(
            {
                "Infrasp": "",
                "BreedingSys": "SC" if index < 102 else "SI",
                "ISI_value": "0.1" if index < 102 else "0.95",
                "ISItype": "ISIfruits",
                "bs.Source": f"study{index % 17}",
                "tnrs_status": "Synonym" if index == 0 else "Accepted",
                "tnrs_family": "Testaceae",
                "tnrs_Sciname": f"Testus {epithet(index)}",
                "tnrs_infrasp": "",
            }
        )
    rows.extend(
        [
            {
                **rows[0],
                "tnrs_Sciname": "Testus outside",
            },
            {
                **rows[1],
                "tnrs_Sciname": f"Testus {epithet(1)}",
                "tnrs_infrasp": "minor",
            },
            {
                **rows[2],
                "tnrs_Sciname": f"Testus {epithet(2)}",
                "tnrs_family": "Wrongaceae",
            },
        ]
    )
    return pd.DataFrame(rows)


def master() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "accepted_species": [f"Testus {epithet(index)}" for index in range(205)],
            "family": ["Testaceae"] * 205,
        }
    )


def test_checkpoint_is_exact_species_direct_and_skips_completed_pairs() -> None:
    baseline = pd.DataFrame(
        {
            "accepted_species": [f"Testus {epithet(3)}"],
            "trait_name": ["self_incompatibility"],
        }
    )
    evidence, audit, selection = build_checkpoint(source_rows(), baseline, master())

    assert len(evidence) == 204
    assert len(audit) == 200
    assert set(evidence["trait_name"]) == {"self_incompatibility"}
    assert set(evidence["axis"]) == {"reproductive_assurance"}
    assert set(evidence["quality"]) == {"high"}
    assert f"Testus {epithet(3)}" not in set(evidence["accepted_species"])
    assert selection["selection_reason"].value_counts().to_dict() == {
        "selected": 204,
        "already_resolved_direct_pair": 1,
        "not_in_strict_island_universe": 1,
        "infraspecific_record": 1,
        "family_conflict": 1,
    }


def test_checkpoint_preserves_synonym_scope_and_original_study_lineage() -> None:
    baseline = pd.DataFrame(columns=["accepted_species", "trait_name"])
    evidence, _, _ = build_checkpoint(source_rows(), baseline, master())
    synonym = evidence.loc[
        evidence["accepted_species"].eq(f"Testus {epithet(0)}")
    ].iloc[0]
    same_study = evidence.loc[evidence["source_lineage"].eq(synonym["source_lineage"])]

    assert synonym["evidence_scope"] == "synonym_direct"
    assert synonym["name_match_method"] == "exact_synonym_to_accepted"
    assert len(same_study) > 1
    assert same_study["source_record_id"].nunique() == len(same_study)


def test_review_attachment_rejects_evidence_edits() -> None:
    baseline = pd.DataFrame(columns=["accepted_species", "trait_name"])
    _, audit, _ = build_checkpoint(source_rows(), baseline, master())
    reviewed = audit.copy()
    reviewed["accepted_correct"] = "true"
    reviewed["cultivar_status"] = "wild_species_database_record"
    reviewed["reviewer"] = "reviewer"
    reviewed["reviewed_at_utc"] = "2026-08-11T00:00:00Z"
    reviewed["audit_reason"] = "accepted after source-row review"
    reviewed["audit_hash"] = "hash"

    attached = apply_reviewed_audit(audit, reviewed)
    assert attached["accepted_correct"].eq("true").all()

    reviewed.loc[0, "normalized_value"] = "SI"
    with pytest.raises(ValueError, match="changed immutable column"):
        apply_reviewed_audit(audit, reviewed)
