import json
from pathlib import Path

import pandas as pd

from island_v2.florml_size_recovery_checkpoint import strict_size_candidate
from island_v2.open_web_common import reviewed_source_package_evidence


def _base(quote: str) -> dict[str, str]:
    return {
        "accepted_species": "Testus exactus",
        "trait_name": "flower_size_class",
        "source_provider": "flora_malesiana",
        "source_excerpt": quote,
        "evidence_scope": "species_direct",
        "name_match_method": "wfo_june_2026_exact_accepted_xml_species_family_consistent",
    }


def test_size_gate_separates_whole_flower_from_other_organs() -> None:
    assert strict_size_candidate(_base("Flowers 5-8 mm long."))[0]
    assert strict_size_candidate(_base("Corolla funnel-shaped, 2 cm long."))[0]
    assert not strict_size_candidate(
        _base("Flowers sessile, subtended by 2 mm long bracts.")
    )[0]
    assert not strict_size_candidate(
        _base("Pedicels of the central flowers 0-0.3 mm long.")
    )[0]
    assert not strict_size_candidate(
        _base("Flowers in a thyrse up to 18 mm long.")
    )[0]
    assert not strict_size_candidate(_base("Male flowers 3 mm long."))[0]


def test_committed_size_checkpoint_passes_trait_review_gate() -> None:
    root = Path(
        "data/v2/staging/traits/open_web_pilot/"
        "florml_size_recovery_checkpoint_20260812"
    )
    evidence = pd.read_csv(
        root / "florml_size_evidence_20260812.csv.gz", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        root / "florml_size_reviewed_audit_200_20260812.csv", dtype=str
    ).fillna("")
    manifest = json.loads(
        (root / "florml_size_manifest_20260812.json").read_text(encoding="utf-8")
    )

    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)

    assert len(evidence) == 321
    assert len(selected) == 320
    assert len(audit) == 200
    assert summary["precision"] == 0.995
    assert summary["cultivar_contamination_rate"] == 0.0
    assert summary["explicitly_rejected_candidate_ids"] == 1
    assert scopes["production_approved"].all()
    assert manifest["guardrails"]["whole_flower_measurement_only"]
