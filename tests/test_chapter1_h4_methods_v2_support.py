from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import yaml


SCRIPT = Path("scripts/audit_chapter1_h4_methods_v2_support.py")
SPEC = importlib.util.spec_from_file_location("h4_support_audit", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

CONTRACT = yaml.safe_load(
    Path("config/chapter1_h4_prospective_temporal_replication_v1.yml").read_text(
        encoding="utf-8"
    )
)


def _traits(n: int = 40) -> pd.DataFrame:
    rows = []
    for i in range(n):
        rows.append(
            {
                "accepted_species": f"Species {i:03d}",
                "autonomous_selfing": i % 2,
                "generalized_form": i % 2,
                "actinomorphic_symmetry": i % 2,
                "shallow_open_tube": i % 2,
            }
        )
    return pd.DataFrame(rows)


def test_status_tiers_do_not_mix_conflicts_into_automatic() -> None:
    methods = pd.DataFrame(
        {
            "doi": ["10/a", "10/b", "10/c"],
            "pmcid": ["", "", ""],
            "target_species": ["Species 000", "Species 001", "Species 002"],
            "candidate_status": [
                "design_candidate_auto",
                "design_candidate_field_unresolved",
                "design_candidate_conflict_unresolved",
            ],
        }
    )
    traits = _traits(3)
    automatic = MODULE.expand_target_species(
        methods,
        traits,
        MODULE.STATUS_TIERS["automatic_only"],
    )
    assert automatic["accepted_species"].tolist() == ["Species 000"]

    reviewable = MODULE.expand_target_species(
        methods,
        traits,
        MODULE.STATUS_TIERS["auto_plus_field_unresolved"],
    )
    assert set(reviewable["accepted_species"]) == {"Species 000", "Species 001"}


def test_support_audit_never_changes_frozen_thresholds() -> None:
    methods_rows = []
    for i in range(40):
        methods_rows.append(
            {
                "doi": f"10/{i:03d}",
                "pmcid": "",
                "target_species": f"Species {i:03d}",
                "candidate_status": "design_candidate_auto",
            }
        )
    result = MODULE.support_audit(
        pd.DataFrame(methods_rows),
        _traits(40),
        CONTRACT,
    )
    assert result["thresholds_changed"] is False
    assert result["outcomes_read"] is False
    assert result["tiers"]["automatic_only"]["H4a_reproductive_assurance"]["evaluable"]
    assert result["tiers"]["automatic_only"]["H4b_accessibility_generalization"]["evaluable"]
