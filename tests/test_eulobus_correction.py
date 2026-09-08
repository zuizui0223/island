import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).parents[1]
SPEC = importlib.util.spec_from_file_location(
    "integration", ROOT / "scripts/integrate_reviewed_restart.py"
)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def sample():
    receipt = json.loads(
        (ROOT / "data/v2/staging/traits/restart_eulobus_correction_20260908.json").read_text(
            encoding="utf-8"
        )
    )
    coverage = pd.DataFrame(
        [
            {
                "accepted_species": "Eulobus californicus",
                "axis": "reproductive_assurance",
                "quality": "medium",
                "trait_names": "self_incompatibility",
                "trait_composition": 'self_incompatibility=["SI"]',
                "source_groups": "reviewed_restart_20260907",
                "source_lineages": "dataset:dryad.cc2fqz6hr:v4",
            }
        ]
    )
    direct = pd.DataFrame(
        [
            {
                "accepted_species": "Eulobus californicus",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "classification": "single_independent_lineage",
                "resolution_status": "resolved",
                "selected_quality": "medium",
                "quality": "medium",
                "state_set": '["SI"]',
                "state_sets": '[["SI"]]',
                "normalized_value": "SI",
                "source_lineages": "dataset:dryad.cc2fqz6hr:v4",
                "source_groups": "reviewed_restart_20260907",
            }
        ]
    )
    ontology = {"traits": {"self_incompatibility": {"allowed_values": ["SI", "SC", "mixed_or_variable"]}}}
    return coverage, direct, receipt, ontology


def test_eulobus_correction_changes_value_not_filled_cell_count():
    coverage, direct, receipt, ontology = sample()
    after, ledger = module.correct_eulobus(coverage, direct, receipt, ontology)
    assert after["quality"].tolist() == ["high"]
    assert after["trait_composition"].tolist() == ['self_incompatibility=["SC"]']
    assert ledger["normalized_value"].tolist() == ["SC"]
    assert ledger["quality"].tolist() == ["high"]
    assert after["quality"].ne("").sum() == coverage["quality"].ne("").sum()


def test_eulobus_correction_fails_if_meyer_precondition_changes():
    coverage, direct, receipt, ontology = sample()
    direct.loc[0, "normalized_value"] = "SC"
    with pytest.raises(ValueError, match="precondition"):
        module.correct_eulobus(coverage, direct, receipt, ontology)


def test_eulobus_correction_fails_if_receipt_is_modified():
    coverage, direct, receipt, ontology = sample()
    receipt["support_basis"] = "weaker claim"
    with pytest.raises(ValueError, match="Unapproved"):
        module.correct_eulobus(coverage, direct, receipt, ontology)
