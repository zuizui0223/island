from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


MODULE_PATH = Path(__file__).parents[1] / "analysis" / "evaluate_chapter1_database_cell_rights.py"
spec = importlib.util.spec_from_file_location("evaluate_chapter1_database_cell_rights", MODULE_PATH)
assert spec is not None and spec.loader is not None
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def _write(path: Path, rows: list[dict[str, object]]) -> Path:
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def test_cell_gate_replaces_derived_lineage_with_upstream_rights(tmp_path: Path) -> None:
    species_axis = _write(
        tmp_path / "species_axis.csv",
        [
            {
                "accepted_species": "Alpha direct",
                "axis": "flower_colour",
                "quality": "high",
                "source_lineages": "dataset:dryad.example:row1",
            },
            {
                "accepted_species": "Beta derived open",
                "axis": "reproductive_assurance",
                "quality": "low",
                "source_lineages": 'validated-low:Beta:self_incompatibility:["SC"]',
            },
            {
                "accepted_species": "Gamma derived blocked",
                "axis": "reproductive_assurance",
                "quality": "low",
                "source_lineages": 'validated-low:Gamma:self_incompatibility:["SC"]',
            },
            {
                "accepted_species": "Delta semantic",
                "axis": "flower_colour",
                "quality": "low",
                "source_lineages": 'validated-low:Delta:flower_primary_color:["red_pink"]',
            },
        ],
    )
    recovery = _write(
        tmp_path / "recovery.csv",
        [
            {
                "accepted_species": "Beta derived open",
                "axis": "reproductive_assurance",
                "derived_lineages": 'validated-low:Beta:self_incompatibility:["SC"]',
                "support_source_families": "dataset:dryad",
                "all_rules_source_provenance_recovered": "true",
                "all_rules_provenance_recovered": "true",
            },
            {
                "accepted_species": "Gamma derived blocked",
                "axis": "reproductive_assurance",
                "derived_lineages": 'validated-low:Gamma:self_incompatibility:["SC"]',
                "support_source_families": "dataset:dryad|database:pladias",
                "all_rules_source_provenance_recovered": "true",
                "all_rules_provenance_recovered": "true",
            },
            {
                "accepted_species": "Delta semantic",
                "axis": "flower_colour",
                "derived_lineages": 'validated-low:Delta:flower_primary_color:["red_pink"]',
                "support_source_families": "dataset:dryad",
                "all_rules_source_provenance_recovered": "true",
                "all_rules_provenance_recovered": "false",
            },
        ],
    )
    policy = tmp_path / "policy.yml"
    policy.write_text(
        """schema_version: 2
default_status: review_required
default_license: null
rules:
  - pattern: '^dataset:dryad$'
    status: redistributable
    license: CC0-1.0
  - pattern: '^database:pladias$'
    status: review_required
    license: null
""",
        encoding="utf-8",
    )

    summary = module.evaluate(species_axis, recovery, policy, tmp_path / "out")

    assert summary["resolved_cells"] == 4
    assert summary["redistributable_cells"] == 2
    assert summary["blocked_cells"] == 2
    assert summary["validated_low_cells"] == 3
    assert summary["validated_low_cells_with_complete_upstream_source_provenance"] == 3
    assert summary["semantic_mismatch_cells"] == 1
    assert summary["distinct_blocker_source_families"] == 1
    assert not summary["release_ready_for_public_zenodo"]
    assert not summary["synthetic_validated_low_family_used_as_rights_decision"]

    cells = pd.read_csv(tmp_path / "out" / "CELL_RELEASE_RIGHTS.csv.gz", dtype=str).fillna("")
    status = dict(zip(cells["accepted_species"], cells["release_status"], strict=True))
    assert status == {
        "Alpha direct": "redistributable",
        "Beta derived open": "redistributable",
        "Gamma derived blocked": "rights_review_required",
        "Delta semantic": "semantic_review_required",
    }

    gamma = cells.loc[cells["accepted_species"].eq("Gamma derived blocked")].iloc[0]
    assert gamma["blocked_source_families"] == "database:pladias"
    assert "derived:validated_low_without_direct_source_lineage" not in gamma["required_source_families"]
