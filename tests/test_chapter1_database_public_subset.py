from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from island_v2.chapter1_database_public_subset import build_public_subset


def _write_inputs(tmp_path: Path) -> tuple[Path, Path]:
    species_axis = pd.DataFrame(
        {
            "accepted_species": ["Alpha one", "Beta two", "Gamma three"],
            "axis": ["flower_colour", "flower_colour", "reproductive_assurance"],
            "quality": ["high", "medium", "unresolved"],
            "trait_composition": ["red", "white", ""],
            "source_lineages": ["dataset:dryad", "domain:blocked.example", ""],
        }
    )
    rights = pd.DataFrame(
        {
            "accepted_species": ["Alpha one", "Beta two"],
            "axis": ["flower_colour", "flower_colour"],
            "quality": ["high", "medium"],
            "release_status": ["redistributable", "rights_review_required"],
        }
    )
    species_path = tmp_path / "species_axis_coverage.csv.gz"
    rights_path = tmp_path / "CELL_RELEASE_RIGHTS.csv.gz"
    species_axis.to_csv(species_path, index=False, compression="gzip")
    rights.to_csv(rights_path, index=False, compression="gzip")
    return species_path, rights_path


def test_public_subset_keeps_only_redistributable_resolved_cells(tmp_path: Path) -> None:
    species_path, rights_path = _write_inputs(tmp_path)
    subset, counts = build_public_subset(species_path, rights_path)
    assert subset[["accepted_species", "axis"]].to_dict("records") == [
        {"accepted_species": "Alpha one", "axis": "flower_colour"}
    ]
    assert counts == {
        "resolved_cells": 2,
        "redistributable_cells": 1,
        "blocked_or_review_cells": 1,
    }


def test_public_subset_fails_if_rights_ledger_omits_resolved_cell(tmp_path: Path) -> None:
    species_path, rights_path = _write_inputs(tmp_path)
    rights = pd.read_csv(rights_path).iloc[:1]
    rights.to_csv(rights_path, index=False, compression="gzip")
    with pytest.raises(ValueError, match="missing 1 resolved"):
        build_public_subset(species_path, rights_path)


def test_public_subset_fails_on_quality_mismatch(tmp_path: Path) -> None:
    species_path, rights_path = _write_inputs(tmp_path)
    rights = pd.read_csv(rights_path)
    rights.loc[0, "quality"] = "low"
    rights.to_csv(rights_path, index=False, compression="gzip")
    with pytest.raises(ValueError, match="quality mismatch"):
        build_public_subset(species_path, rights_path)
