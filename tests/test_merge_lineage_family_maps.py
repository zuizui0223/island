from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


MODULE_PATH = Path(__file__).parents[1] / "analysis" / "merge_lineage_family_maps.py"
spec = importlib.util.spec_from_file_location("merge_lineage_family_maps", MODULE_PATH)
assert spec is not None and spec.loader is not None
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def _write(path: Path, rows: list[dict[str, str]]) -> Path:
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def test_merge_deduplicates_same_assignment(tmp_path: Path) -> None:
    first = _write(
        tmp_path / "first.csv",
        [{"source_lineage": "lineage:a", "source_family": "dataset:a"}],
    )
    second = _write(
        tmp_path / "second.csv",
        [
            {"source_lineage": "lineage:a", "source_family": "dataset:a"},
            {"source_lineage": "lineage:b", "source_family": "database:b"},
        ],
    )

    summary = module.merge_maps(
        (first, second), tmp_path / "merged.csv", tmp_path / "summary.json"
    )
    assert summary["merged_rows"] == 2
    assert summary["duplicate_same_family_rows_removed"] == 1
    assert summary["conflicting_lineages"] == 0
    assert summary["scientific_database_modified"] is False
    assert summary["rights_granted"] is False


def test_merge_rejects_conflicting_assignment(tmp_path: Path) -> None:
    first = _write(
        tmp_path / "first.csv",
        [{"source_lineage": "lineage:a", "source_family": "dataset:a"}],
    )
    second = _write(
        tmp_path / "second.csv",
        [{"source_lineage": "lineage:a", "source_family": "database:b"}],
    )

    try:
        module.merge_maps((first, second), tmp_path / "merged.csv")
    except ValueError as exc:
        assert "conflicting lineage-family assignments" in str(exc)
    else:
        raise AssertionError("conflicting assignment should fail closed")
