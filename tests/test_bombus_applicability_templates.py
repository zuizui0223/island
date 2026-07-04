import csv
from pathlib import Path

from island_v2.bombus_applicability import (
    REQUIRED_ASSIGNMENT_EVIDENCE_COLUMNS,
    REQUIRED_EVIDENCE_COLUMNS,
)


def _header(path: str) -> set[str]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return set(next(csv.reader(handle)))


def test_source_region_evidence_template_covers_all_required_columns():
    header = _header("data/v2/templates/source_region_evidence_template.csv")

    assert REQUIRED_EVIDENCE_COLUMNS.issubset(header)


def test_island_assignment_evidence_template_covers_all_required_columns():
    header = _header("data/v2/templates/island_assignment_evidence_template.csv")

    assert REQUIRED_ASSIGNMENT_EVIDENCE_COLUMNS.issubset(header)
