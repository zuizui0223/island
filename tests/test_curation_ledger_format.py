"""Every hand-curated ledger must parse before a heavy job tries to read it.

These files are written by hand during manual review, so the failure mode is a
free-text field containing a comma or a quote and not being quoted for it. That
is not a hypothetical: the Macrolobium trinitense excerpt -- `The tree had
"inflorescences of delicate, crenellated white flowers" ...` -- was written
unquoted, so pandas split it into an extra field and silently shifted every
column after it. `evidence_scope` came back holding the tail of the sentence and
`name_match_method` came back holding `species_direct`.

It reached CI because nothing parsed the file until the queue-build job
downloaded it, and it stayed red across four heads while more rows were appended
on top. A row that cannot be parsed is worse than a row that is missing, because
the shift is silent: the columns are all populated, just with the wrong values.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

CURATION_DIR = Path("data/v2/curation")
LEDGERS = sorted(CURATION_DIR.glob("*.csv")) if CURATION_DIR.exists() else []

# Columns whose values come from a closed set. A shifted row usually shows up
# here first, holding prose where a token belongs.
CLOSED_VOCABULARY = {
    "evidence_scope": {"species_direct", "synonym_direct", "genus_inferred"},
    "quality": {"high", "medium", "low"},
}


@pytest.mark.skipif(not LEDGERS, reason="no curation ledgers present")
@pytest.mark.parametrize("path", LEDGERS, ids=lambda p: p.name)
def test_every_row_has_the_declared_field_count(path: Path):
    rows = list(csv.reader(path.open(encoding="utf-8")))
    assert rows, f"{path.name} is empty"
    header = rows[0]

    malformed = [
        (number, len(row), row[0] if row else "")
        for number, row in enumerate(rows[1:], start=2)
        if len(row) != len(header)
    ]
    assert not malformed, (
        f"{path.name}: rows do not match the {len(header)}-column header "
        f"{malformed}. A free-text field containing a comma or a quote must be "
        f"quoted, and an inner quote doubled."
    )


@pytest.mark.skipif(not LEDGERS, reason="no curation ledgers present")
@pytest.mark.parametrize("path", LEDGERS, ids=lambda p: p.name)
def test_closed_vocabulary_columns_hold_tokens_not_prose(path: Path):
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    if not rows:
        pytest.skip(f"{path.name} has no data rows")

    for column, permitted in CLOSED_VOCABULARY.items():
        if column not in rows[0]:
            continue
        unexpected = sorted(
            {
                str(row[column]).strip()
                for row in rows
                if str(row.get(column) or "").strip() not in permitted
            }
        )
        assert not unexpected, (
            f"{path.name}: column {column} holds {unexpected}, which is outside "
            f"{sorted(permitted)}. This is the signature of a shifted row: the "
            f"columns are populated, but one field to the left of where they belong."
        )
