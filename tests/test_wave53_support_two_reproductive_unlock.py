from island_v2.wave53_incremental_all_evidence import EXPECTED_NEW_RULES
from island_v2.wave53_reproductive_checkpoint import (
    BASELINE_WAVE52_RUN_ID,
    DIRECT_ROWS,
    EXTERNAL_ROWS,
    IDENTITY_ROWS,
    REJECTED_ROWS,
    REVIEW_ROWS,
)


def test_wave53_contract_constants() -> None:
    assert BASELINE_WAVE52_RUN_ID == 33_470_056_127
    assert (DIRECT_ROWS, EXTERNAL_ROWS, REVIEW_ROWS, IDENTITY_ROWS, REJECTED_ROWS) == (
        1,
        1,
        2,
        2,
        2,
    )
    assert EXPECTED_NEW_RULES == frozenset(
        {
            ("Gomphrena", "reproductive_assurance", "self_incompatibility"),
            ("Marcgravia", "reproductive_assurance", "self_incompatibility"),
        }
    )
