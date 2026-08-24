from pathlib import Path

import pandas as pd

from analysis.prepare_wfo_high_yield_source_package_20260810 import (
    ambiguous_display_quote,
    safe_colour_row,
    safe_display_row,
)
from island_v2.open_web_common import reviewed_source_package_evidence

PACKAGE = Path(
    "data/v2/staging/traits/direct_llm_pilot/"
    "20260810_wfo_high_yield_source_acquisition"
)


def _colour(value: str, quote: str) -> pd.Series:
    return pd.Series({"normalized_value": value, "source_excerpt": quote})


def test_colour_gate_is_complete_and_rejects_nonfloral_context() -> None:
    assert safe_colour_row(_colour("white|red_pink", "Corolla white or pink."))
    assert not safe_colour_row(
        _colour(
            "green_brown_inconspicuous",
            "Male flowers densely brown pubescent; buds 6 mm long.",
        )
    )
    assert not safe_colour_row(_colour("yellow_orange", "Fruit: perianth yellow."))
    assert not safe_colour_row(_colour("white", "Corolla white or pink."))


def test_display_gate_rejects_partial_units_and_sex_specific_counts() -> None:
    assert ambiguous_display_quote(
        "Male inflorescence 30 cm; partial peduncles 2-flowered."
    )
    assert ambiguous_display_quote(
        "In male many-flowered; in female few-flowered."
    )
    assert not ambiguous_display_quote("Inflorescence a terminal compound umbel.")


def test_display_completeness_gate_preserves_every_explicit_state() -> None:
    assert safe_display_row(
        pd.Series(
            {
                "normalized_value": "composite_display|umbel_corymb",
                "source_excerpt": "Capitula arranged in terminal corymbs.",
            }
        )
    )
    assert not safe_display_row(
        pd.Series(
            {
                "normalized_value": "umbel_corymb",
                "source_excerpt": "Capitula arranged in terminal corymbs.",
            }
        )
    )
    assert safe_display_row(
        pd.Series(
            {
                "normalized_value": "raceme_spike_panicle",
                "source_excerpt": "Racemes solitary, including the peduncle.",
            }
        )
    )
    assert safe_display_row(
        pd.Series(
            {
                "normalized_value": "few_flowered|raceme_spike_panicle",
                "source_excerpt": "Panicles narrow; branches 2--4-flowered.",
            }
        )
    )
    assert not safe_display_row(
        pd.Series(
            {
                "normalized_value": "umbel_corymb",
                "source_excerpt": "Racemes corymbiform, flowers rarely solitary.",
            }
        )
    )


def test_committed_package_has_exact_audit_and_trait_specific_gate() -> None:
    evidence = pd.read_csv(
        PACKAGE / "wfo_high_yield_incremental_evidence.csv.gz", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        PACKAGE / "wfo_high_yield_manual_audit_200.csv", dtype=str
    ).fillna("")

    assert len(evidence) == 6395
    assert evidence["source_record_id"].is_unique
    assert set(evidence["trait_name"]) == {
        "floral_form",
        "flower_primary_color",
        "flower_size_class",
        "inflorescence_display",
    }
    assert len(audit) == 200
    assert audit["candidate_id"].is_unique
    assert set(audit["candidate_id"]).issubset(evidence["source_record_id"])

    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)
    approved = set(scopes.loc[scopes["production_approved"], "trait_name"])
    assert approved == {
        "floral_form",
        "flower_primary_color",
        "inflorescence_display",
    }
    assert summary["package_gate_passed"] is True
    assert summary["reviewed"] == 200
    assert summary["accepted_correct"] == 189
    assert summary["precision"] == 0.945
    assert len(selected) == 5467
    assert "flower_size_class" not in set(selected["trait_name"])
