from pathlib import Path

import pandas as pd
import pytest

from analysis.recover_wfo_morphology_measurements_20260808 import (
    classify_measurement_subject,
)
from island_v2.open_web_common import reviewed_source_package_evidence


def _row(quote: str, raw_value: str) -> pd.Series:
    return pd.Series(
        {
            "exact_supporting_quote": quote,
            "raw_value": raw_value,
        }
    )


@pytest.mark.parametrize(
    ("quote", "raw_value", "expected"),
    [
        (
            "corolla white, the tube 10-15 mm long, the lobes spreading;",
            "10-15 mm",
            "tube_depth_class",
        ),
        (
            "Flower: length of the tube of the corolla 4-5 mm.",
            "4-5 mm",
            "tube_depth_class",
        ),
        (
            "Corolla white-pink with a 1 mm long tube;",
            "1 mm",
            "tube_depth_class",
        ),
        (
            "Corolla 3.5-4 cm long, tube 2-3 cm long;",
            "3.5-4 cm",
            "flower_size_class",
        ),
        (
            "Flower: length of the petal 48 mm or more;",
            "48 mm",
            "",
        ),
        (
            "Flowers with the calycine tube 0.35-0.50 mm long;",
            "0.35-0.50 mm",
            "",
        ),
    ],
)
def test_measurement_is_assigned_to_its_trait_specific_subject(
    quote: str,
    raw_value: str,
    expected: str,
) -> None:
    trait, _ = classify_measurement_subject(_row(quote, raw_value))
    assert trait == expected


def test_committed_wfo_morphology_package_is_a_full_census_passing_gate() -> None:
    root = Path("data/v2/staging/traits/direct_llm_pilot/20260808_wfo_morphology_recovery")
    evidence = pd.read_csv(root / "wfo_morphology_recovered_evidence.csv.gz", dtype=str).fillna("")
    audit = pd.read_csv(root / "wfo_morphology_full_census_audit.csv", dtype=str).fillna("")

    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)

    assert summary["reviewed"] == 244
    assert summary["accepted_correct"] == 242
    assert summary["precision"] == pytest.approx(242 / 244)
    assert summary["cultivar_contamination_rate"] == 0.0
    assert summary["package_gate_passed"]
    assert set(scopes.loc[scopes["production_approved"], "trait_name"]) == {
        "flower_size_class",
        "tube_depth_class",
    }
    assert len(selected) == 242
    assert set(selected["source_provider"]) == {
        "flora_north_america",
        "nybg_brittonia",
        "nybg_memoirs",
    }
    assert selected["source_record_id"].is_unique
