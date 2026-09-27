import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SI = ROOT / "submission/chapter1_current/SUPPLEMENTARY_INFORMATION.md"


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def test_supplement_has_complete_s1_s7_structure_and_live_data_paths() -> None:
    text = SI.read_text(encoding="utf-8")
    for idx in range(1, 8):
        assert f"# Appendix S{idx}." in text

    for idx in range(1, 7):
        assert f"Table S{idx}." in text

    refs = sorted(set(re.findall(r"`(results/geography_20260924/[^`]+)`", text)))
    assert refs
    for ref in refs:
        path = ROOT / ref
        assert path.exists(), ref


def test_supplement_h1_values_match_corrected_tables() -> None:
    text = SI.read_text(encoding="utf-8")
    rows = _read_csv(ROOT / "results/geography_20260924/all/beta_binomial_within_omnibus.csv")
    primary = {
        row["context"]: row
        for row in rows
        if row["stratum"] == "all_observed"
    }
    assert len(primary) == 4
    assert all(row["vector_supported"] == "True" for row in primary.values())

    expected = {
        "northern_midlatitude": "3.216 × 10^-10",
        "northern_high_latitude": "2.433 × 10^-5",
        "tropical": "3.498 × 10^-7",
        "southern_extratropical": "2.504 × 10^-18",
    }
    for value in expected.values():
        assert value in text


def test_supplement_h2_accessibility_table_matches_corrected_models() -> None:
    text = SI.read_text(encoding="utf-8")
    all_rows = _read_csv(ROOT / "results/geography_20260924/all/h2_decomposition_models.csv")
    direct_rows = _read_csv(ROOT / "results/geography_20260924/direct/h2_decomposition_models.csv")

    for rows in (all_rows, direct_rows):
        access = [row for row in rows if row["response"] == "generalized_accessible"]
        assert len(access) == 4
        assert all(float(row["distance_estimate"]) > 0 for row in access)

    assert "0.000847" in text
    assert "0.01616" in text
    assert "0.1196" in text


def test_supplement_h3_h4_values_match_corrected_sources() -> None:
    text = SI.read_text(encoding="utf-8")

    h3 = json.loads(
        (ROOT / "results/geography_20260924/h3_original_corrected_comparison.json")
        .read_text(encoding="utf-8")
    )["corrected"]["global_gradient"]
    assert round(float(h3["distance_slope"]), 5) == 0.09191
    assert round(float(h3["distance_slope_se"]), 5) == 0.03806
    assert "0.09191" in text
    assert "0.01575" in text

    h4 = _read_csv(ROOT / "results/geography_20260924/h4_exact_corrected.csv")
    primary = {row["family"]: row for row in h4 if row["analysis"] == "primary"}
    assert round(float(primary["reproductive_assurance"]["estimate"]), 5) == -0.29830
    assert round(float(primary["accessibility_generalization"]["estimate"]), 5) == -0.29566
    assert "-0.29830" in text
    assert "-0.29566" in text


def test_supplement_preserves_inferential_boundaries() -> None:
    text = SI.read_text(encoding="utf-8")
    assert "not mediation" in text
    assert "post-hoc functional triangulation" in text
    assert "does not distinguish species sorting" in text
    assert "before outcome unblinding" in text
