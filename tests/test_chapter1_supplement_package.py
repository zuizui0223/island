import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "submission" / "chapter1_current"
SI = PACKAGE / "SUPPLEMENTARY_INFORMATION.md"
TABLE_DIR = PACKAGE / "supplement"


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def test_supplement_has_complete_s1_s7_structure() -> None:
    text = SI.read_text(encoding="utf-8")
    for idx in range(1, 8):
        assert f"# Appendix S{idx}." in text

    required_tables = (
        "Table S1.",
        "Table S2.",
        "Table S3.",
        "Table S5.",
        "Table S6a.",
        "Table S6b.",
    )
    for label in required_tables:
        assert label in text

    assert "Table S4" in text
    assert "raw-pattern" in text.lower()


def test_supplement_references_live_corrected_paths_and_generated_tables() -> None:
    text = SI.read_text(encoding="utf-8")
    refs = sorted(
        set(
            re.findall(
                r"((?:results/(?:geography_20260924|h1_final_traitwise_t_20261004|h1_final_directional_20261003|h2_complete_selfing_corrected_20261003|h3_offshore_gradient_20261003)|submission/chapter1_current/supplement)/[^\s,;]+)",
                text,
            )
        )
    )
    assert refs
    for ref in refs:
        clean = ref.rstrip(".:)")
        path = ROOT / clean
        assert path.exists(), clean

    for name in (
        "Table_S1_data_summary.csv",
        "Table_S2_H1_traitwise.csv",
        "Table_S3a_H2_all_finite_cluster.csv",
        "Table_S3b_H2_direct_finite_cluster.csv",
        "Table_S5_H3_pollen_limitation.csv",
        "Table_S6a_H4_scores.csv",
        "SUPPLEMENT_TABLES_MANIFEST.json",
    ):
        assert (TABLE_DIR / name).is_file(), name


def test_supplement_h1_values_match_final_traitwise_table() -> None:
    text = SI.read_text(encoding="utf-8")
    h1 = _read_csv(TABLE_DIR / "Table_S2_H1_traitwise.csv")
    assert len(h1) == 112
    assert all(row["optimizer_success"] == "True" for row in h1)

    primary_sc = [
        row
        for row in h1
        if row["evidence_scope"] == "all"
        and row["flora_scope"] == "broad"
        and row["outcome"] == "self_compatibility"
    ]
    assert len(primary_sc) == 4
    assert all(float(row["estimate"]) > 0 for row in primary_sc)
    assert all(float(row["p_two_sided"]) < 0.05 for row in primary_sc)

    assert "0.00569" in text
    assert "0.00671" in text
    assert "0.01812" in text
    assert "0.01123" in text
    assert "Earlier joint/vector H1 tables" in text


def test_supplement_h2_h3_h4_values_match_generated_tables() -> None:
    text = SI.read_text(encoding="utf-8")

    h2_all = _read_csv(TABLE_DIR / "Table_S3a_H2_all_finite_cluster.csv")
    h2_direct = _read_csv(TABLE_DIR / "Table_S3b_H2_direct_finite_cluster.csv")
    access = [
        row for row in h2_all + h2_direct if row["response"] == "generalized_accessible"
    ]
    assert len(access) == 8
    assert all(float(row["distance_estimate"]) > 0 for row in access)
    assert "0.01913" in text
    assert "0.1267" in text

    h3 = _read_csv(TABLE_DIR / "Table_S5_H3_pollen_limitation.csv")
    primary_h3 = next(row for row in h3 if row["analysis"] == "primary")
    assert round(float(primary_h3["estimate"]), 5) == 0.09191
    assert round(float(primary_h3["finite_publication_t_two_sided_p"]), 5) == 0.01594
    assert "0.09191" in text
    assert "0.01594" in text

    h4 = _read_csv(TABLE_DIR / "Table_S6a_H4_scores.csv")
    primary = {row["family"]: row for row in h4 if row["analysis"] == "primary"}
    assert round(float(primary["reproductive_assurance"]["estimate"]), 5) == -0.29830
    assert round(
        float(primary["reproductive_assurance"]["finite_publication_t_two_sided_p"]),
        5,
    ) == 0.00417
    assert round(float(primary["accessibility_generalization"]["estimate"]), 5) == -0.29566
    assert round(
        float(primary["accessibility_generalization"]["finite_publication_t_two_sided_p"]),
        5,
    ) == 0.02334
    assert "-0.29830" in text
    assert "-0.29566" in text


def test_supplement_preserves_claim_and_data_boundaries() -> None:
    text = SI.read_text(encoding="utf-8")
    assert "not mediation" in text
    assert "post-hoc functional triangulation" in text
    assert "does not distinguish species sorting" in text
    assert "before outcome unblinding" in text
    assert "10.5281/zenodo.22704973" in text
    assert "not the complete scientific analysis ledger" in text
    assert "Ecology Letters submission remains blocked" in text
