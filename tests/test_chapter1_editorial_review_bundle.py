import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "scripts/submission/build_chapter1_editorial_review.py"
spec = importlib.util.spec_from_file_location("chapter1_review", SRC)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_current_h1_h4_submission_contract_and_word_limits():
    selector = json.loads((ROOT / "config/chapter1_submission_current.json").read_text())
    assert selector["status"] == "primary_submission_baseline"
    assert selector["manuscript"] == "submission/chapter1_current/MANUSCRIPT.md"
    manuscript = (ROOT / selector["manuscript"]).read_text(encoding="utf-8")
    count = module._word_counts(manuscript)
    assert 100 <= count["abstract"] <= 150
    assert 3000 <= count["main"] <= 5000


def test_review_allowlist_excludes_restricted_trait_source_data():
    approved = set(module.DOC_FILES) | set(module.TABLE_FILES)
    assert "MANUSCRIPT.md" in approved
    assert "SUPPLEMENTARY_INFORMATION.md" in approved
    assert "Table_S5_H3_pollen_limitation.csv" in approved
    forbidden = (
        "ledger", "species_treatment", "raw_trait", "sources",
        "h3_corrected_effect_rows", "occurrence", "evidence_excerpt",
    )
    assert all(not any(term in name.lower() for term in forbidden) for name in approved)
    assert module.REVIEW_WARNING.startswith("EDITORIAL REVIEW DRAFT")
    script = SRC.read_text(encoding="utf-8")
    assert '"contains_restricted_trait_ledger": False' in script
    assert '"ready_for_submission": False' in script
    assert "chapter1-submission-v1" not in script


def test_figure_manifest_is_complete_and_matches_frozen_files():
    figures = ROOT / "submission/chapter1_current/figures"
    meta = json.loads((figures / "FIGURE_MANIFEST.json").read_text())
    assert [x["figure"] for x in meta["figures"]] == list(range(1, 7))
    for row in meta["figures"]:
        assert (figures / row["pdf"]).stat().st_size == row["pdf_bytes"]
        assert (figures / row["svg"]).stat().st_size == row["svg_bytes"]


def test_offshore_slope_kept_as_posthoc_supporting_evidence():
    manuscript = (ROOT / "submission/chapter1_current/MANUSCRIPT.md").read_text()
    assert "β = 0.09191" in manuscript
    assert "β = 0.22031" in manuscript
    assert "cannot be directly compared" in manuscript
    assert "not fewer pollinators or historical decline" in manuscript
    assert "not experimentally assigned" in manuscript
