import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "config/chapter1_h4_prospective_validation_summary_lock.json"
WILD = ROOT / "config/chapter1_h4_prospective_temporal_support_decision_lock.json"
CROP = ROOT / "config/chapter1_h4_pollimcrop_transportability_preflight_result_lock.json"


def test_h4_prospective_validation_summary_preserves_outcome_blinding():
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    wild = json.loads(WILD.read_text(encoding="utf-8"))
    crop = json.loads(CROP.read_text(encoding="utf-8"))

    assert summary["status"] == "prospective_validation_support_gates_not_met"
    assert summary["wild_temporal"]["outcomes_read"] is False
    assert summary["crop_transportability"]["outcomes_read"] is False
    assert summary["integrated_decision"]["any_prospective_outcome_unblinded"] is False

    assert wild["decision"]["outcome_extraction_authorized"] is False
    assert crop["decision"]["outcome_extraction_authorized"] is False
    assert crop["decision"]["admitted_hypotheses"] == []


def test_h4_crop_support_failure_is_balance_limited_not_total_overlap_limited():
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    h4a = summary["crop_transportability"]["H4a_reproductive_assurance"]
    h4b = summary["crop_transportability"]["H4b_accessibility_generalization"]

    assert h4a["matched_species_total"] >= h4a["required_min_species_total"]
    assert h4a["publications_total"] >= h4a["required_min_publications_total"]
    assert h4a["species_state_0"] < h4a["required_min_species_per_state"]

    assert h4b["matched_species_total"] >= h4b["required_min_species_total"]
    assert h4b["publications_total"] >= h4b["required_min_publications_total"]
    assert h4b["score_sd"] >= h4b["required_min_score_sd"]
    assert h4b["low_score_species"] < h4b["required_min_low_score_species"]


def test_no_crop_outcome_result_lock_exists_after_support_failure():
    result_lock = ROOT / "config/chapter1_h4_pollimcrop_transportability_result_lock.json"
    assert not result_lock.exists()
