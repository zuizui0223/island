import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
CANONICAL = ROOT / "config/chapter1_v14_canonical_result_lock.json"
HIERARCHY = ROOT / "config/chapter1_v14_h4_evidence_hierarchy.yml"
SUMMARY = ROOT / "config/chapter1_h4_prospective_validation_summary_lock.json"
MANUSCRIPT = ROOT / "docs/chapter1_manuscript_full_v14_reordered_hypotheses_20260918.md"
FREEZE = ROOT / "docs/chapter1_submission_freeze_v14_20260918.md"
README = ROOT / "README.md"
PIPELINE = ROOT / "docs/PAPER_PIPELINE.md"


def test_postfreeze_validation_does_not_rewrite_canonical_core_lock():
    canonical = json.loads(CANONICAL.read_text(encoding="utf-8"))
    assert canonical["contract"] == "chapter1_v14_canonical_result_lock_v1"
    assert canonical["promotion_provenance"]["ci_run_id"] == 35314955780
    assert canonical["promotion_provenance"]["ci_artifact_id"] == 10535020072
    # The core lock preserves the pre-addendum snapshot by design.
    assert (
        canonical["H4"]["independent_crop_domain"]["status"]
        == "pending_support_gated_execution"
    )


def test_current_h4_hierarchy_uses_completed_support_limited_audit():
    hierarchy = yaml.safe_load(HIERARCHY.read_text(encoding="utf-8"))
    crop = hierarchy["evidence_layers"]["pollimcrop_independent_domain"]
    wild = hierarchy["evidence_layers"]["prospective_post2015_wild"]

    assert wild["status"] == "support_gate_not_met_not_evaluable"
    assert wild["outcomes_unblinded_for_primary_test"] is False

    assert crop["status"] == "support_gate_not_met_not_evaluable"
    assert crop["admitted_hypotheses"] == []
    assert crop["outcome_extraction_authorized"] is False
    assert crop["outcomes_read"] is False
    assert crop["verification_workflow_run_id"] == 35321319394
    assert crop["H4a_support"]["failed_gate"] == "species_state_0"
    assert crop["H4b_support"]["failed_gate"] == "low_score_species"


def test_integrated_prospective_summary_keeps_all_outcomes_blinded():
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    decision = summary["integrated_decision"]

    assert decision["confirmatory_wild_replication_available"] is False
    assert decision["independent_crop_transportability_test_evaluable"] is False
    assert decision["any_prospective_outcome_unblinded"] is False
    assert decision["threshold_relaxation_authorized"] is False
    assert decision["existing_posthoc_H4_discovery_reclassified_as_confirmatory"] is False
    assert decision["biological_null_inferred"] is False


def test_publication_surface_reports_support_failure_not_biological_null():
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    freeze = FREEZE.read_text(encoding="utf-8")
    readme = README.read_text(encoding="utf-8")
    pipeline = PIPELINE.read_text(encoding="utf-8")

    assert "No prospective outcome was unblinded in either route" in manuscript
    assert "support/design failures, not evidence against the biological predictions" in manuscript
    assert "## Post-freeze H4 validation addendum" in freeze
    assert "chapter1_h4_prospective_validation_summary_lock.json" in readme
    assert "chapter1_h4_prospective_validation_summary_lock.json" in pipeline
