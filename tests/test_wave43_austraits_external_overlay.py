from __future__ import annotations

import json
from pathlib import Path

from island_v2 import wave43_austraits_external_overlay as overlay


def test_wave43_overlay_records_incremental_and_combined_external_counts(
    tmp_path: Path, monkeypatch
) -> None:
    audit_dir = tmp_path / "audit"
    audit_dir.mkdir()
    (audit_dir / "all_evidence_trait_coverage_summary.json").write_text(
        json.dumps(
            {
                "source_lineage_audit": {
                    "external_congener_support": {
                        "rows": 5_973,
                        "resolved_species_trait": 5_184,
                        "cell_resolution_classification_counts": {
                            "single_independent_lineage": 5_001,
                            "independent_source_agreement": 125,
                            "true_multistate_variable": 58,
                            "unresolved_direct_conflict": 28,
                        },
                    }
                }
            }
        ),
        encoding="utf-8",
    )
    checkpoint = tmp_path / "checkpoint.json"
    checkpoint.write_text(
        json.dumps(
            {"evidence": {"species_trait": 2_378, "underlying_source_lineages": 5}}
        ),
        encoding="utf-8",
    )
    output_dir = tmp_path / "output"

    def fake_build_overlay(**kwargs):
        assert kwargs["wave_label"] == "wave43"
        assert kwargs["baseline_formal_run_id"] == 33_357_283_159
        kwargs["output_dir"].mkdir()
        return {
            "baseline_formal_run_id": 33_357_283_159,
            "delta": {
                "external_resolved_species_trait": 5_184,
                "external_direct_conflicts": 28,
            },
            "checks": {},
        }

    monkeypatch.setattr(overlay, "build_overlay", fake_build_overlay)
    summary = overlay.build_wave43_overlay(
        baseline_csv=tmp_path / "baseline.csv",
        previous_rule_audit_csv=tmp_path / "rules.csv",
        all_evidence_dir=audit_dir,
        external_evidence_csv=tmp_path / "evidence.csv",
        checkpoint_summary_json=checkpoint,
        output_dir=output_dir,
        expected_species=3,
    )

    assert summary["delta"] == {
        "new_external_input_species_trait": 2_378,
        "combined_external_evidence_rows": 5_973,
        "combined_external_resolved_species_trait": 5_184,
        "combined_external_direct_conflicts": 28,
    }
    assert summary["checks"]["underlying_source_lineage_guard"] is True
    saved = json.loads(
        (output_dir / "wave43_coverage_summary.json").read_text(encoding="utf-8")
    )
    assert saved == summary
