from __future__ import annotations

import json
from pathlib import Path

from island_v2 import wave42_external_bulk_overlay as overlay


def test_wave42_overlay_records_incremental_and_combined_external_counts(
    tmp_path: Path, monkeypatch
) -> None:
    audit_dir = tmp_path / "audit"
    audit_dir.mkdir()
    (audit_dir / "all_evidence_trait_coverage_summary.json").write_text(
        json.dumps(
            {
                "source_lineage_audit": {
                    "external_congener_support": {
                        "rows": 2_981,
                        "resolved_species_trait": 2_806,
                        "cell_resolution_classification_counts": {
                            "single_independent_lineage": 2_785,
                            "independent_source_agreement": 21,
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
            {
                "evidence": {
                    "species_trait": 1_318,
                    "shared_redistribution_guard_lineages": 1,
                }
            }
        ),
        encoding="utf-8",
    )
    output_dir = tmp_path / "output"

    def fake_build_overlay(**kwargs):
        assert kwargs["wave_label"] == "wave42"
        assert kwargs["baseline_formal_run_id"] == 33355297355
        kwargs["output_dir"].mkdir()
        return {
            "baseline_formal_run_id": 33355297355,
            "delta": {
                "external_resolved_species_trait": 2_806,
                "external_direct_conflicts": 28,
            },
            "checks": {},
        }

    monkeypatch.setattr(overlay, "build_overlay", fake_build_overlay)
    summary = overlay.build_wave42_overlay(
        baseline_csv=tmp_path / "baseline.csv",
        previous_rule_audit_csv=tmp_path / "rules.csv",
        all_evidence_dir=audit_dir,
        external_evidence_csv=tmp_path / "evidence.csv",
        checkpoint_summary_json=checkpoint,
        output_dir=output_dir,
        expected_species=3,
    )

    assert summary["baseline_formal_run_id"] == 33355297355
    assert summary["delta"] == {
        "new_external_input_species_trait": 1_318,
        "combined_external_evidence_rows": 2_981,
        "combined_external_resolved_species_trait": 2_806,
        "combined_external_direct_conflicts": 28,
    }
    assert summary["checks"]["shared_meta_database_lineage_guard"] is True
    saved = json.loads(
        (output_dir / "wave42_coverage_summary.json").read_text(encoding="utf-8")
    )
    assert saved == summary
