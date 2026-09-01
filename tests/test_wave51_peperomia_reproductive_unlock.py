from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

import island_v2.wave51_incremental_all_evidence as incremental_module
import island_v2.wave51_reproductive_checkpoint as checkpoint_module
import island_v2.wave51_reproductive_overlay as overlay_module


def _write(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, compression="gzip" if path.suffix == ".gz" else None)


def _packet() -> Path:
    return (
        Path(__file__).resolve().parents[1]
        / "data"
        / "v2"
        / "staging"
        / "traits"
        / "wave51_peperomia_reproductive_unlock"
    )


def test_wave51_packet_separates_direct_autonomy_from_external_compatibility() -> None:
    direct = pd.read_csv(_packet() / "wave51_reviewed_direct_evidence.csv", dtype=str)
    external = pd.read_csv(_packet() / "wave51_external_congener_evidence.csv", dtype=str)
    review = pd.read_csv(_packet() / "wave51_source_review_audit.csv", dtype=str)
    manifest = json.loads((_packet() / "source_manifest.json").read_text(encoding="utf-8"))

    assert len(direct) == 4
    assert len(external) == 2
    assert len(review) == 6
    assert set(direct["trait_name"]) == {"autonomous_selfing_capacity"}
    assert set(external["trait_name"]) == {"self_incompatibility"}
    assert set(external["normalized_value"]) == {"SC"}
    assert direct["source_lineage"].nunique() == 1
    assert external["source_lineage"].nunique() == 2
    assert not pd.concat([direct, external])["trait_name"].isin(
        {"pollen_vector_mode", "reward_type"}
    ).any()
    assert manifest["inference_policy"]["join_key"] == "genus x axis x trait_name"
    assert manifest["inference_policy"]["single_source_lineage_can_unlock_rule"] is False
    manual = manifest["sources"][0]
    assert manual["retrieval_mode"] == "manual_browser_pinned"
    assert manual["manual_receipt"]["binary_available_in_ci"] is False
    assert manual["manual_receipt"]["reviewed_page_numbers"] == [55, 56]


def test_wave51_checkpoint_pins_counts_and_wave50(monkeypatch, tmp_path: Path) -> None:
    captured = {}

    def fake_validate(**kwargs):
        captured.update(kwargs)
        return {"contract": kwargs["contract"]}

    monkeypatch.setattr(checkpoint_module, "validate_reproductive_packet", fake_validate)
    result = checkpoint_module.validate_packet(
        packet_dir=tmp_path,
        target_coverage_csv=tmp_path / "coverage.csv.gz",
        retrieved_source_dir=tmp_path / "sources",
        output_dir=tmp_path / "output",
        output_json=tmp_path / "summary.json",
    )
    assert result["contract"] == "wave51_peperomia_reproductive_checkpoint_v1"
    assert captured["packet_label"] == "wave51"
    assert captured["baseline_formal_run_id"] == 33406678174
    assert captured["expected_direct_rows"] == 4
    assert captured["expected_external_rows"] == 2
    assert captured["baseline_check_label"] == "immediate_wave50_baseline_pinned"


def test_wave51_incremental_declares_trait_rule_and_single_lineage_block(
    monkeypatch, tmp_path: Path
) -> None:
    captured = {}

    def fake_build(**kwargs):
        captured.update(kwargs)
        return {"contract": kwargs["contract"]}

    monkeypatch.setattr(incremental_module, "build_incremental_audit", fake_build)
    kwargs = {
        name: tmp_path / f"{name}.csv.gz"
        for name in (
            "master_csv",
            "ontology_yaml",
            "baseline_coverage_csv",
            "previous_rule_audit_csv",
            "previous_resolved_direct_csv",
            "previous_external_resolved_csv",
            "previous_external_conflicts_csv",
            "previous_rebuilt_low_csv",
            "new_direct_evidence_csv",
            "new_external_evidence_csv",
        )
    }
    result = incremental_module.build_wave51_incremental_audit(
        **kwargs, output_dir=tmp_path / "output"
    )
    assert result["contract"] == "wave51_incremental_all_evidence_touched_rule_rebuild_v1"
    assert captured["expected_new_rules"] == incremental_module.EXPECTED_NEW_RULES
    assert (
        captured["expected_blocked_rules"]
        == incremental_module.EXPECTED_SINGLE_LINEAGE_BLOCKED_RULES
    )
    assert captured["expected_counterexample_rules"] == frozenset()
    assert captured["output_label"] == "wave51"


def test_wave51_overlay_relabels_outputs_and_preserves_manual_receipt(
    tmp_path: Path, monkeypatch
) -> None:
    old_names = (
        "wave45_species_axis_coverage.csv.gz",
        "wave45_resolved_direct_species_trait.csv.gz",
        "wave45_new_validated_low_species_trait.csv.gz",
        "wave45_new_trait_specific_genus_rule_audit.csv.gz",
        "wave45_change_audit.csv.gz",
        "wave45_external_congener_resolved_species_trait.csv.gz",
        "wave45_external_congener_source_conflicts.csv.gz",
    )

    def fake_build(**kwargs):
        output = kwargs["output_dir"]
        output.mkdir(parents=True, exist_ok=True)
        for name in old_names:
            _write(pd.DataFrame([{"value": "ok"}]), output / name)
        return {
            "delta": {
                "resolved_wave45_direct_species_trait": 4,
                "resolved_wave45_direct_species_axis": 4,
            },
            "checks": {},
        }

    monkeypatch.setattr(overlay_module, "build_wave45_overlay", fake_build)
    checkpoint = tmp_path / "checkpoint.json"
    checkpoint.write_text(
        json.dumps(
            {
                "queries": {"formal_search_api_queries": 4},
                "query_cost_usd": 0.0,
                "source_binary_coverage": {"verified_in_run": 2, "manual_receipt_only": 1},
                "limitations": {
                    "manual_source_binary_not_in_artifact": [
                        "figueiredo-1997-peperaceae-thesis"
                    ]
                },
                "checks": {
                    "retrieved_sources_verified": True,
                    "manual_browser_receipts_verified": True,
                    "content_fingerprints_verified": True,
                },
            }
        ),
        encoding="utf-8",
    )
    output = tmp_path / "output"
    summary = overlay_module.build_wave51_overlay(
        baseline_csv=tmp_path / "baseline.csv.gz",
        previous_rule_audit_csv=tmp_path / "rules.csv.gz",
        all_evidence_dir=tmp_path / "audit",
        direct_evidence_csv=tmp_path / "direct.csv.gz",
        external_evidence_csv=tmp_path / "external.csv.gz",
        checkpoint_summary_json=checkpoint,
        output_dir=output,
        expected_species=3,
    )
    assert summary["baseline_formal_run_id"] == 33406678174
    assert summary["formal_wave33_run_id"] == 32932103226
    assert summary["delta"]["resolved_wave51_direct_species_trait"] == 4
    assert summary["source_binary_coverage"]["manual_receipt_only"] == 1
    assert summary["checks"]["reproductive_traits_not_interchanged"] is True
    assert (output / "wave51_species_axis_coverage.csv.gz").is_file()
    assert not list(output.glob("wave45_*"))
