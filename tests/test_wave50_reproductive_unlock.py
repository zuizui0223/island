from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

import island_v2.wave50_incremental_all_evidence as incremental_module
import island_v2.wave50_reproductive_checkpoint as checkpoint_module
import island_v2.wave50_reproductive_overlay as overlay_module


def _write(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, compression="gzip" if path.suffix == ".gz" else None)


def test_wave50_packet_keeps_reproductive_traits_separate() -> None:
    packet = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "v2"
        / "staging"
        / "traits"
        / "wave50_source_grounded_reproductive_unlock"
    )
    direct = pd.read_csv(packet / "wave50_reviewed_direct_evidence.csv", dtype=str)
    review = pd.read_csv(packet / "wave50_source_review_audit.csv", dtype=str)
    manifest = json.loads((packet / "source_manifest.json").read_text(encoding="utf-8"))

    assert len(direct) == len(review) == 5
    assert direct[["accepted_species", "trait_name"]].duplicated().sum() == 0
    assert set(direct["trait_name"]) == {
        "autonomous_selfing_capacity",
        "mating_system",
        "self_incompatibility",
    }
    portulaca = direct.loc[direct["accepted_species"].eq("Portulaca pilosa")].iloc[0]
    assert portulaca["trait_name"] == "self_incompatibility"
    assert portulaca["normalized_value"] == "SC"
    assert not direct["trait_name"].isin({"pollen_vector_mode", "reward_type"}).any()
    assert manifest["inference_policy"] == {
        "join_key": "genus x axis x trait_name",
        "minimum_species": 3,
        "family_inference": False,
        "global_fallback": False,
        "reproductive_traits_interchangeable": False,
        "single_source_lineage_can_unlock_rule": False,
    }
    assert manifest["query_accounting"]["formal_search_api_queries"] == 779
    assert manifest["query_accounting"]["query_cost_usd"] == 0.0


def test_wave50_checkpoint_pins_counts_and_wave49(monkeypatch, tmp_path: Path) -> None:
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
    assert result["contract"] == "wave50_source_grounded_reproductive_checkpoint_v1"
    assert captured["packet_label"] == "wave50"
    assert captured["baseline_formal_run_id"] == 33397671718
    assert captured["expected_direct_rows"] == 5
    assert captured["expected_external_rows"] == 0
    assert captured["baseline_check_label"] == "immediate_wave49_baseline_pinned"


def test_wave50_incremental_declares_two_new_rules_and_two_blocks(
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
    result = incremental_module.build_wave50_incremental_audit(
        **kwargs,
        output_dir=tmp_path / "output",
    )
    assert result["contract"] == "wave50_incremental_all_evidence_touched_rule_rebuild_v1"
    assert captured["expected_new_rules"] == incremental_module.EXPECTED_NEW_RULES
    assert (
        captured["expected_counterexample_rules"]
        == incremental_module.EXPECTED_COUNTEREXAMPLE_RULES
    )
    assert captured["expected_blocked_rules"] == frozenset()
    assert captured["output_label"] == "wave50"


def test_wave50_overlay_relabels_outputs_and_pins_wave49(
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
                "resolved_wave45_direct_species_trait": 5,
                "resolved_wave45_direct_species_axis": 5,
            },
            "checks": {},
        }

    monkeypatch.setattr(overlay_module, "build_wave45_overlay", fake_build)
    checkpoint = tmp_path / "checkpoint.json"
    checkpoint.write_text(
        json.dumps(
            {
                "queries": {"formal_search_api_queries": 779},
                "query_cost_usd": 0.0,
                "checks": {
                    "retrieved_sources_verified": True,
                    "content_fingerprints_verified": True,
                },
            }
        ),
        encoding="utf-8",
    )
    output = tmp_path / "output"
    summary = overlay_module.build_wave50_overlay(
        baseline_csv=tmp_path / "baseline.csv.gz",
        previous_rule_audit_csv=tmp_path / "rules.csv.gz",
        all_evidence_dir=tmp_path / "audit",
        direct_evidence_csv=tmp_path / "direct.csv.gz",
        external_evidence_csv=tmp_path / "external.csv.gz",
        checkpoint_summary_json=checkpoint,
        output_dir=output,
        expected_species=3,
    )
    assert summary["baseline_formal_run_id"] == 33397671718
    assert summary["formal_wave33_run_id"] == 32932103226
    assert summary["delta"]["resolved_wave50_direct_species_trait"] == 5
    assert summary["checks"]["reproductive_traits_not_interchanged"] is True
    assert (output / "wave50_species_axis_coverage.csv.gz").is_file()
    assert not list(output.glob("wave45_*"))
