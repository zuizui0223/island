from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

import island_v2.wave46_corrected_public_web_checkpoint as checkpoint_module
import island_v2.wave46_corrected_public_web_overlay as overlay_module

AXES = (
    "floral_structural_complexity",
    "flower_colour",
    "reproductive_assurance",
)


def test_integrated_workflow_pins_formal_wave33_and_latest_public_web() -> None:
    repo = Path(__file__).resolve().parents[1]
    config = json.loads(
        (repo / "config" / "integrated_trait_coverage_sources.json").read_text(encoding="utf-8")
    )
    public_web = next(
        source for source in config["sources"] if source["source_group"] == "latest_public_web"
    )
    assert public_web == {
        "source_group": "latest_public_web",
        "run_id": 32710232989,
        "artifact_name": "reviewed-open-web-evidence-32710232989",
        "artifact_id": 9514635137,
        "artifact_digest": (
            "sha256:54c2245876b3c2aed51d7d5ee23a8dd106ec436607178051d20b748ca2624967"
        ),
        "role": "latest_species_direct_web",
    }

    workflow = (repo / ".github" / "workflows" / "build-integrated-trait-coverage.yml").read_text(
        encoding="utf-8"
    )
    assert 'default: "32932103226"' in workflow
    assert 'default: "integrated-trait-coverage-32932103226"' in workflow
    assert (
        'default: "sha256:d6405d520414b802fe75b2e450f45b8e74295f95049f20fb7ce61f57471326d1"'
        in workflow
    )
    assert 'default: "32710232989"' in workflow
    assert 'default: "9514635137"' in workflow


def _write(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, compression="gzip" if path.suffix == ".gz" else None)


def _direct_row() -> dict[str, str]:
    return {
        "accepted_species": "Ardisia humilis",
        "axis": "reproductive_assurance",
        "trait_name": "autonomous_selfing_capacity",
        "normalized_value": "autonomous",
        "quality": "high",
        "source_group": "wave46_primary_reproductive_articles",
        "source_provider": "Ye et al. 2023 primary article",
        "source_url": "https://example.test/ardisia",
        "source_record_id": "wave46:ardisia:autonomous",
        "source_citation": "Ye et al. 2023",
        "source_excerpt": "autogamy seed set was 52.29%",
        "evidence_scope": "species_direct",
        "name_match_method": "accepted_name_exact",
        "source_lineage": "doi:10.11931/example",
        "lineage_method": "original_article_doi",
        "source_run_id": "web:retrieved-20260831",
        "source_artifact": "example",
        "source_file": "html_sha256:abc",
        "acceptance_contract": "primary_article_controlled_pollination_v1",
    }


def test_wave46_checkpoint_pins_latest_public_web_and_direct_quote(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(checkpoint_module, "PUBLIC_WEB_ROWS", 1)
    species = ("Ardisia humilis", "Ardisia escallonioides", "Ardisia sieboldii")
    target = pd.DataFrame(
        [
            {"accepted_species": accepted_species, "axis": axis}
            for accepted_species in species
            for axis in AXES
        ]
    )
    target_path = tmp_path / "target.csv.gz"
    _write(target, target_path)
    public_web_path = tmp_path / "public-web.csv.gz"
    _write(
        pd.DataFrame(
            [
                {
                    "accepted_species": "Ardisia escallonioides",
                    "axis": "flower_colour",
                    "trait_name": "flower_primary_color",
                    "source_lineage": "url:https://example.test/flora",
                    "source_run_id": "32710232989",
                    "source_artifact": "reviewed-open-web-evidence-32710232989",
                }
            ]
        ),
        public_web_path,
    )
    packet = tmp_path / "packet"
    packet.mkdir()
    direct = _direct_row()
    pd.DataFrame([direct]).to_csv(packet / "wave46_reviewed_direct_evidence.csv", index=False)
    pd.DataFrame(
        [
            {
                "record_id": direct["source_record_id"],
                "exact_supporting_quote": direct["source_excerpt"],
                "accepted_correct": "true",
            }
        ]
    ).to_csv(packet / "wave46_source_review_audit.csv", index=False)
    pd.DataFrame(
        [
            {
                "accepted_species": "Ardisia humilis",
                "target_universe_status": "in_fixed_106295_species",
                "name_match_method": "accepted_name_exact",
            }
        ]
    ).to_csv(packet / "wave46_identity_audit.csv", index=False)
    (packet / "source_manifest.json").write_text(
        json.dumps(
            {
                "fixed_target_species": 3,
                "formal_wave33_baseline": {"run_id": 32932103226},
                "correct_latest_public_web": {
                    "run_id": 32710232989,
                    "evidence_rows": 1,
                },
                "new_primary_sources": [{"retrieved_content_sha256": "abc", "reviewed_rows": 1}],
                "inference_policy": {
                    "join_key": "genus x axis x trait_name",
                    "family_inference": False,
                    "global_fallback": False,
                    "reproductive_traits_interchangeable": False,
                },
            }
        ),
        encoding="utf-8",
    )
    summary = checkpoint_module.validate_packet(
        packet_dir=packet,
        latest_public_web_csv=public_web_path,
        target_coverage_csv=target_path,
        output_dir=tmp_path / "output",
        output_json=tmp_path / "output" / "summary.json",
        expected_species=3,
    )
    assert summary["evidence"]["latest_public_web_rows"] == 1
    assert summary["evidence"]["new_direct_species_trait"] == 1
    assert summary["checks"]["formal_wave33_baseline_pinned"] is True
    assert summary["checks"]["latest_public_web_artifact_corrected"] is True


def test_wave46_overlay_relabels_artifacts_and_formal_runs(tmp_path: Path, monkeypatch) -> None:
    old_names = (
        "wave45_species_axis_coverage.csv.gz",
        "wave45_resolved_direct_species_trait.csv.gz",
        "wave45_new_validated_low_species_trait.csv.gz",
        "wave45_new_trait_specific_genus_rule_audit.csv.gz",
        "wave45_change_audit.csv.gz",
        "wave45_external_congener_resolved_species_trait.csv.gz",
        "wave45_external_congener_source_conflicts.csv.gz",
    )

    def fake_build_wave45_overlay(**kwargs):
        output_dir = kwargs["output_dir"]
        output_dir.mkdir(parents=True, exist_ok=True)
        for name in old_names:
            _write(pd.DataFrame([{"value": "ok"}]), output_dir / name)
        return {
            "contract": "wave45",
            "baseline_formal_run_id": 1,
            "delta": {
                "resolved_wave45_direct_species_trait": 2,
                "resolved_wave45_direct_species_axis": 1,
            },
            "checks": {},
            "artifact_sha256": {},
        }

    monkeypatch.setattr(overlay_module, "build_wave45_overlay", fake_build_wave45_overlay)
    checkpoint = tmp_path / "checkpoint.json"
    checkpoint.write_text(
        json.dumps(
            {
                "evidence": {
                    "latest_public_web_rows": 118884,
                    "new_direct_species_trait": 2,
                }
            }
        ),
        encoding="utf-8",
    )
    external = tmp_path / "external.csv.gz"
    _write(
        pd.DataFrame([{"accepted_species": "External example", "trait_name": "mating_system"}]),
        external,
    )
    output = tmp_path / "output"
    summary = overlay_module.build_wave46_overlay(
        baseline_csv=tmp_path / "unused-baseline.csv.gz",
        previous_rule_audit_csv=tmp_path / "unused-rules.csv.gz",
        all_evidence_dir=tmp_path / "unused-audit",
        direct_evidence_csv=tmp_path / "unused-direct.csv.gz",
        external_evidence_csv=external,
        checkpoint_summary_json=checkpoint,
        output_dir=output,
        expected_species=3,
    )
    assert summary["baseline_formal_run_id"] == 33370692122
    assert summary["formal_wave33_run_id"] == 32932103226
    assert summary["correct_public_web_run_id"] == 32710232989
    assert summary["checks"]["correct_latest_public_web_artifact_used"] is True
    assert (output / "wave46_species_axis_coverage.csv.gz").is_file()
    assert not list(output.glob("wave45_*"))
