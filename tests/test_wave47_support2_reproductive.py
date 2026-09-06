from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

import island_v2.wave47_support2_reproductive_checkpoint as checkpoint_module
import island_v2.wave47_support2_reproductive_overlay as overlay_module
from island_v2.all_evidence_trait_audit import build_rule_audit
from island_v2.wave47_incremental_all_evidence import build_incremental_audit

AXES = (
    "floral_structural_complexity",
    "flower_colour",
    "reproductive_assurance",
)


def _write(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, compression="gzip" if path.suffix == ".gz" else None)


def _external_rows(source_sha: str = "abc") -> list[dict[str, str]]:
    base = {
        "accepted_species": "Buxus wallichiana",
        "axis": "reproductive_assurance",
        "quality": "high",
        "source_group": "wave47_primary_reproductive_articles",
        "source_provider": "Sharma et al. 2020 primary article",
        "source_url": "https://example.test/buxus.pdf",
        "source_citation": "Sharma et al. 2020",
        "evidence_scope": "external_congener_species_direct",
        "name_match_method": "strict_wfo_gbif_two_backbone",
        "source_lineage": "currentscience:118:7:1021-1022:sharma-et-al-2020",
        "lineage_method": "original_article_bibliographic_identity",
        "source_run_id": "web:retrieved-20260831",
        "source_artifact": "currentscience-118-7-1021",
        "source_file": f"pdf_sha256:{source_sha}",
        "acceptance_contract": "external_congener_species_direct_strict_two_backbone_v1",
    }
    return [
        {
            **base,
            "trait_name": "self_incompatibility",
            "normalized_value": "SC",
            "source_record_id": "wave47:buxus:sc",
            "source_excerpt": "the species is highly self-compatible",
        },
        {
            **base,
            "trait_name": "autonomous_selfing_capacity",
            "normalized_value": "autonomous",
            "source_record_id": "wave47:buxus:autonomous",
            "source_excerpt": "bagging resulted in 100% fruit set",
        },
    ]


def _coverage(species: tuple[str, ...]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "accepted_species": accepted_species,
                "axis": axis,
                "trait_composition": "",
                "trait_names": "",
                "source_groups": "",
                "source_lineages": "",
                "quality": "",
            }
            for accepted_species in species
            for axis in AXES
        ]
    )


def test_wave47_checkpoint_pins_source_and_strict_external_identity(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(checkpoint_module, "SOURCE_ROWS", 2)
    species = ("Buxus alpha", "Buxus beta", "Buxus gamma")
    target_path = tmp_path / "target.csv.gz"
    _write(_coverage(species), target_path)
    source = tmp_path / "source.pdf"
    source.write_bytes(b"frozen-current-science-source")
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()

    packet = tmp_path / "packet"
    packet.mkdir()
    evidence = pd.DataFrame(_external_rows(source_sha))
    evidence.to_csv(packet / "wave47_external_congener_evidence.csv", index=False)
    pd.DataFrame(
        [
            {
                "record_id": row["source_record_id"],
                "exact_supporting_quote": row["source_excerpt"],
                "accepted_correct": "true",
            }
            for row in evidence.to_dict("records")
        ]
    ).to_csv(packet / "wave47_source_review_audit.csv", index=False)
    pd.DataFrame(
        [
            {
                "accepted_species": "Buxus wallichiana",
                "target_universe_status": "external_congener_only",
                "family": "Buxaceae",
                "name_match_method": "strict_wfo_gbif_two_backbone",
                "wfo_match_id": "wfo-0000576661",
                "wfo_status": "accepted",
                "wfo_rank": "species",
                "gbif_usage_key": "3988068",
                "gbif_status": "ACCEPTED",
                "gbif_rank": "SPECIES",
                "gbif_match_type": "EXACT",
                "gbif_confidence": "99",
                "gbif_family": "Buxaceae",
            }
        ]
    ).to_csv(packet / "wave47_identity_audit.csv", index=False)
    (packet / "source_manifest.json").write_text(
        json.dumps(
            {
                "fixed_target_species": 3,
                "formal_wave33_baseline": {"run_id": 32932103226},
                "immediate_formal_baseline": {"run_id": 33376311877},
                "source": {"reviewed_rows": 2, "retrieved_content_sha256": source_sha},
                "inference_policy": {
                    "join_key": "genus x axis x trait_name",
                    "minimum_species": 3,
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
        target_coverage_csv=target_path,
        retrieved_source_file=source,
        output_dir=tmp_path / "output",
        output_json=tmp_path / "output" / "summary.json",
        expected_species=3,
    )
    assert summary["evidence"]["new_external_species_trait"] == 2
    assert summary["checks"]["retrieved_source_digest_verified"] is True
    assert summary["checks"]["external_strict_two_backbone"] is True


def test_wave47_incremental_rebuild_unlocks_only_matching_trait(tmp_path: Path) -> None:
    species = ("Buxus alpha", "Buxus beta", "Buxus gamma")
    master = pd.DataFrame(
        [
            {
                "accepted_species": name,
                "genus": "Buxus",
                "family": "Buxaceae",
                "n_islands": "1",
                "n_records": "1",
            }
            for name in species
        ]
    )
    master_path = tmp_path / "master.csv"
    master.to_csv(master_path, index=False)
    baseline_path = tmp_path / "baseline.csv.gz"
    _write(_coverage(species), baseline_path)

    direct = pd.DataFrame(
        [
            {
                "accepted_species": name,
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "classification": "single_independent_lineage",
                "resolution_status": "resolved",
                "quality": "high",
                "state_set": '["SC"]',
                "normalized_value": "SC",
                "source_groups": "primary",
                "source_lineages": f"doi:{name}",
                "genus": "Buxus",
            }
            for name in species[:2]
        ]
    )
    direct_path = tmp_path / "direct.csv.gz"
    _write(direct, direct_path)
    lineages = pd.DataFrame(
        [
            {
                "accepted_species": row["accepted_species"],
                "genus": "Buxus",
                "axis": row["axis"],
                "trait_name": row["trait_name"],
                "state_set": row["state_set"],
                "source_lineage": row["source_lineages"],
                "source_group": "primary",
                "ontology_valid": True,
                "lineage_internal_conflict": False,
            }
            for row in direct.to_dict("records")
        ]
    )
    old_low = pd.DataFrame(columns=["accepted_species", "genus", "axis", "trait_name", "state_set"])
    previous_rules = build_rule_audit(direct, lineages, old_low)
    previous_rules_path = tmp_path / "rules.csv.gz"
    _write(previous_rules, previous_rules_path)

    external_columns = list(direct.columns)
    external_path = tmp_path / "external.csv.gz"
    _write(pd.DataFrame(columns=external_columns), external_path)
    conflicts_path = tmp_path / "conflicts.csv.gz"
    _write(
        pd.DataFrame(
            columns=[
                "accepted_species",
                "axis",
                "trait_name",
                "classification",
                "resolution_status",
            ]
        ),
        conflicts_path,
    )
    old_low_path = tmp_path / "old-low.csv.gz"
    _write(old_low, old_low_path)
    new_external_path = tmp_path / "new-external.csv.gz"
    _write(pd.DataFrame(_external_rows()), new_external_path)

    summary = build_incremental_audit(
        master_csv=master_path,
        ontology_yaml=Path(__file__).resolve().parents[1] / "config" / "trait_ontology.yml",
        baseline_coverage_csv=baseline_path,
        previous_rule_audit_csv=previous_rules_path,
        previous_resolved_direct_csv=direct_path,
        previous_external_resolved_csv=external_path,
        previous_external_conflicts_csv=conflicts_path,
        previous_rebuilt_low_csv=old_low_path,
        new_external_evidence_csv=new_external_path,
        output_dir=tmp_path / "audit",
        expected_species=3,
    )
    assert summary["new_eligible_rules"] == ["Buxus x self_incompatibility"]
    assert summary["new_rule"] == {
        "n_direct_species": 3,
        "dominance": 1.0,
        "species_loo_accuracy": 1.0,
        "lineage_loo_accuracy": 1.0,
        "inferred_value": "SC",
    }
    assert summary["new_validated_low_species_trait"] == 1
    low = pd.read_csv(tmp_path / "audit" / "rebuilt_all_evidence_validated_low.csv.gz")
    assert low[["accepted_species", "trait_name"]].to_dict("records") == [
        {"accepted_species": "Buxus gamma", "trait_name": "self_incompatibility"}
    ]


def test_wave47_overlay_relabels_outputs_and_pins_runs(tmp_path: Path, monkeypatch) -> None:
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
        output = kwargs["output_dir"]
        output.mkdir(parents=True, exist_ok=True)
        for name in old_names:
            _write(pd.DataFrame([{"value": "ok"}]), output / name)
        return {
            "contract": "wave45",
            "delta": {
                "resolved_wave45_direct_species_trait": 0,
                "resolved_wave45_direct_species_axis": 0,
            },
            "checks": {},
            "artifact_sha256": {},
        }

    monkeypatch.setattr(overlay_module, "build_wave45_overlay", fake_build_wave45_overlay)
    checkpoint = tmp_path / "checkpoint.json"
    checkpoint.write_text(
        json.dumps({"checks": {"retrieved_source_digest_verified": True}}),
        encoding="utf-8",
    )
    output = tmp_path / "output"
    summary = overlay_module.build_wave47_overlay(
        baseline_csv=tmp_path / "unused-baseline.csv.gz",
        previous_rule_audit_csv=tmp_path / "unused-rules.csv.gz",
        all_evidence_dir=tmp_path / "unused-audit",
        direct_evidence_csv=tmp_path / "unused-direct.csv.gz",
        external_evidence_csv=tmp_path / "unused-external.csv.gz",
        checkpoint_summary_json=checkpoint,
        output_dir=output,
        expected_species=3,
    )
    assert summary["formal_wave33_run_id"] == 32932103226
    assert summary["baseline_formal_run_id"] == 33376311877
    assert summary["checks"]["external_source_pdf_digest_verified"] is True
    assert (output / "wave47_species_axis_coverage.csv.gz").is_file()
    assert not list(output.glob("wave45_*"))
