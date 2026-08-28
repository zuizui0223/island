from pathlib import Path

import pandas as pd

from island_v2.wave35_trait_overlay import AXES
from island_v2.wave40_source_recovery_overlay import build_source_recovery_overlay


def _write_csv(frame: pd.DataFrame, path: Path) -> Path:
    frame.to_csv(path, index=False)
    return path


def test_source_recovery_overlay_is_lossless_and_trait_specific(tmp_path: Path) -> None:
    species = ["Alpha one", "Alpha two", "Beta one"]
    baseline = pd.DataFrame(
        [
            {
                "accepted_species": name,
                "axis": axis,
                "trait_composition": "",
                "trait_names": "",
                "source_groups": "",
                "source_lineages": "",
                "quality": "",
            }
            for name in species
            for axis in AXES
        ]
    )
    baseline_path = _write_csv(baseline, tmp_path / "baseline.csv")

    reviewed = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "trait_name": "self_incompatibility",
                "normalized_value": "SC",
                "quality": "high",
                "source_group": "reviewed_database",
                "source_provider": "provider",
                "source_url": "https://example.org/alpha-one",
                "source_record_id": "record-1",
                "source_citation": "Published database",
                "source_excerpt": "Alpha one is self-compatible.",
                "evidence_scope": "species_direct",
                "name_match_method": "exact_accepted_name",
                "source_lineage": "lineage:one",
                "lineage_method": "original_record",
                "source_run_id": "fixture",
                "source_artifact": "fixture-reviewed",
                "acceptance_contract": "reviewed_exact_species_v1",
            }
        ]
    )
    reviewed_path = _write_csv(reviewed, tmp_path / "reviewed.csv")

    all_evidence_dir = tmp_path / "audit"
    all_evidence_dir.mkdir()
    resolved = pd.DataFrame(
        [
            {
                "accepted_species": species_name,
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "resolution_status": "resolved",
                "quality": "high",
                "state_set": '["SC"]',
                "source_groups": "reviewed_database",
                "source_lineages": "lineage:one",
            }
            for species_name in ("Alpha one", "Beta one")
        ]
    )
    resolved.to_csv(
        all_evidence_dir / "resolved_direct_species_trait.csv.gz",
        index=False,
    )
    rule = {
        "setting": "current_min3",
        "eligible": "True",
        "genus": "Alpha",
        "axis": "reproductive_assurance",
        "trait_name": "self_incompatibility",
        "inferred_state_set": '["SC"]',
        "n_direct_species": "3",
        "species_loo_accuracy": "1.0",
        "lineage_loo_accuracy": "1.0",
    }
    pd.DataFrame([rule]).to_csv(
        all_evidence_dir / "trait_specific_genus_rule_audit.csv.gz",
        index=False,
    )
    pd.DataFrame(
        [
            {
                "accepted_species": "Alpha two",
                "genus": "Alpha",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "inferred_state_set": '["SC"]',
                "quality": "low",
                "family_inference_used": "False",
                "global_fallback_used": "False",
                "source_lineage": "rule:Alpha:self_incompatibility",
            }
        ]
    ).to_csv(
        all_evidence_dir / "rebuilt_all_evidence_validated_low.csv.gz",
        index=False,
    )
    prior_path = _write_csv(
        pd.DataFrame(
            columns=[
                "setting",
                "eligible",
                "genus",
                "axis",
                "trait_name",
                "inferred_state_set",
            ]
        ),
        tmp_path / "prior_rules.csv",
    )

    summary = build_source_recovery_overlay(
        baseline_csv=baseline_path,
        all_evidence_dir=all_evidence_dir,
        reviewed_direct_csvs=(reviewed_path,),
        previous_rule_audit_csvs=(prior_path,),
        output_dir=tmp_path / "output",
        expected_species=3,
    )

    assert summary["delta"]["gross_gain_species_axis"] == 2
    assert summary["delta"]["loss_species_axis"] == 0
    assert summary["delta"]["action_counts"] == {
        "direct_fill": 1,
        "validated_low_fill": 1,
    }
    assert summary["new_eligible_rules"] == ["Alpha x self_incompatibility"]
    assert all(summary["checks"].values())
