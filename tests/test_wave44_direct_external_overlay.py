from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from island_v2.wave44_direct_external_overlay import build_wave44_overlay

AXES = (
    "floral_structural_complexity",
    "flower_colour",
    "reproductive_assurance",
)


def _write(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, compression="gzip")


def test_wave44_overlay_adds_new_trait_specific_low_and_upgrades_direct(
    tmp_path: Path,
) -> None:
    species = ("A alpha", "B beta", "B gamma")
    baseline_rows = []
    for accepted_species in species:
        for axis in AXES:
            low = accepted_species == "A alpha" and axis == "reproductive_assurance"
            baseline_rows.append(
                {
                    "accepted_species": accepted_species,
                    "axis": axis,
                    "trait_composition": (
                        'self_incompatibility=["SC"]' if low else ""
                    ),
                    "trait_names": "self_incompatibility" if low else "",
                    "source_groups": "prior_low" if low else "",
                    "source_lineages": "prior:A:self_incompatibility" if low else "",
                    "quality": "low" if low else "",
                }
            )
    baseline = tmp_path / "baseline.csv.gz"
    _write(pd.DataFrame(baseline_rows), baseline)

    rule = {
        "setting": "current_min3",
        "genus": "B",
        "axis": "reproductive_assurance",
        "trait_name": "self_incompatibility",
        "inferred_state_set": '["SC"]',
        "eligible": "True",
        "n_direct_species": "3",
        "species_loo_accuracy": "1.0",
        "lineage_loo_accuracy": "1.0",
    }
    previous_rules = tmp_path / "previous-rules.csv.gz"
    _write(pd.DataFrame([{**rule, "eligible": "False"}]), previous_rules)
    audit_dir = tmp_path / "audit"
    audit_dir.mkdir()
    _write(pd.DataFrame([rule]), audit_dir / "trait_specific_genus_rule_audit.csv.gz")

    resolved = pd.DataFrame(
        [
            {
                "accepted_species": "A alpha",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "resolution_status": "resolved",
                "quality": "high",
                "state_set": '["SC"]',
                "source_groups": "wave44_glopl_2025_reproduction",
                "source_lineages": "glopl:A:self_incompatibility",
            }
        ]
    )
    _write(resolved, audit_dir / "resolved_direct_species_trait.csv.gz")
    low_rows = []
    for accepted_species in ("B beta", "B gamma"):
        low_rows.append(
            {
                **rule,
                "accepted_species": accepted_species,
                "quality": "low",
                "evidence_quality": "low",
                "evidence_scope": "genus_consensus",
                "inference_method": "validated_genus_consensus",
                "state_set": '["SC"]',
                "normalized_value": "SC",
                "family_inference_used": "False",
                "global_fallback_used": "False",
                "source_lineage": "validated-low:B:self_incompatibility",
            }
        )
    _write(
        pd.DataFrame(low_rows),
        audit_dir / "rebuilt_all_evidence_validated_low.csv.gz",
    )
    _write(
        pd.DataFrame([{"accepted_species": "B delta", "trait_name": "self_incompatibility"}]),
        audit_dir / "external_congener_resolved_species_trait.csv.gz",
    )
    _write(
        pd.DataFrame(columns=["classification"]),
        audit_dir / "external_congener_source_lineage_conflicts.csv.gz",
    )
    (audit_dir / "all_evidence_trait_coverage_summary.json").write_text(
        json.dumps(
            {
                "source_lineage_audit": {
                    "external_congener_support": {
                        "entered_confirmatory_direct_coverage": 0
                    }
                }
            }
        ),
        encoding="utf-8",
    )

    direct_csv = tmp_path / "direct.csv.gz"
    external_csv = tmp_path / "external.csv.gz"
    _write(
        pd.DataFrame(
            [{"accepted_species": "A alpha", "trait_name": "self_incompatibility"}]
        ),
        direct_csv,
    )
    _write(
        pd.DataFrame(
            [{"accepted_species": "B delta", "trait_name": "self_incompatibility"}]
        ),
        external_csv,
    )
    checkpoint = tmp_path / "checkpoint.json"
    checkpoint.write_text(
        json.dumps(
            {
                "fixed_target_species": 3,
                "evidence": {"external_species_trait": 1},
                "queries": {"query_cost_usd": 0},
                "query_cost_usd": 0,
            }
        ),
        encoding="utf-8",
    )

    summary = build_wave44_overlay(
        baseline_csv=baseline,
        previous_rule_audit_csv=previous_rules,
        all_evidence_dir=audit_dir,
        direct_evidence_csv=direct_csv,
        external_evidence_csv=external_csv,
        checkpoint_summary_json=checkpoint,
        output_dir=tmp_path / "output",
        expected_species=3,
    )
    assert summary["delta"]["net_gain_species_axis"] == 2
    assert summary["delta"]["action_counts"] == {
        "validated_low_fill": 2,
        "direct_upgrade": 1,
    }
    assert summary["delta"]["new_eligible_genus_trait_rules"] == 1
    assert summary["checks"]["trait_specific_genus_join"] is True
