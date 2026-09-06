import json
from pathlib import Path

import pandas as pd

from island_v2.wave35_trait_overlay import AXES
from island_v2.wave41_external_congener_overlay import build_overlay


def _write_csv(frame: pd.DataFrame, path: Path) -> Path:
    frame.to_csv(path, index=False)
    return path


def test_external_congener_overlay_is_lossless_and_trait_specific(
    tmp_path: Path,
) -> None:
    species = ["Alpha target", "Beta one", "Gamma one"]
    baseline = pd.DataFrame(
        [
            {
                "accepted_species": name,
                "axis": axis,
                "trait_composition": (
                    'self_incompatibility=["SC"]'
                    if name == "Beta one" and axis == "reproductive_assurance"
                    else ""
                ),
                "trait_names": (
                    "self_incompatibility"
                    if name == "Beta one" and axis == "reproductive_assurance"
                    else ""
                ),
                "source_groups": (
                    "baseline" if name == "Beta one" and axis == "reproductive_assurance" else ""
                ),
                "source_lineages": (
                    "paper:beta" if name == "Beta one" and axis == "reproductive_assurance" else ""
                ),
                "quality": (
                    "high" if name == "Beta one" and axis == "reproductive_assurance" else ""
                ),
            }
            for name in species
            for axis in AXES
        ]
    )
    baseline_path = _write_csv(baseline, tmp_path / "baseline.csv")

    rule_columns = [
        "setting",
        "eligible",
        "genus",
        "axis",
        "trait_name",
        "inferred_state_set",
        "n_direct_species",
        "species_loo_accuracy",
        "lineage_loo_accuracy",
    ]
    previous_rule_path = _write_csv(
        pd.DataFrame(
            [
                {
                    "setting": "current_min3",
                    "eligible": "True",
                    "genus": "Beta",
                    "axis": "reproductive_assurance",
                    "trait_name": "self_incompatibility",
                    "inferred_state_set": '["SC"]',
                    "n_direct_species": "3",
                    "species_loo_accuracy": "1.0",
                    "lineage_loo_accuracy": "1.0",
                }
            ],
            columns=rule_columns,
        ),
        tmp_path / "previous-rules.csv",
    )

    audit_dir = tmp_path / "audit"
    audit_dir.mkdir()
    pd.DataFrame(
        [
            {
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
        ],
        columns=rule_columns,
    ).to_csv(audit_dir / "trait_specific_genus_rule_audit.csv.gz", index=False)
    pd.DataFrame(
        [
            {
                "accepted_species": "Alpha target",
                "genus": "Alpha",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "inferred_state_set": '["SC"]',
                "quality": "low",
                "family_inference_used": "False",
                "global_fallback_used": "False",
                "source_lineage": "genus-rule:Alpha:self_incompatibility",
            }
        ]
    ).to_csv(
        audit_dir / "rebuilt_all_evidence_validated_low.csv.gz", index=False
    )
    pd.DataFrame(
        [
            {
                "accepted_species": "Alpha outside",
                "trait_name": "self_incompatibility",
                "classification": "single_independent_lineage",
            }
        ]
    ).to_csv(
        audit_dir / "external_congener_resolved_species_trait.csv.gz", index=False
    )
    pd.DataFrame(
        [
            {
                "accepted_species": "Beta outside",
                "trait_name": "self_incompatibility",
                "classification": "unresolved_direct_conflict",
            }
        ]
    ).to_csv(
        audit_dir / "external_congener_source_lineage_conflicts.csv.gz", index=False
    )
    pd.DataFrame(columns=["accepted_species", "trait_name"]).to_csv(
        audit_dir / "external_congener_source_lineage_duplicates.csv.gz", index=False
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

    external_path = _write_csv(
        pd.DataFrame(
            [
                {
                    "accepted_species": "Alpha outside",
                    "trait_name": "self_incompatibility",
                },
                {
                    "accepted_species": "Beta outside",
                    "trait_name": "self_incompatibility",
                },
            ]
        ),
        tmp_path / "external.csv",
    )
    checkpoint_path = tmp_path / "checkpoint.json"
    checkpoint_path.write_text(
        json.dumps(
            {
                "fixed_target_species": 3,
                "queries": {"wfo": 2, "gbif": 2, "total": 4},
                "query_cost_usd": 0,
            }
        ),
        encoding="utf-8",
    )

    summary = build_overlay(
        baseline_csv=baseline_path,
        previous_rule_audit_csv=previous_rule_path,
        all_evidence_dir=audit_dir,
        external_evidence_csv=external_path,
        checkpoint_summary_json=checkpoint_path,
        output_dir=tmp_path / "output",
        expected_species=3,
    )

    assert summary["delta"]["gross_gain_species_axis"] == 1
    assert summary["delta"]["loss_species_axis"] == 0
    assert summary["delta"]["action_counts"] == {"validated_low_fill": 1}
    assert summary["delta"]["external_direct_conflicts"] == 1
    assert summary["new_eligible_rules"] == ["Alpha x self_incompatibility"]
    assert summary["shadow_invalidated_prior_rules"] == [
        "Beta x self_incompatibility"
    ]
    assert summary["shadow_changed_prior_rule_states"] == []
    assert summary["after"]["quality_counts"] == {
        "high": 1,
        "medium": 0,
        "low": 1,
    }
    result = pd.read_csv(
        tmp_path / "output" / "wave41_species_axis_coverage.csv.gz",
        dtype=str,
    ).fillna("")
    assert "Alpha outside" not in set(result["accepted_species"])
    beta = result.loc[
        result["accepted_species"].eq("Beta one")
        & result["axis"].eq("reproductive_assurance")
    ].iloc[0]
    assert beta["quality"] == "high"
    assert beta["source_lineages"] == "paper:beta"
    assert all(summary["checks"].values())
