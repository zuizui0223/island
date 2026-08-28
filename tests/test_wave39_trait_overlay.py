from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from island_v2.floraweb_synonym_recovery import SOURCE_GROUP
from island_v2.wave35_trait_overlay import AXES
from island_v2.wave39_trait_overlay import build_wave39_lossless_overlay


def test_formal_overlay_keeps_shadow_invalidated_low_and_adds_direct(
    tmp_path: Path,
) -> None:
    species = ["Genus one", "Genus two"]
    baseline = pd.DataFrame(
        [(name, axis, "", "", "", "", "") for name in species for axis in AXES],
        columns=[
            "accepted_species",
            "axis",
            "trait_composition",
            "trait_names",
            "source_groups",
            "source_lineages",
            "quality",
        ],
    )
    baseline.loc[
        (baseline["accepted_species"] == "Genus one")
        & (baseline["axis"] == "reproductive_assurance"),
        ["trait_composition", "trait_names", "source_groups", "source_lineages", "quality"],
    ] = [
        'self_incompatibility=["SI"]',
        "self_incompatibility",
        "prior_low",
        "validated-low:Genus:self_incompatibility",
        "low",
    ]
    baseline_path = tmp_path / "baseline.csv.gz"
    baseline.to_csv(baseline_path, index=False)

    reviewed = pd.DataFrame(
        [{"accepted_species": "Genus two", "trait_name": "self_incompatibility"}]
    )
    reviewed_path = tmp_path / "reviewed.csv.gz"
    reviewed.to_csv(reviewed_path, index=False)
    checkpoint_path = tmp_path / "checkpoint.json"
    checkpoint_path.write_text("{}\n", encoding="utf-8")

    rebuild = tmp_path / "rebuild"
    rebuild.mkdir()
    pd.DataFrame(
        [
            {
                "accepted_species": "Genus two",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "resolution_status": "resolved",
                "quality": "high",
                "state_set": '["SC"]',
                "source_groups": SOURCE_GROUP,
                "source_lineages": "biolflor:test",
            }
        ]
    ).to_csv(rebuild / "resolved_direct_species_trait.csv.gz", index=False)
    low_columns = [
        "accepted_species",
        "genus",
        "axis",
        "trait_name",
        "inferred_state_set",
        "quality",
        "family_inference_used",
        "global_fallback_used",
        "source_lineage",
    ]
    pd.DataFrame(columns=low_columns).to_csv(
        rebuild / "rebuilt_all_evidence_validated_low.csv.gz", index=False
    )
    rule_columns = [
        "setting",
        "genus",
        "axis",
        "trait_name",
        "inferred_state_set",
        "eligible",
        "n_direct_species",
        "species_loo_accuracy",
        "lineage_loo_accuracy",
    ]
    pd.DataFrame(
        [
            {
                "setting": "current_min3",
                "genus": "Genus",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "inferred_state_set": '["SC"]',
                "eligible": "False",
                "n_direct_species": "3",
                "species_loo_accuracy": "0.5",
                "lineage_loo_accuracy": "0.5",
            }
        ],
        columns=rule_columns,
    ).to_csv(rebuild / "trait_specific_genus_rule_audit.csv.gz", index=False)

    prior_rules = tmp_path / "prior_rules.csv.gz"
    pd.DataFrame(
        [
            {
                "setting": "current_min3",
                "genus": "Genus",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "inferred_state_set": '["SI"]',
                "eligible": "True",
                "n_direct_species": "3",
                "species_loo_accuracy": "1.0",
                "lineage_loo_accuracy": "1.0",
            }
        ],
        columns=rule_columns,
    ).to_csv(prior_rules, index=False)

    output = tmp_path / "output"
    summary = build_wave39_lossless_overlay(
        baseline_csv=baseline_path,
        all_evidence_dir=rebuild,
        reviewed_direct_csv=reviewed_path,
        checkpoint_summary_json=checkpoint_path,
        previous_rule_audit_csvs=(prior_rules,),
        output_dir=output,
        expected_species=2,
    )

    assert summary["delta"]["gross_gain_species_axis"] == 1
    assert summary["delta"]["loss_species_axis"] == 0
    assert summary["shadow_invalidated_prior_rules"] == [
        "Genus x self_incompatibility"
    ]
    result = pd.read_csv(output / "wave39_species_axis_coverage.csv.gz").fillna("")
    old_low = result.loc[
        (result["accepted_species"] == "Genus one")
        & (result["axis"] == "reproductive_assurance")
    ].iloc[0]
    new_direct = result.loc[
        (result["accepted_species"] == "Genus two")
        & (result["axis"] == "reproductive_assurance")
    ].iloc[0]
    assert old_low["quality"] == "low"
    assert new_direct["quality"] == "high"
    persisted = json.loads((output / "wave39_coverage_summary.json").read_text())
    assert persisted["checks"]["shadow_invalidations_not_silently_deleted"]
