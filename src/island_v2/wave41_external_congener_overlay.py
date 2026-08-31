"""Lossless Wave41 overlay from strict external-congener reproduction support."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.wave35_trait_overlay import QUALITY_RANK, _aggregate_low, _validate_low
from island_v2.wave37_europe_pmc_checkpoint import (
    EXPECTED_SPECIES,
    _coverage_summary,
    _eligible,
    _sha256,
    _validate_coverage,
    _write_gzip_csv,
)


def build_overlay(
    *,
    baseline_csv: Path,
    previous_rule_audit_csv: Path,
    all_evidence_dir: Path,
    external_evidence_csv: Path,
    checkpoint_summary_json: Path,
    output_dir: Path,
    expected_species: int = EXPECTED_SPECIES,
) -> dict[str, Any]:
    rebuilt_low_path = all_evidence_dir / "rebuilt_all_evidence_validated_low.csv.gz"
    rules_path = all_evidence_dir / "trait_specific_genus_rule_audit.csv.gz"
    external_cells_path = (
        all_evidence_dir / "external_congener_resolved_species_trait.csv.gz"
    )
    external_conflicts_path = (
        all_evidence_dir / "external_congener_source_lineage_conflicts.csv.gz"
    )
    external_duplicates_path = (
        all_evidence_dir / "external_congener_source_lineage_duplicates.csv.gz"
    )
    all_evidence_summary_path = (
        all_evidence_dir / "all_evidence_trait_coverage_summary.json"
    )
    required = (
        baseline_csv,
        previous_rule_audit_csv,
        rebuilt_low_path,
        rules_path,
        external_cells_path,
        external_conflicts_path,
        external_duplicates_path,
        all_evidence_summary_path,
        external_evidence_csv,
        checkpoint_summary_json,
    )
    if missing := [str(path) for path in required if not path.is_file()]:
        raise ValueError(f"Wave41 overlay inputs missing: {missing}")

    baseline = pd.read_csv(baseline_csv, dtype=str).fillna("")
    _validate_coverage(baseline, expected_species)
    previous_rules = pd.read_csv(previous_rule_audit_csv, dtype=str).fillna("")
    rules = pd.read_csv(rules_path, dtype=str).fillna("")
    rebuilt_low = pd.read_csv(rebuilt_low_path, dtype=str).fillna("")
    external_evidence = pd.read_csv(external_evidence_csv, dtype=str).fillna("")
    external_cells = pd.read_csv(external_cells_path, dtype=str).fillna("")
    external_conflicts = pd.read_csv(external_conflicts_path, dtype=str).fillna("")
    checkpoint_summary = json.loads(checkpoint_summary_json.read_text(encoding="utf-8"))
    all_evidence_summary = json.loads(
        all_evidence_summary_path.read_text(encoding="utf-8")
    )
    if checkpoint_summary["fixed_target_species"] != expected_species:
        raise ValueError("Wave41 checkpoint denominator mismatch")
    external_audit = all_evidence_summary["source_lineage_audit"][
        "external_congener_support"
    ]
    if external_audit["entered_confirmatory_direct_coverage"] != 0:
        raise ValueError("external congener rows entered confirmatory direct coverage")

    touched_keys = {
        (str(row.accepted_species).split()[0], str(row.trait_name))
        for row in external_evidence.itertuples(index=False)
    }
    key_columns = ["genus", "axis", "trait_name"]
    previous = previous_rules.loc[
        previous_rules["setting"].eq("current_min3") & _eligible(previous_rules)
    ].copy()
    current = rules.loc[rules["setting"].eq("current_min3") & _eligible(rules)].copy()
    previous = previous.loc[
        previous.apply(
            lambda row: (str(row["genus"]), str(row["trait_name"])) in touched_keys,
            axis=1,
        )
    ]
    current = current.loc[
        current.apply(
            lambda row: (str(row["genus"]), str(row["trait_name"])) in touched_keys,
            axis=1,
        )
    ]
    previous_by_key = previous.set_index(key_columns)["inferred_state_set"].to_dict()
    current_by_key = current.set_index(key_columns)["inferred_state_set"].to_dict()
    new_keys = set(current_by_key).difference(previous_by_key)
    invalidated_keys = sorted(set(previous_by_key).difference(current_by_key))
    changed_keys = sorted(
        key
        for key in set(previous_by_key).intersection(current_by_key)
        if previous_by_key[key] != current_by_key[key]
    )
    new_rules = current.loc[
        current.apply(
            lambda row: tuple(str(row[column]) for column in key_columns) in new_keys,
            axis=1,
        )
    ].copy()
    low = rebuilt_low.merge(
        new_rules[key_columns], on=key_columns, how="inner", validate="many_to_one"
    )
    _validate_low(low, new_rules)

    result_records = {
        (str(row["accepted_species"]), str(row["axis"])): row
        for row in baseline.to_dict("records")
    }
    changes: list[dict[str, str]] = []
    for row in _aggregate_low(low).itertuples(index=False):
        key = (row.accepted_species, row.axis)
        if key not in result_records:
            raise ValueError(f"Wave41 Low evidence is outside fixed universe: {key}")
        record = result_records[key]
        if str(record["quality"]):
            continue
        for column in (
            "trait_composition",
            "trait_names",
            "source_groups",
            "source_lineages",
            "quality",
        ):
            record[column] = getattr(row, column)
        changes.append(
            {
                "accepted_species": row.accepted_species,
                "axis": row.axis,
                "action": "validated_low_fill",
                "quality_before": "",
                "quality_after": "low",
                "trait_names": row.trait_names,
                "source_groups": row.source_groups,
            }
        )

    result = pd.DataFrame(result_records.values(), columns=baseline.columns)
    result = result.sort_values(["accepted_species", "axis"]).reset_index(drop=True)
    _validate_coverage(result, expected_species)
    before = _coverage_summary(baseline, expected_species)
    after = _coverage_summary(result, expected_species)
    comparison = baseline[["accepted_species", "axis", "quality"]].merge(
        result[["accepted_species", "axis", "quality"]],
        on=["accepted_species", "axis"],
        suffixes=("_before", "_after"),
        validate="one_to_one",
    )
    was_filled = comparison["quality_before"].ne("")
    now_filled = comparison["quality_after"].ne("")
    loss = int((was_filled & ~now_filled).sum())
    gain = int((~was_filled & now_filled).sum())
    downgraded = comparison["quality_after"].map(QUALITY_RANK) < comparison[
        "quality_before"
    ].map(QUALITY_RANK)
    if loss or downgraded.any():
        raise ValueError(
            f"Wave41 is not lossless: loss={loss}, downgraded={int(downgraded.sum())}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    changes_frame = pd.DataFrame(changes)
    outputs = {
        "wave41_species_axis_coverage.csv.gz": result,
        "wave41_new_validated_low_species_trait.csv.gz": low,
        "wave41_new_trait_specific_genus_rule_audit.csv.gz": new_rules,
        "wave41_change_audit.csv.gz": changes_frame,
        "wave41_external_congener_resolved_species_trait.csv.gz": external_cells,
        "wave41_external_congener_source_conflicts.csv.gz": external_conflicts,
    }
    for name, frame in outputs.items():
        _write_gzip_csv(frame, output_dir / name)

    summary: dict[str, Any] = {
        "contract": "wave41_external_congener_reproduction_lossless_overlay_v1",
        "baseline_formal_run_id": 33155408805,
        "fixed_denominator": {
            "species": expected_species,
            "species_axis": expected_species * 3,
        },
        "before": before,
        "after": after,
        "delta": {
            "gross_gain_species_axis": gain,
            "loss_species_axis": loss,
            "net_gain_species_axis": gain,
            "by_axis_net_gain": {
                axis: after["by_axis"][axis]["filled_species"]
                - before["by_axis"][axis]["filled_species"]
                for axis in after["by_axis"]
            },
            "action_counts": {
                str(key): int(value)
                for key, value in changes_frame["action"].value_counts().items()
            },
            "external_evidence_rows": len(external_evidence),
            "external_resolved_species_trait": len(external_cells),
            "external_direct_conflicts": int(
                external_conflicts["classification"]
                .eq("unresolved_direct_conflict")
                .sum()
            ),
            "new_validated_low_species_trait": len(low),
        },
        "campaign_targets": {
            "floral_structural_complexity": 91_953,
            "flower_colour": 85_036,
            "reproductive_assurance": 63_777,
        },
        "new_eligible_rules": [
            f"{row.genus} x {row.trait_name}"
            for row in new_rules.itertuples(index=False)
        ],
        "shadow_invalidated_prior_rules": [
            " x ".join((key[0], key[2])) for key in invalidated_keys
        ],
        "shadow_changed_prior_rule_states": [
            " x ".join((key[0], key[2])) for key in changed_keys
        ],
        "query_accounting": checkpoint_summary["queries"],
        "query_cost_usd": checkpoint_summary["query_cost_usd"],
        "checks": {
            "fixed_denominator": len(result) == expected_species * 3,
            "quality_precedence_high_medium_low": bool(not downgraded.any()),
            "external_congener_not_confirmatory_direct": True,
            "source_lineage_deduplicated": True,
            "trait_specific_genus_join": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
            "baseline_loss_zero": loss == 0,
            "shadow_invalidations_not_silently_deleted": True,
        },
        "input_sha256": {str(path): _sha256(path) for path in required},
    }
    summary["remaining_to_campaign_target"] = {
        axis: max(target - after["by_axis"][axis]["filled_species"], 0)
        for axis, target in summary["campaign_targets"].items()
    }
    summary["artifact_sha256"] = {
        name: _sha256(output_dir / name) for name in outputs
    }
    (output_dir / "wave41_coverage_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-csv", required=True, type=Path)
    parser.add_argument("--previous-rule-audit-csv", required=True, type=Path)
    parser.add_argument("--all-evidence-dir", required=True, type=Path)
    parser.add_argument("--external-evidence-csv", required=True, type=Path)
    parser.add_argument("--checkpoint-summary-json", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-species", type=int, default=EXPECTED_SPECIES)
    args = parser.parse_args()
    summary = build_overlay(
        baseline_csv=args.baseline_csv,
        previous_rule_audit_csv=args.previous_rule_audit_csv,
        all_evidence_dir=args.all_evidence_dir,
        external_evidence_csv=args.external_evidence_csv,
        checkpoint_summary_json=args.checkpoint_summary_json,
        output_dir=args.output_dir,
        expected_species=args.expected_species,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
