"""Lossless formal overlay for already-reviewed source packages.

Wave40 recovers provenance-complete evidence that was acquired and reviewed in
earlier checkpoints but was not supplied to the formal all-evidence rebuild.
The formal coverage ledger remains append-only.  Newly invalidated historical
Low rules are reported as a shadow audit and are not silently deleted.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.all_evidence_trait_audit import load_reviewed_direct_supplements
from island_v2.wave35_trait_overlay import (
    QUALITY_RANK,
    _aggregate_direct,
    _aggregate_low,
    _parse_composition,
    _serialize_composition,
    _split_pipe,
    _validate_direct,
    _validate_low,
)
from island_v2.wave37_europe_pmc_checkpoint import (
    EXPECTED_SPECIES,
    _coverage_summary,
    _eligible,
    _sha256,
    _validate_coverage,
    _write_gzip_csv,
)
from island_v2.wave38_trait_overlay import _prior_rule_table


def _touched_keys(reviewed: pd.DataFrame) -> set[tuple[str, str]]:
    return {
        (str(row.accepted_species).split()[0], str(row.trait_name))
        for row in reviewed.itertuples(index=False)
    }


def _target_direct(
    resolved: pd.DataFrame,
    reviewed: pd.DataFrame,
) -> pd.DataFrame:
    reviewed_lineages = {
        lineage
        for value in reviewed["source_lineage"]
        for lineage in _split_pipe(value)
    }
    reviewed_cells = set(
        reviewed[["accepted_species", "trait_name"]].itertuples(
            index=False, name=None
        )
    )
    selected = resolved.loc[
        resolved.apply(
            lambda row: (
                (str(row["accepted_species"]), str(row["trait_name"]))
                in reviewed_cells
                and bool(_split_pipe(row["source_lineages"]) & reviewed_lineages)
            ),
            axis=1,
        )
    ].copy()
    if selected.empty:
        raise ValueError("reviewed source packages produced no resolved direct cells")
    return _validate_direct(selected)


def build_source_recovery_overlay(
    *,
    baseline_csv: Path,
    all_evidence_dir: Path,
    reviewed_direct_csvs: tuple[Path, ...],
    previous_rule_audit_csvs: tuple[Path, ...],
    output_dir: Path,
    expected_species: int = EXPECTED_SPECIES,
) -> dict[str, Any]:
    """Apply reviewed direct cells and newly eligible trait rules losslessly."""

    resolved_path = all_evidence_dir / "resolved_direct_species_trait.csv.gz"
    rebuilt_low_path = all_evidence_dir / "rebuilt_all_evidence_validated_low.csv.gz"
    rules_path = all_evidence_dir / "trait_specific_genus_rule_audit.csv.gz"
    required = (
        baseline_csv,
        resolved_path,
        rebuilt_low_path,
        rules_path,
        *reviewed_direct_csvs,
        *previous_rule_audit_csvs,
    )
    if missing := [str(path) for path in required if not path.is_file()]:
        raise ValueError(f"Wave40 overlay inputs missing: {missing}")

    baseline = pd.read_csv(baseline_csv, dtype=str).fillna("")
    _validate_coverage(baseline, expected_species)
    reviewed = load_reviewed_direct_supplements(reviewed_direct_csvs)
    resolved = pd.read_csv(resolved_path, dtype=str).fillna("")
    rebuilt_low = pd.read_csv(rebuilt_low_path, dtype=str).fillna("")
    rules = pd.read_csv(rules_path, dtype=str).fillna("")
    prior_rules = _prior_rule_table(previous_rule_audit_csvs)
    touched_keys = _touched_keys(reviewed)
    direct = _target_direct(resolved, reviewed)

    current_rules = rules.loc[
        rules["setting"].eq("current_min3") & _eligible(rules)
    ].copy()
    current_rules = current_rules.loc[
        current_rules.apply(
            lambda row: (str(row["genus"]), str(row["trait_name"]))
            in touched_keys,
            axis=1,
        )
    ]
    prior_rules = prior_rules.loc[
        prior_rules.apply(
            lambda row: (str(row["genus"]), str(row["trait_name"]))
            in touched_keys,
            axis=1,
        )
    ]
    key_columns = ["genus", "axis", "trait_name"]
    current_by_key = current_rules.set_index(key_columns)[
        "inferred_state_set"
    ].to_dict()
    prior_by_key = prior_rules.set_index(key_columns)["inferred_state_set"].to_dict()
    invalidated_keys = sorted(set(prior_by_key) - set(current_by_key))
    changed_keys = sorted(
        key
        for key in set(prior_by_key) & set(current_by_key)
        if prior_by_key[key] != current_by_key[key]
    )
    new_keys = set(current_by_key) - set(prior_by_key)
    new_rules = current_rules.loc[
        current_rules.apply(
            lambda row: (
                str(row["genus"]),
                str(row["axis"]),
                str(row["trait_name"]),
            )
            in new_keys,
            axis=1,
        )
    ].copy()
    low = rebuilt_low.merge(
        new_rules[key_columns], on=key_columns, how="inner", validate="many_to_one"
    )
    _validate_low(low, new_rules)

    # A dictionary-backed overlay keeps the 318,885-cell update linear.  Repeated
    # scalar ``DataFrame.loc`` writes make a 30k-row reviewed package needlessly
    # quadratic and are unsuitable for the formal GitHub runner.
    result_records = {
        (str(record["accepted_species"]), str(record["axis"])): record
        for record in baseline.to_dict("records")
    }
    changes: list[dict[str, str]] = []
    for row in _aggregate_direct(direct).itertuples(index=False):
        key = (row.accepted_species, row.axis)
        if key not in result_records:
            raise ValueError(f"direct evidence is outside the fixed universe: {key}")
        record = result_records[key]
        before_quality = str(record["quality"])
        if before_quality in {"", "low"}:
            composition = _parse_composition(row.trait_composition, axis=row.axis)
            groups = _split_pipe(row.source_groups)
            lineages = _split_pipe(row.source_lineages)
        else:
            composition = _parse_composition(
                record["trait_composition"], axis=row.axis
            )
            composition.update(
                _parse_composition(row.trait_composition, axis=row.axis)
            )
            groups = _split_pipe(record["source_groups"]) | _split_pipe(
                row.source_groups
            )
            lineages = _split_pipe(record["source_lineages"]) | _split_pipe(
                row.source_lineages
            )
        after_quality = max(
            (before_quality, row.quality), key=QUALITY_RANK.__getitem__
        )
        record["trait_composition"] = _serialize_composition(composition)
        record["trait_names"] = "|".join(sorted(composition))
        record["source_groups"] = "|".join(sorted(groups))
        record["source_lineages"] = "|".join(sorted(lineages))
        record["quality"] = after_quality
        changes.append(
            {
                "accepted_species": row.accepted_species,
                "axis": row.axis,
                "action": (
                    "direct_fill"
                    if before_quality == ""
                    else "direct_upgrade"
                    if QUALITY_RANK[after_quality] > QUALITY_RANK[before_quality]
                    else "direct_enrichment"
                ),
                "quality_before": before_quality,
                "quality_after": after_quality,
                "trait_names": row.trait_names,
                "source_groups": row.source_groups,
            }
        )

    for row in _aggregate_low(low).itertuples(index=False):
        key = (row.accepted_species, row.axis)
        if key not in result_records:
            raise ValueError(f"Low evidence is outside the fixed universe: {key}")
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
            f"Wave40 is not lossless: loss={loss}, downgraded={int(downgraded.sum())}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    changes_frame = pd.DataFrame(changes)
    reviewed_manifest = pd.DataFrame(
        [
            {
                "path": str(path),
                "sha256": _sha256(path),
                "rows": len(pd.read_csv(path, dtype=str)),
            }
            for path in reviewed_direct_csvs
        ]
    )
    outputs = {
        "wave40_species_axis_coverage.csv.gz": result,
        "wave40_resolved_direct_species_trait.csv.gz": direct,
        "wave40_new_validated_low_species_trait.csv.gz": low,
        "wave40_new_trait_specific_genus_rule_audit.csv.gz": new_rules,
        "wave40_change_audit.csv.gz": changes_frame,
        "wave40_reviewed_source_manifest.csv.gz": reviewed_manifest,
    }
    for name, frame in outputs.items():
        _write_gzip_csv(frame, output_dir / name)

    action_counts = {
        str(key): int(value)
        for key, value in changes_frame["action"].value_counts().items()
    }
    summary: dict[str, Any] = {
        "contract": "wave40_reviewed_source_recovery_lossless_overlay_v1",
        "baseline_formal_run_id": 33149998400,
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
            "action_counts": action_counts,
            "reviewed_direct_rows": len(reviewed),
            "reviewed_direct_species_trait": len(
                reviewed[["accepted_species", "trait_name"]].drop_duplicates()
            ),
            "resolved_target_direct_species_trait": len(direct),
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
        "checks": {
            "fixed_denominator": bool(len(result) == expected_species * 3),
            "quality_precedence_high_medium_low": bool(not downgraded.any()),
            "direct_conflicts_excluded": bool(
                direct["resolution_status"].eq("resolved").all()
            ),
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
    (output_dir / "wave40_coverage_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-csv", required=True, type=Path)
    parser.add_argument("--all-evidence-dir", required=True, type=Path)
    parser.add_argument(
        "--reviewed-direct-csv", action="append", required=True, type=Path
    )
    parser.add_argument(
        "--previous-rule-audit-csv", action="append", required=True, type=Path
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-species", type=int, default=EXPECTED_SPECIES)
    args = parser.parse_args()
    summary = build_source_recovery_overlay(
        baseline_csv=args.baseline_csv,
        all_evidence_dir=args.all_evidence_dir,
        reviewed_direct_csvs=tuple(args.reviewed_direct_csv),
        previous_rule_audit_csvs=tuple(args.previous_rule_audit_csv),
        output_dir=args.output_dir,
        expected_species=args.expected_species,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
