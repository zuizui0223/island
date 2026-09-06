"""Build the lossless Wave45 reviewed direct plus external-congener overlay."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

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

BASELINE_FORMAL_RUN_ID = 33_362_634_097


def _touched_keys(frame: pd.DataFrame) -> set[tuple[str, str]]:
    return {
        (str(row.accepted_species).split()[0], str(row.trait_name))
        for row in frame.itertuples(index=False)
    }


def _packet_direct_rows(
    resolved: pd.DataFrame,
    reviewed_direct: pd.DataFrame,
) -> pd.DataFrame:
    """Keep only resolved cells backed by a lineage in the reviewed packet."""
    packet_lineages: dict[tuple[str, str, str], set[str]] = {}
    for row in reviewed_direct.itertuples(index=False):
        key = (str(row.accepted_species), str(row.axis), str(row.trait_name))
        packet_lineages.setdefault(key, set()).add(str(row.source_lineage))

    def selected(row: pd.Series) -> bool:
        key = (
            str(row["accepted_species"]),
            str(row["axis"]),
            str(row["trait_name"]),
        )
        return bool(packet_lineages.get(key, set()) & _split_pipe(row["source_lineages"]))

    return resolved.loc[resolved.apply(selected, axis=1)].copy()


def build_wave45_overlay(
    *,
    baseline_csv: Path,
    previous_rule_audit_csv: Path,
    all_evidence_dir: Path,
    direct_evidence_csv: Path,
    external_evidence_csv: Path,
    checkpoint_summary_json: Path,
    output_dir: Path,
    expected_species: int = EXPECTED_SPECIES,
) -> dict[str, Any]:
    resolved_path = all_evidence_dir / "resolved_direct_species_trait.csv.gz"
    low_path = all_evidence_dir / "rebuilt_all_evidence_validated_low.csv.gz"
    rules_path = all_evidence_dir / "trait_specific_genus_rule_audit.csv.gz"
    external_cells_path = (
        all_evidence_dir / "external_congener_resolved_species_trait.csv.gz"
    )
    external_conflicts_path = (
        all_evidence_dir / "external_congener_source_lineage_conflicts.csv.gz"
    )
    audit_summary_path = all_evidence_dir / "all_evidence_trait_coverage_summary.json"
    required = (
        baseline_csv,
        previous_rule_audit_csv,
        resolved_path,
        low_path,
        rules_path,
        external_cells_path,
        external_conflicts_path,
        audit_summary_path,
        direct_evidence_csv,
        external_evidence_csv,
        checkpoint_summary_json,
    )
    if missing := [str(path) for path in required if not path.is_file()]:
        raise ValueError(f"Wave45 overlay inputs missing: {missing}")

    baseline = pd.read_csv(baseline_csv, dtype=str).fillna("")
    _validate_coverage(baseline, expected_species)
    previous_rules = pd.read_csv(previous_rule_audit_csv, dtype=str).fillna("")
    rules = pd.read_csv(rules_path, dtype=str).fillna("")
    rebuilt_low = pd.read_csv(low_path, dtype=str).fillna("")
    resolved = pd.read_csv(resolved_path, dtype=str).fillna("")
    reviewed_direct = pd.read_csv(direct_evidence_csv, dtype=str).fillna("")
    external = pd.read_csv(external_evidence_csv, dtype=str).fillna("")
    external_cells = pd.read_csv(external_cells_path, dtype=str).fillna("")
    external_conflicts = pd.read_csv(external_conflicts_path, dtype=str).fillna("")
    checkpoint = json.loads(checkpoint_summary_json.read_text(encoding="utf-8"))
    audit_summary = json.loads(audit_summary_path.read_text(encoding="utf-8"))
    if checkpoint["fixed_target_species"] != expected_species:
        raise ValueError("Wave45 checkpoint denominator mismatch")
    external_audit = audit_summary["source_lineage_audit"][
        "external_congener_support"
    ]
    if external_audit["entered_confirmatory_direct_coverage"] != 0:
        raise ValueError("external congener evidence entered direct coverage")

    direct = _validate_direct(_packet_direct_rows(resolved, reviewed_direct))
    touched = _touched_keys(reviewed_direct) | _touched_keys(external)

    key_columns = ["genus", "axis", "trait_name"]
    previous = previous_rules.loc[
        previous_rules["setting"].eq("current_min3") & _eligible(previous_rules)
    ].copy()
    current = rules.loc[
        rules["setting"].eq("current_min3") & _eligible(rules)
    ].copy()
    previous = previous.loc[
        previous.apply(
            lambda row: (str(row["genus"]), str(row["trait_name"])) in touched,
            axis=1,
        )
    ]
    current = current.loc[
        current.apply(
            lambda row: (str(row["genus"]), str(row["trait_name"])) in touched,
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

    result = baseline.copy().set_index(["accepted_species", "axis"])
    changes: list[dict[str, str]] = []
    for row in _aggregate_direct(direct).itertuples(index=False):
        key = (row.accepted_species, row.axis)
        if key not in result.index:
            raise ValueError(f"Wave45 direct evidence is outside the universe: {key}")
        before_quality = str(result.loc[key, "quality"])
        if before_quality in {"", "low"}:
            composition = _parse_composition(row.trait_composition, axis=row.axis)
            groups = _split_pipe(row.source_groups)
            lineages = _split_pipe(row.source_lineages)
        else:
            composition = _parse_composition(
                result.loc[key, "trait_composition"], axis=row.axis
            )
            composition.update(_parse_composition(row.trait_composition, axis=row.axis))
            groups = _split_pipe(result.loc[key, "source_groups"]) | _split_pipe(
                row.source_groups
            )
            lineages = _split_pipe(
                result.loc[key, "source_lineages"]
            ) | _split_pipe(row.source_lineages)
        after_quality = max(
            (before_quality, row.quality), key=QUALITY_RANK.__getitem__
        )
        result.loc[key, "trait_composition"] = _serialize_composition(composition)
        result.loc[key, "trait_names"] = "|".join(sorted(composition))
        result.loc[key, "source_groups"] = "|".join(sorted(groups))
        result.loc[key, "source_lineages"] = "|".join(sorted(lineages))
        result.loc[key, "quality"] = after_quality
        action = (
            "direct_fill"
            if before_quality == ""
            else "direct_upgrade"
            if QUALITY_RANK[after_quality] > QUALITY_RANK[before_quality]
            else "direct_enrichment"
        )
        changes.append(
            {
                "accepted_species": row.accepted_species,
                "axis": row.axis,
                "action": action,
                "quality_before": before_quality,
                "quality_after": after_quality,
                "trait_names": row.trait_names,
            }
        )

    for row in _aggregate_low(low).itertuples(index=False):
        key = (row.accepted_species, row.axis)
        if key not in result.index:
            raise ValueError(f"Wave45 Low is outside the fixed universe: {key}")
        if str(result.loc[key, "quality"]):
            continue
        for column in (
            "trait_composition",
            "trait_names",
            "source_groups",
            "source_lineages",
            "quality",
        ):
            result.loc[key, column] = getattr(row, column)
        changes.append(
            {
                "accepted_species": row.accepted_species,
                "axis": row.axis,
                "action": "validated_low_fill",
                "quality_before": "",
                "quality_after": "low",
                "trait_names": row.trait_names,
            }
        )

    result = result.reset_index().sort_values(["accepted_species", "axis"])
    result = result.reset_index(drop=True)
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
            f"Wave45 is not lossless: loss={loss}, downgraded={int(downgraded.sum())}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    changes_frame = pd.DataFrame(changes)
    outputs = {
        "wave45_species_axis_coverage.csv.gz": result,
        "wave45_resolved_direct_species_trait.csv.gz": direct,
        "wave45_new_validated_low_species_trait.csv.gz": low,
        "wave45_new_trait_specific_genus_rule_audit.csv.gz": new_rules,
        "wave45_change_audit.csv.gz": changes_frame,
        "wave45_external_congener_resolved_species_trait.csv.gz": external_cells,
        "wave45_external_congener_source_conflicts.csv.gz": external_conflicts,
    }
    for name, frame in outputs.items():
        _write_gzip_csv(frame, output_dir / name)

    action_counts = (
        {
            str(key): int(value)
            for key, value in changes_frame["action"].value_counts().items()
        }
        if not changes_frame.empty
        else {}
    )
    summary: dict[str, Any] = {
        "contract": "wave45_promoted_reproductive_lossless_overlay_v1",
        "baseline_formal_run_id": BASELINE_FORMAL_RUN_ID,
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
            "reviewed_direct_rows": len(reviewed_direct),
            "resolved_wave45_direct_species_trait": len(direct),
            "resolved_wave45_direct_species_axis": int(
                direct[["accepted_species", "axis"]].drop_duplicates().shape[0]
            ),
            "new_external_input_species_trait": checkpoint["evidence"][
                "new_external_species_trait"
            ],
            "new_validated_low_species_trait": len(low),
            "new_eligible_genus_trait_rules": len(new_rules),
            "external_direct_conflicts": int(
                external_conflicts["classification"]
                .eq("unresolved_direct_conflict")
                .sum()
            ),
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
        "query_accounting": checkpoint["queries"],
        "query_cost_usd": checkpoint["query_cost_usd"],
        "checks": {
            "fixed_denominator": len(result) == expected_species * 3,
            "quality_precedence_high_medium_low": bool(not downgraded.any()),
            "direct_conflicts_excluded": True,
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
    (output_dir / "wave45_coverage_summary.json").write_text(
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
    parser.add_argument("--direct-evidence-csv", required=True, type=Path)
    parser.add_argument("--external-evidence-csv", required=True, type=Path)
    parser.add_argument("--checkpoint-summary-json", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-species", type=int, default=EXPECTED_SPECIES)
    args = parser.parse_args()
    summary = build_wave45_overlay(
        baseline_csv=args.baseline_csv,
        previous_rule_audit_csv=args.previous_rule_audit_csv,
        all_evidence_dir=args.all_evidence_dir,
        direct_evidence_csv=args.direct_evidence_csv,
        external_evidence_csv=args.external_evidence_csv,
        checkpoint_summary_json=args.checkpoint_summary_json,
        output_dir=args.output_dir,
        expected_species=args.expected_species,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
