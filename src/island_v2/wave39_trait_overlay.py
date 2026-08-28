"""Scientifically rebased Wave39 direct and trait-specific Low overlay.

Unlike the earlier strictly additive overlays, this module does not leave a
stale Low rule in place when newly integrated direct evidence invalidates it.
Only Low components keyed by a Wave39-touched ``genus x trait_name`` are
removed and rebuilt; direct High/Medium cells are never downgraded.  Losses
are reported explicitly so later acquisition can repair them with direct
evidence instead of hiding counterevidence.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.floraweb_synonym_recovery import SOURCE_GROUP
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
from island_v2.wave38_trait_overlay import (
    _prior_rule_table,
    _touched_keys,
)
from island_v2.wave38_trait_overlay import (
    build_touched_rule_rebuild as _build_touched_rule_rebuild,
)


def _rule_lineage(genus: str, trait: str, state_set: str, prefix: str) -> str:
    return f"{prefix}:{genus}:{trait}:{state_set}"


def touched_rebuild_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--formal-resolved-direct-csv", required=True, type=Path)
    parser.add_argument("--reviewed-direct-csv", required=True, type=Path)
    parser.add_argument(
        "--supplemental-direct-evidence-csv", action="append", default=[], type=Path
    )
    parser.add_argument("--master-csv", required=True, type=Path)
    parser.add_argument("--ontology-yaml", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    summary = _build_touched_rule_rebuild(
        formal_resolved_direct_csv=args.formal_resolved_direct_csv,
        reviewed_direct_csv=args.reviewed_direct_csv,
        supplemental_direct_evidence_csvs=tuple(args.supplemental_direct_evidence_csv),
        master_csv=args.master_csv,
        ontology_yaml=args.ontology_yaml,
        output_dir=args.output_dir,
        source_group=SOURCE_GROUP,
        wave_label="wave39",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


def _rebase_touched_low_components(
    baseline: pd.DataFrame,
    touched_keys: set[tuple[str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    result = baseline.copy()
    audit: list[dict[str, str]] = []
    low_indices = result.index[result["quality"].eq("low")]
    for index in low_indices:
        species = str(result.at[index, "accepted_species"])
        axis = str(result.at[index, "axis"])
        genus = species.split()[0]
        composition = _parse_composition(result.at[index, "trait_composition"], axis=axis)
        removed = sorted(trait for trait in composition if (genus, trait) in touched_keys)
        if not removed:
            continue
        for trait in removed:
            audit.append(
                {
                    "accepted_species": species,
                    "axis": axis,
                    "trait_name": trait,
                    "action": "remove_touched_prior_low_component",
                    "prior_state_set": composition[trait],
                }
            )
            del composition[trait]
        if composition:
            result.at[index, "trait_composition"] = _serialize_composition(composition)
            result.at[index, "trait_names"] = "|".join(sorted(composition))
            result.at[index, "source_groups"] = "retained_prior_trait_specific_validated_low"
            result.at[index, "source_lineages"] = "|".join(
                _rule_lineage(genus, trait, states, "retained-prior-low")
                for trait, states in sorted(composition.items())
            )
        else:
            for column in (
                "trait_composition",
                "trait_names",
                "source_groups",
                "source_lineages",
                "quality",
            ):
                result.at[index, column] = ""
    return result, pd.DataFrame(audit)


def build_wave39_lossless_overlay(
    *,
    baseline_csv: Path,
    all_evidence_dir: Path,
    reviewed_direct_csv: Path,
    checkpoint_summary_json: Path,
    previous_rule_audit_csvs: tuple[Path, ...],
    output_dir: Path,
    expected_species: int = EXPECTED_SPECIES,
) -> dict[str, Any]:
    """Add only new direct evidence and new eligible trait-specific rules.

    The formal secondary ledger is append-only.  Newly observed counterexamples
    are still reported as a shadow invalidation audit, but they do not silently
    rewrite previously published Low cells.  A separate scientifically rebased
    diagnostic is available through :func:`build_wave39_overlay`.
    """

    resolved_path = all_evidence_dir / "resolved_direct_species_trait.csv.gz"
    rebuilt_low_path = all_evidence_dir / "rebuilt_all_evidence_validated_low.csv.gz"
    rules_path = all_evidence_dir / "trait_specific_genus_rule_audit.csv.gz"
    required = (
        baseline_csv,
        reviewed_direct_csv,
        checkpoint_summary_json,
        resolved_path,
        rebuilt_low_path,
        rules_path,
        *previous_rule_audit_csvs,
    )
    if missing := [str(path) for path in required if not path.is_file()]:
        raise ValueError(f"Wave39 overlay inputs missing: {missing}")

    baseline = pd.read_csv(baseline_csv, dtype=str).fillna("")
    _validate_coverage(baseline, expected_species)
    reviewed = pd.read_csv(reviewed_direct_csv, dtype=str).fillna("")
    resolved = pd.read_csv(resolved_path, dtype=str).fillna("")
    rebuilt_low = pd.read_csv(rebuilt_low_path, dtype=str).fillna("")
    rules = pd.read_csv(rules_path, dtype=str).fillna("")
    prior_rules = _prior_rule_table(previous_rule_audit_csvs)
    touched_keys = _touched_keys(reviewed)

    direct = _validate_direct(
        resolved.loc[
            resolved["source_groups"].map(lambda value: SOURCE_GROUP in _split_pipe(value))
        ].copy()
    )
    current_rules = rules.loc[rules["setting"].eq("current_min3") & _eligible(rules)].copy()
    current_rules = current_rules.loc[
        current_rules.apply(
            lambda row: (str(row["genus"]), str(row["trait_name"])) in touched_keys,
            axis=1,
        )
    ]
    prior_rules = prior_rules.loc[
        prior_rules.apply(
            lambda row: (str(row["genus"]), str(row["trait_name"])) in touched_keys,
            axis=1,
        )
    ]
    key_columns = ["genus", "axis", "trait_name"]
    current_by_key = current_rules.set_index(key_columns)["inferred_state_set"].to_dict()
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
    low = low.merge(
        baseline[["accepted_species", "axis", "quality"]],
        on=["accepted_species", "axis"],
        how="left",
        validate="many_to_one",
        suffixes=("", "_baseline"),
    )
    low = low.loc[low["quality_baseline"].eq("")].drop(columns="quality_baseline")
    _validate_low(low, new_rules)

    result = baseline.copy().set_index(["accepted_species", "axis"])
    changes: list[dict[str, str]] = []
    for row in _aggregate_direct(direct).itertuples(index=False):
        key = (row.accepted_species, row.axis)
        if key not in result.index:
            raise ValueError(f"direct evidence is outside the fixed universe: {key}")
        before_quality = str(result.loc[key, "quality"])
        if before_quality in {"", "low"}:
            composition = _parse_composition(row.trait_composition, axis=row.axis)
            groups = _split_pipe(row.source_groups)
            lineages = _split_pipe(row.source_lineages)
        else:
            composition = _parse_composition(result.loc[key, "trait_composition"], axis=row.axis)
            composition.update(_parse_composition(row.trait_composition, axis=row.axis))
            groups = _split_pipe(result.loc[key, "source_groups"]) | _split_pipe(
                row.source_groups
            )
            lineages = _split_pipe(result.loc[key, "source_lineages"]) | _split_pipe(
                row.source_lineages
            )
        after_quality = max((before_quality, row.quality), key=QUALITY_RANK.__getitem__)
        result.loc[key, "trait_composition"] = _serialize_composition(composition)
        result.loc[key, "trait_names"] = "|".join(sorted(composition))
        result.loc[key, "source_groups"] = "|".join(sorted(groups))
        result.loc[key, "source_lineages"] = "|".join(sorted(lineages))
        result.loc[key, "quality"] = after_quality
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
            }
        )
    for row in _aggregate_low(low).itertuples(index=False):
        key = (row.accepted_species, row.axis)
        if str(result.loc[key, "quality"]):
            raise ValueError(f"Low attempts to replace a resolved cell: {key}")
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

    result = result.reset_index().sort_values(["accepted_species", "axis"]).reset_index(drop=True)
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
    if loss:
        raise ValueError(f"Wave39 formal overlay lost {loss} previously filled cells")
    if (
        comparison["quality_after"].map(QUALITY_RANK)
        < comparison["quality_before"].map(QUALITY_RANK)
    ).any():
        raise ValueError("Wave39 formal overlay downgraded an existing quality")

    output_dir.mkdir(parents=True, exist_ok=True)
    changes_frame = pd.DataFrame(changes)
    outputs = {
        "wave39_species_axis_coverage.csv.gz": result,
        "wave39_resolved_direct_species_trait.csv.gz": direct,
        "wave39_new_validated_low_species_trait.csv.gz": low,
        "wave39_new_trait_specific_genus_rule_audit.csv.gz": new_rules,
        "wave39_change_audit.csv.gz": changes_frame,
    }
    for name, frame in outputs.items():
        _write_gzip_csv(frame, output_dir / name)
    summary: dict[str, Any] = {
        "contract": "wave39_lossless_additive_secondary_overlay_v1",
        "formal_wave33_run_id": 32932103226,
        "formal_wave37_run_id": 33143109604,
        "formal_wave38_run_id": 33149235490,
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
            "reviewed_direct_species_trait": len(reviewed),
            "resolved_direct_species_trait": len(direct),
            "new_validated_low_species_trait": len(low),
        },
        "new_eligible_rules": [
            f"{row.genus} x {row.trait_name}" for row in new_rules.itertuples(index=False)
        ],
        "shadow_invalidated_prior_rules": [
            " x ".join((key[0], key[2])) for key in invalidated_keys
        ],
        "shadow_changed_prior_rule_states": [
            " x ".join((key[0], key[2])) for key in changed_keys
        ],
        "checks": {
            "fixed_denominator": True,
            "quality_precedence_high_medium_low": True,
            "direct_conflicts_excluded": True,
            "trait_specific_genus_join": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
            "baseline_loss_zero": loss == 0,
            "shadow_invalidations_not_silently_deleted": True,
        },
        "checkpoint": json.loads(checkpoint_summary_json.read_text(encoding="utf-8")),
        "input_sha256": {path.name: _sha256(path) for path in required},
    }
    summary["artifact_sha256"] = {name: _sha256(output_dir / name) for name in outputs}
    (output_dir / "wave39_coverage_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def build_wave39_overlay(
    *,
    baseline_csv: Path,
    all_evidence_dir: Path,
    reviewed_direct_csv: Path,
    checkpoint_summary_json: Path,
    previous_rule_audit_csvs: tuple[Path, ...],
    output_dir: Path,
    expected_species: int = EXPECTED_SPECIES,
) -> dict[str, Any]:
    resolved_path = all_evidence_dir / "resolved_direct_species_trait.csv.gz"
    rebuilt_low_path = all_evidence_dir / "rebuilt_all_evidence_validated_low.csv.gz"
    rules_path = all_evidence_dir / "trait_specific_genus_rule_audit.csv.gz"
    required = (
        baseline_csv,
        reviewed_direct_csv,
        checkpoint_summary_json,
        resolved_path,
        rebuilt_low_path,
        rules_path,
        *previous_rule_audit_csvs,
    )
    if missing := [str(path) for path in required if not path.is_file()]:
        raise ValueError(f"Wave39 overlay inputs missing: {missing}")

    baseline = pd.read_csv(baseline_csv, dtype=str).fillna("")
    _validate_coverage(baseline, expected_species)
    reviewed = pd.read_csv(reviewed_direct_csv, dtype=str).fillna("")
    resolved = pd.read_csv(resolved_path, dtype=str).fillna("")
    rebuilt_low = pd.read_csv(rebuilt_low_path, dtype=str).fillna("")
    rules = pd.read_csv(rules_path, dtype=str).fillna("")
    prior_rules = _prior_rule_table(previous_rule_audit_csvs)
    touched_keys = _touched_keys(reviewed)

    direct = _validate_direct(
        resolved.loc[
            resolved["source_groups"].map(lambda value: SOURCE_GROUP in _split_pipe(value))
        ].copy()
    )
    current_rules = rules.loc[rules["setting"].eq("current_min3") & _eligible(rules)].copy()
    current_rules = current_rules.loc[
        current_rules.apply(
            lambda row: (str(row["genus"]), str(row["trait_name"])) in touched_keys,
            axis=1,
        )
    ]
    prior_rules = prior_rules.loc[
        prior_rules.apply(
            lambda row: (str(row["genus"]), str(row["trait_name"])) in touched_keys,
            axis=1,
        )
    ]
    key_columns = ["genus", "axis", "trait_name"]
    current_by_key = current_rules.set_index(key_columns)["inferred_state_set"].to_dict()
    prior_by_key = prior_rules.set_index(key_columns)["inferred_state_set"].to_dict()
    invalidated_keys = sorted(set(prior_by_key) - set(current_by_key))
    changed_keys = sorted(
        key
        for key in set(prior_by_key) & set(current_by_key)
        if prior_by_key[key] != current_by_key[key]
    )
    new_keys = sorted(set(current_by_key) - set(prior_by_key))

    current_rule_keys = set(current_by_key)
    low = rebuilt_low.loc[
        rebuilt_low.apply(
            lambda row: (
                (
                    str(row["genus"]),
                    str(row["axis"]),
                    str(row["trait_name"]),
                )
                in current_rule_keys
            ),
            axis=1,
        )
    ].copy()
    # ``rules`` also contains the explicitly diagnostic min-species=2
    # settings.  They must remain in the audit artifact, but must never be
    # admitted to (or validated as part of) the strict current-min3 overlay.
    _validate_low(low, current_rules)

    before = _coverage_summary(baseline, expected_species)
    rebased, low_removal_audit = _rebase_touched_low_components(baseline, touched_keys)
    result = rebased.set_index(["accepted_species", "axis"])
    changes: list[dict[str, str]] = []

    for row in _aggregate_direct(direct).itertuples(index=False):
        key = (row.accepted_species, row.axis)
        if key not in result.index:
            raise ValueError(f"direct evidence is outside the fixed universe: {key}")
        before_quality = str(result.loc[key, "quality"])
        if before_quality in {"", "low"}:
            composition = _parse_composition(row.trait_composition, axis=row.axis)
            groups = _split_pipe(row.source_groups)
            lineages = _split_pipe(row.source_lineages)
        else:
            composition = _parse_composition(result.loc[key, "trait_composition"], axis=row.axis)
            composition.update(_parse_composition(row.trait_composition, axis=row.axis))
            groups = _split_pipe(result.loc[key, "source_groups"]) | _split_pipe(row.source_groups)
            lineages = _split_pipe(result.loc[key, "source_lineages"]) | _split_pipe(
                row.source_lineages
            )
        after_quality = max((before_quality, row.quality), key=QUALITY_RANK.__getitem__)
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
                "quality_before_rebased": before_quality,
                "quality_after": after_quality,
                "trait_names": row.trait_names,
            }
        )

    low_axis = _aggregate_low(low)
    for row in low_axis.itertuples(index=False):
        key = (row.accepted_species, row.axis)
        if key not in result.index:
            raise ValueError(f"Low evidence is outside the fixed universe: {key}")
        before_quality = str(result.loc[key, "quality"])
        if before_quality in {"high", "medium"}:
            continue
        new_composition = _parse_composition(row.trait_composition, axis=row.axis)
        if before_quality == "low":
            composition = _parse_composition(result.loc[key, "trait_composition"], axis=row.axis)
            composition.update(new_composition)
            groups = _split_pipe(result.loc[key, "source_groups"]) | {
                "wave39_trait_specific_validated_low"
            }
            lineages = _split_pipe(result.loc[key, "source_lineages"]) | _split_pipe(
                row.source_lineages
            )
        else:
            composition = new_composition
            groups = {"wave39_trait_specific_validated_low"}
            lineages = _split_pipe(row.source_lineages)
        result.loc[key, "trait_composition"] = _serialize_composition(composition)
        result.loc[key, "trait_names"] = "|".join(sorted(composition))
        result.loc[key, "source_groups"] = "|".join(sorted(groups))
        result.loc[key, "source_lineages"] = "|".join(sorted(lineages))
        result.loc[key, "quality"] = "low"
        changes.append(
            {
                "accepted_species": row.accepted_species,
                "axis": row.axis,
                "action": "validated_low_fill" if before_quality == "" else "validated_low_refresh",
                "quality_before_rebased": before_quality,
                "quality_after": "low",
                "trait_names": row.trait_names,
            }
        )

    result = result.reset_index().sort_values(["accepted_species", "axis"]).reset_index(drop=True)
    _validate_coverage(result, expected_species)
    comparison = baseline[["accepted_species", "axis", "quality"]].merge(
        result[["accepted_species", "axis", "quality"]],
        on=["accepted_species", "axis"],
        suffixes=("_before", "_after"),
        validate="one_to_one",
    )
    before_filled = comparison["quality_before"].ne("")
    after_filled = comparison["quality_after"].ne("")
    gross_gain = int((~before_filled & after_filled).sum())
    gross_loss = int((before_filled & ~after_filled).sum())
    if (
        comparison.loc[comparison["quality_before"].isin(["high", "medium"]), "quality_after"].map(
            QUALITY_RANK
        )
        < comparison.loc[
            comparison["quality_before"].isin(["high", "medium"]), "quality_before"
        ].map(QUALITY_RANK)
    ).any():
        raise ValueError("Wave39 downgraded a prior direct High/Medium cell")
    after = _coverage_summary(result, expected_species)
    changes_frame = pd.DataFrame(changes)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "wave39_species_axis_coverage.csv.gz": result,
        "wave39_resolved_direct_species_trait.csv.gz": direct,
        "wave39_candidate_validated_low_species_trait.csv.gz": low,
        "wave39_trait_specific_genus_rule_audit.csv.gz": rules,
        "wave39_low_rebase_audit.csv.gz": low_removal_audit,
        "wave39_coverage_changes.csv.gz": changes_frame,
    }
    for name, frame in outputs.items():
        _write_gzip_csv(frame, output_dir / name)
    summary: dict[str, Any] = {
        "contract": "wave39_scientifically_rebased_trait_overlay_v1",
        "before": before,
        "after": after,
        "gain": {
            "gross_filled_species_axis": gross_gain,
            "gross_unresolved_species_axis": gross_loss,
            "net_filled_species_axis": gross_gain - gross_loss,
            "by_axis_net": {
                axis: after["by_axis"][axis]["filled_species"]
                - before["by_axis"][axis]["filled_species"]
                for axis in after["by_axis"]
            },
        },
        "direct": {
            "reviewed_species_trait": len(reviewed),
            "resolved_species_trait": len(direct),
            "excluded_conflicted_species_trait": len(reviewed) - len(direct),
        },
        "rules": {
            "new": [" x ".join(key[::2]) for key in new_keys],
            "invalidated": [" x ".join(key[::2]) for key in invalidated_keys],
            "changed": [" x ".join(key[::2]) for key in changed_keys],
            "candidate_low_species_trait": len(low),
        },
        "checks": {
            "fixed_denominator": True,
            "direct_high_medium_not_downgraded": True,
            "stale_touched_low_retained": False,
            "trait_specific_genus_join": True,
            "family_inference": False,
            "global_fallback": False,
            "loss_zero": gross_loss == 0,
        },
        "checkpoint": json.loads(checkpoint_summary_json.read_text(encoding="utf-8")),
        "input_sha256": {path.name: _sha256(path) for path in required},
    }
    summary["artifact_sha256"] = {name: _sha256(output_dir / name) for name in outputs}
    (output_dir / "wave39_coverage_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=("formal-lossless", "rebased-diagnostic"),
        default="formal-lossless",
    )
    parser.add_argument("--baseline-csv", required=True, type=Path)
    parser.add_argument("--all-evidence-dir", required=True, type=Path)
    parser.add_argument("--reviewed-direct-csv", required=True, type=Path)
    parser.add_argument("--checkpoint-summary-json", required=True, type=Path)
    parser.add_argument("--previous-rule-audit-csv", action="append", default=[], type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-species", type=int, default=EXPECTED_SPECIES)
    args = parser.parse_args()
    builder = (
        build_wave39_lossless_overlay
        if args.mode == "formal-lossless"
        else build_wave39_overlay
    )
    summary = builder(
        baseline_csv=args.baseline_csv,
        all_evidence_dir=args.all_evidence_dir,
        reviewed_direct_csv=args.reviewed_direct_csv,
        checkpoint_summary_json=args.checkpoint_summary_json,
        previous_rule_audit_csvs=tuple(args.previous_rule_audit_csv),
        output_dir=args.output_dir,
        expected_species=args.expected_species,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
