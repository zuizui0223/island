"""Recompute Wave48-touched genus/axis/trait rules from formal Wave47 evidence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.all_evidence_trait_audit import (
    apply_genus_rules,
    build_rule_audit,
    dedupe_direct_lineages,
    direct_evidence_from_integrated,
    load_external_congener_support,
    load_ontology,
    resolve_direct_cells,
)
from island_v2.wave37_europe_pmc_checkpoint import EXPECTED_SPECIES, _sha256, _write_gzip_csv

EXPECTED_NEW_RULE = (
    "Pelargonium",
    "reproductive_assurance",
    "self_incompatibility",
)
EXPECTED_BLOCKED_RULES = frozenset(
    {
        ("Spermacoce", "reproductive_assurance", "autonomous_selfing_capacity"),
        ("Torenia", "reproductive_assurance", "autonomous_selfing_capacity"),
        ("Torenia", "reproductive_assurance", "self_incompatibility"),
    }
)
EXPECTED_COUNTEREXAMPLE_RULE = (
    "Callicarpa",
    "reproductive_assurance",
    "self_incompatibility",
)


def _eligible(frame: pd.DataFrame) -> pd.Series:
    return frame["eligible"].astype(str).str.casefold().isin({"true", "1"})


def _reconstruct_lineages(cells: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row in cells.itertuples(index=False):
        lineages = [token for token in str(row.source_lineages).split("|") if token]
        groups = [token for token in str(row.source_groups).split("|") if token]
        source_group = "|".join(sorted(groups))
        for lineage in lineages:
            rows.append(
                {
                    "accepted_species": str(row.accepted_species),
                    "genus": str(row.genus),
                    "axis": str(row.axis),
                    "trait_name": str(row.trait_name),
                    "state_set": str(row.state_set),
                    "source_lineage": lineage,
                    "source_group": source_group,
                    "ontology_valid": True,
                    "lineage_internal_conflict": False,
                }
            )
    return pd.DataFrame(rows)


def _touch_mask(frame: pd.DataFrame, touched: set[tuple[str, str, str]]) -> pd.Series:
    return frame.apply(
        lambda row: (
            str(row["genus"]),
            str(row["axis"]),
            str(row["trait_name"]),
        )
        in touched,
        axis=1,
    )


def _key_set(frame: pd.DataFrame) -> set[tuple[str, str, str]]:
    return set(
        frame[["genus", "axis", "trait_name"]].itertuples(index=False, name=None)
    )


def build_incremental_audit(
    *,
    master_csv: Path,
    ontology_yaml: Path,
    baseline_coverage_csv: Path,
    previous_rule_audit_csv: Path,
    previous_resolved_direct_csv: Path,
    previous_external_resolved_csv: Path,
    previous_external_conflicts_csv: Path,
    previous_rebuilt_low_csv: Path,
    new_direct_evidence_csv: Path,
    new_external_evidence_csv: Path,
    output_dir: Path,
    expected_species: int = EXPECTED_SPECIES,
    expected_direct_rows: int = 12,
    expected_external_rows: int = 1,
    expected_new_rule: tuple[str, str, str] = EXPECTED_NEW_RULE,
    expected_blocked_rules: frozenset[tuple[str, str, str]] = EXPECTED_BLOCKED_RULES,
    expected_counterexample_rule: tuple[str, str, str] = EXPECTED_COUNTEREXAMPLE_RULE,
) -> dict[str, Any]:
    required = (
        master_csv,
        ontology_yaml,
        baseline_coverage_csv,
        previous_rule_audit_csv,
        previous_resolved_direct_csv,
        previous_external_resolved_csv,
        previous_external_conflicts_csv,
        previous_rebuilt_low_csv,
        new_direct_evidence_csv,
        new_external_evidence_csv,
    )
    if missing := [str(path) for path in required if not path.is_file()]:
        raise ValueError(f"Wave48 incremental audit inputs missing: {missing}")

    baseline = pd.read_csv(
        baseline_coverage_csv,
        usecols=["accepted_species", "axis"],
        dtype=str,
    ).fillna("")
    target_species = set(baseline["accepted_species"])
    if len(baseline) != expected_species * 3 or len(target_species) != expected_species:
        raise ValueError("Wave48 baseline species-axis denominator mismatch")
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    master = master.loc[master["accepted_species"].isin(target_species)].copy()
    master = master.drop_duplicates("accepted_species")
    if len(master) != expected_species:
        raise ValueError("Wave48 master does not reproduce the fixed baseline universe")

    previous_rules = pd.read_csv(previous_rule_audit_csv, dtype=str).fillna("")
    previous_direct = pd.read_csv(previous_resolved_direct_csv, dtype=str).fillna("")
    previous_external = pd.read_csv(previous_external_resolved_csv, dtype=str).fillna("")
    previous_conflicts = pd.read_csv(previous_external_conflicts_csv, dtype=str).fillna("")
    old_low = pd.read_csv(previous_rebuilt_low_csv, dtype=str).fillna("")
    ontology = load_ontology(ontology_yaml)

    new_direct_raw = direct_evidence_from_integrated(new_direct_evidence_csv)
    new_external_raw = load_external_congener_support((new_external_evidence_csv,))
    if (
        len(new_direct_raw) != expected_direct_rows
        or len(new_external_raw) != expected_external_rows
    ):
        raise ValueError("Wave48 frozen direct/external packet changed")
    if not new_direct_raw["accepted_species"].isin(target_species).all():
        raise ValueError("Wave48 direct evidence is outside the fixed target")
    if new_external_raw["accepted_species"].isin(target_species).any():
        raise ValueError("Wave48 external evidence entered the fixed target")

    new_direct_lineages, direct_duplicates = dedupe_direct_lineages(
        new_direct_raw, ontology
    )
    new_direct_cells, new_direct_audit = resolve_direct_cells(new_direct_lineages)
    new_external_lineages, external_duplicates = dedupe_direct_lineages(
        new_external_raw, ontology
    )
    new_external_cells, new_external_audit = resolve_direct_cells(new_external_lineages)
    if len(direct_duplicates) or len(external_duplicates):
        raise ValueError("Wave48 contains duplicate source-lineage evidence")
    if (
        len(new_direct_cells) != expected_direct_rows
        or len(new_external_cells) != expected_external_rows
    ):
        raise ValueError("Wave48 evidence contains an unresolved direct conflict")
    if not new_direct_audit["resolution_status"].eq("resolved").all():
        raise ValueError("Wave48 direct evidence failed resolution")
    if not new_external_audit["resolution_status"].eq("resolved").all():
        raise ValueError("Wave48 external evidence failed resolution")

    for frame in (
        new_direct_cells,
        new_external_cells,
        new_direct_lineages,
        new_external_lineages,
    ):
        frame["genus"] = frame["accepted_species"].str.split().str[0]

    previous_direct_keys = set(
        previous_direct[["accepted_species", "trait_name"]].itertuples(
            index=False, name=None
        )
    )
    new_direct_keys = set(
        new_direct_cells[["accepted_species", "trait_name"]].itertuples(
            index=False, name=None
        )
    )
    if previous_direct_keys & new_direct_keys:
        raise ValueError("Wave48 tries to re-ingest a completed direct species-trait")
    previous_external_keys = set(
        previous_external[["accepted_species", "trait_name"]].itertuples(
            index=False, name=None
        )
    )
    new_external_keys = set(
        new_external_cells[["accepted_species", "trait_name"]].itertuples(
            index=False, name=None
        )
    )
    if previous_external_keys & new_external_keys:
        raise ValueError("Wave48 tries to re-ingest a completed external species-trait")

    combined_direct = pd.concat(
        [previous_direct, new_direct_cells], ignore_index=True, sort=False
    ).fillna("")
    combined_external = pd.concat(
        [previous_external, new_external_cells], ignore_index=True, sort=False
    ).fillna("")
    touched = _key_set(new_direct_lineages) | _key_set(new_external_lineages)
    prior_direct_touched = previous_direct.loc[
        _touch_mask(previous_direct, touched)
    ].copy()
    prior_external_touched = previous_external.loc[
        _touch_mask(previous_external, touched)
    ].copy()
    all_cells = pd.concat(
        [combined_direct, combined_external], ignore_index=True, sort=False
    ).fillna("")
    touched_cells = all_cells.loc[_touch_mask(all_cells, touched)].copy()
    reconstructed = _reconstruct_lineages(
        pd.concat(
            [prior_direct_touched, prior_external_touched],
            ignore_index=True,
            sort=False,
        ).fillna("")
    )
    touched_lineages = pd.concat(
        [reconstructed, new_direct_lineages, new_external_lineages],
        ignore_index=True,
        sort=False,
    ).fillna("")
    touched_old_low = old_low.loc[_touch_mask(old_low, touched)].copy()
    recalculated = build_rule_audit(touched_cells, touched_lineages, touched_old_low)
    if recalculated.empty:
        raise ValueError("Wave48 produced no touched rule audit")

    keep = ~_touch_mask(previous_rules, touched)
    rules = pd.concat(
        [previous_rules.loc[keep], recalculated], ignore_index=True, sort=False
    ).fillna("")
    rules = rules.sort_values(
        ["setting", "genus", "axis", "trait_name"], kind="stable"
    ).reset_index(drop=True)

    key_columns = ["genus", "axis", "trait_name"]
    previous_eligible = previous_rules.loc[
        previous_rules["setting"].eq("current_min3") & _eligible(previous_rules)
    ]
    current_eligible = rules.loc[
        rules["setting"].eq("current_min3") & _eligible(rules)
    ]
    previous_key_set = set(
        previous_eligible[key_columns].itertuples(index=False, name=None)
    )
    current_key_set = set(current_eligible[key_columns].itertuples(index=False, name=None))
    newly_eligible_keys = current_key_set - previous_key_set
    invalidated_keys = previous_key_set - current_key_set
    if newly_eligible_keys != {expected_new_rule}:
        raise ValueError(f"Wave48 unexpected new rules: {sorted(newly_eligible_keys)}")
    new_rule = current_eligible.loc[
        current_eligible.apply(
            lambda row: tuple(str(row[column]) for column in key_columns)
            == expected_new_rule,
            axis=1,
        )
    ].copy()

    touched_current = recalculated.loc[recalculated["setting"].eq("current_min3")].copy()
    touched_current_by_key = {
        tuple(str(row[column]) for column in key_columns): row
        for _, row in touched_current.iterrows()
    }
    blocked_failures = []
    for key in expected_blocked_rules:
        row = touched_current_by_key.get(key)
        if (
            row is None
            or str(row["eligible"]).casefold() in {"true", "1"}
            or float(row["lineage_loo_accuracy"]) != 0.0
        ):
            blocked_failures.append(key)
    if blocked_failures:
        raise ValueError(
            f"Wave48 single-lineage rules were not fail-closed: {blocked_failures}"
        )
    callicarpa = touched_current_by_key.get(expected_counterexample_rule)
    if (
        callicarpa is None
        or str(callicarpa["eligible"]).casefold() in {"true", "1"}
        or int(callicarpa["counterexample_species"]) < 1
    ):
        raise ValueError("Wave48 Callicarpa counterexample did not block the rule")

    rebuilt_low = apply_genus_rules(master, all_cells, new_rule, "current_min3")
    if rebuilt_low.empty or set(rebuilt_low["axis"]) != {"reproductive_assurance"}:
        raise ValueError("Wave48 produced no trait-specific reproductive Low")

    combined_conflicts = pd.concat(
        [previous_conflicts, new_external_audit], ignore_index=True, sort=False
    ).fillna("")
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "resolved_direct_species_trait.csv.gz": combined_direct,
        "rebuilt_all_evidence_validated_low.csv.gz": rebuilt_low,
        "trait_specific_genus_rule_audit.csv.gz": rules,
        "external_congener_resolved_species_trait.csv.gz": combined_external,
        "external_congener_source_lineage_conflicts.csv.gz": combined_conflicts,
        "wave48_direct_source_lineage_resolution_audit.csv.gz": new_direct_audit,
    }
    for name, frame in outputs.items():
        _write_gzip_csv(frame, output_dir / name)

    rule_row = new_rule.iloc[0]
    summary: dict[str, Any] = {
        "contract": "wave48_incremental_all_evidence_touched_rule_rebuild_v1",
        "fixed_target_species": expected_species,
        "touched_genus_axis_trait": [
            " x ".join(key) for key in sorted(touched)
        ],
        "new_direct_rows": len(new_direct_raw),
        "new_direct_species_trait": len(new_direct_cells),
        "new_external_rows": len(new_external_raw),
        "new_external_species_trait": len(new_external_cells),
        "new_eligible_rules": [f"{expected_new_rule[0]} x {expected_new_rule[2]}"],
        "invalidated_prior_rules": [" x ".join(key) for key in sorted(invalidated_keys)],
        "counterexample_blocked_rules": [
            f"{expected_counterexample_rule[0]} x {expected_counterexample_rule[2]}"
        ],
        "single_lineage_blocked_rules": [
            " x ".join(key) for key in sorted(expected_blocked_rules)
        ],
        "new_rule": {
            "n_direct_species": int(rule_row["n_direct_species"]),
            "dominance": float(rule_row["dominance"]),
            "species_loo_accuracy": float(rule_row["species_loo_accuracy"]),
            "lineage_loo_accuracy": float(rule_row["lineage_loo_accuracy"]),
            "inferred_value": str(rule_row["inferred_value"]),
        },
        "new_validated_low_species_trait": len(rebuilt_low),
        "source_lineage_audit": {
            "new_direct": {
                "rows": len(new_direct_raw),
                "resolved_species_trait": len(new_direct_cells),
            },
            "external_congener_support": {
                "files": 1,
                "rows": len(new_external_raw),
                "resolved_species_trait": len(new_external_cells),
                "entered_confirmatory_direct_coverage": 0,
            },
        },
        "checks": {
            "fixed_denominator": True,
            "incremental_scope_only_touched_genus_axis_trait": True,
            "new_source_lineages_not_previously_ingested": True,
            "external_congener_not_confirmatory_direct": True,
            "trait_specific_genus_join": True,
            "minimum_species_three": int(rule_row["n_direct_species"]) >= 3,
            "leave_one_source_lineage_out_passed": float(rule_row["lineage_loo_accuracy"])
            == 1.0,
            "single_source_lineage_rules_rejected": True,
            "callicarpa_counterexample_retained": True,
            "reproductive_traits_not_interchanged": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
        },
        "input_sha256": {str(path): _sha256(path) for path in required},
        "artifact_sha256": {name: _sha256(output_dir / name) for name in outputs},
    }
    summary_path = output_dir / "all_evidence_trait_coverage_summary.json"
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master-csv", required=True, type=Path)
    parser.add_argument("--ontology-yaml", required=True, type=Path)
    parser.add_argument("--baseline-coverage-csv", required=True, type=Path)
    parser.add_argument("--previous-rule-audit-csv", required=True, type=Path)
    parser.add_argument("--previous-resolved-direct-csv", required=True, type=Path)
    parser.add_argument("--previous-external-resolved-csv", required=True, type=Path)
    parser.add_argument("--previous-external-conflicts-csv", required=True, type=Path)
    parser.add_argument("--previous-rebuilt-low-csv", required=True, type=Path)
    parser.add_argument("--new-direct-evidence-csv", required=True, type=Path)
    parser.add_argument("--new-external-evidence-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-species", type=int, default=EXPECTED_SPECIES)
    args = parser.parse_args()
    summary = build_incremental_audit(
        master_csv=args.master_csv,
        ontology_yaml=args.ontology_yaml,
        baseline_coverage_csv=args.baseline_coverage_csv,
        previous_rule_audit_csv=args.previous_rule_audit_csv,
        previous_resolved_direct_csv=args.previous_resolved_direct_csv,
        previous_external_resolved_csv=args.previous_external_resolved_csv,
        previous_external_conflicts_csv=args.previous_external_conflicts_csv,
        previous_rebuilt_low_csv=args.previous_rebuilt_low_csv,
        new_direct_evidence_csv=args.new_direct_evidence_csv,
        new_external_evidence_csv=args.new_external_evidence_csv,
        output_dir=args.output_dir,
        expected_species=args.expected_species,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
