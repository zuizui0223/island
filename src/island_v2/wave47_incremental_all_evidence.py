"""Recompute only Wave47-touched genus/trait rules from formal Wave46 evidence."""

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
    load_external_congener_support,
    load_ontology,
    resolve_direct_cells,
)
from island_v2.wave37_europe_pmc_checkpoint import EXPECTED_SPECIES, _sha256, _write_gzip_csv


def _eligible(frame: pd.DataFrame) -> pd.Series:
    return frame["eligible"].astype(str).str.casefold().isin({"true", "1"})


def _reconstruct_lineages(cells: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row in cells.itertuples(index=False):
        lineages = [token for token in str(row.source_lineages).split("|") if token]
        for lineage in lineages:
            rows.append(
                {
                    "accepted_species": str(row.accepted_species),
                    "genus": str(row.genus),
                    "axis": str(row.axis),
                    "trait_name": str(row.trait_name),
                    "state_set": str(row.state_set),
                    "source_lineage": lineage,
                    "source_group": str(row.source_groups),
                    "ontology_valid": True,
                    "lineage_internal_conflict": False,
                }
            )
    return pd.DataFrame(rows)


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
    new_external_evidence_csv: Path,
    output_dir: Path,
    expected_species: int = EXPECTED_SPECIES,
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
        new_external_evidence_csv,
    )
    if missing := [str(path) for path in required if not path.is_file()]:
        raise ValueError(f"Wave47 incremental audit inputs missing: {missing}")

    baseline = pd.read_csv(
        baseline_coverage_csv,
        usecols=["accepted_species", "axis"],
        dtype=str,
    ).fillna("")
    target_species = set(baseline["accepted_species"])
    if len(baseline) != expected_species * 3 or len(target_species) != expected_species:
        raise ValueError("Wave47 baseline species-axis denominator mismatch")
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    master = master.loc[master["accepted_species"].isin(target_species)].copy()
    master = master.drop_duplicates("accepted_species")
    if len(master) != expected_species:
        raise ValueError("Wave47 master does not reproduce the fixed baseline universe")

    previous_rules = pd.read_csv(previous_rule_audit_csv, dtype=str).fillna("")
    direct = pd.read_csv(previous_resolved_direct_csv, dtype=str).fillna("")
    previous_external = pd.read_csv(previous_external_resolved_csv, dtype=str).fillna("")
    previous_conflicts = pd.read_csv(previous_external_conflicts_csv, dtype=str).fillna("")
    old_low = pd.read_csv(previous_rebuilt_low_csv, dtype=str).fillna("")
    ontology = load_ontology(ontology_yaml)
    new_raw = load_external_congener_support((new_external_evidence_csv,))
    if len(new_raw) != 2 or set(new_raw["accepted_species"]) != {"Buxus wallichiana"}:
        raise ValueError("Wave47 new external packet changed")
    if new_raw["accepted_species"].isin(target_species).any():
        raise ValueError("Wave47 external evidence entered the fixed target")

    new_lineages, new_duplicates = dedupe_direct_lineages(new_raw, ontology)
    new_cells, new_cell_audit = resolve_direct_cells(new_lineages)
    if len(new_duplicates) or len(new_cells) != 2:
        raise ValueError("Wave47 external source is duplicated or unresolved")
    if not new_cell_audit["resolution_status"].eq("resolved").all():
        raise ValueError("Wave47 external source has a direct conflict")
    new_cells["genus"] = new_cells["accepted_species"].str.split().str[0]
    new_lineages["genus"] = new_lineages["accepted_species"].str.split().str[0]

    prior_keys = set(
        previous_external[["accepted_species", "trait_name"]].itertuples(index=False, name=None)
    )
    new_keys = set(new_cells[["accepted_species", "trait_name"]].itertuples(index=False, name=None))
    if prior_keys & new_keys:
        raise ValueError("Wave47 tries to re-ingest a completed external species-trait")
    combined_external = pd.concat(
        [previous_external, new_cells], ignore_index=True, sort=False
    ).fillna("")

    touched = set(new_raw[["axis", "trait_name"]].itertuples(index=False, name=None))
    fixed_touched = direct.loc[
        direct.apply(
            lambda row: (
                (str(row["axis"]), str(row["trait_name"])) in touched
                and str(row["genus"]) == "Buxus"
            ),
            axis=1,
        )
    ].copy()
    external_touched = combined_external.loc[
        combined_external.apply(
            lambda row: (
                (str(row["axis"]), str(row["trait_name"])) in touched
                and str(row["genus"]) == "Buxus"
            ),
            axis=1,
        )
    ].copy()
    previous_external_touched = previous_external.loc[
        previous_external.apply(
            lambda row: (
                (str(row["axis"]), str(row["trait_name"])) in touched
                and str(row["genus"]) == "Buxus"
            ),
            axis=1,
        )
    ].copy()
    touched_cells = pd.concat(
        [fixed_touched, external_touched], ignore_index=True, sort=False
    ).fillna("")
    reconstructed = _reconstruct_lineages(
        pd.concat(
            [fixed_touched, previous_external_touched],
            ignore_index=True,
            sort=False,
        ).fillna("")
    )
    # Use the original newly reviewed row lineage rather than reconstructing it.
    touched_lineages = pd.concat(
        [reconstructed, new_lineages], ignore_index=True, sort=False
    ).fillna("")
    touched_old_low = old_low.loc[
        old_low["genus"].eq("Buxus")
        & old_low.apply(
            lambda row: (str(row["axis"]), str(row["trait_name"])) in touched,
            axis=1,
        )
    ].copy()
    recalculated = build_rule_audit(touched_cells, touched_lineages, touched_old_low)
    if recalculated.empty:
        raise ValueError("Wave47 produced no touched rule audit")

    touched_keys = {("Buxus", axis, trait) for axis, trait in touched}
    keep = ~previous_rules.apply(
        lambda row: (
            (
                str(row["genus"]),
                str(row["axis"]),
                str(row["trait_name"]),
            )
            in touched_keys
        ),
        axis=1,
    )
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
    current_eligible = rules.loc[rules["setting"].eq("current_min3") & _eligible(rules)]
    previous_key_set = set(previous_eligible[key_columns].itertuples(index=False, name=None))
    current_key_set = set(current_eligible[key_columns].itertuples(index=False, name=None))
    newly_eligible_keys = current_key_set - previous_key_set
    expected_key = ("Buxus", "reproductive_assurance", "self_incompatibility")
    if newly_eligible_keys != {expected_key}:
        raise ValueError(f"Wave47 unexpected new rules: {sorted(newly_eligible_keys)}")
    new_rule = current_eligible.loc[
        current_eligible.apply(
            lambda row: tuple(str(row[column]) for column in key_columns) == expected_key,
            axis=1,
        )
    ].copy()

    combined_cells = pd.concat([direct, combined_external], ignore_index=True, sort=False).fillna(
        ""
    )
    rebuilt_low = apply_genus_rules(
        master,
        combined_cells,
        new_rule,
        "current_min3",
    )
    if set(rebuilt_low["axis"]) != {"reproductive_assurance"}:
        raise ValueError("Wave47 crossed a trait or axis boundary")

    combined_conflicts = pd.concat(
        [previous_conflicts, new_cell_audit], ignore_index=True, sort=False
    ).fillna("")
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "resolved_direct_species_trait.csv.gz": direct,
        "rebuilt_all_evidence_validated_low.csv.gz": rebuilt_low,
        "trait_specific_genus_rule_audit.csv.gz": rules,
        "external_congener_resolved_species_trait.csv.gz": combined_external,
        "external_congener_source_lineage_conflicts.csv.gz": combined_conflicts,
    }
    for name, frame in outputs.items():
        _write_gzip_csv(frame, output_dir / name)

    rule_row = new_rule.iloc[0]
    summary: dict[str, Any] = {
        "contract": "wave47_incremental_all_evidence_touched_rule_rebuild_v1",
        "fixed_target_species": expected_species,
        "touched_genus": "Buxus",
        "touched_traits": sorted(trait for _, trait in touched),
        "new_external_rows": len(new_raw),
        "new_external_species_trait": len(new_cells),
        "new_eligible_rules": ["Buxus x self_incompatibility"],
        "new_rule": {
            "n_direct_species": int(rule_row["n_direct_species"]),
            "dominance": float(rule_row["dominance"]),
            "species_loo_accuracy": float(rule_row["species_loo_accuracy"]),
            "lineage_loo_accuracy": float(rule_row["lineage_loo_accuracy"]),
            "inferred_value": str(rule_row["inferred_value"]),
        },
        "new_validated_low_species_trait": len(rebuilt_low),
        "source_lineage_audit": {
            "external_congener_support": {
                "files": 1,
                "rows": len(new_raw),
                "resolved_species_trait": len(new_cells),
                "entered_confirmatory_direct_coverage": 0,
            }
        },
        "checks": {
            "fixed_denominator": True,
            "incremental_scope_only_touched_genus_trait": True,
            "new_source_lineage_not_previously_ingested": True,
            "external_congener_not_confirmatory_direct": True,
            "trait_specific_genus_join": True,
            "minimum_species_three": int(rule_row["n_direct_species"]) >= 3,
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
        new_external_evidence_csv=args.new_external_evidence_csv,
        output_dir=args.output_dir,
        expected_species=args.expected_species,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
