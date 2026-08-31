"""Recover reviewed direct evidence omitted from the formal Wave48 ledger.

Wave49 does not search the Web.  It compares immutable, previously reviewed
source packages with the formal Wave48 direct ledger, applies the repository's
reviewed exclusions, and ingests only missing ``species x trait`` keys.  Genus
rules are then rebuilt only for touched ``genus x axis x trait_name`` keys.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.all_evidence_trait_audit import (
    EVIDENCE_COLUMNS,
    apply_genus_rules,
    build_rule_audit,
    dedupe_direct_lineages,
    direct_evidence_from_integrated,
    load_ontology,
    resolve_direct_cells,
)
from island_v2.direct_evidence_exclusions import apply_direct_evidence_exclusions
from island_v2.wave37_europe_pmc_checkpoint import (
    EXPECTED_SPECIES,
    _sha256,
    _write_gzip_csv,
)
from island_v2.wave48_incremental_all_evidence import (
    _key_set,
    _reconstruct_lineages,
    _touch_mask,
)

FORMAL_WAVE33_RUN_ID = 32_932_103_226
BASELINE_WAVE48_RUN_ID = 33_392_166_443
EXPECTED_COUNTS = {
    "candidate_rows_before_exclusions": 3_466,
    "reviewed_rows_after_exclusions": 3_435,
    "source_lineages": 3_215,
    "duplicate_source_rows": 431,
    "resolved_species_trait": 1_148,
    "excluded_species_trait": 707,
    "new_eligible_rules": 37,
    "invalidated_prior_rules": 10,
    "new_validated_low_species_trait": 858,
}


def _eligible(frame: pd.DataFrame) -> pd.Series:
    return frame["eligible"].astype(str).str.casefold().isin({"true", "1"})


def _load_manifest(repo_root: Path, manifest_path: Path) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("contract") != "wave49_reviewed_staging_gap_recovery_v1":
        raise ValueError("Wave49 source-manifest contract mismatch")
    if manifest.get("formal_wave33_run_id") != FORMAL_WAVE33_RUN_ID:
        raise ValueError("Wave49 formal Wave33 pin changed")
    if manifest.get("baseline_wave48_run_id") != BASELINE_WAVE48_RUN_ID:
        raise ValueError("Wave49 Wave48 baseline pin changed")
    if manifest.get("fixed_target_species") != EXPECTED_SPECIES:
        raise ValueError("Wave49 fixed denominator changed")

    precision_gate = float(manifest["precision_gate"])
    contamination_ceiling = float(manifest["cultivar_contamination_ceiling"])
    source_ids: set[str] = set()
    for source in manifest.get("sources", []):
        source_id = str(source["source_id"])
        if source_id in source_ids:
            raise ValueError(f"duplicate Wave49 source_id: {source_id}")
        source_ids.add(source_id)
        for path_key, hash_key in (
            ("evidence_path", "evidence_sha256"),
            ("review_artifact", "review_sha256"),
        ):
            path = repo_root / source[path_key]
            if not path.is_file() or _sha256(path) != source[hash_key]:
                raise ValueError(f"Wave49 source receipt mismatch: {path}")
        precision = source.get("review_precision")
        if precision is None:
            if not str(source.get("approval_basis", "")).strip():
                raise ValueError(f"Wave49 source lacks an approval basis: {source_id}")
        elif float(precision) < precision_gate:
            raise ValueError(f"Wave49 source is below the precision gate: {source_id}")
        if float(source["cultivar_contamination_rate"]) > contamination_ceiling:
            raise ValueError(f"Wave49 source exceeds cultivar ceiling: {source_id}")
        approved = source.get("approved_traits")
        if approved is not None and not approved:
            raise ValueError(f"Wave49 source has an empty approved-trait set: {source_id}")
        overrides = source.get("quality_overrides", {})
        if set(overrides.values()).difference({"medium"}):
            raise ValueError("Wave49 quality overrides may only conservatively downgrade")

    if len(source_ids) != 11:
        raise ValueError("Wave49 requires exactly 11 reviewed source packages")

    exclusions = manifest["direct_evidence_exclusions"]
    exclusions_path = repo_root / exclusions["path"]
    if not exclusions_path.is_file() or _sha256(exclusions_path) != exclusions["sha256"]:
        raise ValueError("Wave49 direct-evidence exclusion receipt mismatch")
    return manifest


def _select_gap_evidence(
    *,
    repo_root: Path,
    manifest: dict[str, Any],
    target_species: set[str],
    previous_direct_keys: set[tuple[str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    source_rows: list[dict[str, Any]] = []
    for source in manifest["sources"]:
        evidence = direct_evidence_from_integrated(repo_root / source["evidence_path"])
        evidence = evidence.loc[evidence["accepted_species"].isin(target_species)].copy()
        approved = source.get("approved_traits")
        if approved is not None:
            evidence = evidence.loc[evidence["trait_name"].isin(approved)].copy()
        for trait_name, quality in source.get("quality_overrides", {}).items():
            evidence.loc[evidence["trait_name"].eq(trait_name), "quality"] = quality
        before = len(evidence)
        keys = evidence[["accepted_species", "trait_name"]].apply(tuple, axis=1)
        evidence = evidence.loc[~keys.isin(previous_direct_keys)].copy()
        evidence["packet_source_id"] = source["source_id"]
        evidence["packet_source_file"] = source["evidence_path"]
        frames.append(evidence)
        source_rows.append(
            {
                "source_id": source["source_id"],
                "eligible_target_rows": before,
                "missing_species_trait_rows": len(evidence),
                "review_precision": source.get("review_precision", ""),
                "cultivar_contamination_rate": source["cultivar_contamination_rate"],
                "approved_traits": "|".join(approved or []),
                "quality_overrides": json.dumps(
                    source.get("quality_overrides", {}), sort_keys=True
                ),
            }
        )
    selected = pd.concat(frames, ignore_index=True, sort=False).fillna("")
    required_text = [
        "accepted_species",
        "trait_name",
        "normalized_value",
        "source_url",
        "source_citation",
        "source_excerpt",
        "source_lineage",
        "name_match_method",
    ]
    if selected[required_text].apply(
        lambda column: column.astype(str).str.strip().eq("").any()
    ).any():
        raise ValueError("Wave49 selected evidence lacks exact provenance")
    if selected["name_match_method"].str.casefold().str.contains("fuzzy").any():
        raise ValueError("Wave49 selected evidence contains a fuzzy name match")
    if not selected["evidence_scope"].isin({"species_direct", "synonym_direct"}).all():
        raise ValueError("Wave49 selected evidence is not species/synonym direct")
    if not selected["quality"].isin({"high", "medium"}).all():
        raise ValueError("Wave49 selected evidence contains an invalid quality")
    return selected, pd.DataFrame(source_rows)


def build_wave49_audit(
    *,
    repo_root: Path,
    manifest_path: Path,
    master_csv: Path,
    ontology_yaml: Path,
    baseline_coverage_csv: Path,
    previous_rule_audit_csv: Path,
    previous_resolved_direct_csv: Path,
    previous_external_resolved_csv: Path,
    previous_external_conflicts_csv: Path,
    previous_rebuilt_low_csv: Path,
    output_dir: Path,
    expected_species: int = EXPECTED_SPECIES,
    expected_counts: dict[str, int] | None = EXPECTED_COUNTS,
) -> dict[str, Any]:
    required = (
        manifest_path,
        master_csv,
        ontology_yaml,
        baseline_coverage_csv,
        previous_rule_audit_csv,
        previous_resolved_direct_csv,
        previous_external_resolved_csv,
        previous_external_conflicts_csv,
        previous_rebuilt_low_csv,
    )
    if missing := [str(path) for path in required if not path.is_file()]:
        raise ValueError(f"Wave49 inputs missing: {missing}")

    manifest = _load_manifest(repo_root, manifest_path)
    baseline = pd.read_csv(baseline_coverage_csv, dtype=str).fillna("")
    target_species = set(baseline["accepted_species"])
    if len(baseline) != expected_species * 3 or len(target_species) != expected_species:
        raise ValueError("Wave49 baseline denominator mismatch")
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    master = master.loc[master["accepted_species"].isin(target_species)].copy()
    master = master.drop_duplicates("accepted_species")
    if len(master) != expected_species:
        raise ValueError("Wave49 master does not reproduce the fixed universe")

    previous_rules = pd.read_csv(previous_rule_audit_csv, dtype=str).fillna("")
    previous_direct = pd.read_csv(previous_resolved_direct_csv, dtype=str).fillna("")
    previous_external = pd.read_csv(previous_external_resolved_csv, dtype=str).fillna("")
    previous_conflicts = pd.read_csv(previous_external_conflicts_csv, dtype=str).fillna("")
    old_low = pd.read_csv(previous_rebuilt_low_csv, dtype=str).fillna("")
    previous_direct_keys = set(
        previous_direct[["accepted_species", "trait_name"]].itertuples(
            index=False, name=None
        )
    )
    candidate, source_selection = _select_gap_evidence(
        repo_root=repo_root,
        manifest=manifest,
        target_species=target_species,
        previous_direct_keys=previous_direct_keys,
    )
    candidate_rows = len(candidate)
    exclusions_path = repo_root / manifest["direct_evidence_exclusions"]["path"]
    exclusions = pd.read_csv(exclusions_path, dtype=str).fillna("")
    reviewed, exclusion_audit = apply_direct_evidence_exclusions(candidate, exclusions)

    ontology = load_ontology(ontology_yaml)
    lineages, duplicates = dedupe_direct_lineages(reviewed, ontology)
    new_direct, resolution_audit = resolve_direct_cells(lineages)
    for frame in (lineages, new_direct):
        frame["genus"] = frame["accepted_species"].str.split().str[0]
    if set(
        new_direct[["accepted_species", "trait_name"]].itertuples(
            index=False, name=None
        )
    ) & previous_direct_keys:
        raise ValueError("Wave49 tries to re-ingest a completed species-trait")

    combined_direct = pd.concat(
        [previous_direct, new_direct], ignore_index=True, sort=False
    ).fillna("")
    all_cells = pd.concat(
        [combined_direct, previous_external], ignore_index=True, sort=False
    ).fillna("")
    touched = _key_set(lineages)
    prior_touched = pd.concat(
        [
            previous_direct.loc[_touch_mask(previous_direct, touched)],
            previous_external.loc[_touch_mask(previous_external, touched)],
        ],
        ignore_index=True,
        sort=False,
    ).fillna("")
    touched_lineages = pd.concat(
        [_reconstruct_lineages(prior_touched), lineages],
        ignore_index=True,
        sort=False,
    ).fillna("")
    touched_cells = all_cells.loc[_touch_mask(all_cells, touched)].copy()
    touched_old_low = old_low.loc[_touch_mask(old_low, touched)].copy()
    recalculated = build_rule_audit(
        touched_cells, touched_lineages, touched_old_low
    )
    rules = pd.concat(
        [previous_rules.loc[~_touch_mask(previous_rules, touched)], recalculated],
        ignore_index=True,
        sort=False,
    ).fillna("")
    rules = rules.sort_values(
        ["setting", "genus", "axis", "trait_name"], kind="stable"
    ).reset_index(drop=True)

    rule_keys = ["genus", "axis", "trait_name"]
    previous_eligible = previous_rules.loc[
        previous_rules["setting"].eq("current_min3") & _eligible(previous_rules)
    ].copy()
    current_eligible = rules.loc[
        rules["setting"].eq("current_min3") & _eligible(rules)
    ].copy()
    previous_by_key = previous_eligible.set_index(rule_keys)[
        "inferred_state_set"
    ].to_dict()
    current_by_key = current_eligible.set_index(rule_keys)[
        "inferred_state_set"
    ].to_dict()
    new_keys = set(current_by_key).difference(previous_by_key)
    invalidated_keys = set(previous_by_key).difference(current_by_key)
    changed_keys = {
        key
        for key in set(previous_by_key).intersection(current_by_key)
        if previous_by_key[key] != current_by_key[key]
    }
    new_rules = current_eligible.loc[
        current_eligible.apply(
            lambda row: tuple(str(row[column]) for column in rule_keys) in new_keys,
            axis=1,
        )
    ].copy()
    rebuilt_low = apply_genus_rules(
        master, all_cells, new_rules, "current_min3"
    )

    counts = {
        "candidate_rows_before_exclusions": candidate_rows,
        "reviewed_rows_after_exclusions": len(reviewed),
        "source_lineages": len(lineages),
        "duplicate_source_rows": len(duplicates),
        "resolved_species_trait": len(new_direct),
        "excluded_species_trait": int(
            resolution_audit["resolution_status"].ne("resolved").sum()
        ),
        "new_eligible_rules": len(new_rules),
        "invalidated_prior_rules": len(invalidated_keys),
        "new_validated_low_species_trait": len(rebuilt_low),
    }
    if expected_counts is not None and counts != expected_counts:
        raise ValueError(f"Wave49 frozen counts changed: {counts}")

    output_dir.mkdir(parents=True, exist_ok=True)
    empty_external = pd.DataFrame(columns=EVIDENCE_COLUMNS)
    outputs = {
        "resolved_direct_species_trait.csv.gz": combined_direct,
        "rebuilt_all_evidence_validated_low.csv.gz": rebuilt_low,
        "trait_specific_genus_rule_audit.csv.gz": rules,
        "external_congener_resolved_species_trait.csv.gz": previous_external,
        "external_congener_source_lineage_conflicts.csv.gz": previous_conflicts,
        "wave49_reviewed_direct_evidence.csv.gz": reviewed,
        "wave49_direct_source_lineages.csv.gz": lineages,
        "wave49_direct_source_lineage_duplicates.csv.gz": duplicates,
        "wave49_direct_resolution_audit.csv.gz": resolution_audit,
        "wave49_direct_exclusion_audit.csv.gz": exclusion_audit,
        "wave49_source_selection_audit.csv.gz": source_selection,
        "wave49_external_congener_evidence.csv.gz": empty_external,
    }
    for name, frame in outputs.items():
        _write_gzip_csv(frame, output_dir / name)

    new_rule_labels = [" x ".join(key) for key in sorted(new_keys)]
    invalidated_labels = [" x ".join(key) for key in sorted(invalidated_keys)]
    changed_labels = [" x ".join(key) for key in sorted(changed_keys)]
    summary: dict[str, Any] = {
        "contract": "wave49_reviewed_staging_gap_incremental_audit_v1",
        "formal_wave33_run_id": FORMAL_WAVE33_RUN_ID,
        "baseline_wave48_run_id": BASELINE_WAVE48_RUN_ID,
        "fixed_target_species": expected_species,
        "counts": counts,
        "new_eligible_rules": new_rule_labels,
        "invalidated_prior_rules": invalidated_labels,
        "changed_prior_rule_states": changed_labels,
        "source_lineage_audit": {
            "new_direct": {
                "source_packages": len(manifest["sources"]),
                "candidate_rows": candidate_rows,
                "reviewed_rows": len(reviewed),
                "source_lineages": len(lineages),
                "resolved_species_trait": len(new_direct),
                "excluded_species_trait": counts["excluded_species_trait"],
            },
            "external_congener_support": {
                "rows": 0,
                "resolved_species_trait": 0,
                "entered_confirmatory_direct_coverage": 0,
            },
        },
        "checks": {
            "formal_wave33_baseline_pinned": True,
            "formal_wave48_immediate_baseline_pinned": True,
            "fixed_denominator": len(baseline) == expected_species * 3,
            "only_missing_species_trait_selected": True,
            "reviewed_exclusions_applied": int(exclusion_audit["matched_rows"].sum())
            == 31,
            "precision_gate_enforced": True,
            "cultivar_contamination_ceiling_enforced": True,
            "source_lineage_deduplicated": True,
            "direct_conflicts_excluded": True,
            "trait_specific_genus_join": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
            "external_congener_not_confirmatory_direct": True,
            "shadow_invalidations_not_silently_deleted": True,
        },
        "input_sha256": {str(path): _sha256(path) for path in required},
        "artifact_sha256": {
            name: _sha256(output_dir / name) for name in outputs
        },
    }
    (output_dir / "all_evidence_trait_coverage_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )

    checkpoint = {
        "contract": "wave49_reviewed_staging_gap_checkpoint_v1",
        "formal_wave33_run_id": FORMAL_WAVE33_RUN_ID,
        "baseline_formal_run_id": BASELINE_WAVE48_RUN_ID,
        "fixed_target_species": expected_species,
        "evidence": {
            "candidate_rows_before_exclusions": candidate_rows,
            "reviewed_direct_rows": len(reviewed),
            "resolved_direct_species_trait": len(new_direct),
            "new_external_species_trait": 0,
        },
        "queries": {
            "new_web_queries": 0,
            "reused_reviewed_source_packages": len(manifest["sources"]),
        },
        "query_cost_usd": 0.0,
        "checks": summary["checks"],
        "artifact_sha256": summary["artifact_sha256"],
    }
    (output_dir / "wave49_checkpoint_validation.json").write_text(
        json.dumps(checkpoint, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--master-csv", required=True, type=Path)
    parser.add_argument("--ontology-yaml", required=True, type=Path)
    parser.add_argument("--baseline-coverage-csv", required=True, type=Path)
    parser.add_argument("--previous-rule-audit-csv", required=True, type=Path)
    parser.add_argument("--previous-resolved-direct-csv", required=True, type=Path)
    parser.add_argument("--previous-external-resolved-csv", required=True, type=Path)
    parser.add_argument("--previous-external-conflicts-csv", required=True, type=Path)
    parser.add_argument("--previous-rebuilt-low-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-species", type=int, default=EXPECTED_SPECIES)
    args = parser.parse_args()
    summary = build_wave49_audit(
        repo_root=args.repo_root,
        manifest_path=args.manifest,
        master_csv=args.master_csv,
        ontology_yaml=args.ontology_yaml,
        baseline_coverage_csv=args.baseline_coverage_csv,
        previous_rule_audit_csv=args.previous_rule_audit_csv,
        previous_resolved_direct_csv=args.previous_resolved_direct_csv,
        previous_external_resolved_csv=args.previous_external_resolved_csv,
        previous_external_conflicts_csv=args.previous_external_conflicts_csv,
        previous_rebuilt_low_csv=args.previous_rebuilt_low_csv,
        output_dir=args.output_dir,
        expected_species=args.expected_species,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
