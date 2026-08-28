"""Lossless Wave38 overlay and trait-specific touched-rule rebuild.

Wave38 promotes only manually audited, source-locked LLM extractions.  The LLM
is never a source: every retained row has an exact quote from a frozen page.
This module reconnects those direct rows to the shared ``genus x trait_name``
Validated Low implementation and applies only new resolved evidence to the
formal Wave37 coverage baseline.
"""

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
    load_ontology,
    load_reviewed_direct_supplements,
    resolve_direct_cells,
)
from island_v2.wave35_trait_overlay import (
    AXES,
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
    EVIDENCE_COLUMNS,
    EXPECTED_SPECIES,
    _coverage_summary,
    _eligible,
    _sha256,
    _text,
    _validate_coverage,
    _write_gzip_csv,
)
from island_v2.wave38_source_locked_llm_checkpoint import SOURCE_GROUP


def _touched_keys(reviewed: pd.DataFrame) -> set[tuple[str, str]]:
    return {
        (str(row.accepted_species).split()[0], str(row.trait_name))
        for row in reviewed.itertuples(index=False)
    }


def _reconstruct_formal_lineages(
    formal_resolved_direct_csv: Path,
    touched_keys: set[tuple[str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    formal = pd.read_csv(formal_resolved_direct_csv, dtype=str).fillna("")
    formal["genus"] = formal["accepted_species"].str.split().str[0]
    formal = formal.loc[
        formal.apply(
            lambda row: (str(row["genus"]), str(row["trait_name"])) in touched_keys,
            axis=1,
        )
        & formal["resolution_status"].eq("resolved")
    ].copy()
    reconstructed: list[dict[str, str]] = []
    for row in formal.to_dict("records"):
        lineages = _split_pipe(row.get("source_lineages", ""))
        if not lineages:
            raise ValueError("formal resolved direct cell lacks source lineage")
        for lineage in sorted(lineages):
            reconstructed.append(
                {
                    "accepted_species": _text(row.get("accepted_species")),
                    "axis": _text(row.get("axis")),
                    "trait_name": _text(row.get("trait_name")),
                    "normalized_value": _text(row.get("state_set")),
                    "quality": _text(row.get("quality")).casefold(),
                    "source_group": _text(row.get("source_groups")),
                    "source_provider": "formal_resolved_direct_lineage",
                    "source_url": "artifact://formal-wave33-resolved-direct",
                    "source_record_id": lineage,
                    "source_citation": "formal resolved direct lineage",
                    "source_excerpt": "",
                    "evidence_scope": "species_direct",
                    "name_match_method": "formal_resolved_species_identity",
                    "source_lineage": lineage,
                    "lineage_method": "expanded_from_formal_resolved_cell",
                    "source_run_id": "32932103226",
                    "source_artifact": "integrated-trait-coverage-32932103226",
                    "source_file": str(formal_resolved_direct_csv),
                    "acceptance_contract": "formal_resolved_direct_lineage_reuse_v1",
                }
            )
    return formal, pd.DataFrame(reconstructed, columns=EVIDENCE_COLUMNS)


def build_touched_rule_rebuild(
    *,
    formal_resolved_direct_csv: Path,
    reviewed_direct_csv: Path,
    supplemental_direct_evidence_csvs: tuple[Path, ...],
    master_csv: Path,
    ontology_yaml: Path,
    output_dir: Path,
    source_group: str = SOURCE_GROUP,
    wave_label: str = "wave38",
) -> dict[str, Any]:
    """Rebuild all Wave38-touched genus x trait rules with shared logic."""
    reviewed = pd.read_csv(reviewed_direct_csv, dtype=str).fillna("")
    touched_keys = _touched_keys(reviewed)
    formal, formal_lineages = _reconstruct_formal_lineages(
        formal_resolved_direct_csv, touched_keys
    )
    supplemental = load_reviewed_direct_supplements(supplemental_direct_evidence_csvs)
    supplemental["genus"] = supplemental["accepted_species"].str.split().str[0]
    supplemental = supplemental.loc[
        supplemental.apply(
            lambda row: (str(row["genus"]), str(row["trait_name"])) in touched_keys,
            axis=1,
        )
    ].drop(columns="genus")
    raw = pd.concat([formal_lineages, supplemental], ignore_index=True).fillna("")
    ontology = load_ontology(ontology_yaml)
    lineages, duplicates = dedupe_direct_lineages(raw, ontology)
    direct_cells, conflicts = resolve_direct_cells(lineages)
    direct_cells["genus"] = direct_cells["accepted_species"].str.split().str[0]
    lineages["genus"] = lineages["accepted_species"].str.split().str[0]
    old_low = pd.DataFrame(columns=["genus", "trait_name", "state_set"])
    rules = build_rule_audit(direct_cells, lineages, old_low)
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    rebuilt_low = apply_genus_rules(master, direct_cells, rules, "current_min3")
    source_touched = direct_cells["source_groups"].map(
        lambda value: source_group in _split_pipe(value)
    )
    wave_direct = direct_cells.loc[source_touched].copy()
    if wave_direct.empty:
        raise ValueError(f"{wave_label} produced no resolved direct cells")

    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "resolved_direct_species_trait.csv.gz": wave_direct,
        "rebuilt_all_evidence_validated_low.csv.gz": rebuilt_low,
        "trait_specific_genus_rule_audit.csv.gz": rules,
        f"{wave_label}_touched_source_lineages.csv.gz": lineages,
        f"{wave_label}_touched_source_lineage_duplicates.csv.gz": duplicates,
        f"{wave_label}_touched_source_lineage_conflicts.csv.gz": conflicts,
    }
    for name, frame in outputs.items():
        _write_gzip_csv(frame, output_dir / name)
    summary = {
        "contract": f"{wave_label}_common_trait_specific_touched_rebuild_v1",
        "formal_resolved_cells_reused": len(formal),
        "formal_source_lineages_reconstructed": len(formal_lineages),
        "supplemental_touched_rows": len(supplemental),
        "touched_genus_trait": len(touched_keys),
        f"resolved_{wave_label}_direct_species_trait": len(wave_direct),
        "candidate_low_before_baseline_filter": len(rebuilt_low),
        "checks": {
            "shared_trait_specific_rule_builder": True,
            "genus_x_trait_name_join": True,
            "lineage_leave_one_out_required": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
        },
        "input_sha256": {
            formal_resolved_direct_csv.name: _sha256(formal_resolved_direct_csv),
            reviewed_direct_csv.name: _sha256(reviewed_direct_csv),
            ontology_yaml.name: _sha256(ontology_yaml),
            **{path.name: _sha256(path) for path in supplemental_direct_evidence_csvs},
        },
        "artifact_sha256": {name: _sha256(output_dir / name) for name in outputs},
    }
    (output_dir / f"{wave_label}_touched_rebuild_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def _prior_rule_table(paths: tuple[Path, ...]) -> pd.DataFrame:
    frames = [pd.read_csv(path, dtype=str).fillna("") for path in paths]
    prior = pd.concat(frames, ignore_index=True)
    prior = prior.loc[prior["setting"].eq("current_min3") & _eligible(prior)].copy()
    keys = ["genus", "axis", "trait_name"]
    contradictions = (
        prior.groupby(keys)["inferred_state_set"].nunique().loc[lambda values: values > 1]
    )
    if len(contradictions):
        raise ValueError("prior eligible rule manifests disagree on inferred state")
    return prior.drop_duplicates(keys, keep="last")


def build_formal_overlay(
    *,
    baseline_csv: Path,
    all_evidence_dir: Path,
    reviewed_direct_csv: Path,
    checkpoint_summary_json: Path,
    previous_rule_audit_csvs: tuple[Path, ...],
    output_dir: Path,
    expected_species: int = EXPECTED_SPECIES,
) -> dict[str, Any]:
    """Apply reviewed Wave38 direct evidence and newly eligible Low losslessly."""
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
        raise ValueError(f"Wave38 overlay inputs missing: {missing}")

    baseline = pd.read_csv(baseline_csv, dtype=str).fillna("")
    _validate_coverage(baseline, expected_species)
    reviewed = pd.read_csv(reviewed_direct_csv, dtype=str).fillna("")
    resolved = pd.read_csv(resolved_path, dtype=str).fillna("")
    rebuilt_low = pd.read_csv(rebuilt_low_path, dtype=str).fillna("")
    rules = pd.read_csv(rules_path, dtype=str).fillna("")
    prior_rules = _prior_rule_table(previous_rule_audit_csvs)
    touched_keys = _touched_keys(reviewed)

    source_touched = resolved["source_groups"].map(
        lambda value: SOURCE_GROUP in _split_pipe(value)
    )
    direct = _validate_direct(resolved.loc[source_touched].copy())
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
    if invalidated_keys or changed_keys:
        raise ValueError(
            "Wave38 would leave stale prior Low rules: "
            f"invalidated={invalidated_keys}, changed={changed_keys}"
        )
    new_keys = set(current_by_key) - set(prior_by_key)
    new_rules = current_rules.loc[
        current_rules.apply(
            lambda row: (
                str(row["genus"]), str(row["axis"]), str(row["trait_name"])
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
            raise ValueError(f"direct evidence outside denominator: {key}")
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
    if not low.empty:
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
                    "source_groups": row.source_groups,
                }
            )

    result = result.reset_index().sort_values(["accepted_species", "axis"])
    _validate_coverage(result, expected_species)
    before = _coverage_summary(baseline, expected_species)
    after = _coverage_summary(result, expected_species)
    changes_frame = pd.DataFrame(changes)
    if not changes_frame.empty:
        changes_frame = changes_frame.sort_values(["accepted_species", "axis", "action"])
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
        raise ValueError(f"Wave38 coverage loss must be zero, observed {loss}")
    if (
        comparison["quality_after"].map(QUALITY_RANK)
        < comparison["quality_before"].map(QUALITY_RANK)
    ).any():
        raise ValueError("Wave38 downgraded an existing quality")

    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "wave38_species_axis_coverage.csv.gz": result,
        "wave38_resolved_direct_species_trait.csv.gz": direct,
        "wave38_candidate_validated_low_species_trait.csv.gz": low,
        "wave38_provider_touched_new_rule_audit.csv.gz": new_rules,
        "wave38_change_audit.csv.gz": changes_frame,
    }
    for name, frame in outputs.items():
        _write_gzip_csv(frame, output_dir / name)
    summary = {
        "contract": "wave38_lossless_source_locked_llm_overlay_v1",
        "formal_wave33_run_id": 32932103226,
        "formal_wave37_run_id": 33143109604,
        "wave37_before": before,
        "wave38_after": after,
        "delta": {
            "gross_gain_species_axis": gain,
            "loss_species_axis": loss,
            "net_gain_species_axis": gain,
            "by_axis_net_gain": {
                axis: after["by_axis"][axis]["filled_species"]
                - before["by_axis"][axis]["filled_species"]
                for axis in AXES
            },
            "action_counts": {
                str(key): int(value)
                for key, value in changes_frame["action"].value_counts().items()
            }
            if not changes_frame.empty
            else {},
            "reviewed_direct_rows": len(reviewed),
            "resolved_direct_species_trait": len(direct),
            "new_validated_low_species_trait": len(low),
        },
        "new_eligible_rules": [
            f"{row.genus} x {row.trait_name}" for row in new_rules.itertuples(index=False)
        ],
        "invalidated_prior_rules": [],
        "changed_prior_rule_states": [],
        "checks": {
            "fixed_denominator": True,
            "quality_precedence_high_medium_low": True,
            "direct_conflicts_excluded": True,
            "trait_specific_genus_join": True,
            "lineage_leave_one_out_required": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
            "baseline_loss_zero": loss == 0,
        },
        "checkpoint": json.loads(checkpoint_summary_json.read_text(encoding="utf-8")),
        "input_sha256": {path.name: _sha256(path) for path in required},
    }
    summary["artifact_sha256"] = {name: _sha256(output_dir / name) for name in outputs}
    (output_dir / "wave38_coverage_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return summary


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
    summary = build_touched_rule_rebuild(
        formal_resolved_direct_csv=args.formal_resolved_direct_csv,
        reviewed_direct_csv=args.reviewed_direct_csv,
        supplemental_direct_evidence_csvs=tuple(args.supplemental_direct_evidence_csv),
        master_csv=args.master_csv,
        ontology_yaml=args.ontology_yaml,
        output_dir=args.output_dir,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


def overlay_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-csv", required=True, type=Path)
    parser.add_argument("--all-evidence-dir", required=True, type=Path)
    parser.add_argument("--reviewed-direct-csv", required=True, type=Path)
    parser.add_argument("--checkpoint-summary-json", required=True, type=Path)
    parser.add_argument("--previous-rule-audit-csv", action="append", default=[], type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-species", type=int, default=EXPECTED_SPECIES)
    args = parser.parse_args()
    summary = build_formal_overlay(
        baseline_csv=args.baseline_csv,
        all_evidence_dir=args.all_evidence_dir,
        reviewed_direct_csv=args.reviewed_direct_csv,
        checkpoint_summary_json=args.checkpoint_summary_json,
        previous_rule_audit_csvs=tuple(args.previous_rule_audit_csv),
        output_dir=args.output_dir,
        expected_species=args.expected_species,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
