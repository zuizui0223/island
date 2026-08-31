"""Validate the frozen Wave43 external AusTraits morphology packet."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

EXPECTED_NAMES = 2_580
EXPECTED_TARGET_SPECIES = 106_295
EXPECTED_MAPPING_REASONS = {
    "accepted_strict_two_backbone_external": 2_188,
    "backbone_family_conflict": 28,
    "backbones_disagree": 119,
    "gbif_not_exact": 17,
    "mapped_into_fixed_target_universe": 83,
    "source_wfo_family_conflict": 47,
    "wfo_ambiguous_exact_name": 31,
    "wfo_no_accepted_species_family_route": 42,
    "wfo_no_exact_species_name": 25,
}
EXPECTED_TRAIT_COUNTS = {
    "flower_primary_color": 2_764,
    "floral_symmetry": 228,
}
EXPECTED_AXIS_BY_TRAIT = {
    "flower_primary_color": "flower_colour",
    "floral_symmetry": "floral_structural_complexity",
}
EXPECTED_SHA256 = {
    "archive": "2995fb9e5eebf9271f9241b5a6c2e2cbb7b45e8eda5ce4db43c14c68f95e2f3f",
    "evidence": "aaefd815251fc6025a1accfeb643db533037a647fd2efce1c22207240d1b5c79",
    "mapping": "6550672944a1cba7aafe2eba46b047e364b58f30c729bc15ca5c629fd477f664",
    "selection": "55577de4b15a6d5e9cba94523211395bbae68ffed211262df2281d013360d9a8",
}
CONTRACT = "external_congener_species_direct_strict_two_backbone_v1"
BINOMIAL = re.compile(r"[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _counts(series: pd.Series) -> dict[str, int]:
    return {str(key): int(value) for key, value in series.value_counts().items()}


def _read_jsonl_gz(path: Path) -> list[dict[str, Any]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def validate_mapping(
    mapping: pd.DataFrame,
    *,
    target_species: set[str],
    expected_names: int = EXPECTED_NAMES,
) -> None:
    required = {
        "source_name",
        "source_family",
        "accepted_species",
        "mapping_reason",
        "wfo_accepted_species",
        "wfo_family",
        "wfo_release",
        "gbif_accepted_species",
        "gbif_family",
        "gbif_match_type",
        "gbif_rank",
        "gbif_kingdom",
    }
    if missing := required.difference(mapping.columns):
        raise ValueError(f"Wave43 mapping is missing columns: {sorted(missing)}")
    if len(mapping) != expected_names or mapping["source_name"].duplicated().any():
        raise ValueError("Wave43 mapping does not contain one row per source name")
    accepted = mapping.loc[
        mapping["mapping_reason"].eq("accepted_strict_two_backbone_external")
    ]
    checks = (
        accepted["accepted_species"].map(
            lambda value: bool(BINOMIAL.fullmatch(value))
        ),
        accepted["accepted_species"].eq(accepted["wfo_accepted_species"]),
        accepted["accepted_species"].eq(accepted["gbif_accepted_species"]),
        accepted["source_family"].eq(accepted["wfo_family"]),
        accepted["wfo_family"].eq(accepted["gbif_family"]),
        accepted["wfo_release"].eq("2026-06"),
        accepted["gbif_match_type"].eq("EXACT"),
        accepted["gbif_rank"].eq("SPECIES"),
        accepted["gbif_kingdom"].eq("Plantae"),
        ~accepted["accepted_species"].isin(target_species),
    )
    if not all(bool(check.all()) for check in checks):
        raise ValueError("Wave43 accepted mapping violates the strict identity gate")
    rejected = mapping.loc[
        ~mapping["mapping_reason"].eq("accepted_strict_two_backbone_external")
    ]
    if rejected["accepted_species"].ne("").any():
        raise ValueError("Wave43 rejected mappings retained an accepted species")


def validate_packet(
    *,
    packet_dir: Path,
    target_coverage_csv: Path,
    output_json: Path,
    expected_species: int = EXPECTED_TARGET_SPECIES,
) -> dict[str, Any]:
    paths = {
        "source_manifest": packet_dir / "source_manifest.json",
        "candidate_summary": packet_dir
        / "wave43_austraits_external_candidate_summary.json",
        "manifest": packet_dir / "wave43_austraits_external_name_manifest.csv",
        "source_rows": packet_dir
        / "wave43_austraits_external_source_rows.csv.gz",
        "wfo_mapping": packet_dir / "wave43_wfo_2026_06_local_mapping.csv.gz",
        "wfo_summary": packet_dir
        / "wave43_wfo_2026_06_local_mapping_summary.json",
        "gbif_snapshot": packet_dir / "wave43_gbif_response_snapshot.jsonl.gz",
        "gbif_summary": packet_dir / "wave43_gbif_query_summary.json",
        "mapping": packet_dir
        / "wave43_austraits_two_backbone_mapping_audit.csv.gz",
        "evidence": packet_dir
        / "wave43_austraits_external_morphology_evidence.csv.gz",
        "selection": packet_dir
        / "wave43_austraits_external_selection_audit.csv.gz",
        "evidence_summary": packet_dir
        / "wave43_austraits_external_evidence_summary.json",
    }
    if missing := [str(path) for path in paths.values() if not path.is_file()]:
        raise ValueError(f"Wave43 packet is incomplete: {missing}")

    target = pd.read_csv(
        target_coverage_csv,
        usecols=["accepted_species", "axis"],
        dtype=str,
    ).fillna("")
    target_species = set(target["accepted_species"])
    if len(target) != expected_species * 3 or len(target_species) != expected_species:
        raise ValueError("Wave43 target coverage denominator mismatch")

    manifest = pd.read_csv(paths["manifest"], dtype=str).fillna("")
    local_wfo = pd.read_csv(paths["wfo_mapping"], dtype=str).fillna("")
    mapping = pd.read_csv(paths["mapping"], dtype=str).fillna("")
    evidence = pd.read_csv(paths["evidence"], dtype=str).fillna("")
    selection = pd.read_csv(paths["selection"], dtype=str).fillna("")
    gbif = _read_jsonl_gz(paths["gbif_snapshot"])

    frames = (manifest, local_wfo, mapping)
    if any(len(frame) != EXPECTED_NAMES for frame in frames):
        raise ValueError("Wave43 name tables have an unexpected row count")
    names = set(manifest["source_name"])
    if len(names) != EXPECTED_NAMES:
        raise ValueError("Wave43 manifest has duplicate source names")
    if set(local_wfo["source_name"]) != names or set(mapping["source_name"]) != names:
        raise ValueError("Wave43 name tables do not cover the same source names")
    gbif_names = [str(record.get("source_name", "")) for record in gbif]
    if len(gbif_names) != EXPECTED_NAMES or set(gbif_names) != names:
        raise ValueError("Wave43 GBIF snapshot does not cover the manifest exactly")

    validate_mapping(mapping, target_species=target_species)
    if _counts(mapping["mapping_reason"]) != EXPECTED_MAPPING_REASONS:
        raise ValueError("Wave43 mapping reason counts changed")

    required_evidence = {
        "accepted_species",
        "axis",
        "trait_name",
        "normalized_value",
        "quality",
        "source_provider",
        "source_url",
        "source_record_id",
        "source_citation",
        "source_excerpt",
        "evidence_scope",
        "name_match_method",
        "source_lineage",
        "lineage_method",
        "source_run_id",
        "source_artifact",
        "source_file",
        "acceptance_contract",
    }
    if missing := required_evidence.difference(evidence.columns):
        raise ValueError(f"Wave43 evidence is missing columns: {sorted(missing)}")
    if evidence[list(required_evidence)].apply(
        lambda column: column.astype(str).str.strip().eq("").any()
    ).any():
        raise ValueError("Wave43 evidence has incomplete provenance")
    if len(evidence) != 2_992 or evidence["source_record_id"].duplicated().any():
        raise ValueError("Wave43 evidence row or record-ID count changed")
    accepted_species = set(
        mapping.loc[
            mapping["mapping_reason"].eq("accepted_strict_two_backbone_external"),
            "accepted_species",
        ]
    )
    axis_matches_trait = pd.Series(
        [
            axis == EXPECTED_AXIS_BY_TRAIT.get(trait, "")
            for trait, axis in zip(
                evidence["trait_name"], evidence["axis"], strict=True
            )
        ],
        index=evidence.index,
    )
    evidence_checks = (
        evidence["accepted_species"].isin(accepted_species),
        ~evidence["accepted_species"].isin(target_species),
        axis_matches_trait,
        evidence["quality"].eq("medium"),
        evidence["evidence_scope"].eq("external_congener_species_direct"),
        evidence["name_match_method"].eq("strict_wfo_gbif_two_backbone"),
        evidence["lineage_method"].eq(
            "normalized_underlying_source_citation_fingerprint"
        ),
        evidence["source_lineage"].str.startswith("austraits-original-source:"),
        evidence["acceptance_contract"].eq(CONTRACT),
    )
    if not all(bool(check.all()) for check in evidence_checks):
        raise ValueError("Wave43 evidence violates its fail-closed contract")
    if _counts(evidence["trait_name"]) != EXPECTED_TRAIT_COUNTS:
        raise ValueError("Wave43 trait evidence counts changed")
    if evidence["source_lineage"].nunique() != 5:
        raise ValueError("Wave43 underlying-source lineage count changed")

    selected = selection.loc[selection["selection_reason"].eq("selected")].copy()
    if len(selected) != len(evidence):
        raise ValueError("Wave43 selected rows do not match evidence rows")
    evidence = evidence.copy()
    evidence["selection_record_id"] = (
        evidence["source_record_id"]
        .str.removeprefix("wave43-external:")
        .str.rsplit(":", n=1)
        .str[0]
    )
    join_keys = [
        "accepted_species",
        "trait_name",
        "normalized_value",
        "selection_record_id",
    ]
    selected = selected.rename(columns={"source_record_id": "selection_record_id"})
    if evidence[join_keys].value_counts().to_dict() != selected[join_keys].value_counts().to_dict():
        raise ValueError("Wave43 selection and evidence rows disagree")

    declared = json.loads(paths["evidence_summary"].read_text(encoding="utf-8"))
    if declared["sha256"] != EXPECTED_SHA256:
        raise ValueError("Wave43 declared source hashes changed")
    if _sha256(paths["mapping"]) != EXPECTED_SHA256["mapping"]:
        raise ValueError("Wave43 mapping hash changed")
    if _sha256(paths["evidence"]) != EXPECTED_SHA256["evidence"]:
        raise ValueError("Wave43 evidence hash changed")
    if _sha256(paths["selection"]) != EXPECTED_SHA256["selection"]:
        raise ValueError("Wave43 selection hash changed")

    queries = json.loads(paths["gbif_summary"].read_text(encoding="utf-8"))
    summary = {
        "contract": "wave43_austraits_external_checkpoint_validation_v1",
        "fixed_target_species": expected_species,
        "queries": queries,
        "query_cost_usd": 0,
        "mapping_reason_counts": _counts(mapping["mapping_reason"]),
        "evidence": {
            "rows": len(evidence),
            "species": int(evidence["accepted_species"].nunique()),
            "species_trait": int(
                evidence[["accepted_species", "trait_name"]]
                .drop_duplicates()
                .shape[0]
            ),
            "by_trait": _counts(evidence["trait_name"]),
            "entered_fixed_target_direct_coverage": 0,
            "underlying_source_lineages": int(evidence["source_lineage"].nunique()),
        },
        "checks": {
            "fixed_denominator": True,
            "exact_wfo_2026_06_and_gbif": True,
            "external_species_only": True,
            "trait_specific_rows": True,
            "complete_provenance": True,
            "underlying_source_lineage_deduplicated": True,
            "strict_axes_unchanged": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
        },
        "input_sha256": {key: _sha256(path) for key, path in paths.items()},
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packet-dir", required=True, type=Path)
    parser.add_argument("--target-coverage-csv", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument(
        "--expected-species", type=int, default=EXPECTED_TARGET_SPECIES
    )
    args = parser.parse_args()
    summary = validate_packet(
        packet_dir=args.packet_dir,
        target_coverage_csv=args.target_coverage_csv,
        output_json=args.output_json,
        expected_species=args.expected_species,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
