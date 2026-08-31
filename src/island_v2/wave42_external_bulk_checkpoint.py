"""Validate the frozen Wave42 external reproductive-evidence packet.

Wave42 uses exact WFO 2026-06 and GBIF species matches for species outside the
fixed island denominator.  The resulting rows may train trait-specific genus
rules, but may never enter species-direct coverage for the target universe.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

EXPECTED_NAMES = 1_671
EXPECTED_TARGET_SPECIES = 106_295
EXPECTED_MAPPING_REASONS = {
    "accepted_strict_two_backbone_external": 1_314,
    "backbone_family_conflict": 19,
    "backbones_disagree": 59,
    "gbif_cache_conflict": 25,
    "gbif_not_exact": 46,
    "mapped_into_fixed_target_universe": 54,
    "wfo_ambiguous_exact_name": 81,
    "wfo_no_accepted_species_family_route": 26,
    "wfo_no_exact_species_name": 47,
}
EXPECTED_PROVIDER_COUNTS = {
    "Ferrer et al. 2024 global self-incompatibility database": 650,
    "Goodwillie, Kalisz & Eckert mixed-mating database": 1,
    "Meyer, Galloway & Eckert 2026 Dryad synthesis": 669,
    "Razanajatovo et al. 2016 Nature Communications dataset": 70,
}
EXPECTED_TRAIT_COUNTS = {
    "autonomous_selfing_capacity": 37,
    "cleistogamy": 1,
    "self_incompatibility": 1_352,
}
SHARED_LINEAGE = (
    "wave42-secondary-reproductive-databases:shared-redistribution-guard"
)
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
        raise ValueError(f"Wave42 mapping is missing columns: {sorted(missing)}")
    if len(mapping) != expected_names or mapping["source_name"].duplicated().any():
        raise ValueError("Wave42 mapping does not contain one row per source name")
    accepted = mapping.loc[
        mapping["mapping_reason"].eq("accepted_strict_two_backbone_external")
    ]
    checks = (
        accepted["accepted_species"].map(lambda value: bool(BINOMIAL.fullmatch(value))),
        accepted["accepted_species"].eq(accepted["wfo_accepted_species"]),
        accepted["accepted_species"].eq(accepted["gbif_accepted_species"]),
        accepted["wfo_family"].eq(accepted["gbif_family"]),
        accepted["wfo_release"].eq("2026-06"),
        accepted["gbif_match_type"].eq("EXACT"),
        accepted["gbif_rank"].eq("SPECIES"),
        accepted["gbif_kingdom"].eq("Plantae"),
        ~accepted["accepted_species"].isin(target_species),
    )
    if not all(bool(check.all()) for check in checks):
        raise ValueError("Wave42 accepted mapping violates the strict identity gate")
    rejected = mapping.loc[
        ~mapping["mapping_reason"].eq("accepted_strict_two_backbone_external")
    ]
    if rejected["accepted_species"].ne("").any():
        raise ValueError("Wave42 rejected mappings retained an accepted species")


def validate_packet(
    *,
    packet_dir: Path,
    target_coverage_csv: Path,
    output_json: Path,
    expected_species: int = EXPECTED_TARGET_SPECIES,
) -> dict[str, Any]:
    paths = {
        "source_manifest": packet_dir / "source_manifest.json",
        "manifest": packet_dir / "wave42_new_name_manifest.csv",
        "wfo_mapping": packet_dir / "wave42_wfo_2026_06_local_mapping.csv.gz",
        "wfo_summary": packet_dir
        / "wave42_wfo_2026_06_local_mapping_summary.json",
        "gbif_snapshot": packet_dir / "wave42_gbif_response_snapshot.jsonl.gz",
        "mapping": packet_dir / "wave42_two_backbone_mapping_audit.csv.gz",
        "mapping_summary": packet_dir / "wave42_two_backbone_mapping_summary.json",
        "evidence": packet_dir / "wave42_external_reproductive_evidence.csv.gz",
        "selection": packet_dir / "wave42_external_selection_audit.csv.gz",
        "evidence_summary": packet_dir / "wave42_external_evidence_summary.json",
    }
    if missing := [str(path) for path in paths.values() if not path.is_file()]:
        raise ValueError(f"Wave42 packet is incomplete: {missing}")

    target = pd.read_csv(
        target_coverage_csv,
        usecols=["accepted_species", "axis"],
        dtype=str,
    ).fillna("")
    target_species = set(target["accepted_species"])
    if len(target) != expected_species * 3 or len(target_species) != expected_species:
        raise ValueError("Wave42 target coverage denominator mismatch")

    manifest = pd.read_csv(paths["manifest"], dtype=str).fillna("")
    local_wfo = pd.read_csv(paths["wfo_mapping"], dtype=str).fillna("")
    mapping = pd.read_csv(paths["mapping"], dtype=str).fillna("")
    evidence = pd.read_csv(paths["evidence"], dtype=str).fillna("")
    selection = pd.read_csv(paths["selection"], dtype=str).fillna("")
    gbif = _read_jsonl_gz(paths["gbif_snapshot"])

    frames = (manifest, local_wfo, mapping)
    if any(len(frame) != EXPECTED_NAMES for frame in frames):
        raise ValueError("Wave42 name tables have an unexpected row count")
    names = set(manifest["source_name"])
    if len(names) != EXPECTED_NAMES:
        raise ValueError("Wave42 manifest has duplicate source names")
    if set(local_wfo["source_name"]) != names or set(mapping["source_name"]) != names:
        raise ValueError("Wave42 name tables do not cover the same source names")
    gbif_names = [str(record.get("source_name", "")) for record in gbif]
    if len(gbif_names) != EXPECTED_NAMES or set(gbif_names) != names:
        raise ValueError("Wave42 GBIF snapshot does not cover the manifest exactly")

    validate_mapping(mapping, target_species=target_species)
    if _counts(mapping["mapping_reason"]) != EXPECTED_MAPPING_REASONS:
        raise ValueError("Wave42 mapping reason counts changed")

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
        "acceptance_contract",
    }
    if missing := required_evidence.difference(evidence.columns):
        raise ValueError(f"Wave42 evidence is missing columns: {sorted(missing)}")
    if evidence[list(required_evidence)].apply(
        lambda column: column.astype(str).str.strip().eq("").any()
    ).any():
        raise ValueError("Wave42 evidence has incomplete provenance")
    if len(evidence) != 1_390 or evidence["source_record_id"].duplicated().any():
        raise ValueError("Wave42 evidence row or record-ID count changed")
    accepted_species = set(
        mapping.loc[
            mapping["mapping_reason"].eq("accepted_strict_two_backbone_external"),
            "accepted_species",
        ]
    )
    evidence_checks = (
        evidence["accepted_species"].isin(accepted_species),
        ~evidence["accepted_species"].isin(target_species),
        evidence["axis"].eq("reproductive_assurance"),
        evidence["quality"].isin({"high", "medium"}),
        evidence["evidence_scope"].eq("external_congener_species_direct"),
        evidence["name_match_method"].eq("strict_wfo_gbif_two_backbone"),
        evidence["source_lineage"].eq(SHARED_LINEAGE),
        evidence["acceptance_contract"].eq(CONTRACT),
    )
    if not all(bool(check.all()) for check in evidence_checks):
        raise ValueError("Wave42 evidence violates its fail-closed contract")
    if _counts(evidence["source_provider"]) != EXPECTED_PROVIDER_COUNTS:
        raise ValueError("Wave42 provider evidence counts changed")
    if _counts(evidence["trait_name"]) != EXPECTED_TRAIT_COUNTS:
        raise ValueError("Wave42 trait evidence counts changed")

    selected = selection.loc[selection["selection_reason"].eq("selected")]
    if len(selected) != len(evidence):
        raise ValueError("Wave42 selected rows do not match evidence rows")
    selected_by_id = selected.set_index("source_record_id")
    evidence_by_id = evidence.set_index("source_record_id")
    if set(selected_by_id.index) != set(evidence_by_id.index):
        raise ValueError("Wave42 selected record IDs do not match evidence")
    if not selected_by_id["accepted_species"].eq(
        evidence_by_id["accepted_species"]
    ).all():
        raise ValueError("Wave42 selection and evidence species disagree")
    source_mapping = mapping.set_index("source_name")
    selected_routes = selected_by_id.join(
        source_mapping[["accepted_species", "wfo_family"]],
        on="source_name",
        rsuffix="_mapping",
        validate="many_to_one",
    )
    if not selected_routes["accepted_species"].eq(
        selected_routes["accepted_species_mapping"]
    ).all() or not selected_routes["source_family"].eq(
        selected_routes["wfo_family"]
    ).all():
        raise ValueError("Wave42 selected rows violate accepted-name or family gates")

    summary = {
        "contract": "wave42_external_bulk_checkpoint_validation_v1",
        "fixed_target_species": expected_species,
        "queries": {
            "local_wfo_2026_06": EXPECTED_NAMES,
            "reused_gbif": 1_453,
            "new_gbif": 64,
            "total_network_requests": 64,
        },
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
            "by_provider": _counts(evidence["source_provider"]),
            "by_trait": _counts(evidence["trait_name"]),
            "entered_fixed_target_direct_coverage": 0,
            "shared_redistribution_guard_lineages": int(
                evidence["source_lineage"].nunique()
            ),
        },
        "checks": {
            "fixed_denominator": True,
            "exact_wfo_2026_06_and_gbif": True,
            "external_species_only": True,
            "trait_specific_rows": True,
            "complete_provenance": True,
            "provider_redistribution_guarded": True,
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
