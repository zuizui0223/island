"""Validate and materialize the frozen Wave45 reviewed evidence packet."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.wave37_europe_pmc_checkpoint import _sha256, _write_gzip_csv

EXPECTED_SPECIES = 106_295
PROMOTED_ROWS = 208
NEW_DIRECT_ROWS = 8
NEW_EXTERNAL_ROWS = 3
REVIEW_ROWS = 11
IDENTITY_ROWS = 11
REPRODUCTIVE_TRAITS = {
    "autonomous_selfing_capacity",
    "cleistogamy",
    "mating_system",
    "self_incompatibility",
}
REQUIRED_EVIDENCE_COLUMNS = {
    "accepted_species",
    "axis",
    "trait_name",
    "normalized_value",
    "quality",
    "source_group",
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


def _validate_evidence(frame: pd.DataFrame, *, label: str) -> None:
    if missing := REQUIRED_EVIDENCE_COLUMNS.difference(frame.columns):
        raise ValueError(f"Wave45 {label} evidence is missing columns: {sorted(missing)}")
    required = sorted(REQUIRED_EVIDENCE_COLUMNS)
    if frame[required].apply(lambda col: col.str.strip().eq("").any()).any():
        raise ValueError(f"Wave45 {label} evidence has incomplete provenance")
    if frame["source_record_id"].duplicated().any():
        raise ValueError(f"Wave45 {label} evidence has duplicate record IDs")
    if not frame["quality"].isin({"high", "medium"}).all():
        raise ValueError(f"Wave45 {label} evidence has an invalid quality")
    reproductive = frame["axis"].eq("reproductive_assurance")
    if not frame.loc[reproductive, "trait_name"].isin(REPRODUCTIVE_TRAITS).all():
        raise ValueError("Wave45 interchanges distinct reproductive traits")
    if frame["trait_name"].isin({"pollen_vector_mode", "reward_type"}).any():
        raise ValueError("Wave45 strict axes contain an independent auxiliary trait")


def validate_packet(
    *,
    packet_dir: Path,
    promoted_direct_csv: Path,
    target_coverage_csv: Path,
    output_dir: Path,
    output_json: Path,
    expected_species: int = EXPECTED_SPECIES,
) -> dict[str, Any]:
    paths = {
        "source_manifest": packet_dir / "source_manifest.json",
        "new_direct": packet_dir / "wave45_reviewed_direct_evidence.csv",
        "new_external": packet_dir / "wave45_external_congener_evidence.csv",
        "review": packet_dir / "wave45_source_review_audit.csv",
        "identity": packet_dir / "wave45_identity_reuse_audit.csv",
        "promoted_direct": promoted_direct_csv,
        "target_coverage": target_coverage_csv,
    }
    if missing := [str(path) for path in paths.values() if not path.is_file()]:
        raise ValueError(f"Wave45 packet is incomplete: {missing}")

    target = pd.read_csv(
        target_coverage_csv,
        usecols=["accepted_species", "axis"],
        dtype=str,
    ).fillna("")
    target_species = set(target["accepted_species"])
    if len(target) != expected_species * 3 or len(target_species) != expected_species:
        raise ValueError("Wave45 target denominator mismatch")

    promoted = pd.read_csv(promoted_direct_csv, dtype=str).fillna("")
    direct = pd.read_csv(paths["new_direct"], dtype=str).fillna("")
    external = pd.read_csv(paths["new_external"], dtype=str).fillna("")
    review = pd.read_csv(paths["review"], dtype=str).fillna("")
    identity = pd.read_csv(paths["identity"], dtype=str).fillna("")
    manifest = json.loads(paths["source_manifest"].read_text(encoding="utf-8"))

    for label, frame in (
        ("promoted", promoted),
        ("direct", direct),
        ("external", external),
    ):
        _validate_evidence(frame, label=label)
    if len(promoted) != PROMOTED_ROWS or len(direct) != NEW_DIRECT_ROWS:
        raise ValueError("Wave45 reviewed direct evidence counts changed")
    if len(external) != NEW_EXTERNAL_ROWS:
        raise ValueError("Wave45 external evidence count changed")
    if not promoted["accepted_species"].isin(target_species).all():
        raise ValueError("Wave45 promoted direct evidence is outside the fixed target")
    if not direct["accepted_species"].isin(target_species).all():
        raise ValueError("Wave45 new direct evidence is outside the fixed target")
    if external["accepted_species"].isin(target_species).any():
        raise ValueError("Wave45 external congener evidence entered the fixed target")
    if not direct["evidence_scope"].eq("species_direct").all():
        raise ValueError("Wave45 new direct evidence has the wrong scope")
    if not external["evidence_scope"].eq("external_congener_species_direct").all():
        raise ValueError("Wave45 external evidence has the wrong scope")
    if not external["name_match_method"].eq("strict_wfo_gbif_two_backbone").all():
        raise ValueError("Wave45 external evidence failed the two-backbone gate")

    if len(review) != REVIEW_ROWS or len(identity) != IDENTITY_ROWS:
        raise ValueError("Wave45 audit row counts changed")
    new_record_ids = set(direct["source_record_id"]) | set(external["source_record_id"])
    if set(review["record_id"]) != new_record_ids:
        raise ValueError("Wave45 source review does not cover every new evidence row")
    if not review["accepted_correct"].str.casefold().eq("true").all():
        raise ValueError("Wave45 review contains an unaccepted evidence row")
    expected_identity = set(direct["accepted_species"]) | set(external["accepted_species"])
    if set(identity["accepted_species"]) != expected_identity:
        raise ValueError("Wave45 identity audit does not cover every new species")
    external_identity = identity["target_universe_status"].eq("external_congener_only")
    strict_identity = (
        identity.loc[external_identity, "name_match_method"]
        .eq("strict_wfo_gbif_two_backbone")
        .all()
        and identity.loc[external_identity, "wfo_match_id"].ne("").all()
        and identity.loc[external_identity, "gbif_usage_key"].ne("").all()
        and identity.loc[external_identity, "family"]
        .eq(identity.loc[external_identity, "gbif_family"])
        .all()
    )
    if not strict_identity:
        raise ValueError("Wave45 external identity audit violates the strict gate")

    promoted_declared = manifest["promoted_reviewed_packet"]
    if (
        manifest["fixed_target_species"] != expected_species
        or promoted_declared["rows"] != PROMOTED_ROWS
        or promoted_declared["sha256"] != _sha256(promoted_direct_csv)
        or manifest["inference_policy"]["join_key"] != "genus x axis x trait_name"
        or manifest["inference_policy"]["family_inference"]
        or manifest["inference_policy"]["global_fallback"]
        or manifest["inference_policy"]["reproductive_traits_interchangeable"]
    ):
        raise ValueError("Wave45 source manifest contract changed")

    combined = pd.concat([promoted, direct], ignore_index=True)
    if combined["source_record_id"].duplicated().any():
        raise ValueError("Wave45 combined direct packet has duplicate record IDs")
    output_dir.mkdir(parents=True, exist_ok=True)
    combined_path = output_dir / "wave45_combined_reviewed_direct_evidence.csv.gz"
    external_path = output_dir / "wave45_external_congener_evidence.csv.gz"
    _write_gzip_csv(combined, combined_path)
    _write_gzip_csv(external, external_path)

    summary: dict[str, Any] = {
        "contract": "wave45_promoted_reproductive_checkpoint_v1",
        "fixed_target_species": expected_species,
        "evidence": {
            "promoted_direct_rows": len(promoted),
            "new_direct_rows": len(direct),
            "combined_direct_rows": len(combined),
            "combined_direct_species": int(combined["accepted_species"].nunique()),
            "combined_direct_species_trait": int(
                combined[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
            ),
            "new_external_rows": len(external),
            "new_external_species": int(external["accepted_species"].nunique()),
            "new_external_species_trait": int(
                external[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
            ),
            "review_rows": len(review),
            "identity_rows": len(identity),
        },
        "queries": {
            "formal_search_api_queries": 0,
            "query_cost_usd": 0,
            "source_documents": len(manifest["new_primary_sources"]),
        },
        "query_cost_usd": 0,
        "checks": {
            "fixed_denominator": True,
            "all_new_rows_reviewed": True,
            "direct_inside_fixed_target": True,
            "external_outside_fixed_target": True,
            "external_strict_two_backbone": True,
            "exact_quote_and_provenance_complete": True,
            "trait_specific_genus_join": True,
            "reproductive_traits_not_interchanged": True,
            "auxiliary_traits_outside_strict_axes": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
        },
        "input_sha256": {label: _sha256(path) for label, path in paths.items()},
        "artifact_sha256": {
            combined_path.name: _sha256(combined_path),
            external_path.name: _sha256(external_path),
        },
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
    parser.add_argument("--promoted-direct-csv", required=True, type=Path)
    parser.add_argument("--target-coverage-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--expected-species", type=int, default=EXPECTED_SPECIES)
    args = parser.parse_args()
    summary = validate_packet(
        packet_dir=args.packet_dir,
        promoted_direct_csv=args.promoted_direct_csv,
        target_coverage_csv=args.target_coverage_csv,
        output_dir=args.output_dir,
        output_json=args.output_json,
        expected_species=args.expected_species,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
