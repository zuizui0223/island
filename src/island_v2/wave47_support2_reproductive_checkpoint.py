"""Validate the frozen Wave47 support-2 reproductive evidence packet."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.all_evidence_trait_audit import EVIDENCE_COLUMNS
from island_v2.wave37_europe_pmc_checkpoint import _sha256, _write_gzip_csv
from island_v2.wave45_promoted_reproductive_checkpoint import _validate_evidence

EXPECTED_SPECIES = 106_295
FORMAL_WAVE33_RUN_ID = 32_932_103_226
BASELINE_WAVE46_RUN_ID = 33_376_311_877
SOURCE_ROWS = 2


def validate_packet(
    *,
    packet_dir: Path,
    target_coverage_csv: Path,
    retrieved_source_file: Path,
    output_dir: Path,
    output_json: Path,
    expected_species: int = EXPECTED_SPECIES,
) -> dict[str, Any]:
    paths = {
        "manifest": packet_dir / "source_manifest.json",
        "external": packet_dir / "wave47_external_congener_evidence.csv",
        "review": packet_dir / "wave47_source_review_audit.csv",
        "identity": packet_dir / "wave47_identity_audit.csv",
        "target_coverage": target_coverage_csv,
        "retrieved_source": retrieved_source_file,
    }
    if missing := [str(path) for path in paths.values() if not path.is_file()]:
        raise ValueError(f"Wave47 packet is incomplete: {missing}")

    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    target = pd.read_csv(
        target_coverage_csv,
        usecols=["accepted_species", "axis"],
        dtype=str,
    ).fillna("")
    target_species = set(target["accepted_species"])
    if len(target) != expected_species * 3 or len(target_species) != expected_species:
        raise ValueError("Wave47 target denominator mismatch")

    external = pd.read_csv(paths["external"], dtype=str).fillna("")
    review = pd.read_csv(paths["review"], dtype=str).fillna("")
    identity = pd.read_csv(paths["identity"], dtype=str).fillna("")
    _validate_evidence(external, label="external")
    if len(external) != SOURCE_ROWS or len(review) != SOURCE_ROWS or len(identity) != 1:
        raise ValueError("Wave47 reviewed evidence counts changed")
    if external["accepted_species"].isin(target_species).any():
        raise ValueError("Wave47 external congener entered the fixed target")
    if not external["evidence_scope"].eq("external_congener_species_direct").all():
        raise ValueError("Wave47 external evidence has the wrong scope")
    if not external["name_match_method"].eq("strict_wfo_gbif_two_backbone").all():
        raise ValueError("Wave47 external evidence failed the two-backbone gate")
    if set(external["source_record_id"]) != set(review["record_id"]):
        raise ValueError("Wave47 review does not cover every evidence row")
    if not review["accepted_correct"].str.casefold().eq("true").all():
        raise ValueError("Wave47 review contains a rejected row")
    merged = external.merge(
        review[["record_id", "exact_supporting_quote"]],
        left_on="source_record_id",
        right_on="record_id",
        validate="one_to_one",
    )
    if not merged["source_excerpt"].eq(merged["exact_supporting_quote"]).all():
        raise ValueError("Wave47 source quote and review quote differ")

    identity_row = identity.iloc[0]
    identity_ok = (
        identity_row["accepted_species"] == "Buxus wallichiana"
        and identity_row["target_universe_status"] == "external_congener_only"
        and identity_row["name_match_method"] == "strict_wfo_gbif_two_backbone"
        and identity_row["wfo_match_id"] == "wfo-0000576661"
        and identity_row["wfo_status"].casefold() == "accepted"
        and identity_row["wfo_rank"].casefold() == "species"
        and identity_row["gbif_usage_key"] == "3988068"
        and identity_row["gbif_status"] == "ACCEPTED"
        and identity_row["gbif_rank"] == "SPECIES"
        and identity_row["gbif_match_type"] == "EXACT"
        and int(identity_row["gbif_confidence"]) >= 95
        and identity_row["family"] == identity_row["gbif_family"] == "Buxaceae"
    )
    if not identity_ok:
        raise ValueError("Wave47 strict WFO/GBIF identity gate failed")

    source = manifest["source"]
    source_sha = _sha256(retrieved_source_file)
    expected_source_file = f"pdf_sha256:{source_sha}"
    inference = manifest["inference_policy"]
    if (
        manifest["fixed_target_species"] != expected_species
        or manifest["formal_wave33_baseline"]["run_id"] != FORMAL_WAVE33_RUN_ID
        or manifest["immediate_formal_baseline"]["run_id"] != BASELINE_WAVE46_RUN_ID
        or source["reviewed_rows"] != SOURCE_ROWS
        or source["retrieved_content_sha256"] != source_sha
        or not external["source_file"].eq(expected_source_file).all()
        or inference["join_key"] != "genus x axis x trait_name"
        or inference["minimum_species"] != 3
        or inference["family_inference"]
        or inference["global_fallback"]
        or inference["reproductive_traits_interchangeable"]
    ):
        raise ValueError("Wave47 source or inference contract changed")

    output_dir.mkdir(parents=True, exist_ok=True)
    external_path = output_dir / "wave47_external_congener_evidence.csv.gz"
    empty_direct_path = output_dir / "wave47_direct_evidence_empty.csv.gz"
    _write_gzip_csv(external, external_path)
    _write_gzip_csv(pd.DataFrame(columns=EVIDENCE_COLUMNS), empty_direct_path)

    summary: dict[str, Any] = {
        "contract": "wave47_support2_reproductive_checkpoint_v1",
        "fixed_target_species": expected_species,
        "formal_wave33_run_id": FORMAL_WAVE33_RUN_ID,
        "baseline_formal_run_id": BASELINE_WAVE46_RUN_ID,
        "evidence": {
            "new_external_rows": len(external),
            "new_external_species": int(external["accepted_species"].nunique()),
            "new_external_species_trait": int(
                external[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
            ),
            "review_rows": len(review),
            "identity_rows": len(identity),
            "new_direct_rows": 0,
        },
        "queries": {
            "formal_search_api_queries": 0,
            "source_pages_retrieved": 1,
            "query_cost_usd": 0,
        },
        "query_cost_usd": 0,
        "checks": {
            "fixed_denominator": True,
            "formal_wave33_baseline_pinned": True,
            "immediate_wave46_baseline_pinned": True,
            "all_new_rows_reviewed": True,
            "retrieved_source_digest_verified": True,
            "exact_quote_and_provenance_complete": True,
            "external_outside_fixed_target": True,
            "external_strict_two_backbone": True,
            "trait_specific_genus_join": True,
            "reproductive_traits_not_interchanged": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
        },
        "input_sha256": {label: _sha256(path) for label, path in paths.items()},
        "artifact_sha256": {
            external_path.name: _sha256(external_path),
            empty_direct_path.name: _sha256(empty_direct_path),
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
    parser.add_argument("--target-coverage-csv", required=True, type=Path)
    parser.add_argument("--retrieved-source-file", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--expected-species", type=int, default=EXPECTED_SPECIES)
    args = parser.parse_args()
    summary = validate_packet(
        packet_dir=args.packet_dir,
        target_coverage_csv=args.target_coverage_csv,
        retrieved_source_file=args.retrieved_source_file,
        output_dir=args.output_dir,
        output_json=args.output_json,
        expected_species=args.expected_species,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
