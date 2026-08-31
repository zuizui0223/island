"""Validate Wave46 sources and materialize the corrected public-Web packet."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.wave37_europe_pmc_checkpoint import _sha256, _write_gzip_csv
from island_v2.wave45_promoted_reproductive_checkpoint import _validate_evidence

EXPECTED_SPECIES = 106_295
PUBLIC_WEB_RUN_ID = 32_710_232_989
PUBLIC_WEB_ROWS = 118_884


def validate_packet(
    *,
    packet_dir: Path,
    latest_public_web_csv: Path,
    target_coverage_csv: Path,
    output_dir: Path,
    output_json: Path,
    expected_species: int = EXPECTED_SPECIES,
) -> dict[str, Any]:
    paths = {
        "manifest": packet_dir / "source_manifest.json",
        "new_direct": packet_dir / "wave46_reviewed_direct_evidence.csv",
        "review": packet_dir / "wave46_source_review_audit.csv",
        "identity": packet_dir / "wave46_identity_audit.csv",
        "latest_public_web": latest_public_web_csv,
        "target_coverage": target_coverage_csv,
    }
    if missing := [str(path) for path in paths.values() if not path.is_file()]:
        raise ValueError(f"Wave46 packet is incomplete: {missing}")

    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    target = pd.read_csv(
        target_coverage_csv,
        usecols=["accepted_species", "axis"],
        dtype=str,
    ).fillna("")
    target_species = set(target["accepted_species"])
    if len(target) != expected_species * 3 or len(target_species) != expected_species:
        raise ValueError("Wave46 target denominator mismatch")

    public_web = pd.read_csv(latest_public_web_csv, dtype=str).fillna("")
    required_public_web = {
        "accepted_species",
        "axis",
        "trait_name",
        "source_lineage",
        "source_run_id",
        "source_artifact",
    }
    if missing := required_public_web.difference(public_web.columns):
        raise ValueError(f"Wave46 public-Web evidence is missing columns: {sorted(missing)}")
    declared_public_web = manifest["correct_latest_public_web"]
    if (
        declared_public_web["run_id"] != PUBLIC_WEB_RUN_ID
        or declared_public_web["evidence_rows"] != PUBLIC_WEB_ROWS
        or len(public_web) != PUBLIC_WEB_ROWS
    ):
        raise ValueError("Wave46 latest public-Web pin or row count changed")
    if (
        public_web[list(required_public_web)]
        .apply(lambda column: column.str.strip().eq("").any())
        .any()
    ):
        raise ValueError("Wave46 latest public-Web evidence has blank routing fields")

    direct = pd.read_csv(paths["new_direct"], dtype=str).fillna("")
    review = pd.read_csv(paths["review"], dtype=str).fillna("")
    identity = pd.read_csv(paths["identity"], dtype=str).fillna("")
    _validate_evidence(direct, label="new direct")
    expected_direct_rows = sum(
        int(source["reviewed_rows"]) for source in manifest["new_primary_sources"]
    )
    if (
        len(direct) != expected_direct_rows
        or len(review) != expected_direct_rows
        or len(identity) != direct["accepted_species"].nunique()
    ):
        raise ValueError("Wave46 reviewed evidence counts changed")
    if not direct["accepted_species"].isin(target_species).all():
        raise ValueError("Wave46 direct evidence is outside the fixed target")
    if not direct["evidence_scope"].eq("species_direct").all():
        raise ValueError("Wave46 direct evidence has the wrong scope")
    if set(direct["source_record_id"]) != set(review["record_id"]):
        raise ValueError("Wave46 review does not cover every direct row")
    if not review["accepted_correct"].str.casefold().eq("true").all():
        raise ValueError("Wave46 review contains a rejected row")
    if (
        not direct.merge(
            review[["record_id", "exact_supporting_quote"]],
            left_on="source_record_id",
            right_on="record_id",
            validate="one_to_one",
        )
        .apply(lambda row: row["source_excerpt"] == row["exact_supporting_quote"], axis=1)
        .all()
    ):
        raise ValueError("Wave46 source quote and review quote differ")
    if (
        set(identity["accepted_species"]) != set(direct["accepted_species"])
        or not identity["target_universe_status"].eq("in_fixed_106295_species").all()
        or not identity["name_match_method"].eq("accepted_name_exact").all()
    ):
        raise ValueError("Wave46 direct identity audit failed")

    source = manifest["new_primary_sources"][0]
    expected_source_file = f"html_sha256:{source['retrieved_content_sha256']}"
    if (
        manifest["fixed_target_species"] != expected_species
        or manifest["inference_policy"]["join_key"] != "genus x axis x trait_name"
        or manifest["inference_policy"]["family_inference"]
        or manifest["inference_policy"]["global_fallback"]
        or manifest["inference_policy"]["reproductive_traits_interchangeable"]
        or not direct["source_file"].eq(expected_source_file).all()
    ):
        raise ValueError("Wave46 source manifest contract changed")

    combined = pd.concat([public_web, direct], ignore_index=True, sort=False).fillna("")
    output_dir.mkdir(parents=True, exist_ok=True)
    combined_path = output_dir / "wave46_corrected_public_web_direct_packet.csv.gz"
    _write_gzip_csv(combined, combined_path)
    summary: dict[str, Any] = {
        "contract": "wave46_corrected_public_web_checkpoint_v1",
        "fixed_target_species": expected_species,
        "evidence": {
            "latest_public_web_rows": len(public_web),
            "latest_public_web_species": int(public_web["accepted_species"].nunique()),
            "latest_public_web_species_trait": int(
                public_web[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
            ),
            "new_direct_rows": len(direct),
            "new_direct_species": int(direct["accepted_species"].nunique()),
            "new_direct_species_trait": int(
                direct[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
            ),
            "combined_direct_rows": len(combined),
            "review_rows": len(review),
            "identity_rows": len(identity),
            "new_external_species_trait": 0,
        },
        "queries": {
            "formal_search_api_queries": 0,
            "direct_source_pages_retrieved": 1,
            "query_cost_usd": 0,
        },
        "query_cost_usd": 0,
        "checks": {
            "fixed_denominator": True,
            "formal_wave33_baseline_pinned": (
                manifest["formal_wave33_baseline"]["run_id"] == 32_932_103_226
            ),
            "latest_public_web_artifact_corrected": True,
            "all_new_rows_reviewed": True,
            "exact_quote_and_provenance_complete": True,
            "trait_specific_genus_join": True,
            "reproductive_traits_not_interchanged": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
        },
        "input_sha256": {label: _sha256(path) for label, path in paths.items()},
        "artifact_sha256": {combined_path.name: _sha256(combined_path)},
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
    parser.add_argument("--latest-public-web-csv", required=True, type=Path)
    parser.add_argument("--target-coverage-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--expected-species", type=int, default=EXPECTED_SPECIES)
    args = parser.parse_args()
    summary = validate_packet(
        packet_dir=args.packet_dir,
        latest_public_web_csv=args.latest_public_web_csv,
        target_coverage_csv=args.target_coverage_csv,
        output_dir=args.output_dir,
        output_json=args.output_json,
        expected_species=args.expected_species,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
