"""Validate the frozen Wave44 FloraWeb/GloPL reproductive packet."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd

EXPECTED_SPECIES = 106_295
EXPECTED_NAMES = 836
EXPECTED_MAPPING_REASONS = {
    "accepted_strict_two_backbone_direct": 85,
    "accepted_strict_two_backbone_external": 518,
    "backbone_family_conflict": 9,
    "backbones_disagree": 73,
    "gbif_not_exact": 24,
    "source_wfo_family_conflict": 7,
    "wfo_ambiguous_exact_name": 57,
    "wfo_no_accepted_species_family_route": 48,
    "wfo_no_exact_species_name": 15,
}
EXPECTED_SHA256 = {
    "direct": "bce976497834bcf5c6eac73242cb21a76730935affe1b14e28a979c5b11eabd5",
    "external": "5d4744670032a61c7e1fa68ffdeb3c8f4ca0cef4d0c8f7452938d9e15421e407",
    "mapping": "85678b32511283e2aa0d0f64e198d9af55423e4cfd2df3d7badd7576bcb251c0",
    "selection": "bc8ff21bf30aab6583861c1a199a76e31baf6c155c7bd9205df2f59fe3da59b0",
}
DIRECT_CONTRACT = "species_direct_strict_wfo_gbif_two_backbone_v1"
EXTERNAL_CONTRACT = "external_congener_species_direct_strict_two_backbone_v1"
TRAITS = {
    "autonomous_selfing_capacity",
    "cleistogamy",
    "self_incompatibility",
}


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


def _validate_identity_gate(
    mapping: pd.DataFrame,
    *,
    target_species: set[str],
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
        raise ValueError(f"Wave44 mapping is missing columns: {sorted(missing)}")
    if len(mapping) != EXPECTED_NAMES or mapping["source_name"].duplicated().any():
        raise ValueError("Wave44 mapping must contain one row per source name")

    accepted = mapping.loc[
        mapping["mapping_reason"].str.startswith("accepted_strict_two_backbone")
    ]
    checks = (
        accepted["accepted_species"].eq(accepted["wfo_accepted_species"]),
        accepted["accepted_species"].eq(accepted["gbif_accepted_species"]),
        accepted["source_family"].eq("")
        | accepted["source_family"].eq(accepted["wfo_family"]),
        accepted["wfo_family"].eq(accepted["gbif_family"]),
        accepted["wfo_release"].eq("2026-06"),
        accepted["gbif_match_type"].eq("EXACT"),
        accepted["gbif_rank"].eq("SPECIES"),
        accepted["gbif_kingdom"].eq("Plantae"),
    )
    if not all(bool(check.all()) for check in checks):
        raise ValueError("Wave44 accepted mapping violates the identity gate")
    direct = accepted["mapping_reason"].eq(
        "accepted_strict_two_backbone_direct"
    )
    if not accepted.loc[direct, "accepted_species"].isin(target_species).all():
        raise ValueError("Wave44 direct mapping is outside the fixed target")
    if accepted.loc[~direct, "accepted_species"].isin(target_species).any():
        raise ValueError("Wave44 external mapping entered the fixed target")
    rejected = mapping.loc[
        ~mapping["mapping_reason"].str.startswith("accepted_strict_two_backbone")
    ]
    if rejected["accepted_species"].ne("").any():
        raise ValueError("Wave44 rejected mappings retained an accepted species")


def _validate_evidence(
    frame: pd.DataFrame,
    *,
    target_species: set[str],
    external: bool,
) -> None:
    required = {
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
    if missing := required.difference(frame.columns):
        raise ValueError(f"Wave44 evidence is missing columns: {sorted(missing)}")
    if frame[list(required)].apply(
        lambda column: column.astype(str).str.strip().eq("").any()
    ).any():
        raise ValueError("Wave44 evidence has incomplete provenance")
    if frame["source_record_id"].duplicated().any():
        raise ValueError("Wave44 evidence has duplicate record IDs")
    checks = (
        frame["axis"].eq("reproductive_assurance"),
        frame["trait_name"].isin(TRAITS),
        frame["quality"].eq("high"),
        frame["name_match_method"].eq("strict_wfo_gbif_two_backbone"),
        frame["source_lineage"].str.startswith(
            ("biolflor:floraweb:", "glopl-burns-trait-compilation:")
        ),
        frame["lineage_method"].isin(
            {
                "underlying_glopl_burns_trait_compilation",
                "underlying_source_accepted_species_trait",
            }
        ),
    )
    if external:
        checks += (
            ~frame["accepted_species"].isin(target_species),
            frame["evidence_scope"].eq("external_congener_species_direct"),
            frame["acceptance_contract"].eq(EXTERNAL_CONTRACT),
        )
    else:
        checks += (
            frame["accepted_species"].isin(target_species),
            frame["evidence_scope"].eq("synonym_direct"),
            frame["acceptance_contract"].eq(DIRECT_CONTRACT),
        )
    if not all(bool(check.all()) for check in checks):
        label = "external" if external else "direct"
        raise ValueError(f"Wave44 {label} evidence violates its contract")


def validate_packet(
    *,
    packet_dir: Path,
    target_coverage_csv: Path,
    output_json: Path,
    expected_species: int = EXPECTED_SPECIES,
) -> dict[str, Any]:
    paths = {
        "source_manifest": packet_dir / "source_manifest.json",
        "candidate_summary": packet_dir
        / "wave44_external_reproductive_candidate_summary.json",
        "manifest": packet_dir / "wave44_external_reproductive_name_manifest.csv",
        "floraweb_rows": packet_dir / "wave44_floraweb_external_source_rows.csv.gz",
        "glopl_rows": packet_dir / "wave44_glopl_external_source_rows.csv.gz",
        "wfo_mapping": packet_dir / "wave44_wfo_2026_06_local_mapping.csv.gz",
        "wfo_summary": packet_dir
        / "wave44_wfo_2026_06_local_mapping_summary.json",
        "gbif_snapshot": packet_dir / "wave44_gbif_response_snapshot.jsonl.gz",
        "gbif_summary": packet_dir / "wave44_gbif_query_summary.json",
        "mapping": packet_dir / "wave44_two_backbone_mapping_audit.csv.gz",
        "direct": packet_dir
        / "wave44_synonym_direct_reproductive_evidence.csv.gz",
        "external": packet_dir
        / "wave44_external_congener_reproductive_evidence.csv.gz",
        "selection": packet_dir
        / "wave44_external_reproductive_selection_audit.csv.gz",
        "evidence_summary": packet_dir
        / "wave44_external_reproductive_evidence_summary.json",
    }
    if missing := [str(path) for path in paths.values() if not path.is_file()]:
        raise ValueError(f"Wave44 packet is incomplete: {missing}")

    target = pd.read_csv(
        target_coverage_csv,
        usecols=["accepted_species", "axis"],
        dtype=str,
    ).fillna("")
    target_species = set(target["accepted_species"])
    if len(target) != expected_species * 3 or len(target_species) != expected_species:
        raise ValueError("Wave44 target denominator mismatch")

    manifest = pd.read_csv(paths["manifest"], dtype=str).fillna("")
    local_wfo = pd.read_csv(paths["wfo_mapping"], dtype=str).fillna("")
    mapping = pd.read_csv(paths["mapping"], dtype=str).fillna("")
    direct = pd.read_csv(paths["direct"], dtype=str).fillna("")
    external = pd.read_csv(paths["external"], dtype=str).fillna("")
    selection = pd.read_csv(paths["selection"], dtype=str).fillna("")
    gbif = _read_jsonl_gz(paths["gbif_snapshot"])

    names = set(manifest["source_name"])
    if len(manifest) != EXPECTED_NAMES or len(names) != EXPECTED_NAMES:
        raise ValueError("Wave44 name manifest count changed")
    if set(local_wfo["source_name"]) != names or set(mapping["source_name"]) != names:
        raise ValueError("Wave44 name tables cover different names")
    gbif_names = [str(record.get("source_name", "")) for record in gbif]
    if len(gbif_names) != EXPECTED_NAMES or set(gbif_names) != names:
        raise ValueError("Wave44 GBIF snapshot does not cover the manifest")

    _validate_identity_gate(mapping, target_species=target_species)
    if _counts(mapping["mapping_reason"]) != EXPECTED_MAPPING_REASONS:
        raise ValueError("Wave44 mapping reason counts changed")
    _validate_evidence(direct, target_species=target_species, external=False)
    _validate_evidence(external, target_species=target_species, external=True)
    if len(direct) != 143 or direct["accepted_species"].nunique() != 80:
        raise ValueError("Wave44 direct evidence counts changed")
    if len(external) != 950 or external["accepted_species"].nunique() != 518:
        raise ValueError("Wave44 external evidence counts changed")
    if external[["accepted_species", "trait_name"]].drop_duplicates().shape[0] != 948:
        raise ValueError("Wave44 external species x trait count changed")

    selected = selection.loc[selection["selection_reason"].eq("selected")]
    if len(selected) != len(direct) + len(external):
        raise ValueError("Wave44 selected rows do not match evidence rows")
    declared = json.loads(paths["evidence_summary"].read_text(encoding="utf-8"))
    source_manifest = json.loads(
        paths["source_manifest"].read_text(encoding="utf-8")
    )
    if (
        source_manifest["fixed_target_species"] != expected_species
        or len(source_manifest["sources"]) != 2
        or source_manifest["inference_policy"]["family_inference"]
        or source_manifest["inference_policy"]["global_fallback"]
        or source_manifest["inference_policy"][
            "external_congener_rows_enter_direct_coverage"
        ]
    ):
        raise ValueError("Wave44 source manifest contract changed")
    observed_hashes = {
        label: _sha256(paths[label])
        for label in ("direct", "external", "mapping", "selection")
    }
    if declared["sha256"] != EXPECTED_SHA256 or observed_hashes != EXPECTED_SHA256:
        raise ValueError("Wave44 frozen packet hashes changed")
    queries = json.loads(paths["gbif_summary"].read_text(encoding="utf-8"))
    if queries != declared["queries"] or queries["query_cost_usd"] != 0:
        raise ValueError("Wave44 query accounting changed")

    summary = {
        "contract": "wave44_external_reproduction_checkpoint_v1",
        "fixed_target_species": expected_species,
        "mapping_reason_counts": _counts(mapping["mapping_reason"]),
        "evidence": {
            "direct_rows": len(direct),
            "direct_species": int(direct["accepted_species"].nunique()),
            "direct_species_trait": int(
                direct[["accepted_species", "trait_name"]]
                .drop_duplicates()
                .shape[0]
            ),
            "external_rows": len(external),
            "external_species": int(external["accepted_species"].nunique()),
            "external_species_trait": int(
                external[["accepted_species", "trait_name"]]
                .drop_duplicates()
                .shape[0]
            ),
            "external_source_lineages": int(external["source_lineage"].nunique()),
        },
        "queries": queries,
        "query_cost_usd": queries["query_cost_usd"],
        "sha256": observed_hashes,
        "checks": {
            "strict_wfo_gbif_two_backbone": True,
            "direct_inside_fixed_target": True,
            "external_outside_fixed_target": True,
            "external_not_confirmatory_direct": True,
            "source_lineage_complete": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
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
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--expected-species", type=int, default=EXPECTED_SPECIES)
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
