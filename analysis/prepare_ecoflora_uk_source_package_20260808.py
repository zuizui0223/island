"""Prepare novel Ecoflora species-direct traits from the archived UK_Flora snapshot.

The original Ecoflora site disallows automated crawling.  This adapter therefore
uses the immutable June 2020 snapshot in ``seancanobrien/UK_Flora`` and retains
the original Ecoflora record URL, field/value excerpt, snapshot commit, and input
hash.  It only emits species-trait pairs absent from the supplied current direct
lineage ledger; the full source inventory remains described in the summary.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

SOURCE_COMMIT = "8e524c994086955bf8f95c5a05213acad1cd1677"
SOURCE_REPOSITORY = "https://github.com/seancanobrien/UK_Flora"
SOURCE_RUN_ID = "ecoflora-uk-flora-snapshot-2020-20260808"
SOURCE_ARTIFACT = f"seancanobrien/UK_Flora@{SOURCE_COMMIT}"
REVIEWER = "codex_gpt5_ecoflora_stratified_record_audit_20260808"
REVIEWED_AT = "2026-08-08T14:15:00Z"

AXIS = {
    "flower_primary_color": "flower_colour",
    "inflorescence_display": "floral_structural_complexity",
    "self_incompatibility": "reproductive_assurance",
    "mating_system": "reproductive_assurance",
    "cleistogamy": "reproductive_assurance",
}
AUDIT_TARGET = {
    "flower_primary_color": 40,
    "inflorescence_display": 70,
    "self_incompatibility": 35,
    "mating_system": 35,
    "cleistogamy": 20,
}

COLOUR_PATTERNS = (
    (r"\bwhite\b|\bcream\b", "white"),
    (r"\bgreen\b|\bbrown\b", "green_brown_inconspicuous"),
    (r"\byellow\b|\borange\b|\bgold(?:en)?\b", "yellow_orange"),
    (r"\bred\b|\bpink\b|\bscarlet\b|\bcrimson\b", "red_pink"),
    (r"\bblue\b|\bpurple\b|\bviolet\b|\blilac\b|\bmauve\b", "blue_purple"),
    (r"\bblack\b", "other_described"),
)
INFLORESCENCE_PATTERNS = (
    (r"\bsolitary\b", "solitary"),
    (r"\bfew[- ]flower", "few_flowered"),
    (
        r"\braceme\b|\bspikes?\b|\bpanicle\b|\bthyrse\b|\bcatkin\b",
        "raceme_spike_panicle",
    ),
    (r"\bumbel\b|\bcorymb\b|\bcyme\b|\bwhorl\b", "umbel_corymb"),
    (r"\bhead\b|\bcapitulum\b", "composite_display"),
    (r"\bbrush\b|\bpuff\b", "brush_puff_display"),
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _states(raw: str, patterns: tuple[tuple[str, str], ...]) -> str:
    return "|".join(
        sorted({state for pattern, state in patterns if re.search(pattern, raw, re.IGNORECASE)})
    )


def _self_incompatibility(raw: str) -> str:
    absent = bool(re.search(r"(^|,\s*)none($|,)", raw, re.IGNORECASE))
    present = bool(re.search(r"present|gametophytic|sporophytic", raw, re.IGNORECASE))
    if absent and present:
        return "mixed_or_variable"
    if absent:
        return "SC"
    if present:
        return "SI"
    return ""


def _cleistogamy(raw: str) -> str:
    if re.search(r"entirely cleistogamous", raw, re.IGNORECASE):
        return "obligate"
    if re.search(
        r"some cleistogamous|sometimes cleistogamous|usually cleistogamous", raw, re.IGNORECASE
    ):
        return "facultative"
    # "not recorded", numeric corruption, and pseudo-cleistogamy are not absence.
    return ""


def _mating_system(raw: str) -> str:
    if re.search(r"apomictic|viviparous|no seed produced", raw, re.IGNORECASE):
        return ""
    mixed = bool(re.search(r"cross and self|cross or automatic self", raw, re.IGNORECASE))
    outcrossing = bool(re.search(r"obligatory cross|normally cross", raw, re.IGNORECASE))
    selfing = bool(re.search(r"normally self", raw, re.IGNORECASE))
    if mixed:
        return "mixed_mating"
    if outcrossing and selfing:
        return ""
    if outcrossing:
        return "predominantly_outcrossing"
    if selfing:
        return "predominantly_selfing"
    return ""


def _candidate_id(plant_no: str, trait_name: str) -> str:
    stable = f"ecoflora:{plant_no}:{trait_name}"
    return "ecoflora-" + hashlib.sha256(stable.encode()).hexdigest()[:20]


def _record(
    row: dict[str, str], trait_name: str, raw_field: str, normalized_value: str
) -> dict[str, str]:
    plant_no = row["id"].strip()
    species = row["accepted_species"]
    source_url = f"http://ecoflora.org.uk/search_ecochars.php?plant_no={plant_no}"
    return {
        "accepted_species": species,
        "axis": AXIS[trait_name],
        "trait_name": trait_name,
        "normalized_value": normalized_value,
        "quality": "medium",
        "source_group": "latest_public_web",
        "source_provider": "ecoflora_via_uk_flora_snapshot",
        "source_url": source_url,
        "source_record_id": _candidate_id(plant_no, trait_name),
        "source_citation": (
            "Ecological Flora of the British Isles, Ecoflora species record; "
            f"archived in {SOURCE_ARTIFACT}"
        ),
        "source_excerpt": f'Ecoflora field "{raw_field}": {row[raw_field].strip()}',
        "evidence_scope": "species_direct",
        "name_match_method": "ecoflora_scientific_name_exact_to_master",
        "source_lineage": f"ecoflora:plant_no:{plant_no}",
        "lineage_method": "original_ecoflora_plant_number",
        "source_run_id": SOURCE_RUN_ID,
        "source_artifact": SOURCE_ARTIFACT,
        "source_file": "data/old_data_scrape.csv",
        "acceptance_contract": "ecoflora_archived_field_value_medium_v1",
    }


def build(
    source_csv: Path,
    master_csv: Path,
    current_coverage_csv: Path,
    current_lineages_csv: Path,
    output_dir: Path,
) -> dict[str, object]:
    source = pd.read_csv(source_csv, sep="|", dtype=str, engine="python").fillna("")
    source["accepted_species"] = (
        source["species"]
        .str.replace("_", " ", regex=False)
        .str.strip()
        .str.replace(r"\s+", " ", regex=True)
        .str.capitalize()
    )
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    master_names = set(master["accepted_species"])
    strict_universe = set(
        pd.read_csv(
            current_coverage_csv,
            usecols=["accepted_species"],
            dtype=str,
        )["accepted_species"]
    )
    if len(strict_universe) != 106_295:
        raise ValueError(
            f"current strict coverage has {len(strict_universe)} species, expected 106295"
        )
    source = source.loc[
        source["accepted_species"].isin(master_names & strict_universe)
    ].copy()
    source = source.loc[
        ~source["accepted_species"].str.contains(
            r"\b(?:x|×)\b|\bhybrid(?:a|us|um)?\b|cultivar|cv\.",
            case=False,
            regex=True,
        )
        & ~source["accepted_species"].isin(
            {"Fragaria ananassa", "Rubus loganobaccus"}
        )
    ].copy()

    extractors = (
        ("flower_primary_color", "Flower colour", lambda value: _states(value, COLOUR_PATTERNS)),
        (
            "inflorescence_display",
            "Inflorescence type",
            lambda value: _states(value, INFLORESCENCE_PATTERNS),
        ),
        ("self_incompatibility", "Incompatibility systems", _self_incompatibility),
        ("mating_system", "Fertilization", _mating_system),
        ("cleistogamy", "Cleistogamy", _cleistogamy),
    )
    records: list[dict[str, str]] = []
    for row in source.to_dict("records"):
        for trait_name, raw_field, extractor in extractors:
            raw = row.get(raw_field, "").strip()
            value = extractor(raw) if raw else ""
            if value:
                records.append(_record(row, trait_name, raw_field, value))
    inventory = pd.DataFrame(records).drop_duplicates(
        ["accepted_species", "trait_name", "normalized_value"]
    )

    current = pd.read_csv(
        current_lineages_csv,
        usecols=["accepted_species", "trait_name"],
        dtype=str,
    ).drop_duplicates()
    current_pairs = set(zip(current["accepted_species"], current["trait_name"], strict=True))
    inventory["novel_species_trait"] = [
        (species, trait) not in current_pairs
        for species, trait in zip(
            inventory["accepted_species"], inventory["trait_name"], strict=True
        )
    ]
    evidence = inventory.loc[inventory["novel_species_trait"]].drop(
        columns="novel_species_trait"
    )
    evidence = evidence.sort_values(
        ["trait_name", "accepted_species", "source_record_id"]
    ).reset_index(drop=True)

    audit_parts = []
    for trait_name, target in AUDIT_TARGET.items():
        group = evidence.loc[evidence["trait_name"].eq(trait_name)]
        audit_parts.append(group.sample(n=min(target, len(group)), random_state=20260808))
    audit = pd.concat(audit_parts, ignore_index=True).sort_values(
        ["trait_name", "accepted_species"]
    )
    audit = audit.rename(
        columns={
            "source_provider": "provider",
            "source_record_id": "candidate_id",
            "source_excerpt": "exact_supporting_quote",
        }
    )
    audit["accepted_correct"] = ""
    audit["cultivar_status"] = "wild_or_unspecified_species_record"
    audit["reviewer"] = REVIEWER
    audit["reviewed_at_utc"] = REVIEWED_AT
    audit["audit_reason"] = ""

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "ecoflora_uk_snapshot_evidence.csv.gz"
    audit_path = output_dir / "ecoflora_uk_snapshot_manual_audit_200.csv"
    summary_path = output_dir / "ecoflora_uk_snapshot_source_package_summary.json"
    evidence.to_csv(evidence_path, index=False)
    audit.to_csv(audit_path, index=False)
    summary: dict[str, object] = {
        "contract": "ecoflora_archived_source_package_v1",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "source_repository": SOURCE_REPOSITORY,
        "source_commit": SOURCE_COMMIT,
        "source_csv_sha256": _sha256(source_csv),
        "source_snapshot_rows": len(source),
        "strict_universe_species": len(strict_universe),
        "full_inventory_evidence_rows": len(inventory),
        "full_inventory_species_trait": len(
            inventory.drop_duplicates(["accepted_species", "trait_name"])
        ),
        "novel_evidence_rows": len(evidence),
        "novel_species": int(evidence["accepted_species"].nunique()),
        "novel_species_trait": len(
            evidence.drop_duplicates(["accepted_species", "trait_name"])
        ),
        "evidence_by_trait": evidence["trait_name"].value_counts().sort_index().to_dict(),
        "audit_rows": len(audit),
        "robots_note": (
            "Original site robots.txt disallows automated crawling; no bulk refetch was performed. "
            "Evidence uses the archived snapshot and retains original record URLs."
        ),
        "quality_reason": (
            "Medium: direct database field/value with stable record ID, but consumed through an "
            "archived third-party snapshot rather than a fresh origin crawl."
        ),
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-csv", type=Path, required=True)
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--current-coverage-csv", type=Path, required=True)
    parser.add_argument("--current-lineages-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(
        json.dumps(
            build(
                args.source_csv,
                args.master_csv,
                args.current_coverage_csv,
                args.current_lineages_csv,
                args.output_dir,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
