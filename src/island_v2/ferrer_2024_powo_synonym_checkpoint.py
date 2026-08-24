"""Promote Ferrer 2024 rows rescued by strict WFO and POWO agreement."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.targeted_support2_wave15_checkpoint import (
    AUDIT_COLUMNS,
    EVIDENCE_COLUMNS,
    _sha,
)
from island_v2.targeted_support2_wave15_checkpoint import (
    _row as _wave15_row,
)
from island_v2.targeted_support2_wave15_checkpoint import (
    build_audit as _wave15_build_audit,
)

SOURCE_GROUP = "ferrer_2024_powo_synonym_checkpoint_20260821"
SOURCE_DIR = Path("data/v2/external/traits/ferrer_2024_two_backbone")
SOURCE_ROWS = SOURCE_DIR / "ferrer_selected_source_rows_20260821.csv"
MAPPINGS = SOURCE_DIR / "powo_wfo_selected_mappings_20260821.csv"
RESPONSES = SOURCE_DIR / "ferrer_two_backbone_responses_20260821.jsonl.gz"
MAPPING_AUDIT = SOURCE_DIR / "ferrer_two_backbone_mapping_audit_20260821.csv.gz"
CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "ferrer_2024_powo_synonym_checkpoint_20260821"
)
PRIOR = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "goodwillie_2005_mixed_mating_checkpoint_20260821"
)
MASTER = Path("data/v2/staging/gbif/collected/island_taxa.csv")
RETRIEVED_AT = "2026-08-21T05:08:53Z"
SOURCE_URL = "https://doi.org/10.1093/aob/mcae056"
WORKBOOK_SHA256 = "8ef2f5ac99780ca19a15b847442272457634064577f8a12fbb9710f5521e5899"
RESPONSES_SHA256 = "8fd6076a5895c4b1844f082126f342f653b38c9376b1bb67eea501132f1c3660"
MAPPING_AUDIT_SHA256 = "c778928d169e28d61bbdb0b49889600a25d91edc6e42f0bf78a9578143b10cfb"


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if _file_sha256(RESPONSES) != RESPONSES_SHA256:
        raise ValueError("Ferrer WFO/GBIF response snapshot hash mismatch")
    if _file_sha256(MAPPING_AUDIT) != MAPPING_AUDIT_SHA256:
        raise ValueError("Ferrer two-backbone mapping audit hash mismatch")
    source = pd.read_csv(SOURCE_ROWS, dtype=str).fillna("")
    mappings = pd.read_csv(MAPPINGS, dtype=str).fillna("")
    master = pd.read_csv(MASTER, dtype=str).fillna("")
    if len(source) != 3 or len(mappings) != 3:
        raise ValueError("Ferrer WFO/POWO checkpoint requires exactly three rows")
    if source["source_scientific_name"].duplicated().any():
        raise ValueError("Ferrer selected source names must be unique")
    if mappings["source_name"].duplicated().any():
        raise ValueError("WFO/POWO mapping names must be unique")

    master_family = dict(zip(master["accepted_species"], master["family"], strict=False))
    for row in mappings.itertuples(index=False):
        if master_family.get(row.accepted_species) != row.source_family:
            raise ValueError(f"fixed-master family mismatch for {row.source_name}")
        if "Synonym of: " + row.accepted_species not in row.powo_visible_result:
            raise ValueError(f"POWO accepted name not explicit for {row.source_name}")
        if not row.wfo_id.startswith("wfo-") or row.wfo_classification_version != "2026-06":
            raise ValueError(f"invalid WFO snapshot identity for {row.source_name}")
    return source, mappings, master


def primary_rows() -> list[dict[str, str]]:
    """Return three synonym-direct Medium self-compatibility records."""
    source, mappings, _ = _load_inputs()
    source_by_name = source.set_index("source_scientific_name").to_dict("index")
    rows: list[dict[str, str]] = []
    for mapping in mappings.to_dict("records"):
        source_name = mapping["source_name"]
        record = source_by_name[source_name]
        if record["standardized_value"] != "SC":
            raise ValueError(f"unexpected Ferrer value for {source_name}")
        lineage = record["source_lineage"]
        resolution = (
            f"wfo:{mapping['wfo_classification_version']}:{mapping['wfo_id']};"
            f"powo:{mapping['powo_source_name_lsid']}:accepted:"
            f"{mapping['powo_accepted_name_lsid']};fixed_master_species_family_exact"
        )
        evidence = _wave15_row(
            mapping["accepted_species"],
            "reproductive_assurance",
            "self_incompatibility",
            record["source_value"],
            record["standardized_value"],
            "medium",
            "Ferrer et al. 2024 global self-incompatibility database",
            SOURCE_URL,
            "A global database of plant self-incompatibility systems",
            (
                "Ferrer et al. (2024), Annals of Botany, DOI 10.1093/aob/mcae056; "
                f"underlying source: {record['source_citation']}"
            ),
            record["evidence_excerpt"],
            lineage,
            "A",
            "peer_reviewed_curated_literature_database",
            "en",
            f'Ferrer 2024 exact synonym "{source_name}" WFO POWO',
            matched_name=source_name,
            scope="synonym_direct",
            name_match_method="synonym_exact_two_backbone",
            name_resolution_lineage=resolution,
            status="published_species_level_record_not_cultivar_limited",
            content_sha256=WORKBOOK_SHA256,
            content_sha256_basis="original_oup_supplementary_workbook_bytes",
        )
        evidence.update(
            {
                "source_group": SOURCE_GROUP,
                "source_record_id": (
                    f"ferrer2024:powo-synonym:{record['source_record_id']}:"
                    "self_incompatibility"
                ),
                "lineage_method": "underlying_literature_citation_not_database_provider",
                "domain": "academic.oup.com",
                "retrieved_at_utc": RETRIEVED_AT,
            }
        )
        rows.append(evidence)
    return rows


def build_audit(evidence: pd.DataFrame) -> pd.DataFrame:
    audit = _wave15_build_audit(evidence)
    audit["reviewer"] = "Codex Ferrer WFO-POWO strict synonym audit"
    audit["reviewed_at_utc"] = RETRIEVED_AT
    audit["decision_reason"] = (
        "Accepted as a cited Ferrer species row only after WFO June 2026 and Kew "
        "POWO/WCVP independently mapped the exact synonym to the same fixed-master "
        "species and the source, POWO target, and master families agreed."
    )
    return audit.loc[:, AUDIT_COLUMNS]


def build_checkpoint(output_dir: Path = CHECKPOINT) -> dict[str, Any]:
    evidence = pd.DataFrame(primary_rows(), columns=EVIDENCE_COLUMNS).sort_values(
        ["accepted_species", "trait_name"], kind="stable"
    )
    evidence = evidence.reset_index(drop=True)
    audit = build_audit(evidence)
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("Ferrer WFO/POWO candidate IDs must be unique")

    prior_evidence = pd.read_csv(
        PRIOR / "combined_curated_evidence_20260821.csv", dtype=str
    ).fillna("")
    prior_audit = pd.read_csv(
        PRIOR / "combined_curated_manual_audit_20260821.csv", dtype=str
    ).fillna("")
    combined_evidence = pd.concat([prior_evidence, evidence], ignore_index=True)
    combined_audit = pd.concat([prior_audit, audit], ignore_index=True)
    if combined_evidence["candidate_id"].duplicated().any():
        raise ValueError("combined evidence candidate IDs must be unique")
    if combined_audit["candidate_id"].duplicated().any():
        raise ValueError("combined audit candidate IDs must be unique")

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "evidence": output_dir / "ferrer_2024_powo_synonym_evidence_20260821.csv",
        "audit": output_dir / "ferrer_2024_powo_synonym_manual_audit_20260821.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260821.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260821.csv",
        "manifest": output_dir / "source_acquisition_manifest_ferrer_2024_powo.json",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    mapping_audit = pd.read_csv(MAPPING_AUDIT, dtype=str).fillna("")
    mapping_counts = mapping_audit["mapping_reason"].value_counts().sort_index().to_dict()
    manifest = {
        "checkpoint": SOURCE_GROUP,
        "built_at_utc": RETRIEVED_AT,
        "baseline_formal_run_id": 32447600643,
        "source_rows": 5686,
        "safe_direct_rows_before_name_matching": 3706,
        "prior_direct_candidate_species": 2427,
        "unmatched_safe_rows": 1181,
        "unique_unmatched_names_queried": 1171,
        "wfo_requests": 1171,
        "reused_gbif_responses": 1171,
        "powo_manual_result_pages_reviewed": 3,
        "search_api_cost_usd": 0.0,
        "accepted_evidence_rows": len(evidence),
        "accepted_species_trait": len(evidence),
        "accepted_source_lineages": int(evidence["source_lineage"].nunique()),
        "wfo_gbif_mapping_counts": mapping_counts,
        "theoretical_rule_candidate": {
            "genus": "Chaetogastra",
            "trait_name": "self_incompatibility",
            "new_value": "SC",
            "prior_support_species": 2,
            "potential_unresolved_cells": 12,
            "formal_gain_claimed": False,
        },
        "source_sha256": {
            "oup_supplementary_workbook": WORKBOOK_SHA256,
            RESPONSES.name: RESPONSES_SHA256,
            MAPPING_AUDIT.name: MAPPING_AUDIT_SHA256,
        },
        "guardrails": {
            "single_backbone_match_accepted": False,
            "fuzzy_match_accepted": False,
            "family_conflict_accepted": False,
            "family_inference": False,
            "global_fallback": False,
            "min_species_two_production": False,
            "cross_trait_substitution": False,
            "genus_axis_only_join": False,
        },
        "output_sha256": {
            key: _sha(path.read_text(encoding="utf-8"))
            for key, path in paths.items()
            if key != "manifest"
        },
        "notes": (
            "The 12-cell rule value is a queue ceiling, not a coverage claim. "
            "Only the formal all-evidence rebuild may report strict or Low gain."
        ),
    }
    paths["manifest"].write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=CHECKPOINT)
    args = parser.parse_args()
    print(json.dumps(build_checkpoint(args.output_dir), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
