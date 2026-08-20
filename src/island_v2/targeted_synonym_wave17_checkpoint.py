"""Build Wave 17 from strict BSdb synonyms and one regional-flora treatment.

Only reviewed species/synonym-direct rows are appended here.  Validated Low is
still rebuilt exclusively by the shared all-evidence implementation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.bsdb_two_backbone_synonym_checkpoint import build_review_audit
from island_v2.targeted_support2_wave15_checkpoint import (
    EVIDENCE_COLUMNS,
    _candidate_id,
    _sha,
)

SOURCE_GROUP = "targeted_synonym_wave17_checkpoint_20260820"
CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "targeted_synonym_wave17_checkpoint_20260820"
)
PRIOR = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "targeted_support2_wave16_checkpoint_20260820"
)
RETRIEVED_AT = "2026-08-20T08:00:00Z"


def regional_rows(path: Path) -> list[dict[str, str]]:
    source = pd.read_csv(path, dtype=str).fillna("")
    if len(source) != 2 or set(source["accepted_species"]) != {
        "Monopsis stellarioides"
    }:
        raise ValueError("Wave17 regional snapshot must contain two Monopsis rows")
    rows: list[dict[str, str]] = []
    for record in source.to_dict("records"):
        rows.append(
            {
                "candidate_id": _candidate_id(
                    record["accepted_species"],
                    record["trait_name"],
                    record["source_lineage"],
                    record["normalized_value"],
                ),
                "accepted_species": record["accepted_species"],
                "axis": record["axis"],
                "trait_name": record["trait_name"],
                "raw_value": record["raw_value"],
                "normalized_value": record["normalized_value"],
                "evidence_quality": "medium",
                "evidence_scope": "species_direct",
                "source_group": SOURCE_GROUP,
                "source_provider": record["source_provider"],
                "source_url": record["source_url"],
                "page_title": record["page_title"],
                "source_citation": record["source_citation"],
                "source_excerpt": record["source_excerpt"],
                "source_record_id": record["source_record_id"],
                "source_lineage": record["source_lineage"],
                "lineage_method": record["lineage_method"],
                "name_resolution_lineage": "fixed_master_accepted_species_family_exact",
                "name_match_method": "accepted_name_exact",
                "matched_page_name": record["accepted_species"],
                "source_tier": "B",
                "source_type": record["source_type"],
                "domain": record["domain"],
                "language": record["language"],
                "wild_cultivated_cultivar_status": record[
                    "wild_cultivated_cultivar_status"
                ],
                "evidence_status": "accepted_individual_source_backed_audit",
                "content_sha256": record["content_sha256"].casefold(),
                "content_sha256_basis": record["content_sha256_basis"],
                "retrieved_at_utc": record["retrieved_at_utc"],
                "query": record["query"],
                "search_rank": "",
                "inference_rule": "",
            }
        )
    return rows


def build_checkpoint(output_dir: Path = CHECKPOINT) -> dict[str, Any]:
    bsdb = pd.read_csv(
        output_dir / "bsdb_two_backbone_evidence_20260820.csv", dtype=str
    ).fillna("")
    if len(bsdb) != 10 or bsdb[["accepted_species", "trait_name"]].drop_duplicates().shape[0] != 8:
        raise ValueError("Wave17 BSdb snapshot must contain 10 rows and 8 species traits")
    bsdb["source_group"] = SOURCE_GROUP
    regional = pd.DataFrame(
        regional_rows(output_dir / "regional_flora_source_rows_wave17.csv"),
        columns=EVIDENCE_COLUMNS,
    )
    evidence = pd.concat([bsdb[EVIDENCE_COLUMNS], regional], ignore_index=True)
    evidence = evidence.sort_values(
        ["accepted_species", "trait_name", "source_lineage"], kind="stable"
    ).reset_index(drop=True)
    if len(evidence) != 12:
        raise ValueError("Wave17 must contain exactly 12 reviewed rows")
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("Wave17 candidate IDs must be unique")
    duplicated_pairs = evidence.loc[
        evidence.duplicated(["accepted_species", "trait_name"], keep=False),
        ["accepted_species", "trait_name"],
    ].drop_duplicates()
    if duplicated_pairs.to_dict("records") != [
        {
            "accepted_species": "Doona cordifolia",
            "trait_name": "self_incompatibility",
        }
    ]:
        raise ValueError("only the reviewed Doona direct conflict may repeat a pair")
    audit = build_review_audit(evidence)

    prior_evidence = pd.read_csv(
        PRIOR / "combined_curated_evidence_20260820.csv", dtype=str
    ).fillna("")
    prior_audit = pd.read_csv(
        PRIOR / "combined_curated_manual_audit_20260820.csv", dtype=str
    ).fillna("")
    combined_evidence = pd.concat([prior_evidence, evidence], ignore_index=True)
    combined_audit = pd.concat([prior_audit, audit], ignore_index=True)
    if combined_evidence["candidate_id"].duplicated().any():
        raise ValueError("combined evidence candidate IDs must be unique")
    if combined_audit["candidate_id"].duplicated().any():
        raise ValueError("combined audit candidate IDs must be unique")

    paths = {
        "evidence": output_dir / "targeted_synonym_wave17_evidence_20260820.csv",
        "audit": output_dir / "targeted_synonym_wave17_manual_audit_20260820.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260820.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260820.csv",
        "manifest": output_dir / "source_acquisition_manifest_wave17.json",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    mapping = pd.read_csv(
        output_dir / "bsdb_two_backbone_mapping_audit_20260820.csv.gz", dtype=str
    ).fillna("")
    selection = pd.read_csv(
        output_dir / "bsdb_two_backbone_selection_audit_20260820.csv", dtype=str
    ).fillna("")
    manifest = {
        "checkpoint": SOURCE_GROUP,
        "built_at_utc": RETRIEVED_AT,
        "baseline_formal_run_id": 32334209990,
        "accepted_evidence_rows": len(evidence),
        "accepted_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "bsdb_acquisition": {
            "source_rows": 3609,
            "unresolved_unique_names_queried": len(mapping),
            "wfo_requests": len(mapping),
            "gbif_requests": len(mapping),
            "total_web_api_requests": 2 * len(mapping),
            "strict_two_backbone_agreements": int(
                mapping["mapping_reason"].eq("accepted_strict_two_backbone").sum()
            ),
            "strict_mapping_rate": round(
                mapping["mapping_reason"].eq("accepted_strict_two_backbone").mean(),
                8,
            ),
            "new_direct_rows": int(selection["selection_reason"].eq("selected").sum()),
            "new_direct_species_trait": int(
                bsdb[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
            ),
            "direct_conflict_species_trait": 1,
            "api_cost_usd": 0.0,
            "frozen_response_jsonl_sha256": (
                "898ebdc9f3bf57958f34a3d368e02fd8c476d08684700d8cfd548c2e0f2cd1d3"
            ),
            "pinned_bsdb_csv_sha256": (
                "f935c8150b3d719aba5f62f14335c1a185304155403bf50db1e3ef1393fc55f4"
            ),
        },
        "regional_flora_resume": {
            "completed_species_pages_not_researched": 1207,
            "new_species_pages_fetched": 333,
            "new_reviewed_direct_rows": 2,
            "new_reviewed_species_trait": 2,
            "source_domain": "zambiaflora.com",
            "search_api_cost_usd": 0.0,
        },
        "recorded_queries": 2766,
        "search_cost_usd": 0.0,
        "guardrails": {
            "both_backbones_same_accepted_species": True,
            "fixed_master_family_exact": True,
            "fuzzy_or_one_backbone_match": False,
            "search_snippet_as_evidence": False,
            "family_inference": False,
            "global_fallback": False,
            "min_species_two_production": False,
            "cross_trait_substitution": False,
            "genus_axis_only_join": False,
            "direct_conflict_forced_resolved": False,
        },
        "output_sha256": {
            key: _sha(path.read_text(encoding="utf-8"))
            for key, path in paths.items()
            if key != "manifest"
        },
        "notes": (
            "Doona cordifolia retains three independent High rows with SC/SI "
            "disagreement and must remain an unresolved direct conflict. The formal "
            "all-evidence rebuild, not this checkpoint, determines strict and Low gains."
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
