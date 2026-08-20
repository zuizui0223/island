"""Build the Wave 16 reviewed direct-evidence checkpoint.

Wave 16 promotes only individually reviewed rows that were already present in
completed acquisition artifacts plus two newly retrieved institutional or
peer-reviewed species pages.  It never creates genus inference itself; the
shared all-evidence rebuild remains the sole producer of Validated Low rows.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.targeted_support2_wave15_checkpoint import (
    AUDIT_COLUMNS,
    EVIDENCE_COLUMNS,
    _candidate_id,
    _sha,
)
from island_v2.targeted_support2_wave15_checkpoint import (
    _row as _wave15_row,
)
from island_v2.targeted_support2_wave15_checkpoint import (
    build_audit as _wave15_build_audit,
)

SOURCE_GROUP = "targeted_support2_wave16_checkpoint_20260820"
CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/targeted_support2_wave16_checkpoint_20260820"
)
PRIOR = Path("data/v2/staging/traits/open_web_pilot/targeted_support2_wave15_checkpoint_20260820")


def _row(*args: str, **kwargs: str) -> dict[str, str]:
    row = _wave15_row(*args, **kwargs)
    row["source_group"] = SOURCE_GROUP
    row["retrieved_at_utc"] = "2026-08-20T12:30:00Z"
    return row


def primary_rows() -> list[dict[str, str]]:
    """Return newly retrieved, exact-species reproductive evidence."""
    return [
        _row(
            "Rhamnus alaternus",
            "reproductive_assurance",
            "mating_system",
            "reproduction type dioecious; flowers unisexual",
            "predominantly_outcrossing",
            "medium",
            "University of Porto Natural History and Science Museum",
            "https://www.mhnc.up.pt/pt/colecoes/colecao-botanica-viva/aderno-bastardo/",
            "Rhamnus alaternus L.",
            "University of Porto living botanical collection species account.",
            "Tipo de reprodução: Dióica; Organização sexual: Unissexuais.",
            "institutional-species-page:up-mhnc:Rhamnus_alaternus",
            "A",
            "university_museum_botanical_species_page",
            "pt",
            '"Rhamnus alaternus" dioecious flowers',
            content_sha256="4763748cc8f47740c773c4a8e91e2f3a262fb5b7166a7ffd59215f396b2fbb19",
            content_sha256_basis="retrieved_original_species_page_html_bytes",
        ),
        _row(
            "Viscum coloratum",
            "reproductive_assurance",
            "mating_system",
            "the plant is dioecious",
            "predominantly_outcrossing",
            "medium",
            "Biomolecules / PubMed Central",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC12292218/",
            "Viscum coloratum: A Review of Botany, Phytochemistry, Pharmacology, Pharmacokinetics and Toxicology",
            "Di et al. 2025. Biomolecules 15:974. DOI 10.3390/biom15070974.",
            "The plant is dioecious.",
            "doi:10.3390/biom15070974",
            "A",
            "peer_reviewed_species_botany_review",
            "en",
            '"Viscum coloratum" dioecious',
            content_sha256="95543993c10ec2a67dcc9fa0f54e27019064e642166dd520310be9ea6f58b37b",
            content_sha256_basis="retrieved_pmc_full_article_html_bytes",
        ),
    ]


def artifact_rows(path: Path) -> list[dict[str, str]]:
    """Convert pinned completed-artifact rows to the common evidence schema."""
    source = pd.read_csv(path, dtype=str).fillna("")
    if len(source) != 25:
        raise ValueError("Wave16 artifact snapshot must contain exactly 25 rows")
    providers = {
        "plants.ces.ncsu.edu": "NC State Extension Plant Toolbox",
        "plantnet.rbgsyd.nsw.gov.au": "PlantNET NSW Flora",
    }
    citations = {
        "plants.ces.ncsu.edu": "NC State Extension Plant Toolbox species page.",
        "plantnet.rbgsyd.nsw.gov.au": (
            "Royal Botanic Garden Sydney. PlantNET NSW Flora species treatment."
        ),
    }
    rows: list[dict[str, str]] = []
    for record in source.to_dict("records"):
        domain = record["domain"]
        quality = "high" if record["source_tier"] == "A" else "medium"
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
                "evidence_quality": quality,
                "evidence_scope": "species_direct",
                "source_group": SOURCE_GROUP,
                "source_provider": providers[domain],
                "source_url": record["source_url"],
                "page_title": record["page_title"],
                "source_citation": citations[domain],
                "source_excerpt": record["supporting_excerpt"],
                "source_record_id": record["candidate_id"],
                "source_lineage": record["source_lineage"],
                "lineage_method": record["lineage_method"],
                "name_resolution_lineage": "master_accepted_name_and_family_exact",
                "name_match_method": "accepted_name_exact",
                "matched_page_name": record["matched_page_name"],
                "source_tier": record["source_tier"],
                "source_type": record["source_type"],
                "domain": domain,
                "language": record["language"],
                "wild_cultivated_cultivar_status": ("species_treatment_not_cultivar_limited"),
                "evidence_status": "accepted_individual_source_backed_audit",
                "content_sha256": record["content_sha256"].casefold(),
                "content_sha256_basis": "retrieved_original_species_page_bytes",
                "retrieved_at_utc": record["retrieved_at_utc"],
                "query": record["query"],
                "search_rank": record["search_rank"],
                "inference_rule": "",
            }
        )
    return rows


def build_audit(evidence: pd.DataFrame) -> pd.DataFrame:
    audit = _wave15_build_audit(evidence)
    audit["reviewer"] = "Codex Wave16 individual source audit"
    audit["reviewed_at_utc"] = "2026-08-20T12:30:00Z"
    return audit.loc[:, AUDIT_COLUMNS]


def build_checkpoint(output_dir: Path = CHECKPOINT) -> dict[str, Any]:
    evidence = pd.DataFrame(
        artifact_rows(output_dir / "source_artifact_rows_wave16.csv") + primary_rows(),
        columns=EVIDENCE_COLUMNS,
    ).sort_values(["accepted_species", "trait_name", "source_lineage"], kind="stable")
    evidence = evidence.reset_index(drop=True)
    audit = build_audit(evidence)
    if len(evidence) != 27 or len(audit) != 27:
        raise ValueError("Wave16 must contain exactly 27 individually reviewed rows")
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("Wave16 candidate IDs must be unique")
    if evidence.duplicated(["accepted_species", "trait_name"]).any():
        raise ValueError("Wave16 species x trait pairs must be unique")

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
        "evidence": output_dir / "targeted_support2_wave16_evidence_20260820.csv",
        "audit": output_dir / "targeted_support2_wave16_manual_audit_20260820.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260820.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260820.csv",
        "manifest": output_dir / "source_acquisition_manifest_wave16.json",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    artifact_manifest = (
        pd.read_csv(output_dir / "source_artifact_rows_wave16.csv", dtype=str)
        .groupby(["source_run_id", "source_artifact_id", "source_artifact_name"])
        .size()
        .rename("accepted_rows")
        .reset_index()
        .to_dict("records")
    )
    targeted_rules = [
        "Yucca|flower_size_class",
        "Planchonella|floral_symmetry",
        "Rhamnus|mating_system",
        "Viscum|mating_system",
    ]
    manifest = {
        "checkpoint": SOURCE_GROUP,
        "built_at_utc": "2026-08-20T12:30:00Z",
        "source_artifacts": artifact_manifest,
        "source_artifact_snapshot_sha256": _sha(
            (output_dir / "source_artifact_rows_wave16.csv").read_text(encoding="utf-8-sig")
        ),
        "accepted_species_trait_rows": len(evidence),
        "accepted_species": int(evidence["accepted_species"].nunique()),
        "recorded_queries": int(evidence["query"].nunique()),
        "formal_search_api_queries": 0,
        "search_cost_usd": 0.0,
        "targeted_support2_rules": targeted_rules,
        "targeted_support2_rule_count": len(targeted_rules),
        "guardrails": {
            "search_snippet_as_evidence": False,
            "family_inference": False,
            "global_fallback": False,
            "min_species_two_production": False,
            "cross_trait_substitution": False,
            "genus_axis_only_join": False,
            "pollen_vector_or_reward_mixed_into_structure": False,
            "component_size_promoted_without_source_audit": False,
            "solitary_spikelet_mapped_to_solitary_flower_display": False,
        },
        "notes": (
            "The formal all-evidence rebuild determines realized direct and Low gains. "
            "Seven PlantNET solitary-spikelet rows and coarse NCSU less-than-one-inch "
            "size rows were deliberately excluded."
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
    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(json.dumps(build_checkpoint(args.output_dir), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
