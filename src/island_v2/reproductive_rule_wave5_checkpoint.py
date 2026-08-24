"""Freeze reviewed third-species evidence for reproductive genus rules.

The seven rows are exact target-master species and exact reproductive traits.
They were selected from support=2 acquisition-queue rules, then checked against
primary papers or official/curated species treatments.  Queue potential is kept
as a planning diagnostic only; this module emits no genus inference.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from island_v2.high_leverage_direct_checkpoint import (
    EVIDENCE_COLUMNS,
    _audit,
    _canonical_file_bytes,
    _evidence_row,
)

CREATED_AT = "2026-08-13T17:15:00Z"
SOURCE_GROUP = "reproductive_rule_wave5_checkpoint_20260814"


SOURCES: tuple[dict[str, object], ...] = (
    {
        "species": "Calanthe mannii",
        "family": "Orchidaceae",
        "trait": "self_incompatibility",
        "value": "SC",
        "quality": "medium",
        "provider": "King and Pantling primary orchid monograph",
        "url": "https://archive.org/details/mobot31753002133210",
        "title": "The Orchids of the Sikkim-Himalaya",
        "citation": "King & Pantling (1898), The Orchids of the Sikkim-Himalaya, p. 167",
        "excerpt": (
            "This is a self-fertile species. Mr. Pantling was unable to obtain "
            "pollinia for figuring, as in all the flowers of the only living "
            "specimen found by him, they had been more or less absorbed by the "
            "stigmas, the clinandrium having apparently become absorbed."
        ),
        "raw_value": "self-fertile; observation limited to the only living specimen",
        "record_id": "king-pantling-1898:p167:calanthe-mannii:self-fertile",
        "lineage": "monograph:king_pantling_1898_orchids_sikkim:calanthe_mannii",
        "lineage_method": "primary_monograph_species_treatment",
        "source_tier": "A",
        "source_type": "primary_botanical_monograph",
        "domain": "archive.org",
        "content_sha256": "1e27f565dd9770f5ca293916dc5a9bec7311280bf2f259be0564542ef8a602c7",
        "content_sha256_basis": "downloaded_mobot_scan_ocr_utf8_bytes",
        "potential": 148,
    },
    {
        "species": "Lotus hispidus",
        "family": "Fabaceae",
        "trait": "autonomous_selfing_capacity",
        "value": "autonomous",
        "quality": "medium",
        "provider": "USDA Agriculture Handbook 223",
        "url": (
            "https://www.govinfo.gov/content/pkg/GOVPUB-A-PURL-gpo20952/pdf/"
            "GOVPUB-A-PURL-gpo20952.pdf"
        ),
        "title": "The Trefoils: Adaptation and Culture",
        "citation": "Henson & Schoth (1962), USDA Agriculture Handbook 223, p. 4",
        "excerpt": "The annual tetraploid L. hispidus Desf. is apparently autogamous.",
        "raw_value": "apparently autogamous",
        "record_id": "usda-ah223:p4:lotus-hispidus:autogamous",
        "lineage": "official_monograph:usda_agriculture_handbook_223:lotus_hispidus",
        "lineage_method": "official_monograph_species_statement",
        "source_tier": "A",
        "source_type": "official_government_monograph",
        "domain": "govinfo.gov",
        "content_sha256": "78e780d6b8d9f15492ac21f0634274b105cbe79ec5c0ff225d02c5c79d0e91c0",
        "content_sha256_basis": "downloaded_official_pdf_bytes",
        "potential": 70,
    },
    {
        "species": "Fritillaria meleagris",
        "family": "Liliaceae",
        "trait": "mating_system",
        "value": "predominantly_outcrossing",
        "quality": "high",
        "provider": "Wiley original pollination study",
        "url": "https://doi.org/10.1111/j.1438-8677.2011.00510.x",
        "title": (
            "Neither protogynous nor obligatory out-crossed: pollination biology "
            "and breeding system of Fritillaria meleagris"
        ),
        "citation": "Zych & Stpiczynska (2012), Plant Biology 14:285-294",
        "excerpt": (
            "Selfing, although rare in natural populations, results in fully "
            "developed seeds. The species' dependence on generally rare "
            "pollinators and largely out-crossed breeding system may accelerate "
            "local extinction."
        ),
        "raw_value": "largely out-crossed; selfing rare but successful",
        "record_id": "doi:10.1111/j.1438-8677.2011.00510.x:abstract:mating-system",
        "lineage": "doi:10.1111/j.1438-8677.2011.00510.x",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_breeding_system_study",
        "domain": "onlinelibrary.wiley.com",
        "content_sha256": "4c1916bd5c7d7059614c8a0e731083ec08cdcedaf5fc326abcd639413e8fee03",
        "content_sha256_basis": "crossref_doi_metadata_abstract_utf8_bytes",
        "potential": 36,
    },
    {
        "species": "Fraxinus pennsylvanica",
        "family": "Oleaceae",
        "trait": "mating_system",
        "value": "predominantly_outcrossing",
        "quality": "medium",
        "provider": "USDA Forest Service FEIS",
        "url": "https://research.fs.usda.gov/feis/species-reviews/frapen",
        "title": "Fraxinus pennsylvanica, green ash",
        "citation": "Gucker (2005), USDA Forest Service Fire Effects Information System",
        "excerpt": (
            "Green ash trees are dioecious [110,124,269]. With male and female "
            "flowers on separate, wind-pollinated trees, outcrossing is mandatory "
            "and the potential for genetic exchange is great."
        ),
        "raw_value": "outcrossing is mandatory",
        "record_id": "usfs-feis:frapen:breeding-system",
        "lineage": "official_review:usfs_feis:frapen:gucker_2005",
        "lineage_method": "official_species_review_citation_lineage",
        "source_tier": "A",
        "source_type": "official_government_species_review",
        "domain": "research.fs.usda.gov",
        "content_sha256": "05d29badd931012b10f5b63086c35a8df8f4e0b9ddc8e723487bc647e0fc446e",
        "content_sha256_basis": "downloaded_official_html_utf8_bytes",
        "potential": 25,
    },
    {
        "species": "Asarum caulescens",
        "family": "Aristolochiaceae",
        "trait": "autonomous_selfing_capacity",
        "value": "autonomous",
        "quality": "high",
        "provider": "Plant Species Biology original bagging study",
        "url": "https://doi.org/10.1111/j.1442-1984.1987.tb00040.x",
        "title": "Self-pollination of Asarum caulescens Maxim. in Japan",
        "citation": "Tanaka & Yahara (1987), Plant Species Biology 2:133-136",
        "excerpt": (
            "Observations in two populations of Asarum caulescens belonging to "
            "sect. Asarum indicate that inbreeding predominates becuase (1) no "
            "effective pollinator was observed, (2) bagged flowers set fruits "
            "with well-swollen seeds, and (3) the behaviour of filaments, "
            "changing from recurved to straight posture, results in direct "
            "deposition of pollen grains on the stigmas."
        ),
        "raw_value": "bagged flowers set fruits with well-swollen seeds",
        "record_id": "doi:10.1111/j.1442-1984.1987.tb00040.x:bagged-flowers",
        "lineage": "doi:10.1111/j.1442-1984.1987.tb00040.x",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_pollinator_exclusion_study",
        "domain": "onlinelibrary.wiley.com",
        "content_sha256": "b3ee93efc197da75470b7662b3ddfe6bfe8b5ffb53c22f1d7cb903d76f373771",
        "content_sha256_basis": "downloaded_official_repository_pdf_bytes",
        "potential": 66,
    },
    {
        "species": "Polygonum aviculare",
        "family": "Polygonaceae",
        "trait": "mating_system",
        "value": "predominantly_selfing",
        "quality": "medium",
        "provider": "PROSEA species treatment via Pl@ntUse",
        "url": "https://plantuse.plantnet.org/en/Polygonum_aviculare_%28PROSEA%29?oldid=335207",
        "title": "Polygonum aviculare (PROSEA)",
        "citation": "PROSEA species treatment; reproductive section cites Yurtseva (1998)",
        "excerpt": (
            "Self-pollination in the open flower is the main pollination type. "
            "Outbreeding seems to be extremely rare."
        ),
        "raw_value": "self-pollination is the main type; outbreeding extremely rare",
        "record_id": "prosea:polygonum-aviculare:reproductive-biology:oldid-335207",
        "lineage": "species_treatment:prosea:polygonum_aviculare",
        "lineage_method": "canonical_species_treatment_revision",
        "source_tier": "A",
        "source_type": "curated_botanical_monograph_species_treatment",
        "domain": "plantuse.plantnet.org",
        "content_sha256": "602408ba3f5a3afc998c9f721deef3b7e69ce1538aec9e50e8e491ab888dfef7",
        "content_sha256_basis": "downloaded_fixed_revision_html_utf8_bytes",
        "potential": 57,
    },
    {
        "species": "Ribes nigrum",
        "family": "Grossulariaceae",
        "trait": "self_incompatibility",
        "value": "SC",
        "quality": "medium",
        "provider": "Pladias Database of the Czech Flora and Vegetation",
        "url": "https://pladias.cz/en/taxon/data/Ribes%20nigrum",
        "title": "Ribes nigrum - Pladias",
        "citation": "Pladias species trait record; pollination syndrome sourced to BIOLFLOR",
        "excerpt": "Pollination syndrome: insect-pollination, selfing",
        "raw_value": "selfing",
        "record_id": "pladias:ribes-nigrum:pollination-syndrome:selfing",
        "lineage": "trait_database:biolflor_pladias:ribes_nigrum:pollination_syndrome",
        "lineage_method": "database_feature_original_lineage",
        "source_tier": "A",
        "source_type": "curated_flora_trait_database",
        "domain": "pladias.cz",
        "content_sha256": "a0321c9fafcae041e4613e0894bcfb2a365308cc25d4b7d6f97ca8bfc18e2d56",
        "content_sha256_basis": "downloaded_species_page_html_utf8_bytes",
        "potential": 61,
    },
)


def reviewed_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for source in SOURCES:
        evidence = _evidence_row(
            species=str(source["species"]),
            trait=str(source["trait"]),
            value=str(source["value"]),
            quality=str(source["quality"]),
            provider=str(source["provider"]),
            url=str(source["url"]),
            title=str(source["title"]),
            citation=str(source["citation"]),
            excerpt=str(source["excerpt"]),
            record_id=str(source["record_id"]),
            lineage=str(source["lineage"]),
            lineage_method=str(source["lineage_method"]),
            source_tier=str(source["source_tier"]),
            source_type=str(source["source_type"]),
            domain=str(source["domain"]),
            content_sha256=str(source["content_sha256"]),
            content_sha256_basis=str(source["content_sha256_basis"]),
            retrieved_at_utc=CREATED_AT,
            raw_value=str(source["raw_value"]),
        )
        evidence["source_group"] = SOURCE_GROUP
        evidence["query"] = "support2_reproductive_third_species_primary_official"
        evidence["name_resolution_lineage"] = "master_accepted_name_exact"
        rows.append(evidence)
    return rows


def _sha256(path: Path) -> str:
    return hashlib.sha256(_canonical_file_bytes(path)).hexdigest()


def build(
    *,
    master_csv: Path,
    prior_curated_evidence_csv: Path,
    prior_curated_audit_csv: Path,
    output_dir: Path,
) -> dict[str, object]:
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    master_identity = master.set_index("accepted_species")["family"].to_dict()
    expected_identity = {str(s["species"]): str(s["family"]) for s in SOURCES}
    missing = sorted(set(expected_identity) - set(master_identity))
    if missing:
        raise ValueError(f"checkpoint species missing from master: {missing}")
    conflicts = {
        species: (master_identity[species], family)
        for species, family in expected_identity.items()
        if master_identity[species] != family
    }
    if conflicts:
        raise ValueError(f"checkpoint family conflicts: {conflicts}")

    evidence = pd.DataFrame(reviewed_rows(), columns=EVIDENCE_COLUMNS).fillna("")
    if evidence[["accepted_species", "trait_name"]].duplicated().any():
        raise ValueError("checkpoint species-trait pairs must be unique")
    if not evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all():
        raise ValueError("every checkpoint row requires a 64-character content hash")

    audit = _audit(evidence)
    audit["reviewer"] = "Codex primary and official reproductive evidence audit"
    audit["reviewed_at_utc"] = CREATED_AT
    audit["decision_reason"] = (
        "Accepted after exact target-master identity, exact reproductive trait, "
        "value polarity, full source text, lineage, provenance and cultivar "
        "scope were checked. Queue potential was not treated as coverage."
    )

    prior_evidence = pd.read_csv(prior_curated_evidence_csv, dtype=str).fillna("")
    prior_audit = pd.read_csv(prior_curated_audit_csv, dtype=str).fillna("")
    owned = prior_evidence["source_group"].eq(SOURCE_GROUP)
    prior_owned_ids = set(prior_evidence.loc[owned, "candidate_id"])
    current_ids = set(evidence["candidate_id"])
    combined_evidence = pd.concat(
        [prior_evidence.loc[~owned], evidence], ignore_index=True
    )
    combined_audit = pd.concat(
        [
            prior_audit.loc[
                ~prior_audit["candidate_id"].isin(prior_owned_ids | current_ids)
            ],
            audit,
        ],
        ignore_index=True,
    )
    for label, frame in (("evidence", combined_evidence), ("audit", combined_audit)):
        if frame["candidate_id"].duplicated().any():
            raise ValueError(f"combined {label} candidate IDs must be unique")

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "evidence": output_dir / "reproductive_rule_wave5_evidence_20260814.csv",
        "audit": output_dir / "reproductive_rule_wave5_manual_audit_20260814.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260814.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260814.csv",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    summary: dict[str, object] = {
        "contract": "reproductive_support2_third_species_checkpoint_v1",
        "created_at_utc": CREATED_AT,
        "accepted_rows": len(evidence),
        "species_trait": len(evidence),
        "quality_counts": evidence["evidence_quality"].value_counts().to_dict(),
        "trait_counts": evidence["trait_name"].value_counts().to_dict(),
        "source_lineages": int(evidence["source_lineage"].nunique()),
        "rule_candidates": len(SOURCES),
        "theoretical_queue_cells_before_formal_validation": sum(
            int(source["potential"]) for source in SOURCES
        ),
        "theoretical_only_not_formal_coverage": True,
        "combined": {
            "evidence_rows": len(combined_evidence),
            "audit_rows": len(combined_audit),
        },
        "guardrails": {
            "trait_specific_records": True,
            "genus_inference_emitted_here": False,
            "family_inference": False,
            "global_fallback": False,
            "n2_formal_inference": False,
            "cross_trait_substitution": False,
            "search_snippet_evidence": False,
            "queue_potential_counted_as_coverage": False,
        },
        "inputs": {
            "master_csv": {"path": str(master_csv), "sha256": _sha256(master_csv)},
            "prior_curated_evidence_csv": {
                "path": str(prior_curated_evidence_csv),
                "sha256": _sha256(prior_curated_evidence_csv),
            },
            "prior_curated_audit_csv": {
                "path": str(prior_curated_audit_csv),
                "sha256": _sha256(prior_curated_audit_csv),
            },
        },
        "files": {},
    }
    readme_path = output_dir / "README.md"
    hash_paths = list(paths.values())
    if readme_path.exists():
        hash_paths.append(readme_path)
    for path in hash_paths:
        summary["files"][path.name] = {
            "sha256": _sha256(path),
            "size_bytes": len(_canonical_file_bytes(path)),
        }
    manifest_path = output_dir / "reproductive_rule_wave5_manifest_20260814.json"
    manifest_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--prior-curated-evidence-csv", type=Path, required=True)
    parser.add_argument("--prior-curated-audit-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    print(json.dumps(build(**vars(parser.parse_args())), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
