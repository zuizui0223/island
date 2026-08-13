"""Freeze independent source lineages that recover high-value genus rules.

Every row is species-direct and was re-read in the retrieved source.  This
checkpoint emits no genus inference; the shared all-evidence implementation
must rebuild genus x trait rules, apply dominance and both leave-one-out gates,
and resolve direct conflicts after the rows are promoted.
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

CREATED_AT = "2026-08-13T00:32:40Z"
SOURCE_GROUP = "independent_lineage_recovery_checkpoint_20260813"


ROWS = [
    {
        "species": "Syzygium jambos",
        "trait": "floral_symmetry",
        "value": "actinomorphic",
        "quality": "high",
        "provider": "Bernice Pauahi Bishop Museum, Plants of Hawai'i",
        "url": (
            "https://plantsofhawaii.org/detail/"
            "%7BF83E4802-A6EC-49F2-98A0-9B18D74CEDE8%7D"
        ),
        "title": "Syzygium jambos (L.) Alston",
        "citation": "Bernice Pauahi Bishop Museum, Plants of Hawai'i species profile",
        "excerpt": "Flowers bisexual (perfect), actinomorphic.",
        "record_id": "bishop-museum:plants-of-hawaii:Syzygium-jambos:flowers",
        "lineage": (
            "provider_treatment:bishop_museum_plants_of_hawaii:Syzygium_jambos"
        ),
        "lineage_method": "institutional_species_treatment_lineage",
        "source_tier": "A",
        "source_type": "museum_species_database",
        "domain": "plantsofhawaii.org",
        "content_sha256": (
            "de172db9edbc2dc2ccb4213ac51a34715058bd4f09b27d00b8aeb42ba016570d"
        ),
        "content_sha256_basis": "retrieved_original_species_page_bytes",
        "raw_value": "actinomorphic",
    },
    {
        "species": "Erigeron canadensis",
        "trait": "autonomous_selfing_capacity",
        "value": "autonomous",
        "quality": "high",
        "provider": "Zelaya et al. 2004, PubMed-indexed primary article",
        "url": (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
            "db=pubmed&id=15502914&retmode=xml"
        ),
        "title": "Inheritance of evolved glyphosate resistance in Conyza canadensis",
        "citation": (
            "Zelaya, Owen & VanGessel (2004), Theoretical and Applied Genetics "
            "110:58-70, DOI 10.1007/s00122-004-1804-8"
        ),
        "excerpt": "The autogamous nature of C. canadensis",
        "record_id": "pubmed:15502914:abstract:autogamous-nature",
        "lineage": "doi:10.1007/s00122-004-1804-8",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_article",
        "domain": "eutils.ncbi.nlm.nih.gov",
        "content_sha256": (
            "200aef6a739523442629d5c615f0b70892790bbd28d48df2d0d633ee50842d49"
        ),
        "content_sha256_basis": "retrieved_pubmed_article_xml_bytes",
        "raw_value": "autogamous",
        "matched_page_name": "Conyza (=Erigeron) canadensis",
        "name_match_method": "exact_synonym",
        "name_resolution_lineage": (
            "primary_article_explicit_equation:Conyza_(=Erigeron)_canadensis"
        ),
    },
    {
        "species": "Aster kantoensis",
        "trait": "mating_system",
        "value": "predominantly_outcrossing",
        "quality": "high",
        "provider": "Oxford Academic, Journal of Heredity",
        "url": "https://oup.silverchair-cdn.com/article-minimal/2186764",
        "title": (
            "Inbreeding depression and outcrossing rate in the endangered "
            "autotetraploid plant Aster kantoensis"
        ),
        "citation": (
            "Inoue, Masuda & Maki (1998), Journal of Heredity 89:559-562, "
            "DOI 10.1093/jhered/89.6.559"
        ),
        "excerpt": (
            "The multilocus estimate of outcrossing rate in the population is "
            "0.883 +/- 0.064. The high value of inbreeding depression is "
            "consistent with the predominant outcrossing of the population."
        ),
        "record_id": "doi:10.1093/jhered/89.6.559:abstract:outcrossing-rate",
        "lineage": "doi:10.1093/jhered/89.6.559",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_article",
        "domain": "oup.silverchair-cdn.com",
        "content_sha256": (
            "b3f5792bcbf764e5914965b044481217c4437670a9e498bd004168ad7fd00d88"
        ),
        "content_sha256_basis": "retrieved_publisher_abstract_page_bytes",
        "raw_value": "multilocus outcrossing rate 0.883 +/- 0.064",
    },
    {
        "species": "Convolvulus sepium",
        "trait": "mating_system",
        "value": "predominantly_outcrossing",
        "quality": "medium",
        "provider": "Natuurtijdschriften original journal repository",
        "url": "https://natuurtijdschriften.nl/pub/551742",
        "title": "Over de middelen tot verspreiding van Calystegia (Convolvulus L.) Sepium R. Br.",
        "citation": "Vuyck (1895), Nederlandsch Kruidkundig Archief",
        "excerpt": "C. sepium est une plante exclusivement xenogame",
        "record_id": "natuurtijdschriften:551742:exclusive-xenogamy",
        "lineage": "vuyck1895:convolvulus_sepium",
        "lineage_method": "original_historical_article_lineage",
        "source_tier": "A",
        "source_type": "historical_primary_natural_history_article",
        "domain": "natuurtijdschriften.nl",
        "content_sha256": (
            "455e00010055eb482b4f0c37f279f183de0ecf6491dfe7182c1dfa0027a41b23"
        ),
        "content_sha256_basis": "retrieved_original_repository_record_bytes",
        "raw_value": "exclusivement xenogame",
        "language": "fr",
    },
    {
        "species": "Lolium temulentum",
        "trait": "autonomous_selfing_capacity",
        "value": "autonomous",
        "quality": "medium",
        "provider": "Aberystwyth University author manuscript",
        "url": (
            "https://pure.aber.ac.uk/ws/portalfiles/portal/29213471/"
            "Thomas_2019_Plants_People_Planet.pdf"
        ),
        "title": "Grass blindness",
        "citation": "Thomas (2019), Plants, People, Planet, DOI 10.1002/ppp3.28",
        "excerpt": (
            "L. temulentum is one of a group of three weedy self-fertile "
            "(autogamous) species."
        ),
        "record_id": "doi:10.1002/ppp3.28:lolium-temulentum:autogamous",
        "lineage": "doi:10.1002/ppp3.28",
        "lineage_method": "peer_reviewed_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_species_review",
        "domain": "pure.aber.ac.uk",
        "content_sha256": (
            "80ddb603bbdabb55e074365374e3090851a8e6cbad779bcbc97b7f1a2a14adcb"
        ),
        "content_sha256_basis": "retrieved_author_manuscript_pdf_bytes",
        "raw_value": "self-fertile (autogamous)",
    },
    {
        "species": "Festuca rubra",
        "trait": "autonomous_selfing_capacity",
        "value": "autonomous",
        "quality": "high",
        "provider": "Universidad Nacional de Rio Negro institutional repository",
        "url": (
            "https://rid.unrn.edu.ar/bitstream/20.500.12049/3523/1/"
            "Gundel%20%282014%29%20Fungal%20endophyte%20mediated%20occurrence%20"
            "of%20seminiferous%20and%20pseudoviviparous%20panicles%20in%20"
            "Festuca%20rubra.pdf"
        ),
        "title": (
            "Fungal endophyte mediated occurrence of seminiferous and "
            "pseudoviviparous panicles in Festuca rubra"
        ),
        "citation": (
            "Gundel et al. (2014), Fungal Diversity 66:69-76, "
            "DOI 10.1007/s13225-014-0290-9"
        ),
        "excerpt": (
            "Two panicles in plants with four or more panicles were individually "
            "closed in waxed bags. In contrast to the presumption that F. rubra "
            "is strictly cross pollinated, 38 out of 664 bagged panicles produced "
            "seeds revealing that self-pollination frequently occur."
        ),
        "record_id": "doi:10.1007/s13225-014-0290-9:bagged-panicles:38-of-664",
        "lineage": "doi:10.1007/s13225-014-0290-9",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_pollinator_exclusion_experiment",
        "domain": "rid.unrn.edu.ar",
        "content_sha256": (
            "fc2a77c48a96a2099f05b6755875ad484f37bc398c23ca589407cb3858de9dc0"
        ),
        "content_sha256_basis": "retrieved_institutional_repository_pdf_bytes",
        "raw_value": "38 of 664 bagged panicles produced seeds",
        "wild_status": "wild_source_populations_common_garden_not_cultivar_limited",
    },
    {
        "species": "Lactuca canadensis",
        "trait": "autonomous_selfing_capacity",
        "value": "autonomous",
        "quality": "medium",
        "provider": "Iowa State University Digital Repository",
        "url": (
            "https://dr.lib.iastate.edu/server/api/core/bitstreams/"
            "f7cd5881-f4d4-4a06-9432-c019cac76a01/content"
        ),
        "title": "Pollinators maintain biodiversity in assembling plant communities",
        "citation": (
            "Soley & Wilsey (2026), Ecology 107:e70369, DOI 10.1002/ecy.70369; "
            "species statement cites Parrish & Bazzaz (1979)"
        ),
        "excerpt": (
            "Three volunteer species were not treated because their seed "
            "production is not dependent on pollinators: Conyza canadensis and "
            "Lactuca canadensis, primarily autogamous annual forbs."
        ),
        "record_id": "doi:10.1002/ecy.70369:lactuca-canadensis:autogamy",
        "lineage": "citation:parrish_bazzaz_1979:lactuca_canadensis_autogamy",
        "lineage_method": "underlying_citation_lineage_not_republisher_url",
        "source_tier": "A",
        "source_type": "peer_reviewed_article_species_statement",
        "domain": "dr.lib.iastate.edu",
        "content_sha256": (
            "a8af7243a64a2ddf1048852ec4ff8a799aa286fc34c3e307fc9fd5841925c3e4"
        ),
        "content_sha256_basis": "retrieved_institutional_repository_pdf_bytes",
        "raw_value": "primarily autogamous",
    },
    {
        "species": "Nephelium lappaceum",
        "trait": "self_incompatibility",
        "value": "SI",
        "quality": "high",
        "provider": "Oxford Academic, Botanical Journal of the Linnean Society",
        "url": "https://oup.silverchair-cdn.com/article-minimal/2660262",
        "title": (
            "Reproductive patterns of selected understorey trees in the "
            "Malaysian rain forest: the sexual species"
        ),
        "citation": (
            "Ha, Sands, Soepadmo & Jong (1988), Botanical Journal of the "
            "Linnean Society 97:295-316, DOI 10.1111/j.1095-8339.1988.tb01585.x"
        ),
        "excerpt": "the self-incompatible androdioecious Nephelium lappaceum",
        "record_id": (
            "doi:10.1111/j.1095-8339.1988.tb01585.x:abstract:nephelium-si"
        ),
        "lineage": "doi:10.1111/j.1095-8339.1988.tb01585.x",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_article",
        "domain": "oup.silverchair-cdn.com",
        "content_sha256": (
            "32006aa19aa73828341a2134ef6a4c32e948d506f39d218597dbe0735be7d040"
        ),
        "content_sha256_basis": "retrieved_publisher_abstract_page_bytes",
        "raw_value": "self-incompatible",
    },
]


def reviewed_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for item in ROWS:
        evidence = _evidence_row(
            species=item["species"],
            trait=item["trait"],
            value=item["value"],
            quality=item["quality"],
            provider=item["provider"],
            url=item["url"],
            title=item["title"],
            citation=item["citation"],
            excerpt=item["excerpt"],
            record_id=item["record_id"],
            lineage=item["lineage"],
            lineage_method=item["lineage_method"],
            source_tier=item["source_tier"],
            source_type=item["source_type"],
            domain=item["domain"],
            content_sha256=item["content_sha256"],
            content_sha256_basis=item["content_sha256_basis"],
            retrieved_at_utc=CREATED_AT,
            raw_value=item["raw_value"],
        )
        evidence["source_group"] = SOURCE_GROUP
        evidence["query"] = "corrected_lineage_loo_high_value_rule_recovery"
        for optional in (
            "language",
            "matched_page_name",
            "name_match_method",
            "name_resolution_lineage",
        ):
            if optional in item:
                evidence[optional] = item[optional]
        if "wild_status" in item:
            evidence["wild_cultivated_cultivar_status"] = item["wild_status"]
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
    evidence = pd.DataFrame(reviewed_rows(), columns=EVIDENCE_COLUMNS).fillna("")
    missing = sorted(set(evidence["accepted_species"]) - set(master["accepted_species"]))
    if missing:
        raise ValueError(f"independent-lineage species missing from master: {missing}")
    if evidence[["accepted_species", "trait_name"]].duplicated().any():
        raise ValueError("independent-lineage species-trait pairs are not unique")
    if not evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all():
        raise ValueError("every independent-lineage row needs a retrieved content hash")
    audit = _audit(evidence)
    audit["reviewer"] = "Codex source-backed independent-lineage evidence audit"
    audit["reviewed_at_utc"] = CREATED_AT
    audit["decision_reason"] = (
        "Accepted species or exact source-stated synonym identity, exact source "
        "statement, trait-specific mapping, source URL and retrieved-content hash "
        "rechecked; no cultivar, family inference, cross-trait substitution or "
        "global fallback."
    )

    prior_evidence = pd.read_csv(prior_curated_evidence_csv, dtype=str).fillna("")
    prior_audit = pd.read_csv(prior_curated_audit_csv, dtype=str).fillna("")
    prior_owned = prior_evidence["source_group"].eq(SOURCE_GROUP)
    prior_owned_ids = set(prior_evidence.loc[prior_owned, "candidate_id"])
    owned_ids = set(evidence["candidate_id"])
    combined_evidence = pd.concat(
        [prior_evidence.loc[~prior_owned], evidence], ignore_index=True
    )
    combined_audit = pd.concat(
        [
            prior_audit.loc[
                ~prior_audit["candidate_id"].isin(prior_owned_ids | owned_ids)
            ],
            audit,
        ],
        ignore_index=True,
    )
    for label, frame in (("evidence", combined_evidence), ("audit", combined_audit)):
        if frame["candidate_id"].duplicated().any():
            raise ValueError(f"combined {label} candidate IDs are not unique")

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "evidence": output_dir / "independent_lineage_recovery_evidence_20260813.csv",
        "audit": output_dir
        / "independent_lineage_recovery_manual_audit_20260813.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260813.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260813.csv",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    summary: dict[str, object] = {
        "contract": "independent_lineage_recovery_checkpoint_v1",
        "created_at_utc": CREATED_AT,
        "accepted_rows": len(evidence),
        "species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "quality_counts": evidence["evidence_quality"].value_counts().to_dict(),
        "trait_counts": evidence["trait_name"].value_counts().to_dict(),
        "source_lineages": int(evidence["source_lineage"].nunique()),
        "theoretical_rule_trait_cells": 1_020,
        "theoretical_only_not_formal_coverage": True,
        "corrective_rows": {
            "Lolium_rigidum_cross_species_exclusion": True,
            "Nephelium_lappaceum_medium_SC_to_high_SI": True,
        },
        "explicitly_rejected": [
            "Calamus javensis|flower_size_class|petal_size_not_whole_flower_size",
            "Cirsium brevistylum|autonomous_selfing_capacity|apomixis_not_excluded",
            "Durio zibethinus|self_incompatibility|cultivar_specific_variability",
        ],
        "combined": {
            "evidence_rows": len(combined_evidence),
            "audit_rows": len(combined_audit),
        },
        "guardrails": {
            "trait_specific_records": True,
            "genus_inference_emitted_here": False,
            "family_inference": False,
            "global_fallback": False,
            "cross_trait_substitution": False,
            "search_snippet_evidence": False,
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
    for path in paths.values():
        summary["files"][path.name] = {
            "sha256": _sha256(path),
            "size_bytes": len(_canonical_file_bytes(path)),
        }
    manifest_path = output_dir / "independent_lineage_recovery_manifest_20260813.json"
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
