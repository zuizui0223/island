"""Build the Wave 23 exact-species dioecy checkpoint.

The checkpoint records only the species-direct mating-system consequence of
dioecy: separate male and female plants require outcrossing.  It deliberately
does not substitute dioecy for self-incompatibility or autonomous selfing.
Validated Low remains the output of the shared trait-specific all-evidence
rebuild, never of this checkpoint.
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
    _sha,
)
from island_v2.targeted_support2_wave15_checkpoint import (
    _row as _wave15_row,
)
from island_v2.targeted_support2_wave15_checkpoint import (
    build_audit as _wave15_build_audit,
)

SOURCE_GROUP = "targeted_dioecy_wave23_checkpoint_20260821"
CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "targeted_dioecy_wave23_checkpoint_20260821"
)
PRIOR = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "ferrer_2024_powo_synonym_checkpoint_20260821"
)
RETRIEVED_AT = "2026-08-21T06:00:00Z"


def _row(*args: str, **kwargs: str) -> dict[str, str]:
    row = _wave15_row(*args, **kwargs)
    row["source_group"] = SOURCE_GROUP
    row["retrieved_at_utc"] = RETRIEVED_AT
    return row


def _dioecy_row(
    species: str,
    quality: str,
    provider: str,
    url: str,
    title: str,
    citation: str,
    excerpt: str,
    lineage: str,
    tier: str,
    source_type: str,
    query: str,
    content_sha256: str,
    content_sha256_basis: str,
) -> dict[str, str]:
    return _row(
        species,
        "reproductive_assurance",
        "mating_system",
        "dioecious; male and female flowers occur on separate plants",
        "predominantly_outcrossing",
        quality,
        provider,
        url,
        title,
        citation,
        excerpt,
        lineage,
        tier,
        source_type,
        "en",
        query,
        content_sha256=content_sha256,
        content_sha256_basis=content_sha256_basis,
    )


def primary_rows() -> list[dict[str, str]]:
    """Return 11 individually reviewed, exact-species dioecy records."""
    rows = [
        _dioecy_row(
            "Ilex opaca",
            "high",
            "USDA Forest Service Silvics",
            "https://www.srs.fs.usda.gov/pubs/misc/ag_654/volume_2/ilex/opaca.htm",
            "Ilex opaca Ait. - American Holly",
            "USDA Forest Service. Silvics of North America, Ilex opaca.",
            (
                "Hollies are dioecious; male (staminate) and female "
                "(pistillate) flowers, similar in appearance, with four to six "
                "small white petals, are produced on separate plants on the "
                "current season's growth. Pollination is accomplished by insects."
            ),
            "official-species-treatment:usda-fs-silvics:Ilex_opaca",
            "A",
            "official_public_agency_species_synthesis",
            '"Ilex opaca" dioecious male female separate plants',
            "e92aa6325c66f885d55ca535b40f4ef9ad99de22b82b899fed8c01f3cafa31a9",
            "retrieved_official_species_page_bytes",
        ),
        _dioecy_row(
            "Ilex crenata",
            "high",
            "American Journal of Botany / PubMed",
            "https://pubmed.ncbi.nlm.nih.gov/26199373/",
            (
                "Sexual dimorphism in floral longevity and flowering synchrony "
                "in relation to pollination and mating success in three "
                "dioecious Ilex species"
            ),
            "Matsuhisa & Ushimaru 2015. DOI 10.3732/ajb.1500073.",
            (
                "We examined how sexual dimorphism in floral longevity and "
                "flowering synchrony affects pollination and mating success in "
                "three dioecious Ilex species (I. pedunculosa, I. serrata, and "
                "I. crenata)."
            ),
            "doi:10.3732/ajb.1500073",
            "A",
            "primary_article_species_level_field_study",
            '"Ilex crenata" dioecious pollination mating success',
            "8feb4220c5396ca3ce16547e024c1b8651b4afddcd961e01b9d80d516b28f3b9",
            "retrieved_official_pubmed_primary_abstract_html_bytes",
        ),
        _dioecy_row(
            "Litsea glutinosa",
            "high",
            "Songklanakarin Journal of Science and Technology",
            "https://sjst.psu.ac.th/index.php/article/1920/PaperPDF/download",
            "Floral biology and pollination ecology of Litsea glutinosa",
            "Bhadra & Bandyopadhyay 2019. DOI 10.14456/sjst-psu.2019.4.",
            (
                "Litsea glutinosa is a semi-evergreen wet season blooming tree "
                "species. It is a dioecious plant characterized by separate "
                "staminate and pistillate trees occurring in a 3:1 ratio."
            ),
            "doi:10.14456/sjst-psu.2019.4",
            "A",
            "primary_article_pollination_ecology",
            '"Litsea glutinosa" dioecious staminate pistillate trees',
            "4227fb870d16c8063a08afcc88d10841e7f10821ea284eb9668f82e932fd8838",
            "retrieved_original_article_pdf_bytes",
        ),
        _dioecy_row(
            "Litsea cubeba",
            "high",
            "Genes / MDPI",
            (
                "https://res.mdpi.com/d_attachment/genes/genes-11-00048/"
                "article_deploy/genes-11-00048.pdf"
            ),
            "Chromosome-level genome assembly of Litsea cubeba",
            "Chen et al. 2020. Genes 11:48. DOI 10.3390/genes11010048.",
            (
                "Litsea cubeba (Lour.) Pers., a popular essential oil plant, "
                "is a dioecious species with degenerative sexual organs in both "
                "male and female individuals."
            ),
            "doi:10.3390/genes11010048",
            "A",
            "primary_article_species_genome_and_sex_determination",
            '"Litsea cubeba" dioecious male female individuals',
            "e1fb64cb0b2f446c3462f59ad1eb0d03c00bcbc2969d32247b73fd7278b2a340",
            "retrieved_original_article_pdf_bytes",
        ),
        _dioecy_row(
            "Nepenthes rafflesiana",
            "high",
            "National Parks Board Singapore Flora & Fauna Web",
            "https://www.nparks.gov.sg/florafaunaweb/flora/1/4/1459",
            "Nepenthes rafflesiana - Flora & Fauna Web",
            "NParks Singapore. Exact-species Flora & Fauna Web treatment.",
            "Plant is diocecious, with male and female flowers found on separate plants.",
            "provider-treatment:nparks:flora:1459",
            "A",
            "official_government_species_database",
            '"Nepenthes rafflesiana" dioecious separate plants',
            "805ea4659302b022971e1385a2e1b06b9ad59a2c5cef297a21d6de272992409b",
            "retrieved_official_species_page_bytes",
        ),
        _dioecy_row(
            "Nepenthes mirabilis",
            "high",
            "PLOS ONE",
            (
                "https://journals.plos.org/plosone/article/file?"
                "id=10.1371/journal.pone.0322885&type=printable"
            ),
            "A high quality genome of the common swamp pitcher plant",
            "Clarke et al. 2025. PLOS ONE. DOI 10.1371/journal.pone.0322885.",
            (
                "Here, PacBio HiFi long-read sequencing was used to assemble a "
                "near chromosome-level genome for Nepenthes mirabilis. Moreover, "
                "whereas most plants are functional hermaphrodites, Nepenthes is "
                "dioecious - that is, individual plants are either male or female."
            ),
            "doi:10.1371/journal.pone.0322885",
            "A",
            "primary_article_species_genome_and_sex_chromosome",
            '"Nepenthes mirabilis" dioecious male female',
            "551f42d6d943a1c0cdf3dd204428ed454ec64135f5b8d52456b6836093c923f5",
            "retrieved_original_article_pdf_bytes",
        ),
        _dioecy_row(
            "Myrsine australis",
            "medium",
            "New Zealand Plant Conservation Network",
            "https://www.nzpcn.org.nz/site/assets/files/0/15/453/myrsine-australis.pdf",
            "Myrsine australis fact sheet",
            (
                "New Zealand Plant Conservation Network species fact sheet; "
                "description based on Allan 1961."
            ),
            (
                "Fixed female and inconstant male flowers on different plants, "
                "1.5-2.5 mm diam."
            ),
            "provider_treatment:nzpcn:myrsine-australis",
            "B",
            "specialist_national_species_profile",
            '"Myrsine australis" flowers different plants',
            "2a85c1453b67899ce738d1557d64300aa79bf8c648048d39cdf829d098d33d08",
            "retrieved_specialist_species_factsheet_pdf_bytes",
        ),
        _dioecy_row(
            "Myrsine divaricata",
            "medium",
            "University of Auckland New Zealand Plants",
            (
                "https://www.nzplants.auckland.ac.nz/en/about/seed-plants-"
                "flowering/primulaceae/myrsine-divaricata.html"
            ),
            "Myrsine divaricata - New Zealand Plants",
            "University of Auckland. Exact-species profile.",
            "Sexuality: unisexual on different plants.",
            "institutional-species-profile:university-auckland-nzplants:Myrsine_divaricata",
            "B",
            "university_exact_species_profile",
            '"Myrsine divaricata" unisexual different plants',
            "b000fe5e5d677bcaaf9bbba0d83fde09f5406117b85acb3620058aa19a5d8088",
            "retrieved_university_species_page_bytes",
        ),
        _dioecy_row(
            "Myrsine umbricola",
            "high",
            "Biota of New Zealand",
            (
                "https://biotanz.landcareresearch.co.nz/scientific-names/"
                "369eb446-eefe-457a-a0f7-7b29a9ae0540"
            ),
            "Myrsine umbricola - Biota of New Zealand",
            (
                "Heenan & de Lange 2004. New Zealand Journal of Botany. "
                "DOI 10.1080/0028825X.2004.9512929."
            ),
            "Evergreen, spreading, dioecious shrub.",
            "doi:10.1080/0028825X.2004.9512929",
            "A",
            "primary_taxonomic_species_description_official_database",
            '"Myrsine umbricola" dioecious shrub',
            "953a749d975b071d171954d1694a506e93f7458131d8ee4cfacfa846233e6a9b",
            "retrieved_official_biodiversity_species_page_bytes",
        ),
        _dioecy_row(
            "Asparagus cochinchinensis",
            "high",
            "PLOS ONE",
            "https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0266376",
            "Genome-wide identification of sex-related genes in Asparagus",
            "Murase et al. 2022. PLOS ONE. DOI 10.1371/journal.pone.0266376.",
            (
                "The dioecious species A. cochinchinensis, A. officinalis, and "
                "A. schoberioides formed a monophyletic group."
            ),
            "doi:10.1371/journal.pone.0266376",
            "A",
            "primary_article_comparative_genomics",
            '"Asparagus cochinchinensis" dioecious',
            "7e55471bd46cdff066495d49216d828194b100e18a3fa67cc6e1503203032098",
            "retrieved_original_article_html_bytes",
        ),
        _dioecy_row(
            "Asparagus schoberioides",
            "high",
            "PLOS ONE",
            "https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0266376",
            "Genome-wide identification of sex-related genes in Asparagus",
            "Murase et al. 2022. PLOS ONE. DOI 10.1371/journal.pone.0266376.",
            (
                "The dioecious species A. cochinchinensis, A. officinalis, and "
                "A. schoberioides formed a monophyletic group."
            ),
            "doi:10.1371/journal.pone.0266376",
            "A",
            "primary_article_comparative_genomics",
            '"Asparagus schoberioides" dioecious',
            "7e55471bd46cdff066495d49216d828194b100e18a3fa67cc6e1503203032098",
            "retrieved_original_article_html_bytes",
        ),
    ]
    return rows


def build_audit(evidence: pd.DataFrame) -> pd.DataFrame:
    audit = _wave15_build_audit(evidence)
    audit["reviewer"] = "Codex Wave23 exact-species dioecy source audit"
    audit["reviewed_at_utc"] = RETRIEVED_AT
    audit["decision_reason"] = (
        "Accepted after exact fixed-master species/family match, original-page "
        "retrieval, quote and source-lineage review, and cultivar screening. "
        "Dioecy is recorded only as mating_system=predominantly_outcrossing; "
        "it is not substituted for self-incompatibility or autonomous selfing."
    )
    return audit.loc[:, AUDIT_COLUMNS]


def build_checkpoint(output_dir: Path = CHECKPOINT) -> dict[str, Any]:
    evidence = pd.DataFrame(primary_rows(), columns=EVIDENCE_COLUMNS).sort_values(
        ["accepted_species", "trait_name", "source_lineage"], kind="stable"
    )
    evidence = evidence.reset_index(drop=True)
    audit = build_audit(evidence)
    if len(evidence) != 11 or len(audit) != 11:
        raise ValueError("Wave23 must contain exactly 11 individually reviewed rows")
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("Wave23 candidate IDs must be unique")
    if evidence.duplicated(["accepted_species", "trait_name"]).any():
        raise ValueError("Wave23 species x trait pairs must be unique")

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
        "evidence": output_dir / "targeted_dioecy_wave23_evidence_20260821.csv",
        "audit": output_dir / "targeted_dioecy_wave23_manual_audit_20260821.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260821.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260821.csv",
        "manifest": output_dir / "source_acquisition_manifest_wave23.json",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    manifest = {
        "checkpoint": SOURCE_GROUP,
        "built_at_utc": RETRIEVED_AT,
        "baseline_formal_run_id": 32449801624,
        "accepted_evidence_rows": len(evidence),
        "accepted_species_trait": len(evidence),
        "accepted_species": int(evidence["accepted_species"].nunique()),
        "accepted_source_lineages": int(evidence["source_lineage"].nunique()),
        "quality_counts": evidence["evidence_quality"].value_counts().to_dict(),
        "recorded_queries": int(evidence["query"].nunique()),
        "formal_search_api_queries": 0,
        "search_cost_usd": 0.0,
        "targeted_support2_rules": [
            "Asparagus|mating_system",
            "Ilex|mating_system",
            "Litsea|mating_system",
            "Myrsine|mating_system",
            "Nepenthes|mating_system",
        ],
        "theoretical_rule_cells_touched": 831,
        "guardrails": {
            "search_snippet_as_evidence": False,
            "family_inference": False,
            "global_fallback": False,
            "min_species_two_production": False,
            "cross_trait_substitution": False,
            "genus_axis_only_join": False,
            "dioecy_mapped_to_self_incompatibility": False,
            "dioecy_mapped_to_autonomous_selfing": False,
            "cultivar_or_hybrid_transferred_to_wild_species": False,
        },
        "output_sha256": {
            key: _sha(path.read_text(encoding="utf-8"))
            for key, path in paths.items()
            if key != "manifest"
        },
        "notes": (
            "The 831 cells are a queue ceiling, not observed coverage gain. "
            "The shared all-evidence rebuild must independently pass dominance, "
            "masked species validation, and source-lineage leave-one-out validation."
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
