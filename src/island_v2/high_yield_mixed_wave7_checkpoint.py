"""Freeze a high-yield, trait-specific direct-evidence checkpoint.

The seven rows target current support=2 genus x trait rules.  Queue potential is
diagnostic only: this module emits species-direct evidence, never genus values.
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

CREATED_AT = "2026-08-14T08:45:00Z"
SOURCE_GROUP = "high_yield_mixed_wave7_checkpoint_20260814"

SOURCES: tuple[dict[str, object], ...] = (
    {
        "species": "Betula pendula",
        "family": "Betulaceae",
        "trait": "mating_system",
        "value": "predominantly_outcrossing",
        "quality": "medium",
        "provider": "EUFORGEN species review",
        "url": "https://www.euforgen.org/species/betula-pendula",
        "title": "Betula pendula - EUFORGEN",
        "citation": "Chaplin (2024), EUFORGEN bibliographic species review",
        "excerpt": (
            "This diversity is attributed to silver birch being an outcrossing "
            "diploid species with its pollen able to travel hundreds of kilometres "
            "and its seeds being easily dispersed."
        ),
        "raw_value": "outcrossing diploid species",
        "record_id": "euforgen:betula-pendula:genetic-diversity:outcrossing",
        "lineage": "official_species_review:euforgen:betula_pendula:chaplin_2024",
        "lineage_method": "official_species_review_citation_lineage",
        "source_tier": "A",
        "source_type": "official_species_bibliographic_review",
        "domain": "euforgen.org",
        "language": "en",
        "content_sha256": "376038c5c15a7d25a7936a80ee5286c1f943eedfe30a4f05664c053ad0c13755",
        "content_sha256_basis": "downloaded_official_html_bytes",
        "potential": 44,
    },
    {
        "species": "Gastrodia kuroshimensis",
        "family": "Orchidaceae",
        "trait": "autonomous_selfing_capacity",
        "value": "autonomous",
        "quality": "high",
        "provider": "Phytotaxa original study",
        "url": "https://phytotaxa.mapress.com/pt/article/view/phytotaxa.278.3.6/8921",
        "title": (
            "Gastrodia kuroshimensis (Orchidaceae), a new mycoheterotrophic "
            "and complete cleistogamous plant from Japan"
        ),
        "citation": "Suetsugu (2016), Phytotaxa 278:265-272, DOI 10.11646/phytotaxa.278.3.6",
        "excerpt": (
            "Careful dissection of the G. kuroshimensis flowers at different "
            "stages revealed that their pollinia fragment into massulae before the "
            "flowers matured. The flowers were also found to have degenerated "
            "rostellums, which allowed the massulae to drop directly onto the "
            "surface of the stigma, confirming the obligate self-pollinating nature "
            "of this species."
        ),
        "raw_value": "pollinia drop directly onto stigma; obligate self-pollination",
        "record_id": "doi:10.11646/phytotaxa.278.3.6:p270:autonomous-mechanism",
        "lineage": "doi:10.11646/phytotaxa.278.3.6",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_reproductive_observation",
        "domain": "phytotaxa.mapress.com",
        "language": "en",
        "content_sha256": "dfddb42a2d414acb6d74d5f0df8f8e47b68429b92db26f4f37c598bfde83ebf3",
        "content_sha256_basis": "downloaded_publisher_pdf_bytes",
        "potential": 41,
    },
    {
        "species": "Goodyera tesselata",
        "family": "Orchidaceae",
        "trait": "self_incompatibility",
        "value": "SC",
        "quality": "medium",
        "provider": "New Jersey Natural Heritage Program rare plant profile",
        "url": (
            "https://www.nj.gov/dep/parksandforests/natural/heritage/docs/"
            "goodyera-tesselata-checkered-rattlesnake-plantain.pdf"
        ),
        "title": "Goodyera tesselata Rare Plant Profile",
        "citation": "Dodds (2022), New Jersey Department of Environmental Protection, p. 5",
        "excerpt": (
            "Kallunki (1976, 1981) found that G. tesselata did not set seed when "
            "insects were excluded. However, she also reported that the species is "
            "highly self-compatible and fruit development is not inhibited when "
            "flowers are fertilized by pollen from another flower on the same stem "
            "or ramet in the same clone."
        ),
        "raw_value": "highly self-compatible; no autonomous seed set",
        "record_id": "njdep:goodyera-tesselata:p5:self-compatible",
        "lineage": "official_species_review:njdep:dodds_2022:goodyera_tesselata",
        "lineage_method": "official_species_review_citing_primary_experiment",
        "source_tier": "A",
        "source_type": "official_government_species_review",
        "domain": "nj.gov",
        "language": "en",
        "content_sha256": "807478bfe88f9afeb475f7fe0c13cac4ae99402ecc60596614542b206edeecad",
        "content_sha256_basis": "exact_crawled_official_pdf_excerpt_utf8_bytes",
        "potential": 53,
    },
    {
        "species": "Puccinellia fasciculata",
        "family": "Poaceae",
        "trait": "self_incompatibility",
        "value": "SC",
        "quality": "medium",
        "provider": "New Jersey Natural Heritage Program rare plant profile",
        "url": (
            "https://dspace.njstatelib.org/server/api/core/bitstreams/"
            "88582153-6a46-43e5-95cb-19b30c4d4edc/content"
        ),
        "title": "Puccinellia fasciculata Rare Plant Profile",
        "citation": "Dodds (2022), New Jersey Department of Environmental Protection, p. 4",
        "excerpt": (
            "Connor (1988) reported self-compatibility in two varieties of P. "
            "fasciculata which are now included as synonyms of the species."
        ),
        "raw_value": "self-compatibility in two synonymized varieties",
        "record_id": "njdep:puccinellia-fasciculata:p4:self-compatible",
        "lineage": "official_species_review:njdep:dodds_2022:puccinellia_fasciculata",
        "lineage_method": "official_species_review_citing_primary_experiment",
        "source_tier": "A",
        "source_type": "official_government_species_review",
        "domain": "dspace.njstatelib.org",
        "language": "en",
        "content_sha256": "a8c18c1cdfc7fa708b55e13c886a92033d7e5bfcee3411e9ccc76046d52c64c2",
        "content_sha256_basis": "downloaded_official_pdf_bytes",
        "potential": 48,
    },
    {
        "species": "Daphne gnidium",
        "family": "Thymelaeaceae",
        "trait": "floral_symmetry",
        "value": "actinomorphic",
        "quality": "medium",
        "provider": "University of Evora Virtual Biodiversity Museum",
        "url": (
            "https://museubiodiversidade.uevora.pt/elenco-de-especies/"
            "biodiversidade-actual/plantas/angiospermicas/daphne-gnidium/"
        ),
        "title": "Daphne gnidium - Museu Virtual Biodiversidade",
        "citation": "University of Evora Virtual Biodiversity Museum species treatment",
        "excerpt": (
            "Flores: 4,5-6,5 mm, actinomorfas, hermafroditas, tetrâmeras e "
            "apétalas, florescendo de Julho a Outubro."
        ),
        "raw_value": "flores actinomorfas",
        "record_id": "uevora-museum:daphne-gnidium:flowers:symmetry",
        "lineage": "university_species_treatment:uevora_museum:daphne_gnidium",
        "lineage_method": "university_species_treatment",
        "source_tier": "A",
        "source_type": "university_virtual_biodiversity_museum",
        "domain": "museubiodiversidade.uevora.pt",
        "language": "pt",
        "content_sha256": "b495f058d173b4c5005111a0d7824fbfa69d9a173fbdee5f62257a52ba0a7c3c",
        "content_sha256_basis": "downloaded_university_html_bytes",
        "potential": 17,
    },
    {
        "species": "Ceodes grandis",
        "family": "Nyctaginaceae",
        "trait": "flower_primary_color",
        "value": "white",
        "quality": "medium",
        "provider": "Conservatoire Botanique National de Mascarin",
        "url": (
            "https://ileseparses.cbnm.org/index.php/plus/telechargements/send/"
            "32-guide-de-gestion-iles-eparses/374-hivert-2021-guide-de-gestion-"
            "de-15-especes-vegetales-menacees-d-europa-v2022-1"
        ),
        "title": "Guide de gestion de 15 espèces végétales menacées d'Europa",
        "citation": "Hivert (2021), CBNM management guide, Pisonia grandis treatment",
        "excerpt": (
            "Pisonia grandis. Synonyme : Ceodes grandis (R. Br.) D.Q. Lu. "
            "Fleurs nombreuses, de 3 mm de long, blanches."
        ),
        "raw_value": "fleurs blanches",
        "record_id": "cbnm-hivert-2021:pisonia-grandis:white-flowers",
        "lineage": "official_species_guide:cbnm:hivert_2021:ceodes_grandis",
        "lineage_method": "official_species_treatment_explicit_synonym",
        "source_tier": "A",
        "source_type": "official_conservation_botanical_guide",
        "domain": "ileseparses.cbnm.org",
        "language": "fr",
        "content_sha256": "96728b92c88ada2bc0ac9a48ab125b42a714cb9a8edbb4f0dd71ec10fe588060",
        "content_sha256_basis": "downloaded_official_pdf_bytes",
        "potential": 16,
        "matched_page_name": "Pisonia grandis",
        "evidence_scope": "synonym_direct",
        "name_match_method": "exact_synonym",
        "name_resolution_lineage": "source_explicit_synonym_to_master_accepted_name",
    },
    {
        "species": "Turnera diffusa",
        "family": "Turneraceae",
        "trait": "floral_symmetry",
        "value": "actinomorphic",
        "quality": "medium",
        "provider": "SEMARNAT non-timber forest resources manual",
        "url": "https://centro.paot.org.mx/documentos/semarnat/Manual_Clima_arido.pdf",
        "title": "Manual de recursos forestales no maderables de clima arido y semiarido",
        "citation": "SEMARNAT technical manual, section 4.4 Damiana, p. 70",
        "excerpt": (
            "Flores bisexuales actinomorfas, solitarias, axilares, de 2 a 12 mm "
            "de largo, sésiles."
        ),
        "raw_value": "flores actinomorfas",
        "record_id": "semarnat-manual:turnera-diffusa:p70:symmetry",
        "lineage": "official_species_manual:semarnat:turnera_diffusa",
        "lineage_method": "official_government_species_manual",
        "source_tier": "A",
        "source_type": "official_government_technical_manual",
        "domain": "centro.paot.org.mx",
        "language": "es",
        "content_sha256": "cf9fed1c32c8b0670d8e7d851058c3e2dcf66feb9be5deb3bded08f592f5bdc9",
        "content_sha256_basis": "exact_crawled_official_pdf_excerpt_utf8_bytes",
        "potential": 13,
    },
)


def reviewed_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for source in SOURCES:
        row = _evidence_row(
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
        row.update(
            {
                "source_group": SOURCE_GROUP,
                "query": "support2_third_species_exact_trait_high_yield",
                "language": str(source["language"]),
                "matched_page_name": str(
                    source.get("matched_page_name", source["species"])
                ),
                "evidence_scope": str(source.get("evidence_scope", "species_direct")),
                "name_match_method": str(
                    source.get("name_match_method", "accepted_name_exact")
                ),
                "name_resolution_lineage": str(
                    source.get(
                        "name_resolution_lineage", "master_accepted_name_exact"
                    )
                ),
            }
        )
        rows.append(row)
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
    identity = master.set_index("accepted_species")["family"].to_dict()
    expected = {str(source["species"]): str(source["family"]) for source in SOURCES}
    missing = sorted(set(expected) - set(identity))
    conflicts = {
        species: (identity.get(species, ""), family)
        for species, family in expected.items()
        if species in identity and identity[species] != family
    }
    if missing or conflicts:
        raise ValueError(f"master identity failure: missing={missing}, conflicts={conflicts}")

    evidence = pd.DataFrame(reviewed_rows(), columns=EVIDENCE_COLUMNS).fillna("")
    if evidence[["accepted_species", "trait_name"]].duplicated().any():
        raise ValueError("checkpoint species-trait pairs must be unique")
    if not evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all():
        raise ValueError("every checkpoint row requires a 64-character content hash")

    audit = _audit(evidence)
    audit["reviewer"] = "Codex high-yield trait-specific source audit"
    audit["reviewed_at_utc"] = CREATED_AT
    audit["decision_reason"] = (
        "Accepted after exact master identity or explicit synonym, exact trait, "
        "value polarity, source text, lineage, provenance and cultivar scope audit."
    )

    prior_evidence = pd.read_csv(prior_curated_evidence_csv, dtype=str).fillna("")
    prior_audit = pd.read_csv(prior_curated_audit_csv, dtype=str).fillna("")
    owned = prior_evidence["source_group"].eq(SOURCE_GROUP)
    prior_ids = set(prior_evidence.loc[owned, "candidate_id"])
    current_ids = set(evidence["candidate_id"])
    combined_evidence = pd.concat(
        [prior_evidence.loc[~owned], evidence], ignore_index=True
    )
    combined_audit = pd.concat(
        [
            prior_audit.loc[
                ~prior_audit["candidate_id"].isin(prior_ids | current_ids)
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
        "evidence": output_dir / "high_yield_mixed_wave7_evidence_20260814.csv",
        "audit": output_dir / "high_yield_mixed_wave7_manual_audit_20260814.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260814.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260814.csv",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    summary: dict[str, object] = {
        "contract": "high_yield_trait_specific_support2_checkpoint_v1",
        "created_at_utc": CREATED_AT,
        "accepted_rows": len(evidence),
        "species_trait": len(evidence),
        "quality_counts": evidence["evidence_quality"].value_counts().to_dict(),
        "trait_counts": evidence["trait_name"].value_counts().to_dict(),
        "language_counts": evidence["language"].value_counts().to_dict(),
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
    readme = output_dir / "README.md"
    hashed = [*paths.values(), *([readme] if readme.exists() else [])]
    for path in hashed:
        summary["files"][path.name] = {
            "sha256": _sha256(path),
            "size_bytes": len(_canonical_file_bytes(path)),
        }
    manifest = output_dir / "high_yield_mixed_wave7_manifest_20260814.json"
    manifest.write_text(
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
