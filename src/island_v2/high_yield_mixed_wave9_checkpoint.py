"""Freeze the ninth high-yield species-direct evidence checkpoint.

The checkpoint targets current support=2 genus x trait_name rules and retains
counterevidence that can invalidate an unsafe rule.  Queue potential is
diagnostic only: this module emits direct records and never emits genus values.
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

CREATED_AT = "2026-08-13T23:42:09Z"
SOURCE_GROUP = "high_yield_mixed_wave9_checkpoint_20260814"

SOURCES: tuple[dict[str, object], ...] = (
    {
        "species": "Berberis microphylla",
        "family": "Berberidaceae",
        "trait": "mating_system",
        "value": "predominantly_outcrossing",
        "quality": "high",
        "provider": "Acta Horticulturae primary article",
        "url": (
            "https://repositorio.unimoron.edu.ar/bitstream/10.34073/188/1/"
            "Effect%20of%20different%20pollination%20treatments%20on%20Berberis.pdf"
        ),
        "title": (
            "Effect of different pollination treatments on Berberis microphylla "
            "G. Forst, a Patagonian barberry"
        ),
        "citation": (
            "Radice and Arena (2019), Acta Horticulturae 1231:75-80, "
            "DOI 10.17660/ActaHortic.2019.1231.13"
        ),
        "excerpt": (
            "While several authors defined B. microphylla as self-fertile, this "
            "species is more likely cross-pollinated and Syrphidae insects are "
            "essential in pollen transport. Experiments developed so far support "
            "this hypothesis."
        ),
        "raw_value": "more likely cross-pollinated; cross treatment best in three years",
        "record_id": "doi:10.17660/actahortic.2019.1231.13:conclusions",
        "lineage": "doi:10.17660/actahortic.2019.1231.13",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_pollination_experiment",
        "domain": "repositorio.unimoron.edu.ar",
        "language": "en",
        "content_sha256": "b75b272ba08dcc3b043973344a379f12b5a425701e791804ce0afc49b66ffc74",
        "content_sha256_basis": "downloaded_official_university_repository_pdf_bytes",
        "potential": 56,
        "role": "rule_candidate",
    },
    {
        "species": "Pectis linifolia",
        "family": "Asteraceae",
        "trait": "autonomous_selfing_capacity",
        "value": "autonomous",
        "quality": "medium",
        "provider": "University of Arizona Herbarium regional flora",
        "url": (
            "https://herbarium.arizona.edu/sites/default/files/2024-08/"
            "Part%2021.%20Asteraceae%20%E2%80%93%20Aster%20Family.pdf"
        ),
        "title": (
            "Ajo Peak to Tinajas Altas: A Flora of Southwestern Arizona, "
            "Part 21 Asteraceae"
        ),
        "citation": "Felger and Rutman (2016), Phytoneuron 2016-77:1-164, p. 115",
        "excerpt": (
            "Heads of 5 small ray florets and 1-3 disk florets, rather "
            "inconspicuous, self-fertile (autogamous) as evidenced by low number "
            "of florets per head, reduction in ray size (1 mm long), and small "
            "anthers (less than 1 mm long) with low pollen production."
        ),
        "raw_value": "self-fertile (autogamous), inferred from floral syndrome",
        "record_id": "felger-rutman-2016:asteraceae:p115:pectis-linifolia-autogamous",
        "lineage": (
            "official_regional_flora:felger_rutman_2016_sw_arizona_asteraceae:"
            "pectis_linifolia"
        ),
        "lineage_method": "official_species_treatment_citation_lineage",
        "source_tier": "A",
        "source_type": "academic_regional_flora_species_treatment",
        "domain": "herbarium.arizona.edu",
        "language": "en",
        "content_sha256": "b3029888a372b95550d959a8411aecce17814d265a37df86c42a6f7b05ed8192",
        "content_sha256_basis": "downloaded_official_university_flora_pdf_bytes",
        "potential": 28,
        "role": "rule_candidate",
    },
    {
        "species": "Callistemon phoeniceus",
        "family": "Myrtaceae",
        "trait": "autonomous_selfing_capacity",
        "value": "mixed_or_variable",
        "quality": "high",
        "provider": "University of Western Australia doctoral thesis",
        "url": (
            "https://api.research-repository.uwa.edu.au/ws/portalfiles/portal/"
            "9786652/THESIS_DOCTOR_OF_PHILOSOPHY_JOHNSON_Bridget_Anne_2015.pdf"
        ),
        "title": (
            "Plant-pollinator networks in a restoration planting, and effects "
            "of non-native plants and nitrogen fertilisation"
        ),
        "citation": "Johnson (2015), University of Western Australia PhD thesis, p. 32",
        "excerpt": (
            "The significant difference in the proportion of viable seeds between "
            "the bagged treatment and control shows evidence of little "
            "self-pollination. However, because the proportions of viable seeds "
            "between the bagged treatment and control are above zero suggests that "
            "the plants are capable of some self-pollination perhaps via wind."
        ),
        "raw_value": "little self-pollination; bagged viable seed above zero; perhaps wind",
        "record_id": "uwa:johnson-2015:thesis:p32:callistemon-phoeniceus-bagged",
        "lineage": "university_thesis:uwa:johnson_2015:plant_pollinator_networks",
        "lineage_method": "original_doctoral_thesis_experiment",
        "source_tier": "A",
        "source_type": "doctoral_thesis_primary_pollination_experiment",
        "domain": "api.research-repository.uwa.edu.au",
        "language": "en",
        "content_sha256": "40d0222ed578c4d992f3e331e031d67d23b548e058d10bd3c72f13e2252aa547",
        "content_sha256_basis": "downloaded_official_university_repository_pdf_bytes",
        "potential": 0,
        "role": "direct_counterexample",
    },
    {
        "species": "Triumfetta rhomboidea",
        "family": "Malvaceae",
        "trait": "autonomous_selfing_capacity",
        "value": "delayed",
        "quality": "high",
        "provider": "Annali di Botanica primary article",
        "url": (
            "https://rosa.uniroma1.it/rosa04/annali_di_botanica/article/"
            "download/13285/13654"
        ),
        "title": "Pollination ecology of Triumfetta rhomboidea (Tiliaceae)",
        "citation": (
            "Solomon Raju and Sandhya Rani (2017), Annali di Botanica 7:33-41, "
            "DOI 10.4462/annbotrm-13285"
        ),
        "excerpt": (
            "The flowers gradually close back by the evening of the same day "
            "during which the stigma curls down and contacts the anthers, the "
            "movement of which eventually ends up in autogamy. The mature buds "
            "bagged and followed for a week showed fruit and seed set suggesting "
            "that the plant is self-compatible and self-pollinating without "
            "vector-mediation."
        ),
        "raw_value": "delayed autonomous selfing; bagged buds set fruit and seed",
        "record_id": "doi:10.4462/annbotrm-13285:p35:delayed-autonomous-selfing",
        "lineage": "doi:10.4462/annbotrm-13285",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_pollination_experiment",
        "domain": "rosa.uniroma1.it",
        "language": "en",
        "content_sha256": "90558dcb62001e4674a2c40b7d4e26c07d5faf8e11d630f820e65bd2428d7408",
        "content_sha256_basis": "downloaded_official_journal_pdf_bytes",
        "potential": 0,
        "role": "direct_enrichment",
    },
    {
        "species": "Triumfetta rhomboidea",
        "family": "Malvaceae",
        "trait": "mating_system",
        "value": "mixed_mating",
        "quality": "high",
        "provider": "Annali di Botanica primary article",
        "url": (
            "https://rosa.uniroma1.it/rosa04/annali_di_botanica/article/"
            "download/13285/13654"
        ),
        "title": "Pollination ecology of Triumfetta rhomboidea (Tiliaceae)",
        "citation": (
            "Solomon Raju and Sandhya Rani (2017), Annali di Botanica 7:33-41, "
            "DOI 10.4462/annbotrm-13285"
        ),
        "excerpt": (
            "In Triumfetta rhomboidea, the floral characteristics such as showy "
            "petals, nectar and pollen, short period of anthesis schedule, "
            "self-compatibility, brief period of stigma receptivity, medium "
            "pollen/ovule ratio and delayed autonomous selfing suggest that it is "
            "facultative autogamous with the option kept open for outcrossing."
        ),
        "raw_value": "facultative autogamous with option open for outcrossing",
        "record_id": "doi:10.4462/annbotrm-13285:abstract:mixed-mating",
        "lineage": "doi:10.4462/annbotrm-13285",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_pollination_experiment",
        "domain": "rosa.uniroma1.it",
        "language": "en",
        "content_sha256": "90558dcb62001e4674a2c40b7d4e26c07d5faf8e11d630f820e65bd2428d7408",
        "content_sha256_basis": "downloaded_official_journal_pdf_bytes",
        "potential": 0,
        "role": "direct_enrichment",
    },
    {
        "species": "Triumfetta rhomboidea",
        "family": "Malvaceae",
        "trait": "floral_symmetry",
        "value": "actinomorphic",
        "quality": "high",
        "provider": "Annali di Botanica primary article",
        "url": (
            "https://rosa.uniroma1.it/rosa04/annali_di_botanica/article/"
            "download/13285/13654"
        ),
        "title": "Pollination ecology of Triumfetta rhomboidea (Tiliaceae)",
        "citation": (
            "Solomon Raju and Sandhya Rani (2017), Annali di Botanica 7:33-41, "
            "DOI 10.4462/annbotrm-13285"
        ),
        "excerpt": (
            "The flowers are small, 6 mm long and 7 mm wide, stellate or "
            "star-shaped, yellow, bisexual and actinomorphic."
        ),
        "raw_value": "actinomorphic",
        "record_id": "doi:10.4462/annbotrm-13285:p35:actinomorphic",
        "lineage": "doi:10.4462/annbotrm-13285",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_flower_morphology",
        "domain": "rosa.uniroma1.it",
        "language": "en",
        "content_sha256": "90558dcb62001e4674a2c40b7d4e26c07d5faf8e11d630f820e65bd2428d7408",
        "content_sha256_basis": "downloaded_official_journal_pdf_bytes",
        "potential": 0,
        "role": "direct_enrichment",
    },
    {
        "species": "Sporobolus alterniflorus",
        "family": "Poaceae",
        "trait": "mating_system",
        "value": "predominantly_outcrossing",
        "quality": "high",
        "provider": "LSU institutional primary-article record",
        "url": "https://repository.lsu.edu/plantsoil_pubs/554/",
        "title": (
            "Mode of pollination, pollen germination, and seed set in smooth "
            "cordgrass (Spartina alterniflora, Poaceae)"
        ),
        "citation": (
            "Fang et al. (2004), International Journal of Plant Sciences "
            "165:395-401, DOI 10.1086/382810"
        ),
        "excerpt": (
            "Average seed set under self- and cross-pollination was 26% and 52%, "
            "respectively. These findings indicate S. alterniflora is largely "
            "cross-pollinated and that protogyny can be exploited to make "
            "controlled hybrids for genetic research and improvement."
        ),
        "raw_value": "self seed set 26%; cross seed set 52%; largely cross-pollinated",
        "record_id": "doi:10.1086/382810:abstract:spartina-alterniflora",
        "lineage": "doi:10.1086/382810",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_pollination_experiment",
        "domain": "repository.lsu.edu",
        "language": "en",
        "content_sha256": "4c0f41836c9c2f4108552a286e24d803ecbdee8fd75511e069194accfb2a2957",
        "content_sha256_basis": "downloaded_official_university_repository_html_bytes",
        "potential": 0,
        "role": "direct_counterexample",
        "matched_page_name": "Spartina alterniflora",
        "name_match_method": "synonym_exact",
        "name_resolution_lineage": "powo_strict_synonym_to_master_accepted",
    },
    {
        "species": "Melaleuca quinquenervia",
        "family": "Myrtaceae",
        "trait": "self_incompatibility",
        "value": "SC",
        "quality": "medium",
        "provider": "USDA Forest Service FEIS",
        "url": "https://research.fs.usda.gov/feis/species-reviews/melqui",
        "title": "Melaleuca quinquenervia, melaleuca",
        "citation": (
            "Munger (2005), Fire Effects Information System species review, "
            "citing Vardaman (1994)"
        ),
        "excerpt": (
            "Field studies in southern Florida showed that pollination within the "
            "same flower resulted in significantly reduced fruit set compared with "
            "pollination between flowers on different trees or pollination between "
            "flowers on the same tree. Melaleuca flowers are perfect. Vardaman "
            "concluded that melaleuca has a mixed breeding system that promotes "
            "outcrossing but allows inbreeding if sufficient outcrossing does not occur."
        ),
        "raw_value": "same-flower pollination sets fewer fruit; inbreeding allowed",
        "record_id": "usfs-feis:melqui:pollination-breeding-system:self-compatible",
        "lineage": "official_species_review:usfs_feis:melqui:munger_2005",
        "lineage_method": "official_species_review_citation_lineage",
        "source_tier": "A",
        "source_type": "official_government_species_synthesis",
        "domain": "research.fs.usda.gov",
        "language": "en",
        "content_sha256": "79d81dc31ba076954879972349f0d807b2e5a67c5cf37357c95a61f8612e5c38",
        "content_sha256_basis": "downloaded_official_government_html_bytes",
        "potential": 0,
        "role": "direct_counterexample",
    },
    {
        "species": "Melaleuca quinquenervia",
        "family": "Myrtaceae",
        "trait": "mating_system",
        "value": "mixed_mating",
        "quality": "medium",
        "provider": "USDA Forest Service FEIS",
        "url": "https://research.fs.usda.gov/feis/species-reviews/melqui",
        "title": "Melaleuca quinquenervia, melaleuca",
        "citation": (
            "Munger (2005), Fire Effects Information System species review, "
            "citing Vardaman (1994)"
        ),
        "excerpt": (
            "Melaleuca flowers are perfect. Vardaman concluded that melaleuca has "
            "a mixed breeding system that promotes outcrossing but allows "
            "inbreeding if sufficient outcrossing does not occur."
        ),
        "raw_value": "mixed breeding system; promotes outcrossing but allows inbreeding",
        "record_id": "usfs-feis:melqui:pollination-breeding-system:mixed-mating",
        "lineage": "official_species_review:usfs_feis:melqui:munger_2005",
        "lineage_method": "official_species_review_citation_lineage",
        "source_tier": "A",
        "source_type": "official_government_species_synthesis",
        "domain": "research.fs.usda.gov",
        "language": "en",
        "content_sha256": "79d81dc31ba076954879972349f0d807b2e5a67c5cf37357c95a61f8612e5c38",
        "content_sha256_basis": "downloaded_official_government_html_bytes",
        "potential": 0,
        "role": "direct_counterexample",
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
                "query": "latest_support2_exact_trait_wave9",
                "language": str(source["language"]),
                "matched_page_name": str(
                    source.get("matched_page_name", source["species"])
                ),
                "evidence_scope": "species_direct",
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
        "Accepted after exact master identity or strict synonym reconciliation, "
        "exact trait, value polarity, source text, lineage, provenance, "
        "counterevidence and cultivar scope audit."
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
        "evidence": output_dir / "high_yield_mixed_wave9_evidence_20260814.csv",
        "audit": output_dir / "high_yield_mixed_wave9_manual_audit_20260814.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260814.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260814.csv",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    role_counts = pd.Series(
        [str(source.get("role", "rule_candidate")) for source in SOURCES]
    ).value_counts()
    summary: dict[str, object] = {
        "contract": "high_yield_trait_specific_support2_checkpoint_v1",
        "created_at_utc": CREATED_AT,
        "accepted_rows": len(evidence),
        "species_trait": len(evidence),
        "quality_counts": evidence["evidence_quality"].value_counts().to_dict(),
        "trait_counts": evidence["trait_name"].value_counts().to_dict(),
        "language_counts": evidence["language"].value_counts().to_dict(),
        "source_lineages": int(evidence["source_lineage"].nunique()),
        "role_counts": role_counts.to_dict(),
        "rule_candidates": int(role_counts.get("rule_candidate", 0)),
        "direct_counterevidence_rows": int(
            role_counts.get("direct_counterexample", 0)
        ),
        "direct_enrichment_rows": int(role_counts.get("direct_enrichment", 0)),
        "pending_rejected_rows": 4,
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
            "counterevidence_preserved": True,
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
    manifest = output_dir / "high_yield_mixed_wave9_manifest_20260814.json"
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
