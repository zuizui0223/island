"""Freeze the eighth high-yield species-direct evidence checkpoint.

The rows target current support=2 genus x trait_name rules.  Queue potential is
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

CREATED_AT = "2026-08-13T22:32:48Z"
SOURCE_GROUP = "high_yield_mixed_wave8_checkpoint_20260814"

SOURCES: tuple[dict[str, object], ...] = (
    {
        "species": "Aegilops tauschii",
        "family": "Poaceae",
        "trait": "mating_system",
        "value": "predominantly_selfing",
        "quality": "high",
        "provider": "Evolution Letters primary article",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11637685/",
        "title": (
            "Mating systems and recombination landscape strongly shape genetic "
            "diversity and selection in wheat relatives"
        ),
        "citation": (
            "Burgarella et al. (2024), Evolution Letters 8:866-880, "
            "DOI 10.1093/evlett/qrae039"
        ),
        "excerpt": (
            "Ae. tauschii and Triticum sp. are known to be prevalently "
            "self-fertilizing species (Highly Selfing), while the remaining "
            "species are self-compatible (SC)."
        ),
        "raw_value": "prevalently self-fertilizing; highly selfing",
        "record_id": "doi:10.1093/evlett/qrae039:figure-1:aegilops-tauschii",
        "lineage": "doi:10.1093/evlett/qrae039",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_population_genomics",
        "domain": "pmc.ncbi.nlm.nih.gov",
        "language": "en",
        "content_sha256": "b543741757bf83b738bedb6507a19f6a9171300349501667d17ecc817e3f7f01",
        "content_sha256_basis": "downloaded_europe_pmc_full_text_xml_bytes",
        "potential": 13,
    },
    {
        "species": "Pimenta pseudocaryophyllus",
        "family": "Myrtaceae",
        "trait": "flower_primary_color",
        "value": "white",
        "quality": "medium",
        "provider": "Brazilian Ministry of Environment and Embrapa",
        "url": "https://www.alice.cnptia.embrapa.br/alice/bitstream/doc/906699/1/6524.pdf",
        "title": "Pimenta pseudocaryophyllus - Craveiro-do-mato",
        "citation": (
            "Ruschel (2011), in Espécies nativas da flora brasileira de valor "
            "econômico atual ou potencial: plantas para o futuro - Região Sul, "
            "MMA, p. 223"
        ),
        "excerpt": (
            "O craveiro-do-mato, Pimenta pseudocaryophyllus (Gomes) Landrum, "
            "trata-se de uma espécie arbórea aromática de 4-10m de altura. "
            "Inflorescências axilares em dicásios simples ou composta com duas "
            "a três flores brancas muito perfumadas."
        ),
        "raw_value": "flores brancas",
        "record_id": "mma-embrapa:plantas-futuro-sul:p223:white-flowers",
        "lineage": (
            "official_species_chapter:mma_embrapa:ruschel_2011:"
            "pimenta_pseudocaryophyllus"
        ),
        "lineage_method": "official_species_chapter_citation_lineage",
        "source_tier": "A",
        "source_type": "official_government_species_synthesis",
        "domain": "alice.cnptia.embrapa.br",
        "language": "pt",
        "content_sha256": "1c1e3290cfe98deeb8c90fd17f5c6c7ca1e853775c60597e006ea904815fdfe2",
        "content_sha256_basis": "downloaded_official_pdf_bytes",
        "potential": 15,
    },
    {
        "species": "Scolopia braunii",
        "family": "Salicaceae",
        "trait": "flower_size_class",
        "value": "small",
        "quality": "medium",
        "provider": "Australian Plants Society NSW",
        "url": (
            "https://resources.austplants.com.au/plant/"
            "scolopia-brauniiflintwood-mountain-cherry-brown-birch-scolopia/"
        ),
        "title": "Scolopia braunii - Flintwood",
        "citation": "Clarke and Howes (2023), Australian Plants Society NSW",
        "excerpt": (
            "In this species, creamy white flowers form on axillary racemes to "
            "about 4 cm long by 2 cm wide, each flower about 4 mm wide, with 4 "
            "petals and sepals."
        ),
        "raw_value": "each flower about 4 mm wide",
        "record_id": "aps-nsw:scolopia-braunii:2023-11-23:flower-4mm",
        "lineage": "curated_regional_page:aps_nsw:scolopia_braunii:2023-11-23",
        "lineage_method": "curated_regional_species_page",
        "source_tier": "B",
        "source_type": "specialist_regional_species_account",
        "domain": "resources.austplants.com.au",
        "language": "en",
        "content_sha256": "3643dea4bd91047edd8d8017c97d3a251197afd748fba42ab7917656aa8161aa",
        "content_sha256_basis": "downloaded_source_html_bytes",
        "potential": 19,
    },
    {
        "species": "Scleromitrion verticillatum",
        "family": "Rubiaceae",
        "trait": "flower_primary_color",
        "value": "white",
        "quality": "medium",
        "provider": "National Taiwan University Herbarium Tai2",
        "url": "https://tai2.ntu.edu.tw/species/516%20055%2006%200",
        "title": "Scleromitrion verticillatum - Tai2",
        "citation": "National Taiwan University Herbarium exact species treatment",
        "excerpt": (
            "Flowers several to many in axillary fascicles, rarely solitary; "
            "bracts and bracteoles linear-lanceolate, 1.5-2 mm long, hirsute; "
            "calyx ca. 4 mm long, hirsute; corolla white, sparsely hairy on upper "
            "inside."
        ),
        "raw_value": "corolla white",
        "record_id": "tai2:516-055-06-200:corolla-white",
        "lineage": "institutional_flora:tai2:516-055-06-200",
        "lineage_method": "academic_herbarium_species_treatment",
        "source_tier": "A",
        "source_type": "academic_herbarium_species_treatment",
        "domain": "tai2.ntu.edu.tw",
        "language": "en",
        "content_sha256": "1e6bf776cd7b4ebaa0e90b4d6c12f9faf9706dac73170a999c734334cbc6b9a1",
        "content_sha256_basis": "downloaded_institutional_html_bytes",
        "potential": 13,
    },
    {
        "species": "Lyonia lucida",
        "family": "Ericaceae",
        "trait": "self_incompatibility",
        "value": "SI",
        "quality": "high",
        "provider": "American Midland Naturalist primary article",
        "url": (
            "https://johnbenning.net/wp-content/uploads/2019/04/"
            "benning-2015-odd-for-an-ericad-nocturnal-pollination-of-"
            "lyonia-lucida-ericaceae.pdf"
        ),
        "title": "Odd for an Ericad: Nocturnal Pollination of Lyonia lucida",
        "citation": (
            "Benning (2015), American Midland Naturalist 174:204-217, "
            "DOI 10.1674/0003-0031-174.2.204"
        ),
        "excerpt": (
            "Flowers excluded from all visitors had an extremely low probability "
            "of fruit set (0.015), as did those flowers hand-pollinated with self "
            "pollen (0.026), suggesting that fetterbush is highly self-sterile and "
            "largely dependent on outcrossing for fruit set. Flowers with "
            "supplemental cross pollen showed drastically higher probability of "
            "fruit set than any other treatment (0.826)."
        ),
        "raw_value": "self 0.026 versus supplemental cross 0.826 fruit-set probability",
        "record_id": "doi:10.1674/0003-0031-174.2.204:fruit-set-experiment",
        "lineage": "doi:10.1674/0003-0031-174.2.204",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_pollination_experiment",
        "domain": "johnbenning.net",
        "language": "en",
        "content_sha256": "0c4b748638eb3cdabbf66f978206bad1984d3f416fac6e23c64f300350732692",
        "content_sha256_basis": "downloaded_author_hosted_publisher_pdf_bytes",
        "potential": 0,
        "role": "direct_counterexample",
    },
    {
        "species": "Anthoxanthum odoratum",
        "family": "Poaceae",
        "trait": "self_incompatibility",
        "value": "SI",
        "quality": "medium",
        "provider": "NatureServe Explorer",
        "url": "https://explorer.natureserve.org/api/data/taxon/ELEMENT_GLOBAL.2.160120",
        "title": "Anthoxanthum odoratum - NatureServe Explorer",
        "citation": "NatureServe Explorer species account, citing Antonovics (1972)",
        "excerpt": "The plants are generally self-incompatible (Antonovics 1972).",
        "raw_value": "generally self-incompatible",
        "record_id": "natureserve:element_global.2.160120:self-incompatible",
        "lineage": "natureserve:element_global.2.160120",
        "lineage_method": "official_species_account_citation_lineage",
        "source_tier": "A",
        "source_type": "official_conservation_species_synthesis",
        "domain": "explorer.natureserve.org",
        "language": "en",
        "content_sha256": "97f33f70d3648b27d701e3460a0eccecb7887f9cfeff92a91cc7a4ac4a2dd55e",
        "content_sha256_basis": "downloaded_natureserve_api_json_bytes",
        "potential": 0,
        "role": "direct_counterexample",
    },
    {
        "species": "Oxytropis campestris",
        "family": "Fabaceae",
        "trait": "self_incompatibility",
        "value": "SC",
        "quality": "medium",
        "provider": "Fabaceae experimental book chapter",
        "url": (
            "https://cannalib.eu/wp-content/uploads/2025/03/"
            "Fabaceae-Classification-nutrient-composition-and-health-benefits-2015.pdf"
        ),
        "title": "Pollination Ecology of the Bulgarian Legumes",
        "citation": (
            "Kozuharova et al. (2015), chapter 3 in Fabaceae: Classification, "
            "Nutrient Composition and Health Benefits, ISBN 978-1-63482-200-8"
        ),
        "excerpt": (
            "The results from our field experiments and P/O tests indicate that "
            "all these species are potentially self-fertile and do not have "
            "self-incompatibility. Self-pollination occurs and seed is set. The "
            "prevention of self-pollination is by means of spacial isolation "
            "between stigma and anthers of the flower. The process is more "
            "expressed in Oxytropis urumovii and less expressed in O. campestris."
        ),
        "raw_value": "self-fertile; no self-incompatibility; self-pollination sets seed",
        "record_id": "isbn:978-1-63482-200-8:chapter3:p92:oxytropis-campestris",
        "lineage": "isbn:978-1-63482-200-8:chapter3",
        "lineage_method": "primary_experimental_book_chapter",
        "source_tier": "A",
        "source_type": "primary_experimental_book_chapter",
        "domain": "cannalib.eu",
        "language": "en",
        "content_sha256": "3e7961dea0c9b432440ce908d97d4b01463345d7c5a4be7ac7a8347e447bfbd8",
        "content_sha256_basis": "downloaded_complete_book_pdf_bytes",
        "potential": 0,
        "role": "direct_conflict",
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
                "query": "latest_support2_third_species_exact_trait_wave8",
                "language": str(source["language"]),
                "matched_page_name": str(source["species"]),
                "evidence_scope": "species_direct",
                "name_match_method": "accepted_name_exact",
                "name_resolution_lineage": "master_accepted_name_exact",
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
        "Accepted after exact master identity, exact trait, value polarity, source "
        "text, lineage, provenance, counterevidence and cultivar scope audit."
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
        "evidence": output_dir / "high_yield_mixed_wave8_evidence_20260814.csv",
        "audit": output_dir / "high_yield_mixed_wave8_manual_audit_20260814.csv",
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
        "rule_candidates": sum(
            1 for source in SOURCES if source.get("role", "rule_candidate") == "rule_candidate"
        ),
        "direct_counterevidence_rows": sum(
            1 for source in SOURCES if source.get("role", "rule_candidate") != "rule_candidate"
        ),
        "pending_rejected_rows": 2,
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
    manifest = output_dir / "high_yield_mixed_wave8_manifest_20260814.json"
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
