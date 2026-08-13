"""Freeze a reviewed reproductive-evidence and counterevidence checkpoint.

The batch deliberately imports positive candidates together with direct
counterexamples.  It emits species-direct records only; the shared all-evidence
builder remains responsible for conflict resolution and genus x trait_name Low
rules.  In particular, compatibility, autonomous selfing, mating system and
cleistogamy are never substituted for one another.
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

CREATED_AT = "2026-08-13T13:10:00Z"
SOURCE_GROUP = "reproductive_counterevidence_checkpoint_20260813"


ROWS = [
    {
        "species": "Abutilon theophrasti",
        "trait": "autonomous_selfing_capacity",
        "value": "autonomous",
        "quality": "medium",
        "provider": "CONABIO official invasive-species assessment",
        "url": (
            "https://sivicoff.cnf.gob.mx/ContenidoPublico/MenuPrincipal/"
            "07Fichas%20tecnicas_OK/02Fichas%20tecnicas/"
            "Fichas%20t%C3%A9cnicas%20CONABIO_especies%20ex%C3%B3ticas/"
            "Fichas%20plantas%20invasoras/A_B/Abutilon%20theophrasti.pdf"
        ),
        "title": "Ponderacion de invasividad de la especie Abutilon theophrasti",
        "citation": (
            "CONABIO (2014), citing Andersen (1988), Outcrossing in Velvetleaf, "
            "DOI 10.1017/S0043174500075470"
        ),
        "excerpt": (
            "A. theophrasti is a self-compatible, autogamous species. [...] "
            "pollination would have occurred before stigmas could be exposed "
            "to pollen from another flower."
        ),
        "record_id": "conabio:abutilon-theophrasti:autogamy",
        "lineage": "citation:andersen_1988_abutilon_theophrasti",
        "lineage_method": "underlying_primary_citation_lineage",
        "source_tier": "A",
        "source_type": "official_agency_species_assessment",
        "domain": "sivicoff.cnf.gob.mx",
        "content_sha256": ("2f420f87db066e4b52638e71620408a48144dcc4475e0bcbee66b86f3db9e406"),
        "content_sha256_basis": "downloaded_official_pdf_bytes",
        "raw_value": "self-compatible, autogamous; prior self-pollination",
    },
    {
        "species": "Allium cernuum",
        "trait": "autonomous_selfing_capacity",
        "value": "autonomous",
        "quality": "high",
        "provider": "Ksiazek et al. primary pollinator-exclusion experiment",
        "url": (
            "https://greenroofresearch.wordpress.com/wp-content/uploads/2016/02/ksiazekal2012.pdf"
        ),
        "title": "An assessment of pollen limitation on Chicago green roofs",
        "citation": (
            "Ksiazek, Fant & Skogen (2012), Landscape and Urban Planning 107:"
            "401-408, DOI 10.1016/j.landurbplan.2012.07.008"
        ),
        "excerpt": (
            "To determine the rate of spontaneous autogamy (self-fertilization), "
            "[...] a pollinator exclusion bag [...] was placed over a single "
            "flower bud or buds on a single inflorescence per plant. [...] "
            "Allium cernuum [...] Seed set attributed to autogamy (%) 11"
        ),
        "record_id": "doi:10.1016/j.landurbplan.2012.07.008:table2:allium-cernuum",
        "lineage": "doi:10.1016/j.landurbplan.2012.07.008",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_pollinator_exclusion_experiment",
        "domain": "greenroofresearch.wordpress.com",
        "content_sha256": ("04595eb890c3e3d9ec96fc9b5d9f3d2a739d6239b2977f73067059043d3d229e"),
        "content_sha256_basis": "downloaded_author_hosted_pdf_bytes",
        "raw_value": "pollinator exclusion; 11% seed set attributed to autogamy",
    },
    {
        "species": "Allium oleraceum",
        "trait": "autonomous_selfing_capacity",
        "value": "absent",
        "quality": "high",
        "provider": "Annales Botanici Fennici primary bagging experiment",
        "url": "https://www.sekj.org/PDF/anb41-free/anb41-001.pdf",
        "title": "Generative reproduction in Allium oleraceum (Alliaceae)",
        "citation": "Astrom & Haeggstrom (2004), Annales Botanici Fennici 41:1-14",
        "excerpt": (
            "The bagging experiment gave the result that not a single flower "
            "developed a capsule in the bagged plants in 1998. [...] In the "
            "bagged plants, not a single flower had developed capsules or could "
            "be seen to have begun development into fruits (straightening of "
            "the pedicel). Our interpretation is that without insect visitors, "
            "capsules and seeds do not develop."
        ),
        "record_id": "astrom-haeggstrom-2004:allium-oleraceum:bagging",
        "lineage": "publication:astrom_haeggstrom_2004_allium_oleraceum",
        "lineage_method": "original_primary_article_lineage",
        "source_tier": "A",
        "source_type": "peer_reviewed_pollinator_exclusion_experiment",
        "domain": "sekj.org",
        "content_sha256": ("e899d75b702b24c71e52c6585f91a4737bf3ccd3c643ea0e689221ed6adc182a"),
        "content_sha256_basis": "verified_exact_supporting_excerpt_utf8_bytes",
        "raw_value": "bagged flowers produced no capsules or seeds",
    },
    {
        "species": "Calanthe mannii",
        "trait": "autonomous_selfing_capacity",
        "value": "autonomous",
        "quality": "medium",
        "provider": "Missouri Botanical Garden historical orchid monograph",
        "url": "https://archive.org/details/mobot31753002133210",
        "title": "The Orchids of the Sikkim-Himalaya",
        "citation": "King & Pantling (1898), The Orchids of the Sikkim-Himalaya, p. 167",
        "excerpt": (
            "This is a self-fertile species. Mr. Pantling was unable to obtain "
            "pollinia for figuring, as in all the flowers of the only living "
            "specimen found by him, they had been more or less absorbed by the "
            "stigmas, the clinandrium having apparently become absorbed."
        ),
        "record_id": "king-pantling-1898:calanthe-mannii:p167-autogamy",
        "lineage": "monograph:king_pantling_1898_orchids_sikkim:calanthe_mannii",
        "lineage_method": "original_monograph_treatment_lineage",
        "source_tier": "A",
        "source_type": "primary_species_monograph_observation",
        "domain": "archive.org",
        "content_sha256": ("1e27f565dd9770f5ca293916dc5a9bec7311280bf2f259be0564542ef8a602c7"),
        "content_sha256_basis": "downloaded_mobot_scan_ocr_text_bytes",
        "raw_value": "self-fertile with spontaneous pollinium-stigma contact",
    },
    {
        "species": "Calanthe striata",
        "trait": "autonomous_selfing_capacity",
        "value": "absent",
        "quality": "high",
        "provider": "Horticulturae controlled-pollination experiment",
        "url": "https://doi.org/10.3390/horticulturae10101025",
        "title": "Flowering Phenology and Mating System of Calanthe sieboldii",
        "citation": "Zhang et al. (2024), Horticulturae 10:1025",
        "excerpt": (
            "the rates of ovary enlargement in the spontaneous autogamy test "
            "(flowers bagged without emasculation) and agamospermy test "
            "(flowers emasculated and bagged) for C. sieboldii were 0% and "
            "0.83%, respectively (n = 120), suggesting the absence of "
            "spontaneous autogamy and the possible occurrence of agamospermy."
        ),
        "record_id": "doi:10.3390/horticulturae10101025:table:spontaneous-autogamy",
        "lineage": "doi:10.3390/horticulturae10101025",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_controlled_pollination_experiment",
        "domain": "doi.org",
        "content_sha256": ("47a931a5e99ee5c8172b87214fc66e314fbb86a79e3e722b238573c8b6c711c4"),
        "content_sha256_basis": "downloaded_publisher_pdf_bytes",
        "raw_value": "bagged without emasculation 0% ovary enlargement, n=120",
        "matched_page_name": "Calanthe sieboldii",
        "name_match_method": "exact_synonym",
        "name_resolution_lineage": (
            "powo:wcvp:urn:lsid:ipni.org:names:621044-1:synonym_of_Calanthe_striata"
        ),
    },
    {
        "species": "Hakea carinata",
        "trait": "mating_system",
        "value": "predominantly_selfing",
        "quality": "high",
        "provider": "University of Adelaide doctoral thesis",
        "url": (
            "https://digital.library.adelaide.edu.au/server/api/core/bitstreams/"
            "7e6f9e21-e413-48c2-8655-9396a292d389/content"
        ),
        "title": "Population genetics of Hakea carinata",
        "citation": (
            "Starr (2001), University of Adelaide PhD thesis; chapter 2 reports "
            "the Starr & Carthew (1998) progeny-array lineage, DOI 10.1071/BT97123"
        ),
        "excerpt": (
            "It survives in very small populations by using a predominantly "
            "selfing mating system [...] mean estimated outcrossing rate over "
            "30 populations, t = 0.111."
        ),
        "record_id": "starr-2001:hakea-carinata:chapter2:mating-system",
        "lineage": "doi:10.1071/BT97123",
        "lineage_method": "underlying_primary_study_doi",
        "source_tier": "A",
        "source_type": "doctoral_thesis_progeny_array_analysis",
        "domain": "digital.library.adelaide.edu.au",
        "content_sha256": ("7f81f49a203ea024d7c399bd56ce80ae548e2abcfc95b7a10a0f823111a9bd02"),
        "content_sha256_basis": "downloaded_official_thesis_bitstream_bytes",
        "raw_value": "predominantly selfing; mean multilocus t=0.111",
    },
    {
        "species": "Liparis liliifolia",
        "trait": "autonomous_selfing_capacity",
        "value": "absent",
        "quality": "medium",
        "provider": "Environment and Climate Change Canada SARA",
        "url": (
            "https://www.canada.ca/en/environment-climate-change/services/"
            "species-risk-public-registry/recovery-strategies/"
            "purple-twayblade-2018.html"
        ),
        "title": "Purple Twayblade recovery strategy",
        "citation": (
            "Environment and Climate Change Canada recovery strategy, citing Whigham et al. (2002)"
        ),
        "excerpt": (
            "Purple Twayblade is incapable of self-fertilization, meaning that "
            "it requires cross-pollination to produce viable seed."
        ),
        "record_id": "canada-sara:liparis-liliifolia:pollinators",
        "lineage": "citation:whigham_et_al_2002:liparis_liliifolia",
        "lineage_method": "underlying_primary_citation_lineage",
        "source_tier": "A",
        "source_type": "official_species_recovery_strategy",
        "domain": "canada.ca",
        "content_sha256": ("8c3018d521fd2658cd5a38a69b2c435b17111e91701c3bd5ff32d3e40aa69a2e"),
        "content_sha256_basis": "verified_exact_supporting_excerpt_utf8_bytes",
        "raw_value": "incapable of self-fertilization; requires cross-pollination",
    },
    {
        "species": "Liparis loeselii",
        "trait": "autonomous_selfing_capacity",
        "value": "autonomous",
        "quality": "medium",
        "provider": "Journal of Ecology Biological Flora synthesis",
        "url": "https://doi.org/10.1111/1365-2745.14086",
        "title": "Biological Flora of Britain and Ireland: Liparis loeselii",
        "citation": (
            "Jacquemyn et al. (2023), Journal of Ecology, "
            "DOI 10.1111/1365-2745.14086, citing Catling (1980), "
            "DOI 10.2307/2484083"
        ),
        "excerpt": (
            "Liparis loeselii is predominantly autogamous, and fruit production "
            "is often the result of rain-facilitated self-pollination."
        ),
        "record_id": "doi:10.1111/1365-2745.14086:liparis-loeselii:autogamy",
        "lineage": "citation:catling_1980_liparis_loeselii_autogamy",
        "lineage_method": "underlying_primary_citation_lineage",
        "source_tier": "A",
        "source_type": "peer_reviewed_species_biological_flora_synthesis",
        "domain": "doi.org",
        "content_sha256": ("fea3f8c296a4203816b76075d9fdc4d41d8506d380bfa6f40ee204f52fc36a30"),
        "content_sha256_basis": "verified_exact_supporting_excerpt_utf8_bytes",
        "raw_value": "predominantly autogamous; rain-facilitated self-pollination",
    },
    {
        "species": "Grevillea robusta",
        "trait": "mating_system",
        "value": "predominantly_outcrossing",
        "quality": "high",
        "provider": "Kalinganire et al. primary mating-system study",
        "url": "https://doi.org/10.1006/anbo.2000.1170",
        "title": "Flowering and fruiting phenology of Grevillea robusta",
        "citation": "Kalinganire et al. (2000), Annals of Botany, DOI 10.1006/anbo.2000.1170",
        "excerpt": "Allogamy is its primary breeding behaviour.",
        "record_id": "doi:10.1006/anbo.2000.1170:grevillea-robusta:allogamy",
        "lineage": "doi:10.1006/anbo.2000.1170",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_mating_system_article",
        "domain": "doi.org",
        "content_sha256": ("804df46c9a088ce9447c3bfcf9d6dc330c1214e2e2ca6d051ff81851540dd748"),
        "content_sha256_basis": "verified_exact_supporting_excerpt_utf8_bytes",
        "raw_value": "allogamy is the primary breeding behaviour",
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
        evidence["query"] = "support2_primary_source_with_counterevidence"
        for optional in (
            "matched_page_name",
            "name_match_method",
            "name_resolution_lineage",
        ):
            if optional in item:
                evidence[optional] = item[optional]
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
        raise ValueError(f"checkpoint species missing from master: {missing}")
    if evidence[["accepted_species", "trait_name"]].duplicated().any():
        raise ValueError("checkpoint species-trait pairs must be unique")
    if not evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all():
        raise ValueError("every checkpoint row requires a 64-character content hash")

    audit = _audit(evidence)
    audit["reviewer"] = "Codex source-backed reproductive counterevidence audit"
    audit["reviewed_at_utc"] = CREATED_AT
    audit["decision_reason"] = (
        "Accepted after target-master identity, exact same-trait statement, "
        "source lineage, value polarity and cultivar status were rechecked. "
        "Positive candidates and direct counterexamples are imported together; "
        "no cross-trait substitution or genus inference is performed here."
    )

    prior_evidence = pd.read_csv(prior_curated_evidence_csv, dtype=str).fillna("")
    prior_audit = pd.read_csv(prior_curated_audit_csv, dtype=str).fillna("")
    owned = prior_evidence["source_group"].eq(SOURCE_GROUP)
    prior_owned_ids = set(prior_evidence.loc[owned, "candidate_id"])
    current_ids = set(evidence["candidate_id"])
    combined_evidence = pd.concat([prior_evidence.loc[~owned], evidence], ignore_index=True)
    combined_audit = pd.concat(
        [
            prior_audit.loc[~prior_audit["candidate_id"].isin(prior_owned_ids | current_ids)],
            audit,
        ],
        ignore_index=True,
    )
    for label, frame in (("evidence", combined_evidence), ("audit", combined_audit)):
        if frame["candidate_id"].duplicated().any():
            raise ValueError(f"combined {label} candidate IDs must be unique")

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "evidence": output_dir / "reproductive_counterevidence_evidence_20260813.csv",
        "audit": output_dir / "reproductive_counterevidence_manual_audit_20260813.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260813.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260813.csv",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    summary: dict[str, object] = {
        "contract": "reproductive_counterevidence_checkpoint_v1",
        "created_at_utc": CREATED_AT,
        "accepted_rows": len(evidence),
        "species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "quality_counts": evidence["evidence_quality"].value_counts().to_dict(),
        "trait_counts": evidence["trait_name"].value_counts().to_dict(),
        "source_lineages": int(evidence["source_lineage"].nunique()),
        "positive_and_counterevidence_together": True,
        "theoretical_clean_abutilon_rule_cells": 55,
        "theoretical_only_not_formal_coverage": True,
        "blocked_nominal_rules": {
            "Allium_autonomous": "counterexample lowers dominance below 0.95",
            "Calanthe_autonomous": "corrected polarity and counterevidence",
            "Liparis_absent": "direct autonomous counterexample",
            "Grevillea_mixed_mating": "predominantly-outcrossing counterexample",
        },
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
    manifest_path = output_dir / "reproductive_counterevidence_manifest_20260813.json"
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
