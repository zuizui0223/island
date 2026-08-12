"""Freeze the second trait-specific genus-rule unlock evidence wave.

The checkpoint contains only individually reviewed species-direct statements.
Reproductive traits remain separate, and the five morphology records are direct
species descriptions.  Genus inference is rebuilt later by the shared
all-evidence implementation using ``genus x trait_name`` rules.
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

CREATED_AT = "2026-08-12T09:35:00Z"
REVIEWER = "Codex source-backed rule-unlock wave-2 audit"
SOURCE_GROUP = "rule_unlock_wave2_checkpoint_20260812"


def _text_sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _row(
    *,
    species: str,
    trait: str,
    value: str,
    raw_value: str,
    excerpt: str,
    quality: str,
    provider: str,
    url: str,
    title: str,
    citation: str,
    record_id: str,
    lineage: str,
    lineage_method: str,
    source_tier: str,
    source_type: str,
    domain: str,
    content_sha256: str,
    content_sha256_basis: str,
    language: str = "en",
    name_resolution_lineage: str = "master_accepted_name_exact",
) -> dict[str, str]:
    row = _evidence_row(
        species=species,
        trait=trait,
        value=value,
        quality=quality,
        provider=provider,
        url=url,
        title=title,
        citation=citation,
        excerpt=excerpt,
        record_id=record_id,
        lineage=lineage,
        lineage_method=lineage_method,
        source_tier=source_tier,
        source_type=source_type,
        domain=domain,
        content_sha256=content_sha256,
        content_sha256_basis=content_sha256_basis,
        retrieved_at_utc=CREATED_AT,
        raw_value=raw_value,
    )
    row["source_group"] = SOURCE_GROUP
    row["language"] = language
    row["name_resolution_lineage"] = name_resolution_lineage
    row["wild_cultivated_cultivar_status"] = (
        "wild_or_species_level_statement_not_cultivar_limited"
    )
    row["query"] = "current_support_2_rule_unlock_original_or_direct_species_source"
    return row


def reviewed_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    myrcia_excerpt = (
        "The treatment group marked and bagged pre-anthesis buds, restricting "
        "access to the flower by floral visitors, but allowing self-pollination. "
        "Table S1: Myrcia guianensis; 9 individuals; 18 branches / 843 marked "
        "flowers; D=0.68; open-pollination fruit-set=20%; self-pollination "
        "fruit-set=5%; biotic pollination dependent."
    )
    rows.append(
        _row(
            species="Myrcia guianensis",
            trait="autonomous_selfing_capacity",
            value="mixed_or_variable",
            raw_value="bagged self fruit-set 5%; open fruit-set 20%; ratio 0.25",
            excerpt=myrcia_excerpt,
            quality="high",
            provider="Martins 2025 UNESP doctoral thesis",
            url=(
                "https://repositorio.unesp.br/server/api/core/bitstreams/"
                "4d986862-3528-4b9a-8ff5-507be198c489/content"
            ),
            title=(
                "Estudo de longa duração em cerrado sensu stricto: entendendo "
                "as mudanças do clima e seus efeitos na fenologia e reprodução "
                "das plantas, serviços ecossistêmicos e a sua percepção pela sociedade"
            ),
            citation=(
                "Martins, Amanda Eburneo (2025), doctoral thesis, Universidade "
                "Estadual Paulista (UNESP), Instituto de Biociências, Rio Claro, "
                "Table S1"
            ),
            record_id="unesp:martins-2025:table-s1:myrcia-guianensis",
            lineage="doi:10.1111/1365-2435.70090",
            lineage_method="underlying_experiment_shared_by_thesis_and_article",
            source_tier="A",
            source_type="doctoral_thesis_field_bagging_experiment",
            domain="repositorio.unesp.br",
            content_sha256=(
                "96c6984dea36d3171c583fc5ae6acd1b13f35f984ef6b4727e26c0144bd0e1a6"
            ),
            content_sha256_basis="downloaded_official_repository_pdf_bytes",
        )
    )

    affine_excerpt = (
        "Melastoma affine is self-compatible but does not produce fruit via "
        "autogamy or apomixis, i.e., pollen vectors are required for fruit set."
    )
    affine_common = {
        "species": "Melastoma affine",
        "provider": "Gross 1993 Biotropica",
        "url": "https://www.jstor.org/stable/2388870",
        "title": (
            "The Breeding System and Pollinators of Melastoma affine "
            "(Melastomataceae); A Pioneer Shrub in Tropical Australia"
        ),
        "citation": (
            "Gross, C. L. (1993), Biotropica 25(4):468-474, "
            "DOI 10.2307/2388870"
        ),
        "lineage": "doi:10.2307/2388870",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_pollination_experiment",
        "domain": "jstor.org",
        "content_sha256": _text_sha256(affine_excerpt),
        "content_sha256_basis": "verified_original_abstract_excerpt_utf8_bytes",
    }
    rows.extend(
        [
            _row(
                trait="autonomous_selfing_capacity",
                value="absent",
                raw_value="does not produce fruit via autogamy or apomixis",
                excerpt=affine_excerpt,
                record_id="doi:10.2307/2388870:abstract:no-autogamy",
                quality="high",
                **affine_common,
            ),
            _row(
                trait="self_incompatibility",
                value="SC",
                raw_value="self-compatible",
                excerpt=affine_excerpt,
                record_id="doi:10.2307/2388870:abstract:self-compatible",
                quality="high",
                **affine_common,
            ),
        ]
    )

    hosta_excerpt = (
        "No fruits matured in Treatment 1 (autonomous selfing), but all flowers "
        "in Treatment 2 (hand-selfing) and Treatment 3 (hand-crossing) matured "
        "into fruits. Mean seed/ovule ratios were 0.464 and 0.450 for Treatments "
        "2 and 3 and did not differ (P=0.53). These results indicate that Hosta "
        "ventricosa is self-compatible, but pollinator visitation is necessary."
    )
    hosta_common = {
        "species": "Hosta ventricosa",
        "provider": "Cao et al. 2015 Journal of Plant Ecology",
        "url": "https://academic.oup.com/jpe/article/8/2/142/885353",
        "title": (
            "Floral sex allocation and reproductive success within "
            "inflorescences of Hosta ventricosa, a pseudogamous apomict"
        ),
        "citation": (
            "Cao et al. (2015), Journal of Plant Ecology 8(2):142-150, "
            "DOI 10.1093/jpe/rtv010"
        ),
        "lineage": "doi:10.1093/jpe/rtv010",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_controlled_pollination_experiment",
        "domain": "academic.oup.com",
        "content_sha256": (
            "eb2632a4c9f769a7a1e002b46aa08705b2966fc08b28e5f28373e9f4f1030574"
        ),
        "content_sha256_basis": "retrieved_publisher_fulltext_html_bytes",
    }
    rows.extend(
        [
            _row(
                trait="autonomous_selfing_capacity",
                value="absent",
                raw_value="no fruits under autonomous-selfing treatment",
                excerpt=hosta_excerpt,
                record_id="doi:10.1093/jpe/rtv010:results:autonomous-selfing",
                quality="high",
                **hosta_common,
            ),
            _row(
                trait="self_incompatibility",
                value="SC",
                raw_value="hand-self and hand-cross fruit and seed set equivalent",
                excerpt=hosta_excerpt,
                record_id="doi:10.1093/jpe/rtv010:results:self-compatible",
                quality="high",
                **hosta_common,
            ),
        ]
    )

    ornith_excerpt = "Ornithogalum thyrsoides was capable of autogamy and selfing."
    ornith_common = {
        "species": "Ornithogalum thyrsoides",
        "provider": "Donaldson et al. 2002 Conservation Biology",
        "url": (
            "https://conbio.onlinelibrary.wiley.com/doi/"
            "10.1046/j.1523-1739.2002.99515.x"
        ),
        "title": (
            "Effects of Habitat Fragmentation on Pollinator Diversity and Plant "
            "Reproductive Success in Renosterveld Shrublands of South Africa"
        ),
        "citation": (
            "Donaldson et al. (2002), Conservation Biology 16(5):1267-1276, "
            "DOI 10.1046/j.1523-1739.2002.99515.x"
        ),
        "lineage": "doi:10.1046/j.1523-1739.2002.99515.x",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_pollination_experiment",
        "domain": "conbio.onlinelibrary.wiley.com",
        "content_sha256": _text_sha256(ornith_excerpt),
        "content_sha256_basis": "verified_original_fulltext_excerpt_utf8_bytes",
    }
    rows.extend(
        [
            _row(
                trait="autonomous_selfing_capacity",
                value="autonomous",
                raw_value="capable of autogamy",
                excerpt=ornith_excerpt,
                record_id="doi:10.1046/j.1523-1739.2002.99515.x:autogamy",
                quality="high",
                **ornith_common,
            ),
            _row(
                trait="self_incompatibility",
                value="SC",
                raw_value="capable of selfing",
                excerpt=ornith_excerpt,
                record_id="doi:10.1046/j.1523-1739.2002.99515.x:selfing",
                quality="high",
                **ornith_common,
            ),
        ]
    )

    malabathricum_excerpt = (
        "Table 3: experimentally selfed flowers set fruit (yellow-anther pollen "
        "66.67%; purple-anther pollen 78.23%), experimentally out-crossed flowers "
        "set fruit (75.56% and 73.68%), and 52 bagged flowers had 0% fruit set. "
        "The authors state that the manipulations demonstrate absence of "
        "self-pollination and agamospermy and a facultative-xenogamy system."
    )
    malabathricum_common = {
        "species": "Melastoma malabathricum",
        "provider": "Lu et al. 2009 Biodiversity Science",
        "url": (
            "https://www.biodiversity-science.net/EN/"
            "10.3724/SP.J.1003.2009.08317"
        ),
        "title": "Division of labor of heteromorphic stamens in Melastoma malabathricum",
        "citation": (
            "Lu, Wu, Wang, Li & Wang (2009), Biodiversity Science 17(2):174-181, "
            "DOI 10.3724/SP.J.1003.2009.08317"
        ),
        "lineage": "doi:10.3724/SP.J.1003.2009.08317",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_controlled_pollination_experiment",
        "domain": "biodiversity-science.net",
        "content_sha256": (
            "403ad953fe72baeb79e35eb99ce3d8e6e23973e891962e2da1455f2e4bb22f9f"
        ),
        "content_sha256_basis": "retrieved_publisher_fulltext_html_bytes",
        "language": "zh+en",
    }
    rows.extend(
        [
            _row(
                trait="autonomous_selfing_capacity",
                value="absent",
                raw_value="52 bagged flowers, 0% fruit set",
                excerpt=malabathricum_excerpt,
                record_id="doi:10.3724/SP.J.1003.2009.08317:table3:bagged",
                quality="high",
                **malabathricum_common,
            ),
            _row(
                trait="self_incompatibility",
                value="SC",
                raw_value="selfed and out-crossed treatments had similar fruit set",
                excerpt=malabathricum_excerpt,
                record_id="doi:10.3724/SP.J.1003.2009.08317:table3:self-compatible",
                quality="high",
                **malabathricum_common,
            ),
        ]
    )

    rows.append(
        _row(
            species="Cornus sericea",
            trait="autonomous_selfing_capacity",
            value="absent",
            raw_value="bagged flowers did not produce fruits",
            excerpt=(
                "Cornus sericea, redosier dogwood. In experiments, redosier dogwood "
                "flowers that were bagged to prevent cross pollination did not produce "
                "fruits, suggesting that successful fruit production depends on cross "
                "pollination [92]."
            ),
            quality="high",
            provider="USDA Forest Service Fire Effects Information System",
            url="https://research.fs.usda.gov/feis/species-reviews/corser",
            title="Cornus sericea, redosier dogwood",
            citation=(
                "Gucker, Corey (2012), Cornus sericea, redosier dogwood, Fire "
                "Effects Information System, USDA Forest Service; controlled "
                "pollination result cited as reference 92"
            ),
            record_id="usfs-feis:cornus-sericea:pollination:experiment-92",
            lineage="usfs-feis:cornus-sericea:experiment-92",
            lineage_method="official_agency_species_review_underlying_experiment",
            source_tier="A",
            source_type="official_agency_species_review_controlled_experiment",
            domain="research.fs.usda.gov",
            content_sha256=(
                "cae0ddd72b246270053621e72307e4ee765ffb3cda95701f1d31c9084b466b19"
            ),
            content_sha256_basis="retrieved_usda_forest_service_page_html_bytes",
        )
    )

    daucus_excerpt = (
        "Daucus montanus. In all enclosed inflorescences of both wild species, "
        "normal seeds were formed in 90% to 100% of the individual flowers, an "
        "indication that they are self-compatible. Because pollinating insects "
        "were excluded, the authors postulated that D. montanus is autogamous as well."
    )
    daucus_common = {
        "species": "Daucus montanus",
        "provider": "Ibañez and Camadro 2015 Botany",
        "url": (
            "https://ri.conicet.gov.ar/bitstream/handle/11336/101699/"
            "CONICET_Digital_Nro.1d9c31b9-4e1f-43ee-8b63-5314c600dc71_A.pdf"
            "?sequence=2"
        ),
        "title": (
            "Reproductive behavior of the wild carrots Daucus pusillus and "
            "Daucus montanus from Argentina"
        ),
        "citation": (
            "Ibañez, M.S. & Camadro, E.L. (2015), Botany 93:279-286, "
            "DOI 10.1139/cjb-2014-0243"
        ),
        "lineage": "doi:10.1139/cjb-2014-0243",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_bagging_experiment",
        "domain": "ri.conicet.gov.ar",
        "content_sha256": _text_sha256(daucus_excerpt),
        "content_sha256_basis": "verified_original_pdf_fulltext_excerpt_utf8_bytes",
    }
    rows.extend(
        [
            _row(
                trait="self_incompatibility",
                value="SC",
                raw_value="90%-100% seeded flowers in enclosed inflorescences",
                excerpt=daucus_excerpt,
                record_id="doi:10.1139/cjb-2014-0243:results:self-compatible",
                quality="high",
                **daucus_common,
            ),
            _row(
                trait="autonomous_selfing_capacity",
                value="autonomous",
                raw_value="bagged inflorescences formed seeds without pollinating insects",
                excerpt=daucus_excerpt,
                record_id="doi:10.1139/cjb-2014-0243:discussion:autogamous",
                quality="high",
                **daucus_common,
            ),
        ]
    )

    rows.append(
        _row(
            species="Drypetes assamica",
            trait="floral_symmetry",
            value="actinomorphic",
            raw_value="floral symmetry: actinomorphic",
            excerpt=(
                "Tree Species Details: Genus Drypetes; Species Drypetes assamica; "
                "floral symmetry: actinomorphic."
            ),
            quality="medium",
            provider="Bangladesh Forest Information System",
            url="https://bfis.bforest.gov.bd/nef/index.php/data/dataSpecies/40",
            title="Forest Emission Factor Database - Drypetes assamica",
            citation=(
                "Bangladesh Forest Information System, Forest Emission Factor "
                "Database species record: Drypetes assamica"
            ),
            record_id="bfis:forest-emission-factor:drypetes-assamica:floral-symmetry",
            lineage="url:https://bfis.bforest.gov.bd/nef/index.php/data/dataSpecies/40",
            lineage_method="canonical_government_database_species_page",
            source_tier="A",
            source_type="government_forest_species_database",
            domain="bfis.bforest.gov.bd",
            content_sha256=(
                "46b78cef8b2d77a9df7b02846c46146a43152f531adb9dac8b8fc80e1c6d2079"
            ),
            content_sha256_basis="retrieved_government_database_html_bytes",
            name_resolution_lineage=(
                "master_accepted_name_exact; historical Euphorbiaceae placement "
                "reconciled to current Putranjivaceae"
            ),
        )
    )

    rows.append(
        _row(
            species="Melastoma malabathricum",
            trait="floral_symmetry",
            value="actinomorphic",
            raw_value="flores ... actinomorfas",
            excerpt=(
                "Melastoma malabathricum L. Las flores se agrupan en inflorescencias "
                "cimosas y presentan características hermafroditas, actinomorfas y "
                "pentámeras."
            ),
            quality="medium",
            provider="El Rincón del Botánico",
            url=(
                "https://www.elrincondelbotanico.es/flora/melastomataceae/"
                "melastoma-malabathricum/"
            ),
            title="Melastoma malabathricum - El Rincón del Botánico",
            citation="El Rincón del Botánico species account: Melastoma malabathricum L.",
            record_id="elrincon:melastoma-malabathricum:floral-symmetry",
            lineage=(
                "url:https://www.elrincondelbotanico.es/flora/melastomataceae/"
                "melastoma-malabathricum"
            ),
            lineage_method="canonical_specialist_botanical_species_page",
            source_tier="C",
            source_type="specialist_personal_botanical_species_description",
            domain="elrincondelbotanico.es",
            content_sha256=(
                "4c5171a48fa6ad86991bd777b4aa468d64a13b0600c41ba44a60088b272c0920"
            ),
            content_sha256_basis="retrieved_source_page_html_bytes",
            language="es",
        )
    )

    rows.append(
        _row(
            species="Palaquium obovatum",
            trait="tube_depth_class",
            value="shallow",
            raw_value="corolla tube very short",
            excerpt=(
                "Palaquium obovatum, King & Gamble. Corolla campanulate; tube "
                "very short; lobes 6, imbricate and twisted."
            ),
            quality="high",
            provider="King and Gamble 1906 Flora of the Malayan Peninsula",
            url=(
                "https://archive.org/download/JournalAsiaticSv74pAsia/"
                "JournalAsiaticSv74pAsia.pdf"
            ),
            title="Materials for a Flora of the Malayan Peninsula, No. 17",
            citation=(
                "King, G. & Gamble, J.S. (1906), Journal of the Asiatic Society "
                "of Bengal, Part II, Natural Science 74(2):190-191"
            ),
            record_id=(
                "archive:JournalAsiaticSv74pAsia:palaquium-obovatum:corolla-tube"
            ),
            lineage="publication:king-gamble-1906-flora-malayan-peninsula",
            lineage_method="original_taxonomic_treatment_publication",
            source_tier="A",
            source_type="primary_taxonomic_monograph_species_description",
            domain="archive.org",
            content_sha256=(
                "6e81c50dee3798feacf92d45f0b4183ec6e56dba8d0a68dc29c24fcddbb38f66"
            ),
            content_sha256_basis="downloaded_internet_archive_original_scan_pdf_bytes",
        )
    )

    rows.append(
        _row(
            species="Alangium salviifolium",
            trait="floral_symmetry",
            value="actinomorphic",
            raw_value="floral symmetry: actinomorphic",
            excerpt=(
                "Tree Species Details: Family Name Alangiaceae; Genus Alangium; "
                "Species Alangium salviifolium; Key character: floral symmetry: "
                "actinomorphic."
            ),
            quality="medium",
            provider="Bangladesh Forest Information System",
            url="https://bfis.bforest.gov.bd/nef/index.php/data/dataSpecies#my43",
            title="Forest Emission Factor Database - Alangium salviifolium",
            citation=(
                "Bangladesh Forest Information System, Forest Emission Factor "
                "Database species record: Alangium salviifolium"
            ),
            record_id=(
                "bfis:forest-emission-factor:alangium-salviifolium:floral-symmetry"
            ),
            lineage="bfis:species-record:alangium-salviifolium",
            lineage_method="canonical_government_database_species_record",
            source_tier="A",
            source_type="government_forest_species_database",
            domain="bfis.bforest.gov.bd",
            content_sha256=(
                "fdfc4f564792bfc59a271b9c05dfa65788e8901500afb11a9e99cd90141bc49b"
            ),
            content_sha256_basis="retrieved_government_database_html_bytes",
            name_resolution_lineage=(
                "master_accepted_name_exact; historical Alangiaceae placement "
                "reconciled to current Cornaceae"
            ),
        )
    )

    rows.append(
        _row(
            species="Sideritis canariensis",
            trait="floral_symmetry",
            value="zygomorphic",
            raw_value="corola tubulosa ... con pequeños labios",
            excerpt=(
                "Sideritis canariensis L. Cada flor presenta una pequeña corola "
                "tubulosa de color blanquecino amarillento, con pequeños labios "
                "de tonalidad marrón rojiza en la marchitez."
            ),
            quality="medium",
            provider="Gobierno de Canarias CanariWiki",
            url=(
                "https://www3.gobiernodecanarias.org/medusa/wiki/"
                "index.php?title=Chajorra_de_monte"
            ),
            title="Chajorra de monte - CanariWiki",
            citation=(
                "Gobierno de Canarias, CanariWiki species account: "
                "Sideritis canariensis L."
            ),
            record_id="canariwiki:sideritis-canariensis:bilabiate-corolla",
            lineage=(
                "url:https://www3.gobiernodecanarias.org/medusa/wiki/"
                "index.php?title=Chajorra_de_monte"
            ),
            lineage_method="canonical_government_species_description_page",
            source_tier="A",
            source_type="government_regional_flora_species_description",
            domain="www3.gobiernodecanarias.org",
            content_sha256=(
                "c415210616906d709e4bc4d5410f45753558dc7d179a66183d43a8c499437e35"
            ),
            content_sha256_basis="retrieved_government_species_page_html_bytes",
            language="es",
        )
    )

    commelineae_excerpt = (
        "Breeding system studies in the selected taxa revealed that they were all "
        "self and cross-compatible species. Except Dictyospermum montanum (2.5%) "
        "and Floscopa scandens (14%), all others showed high percentages of "
        "autogamous fruit set (72.38% in Commelina diffusa, 47.62% in Murdannia "
        "nudiflora and 41.43% in Rhopalephora scaberrima). These species thus have "
        "evolved a system of adaptive modifications which resulted in a mixed "
        "mating system that produces both self and cross seeds."
    )
    commelineae_common = {
        "provider": "Veena 2020 University of Calicut doctoral thesis",
        "url": (
            "https://scholar.uoc.ac.in/bitstreams/"
            "dd6a8802-5ffd-4a70-9078-736a78794e8a/download"
        ),
        "title": (
            "Pollination biology of selected taxa of the tribe Commelineae "
            "(Commelinaceae)"
        ),
        "citation": (
            "Veena, V. (2020), doctoral thesis, Department of Botany, "
            "University of Calicut, pp. 174-176 and 181-183"
        ),
        "lineage": "study:veena-2020:commelineae-controlled-pollination",
        "lineage_method": (
            "underlying_multispecies_thesis_experiment; later species articles "
            "must be reconciled to this study lineage"
        ),
        "source_tier": "A",
        "source_type": "doctoral_thesis_controlled_pollination_experiment",
        "domain": "scholar.uoc.ac.in",
        "content_sha256": (
            "1d40da559fe5e601e05fc9fc7d678b7b4d4c18a6f22462e37b212b00fc2263bd"
        ),
        "content_sha256_basis": "downloaded_official_repository_pdf_bytes",
    }
    commelineae_species = (
        "Commelina diffusa",
        "Dictyospermum montanum",
        "Floscopa scandens",
        "Murdannia nudiflora",
        "Rhopalephora scaberrima",
    )
    for species in commelineae_species:
        slug = species.casefold().replace(" ", "-")
        rows.extend(
            [
                _row(
                    species=species,
                    trait="self_incompatibility",
                    value="SC",
                    raw_value="self- and cross-compatible in controlled pollinations",
                    excerpt=commelineae_excerpt,
                    quality="high",
                    record_id=f"uoc:veena-2020:{slug}:self-compatible",
                    **commelineae_common,
                ),
                _row(
                    species=species,
                    trait="mating_system",
                    value="mixed_mating",
                    raw_value="mixed mating system producing self and cross seeds",
                    excerpt=commelineae_excerpt,
                    quality="high",
                    record_id=f"uoc:veena-2020:{slug}:mixed-mating",
                    **commelineae_common,
                ),
            ]
        )
    return rows


def _review_audit(evidence: pd.DataFrame) -> pd.DataFrame:
    audit = _audit(evidence)
    audit["reviewer"] = REVIEWER
    audit["reviewed_at_utc"] = CREATED_AT
    audit["decision_reason"] = evidence.apply(
        lambda row: (
            "Accepted after exact species identity, source page or original-study "
            "statement, trait-specific ontology, lineage and content fingerprint "
            f"review; mapping retained only for {row['trait_name']}."
        ),
        axis=1,
    )
    return audit


def _sha256(path: Path) -> str:
    return hashlib.sha256(_canonical_file_bytes(path)).hexdigest()


def build(
    *,
    master_csv: Path,
    output_dir: Path,
    prior_curated_evidence_csv: Path | None = None,
    prior_curated_audit_csv: Path | None = None,
) -> dict[str, object]:
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    master_family = master.set_index("accepted_species")["family"].to_dict()
    expected_families = {
        "Myrcia guianensis": "Myrtaceae",
        "Melastoma affine": "Melastomataceae",
        "Hosta ventricosa": "Asparagaceae",
        "Ornithogalum thyrsoides": "Asparagaceae",
        "Melastoma malabathricum": "Melastomataceae",
        "Drypetes assamica": "Putranjivaceae",
        "Cornus sericea": "Cornaceae",
        "Daucus montanus": "Apiaceae",
        "Commelina diffusa": "Commelinaceae",
        "Dictyospermum montanum": "Commelinaceae",
        "Floscopa scandens": "Commelinaceae",
        "Murdannia nudiflora": "Commelinaceae",
        "Rhopalephora scaberrima": "Commelinaceae",
        "Palaquium obovatum": "Sapotaceae",
        "Alangium salviifolium": "Cornaceae",
        "Sideritis canariensis": "Lamiaceae",
    }
    missing = sorted(set(expected_families) - set(master_family))
    if missing:
        raise ValueError(f"reviewed species absent from target master: {missing}")
    conflicts = {
        species: (family, master_family[species])
        for species, family in expected_families.items()
        if master_family[species] != family
    }
    if conflicts:
        raise ValueError(f"family conflicts in reviewed checkpoint: {conflicts}")

    evidence = pd.DataFrame(reviewed_rows(), columns=EVIDENCE_COLUMNS).fillna("")
    evidence = evidence.sort_values(
        ["accepted_species", "trait_name", "candidate_id"]
    ).reset_index(drop=True)
    if len(evidence) != 27:
        raise ValueError(f"expected 27 reviewed trait rows, found {len(evidence)}")
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("wave-2 candidate IDs are not unique")
    audit = _review_audit(evidence)

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "rule_unlock_wave2_evidence_20260812.csv"
    audit_path = output_dir / "rule_unlock_wave2_manual_audit_20260812.csv"
    evidence.to_csv(evidence_path, index=False, lineterminator="\n")
    audit.to_csv(audit_path, index=False, lineterminator="\n")
    outputs = [evidence_path, audit_path]

    combined_evidence: pd.DataFrame | None = None
    combined_audit: pd.DataFrame | None = None
    if prior_curated_evidence_csv is not None:
        if prior_curated_audit_csv is None:
            raise ValueError("prior curated audit is required with prior evidence")
        prior_evidence = pd.read_csv(prior_curated_evidence_csv, dtype=str).fillna("")
        prior_audit = pd.read_csv(prior_curated_audit_csv, dtype=str).fillna("")
        owned = set(evidence["candidate_id"])
        combined_evidence = pd.concat(
            [prior_evidence.loc[~prior_evidence["candidate_id"].isin(owned)], evidence],
            ignore_index=True,
        )
        combined_audit = pd.concat(
            [prior_audit.loc[~prior_audit["candidate_id"].isin(owned)], audit],
            ignore_index=True,
        )
        for name, frame in (("evidence", combined_evidence), ("audit", combined_audit)):
            if frame["candidate_id"].duplicated().any():
                raise ValueError(f"combined {name} candidate IDs are not unique")
        combined_evidence_path = output_dir / "combined_curated_evidence_20260812.csv"
        combined_audit_path = output_dir / "combined_curated_manual_audit_20260812.csv"
        combined_evidence.to_csv(combined_evidence_path, index=False, lineterminator="\n")
        combined_audit.to_csv(combined_audit_path, index=False, lineterminator="\n")
        outputs.extend([combined_evidence_path, combined_audit_path])

    summary: dict[str, object] = {
        "contract": "trait_specific_rule_unlock_wave2_individually_reviewed_v1",
        "created_at_utc": CREATED_AT,
        "new_evidence_rows": len(evidence),
        "new_species": int(evidence["accepted_species"].nunique()),
        "new_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "quality_counts": evidence["evidence_quality"].value_counts().to_dict(),
        "axis_counts": evidence["axis"].value_counts().to_dict(),
        "trait_counts": evidence["trait_name"].value_counts().to_dict(),
        "source_lineages": int(evidence["source_lineage"].nunique()),
        "audit": {
            "reviewed": len(audit),
            "accepted_correct": int(audit["decision"].str.casefold().eq("accept").sum()),
            "precision": float(audit["decision"].str.casefold().eq("accept").mean()),
            "cultivar_contamination_rate": float(
                audit["cultivar_contamination"].str.casefold().eq("true").mean()
            ),
        },
        "guardrails": {
            "trait_specific_records": True,
            "genus_inference_emitted_here": False,
            "family_inference": False,
            "global_fallback": False,
            "n2_formal_inference": False,
            "cross_trait_substitution": False,
            "search_snippet_evidence": False,
        },
        "files": {},
    }
    if combined_evidence is not None and combined_audit is not None:
        summary["combined"] = {
            "evidence_rows": len(combined_evidence),
            "audit_rows": len(combined_audit),
            "species": int(combined_evidence["accepted_species"].nunique()),
            "species_trait": int(
                combined_evidence[["accepted_species", "trait_name"]]
                .drop_duplicates()
                .shape[0]
            ),
        }
    readme = output_dir / "README.md"
    if readme.exists():
        outputs.append(readme)
    for path in outputs:
        summary["files"][path.name] = {
            "sha256": _sha256(path),
            "size_bytes": len(_canonical_file_bytes(path)),
        }
    manifest = output_dir / "rule_unlock_wave2_manifest_20260812.json"
    manifest.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--prior-curated-evidence-csv", type=Path)
    parser.add_argument("--prior-curated-audit-csv", type=Path)
    print(json.dumps(build(**vars(parser.parse_args())), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
