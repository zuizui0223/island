"""Freeze the second trait-specific genus-rule unlock evidence wave.

The checkpoint contains only individually reviewed species-direct statements.
Reproductive traits remain separate, and every morphology record is a direct
species description.  Genus inference is rebuilt later by the shared
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
INDIA_FLORA_REVIEWED_PATH = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "rule_unlock_wave2_checkpoint_20260812/"
    "india_flora_online_reviewed_records_20260812.csv"
)
INDIA_FLORA_FULL_AUDIT_PATH = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "rule_unlock_wave2_checkpoint_20260812/"
    "india_flora_online_full_candidate_audit_20260812.csv"
)
INDIA_FLORA_SYNONYM_REVIEWED_PATH = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "rule_unlock_wave2_checkpoint_20260812/"
    "india_flora_online_synonym_reviewed_records_20260812.csv"
)
INDIA_FLORA_SYNONYM_FULL_AUDIT_PATH = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "rule_unlock_wave2_checkpoint_20260812/"
    "india_flora_online_synonym_full_candidate_audit_20260812.csv"
)
PROTA_REVIEWED_PATH = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "rule_unlock_wave2_checkpoint_20260812/"
    "plantuse_prota_reviewed_records_20260812.csv"
)
PROTA_FULL_AUDIT_PATH = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "rule_unlock_wave2_checkpoint_20260812/"
    "plantuse_prota_full_candidate_audit_20260812.csv"
)

BFIS_BULK_SYMMETRY_ROWS = (
    ("Actephila excelsa", "Euphorbiaceae", "actinomorphic"),
    ("Alchornea tiliifolia", "Euphorbiaceae", "actinomorphic"),
    ("Antidesma acidum", "Euphorbiaceae", "actinomorphic"),
    ("Antidesma ghaesembilla", "Euphorbiaceae", "actinomorphic"),
    ("Antidesma montanum", "Euphorbiaceae", "actinomorphic"),
    ("Antidesma roxburghii", "Euphorbiaceae", "actinomorphic"),
    ("Antidesma velutinum", "Euphorbiaceae", "actinomorphic"),
    ("Aporosa aurea", "Euphorbiaceae", "actinomorphic"),
    ("Baccaurea ramiflora", "Euphorbiaceae", "actinomorphic"),
    ("Barringtonia acutangula", "Lecythidaceae", "actinomorphic"),
    ("Breynia retusa", "Euphorbiaceae", "actinomorphic"),
    ("Bridelia retusa", "Euphorbiaceae", "actinomorphic"),
    ("Bridelia tomentosa", "Euphorbiaceae", "actinomorphic"),
    ("Careya arborea", "Lecythidaceae", "actinomorphic"),
    ("Casearia tomentosa", "Flacourtiaceae", "actinomorphic"),
    ("Chaetocarpus castanocarpus", "Euphorbiaceae", "actinomorphic"),
    ("Cinnamomum iners", "Lauraceae", "actinomorphic"),
    ("Cleidion javanicum", "Euphorbiaceae", "actinomorphic"),
    ("Couroupita guianensis", "Lecythidaceae", "actinomorphic"),
    ("Croton aromaticus", "Euphorbiaceae", "actinomorphic"),
    ("Croton joufra", "Euphorbiaceae", "actinomorphic"),
    ("Croton tiglium", "Euphorbiaceae", "actinomorphic"),
    ("Endospermum chinense", "Euphorbiaceae", "actinomorphic"),
    ("Engelhardtia roxburghiana", "Juglandaceae", "actinomorphic"),
    ("Engelhardtia spicata", "Juglandaceae", "actinomorphic"),
    ("Epiprinus siletianus", "Euphorbiaceae", "actinomorphic"),
    ("Euphorbia antiquorum", "Euphorbiaceae", "actinomorphic"),
    ("Euphorbia cotinifolia", "Euphorbiaceae", "actinomorphic"),
    ("Euphorbia neriifolia", "Euphorbiaceae", "actinomorphic"),
    ("Euphorbia tirucalli", "Euphorbiaceae", "actinomorphic"),
    ("Excoecaria oppositifolia", "Euphorbiaceae", "actinomorphic"),
    ("Falconeria insignis", "Euphorbiaceae", "actinomorphic"),
    ("Flacourtia jangomas", "Flacourtiaceae", "actinomorphic"),
    ("Flueggea leucopyrus", "Euphorbiaceae", "actinomorphic"),
    ("Flueggea virosa", "Euphorbiaceae", "actinomorphic"),
    ("Gomphandra tetrandra", "Icacinaceae", "actinomorphic"),
    ("Homonoia riparia", "Euphorbiaceae", "actinomorphic"),
    ("Jatropha multifida", "Euphorbiaceae", "actinomorphic"),
    ("Macaranga denticulata", "Euphorbiaceae", "actinomorphic"),
    ("Macaranga peltata", "Euphorbiaceae", "actinomorphic"),
    ("Mallotus nudiflorus", "Euphorbiaceae", "actinomorphic"),
    ("Mallotus philippensis", "Euphorbiaceae", "actinomorphic"),
    ("Mallotus repandus", "Euphorbiaceae", "actinomorphic"),
    ("Mallotus tetracoccus", "Euphorbiaceae", "actinomorphic"),
    ("Margaritaria indica", "Euphorbiaceae", "actinomorphic"),
    ("Ostodes paniculata", "Euphorbiaceae", "actinomorphic"),
    ("Shirakiopsis indica", "Euphorbiaceae", "actinomorphic"),
    ("Sophora wightii", "Fabaceae", "zygomorphic"),
    ("Triadica sebifera", "Euphorbiaceae", "actinomorphic"),
)


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
    retrieved_at_utc: str = CREATED_AT,
    cultivar_status: str = "wild_or_species_level_statement_not_cultivar_limited",
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
        retrieved_at_utc=retrieved_at_utc,
        raw_value=raw_value,
    )
    row["source_group"] = SOURCE_GROUP
    row["language"] = language
    row["name_resolution_lineage"] = name_resolution_lineage
    row["wild_cultivated_cultivar_status"] = cultivar_status
    row["query"] = "current_support_2_rule_unlock_original_or_direct_species_source"
    return row


def _bfis_bulk_symmetry_rows() -> list[dict[str, str]]:
    """Return novel master-species records from one immutable BFIS snapshot.

    The source page is a government database index containing 281 separate
    species treatments.  Only exact master-name matches with an explicit
    ``floral symmetry`` field and no pre-existing direct symmetry evidence in
    the formal checkpoint are frozen here.  Historical source-family labels
    are retained in the quote; current family identity is checked separately
    against the target master in :func:`build`.
    """

    rows: list[dict[str, str]] = []
    for species, source_family, value in BFIS_BULK_SYMMETRY_ROWS:
        slug = species.casefold().replace(" ", "-")
        rows.append(
            _row(
                species=species,
                trait="floral_symmetry",
                value=value,
                raw_value=f"floral symmetry: {value}",
                excerpt=(
                    f"Tree Species Details: Family Name: {source_family}; Genus: "
                    f"{species.split()[0]}; Species: {species}; floral symmetry: "
                    f"{value}."
                ),
                quality="medium",
                provider="Bangladesh Forest Information System",
                url=(
                    "https://bfis.bforest.gov.bd/nef/index.php/data/"
                    "dataSpecies/40"
                ),
                title=f"Forest Emission Factor Database - {species}",
                citation=(
                    "Bangladesh Forest Information System, Forest Emission "
                    f"Factor Database species record: {species}"
                ),
                record_id=(
                    f"bfis:forest-emission-factor:{slug}:floral-symmetry"
                ),
                lineage=(
                    "provider_treatment:bangladesh-forest-information-system:"
                    f"{slug}"
                ),
                lineage_method="canonical_government_database_species_record",
                source_tier="A",
                source_type="government_forest_database_species_record",
                domain="bfis.bforest.gov.bd",
                content_sha256=(
                    "80f4d97075ebb66d41b5707b1ae4bc330dd41ea4df567b4e"
                    "a19f63531aa88b4f"
                ),
                content_sha256_basis=(
                    "downloaded_full_government_database_html_bytes"
                ),
                retrieved_at_utc="2026-08-12T12:24:00Z",
                name_resolution_lineage=(
                    "master_accepted_name_exact; source_family_label_retained_"
                    "and_current_master_family_checked"
                ),
            )
        )
    return rows


def _india_flora_online_rows() -> list[dict[str, str]]:
    """Return the 250 accepted records from the full 262-candidate IISc audit."""

    reviewed = pd.read_csv(INDIA_FLORA_REVIEWED_PATH, dtype=str).fillna("")
    required = {
        "candidate_id",
        "accepted_species",
        "page_id",
        "axis",
        "trait_name",
        "normalized_value",
        "supporting_excerpt",
        "source_url",
        "page_content_sha256",
        "cultivar_status",
    }
    missing = required - set(reviewed.columns)
    if missing:
        raise ValueError(f"IISc reviewed records missing columns: {sorted(missing)}")
    if len(reviewed) != 250 or reviewed["candidate_id"].duplicated().any():
        raise ValueError("IISc checkpoint must contain 250 unique accepted candidates")
    if not reviewed["page_content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all():
        raise ValueError("IISc checkpoint contains invalid page-content hashes")

    rows: list[dict[str, str]] = []
    for record in reviewed.to_dict("records"):
        species = record["accepted_species"]
        page_id = record["page_id"]
        trait = record["trait_name"]
        excerpt = f"Key identification features: {record['supporting_excerpt']}"
        row = _row(
            species=species,
            trait=trait,
            value=record["normalized_value"],
            raw_value=record["supporting_excerpt"],
            excerpt=excerpt,
            quality="medium",
            provider="IISc India Flora Online",
            url=record["source_url"],
            title=f"India Flora Online - {species}",
            citation=(
                "Sankara Rao, K. and Deepak Kumar (2026). India Flora Online, "
                f"species treatment: {species}."
            ),
            record_id=f"iisc-india-flora-online:{page_id}:{trait}",
            lineage=f"provider_treatment:iisc-india-flora-online:{page_id}",
            lineage_method="canonical_university_flora_species_treatment",
            source_tier="A",
            source_type="university_herbarium_regional_flora_species_treatment",
            domain="indiaflora-ces.iisc.ac.in",
            content_sha256=record["page_content_sha256"],
            content_sha256_basis="downloaded_complete_species_treatment_html_bytes",
            retrieved_at_utc="2026-08-12T13:02:00Z",
            cultivar_status=record["cultivar_status"],
        )
        row["query"] = "bulk_exact_master_unresolved_morphology_index"
        rows.append(row)
    return rows


def _india_flora_online_synonym_rows() -> list[dict[str, str]]:
    """Return records whose IISc page names agree in both WFO and GBIF."""

    reviewed = pd.read_csv(INDIA_FLORA_SYNONYM_REVIEWED_PATH, dtype=str).fillna("")
    required = {
        "candidate_id",
        "accepted_species",
        "searched_name",
        "page_id",
        "family",
        "wfo_taxon_id",
        "wfo_accepted_usage_id",
        "gbif_usage_key",
        "gbif_accepted_usage_key",
        "trait_name",
        "normalized_value",
        "supporting_excerpt",
        "source_url",
        "page_content_sha256",
        "name_match_method",
        "cultivar_status",
    }
    missing = required - set(reviewed.columns)
    if missing:
        raise ValueError(f"IISc synonym records missing columns: {sorted(missing)}")
    if len(reviewed) != 70 or reviewed["candidate_id"].duplicated().any():
        raise ValueError("IISc synonym checkpoint must contain 70 unique records")
    if not reviewed["name_match_method"].eq("exact_synonym").all():
        raise ValueError("IISc synonym checkpoint lacks two-backbone agreement")
    if not reviewed["page_content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all():
        raise ValueError("IISc synonym checkpoint contains invalid page hashes")

    rows: list[dict[str, str]] = []
    for record in reviewed.to_dict("records"):
        species = record["accepted_species"]
        searched_name = record["searched_name"]
        page_id = record["page_id"]
        trait = record["trait_name"]
        row = _row(
            species=species,
            trait=trait,
            value=record["normalized_value"],
            raw_value=record["supporting_excerpt"],
            excerpt=f"Key identification features: {record['supporting_excerpt']}",
            quality="medium",
            provider="IISc India Flora Online",
            url=record["source_url"],
            title=f"India Flora Online - {searched_name}",
            citation=(
                "Sankara Rao, K. and Deepak Kumar (2026). India Flora Online, "
                f"species treatment: {searched_name}."
            ),
            record_id=f"iisc-india-flora-online:{page_id}:{trait}",
            lineage=f"provider_treatment:iisc-india-flora-online:{page_id}",
            lineage_method="canonical_university_flora_species_treatment",
            source_tier="A",
            source_type="university_herbarium_regional_flora_species_treatment",
            domain="indiaflora-ces.iisc.ac.in",
            content_sha256=record["page_content_sha256"],
            content_sha256_basis="downloaded_complete_species_treatment_html_bytes",
            retrieved_at_utc="2026-08-12T14:29:00Z",
            cultivar_status=record["cultivar_status"],
            name_resolution_lineage=(
                f"wfo_june_2026:{record['wfo_taxon_id']}:"
                f"{record['wfo_accepted_usage_id']};gbif:"
                f"{record['gbif_usage_key']}:{record['gbif_accepted_usage_key']}"
            ),
        )
        row["name_match_method"] = record["name_match_method"]
        row["matched_page_name"] = searched_name
        row["evidence_scope"] = "synonym_direct"
        row["query"] = "bulk_two_backbone_synonym_unresolved_morphology_index"
        rows.append(row)
    return rows


def _prota_monograph_rows() -> list[dict[str, str]]:
    """Return the 448 accepted records from the full PROTA candidate audit."""

    reviewed = pd.read_csv(PROTA_REVIEWED_PATH, dtype=str).fillna("")
    required = {
        "candidate_id",
        "accepted_species",
        "family",
        "page_id",
        "revision_id",
        "revision_timestamp",
        "axis",
        "trait_name",
        "normalized_value",
        "supporting_excerpt",
        "page_title",
        "source_url",
        "citation",
        "page_content_sha256",
        "cultivar_status",
    }
    missing = required - set(reviewed.columns)
    if missing:
        raise ValueError(f"PROTA reviewed records missing columns: {sorted(missing)}")
    if len(reviewed) != 448 or reviewed["candidate_id"].duplicated().any():
        raise ValueError("PROTA checkpoint must contain 448 unique accepted candidates")
    if not reviewed["page_content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all():
        raise ValueError("PROTA checkpoint contains invalid revision-content hashes")
    if reviewed["citation"].str.len().lt(100).any():
        raise ValueError("PROTA checkpoint contains incomplete monograph citations")

    rows: list[dict[str, str]] = []
    for record in reviewed.to_dict("records"):
        species = record["accepted_species"]
        page_id = record["page_id"]
        revision_id = record["revision_id"]
        trait = record["trait_name"]
        lineage = f"provider_treatment:prota:{page_id}"
        if species == "Adenia cissampeloides":
            lineage = "monograph:prota-11-1:adenia-cissampeloides"
        row = _row(
            species=species,
            trait=trait,
            value=record["normalized_value"],
            raw_value=record["supporting_excerpt"],
            excerpt=record["supporting_excerpt"],
            quality="high",
            provider="Plant Resources of Tropical Africa (PROTA)",
            url=record["source_url"],
            title=record["page_title"],
            citation=record["citation"],
            record_id=f"prota-mediawiki:{page_id}:{revision_id}:{trait}",
            lineage=lineage,
            lineage_method="canonical_authored_prota_species_monograph_revision",
            source_tier="A",
            source_type="expert_botanical_monograph_species_description",
            domain="plantuse.plantnet.org",
            content_sha256=record["page_content_sha256"],
            content_sha256_basis="mediawiki_revision_wikitext_utf8_bytes",
            retrieved_at_utc="2026-08-12T14:10:00Z",
            cultivar_status=record["cultivar_status"],
        )
        row["query"] = "bulk_prota_api_exact_master_unresolved"
        rows.append(row)
    return rows


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

    rows.append(
        _row(
            species="Adenia cissampeloides",
            trait="floral_symmetry",
            value="actinomorphic",
            raw_value="flowers regular",
            excerpt=(
                "Adenia cissampeloides (Planch. ex Hook.) Harms. Flowers "
                "unisexual, regular, 5-merous, pale greenish; pedicel "
                "2-10(-15) mm long in male flowers, slightly shorter in female ones."
            ),
            quality="high",
            provider="Plant Resources of Tropical Africa (PROTA)",
            url=(
                "https://plantuse.plantnet.org/en/"
                "Adenia_cissampeloides_%28PROTA%29"
            ),
            title="Adenia cissampeloides (PROTA)",
            citation=(
                "Grace, O.M. & Fowler, D.G. (2007), Adenia cissampeloides, "
                "Plant Resources of Tropical Africa 11(1): Medicinal Plants 1"
            ),
            record_id="prota:adenia-cissampeloides:flowers-regular",
            lineage="monograph:prota-11-1:adenia-cissampeloides",
            lineage_method="original_expert_monograph_species_account",
            source_tier="A",
            source_type="expert_botanical_monograph_species_description",
            domain="plantuse.plantnet.org",
            content_sha256=(
                "74ba4d8972a69f01339a05cf676f2638ae7c8ed66400224e3a57d17df9e44339"
            ),
            content_sha256_basis="retrieved_prota_species_account_html_bytes",
        )
    )

    tristaniopsis_symmetry_excerpt = (
        "Syzygium floribundum, Syzygium smithii and Tristaniopsis laurina "
        "are small to medium-sized trees. The three species conform to the "
        "'general entomophilous' flower structure in which individual flowers "
        "are small with little depth effect in the corolla, normally do not "
        "possess nectar guides, and are actinomorphic in form."
    )
    rows.append(
        _row(
            species="Tristaniopsis laurina",
            trait="floral_symmetry",
            value="actinomorphic",
            raw_value="actinomorphic in form",
            excerpt=tristaniopsis_symmetry_excerpt,
            quality="high",
            provider="Williams and Adam 2019 Cunninghamia",
            url=(
                "https://www.botanicgardens.org.au/sites/default/files/2023-06/"
                "BGD0541_Cunninghamia-%C2%A0WILLIAMS-ADAM-Myrtaceae.pdf"
            ),
            title=(
                "A Preliminary Checklist of Flower-visiting Insects from "
                "Syzygium floribundum, Syzygium smithii and Tristaniopsis laurina"
            ),
            citation=(
                "Williams, G. & Adam, P. (2019), Cunninghamia 19:57-74, p. 58"
            ),
            record_id="cunninghamia-19:tristaniopsis-laurina:actinomorphic",
            lineage="publication:williams-adam-2019-cunninghamia-19-57-74",
            lineage_method="original_peer_reviewed_species_field_study",
            source_tier="A",
            source_type="peer_reviewed_primary_flower_visitor_study",
            domain="www.botanicgardens.org.au",
            content_sha256=(
                "4ce687100a9a10d72dfec9cc5a5d136ad6d0537decc99cab368945940932307b"
            ),
            content_sha256_basis="downloaded_publisher_pdf_bytes",
        )
    )

    rows.append(
        _row(
            species="Tristaniopsis laurina",
            trait="autonomous_selfing_capacity",
            value="absent",
            raw_value="bagged flowers 189; developing fruit 0; open fruit set 19.1%",
            excerpt=(
                "Table 2. Autogamy (automatic self-pollination) and open "
                "pollination results in lowland rainforest species. "
                "Tristaniopsis laurina: 1 plant; bagged flowers 189; developing "
                "fruit 0; fruit set 0%; open flowers 210; developing fruit 40; "
                "fruit set 19.1%."
            ),
            quality="high",
            provider="Adam and Williams 2001 Cunninghamia",
            url=(
                "https://www.botanicgardens.org.au/sites/default/files/2023-06/"
                "Volume-7%281%29-2001-Cun7Ada089-100_0.pdf"
            ),
            title=(
                "Dioecy, self-compatibility and vegetative reproduction in "
                "Australian subtropical rainforest trees and shrubs"
            ),
            citation=(
                "Adam, P. & Williams, G. (2001), Cunninghamia 7(1):89-100, "
                "Table 2, p. 94"
            ),
            record_id="cunninghamia-7-1:table-2:tristaniopsis-laurina:autogamy",
            lineage="publication:adam-williams-2001-cunninghamia-7-89-100",
            lineage_method="original_primary_bagging_experiment",
            source_tier="A",
            source_type="peer_reviewed_primary_field_bagging_experiment",
            domain="www.botanicgardens.org.au",
            content_sha256=(
                "cc63c25c9759692479f98ad2c5bc7707d00022d3ade911762a5486a3f8c431fe"
            ),
            content_sha256_basis="downloaded_publisher_pdf_bytes",
        )
    )

    rows.append(
        _row(
            species="Polycarpaea corymbosa",
            trait="inflorescence_display",
            value="umbel_corymb",
            raw_value="flowers in terminal or axillary corymbose cymes",
            excerpt=(
                "Polycarpaea corymbosa (L.) Lamk. Erect, branched herbs, "
                "tomentose. Leaves fascicled at nodes, 1-2 x 0.1 cm, "
                "linear-lanceolate. Flowers in terminal or axillary corymbose "
                "cymes; sepals silvery white or pink, longer than petals, "
                "petals about 0.12 cm long. Capsules ovoid, 3-valved."
            ),
            quality="high",
            provider="Botanical Survey of India - Flora of Sindhudurg",
            url=(
                "https://bsi.gov.in/uploads/documents/Public_Information/"
                "publication/books/district_flora_latest/Flora%20of%20Sindhudurg.pdf"
            ),
            title="Flora of Sindhudurg",
            citation=(
                "Kulkarni, B.G. (1988), Flora of Sindhudurg, Botanical Survey "
                "of India, p. 33"
            ),
            record_id=(
                "bsi:flora-of-sindhudurg:p33:polycarpaea-corymbosa:"
                "corymbose-cymes"
            ),
            lineage="publication:kulkarni-1988-flora-of-sindhudurg",
            lineage_method="original_government_district_flora_species_treatment",
            source_tier="A",
            source_type="government_flora_species_description",
            domain="bsi.gov.in",
            content_sha256=(
                "26b72112fa4cc1e796e0c2fcadd44fd3345aacc53b6b7d09f93a6578db44b12f"
            ),
            content_sha256_basis="downloaded_official_government_flora_pdf_bytes",
        )
    )

    rows.append(
        _row(
            species="Boronia muelleri",
            trait="floral_form",
            value="open_radial",
            raw_value="pink with star-like flowers",
            excerpt=(
                "At the entrance to the Sydney Region Flora, Boronia muelleri "
                "[Section 191] is of medium size, pink with star-like flowers."
            ),
            quality="medium",
            provider="Australian National Botanic Gardens",
            url=(
                "https://www.anbg.gov.au/gardens/visiting/iftw/iftw-archive/"
                "iftw-2000-10-27.html"
            ),
            title="In Flower This Week - 27 October 2000",
            citation=(
                "Australian National Botanic Gardens, In Flower This Week, "
                "27 October 2000, updated by Murray Fagg"
            ),
            record_id="anbg:iftw-2000-10-27:boronia-muelleri:star-like",
            lineage="anbg:iftw-2000-10-27:boronia-muelleri",
            lineage_method="official_botanic_garden_species_observation_page",
            source_tier="A",
            source_type="official_botanic_garden_species_observation",
            domain="www.anbg.gov.au",
            content_sha256=(
                "18efaaace71baeaecc42ced86eff645b14caeb0f06a33490b1fa2c4abc6ad0b4"
            ),
            content_sha256_basis="retrieved_official_botanic_garden_html_bytes",
        )
    )

    rows.append(
        _row(
            species="Benstonea foetida",
            trait="inflorescence_display",
            value="raceme_spike_panicle",
            raw_value="spadix branched; male flowers in 5-6 cm long catkins",
            excerpt=(
                "Benstonea foetida (Roxb.) Callm. & Buerki. Shrubs with stilt "
                "roots. Leaves spirally arranged; linear-lanceolate, acuminate, "
                "spiny along the margins and midrib below, glabrous, to 1.5 m "
                "long. Spadix branched, spathe ovate, acute. Male flowers in "
                "5-6 cm long catkins; filaments 1 cm long."
            ),
            quality="medium",
            provider="eFlora of India",
            url="https://efloraofindia.com/efi/pandanus-foetidus/",
            title="Benstonea foetida - eFlora of India",
            citation=(
                "eFlora of India species account for Benstonea foetida; "
                "description attributed to Dr. N. Sasidharan, Kerala Forest "
                "Research Institute, via India Biodiversity Portal"
            ),
            record_id="efloraindia:benstonea-foetida:branched-spadix-catkins",
            lineage=(
                "india-biodiversity-portal:dr-n-sasidharan:benstonea-foetida"
            ),
            lineage_method=(
                "credited_expert_species_description_republished_by_efloraindia"
            ),
            source_tier="B",
            source_type="specialist_regional_flora_species_description",
            domain="efloraofindia.com",
            content_sha256=(
                "05f86f645f42189f94225dd02d4a2d890685ebb0f8215a21aa7c58f45335deab"
            ),
            content_sha256_basis="retrieved_source_page_html_bytes",
        )
    )

    rows.append(
        _row(
            species="Pleioluma balansana",
            trait="flower_primary_color",
            value="white",
            raw_value="corolla cream with white-edged corolla lobes",
            excerpt=(
                "Pleioluma balansana (Pierre ex Baill.) Swenson & Munzinger. "
                "Flowers 5-merous, 1-5 per fascicle; pedicel 2-5 mm long, "
                "usually at least 0.7 mm wide, tomentulose. Sepals 1.0-2.0 "
                "(-2.5) mm long with the same indument as the pedicel. Corolla "
                "cream, with white-edged corolla lobes, 2.5-3.5 mm long."
            ),
            quality="high",
            provider="Swenson et al. 2018 Australian Systematic Botany",
            url="https://doi.org/10.1071/SB17040",
            title=(
                "Phylogeny, species delimitation and revision of Pleioluma "
                "(Sapotaceae) in New Caledonia, a frequently gynodioecious genus"
            ),
            citation=(
                "Swenson, Nylander & Munzinger (2018), Australian Systematic "
                "Botany 31:120-165, DOI 10.1071/SB17040, p. 135"
            ),
            record_id="doi:10.1071/SB17040:p135:pleioluma-balansana:corolla-colour",
            lineage="doi:10.1071/SB17040",
            lineage_method="peer_reviewed_primary_taxonomic_revision_doi",
            source_tier="A",
            source_type="peer_reviewed_primary_taxonomic_revision",
            domain="doi.org",
            content_sha256=(
                "1eb462ba75aee764fc20e77be96650b6647aa0de2667d79584b1a079e33a9997"
            ),
            content_sha256_basis="downloaded_author_copy_of_published_pdf_bytes",
        )
    )

    rows.append(
        _row(
            species="Kunzea ericoides",
            trait="floral_form",
            value="open_radial",
            raw_value="tiny white star-shaped flowers",
            excerpt=(
                "Kunzea ericoides, commonly known as Kanuka, is a versatile and "
                "hardy evergreen shrub or small tree native to New Zealand. This "
                "plant is known for its small, aromatic, needle-like leaves and "
                "profusion of tiny, white, star-shaped flowers that bloom from "
                "late spring to summer."
            ),
            quality="medium",
            provider="The Plant Store New Zealand",
            url=(
                "https://www.theplantstore.co.nz/products/nz-natives/"
                "kunzea-ericoides-kanuka/"
            ),
            title="Kunzea ericoides | Kanuka | The Plant Store",
            citation=(
                "The Plant Store New Zealand species page for Kunzea ericoides "
                "(retrieved 2026-08-12)"
            ),
            record_id="theplantstore:kunzea-ericoides:star-shaped-flowers",
            lineage="url:https://www.theplantstore.co.nz/products/nz-natives/kunzea-ericoides-kanuka",
            lineage_method="canonical_species_page_url",
            source_tier="C",
            source_type="specialist_nursery_species_description_morphology_only",
            domain="www.theplantstore.co.nz",
            content_sha256=(
                "689ef8e0e605868cda12c5e77ba6341d9cc633768901bad6c4e9e103b9237f5f"
            ),
            content_sha256_basis="retrieved_species_page_html_bytes",
        )
    )

    rows.append(
        _row(
            species="Corchorus olitorius",
            trait="floral_symmetry",
            value="actinomorphic",
            raw_value="flowers actinomorphic and regular",
            excerpt=(
                "Corchorus olitorius L. Taxonomy: Family Malvaceae; Genus "
                "Corchorus; Species Corchorus olitorius. Flowers: Complete, "
                "pedicellate, bracteate, small, bisexual, dichlamydeous, "
                "actinomorphic, regular, pentamerous, hypogynous, yellow in color."
            ),
            quality="medium",
            provider="Gupta et al. 2023 Journal of Pharmacognosy and Phytochemistry",
            url=(
                "https://www.phytojournal.com/archives/2023/vol12issue2/PartB/"
                "12-2-16-171.pdf"
            ),
            title=(
                "Morphological and pharmacological study of herbal medicine: "
                "Corchorus olitorius L."
            ),
            citation=(
                "Gupta et al. (2023), Journal of Pharmacognosy and "
                "Phytochemistry 12(2):84-87, p. 85"
            ),
            record_id="gupta-2023:p85:corchorus-olitorius:actinomorphic",
            lineage="publication:gupta-et-al-2023-corchorus-olitorius",
            lineage_method="published_species_morphology_article",
            source_tier="B",
            source_type="published_species_morphology_review",
            domain="www.phytojournal.com",
            content_sha256=(
                "09884189ef5cc681ed53663b1e16e8642df9ab4d6b5275f717e857065ae084f5"
            ),
            content_sha256_basis="downloaded_publisher_pdf_bytes",
        )
    )

    rows.append(
        _row(
            species="Jacobaea aquatica",
            trait="floral_form",
            value="composite_head",
            raw_value="yellow flowerheads with 12-15 rays",
            excerpt=(
                "Marsh Ragwort - Jacobaea aquatica. Flowerheads yellow, 25 to "
                "30 mm, with 12 to 15 rays and borne in loose flat topped "
                "clusters. Bracts not black tipped."
            ),
            quality="medium",
            provider="NatureSpot Leicestershire and Rutland",
            url="https://www.naturespot.org/species/marsh-ragwort",
            title="Marsh Ragwort | NatureSpot",
            citation=(
                "NatureSpot species account for Jacobaea aquatica, "
                "Leicestershire and Rutland (retrieved 2026-08-12)"
            ),
            record_id="naturespot:jacobaea-aquatica:yellow-rayed-flowerheads",
            lineage="url:https://www.naturespot.org/species/marsh-ragwort",
            lineage_method="canonical_specialist_regional_species_page_url",
            source_tier="B",
            source_type="specialist_regional_nature_recording_species_account",
            domain="www.naturespot.org",
            content_sha256=(
                "0bf88934cc5a74a1bac7f8e9fe33338e604fe5c65318ba439f570a430cf7f81c"
            ),
            content_sha256_basis="retrieved_species_page_html_bytes",
        )
    )

    rows.append(
        _row(
            species="Carpinus laxiflora",
            trait="inflorescence_display",
            value="raceme_spike_panicle",
            raw_value="pistillate catkins 5-16 cm long",
            excerpt=(
                "2. Carpinus laxiflora (Siebold & Zucc.) Blume, Mus. Bot. "
                "1: 309, 1851. Pistillate catkins 5-16 cm long; peduncle "
                "slender, sparsely pubescent; bracts 14-72 per catkin."
            ),
            quality="high",
            provider="National Institute of Biological Resources Flora of Korea",
            url=(
                "https://www.nibr.go.kr/aiibook/access/ecatalogt.jsp?"
                "Dir=824&callmode=admin&catimage=&eclang=ko&start=86&um=s"
            ),
            title="Flora of Korea, volume 2b - Hamamelidae",
            citation=(
                "Flora of Korea Editorial Committee, Flora of Korea vol. 2b, "
                "Betulaceae, pp. 73-74, National Institute of Biological Resources"
            ),
            record_id="nibr:flora-korea-v2b:p73:carpinus-laxiflora:catkins",
            lineage="nibr:flora-of-korea-v2b:carpinus-laxiflora",
            lineage_method="official_national_flora_species_treatment",
            source_tier="A",
            source_type="official_national_flora_species_treatment",
            domain="www.nibr.go.kr",
            content_sha256=(
                "75825e4c6eb4c82326adcdd4a1c48e0316b34fbd86492978e58cc538cc687e45"
            ),
            content_sha256_basis="downloaded_official_flora_page_086_jpeg_bytes",
        )
    )

    rows.append(
        _row(
            species="Citharexylum myrianthum",
            trait="inflorescence_display",
            value="raceme_spike_panicle",
            raw_value="flowers occur in raceme-like inflorescences",
            excerpt=(
                "Citharexylum myrianthum is an example of cryptic dioecy in "
                "morphologically perfect flowers. Flowers are small, tubular, "
                "white-colored, and crepuscular, with a sweet, pleasant scent "
                "and occur in raceme-like inflorescences."
            ),
            quality="high",
            provider="Rocca and Sazima 2006 Flora via FAO AGRIS",
            url=(
                "https://agris.fao.org/search/en/providers/122535/records/"
                "65de55f70f3e94b9e5cdae7a"
            ),
            title=(
                "The dioecious, sphingophilous species Citharexylum "
                "myrianthum (Verbenaceae): Pollination and visitor diversity"
            ),
            citation=(
                "Rocca & Sazima (2006), Flora 201(6):440-450, "
                "DOI 10.1016/j.flora.2006.02.001"
            ),
            record_id="doi:10.1016/j.flora.2006.02.001:abstract:raceme-like",
            lineage="doi:10.1016/j.flora.2006.02.001",
            lineage_method="original_peer_reviewed_article_doi",
            source_tier="A",
            source_type=(
                "peer_reviewed_primary_article_abstract_official_fao_record"
            ),
            domain="agris.fao.org",
            content_sha256=(
                "478773f5836b53758c3de852f2ac4034a11969799e55c95fe27a6affdde06aa6"
            ),
            content_sha256_basis="retrieved_fao_agris_article_record_html_bytes",
        )
    )

    rows.append(
        _row(
            species="Hakea carinata",
            trait="autonomous_selfing_capacity",
            value="autonomous",
            raw_value=(
                "fruit set after pre-anthesis insect exclusion in 8/10 and "
                "6/10 wild plants; author concludes seed production by autogamy"
            ),
            excerpt=(
                "Flowers from which insects were excluded set fruit on eight "
                "of ten treated plants in Jenkins Scrub population and six of "
                "ten in Humbug Scrub (RT) population. [...] Hakea carinata can "
                "produce seed by autogamy."
            ),
            quality="high",
            provider="University of Adelaide doctoral thesis repository",
            url=(
                "https://digital.library.adelaide.edu.au/server/api/core/bitstreams/"
                "7e6f9e21-e413-48c2-8655-9396a292d389/content"
            ),
            title="Population Genetics of Hakea carinata F. Muell. ex Meissner",
            citation=(
                "Starr, G. (2001), University of Adelaide PhD thesis, "
                "Chapter 5, pp. 84-85"
            ),
            record_id="starr-2001:hakea-carinata:pp84-85:autogamy-test",
            lineage="hdl:2440/21691#chapter5-autogamy-test",
            lineage_method=(
                "persistent_handle_plus_distinct_original_chapter5_experiment"
            ),
            source_tier="A",
            source_type="doctoral_thesis_controlled_insect_exclusion_experiment",
            domain="digital.library.adelaide.edu.au",
            content_sha256=(
                "7f81f49a203ea024d7c399bd56ce80ae548e2abcfc95b7a10a0f823111a9bd02"
            ),
            content_sha256_basis="downloaded_institutional_repository_pdf_bytes",
            retrieved_at_utc="2026-08-12T09:53:30Z",
        )
    )

    rows.append(
        _row(
            species="Celtis sinensis",
            trait="floral_symmetry",
            value="actinomorphic",
            raw_value="Flower symmetry: Radial (Actinomorphic)",
            excerpt=(
                "Celtis sinensis Pers. Flower symmetry: Radial "
                "(Actinomorphic) 多方向對稱的."
            ),
            quality="medium",
            provider="Shiu Ying Hu Herbarium, Chinese University of Hong Kong",
            url=(
                "https://syhuherbarium.sls.cuhk.edu.hk/collections/"
                "factsheet-pro/celtis-sinensis/"
            ),
            title=(
                "Celtis sinensis Pers. | Shiu Ying Hu Herbarium Collections"
            ),
            citation=(
                "Shiu Ying Hu Herbarium Pro-Factsheet, Chinese University "
                "of Hong Kong (retrieved 2026-08-12)"
            ),
            record_id="cuhk-herbarium:celtis-sinensis:floral-symmetry-radial",
            lineage=(
                "provider_treatment:cuhk_shiu_ying_hu_herbarium:Celtis_sinensis"
            ),
            lineage_method="official_university_herbarium_species_factsheet",
            source_tier="B",
            source_type="university_herbarium_species_factsheet",
            domain="syhuherbarium.sls.cuhk.edu.hk",
            content_sha256=(
                "a155ec19fe0a221291917288b43bba13654f7c349e44d0a5910088535316501a"
            ),
            content_sha256_basis="retrieved_species_page_html_utf8_bytes",
            retrieved_at_utc="2026-08-12T09:53:30Z",
        )
    )

    rows.append(
        _row(
            species="Phoenix roebelenii",
            trait="inflorescence_display",
            value="raceme_spike_panicle",
            raw_value="Inflorescences: Type Panicle",
            excerpt=(
                "Phoenix roebelenii O'Brien. Family Arecaceae (Palmae). "
                "Inflorescences: Position Axillary, Interfoliar; Type Panicle."
            ),
            quality="medium",
            provider="Shiu Ying Hu Herbarium, Chinese University of Hong Kong",
            url=(
                "https://syhuherbarium.sls.cuhk.edu.hk/collections/"
                "factsheet-pro/phoenix-roebelenii/"
            ),
            title=(
                "Phoenix roebelenii O'Brien | Shiu Ying Hu Herbarium Collections"
            ),
            citation=(
                "Shiu Ying Hu Herbarium Pro-Factsheet, Chinese University "
                "of Hong Kong (retrieved 2026-08-12)"
            ),
            record_id=(
                "cuhk-herbarium:phoenix-roebelenii:inflorescence-type-panicle"
            ),
            lineage=(
                "provider_treatment:cuhk_shiu_ying_hu_herbarium:Phoenix_roebelenii"
            ),
            lineage_method="official_university_herbarium_species_factsheet",
            source_tier="C",
            source_type="university_herbarium_horticultural_species_factsheet",
            domain="syhuherbarium.sls.cuhk.edu.hk",
            content_sha256=(
                "1ad7fb1b1d22e4917e4289773214edc6919948e98788bf2a3b6c2747143144ca"
            ),
            content_sha256_basis="retrieved_species_page_html_utf8_bytes",
            retrieved_at_utc="2026-08-12T10:20:00Z",
            cultivar_status=(
                "species_level_horticultural_record_not_cultivar_limited"
            ),
        )
    )

    rows.append(
        _row(
            species="Gonystylus confusus",
            trait="flower_primary_color",
            value="yellow_orange",
            raw_value="Flower Colour(s): Yellow / Golden",
            excerpt=(
                "Gonystylus confusus Airy Shaw. Family Name: Thymelaeaceae. "
                "Floral (Angiosperm): Flower Colour(s) Yellow / Golden."
            ),
            quality="high",
            provider="Singapore National Parks Board Flora & Fauna Web",
            url="https://www.nparks.gov.sg/florafaunaweb/flora/4/2/4227",
            title="Gonystylus confusus | NParks Flora & Fauna Web",
            citation=(
                "National Parks Board Singapore, Flora & Fauna Web species "
                "record 4227 (updated 2025-03-26; retrieved 2026-08-12)"
            ),
            record_id="nparks:flora:4227:flower-colour",
            lineage="url:https://www.nparks.gov.sg/florafaunaweb/flora/4/2/4227",
            lineage_method="canonical_government_species_record_url",
            source_tier="A",
            source_type="government_species_trait_database",
            domain="nparks.gov.sg",
            content_sha256=(
                "95b5251a0261187e344a28648ca8de199cd80433610eda1e90b59ff17c3bca78"
            ),
            content_sha256_basis="retrieved_species_page_html_bytes",
            retrieved_at_utc="2026-08-12T10:20:00Z",
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

    turraea_excerpt = (
        "Turraea cadetii A.J.Scott. Famille Meliaceae. Inflorescences : "
        "cymes axillaires de 2 à 3 fleurs. Fleurs : Les pétales sont blancs. "
        "Cette espèce est endémique de la Réunion."
    )
    turraea_common = {
        "species": "Turraea cadetii",
        "provider": "CIRAD Arbres et arbustes de La Réunion",
        "url": (
            "https://arbres-reunion.cirad.fr/especes/meliaceae/"
            "turraea_cadetii_a_j_scott.html"
        ),
        "title": (
            "Arbres et arbustes de La Réunion - Turraea cadetii A.J.Scott / "
            "Meliaceae"
        ),
        "citation": (
            "CIRAD, Arbres et arbustes de La Réunion, species treatment for "
            "Turraea cadetii A.J.Scott"
        ),
        "lineage": "provider_treatment:cirad_reunion:Turraea_cadetii",
        "lineage_method": "institutional_species_treatment_single_page_lineage",
        "source_tier": "A",
        "source_type": "institutional_regional_flora_species_treatment",
        "domain": "arbres-reunion.cirad.fr",
        "content_sha256": (
            "338cd6013d1d7ac887bb9180bfed16d80c446b3c397564b60c2990487802ed67"
        ),
        "content_sha256_basis": "retrieved_cirad_species_page_http_response_bytes",
        "language": "fr",
        "retrieved_at_utc": "2026-08-12T10:45:00Z",
    }
    rows.extend(
        [
            _row(
                trait="inflorescence_display",
                value="umbel_corymb",
                raw_value="cymes axillaires de 2 à 3 fleurs",
                excerpt=turraea_excerpt,
                quality="high",
                record_id="cirad-reunion:turraea-cadetii:axillary-cymes",
                **turraea_common,
            ),
            _row(
                trait="flower_primary_color",
                value="white",
                raw_value="les pétales sont blancs",
                excerpt=turraea_excerpt,
                quality="high",
                record_id="cirad-reunion:turraea-cadetii:white-petals",
                **turraea_common,
            ),
        ]
    )

    mosiera_excerpt = (
        "Mosiera yamaniguensis Bisse ex Urquiola & Z. Acosta, sp. nov. "
        "Flowers on young branches, solitary, axillary, additional 2-4 flowers "
        "per internode extra-axillary. Petals white, 3 mm long and wide, with "
        "small glands and hairy margins."
    )
    mosiera_common = {
        "species": "Mosiera yamaniguensis",
        "provider": "Urquiola Cruz and Acosta Ramos 2008 Willdenowia",
        "url": (
            "https://www.bgbm.org/sites/default/files/documents/"
            "wi38-2Urquiola%2BAcosta.pdf"
        ),
        "title": "Five new species of Mosiera (Myrtaceae) from Cuba",
        "citation": (
            "Urquiola Cruz, A. J. & Acosta Ramos, Z. (2008), Willdenowia "
            "38:533-544, p. 543, DOI 10.3372/wi.38.38213"
        ),
        "lineage": "doi:10.3372/wi.38.38213",
        "lineage_method": "original_peer_reviewed_taxonomic_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_taxonomic_species_description",
        "domain": "bgbm.org",
        "content_sha256": (
            "5e43877d31e0ab0e91be558534b6ab583e4e0089cf06b5eadc82c557a87c0a2c"
        ),
        "content_sha256_basis": "downloaded_bgbm_original_article_pdf_bytes",
        "language": "en+la",
        "retrieved_at_utc": "2026-08-12T10:45:00Z",
    }
    rows.extend(
        [
            _row(
                trait="inflorescence_display",
                value="solitary",
                raw_value="flowers solitary, axillary",
                excerpt=mosiera_excerpt,
                quality="high",
                record_id="doi:10.3372/wi.38.38213:p543:solitary-flowers",
                **mosiera_common,
            ),
            _row(
                trait="flower_primary_color",
                value="white",
                raw_value="petals white",
                excerpt=mosiera_excerpt,
                quality="high",
                record_id="doi:10.3372/wi.38.38213:p543:white-petals",
                **mosiera_common,
            ),
        ]
    )

    rows.append(
        _row(
            species="Dichaetanthera rutenbergiana",
            trait="flower_primary_color",
            value="green_brown_inconspicuous|yellow_orange",
            raw_value="Flores lutei brunneique",
            excerpt=(
                "D. rutenbergiana Baill. ms. [...] Ambalita? Mai 1878, fl. fr. "
                "jun. Flores lutei brunneique. Descriptio a cl. Baillon "
                "communicata."
            ),
            quality="high",
            provider="Vatke 1885 Reliquiae Rutenbergianae VI",
            url=(
                "https://archive.org/details/abhandlungenhera09natu/"
                "page/116/mode/2up"
            ),
            title="Reliquiae Rutenbergianae VI (Botanik, Fortsetzung)",
            citation=(
                "Vatke, W. (1885), Reliquiae Rutenbergianae VI, "
                "Abhandlungen des Naturwissenschaftlichen Vereins zu Bremen "
                "9:115-138, p. 116"
            ),
            record_id="vatke-1885:p116:dichaetanthera-rutenbergiana:flower-colour",
            lineage="publication:vatke-1885:reliquiae-rutenbergianae-vi",
            lineage_method="original_taxonomic_species_description_publication",
            source_tier="A",
            source_type="original_taxonomic_species_description",
            domain="archive.org",
            content_sha256=(
                "bfdb1926b67126cc794c3e67dc6c9723385379faab3ec206e3026072d979b4b7"
            ),
            content_sha256_basis="downloaded_internet_archive_volume_pdf_bytes",
            language="la+de",
            retrieved_at_utc="2026-08-12T10:50:00Z",
        )
    )

    acropogon_excerpt = (
        "Flower color of this new species is intermediate between A. veillonii, "
        "which has yellow flowers with sometimes some red at the base of the "
        "calyx, and A. bullatus, which has red flowers with some yellow at the "
        "edges. Table 1: Flower (calyx) inside color - A. bullatus: Red with "
        "some yellow at the edges; A. mesophilus: Yellow with 2-4 "
        "grayish-purple stripes per lobe; A. veillonii: Yellow with some red at "
        "the base."
    )
    acropogon_common = {
        "provider": "Munzinger and Gâteblé 2017 Phytotaxa",
        "url": (
            "https://horizon.documentation.ird.fr/exl-doc/pleins_textes/"
            "divers17-07/010070247.pdf"
        ),
        "title": (
            "Novitates neocaledonicae VI: Acropogon mesophilus (Malvaceae, "
            "Sterculioideae), a rare and threatened new species from the mesic "
            "forest of New Caledonia"
        ),
        "citation": (
            "Munzinger, J. & Gâteblé, G. (2017), Phytotaxa 307(3):183-190, "
            "p. 188, DOI 10.11646/phytotaxa.307.3.2"
        ),
        "lineage": "doi:10.11646/phytotaxa.307.3.2",
        "lineage_method": "original_peer_reviewed_comparative_taxonomic_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_taxonomic_comparison",
        "domain": "horizon.documentation.ird.fr",
        "content_sha256": (
            "c83b2f236f4a96a2ef20299a8716f135d35ec30d4f4f4a981605042e934150fb"
        ),
        "content_sha256_basis": "downloaded_ird_original_article_pdf_bytes",
        "retrieved_at_utc": "2026-08-12T10:55:00Z",
    }
    for species, value, raw_value in (
        (
            "Acropogon bullatus",
            "red_pink|yellow_orange",
            "red with some yellow at the edges",
        ),
        (
            "Acropogon mesophilus",
            "blue_purple|yellow_orange",
            "yellow with 2-4 grayish-purple stripes per calyx lobe",
        ),
        (
            "Acropogon veillonii",
            "red_pink|yellow_orange",
            "yellow with sometimes some red at the base of the calyx",
        ),
    ):
        slug = species.casefold().replace(" ", "-")
        rows.append(
            _row(
                species=species,
                trait="flower_primary_color",
                value=value,
                raw_value=raw_value,
                excerpt=acropogon_excerpt,
                quality="high",
                record_id=f"doi:10.11646/phytotaxa.307.3.2:p188:{slug}:colour",
                **acropogon_common,
            )
        )

    galapagos_excerpt_template = (
        "Table 1, {species}: Autonomously self-pollinates: Yes; "
        "Self-compatible: Yes; References: {references}. The table footnote "
        "identifies the numbered original studies used for each species."
    )
    galapagos_common = {
        "provider": "Traveset et al. 2013 Annals of Botany",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3489146/",
        "title": (
            "Pollination patterns and plant breeding systems in the "
            "Galapagos: a review"
        ),
        "citation": (
            "Traveset et al. (2013), Annals of Botany 111:391-404, "
            "Table 1, DOI 10.1093/aob/mcs132"
        ),
        "lineage_method": (
            "peer_reviewed_compilation_carrier_linked_to_numbered_original_study"
        ),
        "source_tier": "B",
        "source_type": "peer_reviewed_species_level_breeding_system_compilation",
        "domain": "pmc.ncbi.nlm.nih.gov",
        "content_sha256": (
            "d94f9c9b6d56ed778f09a53eefc248a4d45374c4c169913c4a78ba9760b2f444"
        ),
        "content_sha256_basis": "downloaded_pmc_full_text_html_bytes",
        "retrieved_at_utc": "2026-08-12T11:16:00Z",
    }
    galapagos_rows = (
        (
            "Alternanthera echinocephala",
            "autonomous_selfing_capacity",
            "autonomous",
            "4, 6",
            "bsdb-original:mcmullen_1987",
        ),
        (
            "Bidens pilosa",
            "self_incompatibility",
            "SC",
            "4, 6",
            "bsdb-original:mcmullen_1987",
        ),
        (
            "Pectis tenuifolia",
            "autonomous_selfing_capacity",
            "autonomous",
            "12",
            "study:philipp-et-al-2006:galapagos-lava-desert-network",
        ),
        (
            "Scalesia helleri",
            "autonomous_selfing_capacity",
            "autonomous",
            "2",
            "source:rick-1966:galapagos-plant-animal-relations",
        ),
        (
            "Scalesia aspera",
            "autonomous_selfing_capacity",
            "autonomous",
            "6",
            "study:mcmullen-1990:galapagos-reproductive-biology",
        ),
        (
            "Scalesia baurii",
            "autonomous_selfing_capacity",
            "autonomous",
            "7",
            "study:mcmullen-naranjo-1994:scalesia-baurii-pollination",
        ),
        (
            "Tournefortia pubescens",
            "autonomous_selfing_capacity",
            "autonomous",
            "4, 6",
            "bsdb-original:mcmullen_1987",
        ),
        (
            "Tournefortia pubescens",
            "self_incompatibility",
            "SC",
            "4, 6",
            "bsdb-original:mcmullen_1987",
        ),
        (
            "Cyperus elegans",
            "autonomous_selfing_capacity",
            "autonomous",
            "4, 6",
            "bsdb-original:mcmullen_1987",
        ),
        (
            "Nolana galapagensis",
            "autonomous_selfing_capacity",
            "autonomous",
            "2",
            "source:rick-1966:galapagos-plant-animal-relations",
        ),
        (
            "Nolana galapagensis",
            "self_incompatibility",
            "SC",
            "2",
            "source:rick-1966:galapagos-plant-animal-relations",
        ),
        (
            "Paspalum conjugatum",
            "autonomous_selfing_capacity",
            "autonomous",
            "4, 6",
            "bsdb-original:mcmullen_1987",
        ),
        (
            "Lycium minimum",
            "autonomous_selfing_capacity",
            "autonomous",
            "2",
            "source:rick-1966:galapagos-plant-animal-relations",
        ),
    )
    for species, trait, value, references, lineage in galapagos_rows:
        slug = species.casefold().replace(" ", "-")
        rows.append(
            _row(
                species=species,
                trait=trait,
                value=value,
                raw_value="Yes",
                excerpt=galapagos_excerpt_template.format(
                    species=species, references=references
                ),
                quality="medium",
                record_id=(
                    f"doi:10.1093/aob/mcs132:table1:{slug}:{trait}"
                ),
                lineage=lineage,
                **galapagos_common,
            )
        )

    encyclia_selfing_excerpt = (
        "This is the first paper of a series that reports the results of over "
        "50 years of selfing (autofecundation) of species by Ruben In Orchids "
        "of Miami, Florida. In most cases usually only one plant was available "
        "of each species, therefore outcrossing was not possible. These "
        "included: Encyclia phoenicea (Lindl.) Neumann, Encyclia plicata "
        "(Lindl.) Schltr. Selﬁngs were made several times of each species. "
        "Accurate photographic records were kept."
    )
    encyclia_selfing_common = {
        "provider": "Sauleda 2017 Orquideologia",
        "url": (
            "https://sco.org.co/wp-content/uploads/2018/07/"
            "Orquideologia-34.pdf"
        ),
        "title": (
            "Artificial Self-pollination (Autofecundation) as a Taxonomic "
            "Tool - Encyclia tampensis (Lindl.) Small"
        ),
        "citation": "Sauleda, R. P. (2017), Orquideologia 34(2):181-197",
        "lineage": "study:sauleda-2017:orchid-species-selfing-series",
        "lineage_method": "same_50_year_living_collection_selfing_program",
        "source_tier": "B",
        "source_type": "orchid_society_journal_documented_artificial_selfing",
        "domain": "sco.org.co",
        "content_sha256": (
            "8d788c4b4701b3752c8f87cb67e5407573305c7754f5bbbeafdf405c326f8df8"
        ),
        "content_sha256_basis": "downloaded_society_journal_issue_pdf_bytes",
        "retrieved_at_utc": "2026-08-12T11:35:00Z",
        "cultivar_status": (
            "living_species_accession_statement_not_cultivar_or_hybrid_limited"
        ),
    }
    for species in ("Encyclia phoenicea", "Encyclia plicata"):
        slug = species.casefold().replace(" ", "-")
        rows.append(
            _row(
                species=species,
                trait="self_incompatibility",
                value="SC",
                raw_value="repeated artificial selfings of the species",
                excerpt=encyclia_selfing_excerpt,
                quality="medium",
                record_id=f"sauleda-2017:p194:{slug}:self-compatible",
                **encyclia_selfing_common,
            )
        )

    tampensis_excerpt = (
        "Experiments using pollinator-exclusion bags revealed that this species "
        "is not capable of self-pollination and requires a pollen vector for "
        "seed capsule development."
    )
    rows.append(
        _row(
            species="Encyclia tampensis",
            trait="autonomous_selfing_capacity",
            value="absent",
            raw_value="not capable of spontaneous self-pollination",
            excerpt=tampensis_excerpt,
            quality="high",
            provider="Ray et al. 2019 Florida Entomologist",
            url=(
                "https://www.ars.usda.gov/research/publications/publication/"
                "?seqNo115=349510"
            ),
            title=(
                "Aspects of the pollination biology of Encyclia tampensis, "
                "the commercially exploited butterfly orchid, and Prosthechea "
                "cochleata, the endangered clamshell orchid, in south Florida"
            ),
            citation=(
                "Ray et al. (2019), Florida Entomologist 102(1):154-160, "
                "DOI 10.1653/024.102.0125"
            ),
            record_id="doi:10.1653/024.102.0125:pollinator-exclusion:no-autonomy",
            lineage="doi:10.1653/024.102.0125",
            lineage_method="original_peer_reviewed_pollinator_exclusion_experiment",
            source_tier="A",
            source_type="peer_reviewed_primary_pollinator_exclusion_experiment",
            domain="ars.usda.gov",
            content_sha256=_text_sha256(tampensis_excerpt),
            content_sha256_basis="verified_usda_ars_technical_abstract_excerpt_utf8",
            retrieved_at_utc="2026-08-12T11:28:00Z",
        )
    )

    southern_ocean_common = {
        "provider": "Lord 2015 AoB PLANTS Southern Ocean Islands",
        "url": (
            "https://academic.oup.com/aobpla/article/doi/10.1093/"
            "aobpla/plv095/1800409"
        ),
        "title": (
            "Patterns in floral traits and plant breeding systems on "
            "Southern Ocean Islands"
        ),
        "citation": (
            "Lord, J. M. (2015), AoB PLANTS 7:plv095, Supporting "
            "Information Table S1, DOI 10.1093/aobpla/plv095"
        ),
        "lineage": "compilation:doi:10.1093/aobpla/plv095:table-s1",
        "lineage_method": (
            "single_conservative_compilation_lineage_pending_original_source_retrieval"
        ),
        "source_tier": "B",
        "source_type": "peer_reviewed_species_level_compatibility_compilation",
        "domain": "academic.oup.com",
        "content_sha256": (
            "91222afed21322c51cfef98d083d75e4719cc4198e739594d6894db368291fe8"
        ),
        "content_sha256_basis": "downloaded_official_supplementary_docx_bytes",
        "retrieved_at_utc": "2026-08-12T12:02:00Z",
    }
    southern_ocean_rows = (
        ("Hydrocotyle chamaemorus", "SC", "SC; 5"),
        ("Agoseris coronopifolia", "SC", "SC; 5"),
        ("Anaphalioides bellidioides", "SI", "SI; 32"),
        ("Gamochaeta antarctica", "SC", "SC; 5, 40"),
        ("Gamochaeta malvinensis", "SC", "SC; 5, 40"),
        ("Hieracium patagonicum", "SC", "SC; 5"),
        ("Calceolaria biflora", "SC", "SC; 5"),
        ("Colobanthus affinis", "mixed_or_variable", "SCp; 1, 18"),
        ("Drosera stenopetala", "SC", "SC; 1, 35"),
        ("Gentianella antarctica", "SC", "SC; 11, 14"),
        ("Gentianella antipoda", "SC", "SC; 14"),
        ("Lobelia pratiana", "SC", "SC; 5"),
        ("Plantago barbata", "SC", "SC; 5"),
        ("Veronica benthamii", "SC", "SC; 42"),
        ("Ranunculus biternatus", "SC", "SC; 2, 5, 10"),
        ("Ranunculus maclovianus", "SC", "SC; 5"),
        ("Ranunculus sericocephalus", "SC", "SC; 5"),
        ("Acaena antarctica", "SC", "SC; 5, 16"),
        ("Acaena lucida", "SC", "SC; 5"),
        ("Acaena magellanica", "SC", "SC; 2, 5, 9, 10"),
        ("Acaena tenera", "SC", "SC; 10"),
        ("Rubus geoides", "SC", "SC; 5"),
        ("Saxifraga magellanica", "SC", "SC; 5, 16"),
        ("Phyllachne colensoi", "SC", "SC; 32"),
        ("Viola magellanica", "SC", "SC; 11"),
        ("Juncus scheuchzerioides", "SC", "SC; 2, 5, 10"),
        ("Luzula alopecurus", "SC", "SC; 5"),
        ("Chiloglottis cornuta", "SC", "SC; 14"),
        ("Deschampsia antarctica", "SC", "SC; 2, 5, 10"),
        ("Puccinellia macquariensis", "SC", "SC; 23"),
    )
    for species, value, source_cell in southern_ocean_rows:
        slug = species.casefold().replace(" ", "-")
        excerpt = (
            f"Supporting Information Table S1: {species}; compatibility and "
            f"reference-number cell: {source_cell}. The table defines SC as "
            "fully self-compatible, SCp as partially self-compatible, and SI "
            "as self-incompatible."
        )
        rows.append(
            _row(
                species=species,
                trait="self_incompatibility",
                value=value,
                raw_value=source_cell.split(";")[0],
                excerpt=excerpt,
                quality="medium",
                record_id=f"doi:10.1093/aobpla/plv095:table-s1:{slug}:compatibility",
                **southern_ocean_common,
            )
        )

    anderson_common = {
        "provider": "Anderson et al. 2001 American Journal of Botany",
        "url": (
            "https://ri.conicet.gov.ar/bitstream/handle/11336/38527/"
            "CONICET_Digital_Nro.2a109d53-78cd-4bee-83fb-c7cfd6b3faab_A.pdf"
            "?isAllowed=y&sequence=2"
        ),
        "title": (
            "Breeding system and pollination of selected plants endemic to "
            "Juan Fernandez Islands"
        ),
        "citation": (
            "Anderson et al. (2001), American Journal of Botany 88:220-233, "
            "DOI 10.2307/2657013"
        ),
        "lineage": "doi:10.2307/2657013",
        "lineage_method": "original_peer_reviewed_field_experiment_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_bagging_and_controlled_cross_experiment",
        "domain": "ri.conicet.gov.ar",
        "retrieved_at_utc": "2026-08-12T12:12:00Z",
    }
    anderson_rows = (
        (
            "Berberis corymbosa",
            "self_incompatibility",
            "SI",
            "fruit set 0% (n=6); pollen tubes 0% (n=8)",
            (
                "Table 1 reports self crosses with fruit set 0% (n=6) and "
                "pollen tubes 0% (n=8). In Berberis corymbosa pollen grains "
                "germinate on the stigmata, but pollen tubes do not grow "
                "beyond the stigmatic level. The authors classify it as SI."
            ),
            "p4:berberis-corymbosa:self-incompatible",
        ),
        (
            "Wahlenbergia berteroi",
            "autonomous_selfing_capacity",
            "autonomous",
            "wind-aided mechanism of autogamy",
            (
                "In Wahlenbergia berteroi there is a wind-aided mechanism of "
                "autogamy: pollen on the inner corolla throat is gathered by "
                "the stigmatic lobes when shaken by the ever-present wind. "
                "Table 1 reports self crosses 100% (n=5)."
            ),
            "p5:wahlenbergia-berteroi:wind-aided-autogamy",
        ),
        (
            "Wahlenbergia fernandeziana",
            "autonomous_selfing_capacity",
            "delayed",
            "spontaneous selfing may occur late",
            (
                "In Wahlenbergia fernandeziana spontaneous selfing may occur "
                "late in the lifetime of a flower because the stigma branches "
                "recurve almost 360 degrees, allowing autogamous pollen "
                "deposition; Table 1 reports self crosses 100% (n=5)."
            ),
            "p5:wahlenbergia-fernandeziana:delayed-selfing",
        ),
        (
            "Escallonia callcottiae",
            "autonomous_selfing_capacity",
            "autonomous",
            "100% fruit production from bagged flowers (N=20)",
            (
                "The stigma position and surface promote self-pollen "
                "deposition as soon as the flower opens. The authors observed "
                "100% fruit production from bagged flowers (N=20 flowers)."
            ),
            "p6:escallonia-callcottiae:bagged-autonomy",
        ),
        (
            "Escallonia callcottiae",
            "mating_system",
            "mixed_mating",
            "mixed breeding system with different degrees of xenogamy and autogamy",
            (
                "The authors classify Escallonia callcottiae as a facultative "
                "selfer and state that it has a mixed breeding system with "
                "different degrees of xenogamy and autogamy, including "
                "geitonogamy promoted by hummingbirds."
            ),
            "p6:escallonia-callcottiae:mixed-breeding-system",
        ),
    )
    for species, trait, value, raw_value, excerpt, record_id in anderson_rows:
        rows.append(
            _row(
                species=species,
                trait=trait,
                value=value,
                raw_value=raw_value,
                excerpt=excerpt,
                quality="high",
                record_id=f"doi:10.2307/2657013:{record_id}",
                content_sha256=_text_sha256(excerpt),
                content_sha256_basis="verified_conicet_repository_pdf_excerpt_utf8",
                **anderson_common,
            )
        )

    rows.extend(
        [
            _row(
                species="Pyrostria commersonii",
                trait="flower_primary_color",
                value="white|yellow_orange",
                raw_value="corolle blanc jaunatre",
                excerpt=(
                    "Etablissement : Ravine des cabris. Nom scientifique: "
                    "Pyrostria commersonii J.F. Gmel. Famille: Rubiaceae. "
                    "Fleur : Corolle blanc jaunatre."
                ),
                quality="medium",
                provider="Academie de La Reunion arboretum species sheets",
                url=(
                    "https://pedagogie.ac-reunion.fr/fileadmin/"
                    "ANNEXES-ACADEMIQUES/03-PEDAGOGIE/02-COLLEGE/"
                    "sciences-vie-terre/Fiches-peda/Arboretum/"
                    "RavineDesCabris.pdf"
                ),
                title="Arboretum - Ravine des Cabris",
                citation=(
                    "Academie de La Reunion, Sciences de la vie et de la Terre, "
                    "Arboretum Ravine des Cabris species sheet, p. 2"
                ),
                record_id=(
                    "academie-reunion:ravine-des-cabris:"
                    "pyrostria-commersonii:flower-colour"
                ),
                lineage=(
                    "provider_treatment:academie-reunion:"
                    "ravine-des-cabris:pyrostria-commersonii"
                ),
                lineage_method="official_education_species_sheet",
                source_tier="B",
                source_type="public_education_arboretum_species_sheet",
                domain="pedagogie.ac-reunion.fr",
                content_sha256=(
                    "29228904c032b6983664cae35c5de182d1e15057676f4f8c6169c7f7fb7102e4"
                ),
                content_sha256_basis="downloaded_official_pdf_bytes",
                language="fr",
                retrieved_at_utc="2026-08-12T12:01:49Z",
            ),
            _row(
                species="Quintinia acutifolia",
                trait="flower_primary_color",
                value="white",
                raw_value="showy white flowers",
                excerpt=(
                    "Quintinia acutifolia. Handsome small tree with long glossy "
                    "bronzy-green marbled leaves. In summer it has showy white "
                    "flowers which look great with the bold foliage. Genus "
                    "Quintinia; Species acutifolia; Cultivar [blank]."
                ),
                quality="medium",
                provider="Vibrant Earth New Zealand Plant Nursery",
                url=(
                    "https://www.vibrantearth.nz/catalogue/plantsdetail.php?"
                    "name=QUINTINIA+acutifolia&pid=1281"
                ),
                title="Quintinia acutifolia - Westland Quintinia",
                citation=(
                    "Vibrant Earth catalogue species treatment for Quintinia "
                    "acutifolia; non-cultivar morphology only"
                ),
                record_id="vibrant-earth:pid-1281:flower-colour",
                lineage="provider_treatment:vibrant-earth:pid-1281",
                lineage_method="canonical_nursery_species_page",
                source_tier="C",
                source_type="specialist_native_plant_nursery_species_page",
                domain="vibrantearth.nz",
                content_sha256=(
                    "ebd976ec71a9c55b00fcb3d8c5a0b9ec667d1fce95ecbf899ed99ff2ae566309"
                ),
                content_sha256_basis="verified_original_page_excerpt_utf8_bytes",
                retrieved_at_utc="2026-08-12T12:01:49Z",
                cultivar_status="species_level_horticultural_record_not_cultivar_limited",
            ),
            _row(
                species="Suregada lanceolata",
                trait="floral_symmetry",
                value="actinomorphic",
                raw_value="floral symmetry: actinomorphic",
                excerpt=(
                    "Tree Species Details: Family Name: Euphorbiaceae; Genus: "
                    "Suregada; Species: Suregada lanceolata; floral symmetry: "
                    "actinomorphic."
                ),
                quality="medium",
                provider="Bangladesh Forest Information System",
                url="https://bfis.bforest.gov.bd/nef/index.php/data/dataSpecies/40",
                title="Forest Emission Factor Database - Suregada lanceolata",
                citation=(
                    "Bangladesh Forest Information System, Forest Emission "
                    "Factor Database species record: Suregada lanceolata"
                ),
                record_id=(
                    "bfis:forest-emission-factor:"
                    "suregada-lanceolata:floral-symmetry"
                ),
                lineage=(
                    "provider_treatment:bangladesh-forest-information-system:"
                    "suregada-lanceolata"
                ),
                lineage_method="canonical_government_database_species_record",
                source_tier="A",
                source_type="government_forest_database_species_record",
                domain="bfis.bforest.gov.bd",
                content_sha256=(
                    "a0ea8e5b9481d4c66d75fd24cf045c43986da201e0720061158cf75244e79cf0"
                ),
                content_sha256_basis="verified_original_page_excerpt_utf8_bytes",
                retrieved_at_utc="2026-08-12T12:01:49Z",
            ),
        ]
    )
    rows.extend(_bfis_bulk_symmetry_rows())
    rows.extend(_india_flora_online_rows())
    rows.extend(_india_flora_online_synonym_rows())
    rows.extend(_prota_monograph_rows())
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
        "Actephila excelsa": "Phyllanthaceae",
        "Alchornea tiliifolia": "Euphorbiaceae",
        "Antidesma acidum": "Phyllanthaceae",
        "Antidesma ghaesembilla": "Phyllanthaceae",
        "Antidesma montanum": "Phyllanthaceae",
        "Antidesma roxburghii": "Phyllanthaceae",
        "Antidesma velutinum": "Phyllanthaceae",
        "Aporosa aurea": "Phyllanthaceae",
        "Baccaurea ramiflora": "Phyllanthaceae",
        "Barringtonia acutangula": "Lecythidaceae",
        "Breynia retusa": "Phyllanthaceae",
        "Bridelia retusa": "Phyllanthaceae",
        "Bridelia tomentosa": "Phyllanthaceae",
        "Careya arborea": "Lecythidaceae",
        "Casearia tomentosa": "Salicaceae",
        "Chaetocarpus castanocarpus": "Peraceae",
        "Cinnamomum iners": "Lauraceae",
        "Cleidion javanicum": "Euphorbiaceae",
        "Couroupita guianensis": "Lecythidaceae",
        "Croton aromaticus": "Euphorbiaceae",
        "Croton joufra": "Euphorbiaceae",
        "Croton tiglium": "Euphorbiaceae",
        "Endospermum chinense": "Euphorbiaceae",
        "Engelhardtia roxburghiana": "Juglandaceae",
        "Engelhardtia spicata": "Juglandaceae",
        "Epiprinus siletianus": "Euphorbiaceae",
        "Euphorbia antiquorum": "Euphorbiaceae",
        "Euphorbia cotinifolia": "Euphorbiaceae",
        "Euphorbia neriifolia": "Euphorbiaceae",
        "Euphorbia tirucalli": "Euphorbiaceae",
        "Excoecaria oppositifolia": "Euphorbiaceae",
        "Falconeria insignis": "Euphorbiaceae",
        "Flacourtia jangomas": "Salicaceae",
        "Flueggea leucopyrus": "Phyllanthaceae",
        "Flueggea virosa": "Phyllanthaceae",
        "Gomphandra tetrandra": "Stemonuraceae",
        "Homonoia riparia": "Euphorbiaceae",
        "Jatropha multifida": "Euphorbiaceae",
        "Macaranga denticulata": "Euphorbiaceae",
        "Macaranga peltata": "Euphorbiaceae",
        "Mallotus nudiflorus": "Euphorbiaceae",
        "Mallotus philippensis": "Euphorbiaceae",
        "Mallotus repandus": "Euphorbiaceae",
        "Mallotus tetracoccus": "Euphorbiaceae",
        "Margaritaria indica": "Phyllanthaceae",
        "Ostodes paniculata": "Euphorbiaceae",
        "Shirakiopsis indica": "Euphorbiaceae",
        "Sophora wightii": "Fabaceae",
        "Triadica sebifera": "Euphorbiaceae",
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
        "Adenia cissampeloides": "Passifloraceae",
        "Tristaniopsis laurina": "Myrtaceae",
        "Polycarpaea corymbosa": "Caryophyllaceae",
        "Boronia muelleri": "Rutaceae",
        "Benstonea foetida": "Pandanaceae",
        "Pleioluma balansana": "Sapotaceae",
        "Kunzea ericoides": "Myrtaceae",
        "Corchorus olitorius": "Malvaceae",
        "Jacobaea aquatica": "Asteraceae",
        "Carpinus laxiflora": "Betulaceae",
        "Citharexylum myrianthum": "Verbenaceae",
        "Hakea carinata": "Proteaceae",
        "Celtis sinensis": "Cannabaceae",
        "Phoenix roebelenii": "Arecaceae",
        "Gonystylus confusus": "Thymelaeaceae",
        "Turraea cadetii": "Meliaceae",
        "Mosiera yamaniguensis": "Myrtaceae",
        "Dichaetanthera rutenbergiana": "Melastomataceae",
        "Acropogon bullatus": "Malvaceae",
        "Acropogon mesophilus": "Malvaceae",
        "Acropogon veillonii": "Malvaceae",
        "Alternanthera echinocephala": "Amaranthaceae",
        "Bidens pilosa": "Asteraceae",
        "Pectis tenuifolia": "Asteraceae",
        "Scalesia helleri": "Asteraceae",
        "Scalesia aspera": "Asteraceae",
        "Scalesia baurii": "Asteraceae",
        "Tournefortia pubescens": "Heliotropiaceae",
        "Cyperus elegans": "Cyperaceae",
        "Nolana galapagensis": "Solanaceae",
        "Paspalum conjugatum": "Poaceae",
        "Lycium minimum": "Solanaceae",
        "Encyclia phoenicea": "Orchidaceae",
        "Encyclia plicata": "Orchidaceae",
        "Encyclia tampensis": "Orchidaceae",
        "Hydrocotyle chamaemorus": "Araliaceae",
        "Agoseris coronopifolia": "Asteraceae",
        "Anaphalioides bellidioides": "Asteraceae",
        "Gamochaeta antarctica": "Asteraceae",
        "Gamochaeta malvinensis": "Asteraceae",
        "Hieracium patagonicum": "Asteraceae",
        "Calceolaria biflora": "Calceolariaceae",
        "Colobanthus affinis": "Caryophyllaceae",
        "Drosera stenopetala": "Droseraceae",
        "Gentianella antarctica": "Gentianaceae",
        "Gentianella antipoda": "Gentianaceae",
        "Lobelia pratiana": "Campanulaceae",
        "Plantago barbata": "Plantaginaceae",
        "Veronica benthamii": "Plantaginaceae",
        "Ranunculus biternatus": "Ranunculaceae",
        "Ranunculus maclovianus": "Ranunculaceae",
        "Ranunculus sericocephalus": "Ranunculaceae",
        "Acaena antarctica": "Rosaceae",
        "Acaena lucida": "Rosaceae",
        "Acaena magellanica": "Rosaceae",
        "Acaena tenera": "Rosaceae",
        "Rubus geoides": "Rosaceae",
        "Saxifraga magellanica": "Saxifragaceae",
        "Phyllachne colensoi": "Stylidiaceae",
        "Viola magellanica": "Violaceae",
        "Juncus scheuchzerioides": "Juncaceae",
        "Luzula alopecurus": "Juncaceae",
        "Chiloglottis cornuta": "Orchidaceae",
        "Deschampsia antarctica": "Poaceae",
        "Puccinellia macquariensis": "Poaceae",
        "Berberis corymbosa": "Berberidaceae",
        "Wahlenbergia berteroi": "Campanulaceae",
        "Wahlenbergia fernandeziana": "Campanulaceae",
        "Escallonia callcottiae": "Escalloniaceae",
        "Pyrostria commersonii": "Rubiaceae",
        "Quintinia acutifolia": "Paracryphiaceae",
        "Suregada lanceolata": "Euphorbiaceae",
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

    synonym_reviewed = pd.read_csv(
        INDIA_FLORA_SYNONYM_REVIEWED_PATH, dtype=str
    ).fillna("")
    synonym_family = synonym_reviewed.set_index("accepted_species")["family"].to_dict()
    synonym_missing = sorted(set(synonym_family) - set(master_family))
    if synonym_missing:
        raise ValueError(
            f"IISc synonym species absent from target master: {synonym_missing}"
        )
    synonym_conflicts = {
        species: (family, master_family[species])
        for species, family in synonym_family.items()
        if master_family[species] != family
    }
    if synonym_conflicts:
        raise ValueError(f"IISc synonym family conflicts: {synonym_conflicts}")

    prota_reviewed = pd.read_csv(PROTA_REVIEWED_PATH, dtype=str).fillna("")
    prota_family = prota_reviewed.set_index("accepted_species")["family"].to_dict()
    prota_missing = sorted(set(prota_family) - set(master_family))
    if prota_missing:
        raise ValueError(f"PROTA species absent from target master: {prota_missing}")
    prota_conflicts = {
        species: (family, master_family[species])
        for species, family in prota_family.items()
        if master_family[species] != family
    }
    if prota_conflicts:
        raise ValueError(f"PROTA family conflicts: {prota_conflicts}")

    evidence = pd.DataFrame(reviewed_rows(), columns=EVIDENCE_COLUMNS).fillna("")
    evidence = evidence.sort_values(
        ["accepted_species", "trait_name", "candidate_id"]
    ).reset_index(drop=True)
    if len(evidence) != 922:
        raise ValueError(f"expected 922 reviewed trait rows, found {len(evidence)}")
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("wave-2 candidate IDs are not unique")
    audit = _review_audit(evidence)

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "rule_unlock_wave2_evidence_20260812.csv"
    audit_path = output_dir / "rule_unlock_wave2_manual_audit_20260812.csv"
    evidence.to_csv(evidence_path, index=False, lineterminator="\n")
    audit.to_csv(audit_path, index=False, lineterminator="\n")
    outputs = [
        evidence_path,
        audit_path,
        INDIA_FLORA_REVIEWED_PATH,
        INDIA_FLORA_FULL_AUDIT_PATH,
        INDIA_FLORA_SYNONYM_REVIEWED_PATH,
        INDIA_FLORA_SYNONYM_FULL_AUDIT_PATH,
        PROTA_REVIEWED_PATH,
        PROTA_FULL_AUDIT_PATH,
    ]

    combined_evidence: pd.DataFrame | None = None
    combined_audit: pd.DataFrame | None = None
    if prior_curated_evidence_csv is not None:
        if prior_curated_audit_csv is None:
            raise ValueError("prior curated audit is required with prior evidence")
        prior_evidence = pd.read_csv(prior_curated_evidence_csv, dtype=str).fillna("")
        prior_audit = pd.read_csv(prior_curated_audit_csv, dtype=str).fillna("")
        owned = set(evidence["candidate_id"])
        prior_owned = prior_evidence["source_group"].eq(SOURCE_GROUP)
        prior_owned_ids = set(
            prior_evidence.loc[prior_owned, "candidate_id"].astype(str)
        )
        combined_evidence = pd.concat(
            [prior_evidence.loc[~prior_owned], evidence],
            ignore_index=True,
        )
        combined_audit = pd.concat(
            [
                prior_audit.loc[
                    ~prior_audit["candidate_id"].astype(str).isin(prior_owned_ids | owned)
                ],
                audit,
            ],
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
        "iisc_full_candidate_audit": {
            "reviewed": 262,
            "accepted_correct": 250,
            "precision": 250 / 262,
            "cultivar_contamination_rate": 2 / 262,
        },
        "iisc_two_backbone_synonym_audit": {
            "reviewed": 71,
            "accepted_correct": 70,
            "precision": 70 / 71,
            "cultivar_contamination_rate": 0.0,
            "identity_contract": "exact_species_family_agreement_wfo_june_2026_and_gbif",
        },
        "prota_full_candidate_audit": {
            "reviewed": 469,
            "accepted_correct": 448,
            "precision": 448 / 469,
            "cultivar_contamination_rate": 0.0,
            "exact_master_species_pages_fetched": 1287,
            "non_redirect_pages_with_descriptions": 699,
            "mediawiki_api_calls": 52,
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
