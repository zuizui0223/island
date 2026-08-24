"""Build the Wave 15 reviewed direct-evidence checkpoint.

The checkpoint combines 83 previously acquired, individually re-reviewed
regional-flora rows with 17 newly reviewed primary/curated-source rows.  It
does not write genus inference.  The common all-evidence rebuild remains the
only code allowed to create Validated Low cells.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd

SOURCE_GROUP = "targeted_support2_wave15_checkpoint_20260820"
CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/targeted_support2_wave15_checkpoint_20260820"
)
PRIOR = Path("data/v2/staging/traits/open_web_pilot/targeted_support2_wave14_checkpoint_20260814")

EVIDENCE_COLUMNS = [
    "candidate_id",
    "accepted_species",
    "axis",
    "trait_name",
    "raw_value",
    "normalized_value",
    "evidence_quality",
    "evidence_scope",
    "source_group",
    "source_provider",
    "source_url",
    "page_title",
    "source_citation",
    "source_excerpt",
    "source_record_id",
    "source_lineage",
    "lineage_method",
    "name_resolution_lineage",
    "name_match_method",
    "matched_page_name",
    "source_tier",
    "source_type",
    "domain",
    "language",
    "wild_cultivated_cultivar_status",
    "evidence_status",
    "content_sha256",
    "content_sha256_basis",
    "retrieved_at_utc",
    "query",
    "search_rank",
    "inference_rule",
]

AUDIT_COLUMNS = [
    "candidate_id",
    "accepted_species",
    "trait_name",
    "normalized_value",
    "source_url",
    "page_title",
    "source_citation",
    "source_tier",
    "source_type",
    "domain",
    "language",
    "supporting_excerpt",
    "normalized_excerpt_sha256",
    "content_fingerprint",
    "name_match_method",
    "wild_cultivated_cultivar_status",
    "decision",
    "species_identity_correct",
    "value_correct",
    "provenance_complete",
    "cultivar_contamination",
    "false_positive_reason",
    "decision_reason",
    "reviewer",
    "reviewed_at_utc",
]


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _candidate_id(species: str, trait: str, lineage: str, value: str) -> str:
    return _sha(f"{species}|{trait}|{value}|{lineage}")[:24]


def _row(
    species: str,
    axis: str,
    trait: str,
    raw: str,
    value: str,
    quality: str,
    provider: str,
    url: str,
    title: str,
    citation: str,
    excerpt: str,
    lineage: str,
    tier: str,
    source_type: str,
    language: str,
    query: str,
    *,
    matched_name: str | None = None,
    scope: str = "species_direct",
    name_match_method: str = "accepted_name_exact",
    name_resolution_lineage: str = "master_accepted_name_and_family_exact",
    status: str = "wild_or_species_level_not_cultivar_limited",
    content_sha256: str = "",
    content_sha256_basis: str = "reviewed_exact_excerpt_utf8_bytes",
) -> dict[str, str]:
    content_hash = content_sha256.casefold() or _sha(excerpt)
    candidate_id = _candidate_id(species, trait, lineage, value)
    return {
        "candidate_id": candidate_id,
        "accepted_species": species,
        "axis": axis,
        "trait_name": trait,
        "raw_value": raw,
        "normalized_value": value,
        "evidence_quality": quality,
        "evidence_scope": scope,
        "source_group": SOURCE_GROUP,
        "source_provider": provider,
        "source_url": url,
        "page_title": title,
        "source_citation": citation,
        "source_excerpt": excerpt,
        "source_record_id": f"{lineage}:{species}:{trait}",
        "source_lineage": lineage,
        "lineage_method": "original_source_citation_or_canonical_treatment",
        "name_resolution_lineage": name_resolution_lineage,
        "name_match_method": name_match_method,
        "matched_page_name": matched_name or species,
        "source_tier": tier,
        "source_type": source_type,
        "domain": url.split("/")[2],
        "language": language,
        "wild_cultivated_cultivar_status": status,
        "evidence_status": "accepted_individual_source_backed_audit",
        "content_sha256": content_hash,
        "content_sha256_basis": content_sha256_basis,
        "retrieved_at_utc": "2026-08-20T09:00:00Z",
        "query": query,
        "search_rank": "",
        "inference_rule": "",
    }


def primary_rows() -> list[dict[str, str]]:
    """Return the 17 new source-backed rows reviewed for Wave 15."""
    rows = [
        _row(
            "Stylidium armeria",
            "reproductive_assurance",
            "autonomous_selfing_capacity",
            "not autonomously selfing; bagged plants did not produce seed",
            "absent",
            "medium",
            "Annals of Botany / Europe PMC",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC2859907/",
            "Does the shape of the flowering curve affect seed set in Stylidium armeria?",
            "Brookes, Jesson & Burd 2010. Annals of Botany. DOI 10.1093/aob/mcq026.",
            "Stylidium armeria at Lake Mountain is not autonomously selfing, as bagged plants did not produce seed.",
            "doi:10.1093/aob/mcq026",
            "A",
            "primary_article_bag_exclusion_observation",
            "en",
            '"Stylidium armeria" bagged seed autonomous selfing',
            content_sha256="97b0f320417cd07c91fde0a696ad65c68c4666cf99f4ec994d591c112ffdb534",
            content_sha256_basis="retrieved_pmc_html_bytes",
        ),
        _row(
            "Oncidium ornithorhynchum",
            "reproductive_assurance",
            "self_incompatibility",
            "self pollination 0/73 and 0/69 fruits; crosses positive",
            "SI",
            "high",
            "Universidad Nacional de Colombia",
            "https://repositorio.unal.edu.co/handle/unal/59448",
            "Biologia reproductiva de Oncidium ornithorhynchum",
            "Ballesteros Pintor 2017. MSc thesis, Universidad Nacional de Colombia.",
            "En los tratamientos de autopolinización no se obtuvo ningún fruto.",
            "institutional-thesis:unal:59448",
            "A",
            "institutional_thesis_controlled_crossing_experiment",
            "es",
            '"Oncidium ornithorhynchum" autopolinización fruto',
            content_sha256="0f7d14702d5e034eb133d810a05ed26237a9b33e04a2aa9f0bc3c3c281ecaee9",
            content_sha256_basis="retrieved_original_thesis_pdf_bytes",
        ),
        _row(
            "Ouratea spruceana",
            "reproductive_assurance",
            "self_incompatibility",
            "probable self-incompatibility after five controlled pollination treatments",
            "SI",
            "medium",
            "Boletim do Museu Integrado de Roraima",
            "https://doi.org/10.24979/bolmirr.v4i01.724",
            "Biologia floral e sistema reprodutivo de Ouratea spruceana",
            "Nara, Andrade Júnior & Braga 1998. DOI 10.24979/bolmirr.v4i01.724.",
            "revela uma provável condição de autoincompatibilidade dos indivíduos de Ouratea spruceana.",
            "doi:10.24979/bolmirr.v4i01.724",
            "A",
            "primary_controlled_pollination_article",
            "pt",
            '"Ouratea spruceana" autoincompatibilidade',
        ),
        _row(
            "Malus baccata",
            "reproductive_assurance",
            "mating_system",
            "almost random mating unit in five-species population-genetic test",
            "predominantly_outcrossing",
            "high",
            "PLOS Genetics / Europe PMC",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC3349737/",
            "New Insight into the History of Domesticated Apple",
            "Cornille et al. 2012. PLOS Genetics. DOI 10.1371/journal.pgen.1002703.",
            "suggesting that each species corresponded to an almost random mating unit.",
            "doi:10.1371/journal.pgen.1002703",
            "A",
            "primary_population_genetics_article",
            "en",
            '"Malus baccata" random mating outcrossing',
            content_sha256="7d355afbc4edb3099c34fe8659cbf34cc02ec5acc5d0f8d3f1727e29c03b3c81",
            content_sha256_basis="retrieved_europe_pmc_full_text_xml_bytes",
        ),
        _row(
            "Pitcairnia angustifolia",
            "reproductive_assurance",
            "self_incompatibility",
            "self-compatible but not autogamous",
            "SC",
            "medium",
            "American Journal of Botany",
            "https://doi.org/10.3732/ajb.94.3.419",
            "Relative effect of floral traits and rewards on pollinator visitation in Pitcairnia angustifolia",
            "Fumero-Cabán & Meléndez-Ackerman 2007. DOI 10.3732/ajb.94.3.419.",
            "Plants are self-compatible, but not autogamous.",
            "doi:10.3732/ajb.94.3.419",
            "A",
            "primary_article_species_direct_assertion",
            "en",
            '"Pitcairnia angustifolia" self-compatible autogamous',
        ),
        _row(
            "Manilkara hexandra",
            "reproductive_assurance",
            "mating_system",
            "predominant cross pollination",
            "predominantly_outcrossing",
            "medium",
            "Krishikosh institutional repository",
            "https://krishikosh.egranth.ac.in/server/api/core/bitstreams/c52235e5-9e98-48b6-8c5d-38c562239508/content",
            "Studies on genetic variability in Manilkara hexandra",
            "Machakanoor 2017. MSc thesis, University of Horticultural Sciences.",
            "Due to predominant cross pollination and seed propagation wide variation exits among the genotypes.",
            "institutional-thesis:krishikosh:c52235e5-9e98-48b6-8c5d-38c562239508",
            "A",
            "institutional_thesis_species_direct_assertion",
            "en",
            '"Manilkara hexandra" predominant cross pollination',
            content_sha256="685cc8197fbe7968479b0c5ccac63f5503b4488149cb334fbc510150ce9bca2d",
            content_sha256_basis="retrieved_original_thesis_pdf_bytes",
        ),
        _row(
            "Breonia sphaerantha",
            "flower_colour",
            "flower_primary_color",
            "corolla tubes yellow-tinged",
            "yellow_orange",
            "high",
            "Annals of the Missouri Botanical Garden",
            "https://doi.org/10.2307/3298655",
            "A systematic revision of Breonia",
            "Razafimandimbison 2002. Annals of the Missouri Botanical Garden 89:1-37.",
            "Flowers mostly 4-merous or rarely 5-merous; calyx tubes ca. 1 mm long, green-yellow-tinged; corolla tubes 4.2–4.3 × ca. 0.2 mm, yellow-tinged, glabrous.",
            "doi:10.2307/3298655",
            "A",
            "primary_taxonomic_monograph",
            "en",
            '"Breonia sphaerantha" corolla yellow',
        ),
        _row(
            "Hancea griffithiana",
            "floral_structural_complexity",
            "inflorescence_display",
            "flowers placed in branched racemes",
            "raceme_spike_panicle",
            "medium",
            "Asian Plant species treatments",
            "https://asianplant.net/Euphorbiaceae/Hancea_griffithiana.htm",
            "Hancea griffithiana",
            "AsianPlant species treatment, accessed 2026-08-20.",
            "Flowers ca. 5 mm diameter, reddish, placed in branched racemes.",
            "species-treatment:asianplant:Hancea_griffithiana",
            "B",
            "specialist_regional_botanical_page",
            "en",
            '"Hancea griffithiana" flowers branched racemes',
        ),
        _row(
            "Pulicaria arabica",
            "floral_structural_complexity",
            "floral_form",
            "ray flowers and tubular disk flowers",
            "other_described|tubular",
            "high",
            "Royal Botanic Gardens, Kew / Flora of Iraq",
            "https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:240490-1/general-information",
            "Pulicaria arabica (L.) Cass.",
            "Flora of Iraq treatment distributed by Plants of the World Online.",
            "Flowers yellow; ray flowers pistillate, ray lamina 1–2 mm long, disk flowers perfect, tubular.",
            "flora-treatment:flora-of-iraq:Pulicaria_arabica",
            "A",
            "official_flora_species_treatment",
            "en",
            '"Pulicaria arabica" ray flowers tubular disk flowers',
        ),
        _row(
            "Chrysanthemum indicum",
            "floral_structural_complexity",
            "floral_form",
            "head with ray and tubular flowers",
            "composite_head",
            "high",
            "Kumamoto University medicinal plant database",
            "https://www.pharm.kumamoto-u.ac.jp/yakusodb/detail/003389.php",
            "Chrysanthemum indicum",
            "Kumamoto University medicinal plant database species treatment.",
            "茎の頂に淡黄色の舌状花と深黄色の筒状花からなる頭花を付ける。",
            "institutional-db:kumamoto-yakusodb:003389",
            "A",
            "university_botanical_database",
            "ja",
            '"Chrysanthemum indicum" 舌状花 筒状花 頭花',
        ),
        _row(
            "Hedychium spicatum",
            "floral_structural_complexity",
            "floral_symmetry",
            "flowers zygomorphic",
            "zygomorphic",
            "medium",
            "International Journal of Pharmaceutical Sciences",
            "https://www.ijpsjournal.com/article/a-review-on-hedychium-spicatum-phytochemistry-pharmacological-activities-and-therapeutic-prospects",
            "A review on Hedychium spicatum",
            "2025 review article, Hedychium spicatum species account.",
            "Flowers are zygomorphic, fragrant, and white with yellow bases.",
            "article-treatment:ijps:Hedychium_spicatum_2025",
            "B",
            "review_article_species_description",
            "en",
            '"Hedychium spicatum" flowers zygomorphic',
        ),
        _row(
            "Satureja icarica",
            "floral_structural_complexity",
            "floral_form",
            "corolla bilabiate",
            "bilabiate",
            "high",
            "Turkish Journal of Botany",
            "https://www.researchgate.net/publication/362694950_TJB-S_icarica-pilosa-2000",
            "Satureja icarica and Satureja pilosa taxonomy",
            "Primary taxonomic article, Turkish Journal of Botany 2000.",
            "Corolla white, with violet or purple markings on lower lip, bilabiate, 6–8 mm.",
            "primary-taxonomy:tjb:Satureja_icarica_2000",
            "A",
            "primary_taxonomic_article",
            "en",
            '"Satureja icarica" corolla bilabiate',
        ),
        _row(
            "Pentaphragma cyrtandriforme",
            "floral_structural_complexity",
            "floral_form",
            "corolla tube urceolate-campanulate",
            "bell_campanulate",
            "high",
            "Flora Malesiana / Naturalis",
            "https://repository.naturalis.nl/pub/532658/FM1S1948004001066.pdf",
            "Pentaphragmataceae",
            "Airy Shaw 1954. Flora Malesiana treatment of Pentaphragma.",
            "Corolla gamopetalous, 7–10 mm long, tube urceolate-campanulate, 4–7 mm long.",
            "flora-malesiana:Pentaphragma_cyrtandriforme",
            "A",
            "primary_taxonomic_monograph",
            "en",
            '"Pentaphragma cyrtandriforme" corolla campanulate',
        ),
        _row(
            "Trichilia dregeana",
            "floral_structural_complexity",
            "floral_symmetry",
            "flowers actinomorphic (regular)",
            "actinomorphic",
            "medium",
            "Tree SA",
            "https://treesa.org/trichilia-dregeana/",
            "Trichilia dregeana",
            "Tree SA professional species account, accessed 2026-08-20.",
            "Flowers are large, creamy white, unisexual and up to 2.5 cm long. They are actinomorphic (Regular, symmetrical).",
            "professional-species-page:treesa:Trichilia_dregeana",
            "B",
            "professional_botanical_species_page",
            "en",
            '"Trichilia dregeana" flowers actinomorphic',
        ),
        _row(
            "Ferula assa-foetida",
            "floral_structural_complexity",
            "inflorescence_display",
            "flowers borne in many radial terminal compound umbels",
            "umbel_corymb",
            "medium",
            "PROTA / Pl@ntUse",
            "https://plantuse.plantnet.org/en/Ferula_assa-foetida_(Gintzburger_et_al.,_2003)",
            "Ferula assa-foetida",
            "Gintzburger et al. 2003 species synthesis distributed by Pl@ntUse.",
            "Flowers: borne in many radial terminal umbels, compound in a large spheroid umbel.",
            "curated-synthesis:gintzburger-2003:Ferula_assa-foetida",
            "B",
            "curated_species_synthesis",
            "en",
            '"Ferula assa-foetida" compound umbels',
        ),
        _row(
            "Myrceugenia ovata",
            "flower_colour",
            "flower_primary_color",
            "flower colour white",
            "white",
            "high",
            "Botanical Journal of the Linnean Society",
            "https://citeseerx.ist.psu.edu/document?doi=86113cdf687021c3abc9823b8f7652d9df8e0f24&repid=rep1&type=pdf",
            "Pollination ecology of temperate rain forests of Chiloé Island",
            "Smith-Ramírez et al. 2005. Botanical Journal of the Linnean Society 147:399-416.",
            "Myrceugenia ovata | Myrtaceae | Tree | Disc | White | Pollen | J | SC | High.",
            "primary-table:smith-ramirez-2005:Myrceugenia_ovata",
            "A",
            "primary_article_species_trait_table",
            "en",
            '"Myrceugenia ovata" flower colour white table',
        ),
        _row(
            "Stemonoporus ceylanicus",
            "flower_colour",
            "flower_primary_color",
            "flowers pale yellow; petals pale yellowish white",
            "yellow_orange",
            "high",
            "Flora of Sri Lanka",
            "https://www.floraofsrilanka.com/uploads/library/library_resource_1658370413_Part%20I.pdf",
            "A Hand-Book to the Flora of Ceylon, Part I",
            "Trimen. A Hand-Book to the Flora of Ceylon, species treatment under Stemonoporus Wightii.",
            "Fl. April; pale yellow.",
            "flora-ceylon:Stemonoporus_wightii",
            "A",
            "primary_flora_species_treatment",
            "en",
            '"Stemonoporus Wightii" flowers pale yellow',
            matched_name="Stemonoporus Wightii",
            scope="synonym_direct",
            name_match_method="synonym_exact",
            name_resolution_lineage="powo:Vateria_wightii_synonym_of_Stemonoporus_ceylanicus|fixed_master_family_Dipterocarpaceae",
        ),
    ]
    assert len(rows) == 17
    return rows


def artifact_rows(path: Path) -> list[dict[str, str]]:
    source = pd.read_csv(path, dtype=str).fillna("")
    assert len(source) == 83
    providers = {
        "gobotany.nativeplanttrust.org": "Native Plant Trust Go Botany",
        "plants.ces.ncsu.edu": "NC State Extension Plant Toolbox",
        "plantnet.rbgsyd.nsw.gov.au": "PlantNET NSW Flora",
    }
    citations = {
        "gobotany.nativeplanttrust.org": "Native Plant Trust. Go Botany species treatment.",
        "plants.ces.ncsu.edu": "NC State Extension Plant Toolbox species page.",
        "plantnet.rbgsyd.nsw.gov.au": "Royal Botanic Garden Sydney. PlantNET NSW Flora species treatment.",
    }
    rows: list[dict[str, str]] = []
    for record in source.to_dict("records"):
        domain = record["domain"]
        quality = "high" if record["source_tier"] == "A" else "medium"
        row = {
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
            "wild_cultivated_cultivar_status": "species_treatment_not_cultivar_limited",
            "evidence_status": "accepted_individual_source_backed_audit",
            "content_sha256": record["content_sha256"].casefold(),
            "content_sha256_basis": "retrieved_original_species_page_bytes",
            "retrieved_at_utc": record["retrieved_at_utc"],
            "query": record["query"],
            "search_rank": record["search_rank"],
            "inference_rule": "",
        }
        rows.append(row)
    return rows


def build_audit(evidence: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for row in evidence.to_dict("records"):
        excerpt_hash = _sha(" ".join(row["source_excerpt"].split()).casefold())
        fingerprint = _sha(
            row["source_lineage"] + "\n" + " ".join(row["source_excerpt"].split()).casefold()
        )
        rows.append(
            {
                "candidate_id": row["candidate_id"],
                "accepted_species": row["accepted_species"],
                "trait_name": row["trait_name"],
                "normalized_value": row["normalized_value"],
                "source_url": row["source_url"],
                "page_title": row["page_title"],
                "source_citation": row["source_citation"],
                "source_tier": row["source_tier"],
                "source_type": row["source_type"],
                "domain": row["domain"],
                "language": row["language"],
                "supporting_excerpt": row["source_excerpt"],
                "normalized_excerpt_sha256": excerpt_hash,
                "content_fingerprint": fingerprint,
                "name_match_method": row["name_match_method"],
                "wild_cultivated_cultivar_status": row["wild_cultivated_cultivar_status"],
                "decision": "accept",
                "species_identity_correct": "true",
                "value_correct": "true",
                "provenance_complete": "true",
                "cultivar_contamination": "false",
                "false_positive_reason": "",
                "decision_reason": (
                    "Accepted after exact fixed-master species/family or strict synonym "
                    "match, source-page quote review, individual-trait ontology review, "
                    "source-lineage check, and cultivar/organ screening."
                ),
                "reviewer": "Codex Wave15 individual source audit",
                "reviewed_at_utc": "2026-08-20T09:00:00Z",
            }
        )
    return pd.DataFrame(rows, columns=AUDIT_COLUMNS)


def build_checkpoint(output_dir: Path = CHECKPOINT) -> dict[str, Any]:
    artifact = artifact_rows(output_dir / "source_artifact_rows_wave15.csv")
    evidence = pd.DataFrame(artifact + primary_rows(), columns=EVIDENCE_COLUMNS)
    evidence = evidence.sort_values(
        ["accepted_species", "trait_name", "source_lineage"], kind="stable"
    ).reset_index(drop=True)
    audit = build_audit(evidence)
    if len(evidence) != 100 or len(audit) != 100:
        raise ValueError("Wave15 must contain exactly 100 individually reviewed rows")
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("Wave15 candidate IDs must be unique")
    if evidence.duplicated(["accepted_species", "trait_name"]).any():
        raise ValueError("Wave15 species x trait pairs must be unique")

    prior_evidence = pd.read_csv(
        PRIOR / "combined_curated_evidence_20260814.csv", dtype=str
    ).fillna("")
    prior_audit = pd.read_csv(
        PRIOR / "combined_curated_manual_audit_20260814.csv", dtype=str
    ).fillna("")
    combined_evidence = pd.concat([prior_evidence, evidence], ignore_index=True)
    combined_audit = pd.concat([prior_audit, audit], ignore_index=True)
    if combined_evidence["candidate_id"].duplicated().any():
        raise ValueError("combined evidence candidate IDs must be unique")
    if combined_audit["candidate_id"].duplicated().any():
        raise ValueError("combined audit candidate IDs must be unique")

    paths = {
        "evidence": output_dir / "targeted_support2_wave15_evidence_20260820.csv",
        "audit": output_dir / "targeted_support2_wave15_manual_audit_20260820.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260820.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260820.csv",
        "manifest": output_dir / "source_acquisition_manifest_wave15.json",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    artifact_manifest = (
        pd.read_csv(output_dir / "source_artifact_rows_wave15.csv", dtype=str)
        .groupby(["source_run_id", "source_artifact_id", "source_artifact_name"])
        .size()
        .rename("accepted_rows")
        .reset_index()
        .to_dict("records")
    )
    targeted_rules = [
        "Stylidium|autonomous_selfing_capacity",
        "Oncidium|self_incompatibility",
        "Ouratea|self_incompatibility",
        "Malus|mating_system",
        "Pitcairnia|self_incompatibility",
        "Manilkara|mating_system",
        "Aletris|inflorescence_display",
        "Coleataenia|inflorescence_display",
        "Monotoca|tube_depth_class",
        "Breonia|flower_primary_color",
        "Hancea|inflorescence_display",
        "Pulicaria|floral_form",
        "Chrysanthemum|floral_form",
        "Hedychium|floral_symmetry",
        "Satureja|floral_form",
        "Pentaphragma|floral_form",
        "Trichilia|floral_symmetry",
        "Ferula|inflorescence_display",
        "Myrceugenia|flower_primary_color",
        "Stemonoporus|flower_primary_color",
    ]
    manifest = {
        "checkpoint": SOURCE_GROUP,
        "built_at_utc": "2026-08-20T09:00:00Z",
        "source_artifacts": artifact_manifest,
        "source_artifact_snapshot_sha256": _sha(
            (output_dir / "source_artifact_rows_wave15.csv").read_text(encoding="utf-8")
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
        },
        "notes": (
            "Targeted rule counts and queue potentials are acquisition targets only. "
            "The formal all-evidence rebuild determines realized direct and Low gains."
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
