"""Build the Wave 24 support-two third-species checkpoint.

Every accepted row is an exact-species, source-backed direct observation for
one trait.  The checkpoint never creates genus inference itself; the shared
trait-specific all-evidence rebuild is the only component allowed to create
Validated Low cells.
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
from island_v2.targeted_support2_wave15_checkpoint import _row as _wave15_row
from island_v2.targeted_support2_wave15_checkpoint import (
    build_audit as _wave15_build_audit,
)

SOURCE_GROUP = "targeted_support2_wave24_checkpoint_20260824"
CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "targeted_support2_wave24_checkpoint_20260824"
)
PRIOR = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "targeted_dioecy_wave23_checkpoint_20260821"
)
RETRIEVED_AT = "2026-08-24T06:00:00Z"
FIXED_INTEGRATED_BASELINE_RUN_ID = 31488685866
PRIOR_FORMAL_RUN_ID = 32702160934


def _row(*args: str, **kwargs: str) -> dict[str, str]:
    row = _wave15_row(*args, **kwargs)
    row["source_group"] = SOURCE_GROUP
    row["retrieved_at_utc"] = RETRIEVED_AT
    row["content_sha256_basis"] = "normalized_supporting_excerpt_utf8"
    row["content_sha256"] = _sha(" ".join(row["source_excerpt"].split()))
    return row


# species, axis, trait, raw, value, quality, provider, URL, title, citation,
# excerpt, lineage, tier, source type, language, query, potential cells
RECORDS: list[tuple[Any, ...]] = [
    (
        "Gnidia squarrosa", "floral_structural_complexity",
        "inflorescence_display", "rounded heads of 6-30 flowers",
        "composite_display", "high", "SANBI PlantZAfrica",
        "https://pza.sanbi.org/gnidia-squarrosa", "Gnidia squarrosa",
        "South African National Biodiversity Institute species treatment.",
        "The flowers are borne in rounded heads of 6 - 30 flowers at the tips of the slender branches.",
        "provider-treatment:sanbi-plantzafrica:gnidia-squarrosa", "A",
        "official_biodiversity_species_treatment", "en",
        '"Gnidia squarrosa" flowers heads', 18,
    ),
    (
        "Selenicereus undatus", "floral_structural_complexity",
        "flower_size_class", "flowers 25-35 cm long and 30 cm wide",
        "very_large", "high", "National Parks Board Singapore",
        "https://www.nparks.gov.sg/florafaunaweb/flora/1/4/1419",
        "Selenicereus undatus - Flora & Fauna Web",
        "NParks Singapore exact-species treatment.",
        "Fragrant, whitish to greenish-yellow flowers (25-35 cm long, 30 cm wide).",
        "provider-treatment:nparks:flora:1419", "A",
        "official_government_species_database", "en",
        '"Selenicereus undatus" flower size', 8,
    ),
    (
        "Apium annuum", "floral_structural_complexity",
        "inflorescence_display", "compound umbels", "umbel_corymb", "high",
        "Journal of the Adelaide Botanic Gardens",
        "https://data.environment.sa.gov.au/Content/Publications/JABG01P205_Short.pdf",
        "A new species of Apium L. (Umbelliferae)",
        "Short 1979. Journal of the Adelaide Botanic Gardens 1:205-209.",
        "All species of Apium have compound umbels. Individuals belonging to A. annuum and A. prostratum have either pedunculate or sessile compound umbels.",
        "journal-article:short-1979-apium-annuum", "A",
        "primary_taxonomic_article", "en", '"Apium annuum" compound umbels', 6,
    ),
    (
        "Roystonea regia", "floral_structural_complexity",
        "inflorescence_display", "many-branched panicle",
        "raceme_spike_panicle", "high", "USDA Forest Service",
        "https://research.fs.usda.gov/treesearch/43219", "Roystonea regia",
        "USDA Forest Service species synthesis, Treesearch publication 43219.",
        "The fragrant flowers are borne on a many-branched panicle.",
        "official-species-treatment:usda-fs:roystonea-regia", "A",
        "official_public_agency_species_synthesis", "en",
        '"Roystonea regia" flowers panicle', 7,
    ),
    (
        "Veronicastrum virginicum", "floral_structural_complexity",
        "floral_form", "corolla tube 2-3 times as long as lobes", "tubular",
        "medium", "Missouri Plants", "https://missouriplants.com/Veronicastrum_virginicum_page.html",
        "Veronicastrum virginicum", "Missouri Plants exact-species botanical treatment.",
        "Corollas 4-6 mm long, nearly actinomorphic, 4-lobed, the tube 2-3 times as long as the lobes.",
        "provider-treatment:missouriplants:veronicastrum-virginicum", "B",
        "specialist_regional_flora", "en", '"Veronicastrum virginicum" corolla tube', 7,
    ),
    (
        "Veronicastrum virginicum", "floral_structural_complexity",
        "flower_size_class", "corolla 4-6 mm long", "small", "medium",
        "Missouri Plants", "https://missouriplants.com/Veronicastrum_virginicum_page.html",
        "Veronicastrum virginicum", "Missouri Plants exact-species botanical treatment.",
        "Corollas 4-6 mm long, nearly actinomorphic, 4-lobed, the tube 2-3 times as long as the lobes.",
        "provider-treatment:missouriplants:veronicastrum-virginicum", "B",
        "specialist_regional_flora", "en", '"Veronicastrum virginicum" corolla size', 7,
    ),
    (
        "Pourthiaea benthamiana", "flower_colour", "flower_primary_color",
        "flowers white", "white", "medium", "Taiwan Native Plants",
        "https://kplant.atfr.org.tw/%E5%8F%B0%E7%81%A3%E7%9F%B3%E6%A5%A0/%E5%8F%B0%E7%81%A3%E7%9F%B3%E6%A5%A0.htm",
        "台灣石楠", "Taiwan specialist native-plant exact-species treatment.",
        "花白色，徑 0.7~0.9 公分，多數，呈頂生的聚繖或複繖房花序。",
        "provider-treatment:kplant:pourthiaea-benthamiana", "B",
        "specialist_regional_flora", "zh", '"Pourthiaea benthamiana" 花 白色', 8,
    ),
    (
        "Thottea tomentosa", "floral_structural_complexity",
        "flower_size_class", "perianth 5-8 by 4-7 mm", "small", "high",
        "Flora of Thailand / World Flora Online",
        "https://www.worldfloraonline.org/taxon/wfo-0000411251",
        "Thottea tomentosa", "Flora of Thailand treatment redistributed by WFO.",
        "Perianth cup-shaped, 5-8 by 4-7 mm, 3-lobed.",
        "flora-treatment:flora-of-thailand:thottea-tomentosa", "A",
        "official_flora_species_treatment", "en", '"Thottea tomentosa" perianth 5-8 mm', 11,
    ),
    (
        "Uncarina peltata", "flower_colour", "flower_primary_color",
        "yellow flowers", "yellow_orange", "medium", "San Diego Zoo / Cactus Society of New Mexico",
        "https://www.new-mexico.cactus-society.org/ToThePoint/TTP_Fall_2022.01_compressed-1.pdf",
        "Madagascar collection garden profile", "San Diego Zoo garden profile reprinted in To The Point, Fall 2022.",
        "Uncarina peltata - This was planted in 2015, with beautiful yellow flowers.",
        "institutional-garden-profile:sandiegozoo:uncarina-peltata", "C",
        "institutional_garden_species_profile", "en", '"Uncarina peltata" yellow flowers', 12,
    ),
    (
        "Sabal palmetto", "floral_structural_complexity",
        "inflorescence_display", "four-times-branched panicle",
        "raceme_spike_panicle", "medium", "Botany Brisbane",
        "https://www.botanybrisbane.com/plants/arecaceae/sabal/sabal-palmetto/",
        "Sabal palmetto", "Botany Brisbane exact-species botanical treatment.",
        "Inflorescences, between the leaves are a panicle that is branched 4 times.",
        "provider-treatment:botanybrisbane:sabal-palmetto", "B",
        "specialist_botanical_species_profile", "en", '"Sabal palmetto" inflorescence panicle', 9,
    ),
    (
        "Hypochaeris radicata", "floral_structural_complexity",
        "tube_depth_class", "corolla tube 5-6 mm long", "intermediate", "high",
        "Flora Zambesiaca / Plants of the World Online",
        "https://powo.science.kew.org/taxon/urn%3Alsid%3Aipni.org%3Anames%3A225575-1/general-information",
        "Hypochaeris radicata", "Flora Zambesiaca treatment redistributed by POWO.",
        "Corolla tube 5-6 mm long.",
        "flora-treatment:flora-zambesiaca:hypochaeris-radicata", "A",
        "official_flora_species_treatment", "en", '"Hypochaeris radicata" corolla tube', 8,
    ),
    (
        "Robiquetia bertholdii", "floral_structural_complexity",
        "inflorescence_display", "30-45-flowered racemose inflorescence",
        "raceme_spike_panicle", "medium", "OrchidRoots / OrchidBoard",
        "https://www.orchids.org/grexes/25546", "Robiquetia bertholdii",
        "Orchid species profile with exact scientific name and diagnostic description.",
        "Blooms on a lateral, 3.6 to 8 inch long, 30 to 45 flowered, racemose inflorescence with non-resupinate flowers.",
        "provider-treatment:orchids-org:robiquetia-bertholdii", "B",
        "specialist_taxon_profile", "en", '"Robiquetia bertholdii" racemose inflorescence', 15,
    ),
    (
        "Meistera aculeata", "floral_structural_complexity",
        "inflorescence_display", "oblong spike", "raceme_spike_panicle", "high",
        "Journal of Natural Fibers / Taylor & Francis",
        "https://www.tandfonline.com/doi/full/10.1080/15440478.2025.2568538",
        "Exploring the natural fibers of Meistera aculeata",
        "Peer-reviewed primary species study. DOI 10.1080/15440478.2025.2568538.",
        "M. aculeata is a perennial herb and its inflorescence manifests as an oblong spike.",
        "doi:10.1080/15440478.2025.2568538", "A", "primary_article_species_study", "en",
        '"Meistera aculeata" inflorescence spike', 15,
    ),
    (
        "Codiaeum luzonicum", "flower_colour", "flower_primary_color",
        "flowers white", "white", "high", "Biodiversity Heritage Library",
        "https://upload.wikimedia.org/wikipedia/commons/7/7b/Botanical_publications_of_E.D._Merrill_%28IA_botanicalpublica05merr%29.pdf",
        "Botanical publications of E. D. Merrill", "Merrill protologue of Codiaeum luzonicum.",
        "Codiaeum luzonicum Merrill, sp. nov. Flowers white, the buds globose, the pedicels 5 to 10 mm long.",
        "protologue:merrill:codiaeum-luzonicum", "A", "primary_taxonomic_protologue", "en",
        '"Codiaeum luzonicum" flowers white', 9,
    ),
    (
        "Trichanthecium parvifolium", "floral_structural_complexity",
        "inflorescence_display", "panicle ovate 1-3 cm", "raceme_spike_panicle", "high",
        "Flora of Tropical East Africa / Plants of the World Online",
        "https://powo.science.kew.org/taxon/urn%3Alsid%3Aipni.org%3Anames%3A77119541-1/general-information",
        "Trichanthecium parvifolium", "Flora of Tropical East Africa treatment redistributed by POWO.",
        "Panicle ovate, 1-3 cm long.",
        "flora-treatment:ftea:trichanthecium-parvifolium", "A",
        "official_flora_species_treatment", "en", '"Trichanthecium parvifolium" panicle', 8,
    ),
    (
        "Gymnanthes borneensis", "floral_structural_complexity",
        "inflorescence_display", "axillary groups of racemes", "raceme_spike_panicle", "medium",
        "Asian Plant", "https://www.asianplant.net/Euphorbiaceae/Gymnanthes_borneensis.htm",
        "Gymnanthes borneensis", "Asian Plant exact-species regional treatment.",
        "Flowers ca. 1 mm diameter, yellowish-red, placed in axillary groups of racemes.",
        "provider-treatment:asianplant:gymnanthes-borneensis", "B",
        "specialist_regional_flora", "en", '"Gymnanthes borneensis" flowers racemes', 9,
    ),
    (
        "Bahiopsis chenopodina", "flower_colour", "flower_primary_color",
        "bright yellow ray and disk florets", "yellow_orange", "high",
        "University of California Riverside",
        "https://ezcurralab.ucr.edu/sites/g/files/rcwecm3506/files/2020-05/wilderandfelger2010_martirflora.pdf",
        "Flora of Isla San Pedro Martir", "Wilder & Felger 2010 regional flora.",
        "Bahiopsis chenopodina. Shrubs with slender, brittle stems. Flower heads with bright yellow ray and disk florets.",
        "flora:wilder-felger-2010:bahiopsis-chenopodina", "A", "primary_regional_flora", "en",
        '"Bahiopsis chenopodina" yellow flowers', 6,
    ),
    (
        "Sapium daphnoides", "floral_structural_complexity",
        "inflorescence_display", "simple axillary spikes", "raceme_spike_panicle", "high",
        "Biodiversity Heritage Library",
        "https://upload.wikimedia.org/wikipedia/commons/a/af/Nachrichten_von_der_K._Gesellschaft_der_Wissenschaften_und_der_Georg-Augusts-Universit%C3%A4t_%28IA_nachrichtengott1865%29.pdf",
        "Nachrichten von der K. Gesellschaft der Wissenschaften",
        "Grisebach 1865 protologue of Sapium daphnoides.",
        "Sapium daphnoides Gr. foliis lanceolato-oblongis; spicis axillaribus folio multo brevioribus simplicibus.",
        "protologue:grisebach-1865:sapium-daphnoides", "A", "primary_taxonomic_protologue", "la",
        '"Sapium daphnoides" spicis axillaribus', 13,
    ),
    (
        "Jurinea consanguinea", "floral_structural_complexity",
        "inflorescence_display", "single capitulum", "composite_display", "high",
        "Istanbul University Flora", "https://ibuflora.ibu.edu.tr/tur/jurinea-consanguinea",
        "Jurinea consanguinea", "Istanbul University exact-species flora treatment.",
        "Uzun çiçek saplarının üstünde kapitulum tek.",
        "university-flora:ibu:jurinea-consanguinea", "A", "university_flora_species_treatment", "tr",
        '"Jurinea consanguinea" kapitulum tek', 7,
    ),
    (
        "Reichardia intermedia", "floral_structural_complexity",
        "floral_form", "capitulum with 10-20 flowers", "composite_head", "high",
        "Macaronesian Botany", "https://www.macaronesian.org/assets/files/file-3c45597466126b.pdf",
        "Revision of Reichardia", "Peer-reviewed taxonomic revision of Reichardia.",
        "El numero menor de flores por capitulo lo presentan los individuos de pequeno tamano de R. intermedia (10-20 flores).",
        "journal-revision:reichardia:intermedia", "A", "primary_taxonomic_revision", "es",
        '"Reichardia intermedia" flores capitulo', 6,
    ),
    (
        "Epithema benthamii", "floral_structural_complexity",
        "tube_depth_class", "corolla tube 5-6(-8) mm", "intermediate", "high",
        "Gardens' Bulletin Singapore",
        "https://www.nparks.gov.sg/sbg/research/publications/gardens-bulletin-singapore/-/media/sbg/gardens-bulletin/gbs_67_01_y2015_v67_01/4-4-67-1-159-y2015-v67p1-gbs-pg159.pdf",
        "A revision of Epithema", "Peer-reviewed taxonomic revision in Gardens' Bulletin Singapore 67.",
        "Epithema benthamii: corolla 5-10 mm long; tube (4-)5-6(-8) by 1.3-2.5 mm; lobes 1-2 mm long.",
        "journal-revision:gbs:epithema-benthamii", "A", "primary_taxonomic_revision", "en",
        '"Epithema benthamii" corolla tube', 7,
    ),
    (
        "Cosmianthemum subglabrum", "floral_structural_complexity",
        "tube_depth_class", "corolla tube 6-8 mm", "intermediate", "high",
        "Royal Botanic Garden Edinburgh",
        "https://journals.rbge.org.uk/notes/article/download/2686/2506",
        "A revision of Cosmianthemum", "Peer-reviewed taxonomic revision in Notes RBGE.",
        "Cosmianthemum subglabrum: corolla cream to greenish-white; tube lacking a dorsal pouch 6-8 mm long.",
        "journal-revision:rbge:cosmianthemum-subglabrum", "A", "primary_taxonomic_revision", "en",
        '"Cosmianthemum subglabrum" corolla tube', 7,
    ),
]


REJECTED_COLUMNS = [
    "accepted_species", "queried_trait", "queued_value", "observed_statement",
    "decision", "decision_reason", "source_url", "source_lineage", "review_status",
]

REJECTED: list[tuple[str, ...]] = [
    ("Flindersia australis", "flower_size_class", "large", "petal length only", "reject", "organ measurement is not whole-flower size", "", "", "reviewed_fail_closed"),
    ("Satureja hortensis", "tube_depth_class", "intermediate", "whole corolla length only", "reject", "tube depth was not separately stated", "", "", "reviewed_fail_closed"),
    ("Grammatophyllum speciosum", "flower_primary_color", "yellow_orange", "yellow-green and patterned flowers", "reject", "multistate description does not support queued single state", "", "", "reviewed_fail_closed"),
    ("Citronella moorei", "flower_primary_color", "white", "white and red states described", "reject", "multistate description does not support queued single state", "", "", "reviewed_fail_closed"),
    ("Pachira aquatica", "inflorescence_display", "solitary", "solitary or small groups", "reject", "variable display does not support queued single state", "", "", "reviewed_fail_closed"),
    ("Spathelia terminalioides", "flower_primary_color", "red", "source described white flowers", "reject", "direct contradiction of queued value", "", "", "reviewed_fail_closed"),
    ("Sabal palmetto", "flower_size_class", "very_small", "flowers 5-6 mm", "reject", "measurement contradicts queued size class", "https://www.botanybrisbane.com/plants/arecaceae/sabal/sabal-palmetto/", "provider-treatment:botanybrisbane:sabal-palmetto", "reviewed_fail_closed"),
    ("Stenogyne rugosa", "flower_primary_color", "yellow_orange", "multiple flower colours", "reject", "multistate description does not support queued single state", "", "", "reviewed_fail_closed"),
    ("Leontopodium alpinum", "flower_primary_color", "white", "white refers to bracts; corollas are yellow", "reject", "wrong floral organ", "", "", "reviewed_fail_closed"),
    ("Calytrix tetragona", "floral_form", "funnel", "star-shaped flowers", "reject", "direct contradiction of queued form", "", "", "reviewed_fail_closed"),
    ("Chamaecytisus prolifer", "flower_primary_color", "white", "Australian usage may be a misapplied name", "reject", "strict synonym identity unresolved", "", "", "reviewed_fail_closed"),
    ("Orthoceras annulatum", "floral_form", "tubular", "taxon name has homonym ambiguity", "reject", "strict species identity unresolved", "", "", "reviewed_fail_closed"),
    ("Conyza incana", "inflorescence_display", "composite_display", "name usage taxonomically ambiguous", "reject", "strict species identity unresolved", "", "", "reviewed_fail_closed"),
]


def primary_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for record in RECORDS:
        *fields, _potential_cells = record
        rows.append(_row(*fields))
    return rows


def build_audit(evidence: pd.DataFrame) -> pd.DataFrame:
    audit = _wave15_build_audit(evidence)
    audit["reviewer"] = "Codex Wave24 support-two exact-species source audit"
    audit["reviewed_at_utc"] = RETRIEVED_AT
    audit["decision_reason"] = (
        "Accepted after exact fixed-master species and family match, original-page "
        "retrieval, exact quote and source-lineage review, cultivar screening, and "
        "trait-specific fail-closed ontology mapping."
    )
    return audit.loc[:, AUDIT_COLUMNS]


def build_checkpoint(output_dir: Path = CHECKPOINT) -> dict[str, Any]:
    evidence = pd.DataFrame(primary_rows(), columns=EVIDENCE_COLUMNS).sort_values(
        ["accepted_species", "trait_name", "source_lineage"], kind="stable"
    ).reset_index(drop=True)
    audit = build_audit(evidence)
    rejected = pd.DataFrame(REJECTED, columns=REJECTED_COLUMNS)
    if len(evidence) != 22 or len(audit) != 22:
        raise ValueError("Wave24 must contain exactly 22 individually reviewed rows")
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("Wave24 candidate IDs must be unique")
    if evidence.duplicated(["accepted_species", "trait_name"]).any():
        raise ValueError("Wave24 species x trait pairs must be unique")

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
        "evidence": output_dir / "targeted_support2_wave24_evidence_20260824.csv",
        "audit": output_dir / "targeted_support2_wave24_manual_audit_20260824.csv",
        "rejected": output_dir / "targeted_support2_wave24_rejected_candidates_20260824.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260824.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260824.csv",
        "manifest": output_dir / "source_acquisition_manifest_wave24.json",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    rejected.to_csv(paths["rejected"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    targeted_rules = sorted(
        {
            f"{str(record[0]).split()[0]}|{record[2]}"
            for record in RECORDS
        }
    )
    manifest = {
        "checkpoint": SOURCE_GROUP,
        "built_at_utc": RETRIEVED_AT,
        "fixed_integrated_baseline_run_id": FIXED_INTEGRATED_BASELINE_RUN_ID,
        "baseline_formal_run_id": PRIOR_FORMAL_RUN_ID,
        "accepted_evidence_rows": len(evidence),
        "accepted_species_trait": len(evidence),
        "accepted_species": int(evidence["accepted_species"].nunique()),
        "accepted_source_lineages": int(evidence["source_lineage"].nunique()),
        "quality_counts": evidence["evidence_quality"].value_counts().to_dict(),
        "rejected_reviewed_candidates": len(rejected),
        "recorded_queries": int(evidence["query"].nunique()),
        "formal_search_api_queries": 0,
        "search_cost_usd": 0.0,
        "targeted_support2_rules": targeted_rules,
        "theoretical_rule_cells_touched": sum(int(row[-1]) for row in RECORDS),
        "guardrails": {
            "search_snippet_as_evidence": False,
            "family_inference": False,
            "global_fallback": False,
            "min_species_two_production": False,
            "cross_trait_substitution": False,
            "genus_axis_only_join": False,
            "cultivar_or_hybrid_transferred_to_wild_species": False,
        },
        "output_sha256": {
            key: _sha(path.read_text(encoding="utf-8"))
            for key, path in paths.items()
            if key != "manifest"
        },
        "notes": (
            "Potential cells are a queue ceiling, not observed coverage gain. "
            "The shared rebuild must pass dominance, masked-species validation, "
            "and independent source-lineage leave-one-out validation. Rejected "
            "ambiguous, conflicting, multistate, wrong-organ, or incomplete-trait "
            "candidates remain visible and are never promoted."
        ),
    }
    paths["manifest"].write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=CHECKPOINT)
    args = parser.parse_args()
    print(json.dumps(build_checkpoint(args.output_dir), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
