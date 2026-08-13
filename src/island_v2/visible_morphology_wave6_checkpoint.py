"""Freeze reviewed visible-morphology evidence for support=2 genus rules.

The checkpoint stores species/synonym-direct statements only.  It preserves
explicit multistate descriptions and emits no genus inference; rule potential
is a queue diagnostic until the formal all-evidence rebuild succeeds.
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

CREATED_AT = "2026-08-14T08:00:00Z"
SOURCE_GROUP = "visible_morphology_wave6_checkpoint_20260814"


SOURCES: tuple[dict[str, object], ...] = (
    {
        "species": "Kayea racemosa",
        "family": "Calophyllaceae",
        "trait": "inflorescence_display",
        "value": "raceme_spike_panicle",
        "quality": "medium",
        "provider": "PROSEA species treatment via Pl@ntUse",
        "url": "https://plantuse.plantnet.org/en/index.php?title=Mesua_racemosa_(PROSEA)&oldid=208222",
        "title": "Mesua racemosa (PROSEA)",
        "citation": "PROSEA timber species treatment, fixed revision 208222",
        "excerpt": (
            "Mesua racemosa (Planchon ex Triana & Planchon) Kosterm. Synonyms: "
            "Kayea racemosa Planchon ex Triana & Planchon (1861). A medium-sized "
            "tree up to 24 m tall; flowers 3 together in a short compact raceme."
        ),
        "raw_value": "flowers 3 together in a short compact raceme",
        "record_id": "prosea:mesua-racemosa:oldid-208222:compact-raceme",
        "lineage": "species_treatment:prosea:mesua_racemosa",
        "lineage_method": "canonical_species_treatment_revision",
        "source_tier": "A",
        "source_type": "curated_botanical_monograph_species_treatment",
        "domain": "plantuse.plantnet.org",
        "content_sha256": "66f11674d3d07dad0bbaf27517160cc77b03d80ccdcee5a077d3f08f7f98c246",
        "content_sha256_basis": "downloaded_fixed_revision_html_bytes",
        "potential": 22,
        "evidence_scope": "synonym_direct",
        "name_match_method": "exact_synonym",
        "matched_page_name": "Mesua racemosa",
        "name_resolution_lineage": "prosea_explicit_synonym_kayea_racemosa",
    },
    {
        "species": "Blepharis ciliaris",
        "family": "Acanthaceae",
        "trait": "flower_size_class",
        "value": "medium",
        "quality": "high",
        "provider": "Flora of Pakistan via World Flora Online",
        "url": "https://www.worldfloraonline.org/taxon/wfo-0000566595",
        "title": "Blepharis ciliaris (L.) B.L.Burtt",
        "citation": "Blepharis ciliaris, Flora of Pakistan treatment",
        "excerpt": (
            "Inflorescence a dense, strobilate spike. Flowers purplish-blue, "
            "fading to white, 1.5-2 cm long."
        ),
        "raw_value": "flowers 1.5-2 cm long",
        "record_id": "wfo-0000566595:flora-pakistan:flower-size",
        "lineage": "provider_treatment:flora_of_pakistan:wfo-0000566595",
        "lineage_method": "original_flora_species_treatment",
        "source_tier": "A",
        "source_type": "official_flora_species_treatment",
        "domain": "worldfloraonline.org",
        "content_sha256": "1214c2eccbb88bc165f552edc1aee7798555cf8c8d24c0fefa95dc7faa490c39",
        "content_sha256_basis": "downloaded_wfo_html_bytes",
        "potential": 13,
    },
    {
        "species": "Blepharis ciliaris",
        "family": "Acanthaceae",
        "trait": "inflorescence_display",
        "value": "raceme_spike_panicle",
        "quality": "high",
        "provider": "Flora of Pakistan via World Flora Online",
        "url": "https://www.worldfloraonline.org/taxon/wfo-0000566595",
        "title": "Blepharis ciliaris (L.) B.L.Burtt",
        "citation": "Blepharis ciliaris, Flora of Pakistan treatment",
        "excerpt": (
            "Inflorescence a dense, strobilate spike. Flowers purplish-blue, "
            "fading to white, 1.5-2 cm long."
        ),
        "raw_value": "inflorescence a dense, strobilate spike",
        "record_id": "wfo-0000566595:flora-pakistan:inflorescence",
        "lineage": "provider_treatment:flora_of_pakistan:wfo-0000566595",
        "lineage_method": "original_flora_species_treatment",
        "source_tier": "A",
        "source_type": "official_flora_species_treatment",
        "domain": "worldfloraonline.org",
        "content_sha256": "1214c2eccbb88bc165f552edc1aee7798555cf8c8d24c0fefa95dc7faa490c39",
        "content_sha256_basis": "downloaded_wfo_html_bytes",
        "potential": 13,
    },
    {
        "species": "Hedycarya arborea",
        "family": "Monimiaceae",
        "trait": "flower_size_class",
        "value": "small",
        "quality": "medium",
        "provider": "University of Auckland NZ Plants",
        "url": "https://www.nzplants.auckland.ac.nz/en/about/seed-plants-flowering/monimiaceae/hedycarya-arborea.html",
        "title": "Hedycarya arborea - porokaiwhiri, pigeonwood",
        "citation": "University of Auckland NZ Plants species page",
        "excerpt": (
            "Hedycarya arborea - porokaiwhiri, pigeonwood. Family: Monimiaceae. "
            "Reproductive characteristics: Flower size: 5-12 mm diam."
        ),
        "raw_value": "flower size 5-12 mm diameter",
        "record_id": "auckland-nzplants:hedycarya-arborea:flower-size",
        "lineage": "provider_treatment:university_of_auckland_nzplants:hedycarya_arborea",
        "lineage_method": "university_species_treatment",
        "source_tier": "A",
        "source_type": "university_botanical_species_page",
        "domain": "nzplants.auckland.ac.nz",
        "content_sha256": "7172d92155689464933b43c77b32ac869d586b15957cc22c054c536c7ba5f28d",
        "content_sha256_basis": "downloaded_university_html_bytes",
        "potential": 13,
    },
    {
        "species": "Cassipourea elliptica",
        "family": "Rhizophoraceae",
        "trait": "inflorescence_display",
        "value": "solitary|umbel_corymb",
        "quality": "high",
        "provider": "Flora de Nicaragua via World Flora Online",
        "url": "https://www.worldfloraonline.org/taxon/wfo-4000006916",
        "title": "Cassipourea Aubl. - Flora de Nicaragua treatment",
        "citation": "Barrie, Flora de Nicaragua, Cassipourea elliptica treatment",
        "excerpt": (
            "Cassipourea elliptica (Sw.) Poir. Inflorescencia usualmente de "
            "numerosas flores en fascículos en las axilas de las hojas o "
            "raramente flores solitarias."
        ),
        "raw_value": "usually many flowers in axillary fascicles; rarely solitary flowers",
        "record_id": "wfo-4000006916:flora-nicaragua:cassipourea-elliptica:inflorescence",
        "lineage": "provider_treatment:flora_de_nicaragua:cassipourea_elliptica",
        "lineage_method": "original_flora_species_treatment",
        "source_tier": "A",
        "source_type": "official_flora_species_treatment",
        "domain": "worldfloraonline.org",
        "language": "es",
        "content_sha256": "0f39610b15e505aa8a09b59cfc372f0ae15652866095a4a5c2ce31eee9fed892",
        "content_sha256_basis": "downloaded_wfo_html_bytes",
        "potential": 0,
        "counterexample": True,
    },
    {
        "species": "Scolopia rhamniphylla",
        "family": "Salicaceae",
        "trait": "inflorescence_display",
        "value": "raceme_spike_panicle",
        "quality": "high",
        "provider": "Flora of Tropical East Africa via Kew POWO",
        "url": "https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:365942-1/general-information",
        "title": "Scolopia rhamniphylla Gilg",
        "citation": "Sleumer (1975), Flora of Tropical East Africa",
        "excerpt": (
            "Racemes from foliate and defoliate axils, solitary or rarely in "
            "pairs, rather few- and lax-flowered, sometimes reduced to fascicles, "
            "1-2(-4) cm long."
        ),
        "raw_value": "racemes; sometimes reduced to fascicles",
        "record_id": "powo:365942-1:ftea:inflorescences",
        "lineage": "provider_treatment:flora_tropical_east_africa:scolopia_rhamniphylla",
        "lineage_method": "original_flora_species_treatment",
        "source_tier": "A",
        "source_type": "official_flora_species_treatment",
        "domain": "powo.science.kew.org",
        "content_sha256": "ae2d6fe9da370b6b0c23b1a7ef3839eb11906b0e44bb446cd357dd732822476f",
        "content_sha256_basis": "exact_rendered_kew_treatment_excerpt_utf8_bytes",
        "potential": 0,
    },
    {
        "species": "Hypobathrum frutescens",
        "family": "Rubiaceae",
        "trait": "flower_primary_color",
        "value": "green_brown_inconspicuous|white",
        "quality": "medium",
        "provider": "Cibodas Botanical Garden primary morphology study",
        "url": "https://pasca.unhas.ac.id/ojs/index.php/ijas/article/view/2718/749",
        "title": "Ornamental Plant's Potentials of Indonesian Native Rubiaceae Collected in Cibodas Botanical Garden",
        "citation": "Putri, Junaedi & Hendrian (2021), DOI 10.20956/ijas.v9i1.2718",
        "excerpt": (
            "Hypobathrum frutescens has symmetrical branch architecture, "
            "white inflorescence flowers with pale green corolla."
        ),
        "raw_value": "white flowers with pale green corolla",
        "record_id": "doi:10.20956/ijas.v9i1.2718:hypobathrum-frutescens:flower-colour",
        "lineage": "doi:10.20956/ijas.v9i1.2718",
        "lineage_method": "peer_reviewed_primary_morphology_observation",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_botanical_garden_study",
        "domain": "pasca.unhas.ac.id",
        "content_sha256": "1c5bcce9d83e5aaac9cc717725c7961519edaf3843cbf9b75293fdbbfab8e036",
        "content_sha256_basis": "downloaded_primary_article_pdf_bytes",
        "potential": 21,
    },
    {
        "species": "Pseudognaphalium luteoalbum",
        "family": "Asteraceae",
        "trait": "tube_depth_class",
        "value": "shallow",
        "quality": "high",
        "provider": "Flora of Tropical East Africa via Kew POWO",
        "url": "https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:909623-1/general-information",
        "title": "Pseudognaphalium luteoalbum (L.) Hilliard & B.L.Burtt",
        "citation": "Beentje, Jeffrey & Hind (2005), Flora of Tropical East Africa",
        "excerpt": (
            "Outer florets white or yellow, about 100, tube filiform, 1.6-2.2 mm "
            "long; inner florets 5-15(-30), tube cylindric, 1.3-2 mm long."
        ),
        "raw_value": "outer floret tubes 1.6-2.2 mm; inner floret tubes 1.3-2 mm",
        "record_id": "powo:909623-1:ftea:floret-tubes",
        "lineage": "provider_treatment:flora_tropical_east_africa:pseudognaphalium_luteoalbum",
        "lineage_method": "original_flora_species_treatment",
        "source_tier": "A",
        "source_type": "official_flora_species_treatment",
        "domain": "powo.science.kew.org",
        "content_sha256": "2efa4fa3e25b7ce59910a44115b987b3aa96ae905c2ad67c15744da3a19610ef",
        "content_sha256_basis": "exact_rendered_kew_treatment_excerpt_utf8_bytes",
        "potential": 15,
    },
    {
        "species": "Guettardella caudata",
        "family": "Rubiaceae",
        "trait": "inflorescence_display",
        "value": "umbel_corymb",
        "quality": "high",
        "provider": "Naturalis primary taxonomic revision",
        "url": "https://repository.naturalis.nl/pub/524956/BLUM1984029002017.pdf",
        "title": "A synopsis of Guettardella and the Old World species of Antirhea",
        "citation": "Jansen (1984), Blumea 29:565-588, Guettardella caudata",
        "excerpt": (
            "Guettardella caudata M.E. Jansen, spec. nov. Inflorescence axillary, "
            "solitary, dichotomous, with either 3-8(-11) male flowers or 1 or 2 "
            "female flowers; peduncle 18-42 mm long."
        ),
        "raw_value": "solitary axillary dichotomous inflorescence with 1-11 flowers",
        "record_id": "jansen-1984:guettardella-caudata:inflorescence",
        "lineage": "primary_revision:jansen_1984_guettardella:guettardella_caudata",
        "lineage_method": "primary_taxonomic_revision_species_treatment",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_taxonomic_revision",
        "domain": "repository.naturalis.nl",
        "content_sha256": "f87217c1e5bde5e320e209a75fbb82bccbaf54f50a75340294ad47d80ba094ff",
        "content_sha256_basis": "downloaded_primary_revision_pdf_bytes",
        "potential": 21,
    },
    {
        "species": "Leucosyke capitellata",
        "family": "Urticaceae",
        "trait": "inflorescence_display",
        "value": "umbel_corymb",
        "quality": "medium",
        "provider": "PROSEA species treatment via Pl@ntUse",
        "url": "https://plantuse.plantnet.org/en/index.php?title=Leucosyke_capitellata_(PROSEA)&oldid=333186",
        "title": "Leucosyke capitellata (PROSEA)",
        "citation": "Lemmens, PROSEA species treatment, fixed revision 333186",
        "excerpt": (
            "Leucosyke capitellata (Poir.) Wedd. Inflorescence a pseudo-axillary, "
            "peduncled, globose head 0.5-1 cm in diameter, often 2 heads close "
            "together."
        ),
        "raw_value": "peduncled globose head, often two heads close together",
        "record_id": "prosea:leucosyke-capitellata:oldid-333186:inflorescence",
        "lineage": "species_treatment:prosea:leucosyke_capitellata",
        "lineage_method": "canonical_species_treatment_revision",
        "source_tier": "A",
        "source_type": "curated_botanical_monograph_species_treatment",
        "domain": "plantuse.plantnet.org",
        "content_sha256": "89aa5dee4bc94d1c9fcf382677d2d97e015075b0a99eba7d9f02cd201c858ef7",
        "content_sha256_basis": "downloaded_fixed_revision_html_bytes",
        "potential": 15,
    },
    {
        "species": "Rungia pectinata",
        "family": "Acanthaceae",
        "trait": "inflorescence_display",
        "value": "raceme_spike_panicle",
        "quality": "high",
        "provider": "Flora of China via World Flora Online",
        "url": "https://www.worldfloraonline.org/taxon/wfo-0000401426",
        "title": "Rungia pectinata (L.) Nees",
        "citation": "Rungia pectinata, Flora of China treatment",
        "excerpt": (
            "Spikes axillary or terminal, 0.5-2 cm, 1-sided, solitary or "
            "sometimes 2 or 3 compound. Corolla blue or white, ca. 5 mm."
        ),
        "raw_value": "axillary or terminal spikes, solitary or 2-3 compound",
        "record_id": "wfo-0000401426:flora-china:inflorescence",
        "lineage": "provider_treatment:flora_of_china:rungia_pectinata",
        "lineage_method": "original_flora_species_treatment",
        "source_tier": "A",
        "source_type": "official_flora_species_treatment",
        "domain": "worldfloraonline.org",
        "content_sha256": "a6b6db43ded530dd88c0262219931453cd9720f9326f7f9065192b37e51c3f5e",
        "content_sha256_basis": "downloaded_wfo_html_bytes",
        "potential": 13,
    },
)


def reviewed_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for source in SOURCES:
        evidence = _evidence_row(
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
        evidence["source_group"] = SOURCE_GROUP
        evidence["query"] = "support2_visible_morphology_third_species"
        evidence["language"] = str(source.get("language", "en"))
        for column in (
            "evidence_scope",
            "name_match_method",
            "matched_page_name",
            "name_resolution_lineage",
        ):
            if column in source:
                evidence[column] = str(source[column])
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
    master_identity = master.set_index("accepted_species")["family"].to_dict()
    expected_identity = {str(source["species"]): str(source["family"]) for source in SOURCES}
    missing = sorted(set(expected_identity) - set(master_identity))
    if missing:
        raise ValueError(f"checkpoint species missing from master: {missing}")
    conflicts = {
        species: (master_identity[species], family)
        for species, family in expected_identity.items()
        if master_identity[species] != family
    }
    if conflicts:
        raise ValueError(f"checkpoint family conflicts: {conflicts}")

    evidence = pd.DataFrame(reviewed_rows(), columns=EVIDENCE_COLUMNS).fillna("")
    if evidence[["accepted_species", "trait_name"]].duplicated().any():
        raise ValueError("checkpoint species-trait pairs must be unique")
    if not evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all():
        raise ValueError("every checkpoint row requires a 64-character content hash")

    audit = _audit(evidence)
    audit["reviewer"] = "Codex visible morphology source audit"
    audit["reviewed_at_utc"] = CREATED_AT
    audit["decision_reason"] = (
        "Accepted after exact target-master identity, whole-trait statement, "
        "value state, full source, lineage, provenance and cultivar scope were "
        "checked. Explicit rare or multiple states were retained."
    )

    prior_evidence = pd.read_csv(prior_curated_evidence_csv, dtype=str).fillna("")
    prior_audit = pd.read_csv(prior_curated_audit_csv, dtype=str).fillna("")
    owned = prior_evidence["source_group"].eq(SOURCE_GROUP)
    prior_owned_ids = set(prior_evidence.loc[owned, "candidate_id"])
    current_ids = set(evidence["candidate_id"])
    combined_evidence = pd.concat(
        [prior_evidence.loc[~owned], evidence], ignore_index=True
    )
    combined_audit = pd.concat(
        [
            prior_audit.loc[
                ~prior_audit["candidate_id"].isin(prior_owned_ids | current_ids)
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
        "evidence": output_dir / "visible_morphology_wave6_evidence_20260814.csv",
        "audit": output_dir / "visible_morphology_wave6_manual_audit_20260814.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260814.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260814.csv",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    summary: dict[str, object] = {
        "contract": "visible_morphology_support2_checkpoint_v1",
        "created_at_utc": CREATED_AT,
        "accepted_rows": len(evidence),
        "species": int(evidence["accepted_species"].nunique()),
        "species_trait": len(evidence),
        "quality_counts": evidence["evidence_quality"].value_counts().to_dict(),
        "trait_counts": evidence["trait_name"].value_counts().to_dict(),
        "source_lineages": int(evidence["source_lineage"].nunique()),
        "rule_candidates": int(sum(int(source["potential"]) > 0 for source in SOURCES)),
        "explicit_counterexample_rows": int(sum(bool(source.get("counterexample")) for source in SOURCES)),
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
            "explicit_multistate_retained": True,
            "genus_inference_emitted_here": False,
            "family_inference": False,
            "global_fallback": False,
            "n2_formal_inference": False,
            "cross_trait_substitution": False,
            "search_snippet_evidence": False,
            "partial_organ_size_as_whole_flower_size": False,
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
    readme_path = output_dir / "README.md"
    hash_paths = list(paths.values())
    if readme_path.exists():
        hash_paths.append(readme_path)
    for path in hash_paths:
        summary["files"][path.name] = {
            "sha256": _sha256(path),
            "size_bytes": len(_canonical_file_bytes(path)),
        }
    manifest_path = output_dir / "visible_morphology_wave6_manifest_20260814.json"
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
