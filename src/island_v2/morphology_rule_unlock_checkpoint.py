"""Freeze independent morphology evidence targeted at support-2 genus rules.

The rows in this checkpoint were re-read on their original species pages.  It
does not emit genus inference: the common all-evidence implementation decides
whether each genus x trait rule remains eligible after adding the direct rows.
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

CREATED_AT = "2026-08-12T19:30:00Z"
SOURCE_GROUP = "morphology_rule_unlock_checkpoint_20260813"


ROWS = [
    {
        "species": "Gmelina arborea",
        "trait": "floral_symmetry",
        "value": "zygomorphic",
        "quality": "high",
        "provider": "PROTA / Pl@ntUse",
        "url": "https://plantuse.plantnet.org/en/Gmelina_arborea_(PROTA)",
        "title": "Gmelina arborea (PROTA)",
        "citation": "Adam & Krampah (2005), Gmelina arborea, PROTA treatment",
        "excerpt": "Flowers bisexual, zygomorphic, 5-merous;",
        "record_id": "prota:110130:revision:325985:flowers",
        "lineage": "provider_treatment:prota:110130",
        "source_tier": "A",
        "source_type": "authored_species_monograph_database",
        "domain": "plantuse.plantnet.org",
        "content_sha256": "49db5019226458e88ec4ac4f1684a702c4c2f856dd324b0421558bdd188e4353",
        "raw_value": "zygomorphic",
    },
    {
        "species": "Corylopsis spicata",
        "trait": "floral_form",
        "value": "bell_campanulate",
        "quality": "medium",
        "provider": "NC State Extension Gardener Plant Toolbox",
        "url": "https://plants.ces.ncsu.edu/plants/corylopsis-spicata/",
        "title": "Corylopsis spicata",
        "citation": "NC State Extension Gardener Plant Toolbox species profile",
        "excerpt": "Flower Shape: Bell",
        "record_id": "ncsu:corylopsis-spicata:flower-shape",
        "lineage": "url:https://plants.ces.ncsu.edu/plants/corylopsis-spicata",
        "source_tier": "B",
        "source_type": "university_extension_species_profile",
        "domain": "plants.ces.ncsu.edu",
        "content_sha256": "925ca589ff1488cd479c3af81fc856146d0780ba11f740c460c4236fb0c41d03",
        "raw_value": "Bell",
    },
    {
        "species": "Daphniphyllum macropodum",
        "trait": "inflorescence_display",
        "value": "raceme_spike_panicle",
        "quality": "medium",
        "provider": "NC State Extension Gardener Plant Toolbox",
        "url": "https://plants.ces.ncsu.edu/plants/daphniphyllum-macropodum/",
        "title": "Daphniphyllum macropodum",
        "citation": "NC State Extension Gardener Plant Toolbox species profile",
        "excerpt": "The flowers are apetalous blooms, arranged in racemes, and emerge from the leaf axils.",
        "record_id": "ncsu:daphniphyllum-macropodum:racemes",
        "lineage": "url:https://plants.ces.ncsu.edu/plants/daphniphyllum-macropodum",
        "source_tier": "B",
        "source_type": "university_extension_species_profile",
        "domain": "plants.ces.ncsu.edu",
        "content_sha256": "d5b68e6514bd7d4da377fbf3bf5f00c0aa82df5e878a7b55c9ac3b40483b992c",
        "raw_value": "racemes",
    },
    {
        "species": "Salacca affinis",
        "trait": "inflorescence_display",
        "value": "raceme_spike_panicle",
        "quality": "high",
        "provider": "NParks Flora & Fauna Web",
        "url": "https://www.nparks.gov.sg/florafaunaweb/flora/7/2/7297",
        "title": "Salacca affinis Griff.",
        "citation": "NParks Flora & Fauna Web species treatment",
        "excerpt": "Its inflorescences are erect, with male flower spikes about 2.5-6.4 cm",
        "record_id": "nparks:flora:7297:inflorescence",
        "lineage": "url:https://www.nparks.gov.sg/florafaunaweb/flora/7/2/7297",
        "source_tier": "A",
        "source_type": "government_species_database",
        "domain": "nparks.gov.sg",
        "content_sha256": "000dccf14e3185b8c689536c1c0b473847a28b0d19b956a9fed6a75e180d6210",
        "raw_value": "flower spikes",
    },
    {
        "species": "Paris japonica",
        "trait": "inflorescence_display",
        "value": "solitary",
        "quality": "high",
        "provider": "Royal Botanic Gardens, Kew Species Profiles",
        "url": "https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:539710-1/general-information",
        "title": "Paris japonica (Franch. & Sav.) Franch.",
        "citation": "Royal Botanic Gardens, Kew Species Profile",
        "excerpt": "In Japan P. japonica blooms between May and August, producing a single, stalked, star-like flower composed of up to 10 whitish tepals (2.5-3.5 cm).",
        "record_id": "kew-species-profile:539710-1:single-flower",
        "lineage": "provider_treatment:kew_species_profile:urn:lsid:ipni.org:names:539710-1",
        "source_tier": "A",
        "source_type": "botanic_garden_species_profile",
        "domain": "powo.science.kew.org",
        "content_sha256": "f9f35698cebfd92c36e5f16257db06141f199854dd42ea92975206363b65642b",
        "raw_value": "single flower",
    },
    {
        "species": "Polyspora singaporeana",
        "trait": "flower_size_class",
        "value": "large",
        "quality": "high",
        "provider": "NParks Flora & Fauna Web",
        "url": "https://www.nparks.gov.sg/florafaunaweb/flora/6/7/6729",
        "title": "Polyspora singaporeana",
        "citation": "NParks Flora & Fauna Web species treatment",
        "excerpt": "Its stalked flowers are fleshy, scentless, about 3.8 cm wide",
        "record_id": "nparks:flora:6729:flower-size",
        "lineage": "url:https://www.nparks.gov.sg/florafaunaweb/flora/6/7/6729",
        "source_tier": "A",
        "source_type": "government_species_database",
        "domain": "nparks.gov.sg",
        "content_sha256": "3303e96044dc777b841499f6421ed5cbc10629db3ae943b9cda7ab3ffecaf687",
        "raw_value": "3.8 cm wide",
    },
    {
        "species": "Salacca magnifica",
        "trait": "inflorescence_display",
        "value": "raceme_spike_panicle",
        "quality": "high",
        "provider": "NParks Flora & Fauna Web",
        "url": "https://www.nparks.gov.sg/florafaunaweb/flora/5/1/5126",
        "title": "Salacca magnifica",
        "citation": "NParks Flora & Fauna Web species treatment",
        "excerpt": "The female inflorescences are erect and spike-like",
        "record_id": "nparks:flora:5126:inflorescence",
        "lineage": "url:https://www.nparks.gov.sg/florafaunaweb/flora/5/1/5126",
        "source_tier": "A",
        "source_type": "government_species_database",
        "domain": "nparks.gov.sg",
        "content_sha256": "e9e9e2c1c54b53a38e928fe28f0722ab97ae1d00676c3e23ac0d072350f4184f",
        "raw_value": "spike-like",
    },
    {
        "species": "Paronychia echinulata",
        "trait": "flower_size_class",
        "value": "very_small",
        "quality": "high",
        "provider": "Flora Iberica, Real Jardin Botanico-CSIC",
        "url": "https://www.floraiberica.es/floraiberica/texto/pdfs/02_049_03_Paronychia.pdf",
        "title": "Flora Iberica: Paronychia",
        "citation": "Chaudhri, Paronychia treatment, Flora Iberica volume II",
        "excerpt": "Flores (1,4)1,5-2 mm -excluidas las aristas-",
        "record_id": "flora-iberica:Paronychia_echinulata:flower-size",
        "lineage": "flora_treatment:flora_iberica:Paronychia_echinulata",
        "source_tier": "A",
        "source_type": "official_taxonomic_flora",
        "domain": "floraiberica.es",
        "content_sha256": "fd466aec191b9aed5e7f01b874e539657673e25031d327c56a690010d32f8fab",
        "raw_value": "(1.4)1.5-2 mm",
    },
    {
        "species": "Corylopsis pauciflora",
        "trait": "floral_form",
        "value": "bell_campanulate",
        "quality": "medium",
        "provider": "Royal Horticultural Society",
        "url": "https://www.rhs.org.uk/plants/4501/i-corylopsis-pauciflora-i/details",
        "title": "Corylopsis pauciflora",
        "citation": "Royal Horticultural Society species page",
        "excerpt": "Fragrant, pale yellow, bell-shaped flowers",
        "record_id": "rhs:4501:bell-shaped-flowers",
        "lineage": "url:https://www.rhs.org.uk/plants/4501/i-corylopsis-pauciflora-i/details",
        "source_tier": "B",
        "source_type": "specialist_horticultural_species_profile",
        "domain": "rhs.org.uk",
        "content_sha256": "89066e0f272abd2a3b10084c901492d9930e36f954b327c576314f72530931cd",
        "raw_value": "bell-shaped",
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
            lineage_method="original_species_treatment_or_page_lineage",
            source_tier=item["source_tier"],
            source_type=item["source_type"],
            domain=item["domain"],
            content_sha256=item["content_sha256"],
            content_sha256_basis="retrieved_original_page_or_pdf_bytes",
            retrieved_at_utc="2026-08-12T19:25:00Z",
            raw_value=item["raw_value"],
        )
        evidence["source_group"] = SOURCE_GROUP
        evidence["query"] = "support2_morphology_rule_unlock_original_source"
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
        raise ValueError(f"morphology checkpoint species missing from master: {missing}")
    if evidence[["accepted_species", "trait_name"]].duplicated().any():
        raise ValueError("morphology checkpoint species-trait pairs are not unique")
    audit = _audit(evidence)

    prior_evidence = pd.read_csv(prior_curated_evidence_csv, dtype=str).fillna("")
    prior_audit = pd.read_csv(prior_curated_audit_csv, dtype=str).fillna("")
    prior_owned = prior_evidence["source_group"].eq(SOURCE_GROUP)
    prior_owned_ids = set(prior_evidence.loc[prior_owned, "candidate_id"])
    owned_ids = set(evidence["candidate_id"])
    combined_evidence = pd.concat(
        [prior_evidence.loc[~prior_owned], evidence], ignore_index=True
    )
    combined_audit = pd.concat(
        [
            prior_audit.loc[
                ~prior_audit["candidate_id"].isin(prior_owned_ids | owned_ids)
            ],
            audit,
        ],
        ignore_index=True,
    )
    for label, frame in (("evidence", combined_evidence), ("audit", combined_audit)):
        if frame["candidate_id"].duplicated().any():
            raise ValueError(f"combined {label} candidate IDs are not unique")

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "evidence": output_dir / "morphology_rule_unlock_evidence_20260813.csv",
        "audit": output_dir / "morphology_rule_unlock_manual_audit_20260813.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260813.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260813.csv",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")
    summary: dict[str, object] = {
        "contract": "morphology_support2_rule_unlock_checkpoint_v1",
        "created_at_utc": CREATED_AT,
        "accepted_rows": len(evidence),
        "species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "quality_counts": evidence["evidence_quality"].value_counts().to_dict(),
        "trait_counts": evidence["trait_name"].value_counts().to_dict(),
        "source_lineages": int(evidence["source_lineage"].nunique()),
        "theoretical_unique_rule_cells": 106,
        "theoretical_only_not_formal_coverage": True,
        "excluded_as_duplicate_of_pngtrees_package": ["Gonystylus macrophyllus|flower_size_class"],
        "explicitly_rejected": ["Sloanea dasycarpa|floral_symmetry|family_description"],
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
    manifest_path = output_dir / "morphology_rule_unlock_manifest_20260813.json"
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
