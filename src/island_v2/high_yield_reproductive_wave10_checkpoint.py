"""Freeze the tenth high-yield reproductive evidence checkpoint.

The checkpoint records only reviewed species-direct traits.  Queue potential is
diagnostic: genus inference remains the responsibility of the common
trait-specific all-evidence rebuild and its masked/source-lineage validation.
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

CREATED_AT = "2026-08-14T09:20:00Z"
SOURCE_GROUP = "high_yield_reproductive_wave10_checkpoint_20260814"
OUTPUT_STEM = "high_yield_reproductive_wave10"

CIRSIUM_SOURCE = {
    "species": "Cirsium brevistylum",
    "family": "Asteraceae",
    "trait": "autonomous_selfing_capacity",
    "value": "autonomous",
    "provider": "Biological Invasions primary article",
    "url": "https://doi.org/10.1007/s10530-010-9878-5",
    "title": (
        "Comparing the reproductive success and pollination biology of an "
        "invasive plant to its rare and common native congeners"
    ),
    "citation": (
        "Powell, Krakos and Knight (2011), Biological Invasions 13:905-917, "
        "DOI 10.1007/s10530-010-9878-5"
    ),
    "excerpt": (
        "Native study species include C. fontinale var. fontinale, C. andrewsii, "
        "C. brevistylum, C. occidentale, and C. quercetorum. We compared all "
        "species' reproductive success, insect visitation rate and composition, "
        "autonomous self-pollination, and level of pollen limitation in multiple "
        "populations. The remaining native species set fewer seeds than C. "
        "vulgare without a pollinator."
    ),
    "raw_value": (
        "positive seed set without a pollinator; lower autonomous selfing than "
        "C. vulgare"
    ),
    "record_id": "doi:10.1007/s10530-010-9878-5:cirsium-brevistylum",
    "lineage": "doi:10.1007/s10530-010-9878-5",
    "lineage_method": "original_primary_article_doi",
    "source_type": "peer_reviewed_primary_pollinator_exclusion_experiment",
    "domain": "link.springer.com",
    "content_sha256": (
        "14ca1d67b9e130caea2eb66db25979881386270e747d7ecb444def86bd007bba"
    ),
    "content_sha256_basis": "downloaded_publisher_full_html_bytes",
    "potential": 168,
    "role": "rule_candidate",
}

INGA_SPECIES = ("Inga vera", "Inga striata", "Inga ingoides")
INGA_SOURCE = {
    "family": "Fabaceae",
    "provider": "Universidade Federal de Pernambuco doctoral thesis",
    "url": (
        "https://repositorio.ufpe.br/bitstream/123456789/931/1/arquivo4845_1.pdf"
    ),
    "title": "Biologia reprodutiva de três espécies simpátricas de Inga",
    "citation": (
        "Cruz Neto, Universidade Federal de Pernambuco thesis; experiment later "
        "published as Cruz-Neto et al. (2015), DOI 10.1111/boj.12236"
    ),
    "lineage": "original-study:cruz-neto-ufpe-inga-coimbra",
    "lineage_method": (
        "one_original_field_experiment_shared_by_thesis_and_later_article"
    ),
    "source_type": "doctoral_thesis_primary_controlled_pollination_experiment",
    "domain": "repositorio.ufpe.br",
    "content_sha256": (
        "75f650535c722b9e36fc9f171062d59a95d75601d487011bc10f63ad742a2c5d"
    ),
    "content_sha256_basis": "downloaded_official_university_repository_pdf_bytes",
}


def _inga_rows() -> list[dict[str, object]]:
    autonomy_excerpt = (
        "Autopolinização espontânea 0% (100/0) 0% (100/0) 0% (100/0); "
        "Autopolinização manual 0% (30/0) 0% (30/0) 0% (30/0); "
        "Polinização cruzada 20% (30/6) 6,6% (30/2) 10% (30/3)."
    )
    si_excerpt = (
        "Inga vera, I. striata e I. ingoides são auto-incompatíveis, uma vez "
        "que não houve formação de frutos em nenhuma das três espécies nos "
        "testes de autopolinização manual e espontânea."
    )
    rows: list[dict[str, object]] = []
    for index, species in enumerate(INGA_SPECIES):
        slug = species.casefold().replace(" ", "-")
        common = {
            "species": species,
            **INGA_SOURCE,
            "quality": "high",
            "source_tier": "A",
            "language": "pt",
        }
        rows.append(
            {
                **common,
                "trait": "autonomous_selfing_capacity",
                "value": "absent",
                "excerpt": autonomy_excerpt,
                "raw_value": (
                    "spontaneous self-pollination 0/100; manual self 0/30; "
                    "cross-pollination positive"
                ),
                "record_id": f"ufpe:123456789-931:table3:{slug}:autonomy",
                "potential": 54 if index == 0 else 0,
                "role": "rule_candidate",
            }
        )
        rows.append(
            {
                **common,
                "trait": "self_incompatibility",
                "value": "SI",
                "excerpt": si_excerpt,
                "raw_value": "manual self 0/30; cross-pollination positive",
                "record_id": f"ufpe:123456789-931:table3:{slug}:si",
                "potential": 0,
                "role": "direct_enrichment",
            }
        )
    return rows


SOURCES: tuple[dict[str, object], ...] = (
    {
        **CIRSIUM_SOURCE,
        "quality": "high",
        "source_tier": "A",
        "language": "en",
    },
    *_inga_rows(),
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
                "query": "latest_support2_exact_trait_wave10",
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
        "Accepted after exact master identity, exact trait and value-polarity "
        "review, full-source retrieval, lineage audit, positive controls where "
        "applicable, and cultivar/context screening."
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
        "evidence": output_dir / f"{OUTPUT_STEM}_evidence_20260814.csv",
        "audit": output_dir / f"{OUTPUT_STEM}_manual_audit_20260814.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260814.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260814.csv",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    role_counts = pd.Series([str(source["role"]) for source in SOURCES]).value_counts()
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
        "targeted_genus_trait_rules": 2,
        "pending_rejected_rows": 6,
        "theoretical_queue_cells_before_formal_validation": sum(
            int(source["potential"]) for source in SOURCES
        ),
        "theoretical_only_not_formal_coverage": True,
        "expected_new_direct_species_axis_before_formal_rebuild": 2,
        "combined": {
            "evidence_rows": len(combined_evidence),
            "audit_rows": len(combined_audit),
        },
        "guardrails": {
            "trait_specific_records": True,
            "autonomy_and_self_incompatibility_separate": True,
            "genus_inference_emitted_here": False,
            "family_inference": False,
            "global_fallback": False,
            "n2_formal_inference": False,
            "cross_trait_substitution": False,
            "search_snippet_evidence": False,
            "queue_potential_counted_as_coverage": False,
            "shared_inga_study_counted_as_one_lineage": True,
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
    for path in [*paths.values(), *([readme] if readme.exists() else [])]:
        summary["files"][path.name] = {
            "sha256": _sha256(path),
            "size_bytes": len(_canonical_file_bytes(path)),
        }
    manifest = output_dir / f"{OUTPUT_STEM}_manifest_20260814.json"
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
