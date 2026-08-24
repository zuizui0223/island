"""Freeze a hand-pollination checkpoint for reproductive rule information gain.

The source article tested autogamous, geitonogamous and allogamous hand
pollinations and evaluated fruiting, seed viability and germination.  Every
row below is an exact target-master species and an explicit compatibility
result from that one study.  The rows intentionally share one DOI lineage;
they must not be treated as independent publications.

This module emits species-direct ``self_incompatibility`` evidence only.  It
does not infer autonomous selfing or mating system from compatibility, and it
does not create genus, family, global or n=2 inference.
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

CREATED_AT = "2026-08-14T06:30:00Z"
SOURCE_GROUP = "reproductive_rule_wave4_checkpoint_20260814"
SOURCE_LINEAGE = "doi:10.1093/botlinnean/boae032"
SOURCE_URL = "https://doi.org/10.1093/botlinnean/boae032"
SOURCE_TITLE = (
    "Seed quality and germination performance increase with cross-pollination "
    "in members of subtribe Orchidinae (Orchidaceae)"
)
SOURCE_CITATION = (
    "Bazzicalupo, Masullo, Duffy, Fay & Calevo (2025), Botanical Journal of "
    "the Linnean Society 207:1-9, DOI 10.1093/botlinnean/boae032"
)
SOURCE_EXCERPT = (
    "Our germination rate experiment revealed that A. papilionacea, "
    "H. adriaticum, Ophrys apifera, Ophrys bertolonii, Orchis patens subsp. "
    "brevicornis, Orchis provincialis, and S. vomeracea are self-compatible."
)
SPECIES = (
    "Anacamptis papilionacea",
    "Himantoglossum adriaticum",
    "Ophrys apifera",
    "Ophrys bertolonii",
    "Orchis provincialis",
    "Serapias vomeracea",
)


def reviewed_rows() -> list[dict[str, str]]:
    excerpt_sha256 = hashlib.sha256(SOURCE_EXCERPT.encode()).hexdigest()
    rows: list[dict[str, str]] = []
    for species in SPECIES:
        evidence = _evidence_row(
            species=species,
            trait="self_incompatibility",
            value="SC",
            quality="high",
            provider="Oxford Academic original controlled-pollination article",
            url=SOURCE_URL,
            title=SOURCE_TITLE,
            citation=SOURCE_CITATION,
            excerpt=SOURCE_EXCERPT,
            record_id=(
                "doi:10.1093/botlinnean/boae032:germination-self-compatibility:"
                + species.casefold().replace(" ", "-")
            ),
            lineage=SOURCE_LINEAGE,
            lineage_method="original_primary_article_doi",
            source_tier="A",
            source_type="peer_reviewed_controlled_pollination_experiment",
            domain="academic.oup.com",
            content_sha256=excerpt_sha256,
            content_sha256_basis="verified_publisher_fulltext_excerpt_utf8_bytes",
            retrieved_at_utc=CREATED_AT,
            raw_value="confirmed as self-compatible by germination-rate experiment",
        )
        evidence["source_group"] = SOURCE_GROUP
        evidence["query"] = "support2_self_compatibility_primary_hand_pollination"
        evidence["name_resolution_lineage"] = "master_accepted_name_exact"
        evidence["wild_cultivated_cultivar_status"] = (
            "wild_in_situ_northwestern_italy_not_cultivar_limited"
        )
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
    missing = sorted(set(SPECIES) - set(master_identity))
    if missing:
        raise ValueError(f"checkpoint species missing from master: {missing}")
    family_conflicts = {
        species: master_identity[species]
        for species in SPECIES
        if master_identity[species] != "Orchidaceae"
    }
    if family_conflicts:
        raise ValueError(f"checkpoint family conflicts: {family_conflicts}")

    evidence = pd.DataFrame(reviewed_rows(), columns=EVIDENCE_COLUMNS).fillna("")
    if evidence[["accepted_species", "trait_name"]].duplicated().any():
        raise ValueError("checkpoint species-trait pairs must be unique")
    if not evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all():
        raise ValueError("every checkpoint row requires a 64-character content hash")

    audit = _audit(evidence)
    audit["reviewer"] = "Codex primary reproductive evidence audit"
    audit["reviewed_at_utc"] = CREATED_AT
    audit["decision_reason"] = (
        "Accepted after exact target-master species identity, in-situ plant "
        "material, hand-pollination methods, germination-based compatibility, "
        "source lineage, value polarity and cultivar status were rechecked."
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
        "evidence": output_dir / "reproductive_rule_wave4_evidence_20260814.csv",
        "audit": output_dir / "reproductive_rule_wave4_manual_audit_20260814.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260814.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260814.csv",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    summary: dict[str, object] = {
        "contract": "reproductive_primary_hand_pollination_checkpoint_v1",
        "created_at_utc": CREATED_AT,
        "accepted_rows": len(evidence),
        "species_trait": len(evidence),
        "quality_counts": evidence["evidence_quality"].value_counts().to_dict(),
        "trait_counts": evidence["trait_name"].value_counts().to_dict(),
        "source_lineages": int(evidence["source_lineage"].nunique()),
        "theoretical_ophrys_rule_cells_before_formal_validation": 52,
        "theoretical_only_not_formal_coverage": True,
        "excluded_infraspecific_result": "Orchis patens subsp. brevicornis",
        "combined": {
            "evidence_rows": len(combined_evidence),
            "audit_rows": len(combined_audit),
        },
        "guardrails": {
            "trait_specific_records": True,
            "shared_paper_counted_as_one_lineage": True,
            "genus_inference_emitted_here": False,
            "family_inference": False,
            "global_fallback": False,
            "n2_formal_inference": False,
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
    readme_path = output_dir / "README.md"
    hash_paths = list(paths.values())
    if readme_path.exists():
        hash_paths.append(readme_path)
    for path in hash_paths:
        summary["files"][path.name] = {
            "sha256": _sha256(path),
            "size_bytes": len(_canonical_file_bytes(path)),
        }
    manifest_path = output_dir / "reproductive_rule_wave4_manifest_20260814.json"
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
