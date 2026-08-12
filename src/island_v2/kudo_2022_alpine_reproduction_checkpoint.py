"""Freeze species-level reproductive evidence from Kudo's alpine dataset.

The source is a single peer-reviewed community dataset.  Its Table 1 reports
self-compatibility indices and a trait-specific mating-system classification
for 46 taxa.  This checkpoint keeps the original table row as the supporting
excerpt, rejects unknown and infraspecific rows, and uses a previously frozen
two-backbone synonym only for ``Loiseleuria procumbens``.
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

CREATED_AT = "2026-08-12T17:50:00Z"
REVIEWER = "Codex Kudo 2022 species-table audit"
SOURCE_GROUP = "kudo_2022_alpine_reproduction_checkpoint_20260813"
SOURCE_URL = (
    "https://esj-journals.onlinelibrary.wiley.com/doi/"
    "10.1111/1440-1703.12314"
)
SOURCE_LINEAGE = "doi:10.1111/1440-1703.12314"
TABLE_SHA256 = "24737a195d4aca322275119f53b949a9adef070edef7c64a176afffe0e0040fe"

# Species, family, pollinator, populations, observation years, SCI, category.
# Values reproduce the rendered Table 1, including source spelling.
TABLE_ROWS = [
    ("Bupleurum triradiatum", "Apiaceae", "Fly", "1", "1", "Unknown", "Unknown"),
    ("Peucedanum multivittatum", "Apiaceae", "Fly", "8", "3–13", "0.18", "Obligate outcrosser"),
    ("Tilingia ajanensis", "Apiaceae", "Fly", "4", "1", "Unknown", "Unknown"),
    ("Arnica unalaschcensis", "Asteraceae", "Fly", "2", "1", "0.01", "Self-incompatibility"),
    ("Saussurea riederi ssp. yezoensis", "Asteraceae", "Bee", "1", "1", "0.04", "Self-incompatibility"),
    ("Saussurea yanagisawae", "Asteraceae", "Bee", "1", "3", "1.00", "Autonomous selfer"),
    ("Senecio kawakamii", "Asteraceae", "Fly", "1", "1", "Unknown", "Unknown"),
    ("Solidago virgaurea ssp. leiocarpa", "Asteraceae", "Fly", "4", "1–3", "0.00", "Obligate outcrosser"),
    ("Campanula chamissonis", "Campanulaceae", "Bee", "1", "3", "0.00", "Self-incompatibility"),
    ("Weigela middendorffiana", "Caprifoliaceae", "Bee", "2", "3", "0.00", "Self-incompatibility"),
    ("Parnassia palustris", "Celastraceae", "Fly", "2", "1–2", "0.00", "Obligate outcrosser"),
    ("Diapensia lapponica", "Diapensiaceae", "Mixture", "2", "3–6", "0.00", "Self-incompatibility"),
    ("Arcterica nana", "Ericaceae", "Bee", "2", "2–4", "0.00", "Self-incompatibility"),
    ("Arctous alpinus", "Ericaceae", "Bee", "2", "7–8", "0.00", "Self-incompatibility"),
    ("Bryanthus gmelinii", "Ericaceae", "Mixture", "3", "1–9", "0.07", "Obligate outcrosser"),
    ("Harrimanella stelleriana", "Ericaceae", "Mixture", "1", "2", "0.00", "Obligate outcrosser"),
    ("Loiseleuria procumbens", "Ericaceae", "Mixture", "3", "2–9", "0.00", "Self-incompatibility"),
    ("Phyllodoce aleutica", "Ericaceae", "Bee", "7", "1–4", "0.47", "Mixed mating"),
    ("Phyllodoce caerulea", "Ericaceae", "Bee", "1", "4", "0.25", "Mixed mating"),
    ("Phyllodoce caerulea var. yezoensis", "Ericaceae", "Bee", "4", "1–2", "0.27", "Mixed mating"),
    ("Rhododendron aureum", "Ericaceae", "Bee", "7", "2–16", "0.12", "Obligate outcrosser"),
    ("Rhododendron camtschaticum", "Ericaceae", "Bee", "1", "9", "0.03", "Self-incompatibility"),
    ("Rhododendron diversipilosum", "Ericaceae", "Fly", "1", "2", "0.07", "Obligate outcrosser"),
    ("Rhododendron subarcticum", "Ericaceae", "Fly", "2", "1–6", "0.03", "Self-incompatibility"),
    ("Vaccinium ovalifolium", "Ericaceae", "Bee", "2", "2", "0.32", "Mixed mating"),
    ("Vaccinium uliginosum", "Ericaceae", "Bee", "2", "5–6", "0.00", "Obligate ourcrosser"),
    ("Vaccinium vitis-idaea", "Ericaceae", "Bee", "2", "7–8", "0.04", "Self-incompatibility"),
    ("Gentiana nipponica", "Gentianaceae", "Mixture", "4", "1–4", "0.25", "Mixed mating"),
    ("Geranium erianthum", "Geraniaceae", "Mixture", "2", "1", "0.00", "Obligate outcrosser"),
    ("Hypericum kamtschaticum", "Hypericaceae", "Fly", "1", "1", "0.00", "Obligate outcrosser"),
    ("Veratrum album ssp. oxysepalum", "Melanthiaceae", "Mixture", "1", "2", "0.00", "Self-incompatibility"),
    ("Nephrophyllidium crista-galli", "Menyanthaceae", "Fly", "3", "1–2", "0.00", "Self-incompatibility"),
    ("Platanthera chorisiana", "Orchidaceae", "Fly", "1", "1", "Unknown", "Unknown"),
    ("Platanthera tipuloides var. sororia", "Orchidaceae", "Moth", "2", "1", "0.00", "Self-incompatibility"),
    ("Pedicularis chamissonis", "Orobanchaceae", "Bee", "5", "1–6", "0.12", "Obligate outcrosser"),
    ("Pennelianthus frutescens", "Plantaginaceae", "Bee", "1", "3", "0.00", "Self-incompatibility"),
    ("Veronica stelleri var. longistyla", "Plantaginaceae", "Mixture", "2", "3", "0.02", "Obligate outcrosser"),
    ("Primura cuneifolia", "Primulaceae", "Bee", "6", "1–4", "0.00", "Self-incompatibility"),
    ("Anemone narcissiflora ssp. sachalinensis", "Ranunculaceae", "Fly", "2", "3–6", "0.13", "Obligate outcrosser"),
    ("Anemone yezoensis", "Ranunculaceae", "Fly", "1", "3", "0.11", "Obligate outcrosser"),
    ("Trollius riederianus", "Ranunculaceae", "Mixture", "1", "1", "0.06", "Obligate outcrosser"),
    ("Potentilla matsumurae", "Rosaceae", "Fly", "7", "1–5", "0.04", "Obligate outcrosser"),
    ("Sieversia pentapetala", "Rosaceae", "Mixture", "5", "2–3", "0.05", "Self-incompatibility"),
    ("Spiraea betulifolia var. aemiliana", "Rosaceae", "Fly", "2", "2–3", "0.01", "Self-incompatibility"),
    ("Patrinia sibirica", "Valerianaceae", "Fly", "1", "3", "0.05", "Obligate outcrosser"),
    ("Hemerocallis middendorffii", "Xanthorrhoeaceae", "Butterfly", "1", "1", "Unknown", "Unknown"),
]

NORMALIZATION = {
    "Self-incompatibility": ("self_incompatibility", "SI"),
    "Autonomous selfer": ("autonomous_selfing_capacity", "autonomous"),
    "Mixed mating": ("mating_system", "mixed_mating"),
    "Obligate outcrosser": ("mating_system", "predominantly_outcrossing"),
    "Obligate ourcrosser": ("mating_system", "predominantly_outcrossing"),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(_canonical_file_bytes(path)).hexdigest()


def _row_excerpt(row: tuple[str, ...]) -> str:
    return "\t".join(row)


def _selected(
    master_family: dict[str, str], synonym_csv: Path
) -> tuple[list[dict[str, str]], pd.DataFrame]:
    synonyms = pd.read_csv(synonym_csv, dtype=str).fillna("")
    synonym = synonyms.loc[
        synonyms["submitted_name"].eq("Loiseleuria procumbens")
        & synonyms["accepted_species"].eq("Kalmia procumbens")
        & synonyms["resolution_method"].eq("strict_exact_two_backbone_agreement")
    ]
    if len(synonym) != 1:
        raise ValueError("missing frozen two-backbone Loiseleuria synonym")

    evidence: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []
    for row in TABLE_ROWS:
        submitted, source_family, pollinator, populations, years, sci, category = row
        accepted = "Kalmia procumbens" if submitted == "Loiseleuria procumbens" else submitted
        method = (
            "exact_synonym"
            if submitted == "Loiseleuria procumbens"
            else "accepted_name_exact"
        )
        reason = "accepted"
        if category == "Unknown":
            reason = "no_trait_statement"
        elif " ssp. " in submitted or " var. " in submitted:
            reason = "rank_below_species_rejected"
        elif submitted == "Primura cuneifolia":
            reason = "source_name_typo_not_fuzzily_matched"
        elif accepted not in master_family:
            reason = "no_strict_target_master_match"
        elif master_family[accepted] != source_family:
            reason = "family_conflict_rejected"

        audit_rows.append(
            {
                "submitted_name": submitted,
                "accepted_species": accepted if reason == "accepted" else "",
                "source_family": source_family,
                "master_family": master_family.get(accepted, ""),
                "pollinator": pollinator,
                "populations": populations,
                "observation_years": years,
                "self_compatibility_index": sci,
                "source_category": category,
                "decision": "accept" if reason == "accepted" else "reject",
                "decision_reason": reason,
                "name_match_method": method if reason == "accepted" else "",
                "source_excerpt": _row_excerpt(row),
            }
        )
        if reason != "accepted":
            continue

        trait, value = NORMALIZATION[category]
        evidence_row = _evidence_row(
            species=accepted,
            trait=trait,
            value=value,
            raw_value=f"SCI={sci}; {category}",
            quality="high",
            provider="Kudo 2022 Ecological Research alpine pollination dataset",
            url=SOURCE_URL,
            title=(
                "Outcrossing syndrome in alpine plants: Implications for "
                "flowering phenology and pollination success"
            ),
            citation=(
                "Kudo, G. (2022), Ecological Research 37:288-300, Table 1, "
                "DOI 10.1111/1440-1703.12314"
            ),
            excerpt=_row_excerpt(row),
            record_id=(
                "doi:10.1111/1440-1703.12314:table1:"
                + submitted.casefold().replace(" ", "-")
            ),
            lineage=SOURCE_LINEAGE,
            lineage_method="original_peer_reviewed_community_experiment_doi",
            source_tier="A",
            source_type="peer_reviewed_species_level_controlled_pollination_dataset",
            domain="esj-journals.onlinelibrary.wiley.com",
            content_sha256=TABLE_SHA256,
            content_sha256_basis="rendered_table_1_inner_text_utf8",
            retrieved_at_utc=CREATED_AT,
        )
        evidence_row["source_group"] = SOURCE_GROUP
        evidence_row["matched_page_name"] = submitted
        evidence_row["name_match_method"] = method
        evidence_row["evidence_scope"] = (
            "synonym_direct"
            if submitted == "Loiseleuria procumbens"
            else "species_direct"
        )
        evidence_row["name_resolution_lineage"] = (
            "gbif_wfo_two_backbone_snapshot"
            if method == "exact_synonym"
            else "master_accepted_name_exact"
        )
        evidence_row["query"] = "peer_reviewed_community_table_reproductive_acquisition"
        evidence.append(evidence_row)

    return evidence, pd.DataFrame(audit_rows)


def reviewed_rows(master_csv: Path, synonym_csv: Path) -> list[dict[str, str]]:
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    master_family = master.set_index("accepted_species")["family"].to_dict()
    rows, _ = _selected(master_family, synonym_csv)
    return rows


def build(
    *,
    master_csv: Path,
    synonym_csv: Path,
    output_dir: Path,
    prior_curated_evidence_csv: Path,
    prior_curated_audit_csv: Path,
) -> dict[str, object]:
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    master_family = master.set_index("accepted_species")["family"].to_dict()
    rows, selection_audit = _selected(master_family, synonym_csv)
    evidence = pd.DataFrame(rows, columns=EVIDENCE_COLUMNS).fillna("")
    evidence = evidence.sort_values(
        ["accepted_species", "trait_name", "candidate_id"]
    ).reset_index(drop=True)
    if len(TABLE_ROWS) != 46 or len(selection_audit) != 46:
        raise ValueError("Kudo Table 1 must contain exactly 46 audited rows")
    if len(evidence) != 23:
        raise ValueError(f"expected 23 accepted species-level rows, found {len(evidence)}")
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("Kudo candidate IDs are not unique")

    audit = _audit(evidence)
    audit["reviewer"] = REVIEWER
    audit["reviewed_at_utc"] = CREATED_AT
    audit["decision_reason"] = evidence["trait_name"].map(
        lambda trait: (
            "Accepted from peer-reviewed Table 1 after exact species or frozen "
            "two-backbone synonym matching; the source category maps only to "
            f"{trait}."
        )
    )

    prior_evidence = pd.read_csv(prior_curated_evidence_csv, dtype=str).fillna("")
    prior_audit = pd.read_csv(prior_curated_audit_csv, dtype=str).fillna("")
    prior_owned = prior_evidence["source_group"].eq(SOURCE_GROUP)
    prior_owned_ids = set(prior_evidence.loc[prior_owned, "candidate_id"].astype(str))
    owned = set(evidence["candidate_id"])
    combined_evidence = pd.concat(
        [prior_evidence.loc[~prior_owned], evidence], ignore_index=True
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

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "readme": output_dir / "README.md",
        "selection_audit": output_dir / "kudo_2022_full_table_audit_20260813.csv",
        "evidence": output_dir / "kudo_2022_evidence_20260813.csv",
        "audit": output_dir / "kudo_2022_manual_audit_20260813.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260813.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260813.csv",
    }
    if not paths["readme"].exists():
        raise ValueError(f"checkpoint README is missing: {paths['readme']}")
    selection_audit.to_csv(paths["selection_audit"], index=False, lineterminator="\n")
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    summary: dict[str, object] = {
        "contract": "kudo_2022_trait_specific_species_table_checkpoint_v1",
        "created_at_utc": CREATED_AT,
        "source_table_rows": len(TABLE_ROWS),
        "accepted_rows": len(evidence),
        "rejected_rows": int(selection_audit["decision"].eq("reject").sum()),
        "selection_reasons": selection_audit["decision_reason"].value_counts().to_dict(),
        "new_species": int(evidence["accepted_species"].nunique()),
        "new_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "quality_counts": evidence["evidence_quality"].value_counts().to_dict(),
        "trait_counts": evidence["trait_name"].value_counts().to_dict(),
        "source_lineages": int(evidence["source_lineage"].nunique()),
        "source_table_sha256": TABLE_SHA256,
        "inputs": {
            "master_csv": {"path": str(master_csv), "sha256": _sha256(master_csv)},
            "synonym_csv": {"path": str(synonym_csv), "sha256": _sha256(synonym_csv)},
            "prior_curated_evidence_csv": {
                "path": str(prior_curated_evidence_csv),
                "sha256": _sha256(prior_curated_evidence_csv),
            },
            "prior_curated_audit_csv": {
                "path": str(prior_curated_audit_csv),
                "sha256": _sha256(prior_curated_audit_csv),
            },
        },
        "combined": {
            "evidence_rows": len(combined_evidence),
            "audit_rows": len(combined_audit),
        },
        "guardrails": {
            "trait_specific_records": True,
            "family_conflicts_rejected": True,
            "infraspecific_rows_rejected": True,
            "fuzzy_name_matching": False,
            "genus_inference_emitted_here": False,
            "family_inference": False,
            "global_fallback": False,
            "cross_trait_substitution": False,
        },
        "files": {},
    }
    for path in paths.values():
        summary["files"][path.name] = {
            "sha256": _sha256(path),
            "size_bytes": len(_canonical_file_bytes(path)),
        }
    manifest = output_dir / "kudo_2022_manifest_20260813.json"
    manifest.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--synonym-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--prior-curated-evidence-csv", type=Path, required=True)
    parser.add_argument("--prior-curated-audit-csv", type=Path, required=True)
    print(json.dumps(build(**vars(parser.parse_args())), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
