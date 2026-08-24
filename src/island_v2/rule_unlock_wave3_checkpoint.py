"""Freeze a small, high-information species-direct evidence checkpoint.

Each record was selected from the current ``support == 2`` acquisition queue
or upgrades an existing Validated Low value.  This module only records direct
trait statements.  The shared all-evidence rebuild remains responsible for
lineage deduplication, conflict handling, genus x trait-name validation, and
strict three-axis coverage.
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

CREATED_AT = "2026-08-13T03:40:00Z"
REVIEWER = "Codex source-backed rule-unlock wave-3 audit"
SOURCE_GROUP = "rule_unlock_wave3_checkpoint_20260813"


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
) -> dict[str, str]:
    row = _evidence_row(
        species=species,
        trait=trait,
        value=value,
        raw_value=raw_value,
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
    )
    row["source_group"] = SOURCE_GROUP
    row["language"] = language
    row["name_resolution_lineage"] = "master_accepted_name_exact"
    row["wild_cultivated_cultivar_status"] = (
        "wild_or_species_level_statement_not_cultivar_limited"
    )
    row["query"] = "current_support_2_rule_unlock_exact_species_acquisition"
    return row


def reviewed_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    bidens_excerpt = (
        "The results from four treatments (open pollination, bagged capitulum, "
        "decapitated capitulum, and pollen supplement) indicated that, with the "
        "exception of Bidens pilosa var. radiata, all the other five Bidens "
        "taxa; namely, B. pilosa var. pilosa, B. pilosa var. minor, B. frondosa, "
        "B. bipinnata, and B. biternata, were self-fertile with a high seed-set "
        "in bagged capitula, and a high autofertility index."
    )
    bidens_common = {
        "trait": "autonomous_selfing_capacity",
        "value": "autonomous",
        "raw_value": "high seed-set in bagged capitula and high autofertility index",
        "excerpt": bidens_excerpt,
        "quality": "high",
        "provider": "Hao et al. 2018 Weed Biology and Management",
        "url": (
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
            "query=EXT_ID:IND605920380&format=json&resultType=core"
        ),
        "title": "Breeding systems and seed production for six weedy taxa of Bidens",
        "citation": (
            "Hao, Bhattacharya, Ma & Wang (2018), Weed Biology and Management "
            "18:41-49, DOI 10.1111/wbm.12142"
        ),
        "lineage": "doi:10.1111/wbm.12142",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "peer_reviewed_primary_bagging_experiment",
        "domain": "ebi.ac.uk",
        "content_sha256": (
            "ecddf41dbebe2afc31aea1328265e2f814bbabe7e572b29defe650be3443a8e0"
        ),
        "content_sha256_basis": "downloaded_europe_pmc_core_api_json_bytes",
    }
    for species in ("Bidens bipinnata", "Bidens biternata", "Bidens frondosa"):
        slug = species.casefold().replace(" ", "-")
        rows.append(
            _row(
                species=species,
                record_id=f"doi:10.1111/wbm.12142:abstract:{slug}:autofertility",
                **bidens_common,
            )
        )

    rows.extend(
        [
            _row(
                species="Coccothrinax barbadensis",
                trait="flower_size_class",
                value="small",
                raw_value="flowers about 2 mm in diameter",
                excerpt=(
                    "Květenství je kratší než listy, vzpřímené nebo obloukovitě "
                    "skloněné, 25–45 cm dlouhé; květy jsou nažloutlé, redukované, "
                    "dosahují asi 2 mm v průměru, tyčinek je 9–12."
                ),
                quality="medium",
                provider="BOTANY.cz expert botanical atlas",
                url="https://botany.cz/cs/coccothrinax-barbadensis/",
                title="COCCOTHRINAX BARBADENSIS (Lodd. ex Mart.) Becc.",
                citation=(
                    "Hoskovec, L. (2023), BOTANY.cz species treatment for "
                    "Coccothrinax barbadensis"
                ),
                record_id="botany-cz:coccothrinax-barbadensis:flower-size",
                lineage="provider_treatment:botany.cz:coccothrinax-barbadensis",
                lineage_method="canonical_expert_species_treatment",
                source_tier="B",
                source_type="expert_authored_botanical_atlas_species_treatment",
                domain="botany.cz",
                content_sha256=(
                    "c55c8eb81b0321ec932b0e81ab8b3de1c3ec981e793cde90cce2f53aabb8ef57"
                ),
                content_sha256_basis="downloaded_complete_species_page_html_bytes",
                language="cs",
            ),
            _row(
                species="Coccothrinax barbadensis",
                trait="inflorescence_display",
                value="raceme_spike_panicle",
                raw_value="inflorescence branched to two orders",
                excerpt=(
                    "Flowers and fruits: Inflorescence is shorter than the leaves, "
                    "to 1.5 m long, branched to two orders, with up to 10 primary "
                    "branches."
                ),
                quality="medium",
                provider="Palmpedia archived palm species treatment",
                url=(
                    "https://web.archive.org/web/20200119123253id_/"
                    "http://www.palmpedia.net:80/wiki/Coccothrinax_barbadensis"
                ),
                title="Coccothrinax barbadensis - Palmpedia",
                citation="Palmpedia species treatment, archived 19 January 2020",
                record_id="palmpedia:Coccothrinax_barbadensis:inflorescence-display",
                lineage="provider_treatment:palmpedia:Coccothrinax_barbadensis",
                lineage_method="canonical_archived_species_treatment",
                source_tier="C",
                source_type="specialist_palm_encyclopedia_species_treatment",
                domain="palmpedia.net",
                content_sha256=(
                    "e76f017290fb712a050536c9044af474003d225ca519a216a39aa73dddf7c89e"
                ),
                content_sha256_basis="downloaded_wayback_species_page_html_bytes",
            ),
            _row(
                species="Bambusa tulda",
                trait="inflorescence_display",
                value="raceme_spike_panicle",
                raw_value="leafless panicle",
                excerpt=(
                    "Inflorescence a leafless panicle with a branching pattern "
                    "similar to that of a vegetative culm."
                ),
                quality="high",
                provider="Bhattacharya et al. 2006 Annals of Botany",
                url="https://pmc.ncbi.nlm.nih.gov/articles/PMC2803566/",
                title=(
                    "Morphological and Molecular Characterization of Bambusa "
                    "tulda with a Note on Flowering"
                ),
                citation=(
                    "Bhattacharya, Das, Bar & Pal (2006), Annals of Botany "
                    "98:529-535, DOI 10.1093/aob/mcl143"
                ),
                record_id="doi:10.1093/aob/mcl143:results:inflorescence",
                lineage="doi:10.1093/aob/mcl143",
                lineage_method="original_primary_article_doi",
                source_tier="A",
                source_type="peer_reviewed_primary_morphological_study",
                domain="pmc.ncbi.nlm.nih.gov",
                content_sha256=(
                    "2670b702fa42594d8c82b55f755b3e60a7500d16c15e2c4a52cd069eb981813d"
                ),
                content_sha256_basis="downloaded_pmc_fulltext_html_bytes",
            ),
            _row(
                species="Bulbophyllum tseanum",
                trait="self_incompatibility",
                value="SI",
                raw_value="self-incompatible",
                excerpt=(
                    "Attempts at hand-pollinating this species in the nursery have "
                    "demonstrated that it is self-incompatible – that is, unable "
                    "to produce seeds through self-pollination."
                ),
                quality="medium",
                provider="Kadoorie Farm and Botanic Garden",
                url="https://www.kfbg.org/images/ar/report/kfbg-annual-report-2020.pdf",
                title="Kadoorie Farm and Botanic Garden Annual Report 2020",
                citation=(
                    "Kadoorie Farm and Botanic Garden (2020), Annual Report, p. 31"
                ),
                record_id="kfbg:annual-report-2020:p31:bulbophyllum-tseanum:si",
                lineage=(
                    "provider_report:kfbg:annual-report-2020:bulbophyllum-tseanum"
                ),
                lineage_method="official_botanic_garden_report_species_experiment",
                source_tier="A",
                source_type="botanic_garden_hand_pollination_report",
                domain="kfbg.org",
                content_sha256=(
                    "be9acfa19a33c7545143d991c09d395b7f36c9b8a8b810e1c0c4516f55fba6f6"
                ),
                content_sha256_basis="downloaded_official_annual_report_pdf_bytes",
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
            "Accepted after exact target-master species identity, original-page "
            "quote, trait-specific ontology, cultivar scope, lineage and content "
            f"fingerprint review; mapping retained only for {row['trait_name']}."
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
    prior_curated_evidence_csv: Path,
    prior_curated_audit_csv: Path,
) -> dict[str, object]:
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    master_family = master.set_index("accepted_species")["family"].to_dict()
    expected_families = {
        "Bambusa tulda": "Poaceae",
        "Bidens bipinnata": "Asteraceae",
        "Bidens biternata": "Asteraceae",
        "Bidens frondosa": "Asteraceae",
        "Bulbophyllum tseanum": "Orchidaceae",
        "Coccothrinax barbadensis": "Arecaceae",
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
    if len(evidence) != 7:
        raise ValueError(f"expected 7 reviewed trait rows, found {len(evidence)}")
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("wave-3 candidate IDs are not unique")
    audit = _review_audit(evidence)

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
        "evidence": output_dir / "rule_unlock_wave3_evidence_20260813.csv",
        "audit": output_dir / "rule_unlock_wave3_manual_audit_20260813.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260813.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260813.csv",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(
        paths["combined_evidence"], index=False, lineterminator="\n"
    )
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    summary: dict[str, object] = {
        "contract": "trait_specific_rule_unlock_wave3_individually_reviewed_v1",
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
        "inputs": {
            "master_csv": {
                "path": str(master_csv),
                "sha256": _sha256(master_csv),
                "size_bytes": len(_canonical_file_bytes(master_csv)),
            },
            "prior_curated_evidence_csv": {
                "path": str(prior_curated_evidence_csv),
                "sha256": _sha256(prior_curated_evidence_csv),
                "size_bytes": len(_canonical_file_bytes(prior_curated_evidence_csv)),
            },
            "prior_curated_audit_csv": {
                "path": str(prior_curated_audit_csv),
                "sha256": _sha256(prior_curated_audit_csv),
                "size_bytes": len(_canonical_file_bytes(prior_curated_audit_csv)),
            },
        },
        "combined": {
            "evidence_rows": len(combined_evidence),
            "audit_rows": len(combined_audit),
        },
        "pre_formal_targets": {
            "candidate_genus_trait_rules": 4,
            "maximum_distinct_species_axis_unlocked": 135,
            "validated_low_to_direct_candidates": 1,
            "status": "theoretical_until_shared_all_evidence_rebuild",
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
    for path in paths.values():
        summary["files"][path.name] = {
            "sha256": _sha256(path),
            "size_bytes": len(_canonical_file_bytes(path)),
        }
    manifest = output_dir / "rule_unlock_wave3_manifest_20260813.json"
    manifest.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--prior-curated-evidence-csv", type=Path, required=True)
    parser.add_argument("--prior-curated-audit-csv", type=Path, required=True)
    print(json.dumps(build(**vars(parser.parse_args())), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
