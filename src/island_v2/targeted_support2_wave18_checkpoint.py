"""Build the Wave 18 reviewed reproductive-evidence checkpoint.

Wave 18 adds only species-direct rows with an exact fixed-master identity,
an exact source quotation, and an auditable source lineage.  It targets four
support=2 genus x trait rules but never creates genus inference itself; the
shared all-evidence rebuild remains the sole producer of Validated Low rows.
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
from island_v2.targeted_support2_wave15_checkpoint import (
    _row as _wave15_row,
)
from island_v2.targeted_support2_wave15_checkpoint import (
    build_audit as _wave15_build_audit,
)

SOURCE_GROUP = "targeted_support2_wave18_checkpoint_20260820"
CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/targeted_support2_wave18_checkpoint_20260820"
)
PRIOR = Path("data/v2/staging/traits/open_web_pilot/targeted_synonym_wave17_checkpoint_20260820")
RETRIEVED_AT = "2026-08-20T16:30:00Z"


def _row(*args: str, **kwargs: str) -> dict[str, str]:
    row = _wave15_row(*args, **kwargs)
    row["source_group"] = SOURCE_GROUP
    row["retrieved_at_utc"] = RETRIEVED_AT
    return row


def primary_rows() -> list[dict[str, str]]:
    """Return individually reviewed, exact-species reproductive evidence."""
    spathoglottis_excerpt = (
        "The much-reduced lip of S. aurea and S. microchilina, the "
        "underdeveloped rostellum, and the swollen ovary with developed pollen "
        "tubes during the floral bud stage are evidence of a self-pollination "
        "strategy. In this study, each of the individuals of S. microchilina "
        "from Borneo was observed to be completely cleistogamous or occasionally "
        "geitonogamous. The flower does not open fully, and fruit set percentage "
        "is very high; a sign of a successful self-pollination strategy. Likewise, "
        "similar observation was noticed among individuals of S. aurea from "
        "different populations in Peninsular Malaysia. A complete cleistogamous "
        "form of S. aurea is usually stunted or smaller in size, in comparison to "
        "the geitonogamous or insect-pollinated S. aurea."
    )
    echinochloa_excerpt = (
        "Late season florets are cleistogamous (not opening). It is primarily "
        "self-pollinating and self-compatible. Some degree of outcrossing "
        "recorded which was facilitated by wind pollination."
    )
    rows = [
        _row(
            "Alternanthera pungens",
            "reproductive_assurance",
            "self_incompatibility",
            "ISI fruit=1.2273 SC; seed=1.2273 SC; composite NAG/PSSP/PX/SC",
            "SC",
            "high",
            "Open Journal of Ecology",
            "https://www.scirp.org/pdf/oje_2022071815180097.pdf",
            "Sexual and Breeding Systems in a Xerophytic Shrubland",
            "Ramirez 2022. Open Journal of Ecology 12:434-482. DOI 10.4236/oje.2022.127025.",
            "Alternanthera pungens 0.2171 PSSP 0.2171 PSSP PSSP 0.2664 PX 0.2664 PX PX 1.2273 SC 1.2273 SC SC; NAG, PSSP, PX, SC.",
            "doi:10.4236/oje.2022.127025",
            "A",
            "primary_article_controlled_pollination_table",
            "en",
            '"Alternanthera pungens" self compatibility breeding system',
            content_sha256=("d3f46146de8bb7d8d6128050f20e76d88b2427ae426930035532b4a2d97b804b"),
            content_sha256_basis="retrieved_original_article_pdf_bytes",
        ),
        _row(
            "Lolium canariense",
            "reproductive_assurance",
            "autonomous_selfing_capacity",
            "preferentially autogamous; outcrossing rate 0.08-0.16",
            "autonomous",
            "medium",
            "University of Kentucky International Grassland Congress Proceedings",
            "https://uknowledge.uky.edu/igc/1993/session9/4/",
            "Lolium canariense: An Endemic Ryegrass from Canary Islands",
            "Oliveira, Charmet, Balfourier & Ravel 1993. IGC Proceedings.",
            "Preliminary results show that L. canariense is preferentially autogamous, with an outcrossing rate of 0.08-0.16.",
            "url:https://uknowledge.uky.edu/igc/1993/session9/4/",
            "B",
            "university_hosted_population_genetics_proceedings",
            "en",
            '"Lolium canariense" autogamous outcrossing rate',
            content_sha256=("c61ebcfbd939ec555e9ed3bb55334c790e2827e3cf298bd5d1e77b0ed16dfbc2"),
            content_sha256_basis="retrieved_university_repository_html_bytes",
        ),
    ]

    spathoglottis_common = {
        "quality": "medium",
        "provider": "Forests / MDPI",
        "url": "https://www.mdpi.com/1999-4907/14/5/940",
        "title": (
            "Morphological Systematics of Spathoglottis Blume "
            "(Orchidaceae: Collabieae) in Peninsular Malaysia and Borneo"
        ),
        "citation": "Nordin et al. 2023. Forests 14:940. DOI 10.3390/f14050940.",
        "excerpt": spathoglottis_excerpt,
        "lineage": "doi:10.3390/f14050940",
        "tier": "A",
        "source_type": "peer_reviewed_morphological_systematics_observation",
        "language": "en",
        "content_sha256": ("8d9854d45b8975df7a0c56fdaa3b58e736d56ff962b8114db68ecc540446a1a8"),
    }
    for species in ("Spathoglottis aurea", "Spathoglottis microchilina"):
        rows.extend(
            [
                _row(
                    species,
                    "reproductive_assurance",
                    "autonomous_selfing_capacity",
                    "self-pollination in floral bud; developed pollen tubes and high fruit set",
                    "autonomous",
                    spathoglottis_common["quality"],
                    spathoglottis_common["provider"],
                    spathoglottis_common["url"],
                    spathoglottis_common["title"],
                    spathoglottis_common["citation"],
                    spathoglottis_common["excerpt"],
                    spathoglottis_common["lineage"],
                    spathoglottis_common["tier"],
                    spathoglottis_common["source_type"],
                    spathoglottis_common["language"],
                    f'"{species}" self-pollination cleistogamous',
                    content_sha256=spathoglottis_common["content_sha256"],
                    content_sha256_basis="retrieved_original_article_pdf_bytes",
                ),
                _row(
                    species,
                    "reproductive_assurance",
                    "cleistogamy",
                    "cleistogamous forms with occasional geitonogamous or insect-pollinated forms",
                    "facultative",
                    spathoglottis_common["quality"],
                    spathoglottis_common["provider"],
                    spathoglottis_common["url"],
                    spathoglottis_common["title"],
                    spathoglottis_common["citation"],
                    spathoglottis_common["excerpt"],
                    spathoglottis_common["lineage"],
                    spathoglottis_common["tier"],
                    spathoglottis_common["source_type"],
                    spathoglottis_common["language"],
                    f'"{species}" cleistogamous',
                    content_sha256=spathoglottis_common["content_sha256"],
                    content_sha256_basis="retrieved_original_article_pdf_bytes",
                ),
            ]
        )

    echinochloa_common = {
        "quality": "medium",
        "provider": "International Journal of Current Research",
        "url": "https://journalcra.com/sites/default/files/issue-pdf/45982.pdf",
        "title": (
            "Origin, domestication, taxonomy, botanical description, genetics and "
            "cytogenetics, genetic diversity and breeding of barnyard millet"
        ),
        "citation": (
            "Swamy 2023. International Journal of Current Research 15(9):25839-25864. "
            "DOI 10.24941/ijcr.45982.09.2023."
        ),
        "excerpt": echinochloa_excerpt,
        "lineage": "doi:10.24941/ijcr.45982.09.2023",
        "tier": "B",
        "source_type": "peer_reviewed_crop_species_review",
        "language": "en",
        "content_sha256": ("da4b017aa5c870104ef3037dc7d63a0ee247b89097c8f4df414db99c1ad3bbe9"),
        "status": "cultivated_species_level_not_named_cultivar",
    }
    for trait, raw, value in (
        ("self_incompatibility", "self-compatible", "SC"),
        ("mating_system", "primarily self-pollinating", "predominantly_selfing"),
        ("cleistogamy", "late season florets are cleistogamous", "facultative"),
    ):
        rows.append(
            _row(
                "Echinochloa frumentacea",
                "reproductive_assurance",
                trait,
                raw,
                value,
                echinochloa_common["quality"],
                echinochloa_common["provider"],
                echinochloa_common["url"],
                echinochloa_common["title"],
                echinochloa_common["citation"],
                echinochloa_common["excerpt"],
                echinochloa_common["lineage"],
                echinochloa_common["tier"],
                echinochloa_common["source_type"],
                echinochloa_common["language"],
                f'"Echinochloa frumentacea" {raw}',
                status=echinochloa_common["status"],
                content_sha256=echinochloa_common["content_sha256"],
                content_sha256_basis="retrieved_original_review_pdf_bytes",
            )
        )
    if len(rows) != 9:
        raise AssertionError("Wave18 must define exactly nine direct rows")
    return rows


def build_audit(evidence: pd.DataFrame) -> pd.DataFrame:
    audit = _wave15_build_audit(evidence)
    audit["reviewer"] = "Codex Wave18 individual source audit"
    audit["reviewed_at_utc"] = RETRIEVED_AT
    return audit.loc[:, AUDIT_COLUMNS]


def build_checkpoint(output_dir: Path = CHECKPOINT) -> dict[str, Any]:
    evidence = pd.DataFrame(primary_rows(), columns=EVIDENCE_COLUMNS).sort_values(
        ["accepted_species", "trait_name", "source_lineage"], kind="stable"
    )
    evidence = evidence.reset_index(drop=True)
    audit = build_audit(evidence)
    if len(evidence) != 9 or len(audit) != 9:
        raise ValueError("Wave18 must contain exactly nine individually reviewed rows")
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("Wave18 candidate IDs must be unique")
    if evidence.duplicated(["accepted_species", "trait_name"]).any():
        raise ValueError("Wave18 species x trait pairs must be unique")

    prior_evidence = pd.read_csv(
        PRIOR / "combined_curated_evidence_20260820.csv", dtype=str
    ).fillna("")
    prior_audit = pd.read_csv(
        PRIOR / "combined_curated_manual_audit_20260820.csv", dtype=str
    ).fillna("")
    combined_evidence = pd.concat([prior_evidence, evidence], ignore_index=True)
    combined_audit = pd.concat([prior_audit, audit], ignore_index=True)
    if combined_evidence["candidate_id"].duplicated().any():
        raise ValueError("combined evidence candidate IDs must be unique")
    if combined_audit["candidate_id"].duplicated().any():
        raise ValueError("combined audit candidate IDs must be unique")

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "evidence": output_dir / "targeted_support2_wave18_evidence_20260820.csv",
        "audit": output_dir / "targeted_support2_wave18_manual_audit_20260820.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260820.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260820.csv",
        "manifest": output_dir / "source_acquisition_manifest_wave18.json",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    manifest = {
        "checkpoint": SOURCE_GROUP,
        "built_at_utc": RETRIEVED_AT,
        "baseline_formal_run_id": 32339543372,
        "accepted_evidence_rows": len(evidence),
        "accepted_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "accepted_species": int(evidence["accepted_species"].nunique()),
        "recorded_queries": int(evidence["query"].nunique()),
        "formal_search_api_queries": 0,
        "search_cost_usd": 0.0,
        "targeted_support2_rules": [
            "Alternanthera|self_incompatibility",
            "Echinochloa|self_incompatibility",
            "Lolium|autonomous_selfing_capacity",
            "Spathoglottis|autonomous_selfing_capacity",
        ],
        "theoretical_rule_cells_touched": 118,
        "guardrails": {
            "search_snippet_as_evidence": False,
            "family_inference": False,
            "global_fallback": False,
            "min_species_two_production": False,
            "cross_trait_substitution": False,
            "genus_axis_only_join": False,
            "pollen_vector_or_reward_mixed_into_structure": False,
            "self_fertility_silently_mapped_to_autonomous_selfing": False,
            "paywall_or_access_control_bypassed": False,
            "ambiguous_peperomia_mechanism_promoted": False,
        },
        "output_sha256": {
            key: _sha(path.read_text(encoding="utf-8"))
            for key, path in paths.items()
            if key != "manifest"
        },
        "notes": (
            "The 118 cells are a queue ceiling, not observed coverage gain. Formal "
            "all-evidence rebuilding and species/lineage leave-one-out validation "
            "determine realized direct and Validated Low changes. Four Peperomia "
            "bagging results were retained in the research audit only because the "
            "source could not separate autogamy, wind pollination, and agamospermy."
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
    print(json.dumps(build_checkpoint(args.output_dir), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
