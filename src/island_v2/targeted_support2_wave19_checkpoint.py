"""Build the Wave 19 reviewed support-two evidence checkpoint.

Wave 19 adds exact-species evidence for three current support=2 genus x trait
rules plus two additional traits stated by the same Pulicaria species
treatment.  It never creates genus inference itself; the shared all-evidence
rebuild remains the sole producer of Validated Low rows.
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

SOURCE_GROUP = "targeted_support2_wave19_checkpoint_20260820"
CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/targeted_support2_wave19_checkpoint_20260820"
)
PRIOR = Path(
    "data/v2/staging/traits/open_web_pilot/targeted_support2_wave18_checkpoint_20260820"
)
RETRIEVED_AT = "2026-08-20T18:00:00Z"


def _row(*args: str, **kwargs: str) -> dict[str, str]:
    row = _wave15_row(*args, **kwargs)
    row["source_group"] = SOURCE_GROUP
    row["retrieved_at_utc"] = RETRIEVED_AT
    return row


def primary_rows() -> list[dict[str, str]]:
    """Return five individually reviewed exact-species evidence rows."""
    pulicaria_oman_excerpt = (
        "Capitula radiate, at the end of the stems in a loosely corymbiform "
        "synflorescence, rarely single. Flowers yellow; ray flowers pistillate, "
        "ray lamina 1-2 mm long, disk flowers perfect, tubular."
    )
    pulicaria_iraq_excerpt = (
        "Female florets 2-4 mm long, disc florets 2-3 mm long, lobes glandular "
        "and papillose."
    )
    lilium_excerpt = (
        "It has been characterized as a butterfly pollinated obligate "
        "outcrosser (Edwards and Jordan, 1992)."
    )

    rows = [
        _row(
            "Lilium philadelphicum",
            "reproductive_assurance",
            "mating_system",
            "obligate outcrosser",
            "predominantly_outcrossing",
            "medium",
            "USDA Forest Service Treesearch / American Midland Naturalist",
            "https://research.fs.usda.gov/download/treesearch/34945.pdf",
            "Conservation Genetics of Remnant Lilium philadelphicum Populations",
            (
                "Horning & Webster 2009. American Midland Naturalist 161:286-300; "
                "characterization attributed to Edwards & Jordan 1992."
            ),
            lilium_excerpt,
            "cited_original:edwards-jordan-1992:lilium-philadelphicum",
            "A",
            "peer_reviewed_article_species_background_with_original_citation",
            "en",
            '"Lilium philadelphicum" obligate outcrosser',
            content_sha256=(
                "e53cddeb02d75a1b1af7c4e86e44e62390c32ad2b3ee45297a38c105a07a550a"
            ),
            content_sha256_basis="retrieved_usda_hosted_article_pdf_bytes",
        ),
    ]

    pulicaria_common = {
        "provider": "Kew Plants of the World Online / Flora of Oman",
        "url": (
            "https://powo.science.kew.org/taxon/"
            "urn%3Alsid%3Aipni.org%3Anames%3A240490-1/general-information"
        ),
        "title": "Pulicaria arabica in Plants of the World Online",
        "citation": (
            "Ghazanfar 2015. Flora of Oman 3: Loganiaceae-Asteraceae; "
            "Kew Plants of the World Online accepted-species treatment."
        ),
        "excerpt": pulicaria_oman_excerpt,
        "lineage": "flora-treatment:flora-of-oman:Pulicaria_arabica",
        "tier": "A",
        "source_type": "official_flora_exact_species_treatment",
        "language": "en",
    }
    for trait, raw, value, axis in (
        (
            "floral_symmetry",
            "radiate capitulum with ray flowers and tubular disk flowers",
            "actinomorphic|zygomorphic",
            "floral_structural_complexity",
        ),
        (
            "inflorescence_display",
            "capitula in a loosely corymbiform synflorescence",
            "composite_display|umbel_corymb",
            "floral_structural_complexity",
        ),
        (
            "flower_primary_color",
            "flowers yellow",
            "yellow_orange",
            "flower_colour",
        ),
    ):
        rows.append(
            _row(
                "Pulicaria arabica",
                axis,
                trait,
                raw,
                value,
                "high",
                pulicaria_common["provider"],
                pulicaria_common["url"],
                pulicaria_common["title"],
                pulicaria_common["citation"],
                pulicaria_common["excerpt"],
                pulicaria_common["lineage"],
                pulicaria_common["tier"],
                pulicaria_common["source_type"],
                pulicaria_common["language"],
                f'"Pulicaria arabica" {raw}',
            )
        )
    rows.append(
        _row(
            "Pulicaria arabica",
            "floral_structural_complexity",
            "flower_size_class",
            "female florets 2-4 mm; disc florets 2-3 mm",
            "small",
            "high",
            "Kew Plants of the World Online / Flora of Iraq",
            pulicaria_common["url"],
            pulicaria_common["title"],
            (
                "Ghazanfar, Edmondson & Hind 2019. Flora of Iraq 6: Compositae; "
                "Kew Plants of the World Online accepted-species treatment."
            ),
            pulicaria_iraq_excerpt,
            "flora-treatment:flora-of-iraq:Pulicaria_arabica",
            "A",
            "official_flora_exact_species_treatment",
            "en",
            '"Pulicaria arabica" floret length',
        )
    )
    if len(rows) != 5:
        raise AssertionError("Wave19 must define exactly five direct rows")
    return rows


def build_audit(evidence: pd.DataFrame) -> pd.DataFrame:
    audit = _wave15_build_audit(evidence)
    audit["reviewer"] = "Codex Wave19 individual source audit"
    audit["reviewed_at_utc"] = RETRIEVED_AT
    pulicaria_symmetry = (
        audit["accepted_species"].eq("Pulicaria arabica")
        & audit["trait_name"].eq("floral_symmetry")
    )
    audit.loc[pulicaria_symmetry, "decision_reason"] = (
        "Accepted after exact fixed-master identity review: the species-level "
        "official flora treatment explicitly describes a radiate capitulum with "
        "ray flowers and tubular disk flowers, preserving the multistate symmetry."
    )
    lilium = audit["accepted_species"].eq("Lilium philadelphicum")
    audit.loc[lilium, "decision_reason"] = (
        "Accepted as Medium after exact fixed-master identity and article quote "
        "review; the characterization is attributed to Edwards & Jordan 1992, so "
        "that cited original is preserved as lineage rather than counted as a new "
        "independent experiment."
    )
    return audit.loc[:, AUDIT_COLUMNS]


def build_checkpoint(output_dir: Path = CHECKPOINT) -> dict[str, Any]:
    evidence = pd.DataFrame(primary_rows(), columns=EVIDENCE_COLUMNS).sort_values(
        ["accepted_species", "trait_name", "source_lineage"], kind="stable"
    )
    evidence = evidence.reset_index(drop=True)
    audit = build_audit(evidence)
    if len(evidence) != 5 or len(audit) != 5:
        raise ValueError("Wave19 must contain exactly five individually reviewed rows")
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("Wave19 candidate IDs must be unique")
    if evidence.duplicated(["accepted_species", "trait_name"]).any():
        raise ValueError("Wave19 species x trait pairs must be unique")

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
        "evidence": output_dir / "targeted_support2_wave19_evidence_20260820.csv",
        "audit": output_dir / "targeted_support2_wave19_manual_audit_20260820.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260820.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260820.csv",
        "manifest": output_dir / "source_acquisition_manifest_wave19.json",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    manifest = {
        "checkpoint": SOURCE_GROUP,
        "built_at_utc": RETRIEVED_AT,
        "baseline_formal_run_id": 32344088168,
        "accepted_evidence_rows": len(evidence),
        "accepted_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "accepted_species": int(evidence["accepted_species"].nunique()),
        "recorded_queries": int(evidence["query"].nunique()),
        "formal_search_api_queries": 0,
        "search_cost_usd": 0.0,
        "targeted_support2_rules": [
            "Lilium|mating_system",
            "Pulicaria|floral_symmetry",
            "Pulicaria|flower_size_class",
        ],
        "theoretical_rule_cells_touched": 67,
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
            "cited_characterization_counted_as_independent_replication": False,
        },
        "output_sha256": {
            key: _sha(path.read_text(encoding="utf-8"))
            for key, path in paths.items()
            if key != "manifest"
        },
        "notes": (
            "The 67 cells are a queue ceiling, not observed coverage gain. Formal "
            "all-evidence rebuilding and species/lineage leave-one-out validation "
            "determine realized direct and Validated Low changes. Pulicaria colour "
            "and display rows are direct additions but are not counted in the rule "
            "ceiling."
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
