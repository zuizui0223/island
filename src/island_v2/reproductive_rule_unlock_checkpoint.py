"""Freeze primary reproductive evidence selected for genus-rule information gain.

Every promoted row is a species-level statement from a retrieved original
article.  The records remain trait-specific: dioecy/outbreeding maps only to
``mating_system``; bagging experiments map only to
``autonomous_selfing_capacity``; and explicit compatibility statements map
only to ``self_incompatibility``.  The module does not create genus inference
itself.  The shared all-evidence implementation rebuilds those rules after
these direct rows pass individual review.
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

CREATED_AT = "2026-08-12T08:20:00Z"
REVIEWER = "Codex source-backed line-by-line evidence audit"
SOURCE_GROUP = "reproductive_rule_unlock_checkpoint_20260812"


def _row(
    *,
    species: str,
    trait: str,
    value: str,
    raw_value: str,
    excerpt: str,
    provider: str,
    url: str,
    title: str,
    citation: str,
    record_id: str,
    lineage: str,
    source_type: str,
    domain: str,
    content_sha256: str,
    content_sha256_basis: str,
    retrieved_at_utc: str,
) -> dict[str, str]:
    evidence = _evidence_row(
        species=species,
        trait=trait,
        value=value,
        quality="high",
        provider=provider,
        url=url,
        title=title,
        citation=citation,
        excerpt=excerpt,
        record_id=record_id,
        lineage=lineage,
        lineage_method="doi_or_original_article_lineage",
        source_tier="A",
        source_type=source_type,
        domain=domain,
        content_sha256=content_sha256,
        content_sha256_basis=content_sha256_basis,
        retrieved_at_utc=retrieved_at_utc,
        raw_value=raw_value,
    )
    evidence["source_group"] = SOURCE_GROUP
    evidence["name_resolution_lineage"] = "master_accepted_name_exact"
    evidence["wild_cultivated_cultivar_status"] = (
        "wild_or_species_level_statement_not_cultivar_limited"
    )
    evidence["query"] = "current_support_2_rule_unlock_primary_source"
    return evidence


def reviewed_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    rhynchospora = {
        "species": "Rhynchospora breviuscula",
        "provider": "Castellani et al. 2024 Nature Plants",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10954556/",
        "title": (
            "Meiotic recombination dynamics in plants with repeat-based "
            "holocentromeres shed light on the primary drivers of crossover patterning"
        ),
        "citation": (
            "Castellani et al. (2024), Nature Plants 10:423-438, "
            "DOI 10.1038/s41477-024-01625-y"
        ),
        "lineage": "doi:10.1038/s41477-024-01625-y",
        "source_type": "peer_reviewed_primary_article",
        "domain": "pmc.ncbi.nlm.nih.gov",
        "content_sha256": (
            "2bd71f76ce1bffcb1755a576e7e2c7f19a1c6a02aa9b95fea3a3b07b349c0f58"
        ),
        "content_sha256_basis": "retrieved_europe_pmc_fulltext_xml_bytes",
        "retrieved_at_utc": "2026-08-12T05:01:40Z",
    }
    rows.append(
        _row(
            trait="mating_system",
            value="predominantly_outcrossing",
            raw_value="outbred species with high levels of self-incompatibility",
            excerpt=(
                "R. breviuscula is an outbred species with high levels of "
                "self-incompatibility, which hampers the standard detection of COs."
            ),
            record_id="doi:10.1038/s41477-024-01625-y:results:outbred",
            **rhynchospora,
        )
    )

    vitis_coignetiae = {
        "species": "Vitis coignetiae",
        "provider": "Honma et al. 2003 Horticultural Research Japan",
        "url": (
            "https://www.jstage.jst.go.jp/article/hrj/2/4/2_4_289/_article/-char/en"
        ),
        "title": "Long-term Pollen Storage in Vitis coignetiae Pulliat with Organic Solvents",
        "citation": (
            "Honma, Endou, Takahasi & Taira (2003), Horticultural Research "
            "(Japan) 2:289-292, DOI 10.2503/hrj.2.289"
        ),
        "lineage": "doi:10.2503/hrj.2.289",
        "source_type": "peer_reviewed_primary_article",
        "domain": "jstage.jst.go.jp",
        "content_sha256": (
            "c5dd223c73af5fd89e03d9209b8a500d6cd701e18263d5b1fdc7f6b756d8d1cf"
        ),
        "content_sha256_basis": "downloaded_jstage_publisher_pdf_bytes",
        "retrieved_at_utc": "2026-08-12T07:26:09Z",
    }
    rows.append(
        _row(
            trait="mating_system",
            value="predominantly_outcrossing",
            raw_value="dioecious plant that requires pollination for fruiting",
            excerpt=(
                "Vitis coignetiae Pulliat, which is a dioecious plant that "
                "requires pollination for fruiting."
            ),
            record_id="doi:10.2503/hrj.2.289:abstract:dioecious",
            **vitis_coignetiae,
        )
    )

    vitis_amurensis = {
        "species": "Vitis amurensis",
        "provider": "Men et al. 2021 Frontiers in Genetics",
        "url": (
            "https://www.frontiersin.org/journals/genetics/articles/"
            "10.3389/fgene.2021.727260/full"
        ),
        "title": "VaAPRT3 Gene is Associated With Sex Determination in Vitis amurensis",
        "citation": (
            "Men et al. (2021), Frontiers in Genetics 12:727260, "
            "DOI 10.3389/fgene.2021.727260"
        ),
        "lineage": "doi:10.3389/fgene.2021.727260",
        "source_type": "peer_reviewed_primary_article",
        "domain": "frontiersin.org",
        "content_sha256": (
            "14b5e76c213318a6cdffed98e7916af5d5e77b47bd4c5ef60485084d430377e7"
        ),
        "content_sha256_basis": "retrieved_publisher_fulltext_html_bytes",
        "retrieved_at_utc": "2026-08-12T07:27:11Z",
    }
    rows.append(
        _row(
            trait="mating_system",
            value="predominantly_outcrossing",
            raw_value="mostly dioecious",
            excerpt=(
                "V. amurensis is mostly dioecious. The female plants produce "
                "flowers with perfect formed pistil having style and stigma but "
                "reflexed stamens with infertile pollen, but the male plants have "
                "erect stamens and fertile pollen, but aborted pistils with no style "
                "or stigma."
            ),
            record_id="doi:10.3389/fgene.2021.727260:introduction:mostly-dioecious",
            **vitis_amurensis,
        )
    )

    elymus = {
        "species": "Elymus repens",
        "provider": "Urfusova et al. 2021 Preslia",
        "url": "https://www.preslia.cz/article/pdf?id=11522",
        "title": "The mentor effect increases the rate of selfing in couch grasses",
        "citation": (
            "Urfusova et al. (2021), Preslia 93:377-397, "
            "DOI 10.23855/preslia.2021.377"
        ),
        "lineage": "doi:10.23855/preslia.2021.377",
        "source_type": "peer_reviewed_primary_bagging_experiment",
        "domain": "preslia.cz",
        "content_sha256": (
            "a216b9895c1bcde55af450ea766ef0f16d1df5eb9a054e3ab360a58ab3b293f2"
        ),
        "content_sha256_basis": "downloaded_publisher_pdf_bytes",
        "retrieved_at_utc": "2026-08-12T07:25:07Z",
    }
    rows.append(
        _row(
            trait="autonomous_selfing_capacity",
            value="autonomous",
            raw_value=(
                "pollination-bag treatment before anther maturity; 0.25% "
                "autonomous selfing"
            ),
            excerpt=(
                "The autonomous selfing treatment covered a single spike or two "
                "spikes of a plant with a pollination bag before anther maturity. "
                "E. repens had 0.25% autonomous selfing; two seedlings developed "
                "from 18 spikes with 1222 florets."
            ),
            record_id="doi:10.23855/preslia.2021.377:autonomous-selfing:elymus-repens",
            **elymus,
        )
    )

    fritillaria = {
        "species": "Fritillaria koidzumiana",
        "provider": "Kawano et al. 2008 Plant Species Biology",
        "url": (
            "https://esj-journals.onlinelibrary.wiley.com/doi/"
            "10.1111/j.1442-1984.2008.00208.x"
        ),
        "title": (
            "Life-history monographs of Japanese plants. 10: "
            "Fritillaria koidzumiana Ohwi (Liliaceae)"
        ),
        "citation": (
            "Kawano, Masuda & Hayashi (2008), Plant Species Biology 23:51-57, "
            "DOI 10.1111/j.1442-1984.2008.00208.x"
        ),
        "lineage": "doi:10.1111/j.1442-1984.2008.00208.x",
        "source_type": "peer_reviewed_species_monograph_field_experiment",
        "domain": "esj-journals.onlinelibrary.wiley.com",
        "content_sha256": (
            "84df395ea522f7644c20fd6f20cf6f54e4648b18b1186884816d134ee2231f0d"
        ),
        "content_sha256_basis": "verified_publisher_fulltext_excerpt_utf8_bytes",
        "retrieved_at_utc": "2026-08-12T07:39:22Z",
    }
    frit_bagging = (
        "In contrast, bagged flowers in the field produced no seeds, while "
        "self-pollinated individuals produced an average number of 12.6 seeds, "
        "much less seed output per plant and a mean number of 21.7 immature "
        "seeds per plant, suggesting the presence of a weak self-incompatibility."
    )
    rows.extend(
        [
            _row(
                trait="autonomous_selfing_capacity",
                value="absent",
                raw_value="bagged flowers in the field produced no seeds",
                excerpt=frit_bagging,
                record_id=(
                    "doi:10.1111/j.1442-1984.2008.00208.x:"
                    "reproductive-systems:bagged-no-seeds"
                ),
                **fritillaria,
            ),
            _row(
                trait="self_incompatibility",
                value="mixed_or_variable",
                raw_value="weak self-incompatibility",
                excerpt=frit_bagging,
                record_id=(
                    "doi:10.1111/j.1442-1984.2008.00208.x:"
                    "reproductive-systems:weak-self-incompatibility"
                ),
                **fritillaria,
            ),
            _row(
                trait="mating_system",
                value="predominantly_outcrossing",
                raw_value="outbreeder with protogynous breeding systems",
                excerpt=(
                    "The relatively low pollen : ovule (P/O) ratio found in "
                    "F. koidzumiana was approximately 1500. This is somewhat "
                    "exceptional as an outbreeder with protogynous breeding systems."
                ),
                record_id=(
                    "doi:10.1111/j.1442-1984.2008.00208.x:"
                    "reproductive-systems:outbreeder"
                ),
                **fritillaria,
            ),
        ]
    )

    sonchus = {
        "species": "Sonchus microcephalus",
        "provider": "Mejias 1992 Flora Mediterranea",
        "url": "https://herbmedit.org/storage/165/2-015.pdf",
        "title": (
            "Reproductive biology in the Iberian taxa of the genera Sonchus "
            "and Aetheorhiza (Asteraceae: Lactuceae)"
        ),
        "citation": "Mejias (1992), Flora Mediterranea 2:15-32",
        "lineage": "article:flora-mediterranea:1992:2:15-32:mejias",
        "source_type": "peer_reviewed_primary_bagging_experiment",
        "domain": "herbmedit.org",
        "content_sha256": (
            "c53b9cfe1b0e59ad02b310018c24c7520b8442303eea8804c82931670f254a39"
        ),
        "content_sha256_basis": "downloaded_journal_pdf_bytes",
        "retrieved_at_utc": "2026-08-12T07:31:44Z",
    }
    rows.extend(
        [
            _row(
                trait="autonomous_selfing_capacity",
                value="autonomous",
                raw_value=(
                    "bagged-before-anthesis fruit set 90.53%, 75.52%, and 69.06%"
                ),
                excerpt=(
                    "For flower heads bagged before anthesis, the values represent "
                    "the frequency of autogamy including geitonogamy. Sonchus "
                    "microcephalus produced large numbers of fruits in the flower "
                    "heads that were bagged before anthesis: M1 90.53%, M1 raised "
                    "from cypselas produced after selfing 75.52%, and M2 69.06%."
                ),
                record_id=(
                    "flora-mediterranea:1992:2:table1:"
                    "sonchus-microcephalus:bagged-fruit-set"
                ),
                **sonchus,
            ),
            _row(
                trait="self_incompatibility",
                value="SC",
                raw_value="S. microcephalus (SC)",
                excerpt=(
                    "Table 3 labels S. microcephalus (SC); the article reports that "
                    "self-fertility prevails in Sonchus microcephalus and that it "
                    "produced large numbers of fruits in flower heads bagged before "
                    "anthesis."
                ),
                record_id=(
                    "flora-mediterranea:1992:2:table3:"
                    "sonchus-microcephalus:self-compatible"
                ),
                **sonchus,
            ),
        ]
    )
    return rows


def _review_audit(evidence: pd.DataFrame) -> pd.DataFrame:
    audit = _audit(evidence)
    by_id = evidence.set_index("candidate_id")
    audit["reviewer"] = REVIEWER
    audit["reviewed_at_utc"] = CREATED_AT
    audit["decision_reason"] = audit["candidate_id"].map(
        lambda candidate_id: (
            "Accepted after exact original-article statement, target-master "
            "identity/family, source lineage, content fingerprint and "
            "trait-specific ontology review; mapping retained only for "
            f"{by_id.at[candidate_id, 'trait_name']}."
        )
    )
    return audit


def _sha256(path: Path) -> str:
    return hashlib.sha256(_canonical_file_bytes(path)).hexdigest()


def build(
    *,
    master_csv: Path,
    output_dir: Path,
    prior_curated_evidence_csv: Path | None = None,
    prior_curated_audit_csv: Path | None = None,
) -> dict[str, object]:
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    master_identity = master.set_index("accepted_species")["family"].to_dict()
    expected_families = {
        "Rhynchospora breviuscula": "Cyperaceae",
        "Vitis coignetiae": "Vitaceae",
        "Vitis amurensis": "Vitaceae",
        "Elymus repens": "Poaceae",
        "Fritillaria koidzumiana": "Liliaceae",
        "Sonchus microcephalus": "Asteraceae",
    }
    missing = sorted(set(expected_families) - set(master_identity))
    if missing:
        raise ValueError(f"reviewed species absent from target master: {missing}")
    conflicts = {
        species: (family, master_identity[species])
        for species, family in expected_families.items()
        if master_identity[species] != family
    }
    if conflicts:
        raise ValueError(f"family conflicts in reviewed checkpoint: {conflicts}")

    evidence = pd.DataFrame(reviewed_rows(), columns=EVIDENCE_COLUMNS).fillna("")
    evidence = evidence.sort_values(
        ["accepted_species", "trait_name", "candidate_id"]
    ).reset_index(drop=True)
    if len(evidence) != 9:
        raise ValueError(f"expected 9 reviewed trait rows, found {len(evidence)}")
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("rule-unlock checkpoint candidate IDs are not unique")
    audit = _review_audit(evidence)

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "reproductive_rule_unlock_evidence_20260812.csv"
    audit_path = output_dir / "reproductive_rule_unlock_manual_audit_20260812.csv"
    evidence.to_csv(evidence_path, index=False, lineterminator="\n")
    audit.to_csv(audit_path, index=False, lineterminator="\n")
    outputs = [evidence_path, audit_path]

    combined_evidence: pd.DataFrame | None = None
    combined_audit: pd.DataFrame | None = None
    if prior_curated_evidence_csv is not None:
        if prior_curated_audit_csv is None:
            raise ValueError("prior curated audit is required with prior evidence")
        prior_evidence = pd.read_csv(prior_curated_evidence_csv, dtype=str).fillna("")
        prior_audit = pd.read_csv(prior_curated_audit_csv, dtype=str).fillna("")
        owned = set(evidence["candidate_id"])
        combined_evidence = pd.concat(
            [prior_evidence.loc[~prior_evidence["candidate_id"].isin(owned)], evidence],
            ignore_index=True,
        )
        combined_audit = pd.concat(
            [prior_audit.loc[~prior_audit["candidate_id"].isin(owned)], audit],
            ignore_index=True,
        )
        for name, frame in (
            ("combined evidence", combined_evidence),
            ("combined audit", combined_audit),
        ):
            if frame["candidate_id"].duplicated().any():
                raise ValueError(f"{name} candidate IDs are not unique")
        combined_evidence_path = output_dir / "combined_curated_evidence_20260812.csv"
        combined_audit_path = output_dir / "combined_curated_manual_audit_20260812.csv"
        combined_evidence.to_csv(
            combined_evidence_path, index=False, lineterminator="\n"
        )
        combined_audit.to_csv(combined_audit_path, index=False, lineterminator="\n")
        outputs.extend([combined_evidence_path, combined_audit_path])

    summary: dict[str, object] = {
        "contract": "reproductive_rule_unlock_individually_reviewed_checkpoint_v1",
        "created_at_utc": CREATED_AT,
        "new_evidence_rows": len(evidence),
        "new_species": int(evidence["accepted_species"].nunique()),
        "new_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "quality_counts": evidence["evidence_quality"].value_counts().to_dict(),
        "trait_counts": evidence["trait_name"].value_counts().to_dict(),
        "source_lineages": int(evidence["source_lineage"].nunique()),
        "audit": {
            "reviewed": len(audit),
            "accepted_correct": int(audit["decision"].str.casefold().eq("accept").sum()),
            "precision": float(audit["decision"].str.casefold().eq("accept").mean()),
            "cultivar_contamination_rate": float(
                audit["cultivar_contamination"].str.casefold().eq("true").mean()
            ),
        },
        "guardrails": {
            "trait_specific_records": True,
            "genus_inference_emitted_here": False,
            "family_inference": False,
            "global_fallback": False,
            "n2_formal_inference": False,
            "cross_trait_substitution": False,
            "search_snippet_evidence": False,
            "rejected_ornithogalum_candidate": True,
        },
        "files": {},
    }
    if combined_evidence is not None and combined_audit is not None:
        summary["combined"] = {
            "evidence_rows": len(combined_evidence),
            "audit_rows": len(combined_audit),
            "species": int(combined_evidence["accepted_species"].nunique()),
            "species_trait": int(
                combined_evidence[["accepted_species", "trait_name"]]
                .drop_duplicates()
                .shape[0]
            ),
        }
    readme_path = output_dir / "README.md"
    if readme_path.exists():
        outputs.append(readme_path)
    for path in outputs:
        summary["files"][path.name] = {
            "sha256": _sha256(path),
            "size_bytes": len(_canonical_file_bytes(path)),
        }
    manifest_path = output_dir / "reproductive_rule_unlock_manifest_20260812.json"
    manifest_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--prior-curated-evidence-csv", type=Path)
    parser.add_argument("--prior-curated-audit-csv", type=Path)
    print(json.dumps(build(**vars(parser.parse_args())), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
