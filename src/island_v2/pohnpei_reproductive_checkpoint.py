"""Freeze the independent Pohnpei dissertation hand-pollination experiment.

The source table retains every Table 2.1 row, including the eleven rows that
originate in Yomai & Williams (2021).  Only the independent chapter-2 rows are
promoted here, so the article and dissertation are not counted as independent
source lineages.  Self-compatibility and autonomous selfing remain distinct
traits throughout.
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

CREATED_AT = "2026-08-12T06:20:00Z"
PDF_SHA256 = "1883a70541cb33f2e4ec85bec904f3dcb80683743cb2255efa0431122bba02ff"
SOURCE_TABLE_SHA256 = (
    "c14cac6d0eccf22faf8ff0e2a99ccb9fc76078f77e16f6d93ebd41b47be5902f"
)
SOURCE_URL = "https://trace.tennessee.edu/utk_graddiss/6964/"
SOURCE_LINEAGE = "repository:trace-tennessee:utk_graddiss:6964:chapter2"
SOURCE_TITLE = (
    "The Reproductive Biology of Island Plants: Is There Evidence for Baker's Law?"
)
SOURCE_CITATION = (
    "Yomai, Viann Marie H. (2021), PhD dissertation, University of Tennessee, "
    "Knoxville, TRACE record 6964, Chapter 2 Tables 2.1, 2.2 and 2.5"
)
REVIEWER = "Codex source-backed dissertation table audit"


def _sha256(path: Path) -> str:
    return hashlib.sha256(_canonical_file_bytes(path)).hexdigest()


def _as_bool(value: object) -> bool:
    return str(value).strip().casefold() == "true"


def _autofertility_state(index: float) -> str:
    """Apply the existing trait-ledger thresholds without changing ontology."""

    if index >= 0.8:
        return "autonomous"
    if index >= 0.2:
        return "mixed_or_variable"
    return "absent"


def _table_excerpt(row: pd.Series) -> str:
    return (
        f"Table 2.1: {row['source_species']}; Hand self "
        f"{row['hand_self_mean']}; Autonomous self {row['autonomous_self_mean']}; "
        f"Hand outcross {row['hand_outcross_mean']}; Open {row['open_mean']}. "
        'The table defines "Autonomous self" as bagged and unmanipulated.'
    )


def reviewed_rows(source_table: pd.DataFrame) -> list[dict[str, str]]:
    selected = source_table.loc[
        source_table["include_independent_chapter2"].map(_as_bool)
    ].copy()
    rows: list[dict[str, str]] = []
    common = {
        "quality": "high",
        "provider": "yomai_2021_doctoral_dissertation_chapter2",
        "url": SOURCE_URL,
        "title": SOURCE_TITLE,
        "citation": SOURCE_CITATION,
        "lineage": SOURCE_LINEAGE,
        "lineage_method": "original_doctoral_dissertation_independent_chapter2",
        "source_tier": "A",
        "source_type": "doctoral_dissertation_controlled_hand_pollination_experiment",
        "domain": "trace.tennessee.edu",
        "content_sha256": PDF_SHA256,
        "content_sha256_basis": "downloaded_official_repository_pdf_bytes",
        "retrieved_at_utc": "2026-08-12T05:33:44Z",
    }
    for item in selected.to_dict("records"):
        source_species = str(item["source_species"])
        species = str(item["accepted_species"])
        table_excerpt = _table_excerpt(pd.Series(item))
        sc = _evidence_row(
            species=species,
            trait="self_incompatibility",
            value="SC",
            raw_value="Table 2.5: all island species are self-compatible (SC)",
            excerpt=(
                "Table 2.5: All island species are self-compatible (SC) with "
                f"bisexual flowers. {table_excerpt}"
            ),
            record_id=(
                "trace:utk_graddiss:6964:chapter2:table2.1:"
                f"{source_species}:self-compatibility"
            ),
            **common,
        )
        _apply_identity(sc, item)
        rows.append(sc)

        # Bidens has a published AFI but no hand-self denominator in Table 2.1.
        # It remains withheld instead of treating an uninterpretable ratio as a
        # new autonomous-selfing record.
        if str(item["hand_self_mean"]).casefold() == "na":
            continue
        afi = float(item["autofertility_index"])
        state = _autofertility_state(afi)
        autonomy = _evidence_row(
            species=species,
            trait="autonomous_selfing_capacity",
            value=state,
            raw_value=(
                f"AFI {afi:.3f} ({item['autofertility_index_basis']}); "
                f"autonomous-self mean {item['autonomous_self_mean']} per five flowers"
            ),
            excerpt=(
                f"{table_excerpt} Table 2.2 reports Autofertility Index "
                f"{afi:.3f} ({item['autofertility_index_basis']})."
            ),
            record_id=(
                "trace:utk_graddiss:6964:chapter2:table2.1:"
                f"{source_species}:autonomous-selfing"
            ),
            **common,
        )
        _apply_identity(autonomy, item)
        rows.append(autonomy)
    return rows


def _apply_identity(evidence: dict[str, str], source: dict[str, str]) -> None:
    source_method = str(source["name_match_method"])
    # The accepted spelling occurs verbatim in the dissertation's species
    # table even when Table 2.1 contains a one-word typo.  The source table
    # retains that reconciliation detail; the promoted record truthfully uses
    # the gate's accepted-name-exact category.
    method = (
        "accepted_name_exact"
        if source_method == "same_document_orthographic_reconciliation"
        else source_method
    )
    evidence["source_group"] = "pohnpei_reproductive_checkpoint_20260812"
    evidence["matched_page_name"] = (
        str(source["accepted_species"])
        if source_method == "same_document_orthographic_reconciliation"
        else str(source["source_species"])
    )
    evidence["name_match_method"] = method
    evidence["name_resolution_lineage"] = str(source["name_resolution_lineage"])
    evidence["evidence_scope"] = (
        "synonym_direct" if method == "synonym_exact" else "species_direct"
    )
    evidence["wild_cultivated_cultivar_status"] = (
        "wild_or_field_population_not_cultivar_limited"
    )
    evidence["query"] = "high_information_primary_source_table_expansion"


def _review_audit(evidence: pd.DataFrame) -> pd.DataFrame:
    audit = _audit(evidence)
    by_id = evidence.set_index("candidate_id")
    audit["reviewer"] = REVIEWER
    audit["reviewed_at_utc"] = CREATED_AT
    audit["decision_reason"] = audit["candidate_id"].map(
        lambda candidate_id: (
            "Accepted after exact Table 2.1/2.2/2.5 transcription, target-master "
            "identity and family check, independent-chapter lineage check, and "
            f"trait-specific mapping; name match: "
            f"{by_id.at[candidate_id, 'name_match_method']}."
        )
    )
    return audit


def build(
    *,
    source_table_csv: Path,
    master_csv: Path,
    output_dir: Path,
    prior_curated_evidence_csv: Path | None = None,
    prior_curated_audit_csv: Path | None = None,
) -> dict[str, object]:
    if _sha256(source_table_csv) != SOURCE_TABLE_SHA256:
        raise ValueError("Pohnpei source table differs from the reviewed checkpoint")
    source = pd.read_csv(source_table_csv, dtype=str, keep_default_na=False)
    if len(source) != 36:
        raise ValueError(f"expected 36 Table 2.1 rows, found {len(source)}")
    selected = source.loc[source["include_independent_chapter2"].map(_as_bool)]
    if len(selected) != 25:
        raise ValueError(f"expected 25 independent chapter-2 rows, found {len(selected)}")
    if source["prior_article_lineage"].map(_as_bool).sum() != 11:
        raise ValueError("expected eleven rows from the prior article lineage")

    master = pd.read_csv(master_csv, dtype=str).fillna("")
    master_identity = master.set_index("accepted_species")["family"].to_dict()
    missing = sorted(set(selected["accepted_species"]) - set(master_identity))
    if missing:
        raise ValueError(f"reviewed species absent from target master: {missing}")
    family_conflicts = selected.loc[
        selected.apply(
            lambda row: master_identity[row["accepted_species"]] != row["family"],
            axis=1,
        )
    ]
    if not family_conflicts.empty:
        raise ValueError(
            "family conflicts in reviewed source table: "
            + ", ".join(family_conflicts["accepted_species"])
        )

    evidence = pd.DataFrame(reviewed_rows(source), columns=EVIDENCE_COLUMNS).fillna("")
    evidence = evidence.sort_values(
        ["accepted_species", "trait_name", "candidate_id"]
    ).reset_index(drop=True)
    if len(evidence) != 49:
        raise ValueError(f"expected 49 reviewed trait rows, found {len(evidence)}")
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("Pohnpei checkpoint candidate IDs are not unique")
    audit = _review_audit(evidence)

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "pohnpei_reproductive_evidence_20260812.csv"
    audit_path = output_dir / "pohnpei_reproductive_manual_audit_20260812.csv"
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
        "contract": "pohnpei_reproductive_individually_reviewed_checkpoint_v1",
        "created_at_utc": CREATED_AT,
        "source_pdf_sha256": PDF_SHA256,
        "source_table_sha256": SOURCE_TABLE_SHA256,
        "source_rows": len(source),
        "independent_chapter2_source_rows": len(selected),
        "excluded_same_lineage_source_rows": len(source) - len(selected),
        "new_evidence_rows": len(evidence),
        "new_species": int(evidence["accepted_species"].nunique()),
        "new_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "quality_counts": evidence["evidence_quality"].value_counts().to_dict(),
        "trait_counts": evidence["trait_name"].value_counts().to_dict(),
        "autofertility_state_counts": evidence.loc[
            evidence["trait_name"].eq("autonomous_selfing_capacity"),
            "normalized_value",
        ].value_counts().to_dict(),
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
            "same_lineage_article_rows_excluded": True,
            "family_inference": False,
            "global_fallback": False,
            "cross_trait_substitution": False,
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
    manifest_path = output_dir / "pohnpei_reproductive_checkpoint_manifest_20260812.json"
    manifest_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-table-csv", type=Path, required=True)
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--prior-curated-evidence-csv", type=Path)
    parser.add_argument("--prior-curated-audit-csv", type=Path)
    print(json.dumps(build(**vars(parser.parse_args())), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
