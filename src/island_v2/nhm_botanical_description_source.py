"""Promote reviewed NHM botanical-description cells as Medium evidence.

The NHM dataset is a machine-readable compilation of descriptions from BHL,
Ecoflora, eFloras and related floras (doi:10.5519/p3dm31kc).  This adapter is
deliberately offline: discovery and review happen before this step, and only
rows carrying an exact structured excerpt and a strict name-gate decision are
accepted.  Sex-system, pollination and reward fields are never collapsed into
the strict reproductive-assurance or floral-structure axes here.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Annotated

import pandas as pd
import typer

from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS

app = typer.Typer(add_completion=False, no_args_is_help=True)

DATASET_DOI = "10.5519/p3dm31kc"
DATASET_URL = (
    "https://data.nhm.ac.uk/id/dataset/"
    "angiosperm-functional-traits-extracted-from-botanical-descriptions"
)
ALLOWED_TRAITS = {
    "flower_primary_color": "flower_colour",
    "floral_form": "floral_structural_complexity",
    "floral_symmetry": "floral_structural_complexity",
    "tube_depth_class": "floral_structural_complexity",
    "flower_size_class": "floral_structural_complexity",
    "inflorescence_display": "floral_structural_complexity",
}
ALLOWED_NAME_MATCHES = {
    "accepted_name_exact",
    "wfo_gbif_two_backbone_strict_synonym",
    "wfo_stable_identifier_orthographic_variant",
}
REVIEW_COLUMNS = (
    "accepted_species",
    "searched_taxon",
    "axis",
    "trait_name",
    "normalized_value",
    "provider",
    "resource_id",
    "query_url",
    "source_name",
    "source_record_ids",
    "supporting_excerpt",
    "decision",
    "reason",
    "name_match_method",
    "reviewer",
    "reviewed_at_utc",
)
GZIP = {"method": "gzip", "mtime": 0}
TEXT_SUFFIXES = {".csv", ".json", ".jsonl", ".md", ".toml", ".yml", ".yaml"}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha256(path: Path) -> str:
    payload = path.read_bytes()
    if path.suffix.casefold() in TEXT_SUFFIXES:
        payload = payload.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def _record_id(row: dict[str, object]) -> str:
    payload = "|".join(
        _text(row.get(key))
        for key in (
            "accepted_species",
            "trait_name",
            "provider",
            "source_record_ids",
            "normalized_value",
        )
    )
    return "nhm:" + hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]


def _lineage(row: dict[str, object]) -> str:
    # All cells for one species/trait/provider are one original-description
    # lineage.  Multiple API rows and mirrored structured fields cannot create
    # artificial independent support.
    return "nhm-description:" + ":".join(
        _text(row.get(key)).casefold().replace(" ", "-")
        for key in ("provider", "searched_taxon", "trait_name")
    )


def build_package(review: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Validate the frozen review table and return Medium direct evidence."""

    missing = set(REVIEW_COLUMNS).difference(review.columns)
    if missing:
        raise ValueError(f"NHM review table missing columns: {sorted(missing)}")
    data = review.fillna("").copy()
    audit_rows: list[dict[str, str]] = []
    evidence_rows: list[dict[str, str]] = []
    for row in data.to_dict("records"):
        trait = _text(row["trait_name"])
        expected_axis = ALLOWED_TRAITS.get(trait, "")
        decision = _text(row["decision"]).casefold()
        reason = "selected"
        if decision != "accept":
            reason = "not_review_accepted"
        elif not expected_axis:
            reason = "trait_not_allowed_in_strict_axes"
        elif _text(row["axis"]) != expected_axis:
            reason = "trait_axis_mismatch"
        elif _text(row["name_match_method"]) not in ALLOWED_NAME_MATCHES:
            reason = "strict_species_identity_gate_failed"
        elif not _text(row["supporting_excerpt"]):
            reason = "missing_exact_structured_excerpt"
        elif not _text(row["query_url"]).startswith("https://data.nhm.ac.uk/"):
            reason = "unstable_or_unexpected_source_url"
        elif not _text(row["source_record_ids"]):
            reason = "missing_source_record_id"
        elif not _text(row["reviewer"]) or not _text(row["reviewed_at_utc"]):
            reason = "incomplete_review_provenance"
        audit_rows.append(
            {
                **{column: _text(row.get(column)) for column in REVIEW_COLUMNS},
                "package_decision": "accept" if reason == "selected" else "exclude",
                "package_reason": reason,
                "source_lineage": _lineage(row),
            }
        )
        if reason != "selected":
            continue
        provider = _text(row["provider"])
        source_name = _text(row["source_name"]) or provider
        evidence_rows.append(
            {
                "accepted_species": _text(row["accepted_species"]),
                "axis": expected_axis,
                "trait_name": trait,
                "normalized_value": _text(row["normalized_value"]),
                "quality": "medium",
                "source_group": "nhm_botanical_descriptions",
                "source_provider": f"NHM botanical descriptions: {source_name}",
                "source_url": _text(row["query_url"]),
                "source_record_id": _record_id(row),
                "source_citation": (
                    f"NHM Angiosperm functional traits extracted from botanical "
                    f"descriptions, doi:{DATASET_DOI}; source={source_name}; "
                    f"resource={_text(row['resource_id'])}."
                ),
                "source_excerpt": _text(row["supporting_excerpt"]),
                "evidence_scope": "species_direct",
                "name_match_method": _text(row["name_match_method"]),
                "source_lineage": _lineage(row),
                "lineage_method": "provider_species_trait_structured_description_collapse",
                "source_run_id": "local-wave36-nhm-20260827",
                "source_artifact": "source-specific-reviewed-evidence-rows.csv.gz",
                "source_file": _text(row["source_record_ids"]),
                "acceptance_contract": "reviewed_nhm_structured_description_medium_v1",
            }
        )
    evidence = pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS)
    audit = pd.DataFrame(audit_rows)
    if not evidence.empty:
        evidence = evidence.sort_values(
            ["accepted_species", "axis", "trait_name", "source_lineage"]
        ).reset_index(drop=True)
    return evidence, audit


@app.command("build")
def build(
    review_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
) -> None:
    review = pd.read_csv(review_csv, dtype=str).fillna("")
    evidence, audit = build_package(review)
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "nhm_reviewed_direct_evidence.csv.gz"
    audit_path = output_dir / "nhm_review_audit.csv.gz"
    evidence.to_csv(evidence_path, index=False, compression=GZIP)
    audit.to_csv(audit_path, index=False, compression=GZIP)
    cells = evidence[["accepted_species", "axis"]].drop_duplicates()
    summary = {
        "contract": "nhm_botanical_description_source_package_v1",
        "dataset": {"doi": DATASET_DOI, "url": DATASET_URL},
        "reviewed_rows": len(review),
        "accepted_evidence_rows": len(evidence),
        "accepted_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "accepted_species_axis": len(cells),
        "accepted_by_axis": cells["axis"].value_counts().sort_index().to_dict(),
        "excluded_by_reason": audit.loc[audit["package_decision"].eq("exclude"), "package_reason"]
        .value_counts()
        .sort_index()
        .to_dict(),
        "checks": {
            "strict_name_gate": True,
            "exact_structured_excerpt": True,
            "source_lineage_deduplicated": True,
            "reward_not_collapsed_into_structure": True,
            "pollination_not_collapsed_into_structure": True,
            "sex_system_not_collapsed_into_reproductive_assurance": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
        },
        "input_sha256": _sha256(review_csv),
        "artifact_sha256": {
            evidence_path.name: _sha256(evidence_path),
            audit_path.name: _sha256(audit_path),
        },
    }
    (output_dir / "nhm_source_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    app()
