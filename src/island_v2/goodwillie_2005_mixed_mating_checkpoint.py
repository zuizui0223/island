"""Build reviewed reproductive evidence from the Goodwillie mixed-mating data.

The source workbook records multilocus outcrossing rates and an explicit
``Mode of selfing`` field.  This checkpoint deliberately does not reinterpret
outcrossing rates, ``none``, geitonogamy, or biparental inbreeding as another
trait.  It accepts only exact natural-population records that explicitly mark
autonomous selfing (``auton``) or cleistogamy (``cleis``).
"""

from __future__ import annotations

import argparse
import hashlib
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

SOURCE_GROUP = "goodwillie_2005_mixed_mating_checkpoint_20260821"
SOURCE_DIR = Path("data/v2/external/traits/goodwillie_2005")
SOURCE_CSV = SOURCE_DIR / "mixed_mating_outcrossing_rate_database_normalized.csv"
SOURCE_XLS = SOURCE_DIR / "mixed_mating_outcrossing_rate_database.xls"
SOURCE_README = SOURCE_DIR / "README.txt"
CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/goodwillie_2005_mixed_mating_checkpoint_20260821"
)
PRIOR = Path("data/v2/staging/traits/open_web_pilot/targeted_support2_wave20_checkpoint_20260821")
RETRIEVED_AT = "2026-08-21T04:15:00Z"
SOURCE_URL = "https://doi.org/10.5061/dryad.292q34fp"
SOURCE_XLS_SHA256 = "9ec35a0d56bc695c5792f5037cfd325f465a01737c7666f3eb042cc587f0fdc9"
SOURCE_README_SHA256 = "49f8283f2f52cee33d3e9fa429df24127572f53c2bddb7b636915c19c384fb9b"
SOURCE_CSV_SHA256 = "7e2a72856f85f9f7db25b9d3888248b575f4502a0d1713407592fb1c95d5fbff"

# Excel row, accepted species, expected family, trait, allowed source token,
# normalized value.  The row number makes the selection reproducible against
# the immutable original workbook and its committed normalized representation.
SELECTIONS = (
    (20, "Senecio squalidus", "Asteraceae", "autonomous_selfing_capacity", "auton", "autonomous"),
    (24, "Crepis sancta", "Asteraceae", "autonomous_selfing_capacity", "auton", "autonomous"),
    (56, "Begonia hirsuta", "Begoniaceae", "autonomous_selfing_capacity", "auton", "autonomous"),
    (57, "Begonia semiovata", "Begoniaceae", "autonomous_selfing_capacity", "auton", "autonomous"),
    (355, "Hordeum spontaneum", "Poaceae", "cleistogamy", "cleis", "facultative"),
    (359, "Avena barbata", "Poaceae", "cleistogamy", "cleis", "facultative"),
    (366, "Hordeum spontaneum", "Poaceae", "autonomous_selfing_capacity", "auton", "autonomous"),
    (367, "Deschampsia cespitosa", "Poaceae", "autonomous_selfing_capacity", "auton", "autonomous"),
    (
        403,
        "Grevillea barklyana",
        "Proteaceae",
        "autonomous_selfing_capacity",
        "auton",
        "autonomous",
    ),
)


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _text_sha256(path: Path) -> str:
    canonical = path.read_text(encoding="utf-8").replace("\r\n", "\n")
    return _sha(canonical)


def _normalized_species(value: str) -> str:
    return " ".join(value.split())


def _lineage(citation: str) -> str:
    normalized = " ".join(citation.casefold().split())
    return f"citation:{_sha(normalized)[:20]}"


def _source_rows() -> pd.DataFrame:
    if _file_sha256(SOURCE_XLS) != SOURCE_XLS_SHA256:
        raise ValueError(f"Goodwillie source hash mismatch: {SOURCE_XLS}")
    for path, expected in (
        (SOURCE_README, SOURCE_README_SHA256),
        (SOURCE_CSV, SOURCE_CSV_SHA256),
    ):
        if _text_sha256(path) != expected:
            raise ValueError(f"Goodwillie canonical text hash mismatch: {path}")
    rows = pd.read_csv(SOURCE_CSV, dtype=str).fillna("")
    if len(rows) != 469 or rows["source_row_number"].duplicated().any():
        raise ValueError("Goodwillie normalized source must contain 469 unique rows")
    return rows.set_index("source_row_number", drop=False)


def primary_rows() -> list[dict[str, str]]:
    """Return nine exact-master, natural-population reproductive records."""
    source = _source_rows()
    rows: list[dict[str, str]] = []
    for row_number, species, family, trait, token, value in SELECTIONS:
        record = source.loc[str(row_number)]
        source_species = _normalized_species(record["Genus species"])
        source_tokens = {
            part.strip().casefold() for part in record["Mech of selfing"].split(",") if part.strip()
        }
        if source_species != species or record["Family"] != family:
            raise ValueError(f"identity mismatch at Goodwillie row {row_number}")
        if record["PopType (0=natural, 1= experimental, 2=seed orchard, 3=agricultural)"] != "0":
            raise ValueError(f"non-natural population at Goodwillie row {row_number}")
        if token not in source_tokens:
            raise ValueError(f"missing explicit {token} token at row {row_number}")

        citation = record["Reference.1"].strip()
        lineage = _lineage(citation)
        excerpt = (
            "Goodwillie et al. Dryad structured record: "
            f"Genus species={species}; Mean-tm={record['Mean-tm']}; "
            f"Mech of selfing={record['Mech of selfing']}; PopType=0; "
            f"Where-t={record['Where-t']}; underlying reference={citation}"
        )
        evidence = _wave15_row(
            species,
            "reproductive_assurance",
            trait,
            f"Mech of selfing={record['Mech of selfing']}",
            value,
            "high",
            "Goodwillie, Kalisz & Eckert 2005/2012 Dryad mixed-mating database",
            SOURCE_URL,
            "The evolutionary enigma of mixed mating systems in plants: Dryad dataset",
            (
                "Goodwillie, Kalisz & Eckert (2005), Annual Review of Ecology, "
                "Evolution, and Systematics, DOI "
                "10.1146/annurev.ecolsys.36.091704.175539; underlying study: "
                f"{citation}"
            ),
            excerpt,
            lineage,
            "A",
            "published_research_database",
            "en",
            "Goodwillie 2005 Dryad mixed mating mode of selfing",
            status="wild_natural_population",
            content_sha256=SOURCE_XLS_SHA256,
            content_sha256_basis="original_dryad_xls_file_bytes",
        )
        evidence.update(
            {
                "source_group": SOURCE_GROUP,
                "source_record_id": f"goodwillie2005:row:{row_number}:{trait}",
                "lineage_method": ("normalized_underlying_study_citation_not_dryad_redistributor"),
                "domain": "datadryad.org",
                "retrieved_at_utc": RETRIEVED_AT,
            }
        )
        rows.append(evidence)
    return rows


def build_audit(evidence: pd.DataFrame) -> pd.DataFrame:
    audit = _wave15_build_audit(evidence)
    audit["reviewer"] = "Codex Goodwillie source-scale audit"
    audit["reviewed_at_utc"] = RETRIEVED_AT
    audit["decision_reason"] = (
        "Accepted from an exact fixed-master species row in the CC0 Dryad "
        "database: the Mode of selfing field is explicit, PopType is natural, "
        "and the underlying study is retained as the source lineage."
    )
    return audit.loc[:, AUDIT_COLUMNS]


def build_checkpoint(output_dir: Path = CHECKPOINT) -> dict[str, Any]:
    evidence = pd.DataFrame(primary_rows(), columns=EVIDENCE_COLUMNS).sort_values(
        ["accepted_species", "trait_name", "source_lineage"], kind="stable"
    )
    evidence = evidence.reset_index(drop=True)
    audit = build_audit(evidence)
    if len(evidence) != 9 or len(audit) != 9:
        raise ValueError("Goodwillie checkpoint must contain exactly nine rows")
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("Goodwillie candidate IDs must be unique")
    if evidence.duplicated(["accepted_species", "trait_name"]).any():
        raise ValueError("Goodwillie species x trait pairs must be unique")

    prior_evidence = pd.read_csv(
        PRIOR / "combined_curated_evidence_20260821.csv", dtype=str
    ).fillna("")
    prior_audit = pd.read_csv(
        PRIOR / "combined_curated_manual_audit_20260821.csv", dtype=str
    ).fillna("")
    combined_evidence = pd.concat([prior_evidence, evidence], ignore_index=True)
    combined_audit = pd.concat([prior_audit, audit], ignore_index=True)
    if combined_evidence["candidate_id"].duplicated().any():
        raise ValueError("combined evidence candidate IDs must be unique")
    if combined_audit["candidate_id"].duplicated().any():
        raise ValueError("combined audit candidate IDs must be unique")

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "evidence": output_dir / "goodwillie_2005_mode_of_selfing_evidence_20260821.csv",
        "audit": output_dir / "goodwillie_2005_mode_of_selfing_manual_audit_20260821.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260821.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260821.csv",
        "manifest": output_dir / "source_acquisition_manifest_goodwillie_2005.json",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    manifest = {
        "checkpoint": SOURCE_GROUP,
        "built_at_utc": RETRIEVED_AT,
        "baseline_formal_run_id": 32445414491,
        "dryad_dataset_doi": "10.5061/dryad.292q34fp",
        "dryad_file_id": 27237,
        "license": "CC0-1.0",
        "source_rows": len(_source_rows()),
        "exact_master_mating_system_rows": 243,
        "novel_mating_system_rows": 0,
        "accepted_evidence_rows": len(evidence),
        "accepted_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "accepted_species": int(evidence["accepted_species"].nunique()),
        "accepted_source_lineages": int(evidence["source_lineage"].nunique()),
        "accepted_by_trait": evidence["trait_name"].value_counts().sort_index().to_dict(),
        "recorded_queries": int(evidence["query"].nunique()),
        "formal_search_api_queries": 0,
        "search_cost_usd": 0.0,
        "local_formal_simulation": {
            "strict_species_axis_before": 159443,
            "strict_species_axis_after": 160103,
            "net_change": 660,
            "reproductive_assurance_before": 30398,
            "reproductive_assurance_after": 31058,
            "new_eligible_rule": "Begonia|autonomous_selfing_capacity",
        },
        "source_sha256": {
            SOURCE_XLS.name: SOURCE_XLS_SHA256,
            SOURCE_README.name: SOURCE_README_SHA256,
            SOURCE_CSV.name: SOURCE_CSV_SHA256,
        },
        "guardrails": {
            "outcrossing_rate_reinterpreted_as_other_trait": False,
            "none_mapped_to_absent_autonomous_selfing": False,
            "geitonogamy_mapped_to_autonomous_selfing": False,
            "family_inference": False,
            "global_fallback": False,
            "min_species_two_production": False,
            "cross_trait_substitution": False,
            "genus_axis_only_join": False,
        },
        "output_sha256": {
            key: _sha(path.read_text(encoding="utf-8"))
            for key, path in paths.items()
            if key != "manifest"
        },
        "notes": (
            "The +660 value is a local dry-run expectation, not a formal artifact "
            "claim. The GitHub all-evidence rebuild remains authoritative. The "
            "243 exact-master mating-system rows were already represented by prior "
            "lineages and were not duplicated."
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
