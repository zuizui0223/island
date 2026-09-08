"""Extract strict SI/SC cells from Ackerman orchid database v1.4.0 (2025 update).

The 2025 Zenodo release is the same Ackerman et al. global orchid pollination
resource previously reviewed for the 2024 source batch, with the same creators,
research scope and structured breeding-system columns, updated through 2025.
Only explicit reference-backed SI/SC states are mapped to the strict
``self_incompatibility`` trait.  Mixed pollination, autonomous
selfing/agamospermy and generic evidence-for-selfing fields remain held out.

The script emits a reviewed source packet only for cells that are still empty in
the supplied cumulative checkpoint.  Formal coverage changes happen later in
the cumulative source-batch integrator.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

AXIS = "reproductive_assurance"
SOURCE_RECORD = 18208949
SOURCE_DOI = "10.5281/zenodo.18208949"
SOURCE_URL = "https://zenodo.org/api/records/18208949/files/Pollination%20List%20Thru%202025.xlsx/content"
SOURCE_TITLE = "Beyond the various contrivances by which orchids are pollinated: global patterns in orchid pollination biology"
SOURCE_VERSION = "1.4.0"
SOURCE_ARTICLE_URL = "https://academic.oup.com/botlinnean/article/202/3/295/7076252"
EXPECTED_CREATORS = {
    "James D. Ackerman", "Ryan D. Phillips", "Raymond L. Tremblay", "Adam Karremans",
    "Noushka Reiter", "Craig I. Peter", "Diego Bogarín", "Oscar A. Pérez-Escobar", "Hong Liu",
}
BINOMIAL = re.compile(r"^[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+$")
REVIEW_STATUS = "source_methodology_reviewed_reference_backed_strict_direct"


def text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def present(value: object) -> bool:
    return text(value).casefold() not in {"", "0", "na", "n/a", "none", "?", "nan"}


def lineage(refs: str) -> str:
    return "orchid-2025-reference-set:" + hashlib.sha256(refs.casefold().encode("utf-8")).hexdigest()[:24]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--workbook", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--zenodo-metadata", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    meta = json.loads(args.zenodo_metadata.read_text(encoding="utf-8"))
    md = meta.get("metadata", {})
    title = text(md.get("title"))
    version = text(md.get("version"))
    doi = text(md.get("doi")) or text(meta.get("doi"))
    creators = {text(row.get("name")) for row in md.get("creators", []) if isinstance(row, dict)}
    if title != SOURCE_TITLE or version != SOURCE_VERSION or doi != SOURCE_DOI or creators != EXPECTED_CREATORS:
        raise ValueError(
            f"Zenodo source identity drift: title={title!r}, version={version!r}, doi={doi!r}, creators={sorted(creators)}"
        )
    files = {text(row.get("key")): row for row in meta.get("files", []) if isinstance(row, dict)}
    source_file = files.get("Pollination List Thru 2025.xlsx")
    if not source_file:
        raise ValueError("Zenodo record lacks expected 2025 workbook")
    expected_md5 = text(source_file.get("checksum")).removeprefix("md5:")
    actual_md5 = hashlib.md5(args.workbook.read_bytes()).hexdigest()  # noqa: S324 - source integrity identifier
    if expected_md5 and actual_md5 != expected_md5:
        raise ValueError(f"workbook MD5 mismatch: {actual_md5} != {expected_md5}")

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    rep = coverage.loc[coverage["axis"].eq(AXIS), ["accepted_species", "quality"]].drop_duplicates()
    if len(rep) != 106_295 or rep["accepted_species"].duplicated().any():
        raise ValueError(f"expected fixed reproductive universe, got {len(rep)}")
    fixed = set(rep["accepted_species"])
    unresolved = set(rep.loc[rep["quality"].eq(""), "accepted_species"])

    book = pd.ExcelFile(args.workbook)
    if "species" not in book.sheet_names:
        raise ValueError(f"expected species sheet; sheets={book.sheet_names}")
    source = pd.read_excel(args.workbook, sheet_name="species", dtype=object).fillna("")
    required = {
        "genus", "species", "SI", "SC", "Mixed mating",
        "autonomous selfing/agagamospermy", "evidence for selfing", "references",
    }
    missing = required.difference(source.columns)
    if missing:
        raise ValueError(f"2025 orchid schema no longer matches reviewed contract: {sorted(missing)}")

    genus = source["genus"].map(text).replace("", pd.NA).ffill().fillna("")
    epithet = source["species"].map(text)
    source["source_species_name"] = (genus + " " + epithet).str.strip()
    source["source_excel_row"] = source.index + 2

    selected: list[dict[str, object]] = []
    holdouts: list[dict[str, object]] = []
    all_si_sc = 0
    fixed_si_sc = 0
    for raw in source.to_dict("records"):
        species = text(raw["source_species_name"])
        if not BINOMIAL.fullmatch(species):
            continue
        si = present(raw["SI"])
        sc = present(raw["SC"])
        if not (si or sc):
            continue
        all_si_sc += 1
        if species not in fixed:
            continue
        fixed_si_sc += 1
        refs = text(raw["references"])
        value = "mixed_or_variable" if si and sc else ("SI" if si else "SC")
        base = {
            "accepted_species": species,
            "source_species_name": species,
            "axis": AXIS,
            "trait_name": "self_incompatibility",
            "source_excel_row": int(raw["source_excel_row"]),
            "source_reference_raw": refs,
            "source_url": SOURCE_URL,
            "source_article_url": SOURCE_ARTICLE_URL,
            "source_lineage": lineage(refs) if refs else f"orchid-2025-row-without-reference:{int(raw['source_excel_row'])}",
            "name_match_method": "exact_fixed_universe_name_after_source_genus_filldown",
            "source_column": "SI|SC",
            "source_raw_value": f"SI={text(raw['SI'])}; SC={text(raw['SC'])}",
            "normalized_value": value,
            "quality": "high",
            "promotion_allowed": "false",
            "genus_rule_training_allowed": "false",
            "source_version": SOURCE_VERSION,
            "source_record": SOURCE_DOI,
        }
        if species not in unresolved:
            holdouts.append({**base, "holdout_reason": "reproductive_axis_already_resolved"})
            continue
        if not refs:
            holdouts.append({**base, "holdout_reason": "strict_si_sc_state_without_row_reference"})
            continue
        selected.append({
            **base,
            "review_status": REVIEW_STATUS,
            "acceptance_basis": (
                "Ackerman global orchid pollination database v1.4.0 is the 2025 update of the previously reviewed source; "
                "same creators, title/resource lineage and SI/SC schema; exact species; reference-backed explicit SI/SC"
            ),
        })

    packet = pd.DataFrame(selected)
    if packet.empty:
        packet = pd.DataFrame(columns=[
            "accepted_species", "source_species_name", "axis", "trait_name", "source_excel_row",
            "source_reference_raw", "source_url", "source_article_url", "source_lineage",
            "name_match_method", "source_column", "source_raw_value", "normalized_value", "quality",
            "promotion_allowed", "genus_rule_training_allowed", "source_version", "source_record",
            "review_status", "acceptance_basis",
        ])
    if packet.duplicated(["accepted_species", "trait_name"]).any():
        dup = packet.loc[packet.duplicated(["accepted_species", "trait_name"], keep=False), "accepted_species"].tolist()
        raise ValueError(f"duplicate current-gap species-trait rows: {dup}")
    packet.to_csv(
        args.output / "orchid_2025_unresolved_si_sc_batch.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    pd.DataFrame(holdouts).to_csv(
        args.output / "orchid_2025_si_sc_holdouts.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    summary = {
        "contract": "ackerman_orchid_2025_structured_breeding_v1",
        "source_record": SOURCE_RECORD,
        "source_doi": SOURCE_DOI,
        "source_version": SOURCE_VERSION,
        "source_title": SOURCE_TITLE,
        "source_workbook_md5": actual_md5,
        "source_rows": int(len(source)),
        "si_sc_source_species_rows": all_si_sc,
        "si_sc_exact_fixed_species_rows": fixed_si_sc,
        "reviewed_current_gap_rows": int(len(packet)),
        "reviewed_current_gap_species": int(packet["accepted_species"].nunique()) if len(packet) else 0,
        "reviewed_current_gap_state_counts": packet["normalized_value"].value_counts().sort_index().astype(int).to_dict() if len(packet) else {},
        "reviewed_current_gap_species_names": sorted(packet["accepted_species"].tolist()) if len(packet) else [],
        "source_identity_equivalence": {
            "same_database_title": True,
            "same_creator_set": True,
            "same_required_si_sc_schema": True,
            "version": SOURCE_VERSION,
            "update_scope": "literature and pollination list updated through 2025",
        },
        "excluded_cross_trait_mappings": {
            "Mixed mating": "not mapped to population-genetic mating_system",
            "autonomous selfing/agagamospermy": "not mapped because sexual selfing and apomixis are conflated",
            "evidence for selfing": "evidence field, not a strict trait state",
        },
        "formal_gain": 0,
        "promotion_allowed": False,
        "genus_rule_training_allowed": False,
        "integration_deferred_to_source_batch": True,
    }
    (args.output / "orchid_2025_structured_breeding_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
