"""Extract strict reproductive candidates from the complete Ackerman orchid table.

This second source-scale stage is intentionally schema-specific after the whole
workbook inventory identified explicit `SI`, `SC`, and `Mixed mating` columns.
The `autonomous selfing/agagamospermy` field is retained only as a holdout
because it conflates sexual selfing with apomixis.  Likewise, `evidence for
selfing` is not substituted for any stricter reproductive trait.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

AXIS = "reproductive_assurance"
SOURCE_DOI = "10.5281/zenodo.14601785"
SOURCE_URL = "https://zenodo.org/records/14601785/files/Pollination%20List%20Thru%202024.xlsx?download=1"
BINOMIAL = re.compile(r"^[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+$")


def text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def digest(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()[:24]


def present(value: object) -> bool:
    return text(value).casefold() not in {"", "0", "na", "n/a", "none", "?", "nan"}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--workbook", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    rep = coverage.loc[coverage["axis"].eq(AXIS)].copy()
    if len(rep) != 106_295 or rep["accepted_species"].nunique() != 106_295:
        raise ValueError("expected one reproductive row for each fixed species")
    fixed = set(rep["accepted_species"])
    unresolved = set(rep.loc[rep["quality"].eq(""), "accepted_species"])

    source = pd.read_excel(args.workbook, sheet_name="species", dtype=object).fillna("")
    required = {
        "genus",
        "species",
        "SI",
        "SC",
        "Mixed mating",
        "autonomous selfing/agagamospermy",
        "evidence for selfing",
        "references",
    }
    missing = required.difference(source.columns)
    if missing:
        raise ValueError(f"orchid source missing structured columns: {sorted(missing)}")

    # The published table prints the genus only on the first row of each genus.
    # Forward fill is therefore part of reconstructing the source's row key, not
    # a taxonomic inference.
    genus = source["genus"].map(text).replace("", pd.NA).ffill().fillna("")
    epithet = source["species"].map(text)
    source["source_species_name"] = (genus + " " + epithet).str.strip()
    source["source_excel_row"] = source.index + 2

    selected: list[dict] = []
    holdouts: list[dict] = []
    matched_rows = 0
    unmatched_structured = 0

    for row in source.to_dict("records"):
        species = text(row["source_species_name"])
        if not BINOMIAL.fullmatch(species):
            continue
        in_fixed = species in fixed
        structured_present = any(
            present(row[col])
            for col in (
                "SI",
                "SC",
                "Mixed mating",
                "autonomous selfing/agagamospermy",
                "evidence for selfing",
            )
        )
        if structured_present and not in_fixed:
            unmatched_structured += 1
        if not in_fixed:
            continue
        matched_rows += 1
        refs = text(row["references"])
        base = {
            "accepted_species": species,
            "source_species_name": species,
            "axis": AXIS,
            "source_excel_row": int(row["source_excel_row"]),
            "source_reference_raw": refs,
            "source_url": SOURCE_URL,
            "source_lineage": (
                f"orchid-reference-set:{digest(refs.casefold())}"
                if refs
                else f"orchid-row-without-reference:{int(row['source_excel_row'])}"
            ),
            "name_match_method": "exact_fixed_universe_name_after_source_genus_filldown",
            "unresolved_reproductive_before_source": str(species in unresolved).lower(),
            "promotion_allowed": "false",
            "genus_rule_training_allowed": "false",
        }

        si = present(row["SI"])
        sc = present(row["SC"])
        if si or sc:
            value = "mixed_or_variable" if si and sc else ("SI" if si else "SC")
            selected.append(
                {
                    **base,
                    "trait_name": "self_incompatibility",
                    "normalized_value": value,
                    "quality": "unreviewed",
                    "source_column": "SI|SC",
                    "source_raw_value": f"SI={text(row['SI'])}; SC={text(row['SC'])}",
                    "review_status": (
                        "structured_direct_with_reference_needs_batch_review"
                        if refs
                        else "structured_direct_missing_reference_holdout"
                    ),
                }
            )
        if present(row["Mixed mating"]):
            selected.append(
                {
                    **base,
                    "trait_name": "mating_system",
                    "normalized_value": "mixed_mating",
                    "quality": "unreviewed",
                    "source_column": "Mixed mating",
                    "source_raw_value": text(row["Mixed mating"]),
                    "review_status": (
                        "structured_direct_with_reference_needs_batch_review"
                        if refs
                        else "structured_direct_missing_reference_holdout"
                    ),
                }
            )

        for col, reason in (
            (
                "autonomous selfing/agagamospermy",
                "conflates_autonomous_selfing_and_agamospermy_no_strict_mapping",
            ),
            (
                "evidence for selfing",
                "selfing_evidence_does_not_identify_autonomy_or_mating_system",
            ),
        ):
            if present(row[col]):
                holdouts.append(
                    {
                        **base,
                        "source_column": col,
                        "source_raw_value": text(row[col]),
                        "holdout_reason": reason,
                    }
                )

    candidates = pd.DataFrame(selected)
    withheld = pd.DataFrame(holdouts)
    candidates.to_csv(
        args.output / "orchid_structured_reproductive_review_queue.csv.gz",
        index=False,
        compression="gzip",
    )
    withheld.to_csv(
        args.output / "orchid_reproductive_holdouts.csv.gz",
        index=False,
        compression="gzip",
    )

    unresolved_candidates = candidates.loc[
        candidates["unresolved_reproductive_before_source"].eq("true")
    ] if not candidates.empty else candidates
    unresolved_referenced = unresolved_candidates.loc[
        unresolved_candidates["source_reference_raw"].ne("")
    ] if not unresolved_candidates.empty else unresolved_candidates
    summary = {
        "contract": "ackerman_orchid_structured_breeding_v1",
        "source_doi": SOURCE_DOI,
        "source_rows": int(len(source)),
        "fixed_universe_species_rows_after_genus_filldown": matched_rows,
        "structured_rows_outside_exact_fixed_names": unmatched_structured,
        "strict_candidate_rows": int(len(candidates)),
        "strict_candidate_species": int(candidates["accepted_species"].nunique()) if not candidates.empty else 0,
        "unresolved_reproductive_candidate_rows": int(len(unresolved_candidates)),
        "unresolved_reproductive_candidate_species": int(unresolved_candidates["accepted_species"].nunique()) if not unresolved_candidates.empty else 0,
        "unresolved_reproductive_candidate_species_with_reference": int(unresolved_referenced["accepted_species"].nunique()) if not unresolved_referenced.empty else 0,
        "holdout_rows": int(len(withheld)),
        "candidate_counts_by_trait_value": {
            f"{trait}:{value}": int(n)
            for (trait, value), n in (
                candidates.groupby(["trait_name", "normalized_value"]).size().items()
                if not candidates.empty
                else []
            )
        },
        "unresolved_counts_by_trait_value": {
            f"{trait}:{value}": int(n)
            for (trait, value), n in (
                unresolved_candidates.groupby(["trait_name", "normalized_value"]).size().items()
                if not unresolved_candidates.empty
                else []
            )
        },
        "formal_gain": 0,
        "promotion_allowed": False,
        "genus_rule_training_allowed": False,
        "integration_deferred_to_source_batch": True,
    }
    (args.output / "orchid_structured_breeding_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
