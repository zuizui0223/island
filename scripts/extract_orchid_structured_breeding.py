"""Extract strict reproductive candidates from the complete Ackerman orchid table.

The whole-source inventory identified explicit `SI`, `SC`, `Mixed mating`,
`autonomous selfing/agagamospermy`, and `evidence for selfing` columns.  Only
SI/SC are promoted into the reviewed source-scale packet here.

Ackerman et al. define SI/SC from controlled hand-pollination comparisons:
self-pollination with 0% fruit set or <5% seed set is SI, otherwise SC.  That
source-level definition is sufficiently specific for the strict
`self_incompatibility` ontology when the source row also retains a literature
reference.  By contrast, the database's `Mixed mating` column means a mixed
*pollination system* (chasmogamy plus autonomous pollination), not a population-
genetic mixed-mating estimate, so it is explicitly held out.  The composite
`autonomous selfing/agagamospermy` field is also held out because it conflates
sexual autonomous selfing with apomixis; `evidence for selfing` is provenance
about that composite state and is not substituted for a strict trait.

This script completes one source-scale packet; it does not itself change formal
coverage.  Candidates remain promotion_allowed=false until the packet is
combined with other completed source packets in a manual batch integration.
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
SOURCE_ARTICLE_URL = "https://academic.oup.com/botlinnean/article/202/3/295/7076252"
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

    # The source workbook prints the genus only on the first row of each genus.
    # Forward fill reconstructs the source row key; it is not taxonomic inference.
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
            "source_article_url": SOURCE_ARTICLE_URL,
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
            if refs:
                selected.append(
                    {
                        **base,
                        "trait_name": "self_incompatibility",
                        "normalized_value": value,
                        "quality": "high",
                        "source_column": "SI|SC",
                        "source_raw_value": f"SI={text(row['SI'])}; SC={text(row['SC'])}",
                        "review_status": "source_methodology_reviewed_reference_backed_strict_direct",
                        "acceptance_basis": (
                            "Ackerman et al. species SI/SC scoring is defined from controlled hand-"
                            "pollination comparisons; exact fixed-universe species; source references retained"
                        ),
                    }
                )
            else:
                holdouts.append(
                    {
                        **base,
                        "source_column": "SI|SC",
                        "source_raw_value": f"SI={text(row['SI'])}; SC={text(row['SC'])}",
                        "holdout_reason": "strict_si_sc_state_without_row_reference",
                    }
                )

        for col, reason in (
            (
                "Mixed mating",
                "source_mixed_pollination_system_not_population_genetic_mating_system",
            ),
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
        compression={"method": "gzip", "mtime": 0},
    )
    withheld.to_csv(
        args.output / "orchid_reproductive_holdouts.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    unresolved_candidates = candidates.loc[
        candidates["unresolved_reproductive_before_source"].eq("true")
    ] if not candidates.empty else candidates
    unresolved_candidates.to_csv(
        args.output / "orchid_unresolved_si_sc_batch.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    summary = {
        "contract": "ackerman_orchid_structured_breeding_v2",
        "source_doi": SOURCE_DOI,
        "source_article_url": SOURCE_ARTICLE_URL,
        "source_scale_complete": True,
        "source_rows": int(len(source)),
        "fixed_universe_species_rows_after_genus_filldown": matched_rows,
        "structured_rows_outside_exact_fixed_names": unmatched_structured,
        "strict_reference_backed_si_sc_rows": int(len(candidates)),
        "strict_reference_backed_si_sc_species": int(candidates["accepted_species"].nunique()) if not candidates.empty else 0,
        "unresolved_reproductive_si_sc_rows": int(len(unresolved_candidates)),
        "unresolved_reproductive_si_sc_species": int(unresolved_candidates["accepted_species"].nunique()) if not unresolved_candidates.empty else 0,
        "holdout_rows": int(len(withheld)),
        "candidate_counts_by_value": {
            str(value): int(n)
            for value, n in (
                candidates["normalized_value"].value_counts().items()
                if not candidates.empty else []
            )
        },
        "unresolved_counts_by_value": {
            str(value): int(n)
            for value, n in (
                unresolved_candidates["normalized_value"].value_counts().items()
                if not unresolved_candidates.empty else []
            )
        },
        "excluded_cross_trait_mappings": {
            "Mixed mating": "held out; mixed pollination system is not strict mating_system",
            "autonomous selfing/agagamospermy": "held out; autonomous sexual selfing and apomixis conflated",
            "evidence for selfing": "held out; evidence grade is not a trait state",
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
