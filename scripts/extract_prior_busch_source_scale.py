"""Source-scale audit of the Prior & Busch rows retained in Meyer et al.'s input workbook.

This script deliberately does not promote evidence. Meyer et al.'s public analysis
sets the Prior & Busch `mating_system` field to NA before its final SI/SC export,
so this lane recovers the complete underlying row block and measures its overlap
with the current unresolved reproductive-assurance universe.

`quant` is preserved verbatim and numerically summarized. A provisional
mating-system bin is emitted only under the explicit conditional label
`if_quant_is_selfing_rate`; promotion remains false until the original
Prior & Busch variable semantics and species-level aggregation contract are
verified.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

BINOMIAL = re.compile(r"^[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+$")
AXIS = "reproductive_assurance"
SOURCE_REPO = "elenacicada/interesting_flowers"
SOURCE_COMMIT = "3d57315fe77097fe582025b884917550139801e8"
SOURCE_BLOB = "6c1f7de0f3e5a6c62d1fca2647d7adf33a19c380"
SOURCE_ARTIFACT = "Input_Data_MS.xlsx"
SOURCE_DOI = "10.5061/dryad.cc2fqz6hr"
UNDERLYING_DOI = "10.5061/dryad.3j9kd51jw"
UNDERLYING_ARTICLE_DOI = "10.1002/ajb2.1766"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def locate_prior_busch_sheet(path: Path) -> tuple[str, pd.DataFrame]:
    book = pd.ExcelFile(path)
    hits: list[tuple[str, pd.DataFrame]] = []
    for sheet in book.sheet_names:
        frame = pd.read_excel(path, sheet_name=sheet, dtype=object).fillna("")
        lookup = {str(c).strip().casefold(): c for c in frame.columns}
        if "ref" not in lookup or "genus_species" not in lookup:
            continue
        ref_col = lookup["ref"]
        subset = frame.loc[
            frame[ref_col].map(text).str.casefold().str.contains("prior and busch", regex=False)
        ].copy()
        if not subset.empty:
            hits.append((sheet, subset))
    if len(hits) != 1:
        raise ValueError(f"expected exactly one Prior & Busch source sheet, found {[h[0] for h in hits]}")
    return hits[0]


def provisional_selfing_bin(value: float) -> str:
    # Conditional only: valid iff source semantics confirm quant = selfing rate.
    if value < 0.2:
        return "predominantly_outcrossing"
    if value > 0.8:
        return "predominantly_selfing"
    return "mixed_mating"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-xlsx", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    sheet, source = locate_prior_busch_sheet(args.input_xlsx)
    col = {str(c).strip().casefold(): c for c in source.columns}
    if "quant" not in col:
        raise ValueError(f"Prior & Busch block lacks quant; columns={list(source.columns)}")
    species_col = col["genus_species"]
    quant_col = col["quant"]
    ref_col = col["ref"]

    source = source.copy()
    source["source_sheet"] = sheet
    source["source_row_in_sheet"] = source.index + 2
    source["accepted_species_candidate"] = source[species_col].map(text).str.replace("_", " ", regex=False)
    source["quant_raw"] = source[quant_col].map(text)
    source["quant_numeric"] = pd.to_numeric(source[quant_col], errors="coerce")
    source["exact_binomial"] = source["accepted_species_candidate"].map(lambda x: bool(BINOMIAL.fullmatch(x)))
    source["quant_in_unit_interval"] = source["quant_numeric"].between(0, 1, inclusive="both")
    source["provisional_mating_system_if_quant_is_selfing_rate"] = source["quant_numeric"].map(
        lambda x: provisional_selfing_bin(float(x)) if pd.notna(x) and 0 <= float(x) <= 1 else ""
    )
    source["trait_name"] = "mating_system"
    source["axis"] = AXIS
    source["source_lineage"] = source.apply(
        lambda r: f"prior-busch-2021:{text(r.get('citation', '')) or text(r.get('source', '')) or int(r['source_row_in_sheet'])}",
        axis=1,
    )
    source["evidence_scope"] = "source_scale_population_quantitative_candidate"
    source["normalized_value"] = ""
    source["evidence_quality"] = "unreviewed"
    source["promotion_allowed"] = False
    source["genus_rule_training_allowed"] = False
    source["review_status"] = "needs_original_variable_semantics_and_species_aggregation_review"

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    rep = coverage.loc[coverage["axis"].eq(AXIS), ["accepted_species", "quality"]].drop_duplicates()
    if len(rep) != 106295 or rep["accepted_species"].duplicated().any():
        raise ValueError(f"expected fixed 106295-species reproductive coverage, got {len(rep)}")
    quality = dict(zip(rep["accepted_species"], rep["quality"], strict=True))
    universe = set(quality)
    source["in_fixed_universe_exact"] = source["accepted_species_candidate"].isin(universe)
    source["current_reproductive_quality"] = source["accepted_species_candidate"].map(quality).fillna("")
    source["currently_unresolved_reproductive"] = source["in_fixed_universe_exact"] & source["current_reproductive_quality"].eq("")

    # A species-level candidate is only considered stable if every available
    # population/source row falls in the same provisional bin. This is still not
    # promotion-ready because quant semantics are intentionally unconfirmed here.
    eligible_rows = source.loc[
        source["currently_unresolved_reproductive"]
        & source["exact_binomial"]
        & source["quant_in_unit_interval"]
    ].copy()
    species_rows: list[dict[str, object]] = []
    for species, group in eligible_rows.groupby("accepted_species_candidate", sort=True):
        states = sorted(set(group["provisional_mating_system_if_quant_is_selfing_rate"]) - {""})
        values = sorted(float(v) for v in group["quant_numeric"].dropna())
        species_rows.append({
            "accepted_species": species,
            "axis": AXIS,
            "trait_name": "mating_system",
            "source_rows": len(group),
            "quant_min": min(values) if values else "",
            "quant_max": max(values) if values else "",
            "quant_mean": sum(values) / len(values) if values else "",
            "provisional_state_set_if_quant_is_selfing_rate": "|".join(states),
            "stable_single_state": len(states) == 1,
            "provisional_normalized_value": states[0] if len(states) == 1 else "",
            "promotion_allowed": False,
            "genus_rule_training_allowed": False,
            "review_status": "needs_original_variable_semantics_and_species_aggregation_review",
        })
    species = pd.DataFrame(species_rows)

    source.to_csv(args.output / "prior_busch_full_source_rows.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    species.to_csv(args.output / "prior_busch_unresolved_species_candidates.csv", index=False)

    numeric = source["quant_numeric"].dropna()
    summary = {
        "contract": "prior_busch_2021_source_scale_audit_v1",
        "source_repo": SOURCE_REPO,
        "source_commit": SOURCE_COMMIT,
        "source_blob": SOURCE_BLOB,
        "source_artifact": SOURCE_ARTIFACT,
        "source_sha256": sha256(args.input_xlsx),
        "source_doi": SOURCE_DOI,
        "underlying_dataset_doi": UNDERLYING_DOI,
        "underlying_article_doi": UNDERLYING_ARTICLE_DOI,
        "source_sheet": sheet,
        "source_columns": [str(c) for c in source.columns],
        "prior_busch_rows": int(len(source)),
        "raw_unique_taxon_labels": int(source["accepted_species_candidate"].nunique()),
        "exact_binomial_rows": int(source["exact_binomial"].sum()),
        "exact_fixed_universe_species": int(source.loc[source["in_fixed_universe_exact"], "accepted_species_candidate"].nunique()),
        "exact_current_unresolved_reproductive_species": int(source.loc[source["currently_unresolved_reproductive"], "accepted_species_candidate"].nunique()),
        "numeric_quant_rows": int(numeric.size),
        "quant_in_unit_interval_rows": int(source["quant_in_unit_interval"].sum()),
        "quant_min": float(numeric.min()) if len(numeric) else None,
        "quant_max": float(numeric.max()) if len(numeric) else None,
        "unresolved_species_with_numeric_quant": int(species["accepted_species"].nunique()) if len(species) else 0,
        "stable_single_state_unresolved_species_if_quant_is_selfing_rate": int(species["stable_single_state"].sum()) if len(species) else 0,
        "variable_state_unresolved_species_if_quant_is_selfing_rate": int((~species["stable_single_state"]).sum()) if len(species) else 0,
        "formal_gain": 0,
        "promotion_allowed": False,
        "semantics_gate": "verify Prior & Busch quant is selfing rate and define species aggregation before any mapping",
    }
    (args.output / "prior_busch_source_scale_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
