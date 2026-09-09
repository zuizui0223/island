"""Measure current reproductive ROI of the complete Goodwillie Dryad dataset.

This is a source-selection audit, not a promotion step.  It reads all 469 source
rows and measures exact overlap with the current recoverable public
reproductive-assurance gaps after the Orchid source batch.

`Mean-tm` is retained as the published multilocus outcrossing estimate.  For ROI
screening only, the script also reports the same <0.2 / 0.2-0.8 / >0.8 bins used
by the repository's already-reviewed Moeller source package.  Those diagnostic
bins are NOT emitted as strict evidence here: Goodwillie's source-level
measurement/aggregation contract must be reviewed before any `mating_system`
promotion.  Explicit `auton` and `cleis` tokens are reported separately because
they identify different reproductive traits and must never be substituted for
Mean-tm.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

AXIS = "reproductive_assurance"
BINOMIAL = re.compile(r"^[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+$")
POP = "PopType (0=natural, 1= experimental, 2=seed orchard, 3=agricultural)"


def text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def tm_class(value: float) -> str:
    if value < 0.2:
        return "predominantly_selfing"
    if value <= 0.8:
        return "mixed_mating"
    return "predominantly_outcrossing"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    source = pd.read_csv(args.source, dtype=str).fillna("")
    required = {"source_row_number", "Genus species", POP, "Mean-tm", "Mech of selfing", "Reference.1"}
    missing = required.difference(source.columns)
    if missing:
        raise ValueError(f"Goodwillie source missing columns: {sorted(missing)}")
    if len(source) != 469:
        raise ValueError(f"expected 469 Goodwillie rows, got {len(source)}")

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    rep = coverage.loc[coverage["axis"].eq(AXIS), ["accepted_species", "quality"]].drop_duplicates()
    if len(rep) != 106_295 or rep["accepted_species"].duplicated().any():
        raise ValueError("coverage is not the fixed 106295-species reproductive universe")
    quality = dict(zip(rep["accepted_species"], rep["quality"], strict=True))
    fixed = set(quality)

    work = source.copy()
    work["source_species_name"] = work["Genus species"].map(text)
    work["exact_binomial"] = work["source_species_name"].map(lambda x: bool(BINOMIAL.fullmatch(x)))
    work["accepted_species"] = work["source_species_name"]
    work["in_fixed_universe_exact"] = work["accepted_species"].isin(fixed)
    work["current_reproductive_quality"] = work["accepted_species"].map(quality).fillna("")
    work["currently_unresolved_reproductive"] = (
        work["in_fixed_universe_exact"] & work["current_reproductive_quality"].eq("")
    )
    work["natural_population"] = work[POP].map(text).eq("0")
    work["underlying_reference_present"] = work["Reference.1"].map(text).ne("")
    work["mean_tm"] = pd.to_numeric(work["Mean-tm"], errors="coerce")
    work["mean_tm_valid"] = work["mean_tm"].between(0, 1, inclusive="both")
    work["diagnostic_tm_class"] = ""
    valid_tm = work["mean_tm_valid"]
    work.loc[valid_tm, "diagnostic_tm_class"] = work.loc[valid_tm, "mean_tm"].map(tm_class)
    work["selfing_tokens"] = work["Mech of selfing"].map(text).str.casefold()
    work["explicit_autonomous"] = work["selfing_tokens"].map(
        lambda x: "auton" in {part.strip() for part in x.split(",") if part.strip()}
    )
    work["explicit_cleistogamy"] = work["selfing_tokens"].map(
        lambda x: "cleis" in {part.strip() for part in x.split(",") if part.strip()}
    )

    roi_rows = work.loc[
        work["currently_unresolved_reproductive"]
        & work["natural_population"]
        & work["underlying_reference_present"]
    ].copy()
    tm_rows = roi_rows.loc[roi_rows["mean_tm_valid"]].copy()

    # Multiple study rows for one species are kept for provenance.  For source
    # selection, classify a species as concordant only when all natural-study
    # rows fall in the same diagnostic bin.  Conflicts are explicitly counted.
    species_tm = []
    for species, group in tm_rows.groupby("accepted_species"):
        classes = sorted(set(group["diagnostic_tm_class"]) - {""})
        species_tm.append(
            {
                "accepted_species": species,
                "n_natural_study_rows": int(len(group)),
                "n_underlying_references": int(group["Reference.1"].map(text).nunique()),
                "diagnostic_classes": "|".join(classes),
                "diagnostic_class_count": len(classes),
                "diagnostic_concordant": len(classes) == 1,
                "diagnostic_candidate_class": classes[0] if len(classes) == 1 else "",
                "strict_promotion_allowed": False,
                "review_needed": "source_level_mean_tm_to_mating_system_contract",
            }
        )
    species_tm_frame = pd.DataFrame(species_tm)

    explicit = roi_rows.loc[roi_rows["explicit_autonomous"] | roi_rows["explicit_cleistogamy"]].copy()
    explicit_rows = []
    for row in explicit.to_dict("records"):
        if bool(row["explicit_autonomous"]):
            explicit_rows.append(
                {
                    "accepted_species": row["accepted_species"],
                    "trait_name": "autonomous_selfing_capacity",
                    "source_token": "auton",
                    "normalized_value": "autonomous",
                    "source_row_number": row["source_row_number"],
                    "underlying_reference": text(row["Reference.1"]),
                    "strict_promotion_allowed": False,
                    "review_needed": "current_cell_and_duplicate_lineage_review",
                }
            )
        if bool(row["explicit_cleistogamy"]):
            explicit_rows.append(
                {
                    "accepted_species": row["accepted_species"],
                    "trait_name": "cleistogamy",
                    "source_token": "cleis",
                    "normalized_value": "facultative",
                    "source_row_number": row["source_row_number"],
                    "underlying_reference": text(row["Reference.1"]),
                    "strict_promotion_allowed": False,
                    "review_needed": "current_cell_and_duplicate_lineage_review",
                }
            )
    explicit_frame = pd.DataFrame(explicit_rows)

    work.to_csv(
        args.output / "goodwillie_full_source_roi_audit.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    species_tm_frame.to_csv(args.output / "goodwillie_unresolved_tm_species_roi.csv", index=False)
    explicit_frame.to_csv(args.output / "goodwillie_unresolved_explicit_mechanism_roi.csv", index=False)

    concordant = species_tm_frame.loc[species_tm_frame["diagnostic_concordant"]].copy() if not species_tm_frame.empty else species_tm_frame
    class_counts = (
        concordant["diagnostic_candidate_class"].value_counts().sort_index().astype(int).to_dict()
        if not concordant.empty else {}
    )
    summary = {
        "contract": "goodwillie_source_scale_roi_v1",
        "source_rows": int(len(source)),
        "current_reproductive_filled": int(rep["quality"].ne("").sum()),
        "current_reproductive_unresolved": int(rep["quality"].eq("").sum()),
        "exact_fixed_universe_source_species": int(work.loc[work["in_fixed_universe_exact"], "accepted_species"].nunique()),
        "exact_current_unresolved_source_species": int(work.loc[work["currently_unresolved_reproductive"], "accepted_species"].nunique()),
        "unresolved_natural_reference_backed_species": int(roi_rows["accepted_species"].nunique()),
        "unresolved_natural_reference_backed_mean_tm_rows": int(len(tm_rows)),
        "unresolved_natural_reference_backed_mean_tm_species": int(tm_rows["accepted_species"].nunique()),
        "diagnostic_concordant_mean_tm_species": int(len(concordant)),
        "diagnostic_conflicting_mean_tm_species": int((species_tm_frame["diagnostic_class_count"] > 1).sum()) if not species_tm_frame.empty else 0,
        "diagnostic_class_counts": class_counts,
        "unresolved_explicit_autonomous_or_cleistogamy_rows": int(len(explicit_frame)),
        "unresolved_explicit_autonomous_or_cleistogamy_species": int(explicit_frame["accepted_species"].nunique()) if not explicit_frame.empty else 0,
        "formal_gain": 0,
        "mean_tm_strict_promotion_allowed": False,
        "thresholds_used_only_for_roi_diagnostic": {"predominantly_selfing": "tm<0.2", "mixed_mating": "0.2<=tm<=0.8", "predominantly_outcrossing": "tm>0.8"},
        "next_gate": "review Goodwillie source-level Mean-tm measurement/aggregation contract before any strict mating_system mapping",
    }
    (args.output / "goodwillie_source_scale_roi_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
