"""Audit remaining strict reproductive Validated-Low headroom without relaxing rules.

The recovered Wave53 rule audit is the last hash-verified all-evidence rule
rebuild.  Restart/source-scale direct packets are deliberately marked
`genus_rule_training_allowed=false`, so they cannot create new genus rules in
this restart lane.  This audit therefore asks a narrower exact question: among
rules already eligible under `current_min3`, which still have unresolved valid
species-axis cells after the current batch?

No min2/relaxed setting, family inference, axis-only join, or global fallback is
used.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

AXIS = "reproductive_assurance"
BINOMIAL = re.compile(r"^[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+$")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--rule-audit", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    rules = pd.read_csv(args.rule_audit, dtype=str).fillna("")
    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    eligible = rules.loc[
        rules["setting"].eq("current_min3")
        & rules["axis"].eq(AXIS)
        & rules["eligible"].str.casefold().eq("true")
        & rules["diagnostic_only"].str.casefold().ne("true")
    ].copy()
    if eligible.empty:
        raise ValueError("no current_min3 reproductive rules found")

    rep = coverage.loc[coverage["axis"].eq(AXIS)].copy()
    if len(rep) != 106_295 or rep["accepted_species"].nunique() != 106_295:
        raise ValueError("coverage is not the fixed 106295-species reproductive universe")
    rep["genus"] = rep["accepted_species"].str.split().str[0]
    unresolved = rep.loc[rep["quality"].eq("")].copy()
    unresolved["exact_binomial"] = unresolved["accepted_species"].map(
        lambda value: bool(BINOMIAL.fullmatch(value))
    )

    candidate = unresolved.merge(
        eligible[
            [
                "genus",
                "trait_name",
                "inferred_value",
                "n_direct_species",
                "dominance",
                "species_loo_accuracy",
                "lineage_loo_accuracy",
                "required_min_species",
                "required_dominance",
                "required_masked_accuracy",
            ]
        ],
        on="genus",
        how="inner",
        validate="many_to_many",
    )
    valid = candidate.loc[candidate["exact_binomial"]].copy()
    candidate.to_csv(args.output / "strict_low_all_remaining_candidates.csv", index=False)
    valid.to_csv(args.output / "strict_low_valid_binomial_candidates.csv", index=False)

    summary = {
        "contract": "reproductive_strict_low_headroom_v1",
        "rule_setting": "current_min3",
        "thresholds_relaxed": False,
        "eligible_reproductive_genus_trait_rules": int(len(eligible)),
        "eligible_reproductive_genera": int(eligible["genus"].nunique()),
        "remaining_unresolved_reproductive_cells": int(len(unresolved)),
        "unresolved_cells_in_eligible_genera_all_labels": int(candidate["accepted_species"].nunique()),
        "valid_binomial_cells_in_eligible_genera": int(valid["accepted_species"].nunique()),
        "valid_binomial_species": sorted(valid["accepted_species"].unique().tolist()),
        "family_inference": False,
        "global_fallback": False,
        "axis_only_rules": False,
        "new_source_packets_train_rules": False,
        "claim_limit": (
            "Headroom under the last hash-verified all-evidence current_min3 rule audit; "
            "restart/source packets marked genus_rule_training_allowed=false do not create new rules."
        ),
    }
    (args.output / "strict_low_headroom_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
