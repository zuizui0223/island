from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import pandas as pd

AXIS = "reproductive_assurance"
MIN_SPECIES = 3
MIN_DOMINANCE = 0.95
MIN_MASKED_ACCURACY = 0.85
EXCLUDED_SOURCE_GROUP = "reviewed_restart_20260907"


def _mode(values: list[str]) -> str | None:
    if not values:
        return None
    counts = Counter(values)
    top = counts.most_common()
    if len(top) > 1 and top[0][1] == top[1][1]:
        return None
    return top[0][0]


def _species_loo(rows: pd.DataFrame) -> tuple[int, int, float]:
    species_state = dict(zip(rows["accepted_species"], rows["normalized_value"], strict=True))
    correct = 0
    tested = 0
    for species, state in species_state.items():
        train = [v for s, v in species_state.items() if s != species]
        pred = _mode(train)
        if pred is None:
            continue
        tested += 1
        correct += int(pred == state)
    return tested, correct, (correct / tested if tested else 0.0)


def _split_lineages(value: str) -> set[str]:
    return {token.strip() for token in str(value).split("|") if token.strip()}


def _lineage_loo_proxy(rows: pd.DataFrame) -> tuple[int, int, float]:
    records = [
        (row.accepted_species, row.normalized_value, _split_lineages(row.source_lineages))
        for row in rows.itertuples(index=False)
    ]
    lineages = sorted({lineage for _, _, ls in records for lineage in ls})
    tested = 0
    correct = 0
    for lineage in lineages:
        held = [(s, v) for s, v, ls in records if lineage in ls and len(ls) == 1]
        if not held:
            continue
        train = [(s, v) for s, v, ls in records if lineage not in ls or len(ls) > 1]
        pred = _mode([v for _, v in train])
        if pred is None:
            continue
        tested += len(held)
        correct += sum(int(v == pred) for _, v in held)
    return tested, correct, (correct / tested if tested else 0.0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--direct", type=Path, required=True)
    parser.add_argument("--coverage", type=Path, required=True)
    parser.add_argument("--baseline-rule-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    direct = pd.read_csv(args.direct, dtype=str).fillna("")
    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    baseline = pd.read_csv(args.baseline_rule_audit, dtype=str).fillna("")

    direct = direct.loc[
        direct["axis"].eq(AXIS)
        & direct["resolution_status"].eq("resolved")
        & direct["quality"].isin({"high", "medium"})
        & direct["normalized_value"].ne("")
    ].copy()
    # PR151's new reviewed rows explicitly prohibit genus-rule training. They are
    # valid direct evidence but must not leak into a Low rule through this audit.
    direct = direct.loc[~direct["source_groups"].str.contains(EXCLUDED_SOURCE_GROUP, regex=False)].copy()
    direct["genus"] = direct["accepted_species"].str.split().str[0]

    unresolved = coverage.loc[
        coverage["axis"].eq(AXIS) & coverage["quality"].eq(""), "accepted_species"
    ].drop_duplicates()
    unresolved_by_genus = unresolved.str.split().str[0].value_counts().to_dict()

    baseline_current = baseline.loc[
        baseline["setting"].eq("current_min3")
        & baseline["axis"].eq(AXIS)
        & baseline["eligible"].astype(str).str.casefold().eq("true")
    ].copy()
    baseline_keys = set(zip(baseline_current["genus"], baseline_current["trait_name"], strict=True))

    rows: list[dict[str, object]] = []
    for (genus, trait), group in direct.groupby(["genus", "trait_name"], sort=True):
        group = group.drop_duplicates("accepted_species")
        n_species = group["accepted_species"].nunique()
        if n_species < MIN_SPECIES:
            continue
        counts = group["normalized_value"].value_counts()
        dominant = str(counts.index[0])
        dominant_n = int(counts.iloc[0])
        dominance = dominant_n / n_species
        species_n, species_correct, species_acc = _species_loo(group)
        lineage_n, lineage_correct, lineage_acc = _lineage_loo_proxy(group)
        eligible = (
            n_species >= MIN_SPECIES
            and dominance >= MIN_DOMINANCE
            and species_n > 0
            and species_acc >= MIN_MASKED_ACCURACY
            and lineage_n > 0
            and lineage_acc >= MIN_MASKED_ACCURACY
        )
        key = (genus, trait)
        rows.append(
            {
                "genus": genus,
                "trait_name": trait,
                "n_direct_species": n_species,
                "dominant_state": dominant,
                "counterexample_species": n_species - dominant_n,
                "dominance": round(dominance, 6),
                "species_loo_n": species_n,
                "species_loo_correct": species_correct,
                "species_loo_accuracy": round(species_acc, 6),
                "lineage_loo_proxy_n": lineage_n,
                "lineage_loo_proxy_correct": lineage_correct,
                "lineage_loo_proxy_accuracy": round(lineage_acc, 6),
                "current_threshold_candidate": eligible,
                "already_eligible_in_recovered_wave53_audit": key in baseline_keys,
                "current_unresolved_species_axis_in_genus": int(unresolved_by_genus.get(genus, 0)),
                "potential_new_low_species_axis_upper_bound": (
                    int(unresolved_by_genus.get(genus, 0)) if eligible and key not in baseline_keys else 0
                ),
                "promotion_allowed": False,
                "audit_note": (
                    "Resolved-direct proxy only; final promotion requires the shared raw-lineage auditor. "
                    "PR151 restart rows are excluded because all have genus_rule_training_allowed=false."
                ),
            }
        )

    audit = pd.DataFrame(rows)
    if audit.empty:
        audit = pd.DataFrame(
            columns=[
                "genus", "trait_name", "n_direct_species", "dominant_state",
                "counterexample_species", "dominance", "species_loo_n",
                "species_loo_correct", "species_loo_accuracy", "lineage_loo_proxy_n",
                "lineage_loo_proxy_correct", "lineage_loo_proxy_accuracy",
                "current_threshold_candidate", "already_eligible_in_recovered_wave53_audit",
                "current_unresolved_species_axis_in_genus", "potential_new_low_species_axis_upper_bound",
                "promotion_allowed", "audit_note",
            ]
        )
    audit = audit.sort_values(
        ["potential_new_low_species_axis_upper_bound", "n_direct_species", "genus", "trait_name"],
        ascending=[False, False, True, True],
    )
    audit.to_csv(args.output / "reproductive_low_capacity_audit.csv", index=False)
    candidate = audit.loc[
        audit["current_threshold_candidate"].astype(str).str.casefold().eq("true")
        & ~audit["already_eligible_in_recovered_wave53_audit"].astype(str).str.casefold().eq("true")
    ].copy()
    candidate.to_csv(args.output / "reproductive_new_low_rule_candidates.csv", index=False)

    summary = {
        "contract": "reproductive_current_threshold_low_capacity_audit_v1",
        "thresholds": {
            "min_species": MIN_SPECIES,
            "dominance": MIN_DOMINANCE,
            "masked_accuracy": MIN_MASKED_ACCURACY,
        },
        "excluded_pr151_restart_rule_training": True,
        "reproductive_direct_rows_audited": len(direct),
        "recovered_wave53_current_min3_eligible_rules": len(baseline_keys),
        "current_threshold_proxy_candidates": int(
            audit["current_threshold_candidate"].astype(str).str.casefold().eq("true").sum()
        ),
        "new_rule_candidates_vs_recovered_wave53": len(candidate),
        "potential_new_unresolved_species_axis_upper_bound": int(
            candidate["potential_new_low_species_axis_upper_bound"].sum()
        ) if len(candidate) else 0,
        "formal_gain": 0,
        "promotion_allowed": False,
    }
    (args.output / "reproductive_low_capacity_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
