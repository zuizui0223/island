"""Finalize reviewed open-Web evidence with trait-specific Low rebuilding.

Audit errors remain in the precision denominator. Genus rules are applied via
``genus x axis x trait_name``. Reward and pollen-vector evidence stay as
separate outcomes and never fill the floral-structure axis.
"""

# Typer command defaults intentionally use Option metadata.
# ruff: noqa: B008

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd
import typer

from island_v2.open_web_evidence import promote_reviewed
from island_v2.open_web_pilot import (
    DIRECT_QUALITY,
    _bool,
    _production_approved_traits,
    _read_baseline,
    _text,
    _validate_denominator,
)

CONTRACT = "open_web_review_finalization_v2"
app = typer.Typer(add_completion=False, no_args_is_help=True)

STRICT_AXIS_TRAITS = {
    "flower_colour": ("flower_primary_color",),
    "floral_structural_complexity": (
        "floral_form",
        "floral_symmetry",
        "tube_depth_class",
        "flower_size_class",
        "inflorescence_display",
    ),
    "reproductive_assurance": (
        "self_incompatibility",
        "autonomous_selfing_capacity",
        "mating_system",
        "cleistogamy",
    ),
}
STRICT_TRAIT_AXIS = {
    trait: axis for axis, traits in STRICT_AXIS_TRAITS.items() for trait in traits
}
DOMINANCE = {
    "flower_colour": 0.90,
    "floral_structural_complexity": 0.80,
    "reproductive_assurance": 0.95,
}
MASKED = {
    "flower_colour": 0.80,
    "floral_structural_complexity": 0.75,
    "reproductive_assurance": 0.85,
}
RULE_COLUMNS = [
    "genus",
    "axis",
    "trait_name",
    "inferred_value",
    "n_direct_species",
    "dominant_species",
    "counterexample_species",
    "dominance",
    "required_dominance",
    "masked_n",
    "masked_correct",
    "masked_accuracy",
    "required_masked_accuracy",
    "eligible",
]
LOW_COLUMNS = [
    "accepted_species",
    "genus",
    "axis",
    "trait_name",
    *[column for column in RULE_COLUMNS if column not in {"genus", "axis", "trait_name"}],
    "normalized_value",
    "quality",
    "evidence_scope",
    "inference_method",
    "family_inference_used",
    "global_fallback_used",
]


def _now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _empty_rules() -> pd.DataFrame:
    return pd.DataFrame(columns=RULE_COLUMNS)


def _empty_low() -> pd.DataFrame:
    return pd.DataFrame(columns=LOW_COLUMNS)


def _winning_direct(
    evidence: pd.DataFrame,
    promoted: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    columns = [
        "accepted_species",
        "trait_name",
        "normalized_value",
        "quality",
        "source_lineage",
    ]
    direct = evidence.loc[
        evidence["quality"].isin(DIRECT_QUALITY)
        & evidence["evidence_scope"].isin({"species_direct", "synonym_direct"})
        & evidence["trait_name"].isin(STRICT_TRAIT_AXIS),
        columns,
    ].copy()
    if not promoted.empty:
        added = promoted[
            [
                "accepted_species",
                "trait_name",
                "normalized_value",
                "confidence",
                "source_lineage",
            ]
        ].rename(columns={"confidence": "quality"})
        direct = pd.concat(
            [direct, added.loc[added["trait_name"].isin(STRICT_TRAIT_AXIS)]],
            ignore_index=True,
        )
    if direct.empty:
        return pd.DataFrame(columns=columns), pd.DataFrame(
            columns=["accepted_species", "trait_name", "normalized_value"]
        )

    direct["quality_rank"] = direct["quality"].map({"medium": 1, "high": 2}).fillna(0)
    best_rank = direct.groupby(["accepted_species", "trait_name"])[
        "quality_rank"
    ].transform("max")
    winning = direct.loc[direct["quality_rank"].eq(best_rank)].drop_duplicates(
        [
            "accepted_species",
            "trait_name",
            "normalized_value",
            "source_lineage",
        ]
    )
    values = (
        winning.groupby(["accepted_species", "trait_name"])["normalized_value"]
        .agg(lambda items: sorted({_text(item) for item in items if _text(item)}))
        .reset_index()
    )
    conflicts = values.loc[values["normalized_value"].map(len).gt(1)].copy()
    unambiguous_keys = values.loc[
        values["normalized_value"].map(len).eq(1),
        ["accepted_species", "trait_name"],
    ]
    unambiguous = winning.merge(
        unambiguous_keys,
        on=["accepted_species", "trait_name"],
        how="inner",
    )
    unambiguous = unambiguous.loc[
        ~unambiguous["normalized_value"].isin(
            {"multicolored_variable", "multistate_variable"}
        )
    ].copy()
    return unambiguous, conflicts


def _genus_rules(direct: pd.DataFrame) -> pd.DataFrame:
    if direct.empty:
        return _empty_rules()
    votes = direct.drop_duplicates(
        ["accepted_species", "trait_name", "normalized_value"]
    ).copy()
    votes["genus"] = votes["accepted_species"].str.split().str[0]
    species_values = (
        votes.groupby(["genus", "trait_name", "accepted_species"])[
            "normalized_value"
        ]
        .agg(lambda values: sorted(set(values)))
        .reset_index()
    )
    species_values = species_values.loc[
        species_values["normalized_value"].map(len).eq(1)
    ].copy()
    species_values["value"] = species_values["normalized_value"].str[0]
    rows: list[dict[str, object]] = []
    for (genus, trait), group in species_values.groupby(["genus", "trait_name"]):
        axis = STRICT_TRAIT_AXIS[trait]
        counts = group["value"].value_counts()
        n_species = int(group["accepted_species"].nunique())
        inferred_value = _text(counts.index[0])
        dominant_species = int(counts.iloc[0])
        dominance = dominant_species / n_species
        outcomes: list[bool] = []
        for index, held in group.iterrows():
            training = group.drop(index=index)
            if training.empty:
                continue
            prediction = _text(training["value"].value_counts().index[0])
            outcomes.append(prediction == _text(held["value"]))
        masked_accuracy = sum(outcomes) / len(outcomes) if outcomes else 0.0
        rows.append(
            {
                "genus": genus,
                "axis": axis,
                "trait_name": trait,
                "inferred_value": inferred_value,
                "n_direct_species": n_species,
                "dominant_species": dominant_species,
                "counterexample_species": n_species - dominant_species,
                "dominance": dominance,
                "required_dominance": DOMINANCE[axis],
                "masked_n": len(outcomes),
                "masked_correct": int(sum(outcomes)),
                "masked_accuracy": masked_accuracy,
                "required_masked_accuracy": MASKED[axis],
                "eligible": bool(
                    n_species >= 3
                    and dominance >= DOMINANCE[axis]
                    and masked_accuracy >= MASKED[axis]
                ),
            }
        )
    return pd.DataFrame(rows, columns=RULE_COLUMNS)


def rebuild_validated_low(
    *,
    coverage: pd.DataFrame,
    evidence: pd.DataFrame,
    promoted: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    direct, conflicts = _winning_direct(evidence, promoted)
    rules = _genus_rules(direct)
    direct_axes = {
        (species, STRICT_TRAIT_AXIS[trait])
        for species, trait in zip(direct["accepted_species"], direct["trait_name"])
    }
    unresolved_axes = coverage.loc[
        [
            (species, axis) not in direct_axes
            for species, axis in zip(
                coverage["accepted_species"], coverage["axis"]
            )
        ],
        ["accepted_species", "axis"],
    ]
    task_rows = [
        {
            "accepted_species": _text(row.accepted_species),
            "genus": _text(row.accepted_species).split()[0],
            "axis": _text(row.axis),
            "trait_name": trait,
        }
        for row in unresolved_axes.itertuples(index=False)
        for trait in STRICT_AXIS_TRAITS[_text(row.axis)]
    ]
    tasks = pd.DataFrame(
        task_rows,
        columns=["accepted_species", "genus", "axis", "trait_name"],
    )
    eligible = rules.loc[rules["eligible"]].copy() if not rules.empty else _empty_rules()
    if tasks.empty or eligible.empty:
        low = _empty_low()
    else:
        low = tasks.merge(
            eligible,
            on=["genus", "axis", "trait_name"],
            how="inner",
        )
        low["normalized_value"] = low["inferred_value"]
        low["quality"] = "low"
        low["evidence_scope"] = "genus_consensus"
        low["inference_method"] = "validated_genus_consensus"
        low["family_inference_used"] = False
        low["global_fallback_used"] = False
        low = low.drop_duplicates(
            ["accepted_species", "trait_name", "normalized_value"]
        ).reindex(columns=LOW_COLUMNS)

    old_low_axes = set(
        zip(
            coverage.loc[coverage["after_quality"].eq("low"), "accepted_species"],
            coverage.loc[coverage["after_quality"].eq("low"), "axis"],
        )
    )
    rebuilt_low_axes = set(zip(low["accepted_species"], low["axis"]))
    report: dict[str, object] = {
        "contract": "trait_specific_all_direct_validated_low_v2",
        "quality_precedence": ["high", "medium"],
        "trait_specific_join": True,
        "family_inference": False,
        "global_fallback": False,
        "direct_species_trait_conflicts_excluded": len(conflicts),
        "eligible_genus_trait_rules": int(rules["eligible"].sum()) if not rules.empty else 0,
        "rebuilt_validated_low_species_trait": int(
            low.drop_duplicates(["accepted_species", "trait_name"]).shape[0]
        ),
        "rebuilt_validated_low_species_axis": len(rebuilt_low_axes),
        "validated_low_added": len(rebuilt_low_axes - old_low_axes),
        "validated_low_invalidated": len(old_low_axes - rebuilt_low_axes),
        "validated_low_retained": len(old_low_axes & rebuilt_low_axes),
    }
    return rules, low, conflicts, report


def apply_review(
    *,
    baseline_dir: Path,
    candidates_csv: Path,
    audit_csv: Path,
    output_dir: Path,
) -> dict[str, object]:
    summary, coverage, evidence, _ = _read_baseline(baseline_dir)
    _validate_denominator(summary, coverage)
    candidates = pd.read_csv(candidates_csv, dtype=str).fillna("")
    audits = pd.read_csv(audit_csv, dtype=str).fillna("")
    required = {
        "candidate_id",
        "decision",
        "species_identity_correct",
        "value_correct",
        "provenance_complete",
        "cultivar_contamination",
    }
    missing = required.difference(audits.columns)
    if missing:
        raise ValueError(f"audit is missing columns: {sorted(missing)}")
    reviewed = audits.loc[
        audits["decision"].str.casefold().isin({"accept", "reject"})
    ].copy()
    good_mask = (
        reviewed["decision"].str.casefold().eq("accept")
        & reviewed["species_identity_correct"].map(_bool)
        & reviewed["value_correct"].map(_bool)
        & reviewed["provenance_complete"].map(_bool)
        & ~reviewed["cultivar_contamination"].map(_bool)
    )
    accepted_audit = reviewed.loc[good_mask]
    promoted = pd.DataFrame(
        promote_reviewed(
            candidates.to_dict("records"),
            accepted_audit.assign(decision="accept").to_dict("records"),
        )
    )
    precision = len(accepted_audit) / len(reviewed) if len(reviewed) else None
    cultivar_rate = (
        float(reviewed["cultivar_contamination"].map(_bool).mean())
        if len(reviewed)
        else None
    )
    by_trait: dict[str, dict[str, object]] = {}
    for trait, group in reviewed.groupby("trait_name"):
        good = group.loc[good_mask.loc[group.index]]
        by_trait[_text(trait)] = {
            "audited": len(group),
            "accepted_correct": len(good),
            "precision": len(good) / len(group) if len(group) else None,
        }
    approved_traits = _production_approved_traits(by_trait)
    rules, low, conflicts, low_report = rebuild_validated_low(
        coverage=coverage,
        evidence=evidence,
        promoted=promoted,
    )
    baseline_direct_axes = set(
        zip(
            coverage.loc[
                coverage["after_quality"].isin(DIRECT_QUALITY), "accepted_species"
            ],
            coverage.loc[coverage["after_quality"].isin(DIRECT_QUALITY), "axis"],
        )
    )
    promoted_axes = {
        (row.accepted_species, STRICT_TRAIT_AXIS[row.trait_name])
        for row in promoted.itertuples(index=False)
        if row.trait_name in STRICT_TRAIT_AXIS
    }
    low_axes = set(zip(low["accepted_species"], low["axis"]))
    final_axes = baseline_direct_axes | promoted_axes | low_axes
    baseline_axes = set(
        zip(
            coverage.loc[coverage["after_quality"].ne(""), "accepted_species"],
            coverage.loc[coverage["after_quality"].ne(""), "axis"],
        )
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    promoted.to_csv(output_dir / "accepted_open_web_evidence.csv", index=False)
    rules.to_csv(output_dir / "rebuilt_validated_genus_rules.csv.gz", index=False)
    low.to_csv(output_dir / "rebuilt_validated_low_evidence.csv.gz", index=False)
    conflicts.to_csv(
        output_dir / "direct_species_trait_conflicts_excluded.csv.gz", index=False
    )
    report: dict[str, object] = {
        "contract": CONTRACT,
        "completed_at_utc": _now(),
        "audited_pages_or_candidates": len(reviewed),
        "accepted_correct": len(accepted_audit),
        "precision_denominator": len(reviewed),
        "precision": precision,
        "cultivar_contamination_rate": cultivar_rate,
        "by_trait": by_trait,
        "production_approved_traits": approved_traits,
        "new_direct_species_trait": (
            int(promoted.drop_duplicates(["accepted_species", "trait_name"]).shape[0])
            if not promoted.empty
            else 0
        ),
        "new_direct_species_axis": len(promoted_axes - baseline_direct_axes),
        "validated_low_rebuild": low_report,
        "strict_filled_species_axis_before": len(baseline_axes),
        "strict_filled_species_axis_after": len(final_axes),
        "strict_coverage_before": len(baseline_axes) / 318_885,
        "strict_coverage_after": len(final_axes) / 318_885,
        "strict_coverage_increment": (len(final_axes) - len(baseline_axes)) / 318_885,
        "production_approved": bool(
            len(reviewed) >= 100
            and precision is not None
            and precision >= 0.90
            and (cultivar_rate or 0.0) <= 0.02
            and approved_traits
        ),
        "family_inference": False,
        "global_fallback": False,
    }
    (output_dir / "manual_audit_gain_report_v2.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report


@app.command("apply-audit")
def apply_audit_command(
    baseline_dir: Path = typer.Option(..., exists=True, file_okay=False),
    candidates_csv: Path = typer.Option(..., exists=True),
    audit_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    report = apply_review(
        baseline_dir=baseline_dir,
        candidates_csv=candidates_csv,
        audit_csv=audit_csv,
        output_dir=output_dir,
    )
    typer.echo(json.dumps(report, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    app()
