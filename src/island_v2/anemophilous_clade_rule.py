"""Validate declared anemophilous clade rules as a separate, reportable tier.

The strict coverage contract forbids family inference, and that prohibition is
correct as a default: most families carry real floral variation, so a family
modal value is not evidence about any particular species.

A small number of clades are different. Poaceae, Cyperaceae, Juncaceae and
*Plantago* are anemophilous throughout, and the floral consequences of wind
pollination -- reduced, inconspicuous perianth -- are near-invariant within
them. This module tests that claim instead of assuming it, using the same
machinery the genus Validated Low track already uses (support, dominance,
leave-one-out), with three additional constraints:

1. **Declared clades only.** A clade is eligible only if it is named in the
   configuration. There is no general family fallback, so a passing rule can
   never widen into inference over families nobody vetted.
2. **Stricter gates than genus.** A family is a broader clade than a genus, so
   it must clear a higher support floor and a higher dominance floor to say the
   same thing.
3. **Both leave-one-out passes.** Species LOO and source-lineage LOO must both
   hold, so a rule cannot rest on one over-sampled species or one prolific
   source.

Output is tagged `family_validated` and is never merged into the strict
numerator. Whether to admit the tier at all is a submission-contract decision
made by a human; this module only reports what would pass.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

TIER = "family_validated"
DIRECT_SCOPES = {"species_direct", "synonym_direct"}

REQUIRED_COLUMNS = {"accepted_species", "trait_name", "normalized_value", "evidence_scope"}


@dataclass(frozen=True)
class CladeThresholds:
    """Gates a declared clade rule must clear. Deliberately stricter than genus."""

    min_species: int = 20
    min_dominance: float = 0.95
    min_species_loo_accuracy: float = 0.95
    min_lineage_loo_accuracy: float = 0.95
    min_source_lineages: int = 3


DEFAULT_THRESHOLDS = CladeThresholds()

RULE_COLUMNS = [
    "clade_rank",
    "clade",
    "trait_name",
    "inferred_value",
    "n_direct_species",
    "dominant_species",
    "counterexample_species",
    "dominance",
    "n_source_lineages",
    "species_loo_n",
    "species_loo_accuracy",
    "lineage_loo_n",
    "lineage_loo_accuracy",
    "eligible",
    "ineligible_reason",
    "tier",
]


@app.callback()
def main() -> None:
    """Validate declared anemophilous clade rules and report them as a separate tier."""


def _text(value: object) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return " ".join(str(value).strip().split())


def load_declared_clades(config_path: Path) -> dict[str, set[str]]:
    """Read the declared anemophilous clades. Nothing outside this list is eligible."""
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    scope_classes = (config or {}).get("scope_classes") or {}
    families = set(scope_classes.get("anemophilous_families") or ())
    genera = set(scope_classes.get("anemophilous_genera") or ())
    if not families and not genera:
        raise typer.BadParameter(f"{config_path} declares no anemophilous clades")
    return {"family": families, "genus": genera}


def direct_rows(evidence: pd.DataFrame) -> pd.DataFrame:
    """Normalise the ledger and keep only species/synonym-direct rows.

    Kept separate from clade selection so the "already has direct evidence"
    guard can consult the whole ledger, not just the declared clades.
    """
    missing = REQUIRED_COLUMNS - set(evidence.columns)
    if missing:
        raise typer.BadParameter(f"evidence is missing columns: {sorted(missing)}")

    frame = evidence.copy()
    for column in ("accepted_species", "trait_name", "normalized_value", "evidence_scope"):
        frame[column] = frame[column].map(_text)
    if "family" in frame.columns:
        frame["family"] = frame["family"].map(_text)
    else:
        frame["family"] = ""
    # Source lineage is what makes two records independent. Without it a single
    # redistributed source could satisfy the support floor on its own.
    if "source_lineage" in frame.columns:
        frame["source_lineage"] = frame["source_lineage"].map(_text)
    else:
        frame["source_lineage"] = ""

    frame = frame.loc[frame["evidence_scope"].str.lower().isin(DIRECT_SCOPES)]
    frame = frame.loc[frame["accepted_species"].ne("") & frame["normalized_value"].ne("")]
    frame["genus"] = frame["accepted_species"].str.split().str[0]
    return frame


def _prepare(evidence: pd.DataFrame, declared: dict[str, set[str]]) -> pd.DataFrame:
    """Restrict the direct ledger to the declared clades, tagged by clade rank."""
    frame = direct_rows(evidence)

    rows = []
    for rank, column in (("family", "family"), ("genus", "genus")):
        allowed = declared.get(rank) or set()
        if not allowed:
            continue
        subset = frame.loc[frame[column].isin(allowed)].copy()
        subset["clade_rank"] = rank
        subset["clade"] = subset[column]
        rows.append(subset)

    if not rows:
        return frame.iloc[0:0].assign(clade_rank="", clade="")
    return pd.concat(rows, ignore_index=True)


def _unambiguous_species_values(group: pd.DataFrame) -> pd.Series:
    """One value per species, dropping species whose own evidence disagrees."""
    per_species = group.groupby("accepted_species")["normalized_value"].agg(
        lambda values: sorted(set(values))
    )
    return per_species.loc[per_species.map(len).eq(1)].map(lambda values: values[0])


def _loo_accuracy(labels: pd.Series) -> tuple[int, float]:
    """Leave-one-out: does the modal value survive dropping each holdout in turn?"""
    outcomes: list[bool] = []
    for key, truth in labels.items():
        training = labels.drop(index=key)
        if training.empty:
            continue
        outcomes.append(str(training.value_counts().index[0]) == str(truth))
    if not outcomes:
        return 0, 0.0
    return len(outcomes), sum(outcomes) / len(outcomes)


def build_clade_rules(
    evidence: pd.DataFrame,
    declared: dict[str, set[str]],
    thresholds: CladeThresholds = DEFAULT_THRESHOLDS,
) -> pd.DataFrame:
    """Score every declared `clade x trait` and mark which rules would pass."""
    frame = _prepare(evidence, declared)
    if frame.empty:
        return pd.DataFrame(columns=RULE_COLUMNS)

    rows: list[dict[str, Any]] = []
    for (clade_rank, clade, trait), group in frame.groupby(
        ["clade_rank", "clade", "trait_name"], sort=True
    ):
        labels = _unambiguous_species_values(group)
        n_species = len(labels)
        if n_species == 0:
            continue

        counts = labels.value_counts()
        dominant_value = str(counts.index[0])
        dominant_species = int(counts.iloc[0])
        dominance = dominant_species / n_species
        counterexamples = n_species - dominant_species

        lineages = {lineage for lineage in group["source_lineage"] if lineage}
        n_lineages = len(lineages)

        species_loo_n, species_loo_accuracy = _loo_accuracy(labels)

        # Lineage LOO collapses each lineage to its own modal value first, so a
        # single prolific source counts once rather than once per record.
        lineage_labels = (
            group.loc[group["source_lineage"].ne("")]
            .groupby("source_lineage")["normalized_value"]
            .agg(lambda values: str(pd.Series(list(values)).value_counts().index[0]))
        )
        lineage_loo_n, lineage_loo_accuracy = _loo_accuracy(lineage_labels)

        reasons: list[str] = []
        if n_species < thresholds.min_species:
            reasons.append(f"support {n_species} < {thresholds.min_species}")
        if dominance < thresholds.min_dominance:
            reasons.append(f"dominance {dominance:.4f} < {thresholds.min_dominance}")
        if n_lineages < thresholds.min_source_lineages:
            reasons.append(f"lineages {n_lineages} < {thresholds.min_source_lineages}")
        if species_loo_accuracy < thresholds.min_species_loo_accuracy:
            reasons.append(
                f"species LOO {species_loo_accuracy:.4f} < {thresholds.min_species_loo_accuracy}"
            )
        if lineage_loo_accuracy < thresholds.min_lineage_loo_accuracy:
            reasons.append(
                f"lineage LOO {lineage_loo_accuracy:.4f} < {thresholds.min_lineage_loo_accuracy}"
            )

        rows.append(
            {
                "clade_rank": clade_rank,
                "clade": clade,
                "trait_name": trait,
                "inferred_value": dominant_value,
                "n_direct_species": n_species,
                "dominant_species": dominant_species,
                "counterexample_species": counterexamples,
                "dominance": dominance,
                "n_source_lineages": n_lineages,
                "species_loo_n": species_loo_n,
                "species_loo_accuracy": species_loo_accuracy,
                "lineage_loo_n": lineage_loo_n,
                "lineage_loo_accuracy": lineage_loo_accuracy,
                "eligible": not reasons,
                "ineligible_reason": "; ".join(reasons),
                "tier": TIER,
            }
        )

    return pd.DataFrame(rows, columns=RULE_COLUMNS)


def apply_clade_rules(
    rules: pd.DataFrame, unresolved: pd.DataFrame, evidence: pd.DataFrame
) -> pd.DataFrame:
    """Emit tier rows for unresolved species covered by an eligible clade rule.

    A species that already carries direct evidence for the trait is never
    overwritten, and only eligible rules are applied.
    """
    columns = [
        "accepted_species",
        "trait_name",
        "inferred_value",
        "clade_rank",
        "clade",
        "dominance",
        "tier",
    ]
    if rules.empty or unresolved.empty:
        return pd.DataFrame(columns=columns)

    eligible = rules.loc[rules["eligible"].astype(bool)]
    if eligible.empty:
        return pd.DataFrame(columns=columns)

    targets = unresolved.copy()
    for column in ("accepted_species", "trait_name"):
        if column not in targets.columns:
            raise typer.BadParameter(f"unresolved table is missing column: {column}")
        targets[column] = targets[column].map(_text)
    if "family" in targets.columns:
        targets["family"] = targets["family"].map(_text)
    else:
        targets["family"] = ""
    targets["genus"] = targets["accepted_species"].str.split().str[0]

    direct = direct_rows(evidence)
    already = set()
    if not direct.empty:
        already = set(zip(direct["accepted_species"], direct["trait_name"], strict=True))

    rows: list[dict[str, Any]] = []
    for _, rule in eligible.iterrows():
        column = "family" if rule["clade_rank"] == "family" else "genus"
        subset = targets.loc[
            targets[column].eq(rule["clade"]) & targets["trait_name"].eq(rule["trait_name"])
        ]
        for _, target in subset.iterrows():
            key = (target["accepted_species"], target["trait_name"])
            if key in already:
                continue
            rows.append(
                {
                    "accepted_species": target["accepted_species"],
                    "trait_name": target["trait_name"],
                    "inferred_value": rule["inferred_value"],
                    "clade_rank": rule["clade_rank"],
                    "clade": rule["clade"],
                    "dominance": rule["dominance"],
                    "tier": TIER,
                }
            )

    return pd.DataFrame(rows, columns=columns).drop_duplicates(
        ["accepted_species", "trait_name"]
    )


@app.command("build")
def build(
    evidence_csv: Path = typer.Option(..., "--evidence-csv", exists=True),
    output_dir: Path = typer.Option(..., "--output-dir"),
    config_path: Path = typer.Option(
        Path("config/island_weighted_acquisition.yml"), "--config", exists=True
    ),
    unresolved_csv: Path | None = typer.Option(None, "--unresolved-csv"),
) -> None:
    """Score declared clade rules and, optionally, apply the eligible ones."""
    declared = load_declared_clades(config_path)
    evidence = pd.read_csv(evidence_csv)
    rules = build_clade_rules(evidence, declared)

    output_dir.mkdir(parents=True, exist_ok=True)
    rules.to_csv(output_dir / "anemophilous_clade_rules.csv", index=False)

    summary: dict[str, Any] = {
        "version": "1.0",
        "tier": TIER,
        "n_declared_families": len(declared["family"]),
        "n_declared_genera": len(declared["genus"]),
        "n_rules_scored": int(len(rules)),
        "n_rules_eligible": int(rules["eligible"].sum()) if len(rules) else 0,
        "thresholds": {
            "min_species": DEFAULT_THRESHOLDS.min_species,
            "min_dominance": DEFAULT_THRESHOLDS.min_dominance,
            "min_source_lineages": DEFAULT_THRESHOLDS.min_source_lineages,
            "min_species_loo_accuracy": DEFAULT_THRESHOLDS.min_species_loo_accuracy,
            "min_lineage_loo_accuracy": DEFAULT_THRESHOLDS.min_lineage_loo_accuracy,
        },
        "interpretation": (
            f"Rows are tier '{TIER}' and are never part of the strict numerator. "
            "Only clades declared in the configuration are eligible; there is no "
            "general family fallback. Admitting this tier is a submission-contract "
            "decision for a human reviewer."
        ),
    }

    if unresolved_csv is not None:
        unresolved = pd.read_csv(unresolved_csv)
        applied = apply_clade_rules(rules, unresolved, evidence)
        applied.to_csv(output_dir / "anemophilous_clade_tier_rows.csv.gz", index=False)
        summary["n_tier_rows_emitted"] = int(len(applied))

    (output_dir / "anemophilous_clade_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    app()
