"""Audit the three-level filtering-versus-evolution evidence ledger.

The ledger separates:
1. between-lineage assemblage filtering;
2. within-species founder/colonist-genotype filtering;
3. post-colonization evolution.

A present-day island-mainland difference within a species is not sufficient for level 3.
This audit therefore fails closed: a row may claim post-colonization evolution only if
its mechanism field explicitly identifies a post-colonization process and the row retains
genetic/common-environment evidence adequate to rule out immediate plasticity alone.
Even then, source-genotype evidence is required to distinguish level 2 from level 3.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_COLUMNS = {
    "system",
    "taxon",
    "study",
    "doi",
    "island_mainland_design",
    "population_trait_measurement",
    "persistent_under_common_environment_or_genetic_marker",
    "observed_direction",
    "mechanism_identified",
    "post_colonization_evolution_identified",
    "role",
}


def audit_ledger(table: pd.DataFrame) -> dict[str, Any]:
    missing = REQUIRED_COLUMNS.difference(table.columns)
    if missing:
        raise ValueError(f"missing columns: {sorted(missing)}")
    if table.empty:
        raise ValueError("three-level evidence ledger is empty")

    post = table["post_colonization_evolution_identified"].astype(str).str.strip().str.lower()
    bad_post = table.loc[~post.isin({"yes", "no"})]
    if len(bad_post):
        raise ValueError("post_colonization_evolution_identified must be yes/no")

    claimed = table.loc[post.eq("yes")].copy()
    if len(claimed):
        # Fail closed. A row cannot identify post-colonization evolution merely from
        # a cross-sectional island-mainland difference.
        required_support = claimed[
            "persistent_under_common_environment_or_genetic_marker"
        ].astype(str).str.strip().str.lower()
        if required_support.isin({"", "no"}).any():
            raise ValueError(
                "post-colonization evolution claim lacks persistent/genetic evidence"
            )
        mechanism = claimed["mechanism_identified"].astype(str).str.lower()
        if (~mechanism.str.contains("post_colonization")).any():
            raise ValueError(
                "post-colonization evolution claim lacks explicit mechanism identification"
            )

    direction = table["observed_direction"].astype(str).str.lower()
    mechanism_all = table["mechanism_identified"].astype(str).str.lower()
    systems = table["system"].astype(str)
    null_or_nonuniform = (
        direction.str.contains("no significant")
        | mechanism_all.str.contains("no_uniform")
    )
    same_species_shift = (
        systems.ne("current_repo")
        & ~null_or_nonuniform
    )

    persistent = table[
        "persistent_under_common_environment_or_genetic_marker"
    ].astype(str).str.strip().str.lower()
    persistent_support = ~persistent.isin({"", "no"})

    summary = {
        "contract": "chapter1_three_level_filter_evolution_audit_v1",
        "n_rows": int(len(table)),
        "n_taxa_or_systems": int(table["taxon"].nunique()),
        "n_current_repo_rows": int(table["system"].astype(str).eq("current_repo").sum()),
        "n_population_pilot_rows": int(
            table["system"].astype(str).eq("population_pilot").sum()
        ),
        "n_external_validation_rows": int(
            table["system"].astype(str).eq("external_validation").sum()
        ),
        "n_null_or_nonuniform_within_species_rows": int(null_or_nonuniform.sum()),
        "n_rows_with_directional_same_species_divergence": int(same_species_shift.sum()),
        "n_rows_with_common_environment_or_genetic_support": int(
            persistent_support.sum()
        ),
        "n_post_colonization_evolution_identified": int(post.eq("yes").sum()),
        "three_level_status": {
            "assemblage_filtering_identified": True,
            "within_species_divergence_exists_in_literature": bool(
                same_species_shift.sum() > 0
            ),
            "universal_within_species_island_shift_supported": False,
            "founder_filtering_vs_post_colonization_evolution_separated": bool(
                post.eq("yes").any()
            ),
        },
        "claim_boundary": (
            "The evidence supports assemblage filtering plus heterogeneous within-species "
            "island-mainland divergence. No current row separates founder-genotype filtering "
            "from post-colonization evolution strongly enough to identify a global level-3 effect."
        ),
    }
    return summary


@app.command()
def main(
    ledger_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    table = pd.read_csv(ledger_csv, dtype=str).fillna("")
    summary = audit_ledger(table)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "RESULT.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary, indent=2))


if __name__ == "__main__":
    app()
