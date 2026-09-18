"""Build the frozen trait-state table for prospective H4 replication.

This only recodes the already frozen v13 high/medium species-direct ledger using the
unchanged parent H5 trait mappings. It reads no pollen-limitation outcomes.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import typer
import yaml

from island_v2.chapter1_h5_glopl_floral_architecture_moderation import (
    build_trait_assignments as build_architecture_assignments,
)
from island_v2.chapter1_h5_glopl_reproductive_assurance_moderation import (
    build_trait_assignments as build_reproductive_assignments,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

RA_CONFIG = Path("config/chapter1_h5_glopl_reproductive_assurance_moderation_v1.yml")
ARCH_CONFIG = Path("config/chapter1_h5_glopl_floral_architecture_moderation_v1.yml")
TARGET_TRAITS = (
    "autonomous_selfing",
    "generalized_form",
    "actinomorphic_symmetry",
    "shallow_open_tube",
)


def _load_yaml(path: Path) -> dict:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise typer.BadParameter(f"invalid config: {path}")
    return value


def build_state_table(ledger: pd.DataFrame) -> pd.DataFrame:
    ra = build_reproductive_assignments(ledger, _load_yaml(RA_CONFIG))
    architecture = build_architecture_assignments(ledger, _load_yaml(ARCH_CONFIG))
    long = pd.concat([ra, architecture], ignore_index=True)
    long = long.loc[long["trait"].isin(TARGET_TRAITS)].copy()
    if long.empty:
        return pd.DataFrame(columns=["accepted_species", *TARGET_TRAITS])

    conflicts = long.groupby(["species_key", "trait"])["trait_state"].nunique()
    if bool((conflicts > 1).any()):
        raise typer.BadParameter("conflicting frozen trait states after parent recoding")

    names = (
        long.sort_values(["species_key", "accepted_species"])
        .drop_duplicates("species_key")[["species_key", "accepted_species"]]
    )
    wide = (
        long.pivot_table(
            index="species_key",
            columns="trait",
            values="trait_state",
            aggfunc="first",
        )
        .reset_index()
    )
    out = names.merge(wide, on="species_key", how="left", validate="one_to_one")
    for trait in TARGET_TRAITS:
        if trait not in out:
            out[trait] = pd.NA
    return (
        out[["accepted_species", *TARGET_TRAITS]]
        .sort_values("accepted_species")
        .reset_index(drop=True)
    )


@app.command("run")
def run(
    ledger_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    output_csv: Path = typer.Option(...),
) -> None:
    ledger = pd.read_csv(ledger_csv, dtype=str).fillna("")
    states = build_state_table(ledger)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    states.to_csv(output_csv, index=False)
    typer.echo(
        f"species={len(states)} "
        + " ".join(
            f"{trait}={int(states[trait].notna().sum())}" for trait in TARGET_TRAITS
        )
    )


if __name__ == "__main__":
    app()
