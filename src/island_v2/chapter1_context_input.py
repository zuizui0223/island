"""Build the canonical Bombus-free Chapter 1 island trait-composition table.

The output is long-form island x floristic-status stratum x trait x category counts.
Only unambiguous direct species-level trait evidence is admitted. Pollinator guild,
Bombus occurrence, and pollination-syndrome interpretation are deliberately absent.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

from island_v2.flora_status_support import DIRECT_SCOPES, STRATA, _text, attach_floristic_status, stratum_mask

app = typer.Typer(add_completion=False, no_args_is_help=True)

PRIMARY_TRAITS = (
    "flower_primary_color",
    "floral_form",
    "self_incompatibility",
)


def collapse_direct_trait_states(
    direct_evidence: pd.DataFrame,
    *,
    traits: tuple[str, ...] = PRIMARY_TRAITS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {"accepted_species", "trait_name", "normalized_value"}
    missing = required - set(direct_evidence.columns)
    if missing:
        raise typer.BadParameter(f"direct evidence missing columns: {sorted(missing)}")

    frame = direct_evidence.copy()
    frame["accepted_species"] = _text(frame["accepted_species"])
    frame["trait_name"] = _text(frame["trait_name"])
    frame["normalized_value"] = _text(frame["normalized_value"])
    frame = frame.loc[
        frame["trait_name"].isin(traits)
        & frame["accepted_species"].ne("")
        & frame["normalized_value"].ne("")
    ].copy()
    if "evidence_scope" in frame.columns:
        frame = frame.loc[_text(frame["evidence_scope"]).str.lower().isin(DIRECT_SCOPES)].copy()
    if "resolution_status" in frame.columns:
        frame = frame.loc[_text(frame["resolution_status"]).str.lower().eq("resolved")].copy()

    rows = []
    audit = []
    for (species, trait), group in frame.groupby(["accepted_species", "trait_name"], sort=True):
        values = sorted(set(group["normalized_value"]))
        # Multi-state strings are retained only when every direct record agrees on the
        # same exact normalized state. Cross-record conflicts are excluded from primary.
        resolved = len(values) == 1
        audit.append({
            "accepted_species": species,
            "trait_name": trait,
            "n_direct_rows": int(len(group)),
            "n_distinct_values": int(len(values)),
            "resolved_for_primary": resolved,
            "values": "|".join(values),
        })
        if resolved:
            rows.append({
                "accepted_species": species,
                "trait_name": trait,
                "trait_state": values[0],
            })
    return pd.DataFrame(rows), pd.DataFrame(audit)


def build_island_trait_composition(
    island_species: pd.DataFrame,
    status_ledger: pd.DataFrame,
    direct_evidence: pd.DataFrame,
    *,
    strata: tuple[str, ...] = ("all_native", "native_nonendemic", "endemic"),
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    invalid = set(strata) - set(STRATA)
    if invalid:
        raise typer.BadParameter(f"unknown strata: {sorted(invalid)}")
    status_flora = attach_floristic_status(island_species, status_ledger)
    states, evidence_audit = collapse_direct_trait_states(direct_evidence)

    long_rows: list[pd.DataFrame] = []
    for stratum in strata:
        subset = status_flora.loc[
            stratum_mask(status_flora, stratum), ["island_id", "accepted_species"]
        ].copy()
        joined = subset.merge(states, on="accepted_species", how="inner", validate="many_to_many")
        if joined.empty:
            continue
        category_counts = joined.groupby(
            ["island_id", "trait_name", "trait_state"], as_index=False
        ).agg(successes=("accepted_species", "nunique"))
        trials = joined.groupby(["island_id", "trait_name"], as_index=False).agg(
            trials=("accepted_species", "nunique")
        )
        out = category_counts.merge(trials, on=["island_id", "trait_name"], how="left", validate="many_to_one")
        out["share"] = out["successes"] / out["trials"]
        out["stratum"] = stratum
        long_rows.append(out)

    composition = pd.concat(long_rows, ignore_index=True) if long_rows else pd.DataFrame()
    if not composition.empty:
        composition = composition.sort_values(
            ["stratum", "trait_name", "trait_state", "island_id"]
        ).reset_index(drop=True)
    return composition, status_flora, evidence_audit


@app.command("run")
def run(
    island_species_csv: Path = typer.Option(..., exists=True),
    status_ledger_csv: Path = typer.Option(..., exists=True),
    direct_evidence_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    composition, status_flora, evidence_audit = build_island_trait_composition(
        pd.read_csv(island_species_csv),
        pd.read_csv(status_ledger_csv),
        pd.read_csv(direct_evidence_csv),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    composition.to_csv(output_dir / "chapter1_trait_composition_long.csv.gz", index=False)
    status_flora.to_csv(output_dir / "chapter1_status_flora.csv.gz", index=False)
    evidence_audit.to_csv(output_dir / "chapter1_direct_trait_state_audit.csv.gz", index=False)
    manifest = {
        "contract": "chapter1_context_input_v1",
        "primary_traits": list(PRIMARY_TRAITS),
        "primary_strata": ["all_native", "native_nonendemic", "endemic"],
        "pollinator_predictors_included": False,
        "n_composition_rows": int(len(composition)),
        "n_trait_state_audit_rows": int(len(evidence_audit)),
        "n_primary_resolved_trait_states": int(evidence_audit["resolved_for_primary"].sum()) if not evidence_audit.empty else 0,
    }
    (output_dir / "chapter1_context_input_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
