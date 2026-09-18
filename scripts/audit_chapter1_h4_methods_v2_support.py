"""Outcome-blind support audit for H4 Methods v2 candidates.

This is diagnostic only. It never changes the frozen H4a/H4b support thresholds and
never reads reproductive outcomes. It reports support for:
- automatic design candidates only;
- automatic + field-unresolved candidates (reviewable upper bound);
- all core-design candidates including exclusion conflicts (diagnostic ceiling).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.chapter1_h4_prospective_temporal_replication import (
    _h4a_support,
    _h4b_support,
    prepare_frozen_trait_states,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

STATUS_TIERS = {
    "automatic_only": {"design_candidate_auto"},
    "auto_plus_field_unresolved": {
        "design_candidate_auto",
        "design_candidate_field_unresolved",
    },
    "all_core_design_diagnostic": {
        "design_candidate_auto",
        "design_candidate_field_unresolved",
        "design_candidate_conflict_unresolved",
    },
}


def _load_contract(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    expected = "chapter1_h4_prospective_temporal_replication_v1"
    if not isinstance(config, dict) or config.get("contract") != expected:
        raise typer.BadParameter("unexpected prospective H4 contract")
    return config


def expand_target_species(
    methods: pd.DataFrame,
    traits: pd.DataFrame,
    statuses: set[str],
) -> pd.DataFrame:
    required = {"doi", "pmcid", "target_species", "candidate_status"}
    if missing := required - set(methods.columns):
        raise typer.BadParameter(f"Methods v2 table missing columns: {sorted(missing)}")

    frozen = prepare_frozen_trait_states(traits, _load_contract(
        Path("config/chapter1_h4_prospective_temporal_replication_v1.yml")
    ))
    rows: list[dict[str, str]] = []
    admitted = methods.loc[methods["candidate_status"].astype(str).isin(statuses)].copy()
    for row in admitted.itertuples(index=False):
        study_key = str(row.doi or "").strip().casefold() or str(row.pmcid or "").strip()
        if not study_key:
            continue
        for species in sorted(
            {
                value.strip()
                for value in str(row.target_species or "").split("|")
                if value.strip()
            }
        ):
            rows.append(
                {
                    "study_key": study_key,
                    "accepted_species": species,
                }
            )
    base = pd.DataFrame(rows, columns=["study_key", "accepted_species"]).drop_duplicates()
    if base.empty:
        return base.merge(frozen, on="accepted_species", how="left")
    return base.merge(
        frozen,
        on="accepted_species",
        how="left",
        validate="many_to_one",
    )


def support_audit(
    methods: pd.DataFrame,
    traits: pd.DataFrame,
    contract: dict[str, Any],
) -> dict[str, Any]:
    tiers: dict[str, Any] = {}
    for name, statuses in STATUS_TIERS.items():
        matched = expand_target_species(methods, traits, statuses)
        tiers[name] = {
            "candidate_statuses": sorted(statuses),
            "n_publications": int(matched["study_key"].nunique()) if not matched.empty else 0,
            "n_target_species": int(matched["accepted_species"].nunique()) if not matched.empty else 0,
            "H4a_reproductive_assurance": _h4a_support(matched, contract)
            if not matched.empty
            else {"evaluable": False, "matched_species_total": 0, "publications_total": 0},
            "H4b_accessibility_generalization": _h4b_support(matched, contract)
            if not matched.empty
            else {"evaluable": False, "matched_species_total": 0, "publications_total": 0},
        }
    return {
        "inferential_role": "outcome_blind_support_diagnostic_not_support_lock",
        "thresholds_changed": False,
        "outcomes_read": False,
        "tiers": tiers,
    }


@app.command("run")
def run(
    methods_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    trait_states_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    contract_path: Path = typer.Option(..., exists=True, dir_okay=False),
    output_json: Path = typer.Option(...),
) -> None:
    methods = pd.read_csv(methods_csv).fillna("")
    traits = pd.read_csv(trait_states_csv)
    contract = _load_contract(contract_path)
    result = support_audit(methods, traits, contract)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    typer.echo(json.dumps(result, indent=2))


if __name__ == "__main__":
    app()
