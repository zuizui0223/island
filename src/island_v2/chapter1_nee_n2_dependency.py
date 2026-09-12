"""Prospective validation for N2 lineage functional dependency.

This module validates independently measured plant-lineage x pollination-channel
functional dependency before any N2 genus-entry outcome is inspected. It also
implements the predeclared support gate once source-specific genus dependency and
N1 channel states exist. It does not estimate N2 effects.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

LEDGER_REQUIRED = {
    "plant_species",
    "accepted_genus",
    "channel_id",
    "study_id",
    "evidence_design",
    "evidence_tier",
    "dependency_estimate",
    "dependency_lower",
    "dependency_upper",
    "source_citation",
    "source_url",
    "review_status",
    "evidence_origin",
}

GENUS_DEPENDENCY_REQUIRED = {
    "source_region_id",
    "accepted_genus",
    "channel_id",
    "dependency_mean",
    "dependency_sd",
    "n_species_evidenced",
    "evidence_tiers",
    "primary_dependency_evaluable",
}

CHANNEL_STATE_REQUIRED = {"island_id", "channel_id", "channel_state"}

PROHIBITED_ORIGIN_TOKENS = {
    "flower_color",
    "flower_colour",
    "flower_shape",
    "floral_architecture",
    "pollination_syndrome",
    "pollination_guild",
    "pollination_notes",
    "island_trait",
    "island_genus_entry",
    "channel_retention",
    "channel_disruption",
    "self_compatibility",
    "mating_system",
    "autonomous_selfing",
}


def load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("N2 dependency config must be a mapping")
    if payload.get("contract") != "chapter1_nee_n2_dependency_v1":
        raise typer.BadParameter("unexpected N2 dependency contract")
    return payload


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _require(table: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = columns.difference(table.columns)
    if missing:
        raise ValueError(f"{label} missing columns: {sorted(missing)}")


def _bool_series(series: pd.Series) -> pd.Series:
    text = _text(series).str.lower()
    mapping = {
        "true": True,
        "1": True,
        "yes": True,
        "false": False,
        "0": False,
        "no": False,
    }
    invalid = sorted(set(text).difference(mapping))
    if invalid:
        raise ValueError(f"invalid boolean values: {invalid}")
    return text.map(mapping).astype(bool)


def _tier_set(value: object) -> set[str]:
    return {token.strip() for token in str(value).split("|") if token.strip()}


def validate_dependency_ledger(
    ledger: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Validate direct dependency evidence and return primary-eligible rows only."""
    _require(ledger, LEDGER_REQUIRED, "N2 dependency ledger")
    work = ledger[list(LEDGER_REQUIRED)].copy()
    text_columns = LEDGER_REQUIRED.difference(
        {"dependency_estimate", "dependency_lower", "dependency_upper"}
    )
    for column in text_columns:
        work[column] = _text(work[column])

    key = ["plant_species", "channel_id", "study_id"]
    if work[key].duplicated().any():
        raise ValueError("duplicate study x plant_species x channel_id dependency evidence")

    channels = set(config["channels"])
    invalid_channels = sorted(set(work["channel_id"]).difference(channels))
    if invalid_channels:
        raise ValueError(f"unregistered N2 channels: {invalid_channels}")

    allowed_tiers = set(config["dependency_ledger_schema"]["accepted_primary_tiers"]) | set(
        config["dependency_ledger_schema"]["sensitivity_only_tiers"]
    )
    invalid_tiers = sorted(set(work["evidence_tier"]).difference(allowed_tiers))
    if invalid_tiers:
        raise ValueError(f"invalid dependency evidence tiers: {invalid_tiers}")

    allowed_review = set(config["dependency_ledger_schema"]["allowed_review_status"])
    invalid_review = sorted(set(work["review_status"]).difference(allowed_review))
    if invalid_review:
        raise ValueError(f"invalid dependency review status: {invalid_review}")

    accepted_designs = set(config["primary_evidence"]["accepted_designs"])
    direct_tiers = set(config["dependency_ledger_schema"]["accepted_primary_tiers"])
    candidate_primary = work["evidence_tier"].isin(direct_tiers)
    wrong_design = candidate_primary & ~work["evidence_design"].isin(accepted_designs)
    if wrong_design.any():
        bad = sorted(work.loc[wrong_design, "evidence_design"].unique().tolist())
        raise ValueError(f"primary dependency evidence uses unregistered designs: {bad}")

    origins = work["evidence_origin"].str.casefold()
    forbidden_rows = origins.map(
        lambda value: any(token in value for token in PROHIBITED_ORIGIN_TOKENS)
    )
    if (candidate_primary & forbidden_rows).any():
        bad = sorted(
            work.loc[candidate_primary & forbidden_rows, "evidence_origin"].unique().tolist()
        )
        raise ValueError(f"prohibited primary dependency evidence origin: {bad}")

    for column in ["dependency_estimate", "dependency_lower", "dependency_upper"]:
        work[column] = pd.to_numeric(work[column], errors="coerce")

    primary = work.loc[candidate_primary & work["review_status"].eq("accepted")].copy()
    if not primary.empty:
        if primary["dependency_estimate"].isna().any():
            raise ValueError("accepted D1/D2 evidence requires a quantitative dependency estimate")
        lo = float(config["dependency_ledger_schema"]["dependency_estimate_range"]["min"])
        hi = float(config["dependency_ledger_schema"]["dependency_estimate_range"]["max"])
        for column in ["dependency_estimate", "dependency_lower", "dependency_upper"]:
            finite = primary[column].dropna()
            if ((finite < lo) | (finite > hi)).any():
                raise ValueError(f"{column} must remain within [{lo}, {hi}]")
        bounded = primary["dependency_lower"].notna() & primary["dependency_upper"].notna()
        invalid_bounds = bounded & (
            (primary["dependency_lower"] > primary["dependency_estimate"])
            | (primary["dependency_estimate"] > primary["dependency_upper"])
        )
        if invalid_bounds.any():
            raise ValueError("dependency interval must contain the point estimate")
        if primary["accepted_genus"].eq("").any() or primary["plant_species"].eq("").any():
            raise ValueError("accepted primary dependency evidence requires species and genus")
        if primary["source_citation"].eq("").any():
            raise ValueError("accepted primary dependency evidence requires a source citation")

    per_channel = []
    for channel in config["channels"]:
        rows = primary.loc[primary["channel_id"].eq(channel)]
        per_channel.append(
            {
                "channel_id": channel,
                "n_primary_evidence_rows": int(len(rows)),
                "n_primary_species": int(rows["plant_species"].nunique()),
                "n_candidate_genera": int(rows["accepted_genus"].nunique()),
            }
        )
    report = {
        "contract": config["contract"],
        "status": "dependency_ledger_validated_pre_N2",
        "n_input_rows": int(len(work)),
        "n_primary_accepted_rows": int(len(primary)),
        "n_D3_rows": int(work["evidence_tier"].eq("D3").sum()),
        "channels": per_channel,
        "pollination_guild_used_as_primary_dependency": False,
        "baker_rival_relabelled_as_dependency": False,
        "N2_fitted": False,
    }
    return primary.reset_index(drop=True), report


def qualify_n2_support(
    genus_dependency: pd.DataFrame,
    channel_states: pd.DataFrame,
    config: dict[str, Any],
) -> dict[str, Any]:
    """Apply only the frozen design-support gate; no genus-entry outcome is read."""
    _require(genus_dependency, GENUS_DEPENDENCY_REQUIRED, "source-specific genus dependency")
    _require(channel_states, CHANNEL_STATE_REQUIRED, "N1 channel states")

    dep = genus_dependency.copy()
    for column in ["source_region_id", "accepted_genus", "channel_id", "evidence_tiers"]:
        dep[column] = _text(dep[column])
    dep["dependency_mean"] = pd.to_numeric(dep["dependency_mean"], errors="coerce")
    dep["dependency_sd"] = pd.to_numeric(dep["dependency_sd"], errors="coerce")
    dep["n_species_evidenced"] = pd.to_numeric(dep["n_species_evidenced"], errors="coerce")
    dep["primary_dependency_evaluable"] = _bool_series(dep["primary_dependency_evaluable"])

    invalid_channels = sorted(set(dep["channel_id"]).difference(set(config["channels"])))
    if invalid_channels:
        raise ValueError(f"unregistered dependency channels: {invalid_channels}")

    direct_tiers = set(config["dependency_ledger_schema"]["accepted_primary_tiers"])
    tier_sets = dep["evidence_tiers"].map(_tier_set)
    tier_ok = tier_sets.map(lambda tiers: bool(tiers) and tiers.issubset(direct_tiers))
    claimed_evaluable_with_D3 = dep["primary_dependency_evaluable"] & ~tier_ok
    if claimed_evaluable_with_D3.any():
        bad = sorted(dep.loc[claimed_evaluable_with_D3, "evidence_tiers"].unique().tolist())
        raise ValueError(f"D3/non-primary provenance cannot be primary dependency-evaluable: {bad}")

    single_species = dep["n_species_evidenced"].eq(1)
    exact_single_species = single_species & dep["dependency_sd"].fillna(0).le(0)
    if (dep["primary_dependency_evaluable"] & exact_single_species).any():
        raise ValueError("single-species genus dependency must retain non-zero uncertainty")

    evaluable = dep.loc[
        dep["primary_dependency_evaluable"]
        & tier_ok
        & dep["dependency_mean"].between(0.0, 1.0)
        & dep["dependency_sd"].ge(0.0)
        & dep["dependency_sd"].notna()
        & dep["n_species_evidenced"].ge(1)
        & dep["accepted_genus"].ne("")
        & dep["source_region_id"].ne("")
    ].copy()

    states = channel_states.copy()
    for column in CHANNEL_STATE_REQUIRED:
        states[column] = _text(states[column])
    invalid_states = sorted(
        set(states["channel_state"]).difference(
            {"retained", "disrupted", "unresolved", "structurally_absent"}
        )
    )
    if invalid_states:
        raise ValueError(f"invalid channel states: {invalid_states}")

    gate = config["support_qualification_before_N2_fit"]
    min_genera = int(gate["minimum_dependency_resolved_genera_per_eligible_channel"])
    min_retained = int(gate["minimum_islands_with_retained_state_per_eligible_channel"])
    min_disrupted = int(gate["minimum_islands_with_disrupted_state_per_eligible_channel"])

    channel_rows: list[dict[str, Any]] = []
    eligible_channels: list[str] = []
    for channel in config["channels"]:
        d = evaluable.loc[evaluable["channel_id"].eq(channel)]
        s = states.loc[states["channel_id"].eq(channel)]
        n_genera = int(d["accepted_genus"].nunique())
        n_retained = int(s.loc[s["channel_state"].eq("retained"), "island_id"].nunique())
        n_disrupted = int(s.loc[s["channel_state"].eq("disrupted"), "island_id"].nunique())
        eligible = (
            n_genera >= min_genera
            and n_retained >= min_retained
            and n_disrupted >= min_disrupted
        )
        if eligible:
            eligible_channels.append(channel)
        channel_rows.append(
            {
                "channel_id": channel,
                "n_dependency_resolved_genera": n_genera,
                "n_retained_islands": n_retained,
                "n_disrupted_islands": n_disrupted,
                "eligible_for_primary_N2": eligible,
            }
        )

    primary_support = evaluable.loc[evaluable["channel_id"].isin(eligible_channels)]
    total_genera = int(primary_support["accepted_genus"].nunique())
    passed = (
        len(eligible_channels) >= int(gate["minimum_confirmatory_channels"])
        and total_genera >= int(gate["minimum_unique_dependency_resolved_genera_total"])
    )
    return {
        "contract": config["contract"],
        "status": "N2_support_qualification_only",
        "passed": bool(passed),
        "eligible_channels": eligible_channels,
        "n_eligible_channels": int(len(eligible_channels)),
        "n_unique_dependency_resolved_genera_total": total_genera,
        "channels": channel_rows,
        "genus_entry_outcome_read": False,
        "D3_can_satisfy_primary_support": False,
        "failure_action": gate["qualification_failure_action"],
        "N2_fitted": False,
    }


@app.command("validate-ledger")
def validate_ledger_command(
    ledger_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/chapter1_nee_n2_dependency.yml")),
) -> None:
    config = load_config(config_path)
    ledger = pd.read_csv(ledger_csv, dtype=str).fillna("")
    primary, report = validate_dependency_ledger(ledger, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    primary.to_csv(output_dir / "n2_primary_dependency_evidence.csv", index=False)
    (output_dir / "n2_dependency_ledger_receipt.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(report))


@app.command("qualify-support")
def qualify_support_command(
    genus_dependency_csv: Path = typer.Option(..., exists=True),
    channel_states_csv: Path = typer.Option(..., exists=True),
    output_json: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/chapter1_nee_n2_dependency.yml")),
) -> None:
    config = load_config(config_path)
    dependency = pd.read_csv(genus_dependency_csv, dtype=str).fillna("")
    states = pd.read_csv(channel_states_csv, dtype=str).fillna("")
    report = qualify_n2_support(dependency, states, config)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    typer.echo(json.dumps(report))


if __name__ == "__main__":
    app()
