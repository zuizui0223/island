"""Rank trait-acquisition targets by the island composition they actually unlock.

The integrated coverage ledger scores progress on a fixed
`106,295 accepted species x 3 axes` denominator, which weights every species
equally. The confirmatory models do not read species: they read island-level
composition built from the realised flora of each island. A species occupying
993 islands and a single-island endemic are one cell each in the coverage
denominator and differ by three orders of magnitude in the analysis.

This module makes that difference measurable. It reports, for the frozen flora
universe:

- `flora mass` -- total `island x species` memberships, the quantity island
  composition is actually averaged over;
- per-species acquisition value, the share of flora mass a species carries;
- island readiness under any candidate coverage set; and
- the planning frontier, comparing island-weighted targeting against the
  unweighted species-count strategy at equal species budgets.

Nothing here accepts, promotes or infers a trait value. It ranks what to go and
find, and measures what a set of already-covered species would deliver.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_MASTER_COLUMNS = {"accepted_species", "genus", "family", "n_islands"}
REQUIRED_ISLAND_COLUMNS = {"island_id", "species"}

BIOTIC = "biotic_floral_evidence_required"
ANEMOPHILOUS = "anemophilous_candidate"
NON_ANGIOSPERM = "non_angiosperm_out_of_floral_scope"


@app.callback()
def main() -> None:
    """Rank and evaluate trait acquisition by island-composition value."""


def load_config(path: Path) -> dict[str, Any]:
    """Load and minimally validate the versioned acquisition-value configuration."""
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    required = {"inputs", "readiness", "frontier_targets", "scope_classes"}
    if not isinstance(config, dict) or not required.issubset(config):
        raise typer.BadParameter(
            "config must contain inputs, readiness, frontier_targets and scope_classes"
        )
    readiness = config["readiness"]
    if not {"min_covered_species", "min_covered_fraction"}.issubset(readiness):
        raise typer.BadParameter(
            "readiness must define min_covered_species and min_covered_fraction"
        )
    return config


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def scope_class(family: str, genus: str, scope_classes: dict[str, Any]) -> str:
    """Label how floral evidence for a species must be obtained.

    The label is an acquisition-planning class, not a trait value. An
    `anemophilous_candidate` still needs reviewed source-backed evidence or a
    separately validated and separately reported clade rule.
    """
    if family in set(scope_classes.get("non_angiosperm_families") or ()):
        return NON_ANGIOSPERM
    if family in set(scope_classes.get("anemophilous_families") or ()):
        return ANEMOPHILOUS
    if genus in set(scope_classes.get("anemophilous_genera") or ()):
        return ANEMOPHILOUS
    return BIOTIC


def build_species_value(master: pd.DataFrame, scope_classes: dict[str, Any]) -> pd.DataFrame:
    """Rank every master species by the share of island flora mass it carries."""
    missing = REQUIRED_MASTER_COLUMNS - set(master.columns)
    if missing:
        raise typer.BadParameter(f"master taxa file is missing columns: {sorted(missing)}")

    frame = master.loc[:, sorted(REQUIRED_MASTER_COLUMNS)].copy()
    for column in ("accepted_species", "genus", "family"):
        frame[column] = _text(frame[column])
    frame["n_islands"] = pd.to_numeric(frame["n_islands"], errors="coerce").fillna(0).astype(int)
    frame = frame.loc[frame["accepted_species"].ne("")]

    total_mass = int(frame["n_islands"].sum())
    if total_mass <= 0:
        raise typer.BadParameter("master taxa file carries no island memberships")

    # Ties broken by name so the ranking is reproducible across runs.
    frame = frame.sort_values(
        ["n_islands", "accepted_species"], ascending=[False, True]
    ).reset_index(drop=True)
    frame["priority_rank"] = frame.index + 1
    frame["flora_mass_share"] = frame["n_islands"] / total_mass
    frame["cumulative_flora_mass_share"] = frame["flora_mass_share"].cumsum()
    frame["scope_class"] = [
        scope_class(family, genus, scope_classes)
        for family, genus in zip(frame["family"], frame["genus"], strict=True)
    ]
    return frame


def island_readiness(
    island_species: pd.DataFrame,
    covered: set[str],
    min_covered_species: int,
    min_covered_fraction: float,
) -> pd.DataFrame:
    """Score every island by how much of its realised flora the covered set reaches."""
    missing = REQUIRED_ISLAND_COLUMNS - set(island_species.columns)
    if missing:
        raise typer.BadParameter(f"island species file is missing columns: {sorted(missing)}")

    frame = island_species.loc[:, ["island_id", "species"]].copy()
    frame["island_id"] = _text(frame["island_id"])
    frame["species"] = _text(frame["species"])
    frame = frame.loc[frame["island_id"].ne("") & frame["species"].ne("")].drop_duplicates()
    frame["is_covered"] = frame["species"].isin(covered)

    readiness = (
        frame.groupby("island_id", as_index=False)
        .agg(n_flora_species=("species", "size"), n_covered_species=("is_covered", "sum"))
        .sort_values("island_id")
        .reset_index(drop=True)
    )
    readiness["n_covered_species"] = readiness["n_covered_species"].astype(int)
    readiness["covered_fraction"] = readiness["n_covered_species"] / readiness["n_flora_species"]
    readiness["meets_min_covered_species"] = readiness["n_covered_species"] >= min_covered_species
    readiness["meets_min_covered_fraction"] = readiness["covered_fraction"] >= min_covered_fraction
    readiness["analysis_ready"] = (
        readiness["meets_min_covered_species"] & readiness["meets_min_covered_fraction"]
    )
    return readiness


def _readiness_summary(readiness: pd.DataFrame, n_species: int) -> dict[str, Any]:
    return {
        "n_target_species": int(n_species),
        "n_islands": len(readiness),
        "n_islands_meeting_min_covered_species": int(readiness["meets_min_covered_species"].sum()),
        "n_islands_meeting_min_covered_fraction": int(
            readiness["meets_min_covered_fraction"].sum()
        ),
        "n_islands_analysis_ready": int(readiness["analysis_ready"].sum()),
        "median_island_covered_fraction": float(readiness["covered_fraction"].median()),
    }


def build_frontier(
    species_value: pd.DataFrame,
    island_species: pd.DataFrame,
    targets: list[int],
    min_covered_species: int,
    min_covered_fraction: float,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Compare island-weighted targeting with the unweighted strategy per budget.

    The unweighted comparison is deterministic rather than sampled: it takes the
    same budget in master name order, which is how a source-driven campaign
    reaches species -- by whatever the next source happens to contain, not by
    island footprint.
    """
    ranked = species_value["accepted_species"].tolist()
    unweighted = sorted(ranked)

    records: list[dict[str, Any]] = []
    for budget in targets:
        budget = int(budget)
        if budget <= 0:
            continue
        for strategy, order in (
            ("island_weighted", ranked),
            ("unweighted_species_count", unweighted),
        ):
            selected = set(order[:budget])
            readiness = island_readiness(
                island_species, selected, min_covered_species, min_covered_fraction
            )
            summary = _readiness_summary(readiness, len(selected))
            summary["strategy"] = strategy
            summary["flora_mass_share"] = float(
                species_value.loc[
                    species_value["accepted_species"].isin(selected), "flora_mass_share"
                ].sum()
            )
            records.append(summary)

    frontier = pd.DataFrame.from_records(records)
    frontier = frontier.loc[
        :,
        [
            "strategy",
            "n_target_species",
            "flora_mass_share",
            "median_island_covered_fraction",
            "n_islands_meeting_min_covered_species",
            "n_islands_meeting_min_covered_fraction",
            "n_islands_analysis_ready",
            "n_islands",
        ],
    ].sort_values(["n_target_species", "strategy"]).reset_index(drop=True)

    weighted = frontier.loc[frontier["strategy"] == "island_weighted"]
    plain = frontier.loc[frontier["strategy"] == "unweighted_species_count"]
    summary = {
        "version": "1.0",
        "n_master_species": len(species_value),
        "total_flora_mass": int(species_value["n_islands"].sum()),
        "n_islands": int(frontier["n_islands"].max()) if len(frontier) else 0,
        "readiness_gate": {
            "min_covered_species": int(min_covered_species),
            "min_covered_fraction": float(min_covered_fraction),
        },
        "best_island_weighted_ready": (
            int(weighted["n_islands_analysis_ready"].max()) if len(weighted) else 0
        ),
        "best_unweighted_ready": (
            int(plain["n_islands_analysis_ready"].max()) if len(plain) else 0
        ),
        "interpretation": (
            "Analysis-ready counts are what a covered species set would deliver to "
            "island composition models. They are not trait values, not accepted "
            "evidence, and not biological absences."
        ),
    }
    return frontier, summary


def scope_breakdown(species_value: pd.DataFrame, budget: int) -> pd.DataFrame:
    """Split a target budget by how its floral evidence has to be obtained."""
    head = species_value.loc[species_value["priority_rank"] <= int(budget)]
    breakdown = (
        head.groupby("scope_class", as_index=False)
        .agg(n_species=("accepted_species", "size"), flora_mass_share=("flora_mass_share", "sum"))
        .sort_values("flora_mass_share", ascending=False)
        .reset_index(drop=True)
    )
    breakdown.insert(0, "target_budget", int(budget))
    return breakdown


DIRECT_SCOPES = {"species_direct", "synonym_direct"}
VALIDATED_LOW_SCOPES = {"validated_low", "genus_validated_low", "genus_consensus"}


def _axis_of(trait_name: str, axes: dict[str, Any]) -> str:
    for axis, traits in axes.items():
        if trait_name in set(traits):
            return axis
    return ""


def covered_species_by_axis(
    evidence: pd.DataFrame, axes: dict[str, Any], include_validated_low: bool
) -> dict[str, set[str]]:
    """Resolve which species carry usable evidence for each strict axis.

    Rows are kept only when the evidence scope is species/synonym-direct, or --
    when the tier explicitly asks for it -- validated genus Low. Any other scope
    (family inference, global fallback, unreviewed candidate) is dropped, so a
    tier can never quietly widen into prohibited inference.
    """
    required = {"accepted_species", "trait_name", "evidence_scope"}
    missing = required - set(evidence.columns)
    if missing:
        raise typer.BadParameter(f"evidence ledger is missing columns: {sorted(missing)}")

    frame = evidence.loc[:, sorted(required)].copy()
    for column in frame.columns:
        frame[column] = _text(frame[column])

    allowed = set(DIRECT_SCOPES)
    if include_validated_low:
        allowed |= VALIDATED_LOW_SCOPES
    frame = frame.loc[frame["evidence_scope"].str.lower().isin(allowed)]

    frame["axis"] = [_axis_of(trait, axes) for trait in frame["trait_name"]]
    frame = frame.loc[frame["axis"].ne("") & frame["accepted_species"].ne("")]

    return {
        axis: set(group["accepted_species"])
        for axis, group in frame.groupby("axis", sort=True)
    }


def evaluate_tiers(
    evidence: pd.DataFrame,
    island_species: pd.DataFrame,
    axes: dict[str, Any],
    tiers: dict[str, Any],
    min_covered_species: int,
    min_covered_fraction: float,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Report analysis-ready islands per axis for every declared evidence tier."""
    records: list[dict[str, Any]] = []
    for tier_name, tier in sorted(tiers.items()):
        include_low = bool(tier.get("include_validated_low", False))
        by_axis = covered_species_by_axis(evidence, axes, include_low)
        for axis in sorted(axes):
            covered = by_axis.get(axis, set())
            readiness = island_readiness(
                island_species, covered, min_covered_species, min_covered_fraction
            )
            record = _readiness_summary(readiness, len(covered))
            record["tier"] = tier_name
            record["role"] = str(tier.get("role", ""))
            record["axis"] = axis
            records.append(record)

    table = pd.DataFrame.from_records(records)
    table = table.loc[
        :,
        [
            "tier",
            "role",
            "axis",
            "n_target_species",
            "median_island_covered_fraction",
            "n_islands_analysis_ready",
            "n_islands",
        ],
    ].rename(columns={"n_target_species": "n_covered_species"})
    table = table.sort_values(["role", "tier", "axis"]).reset_index(drop=True)

    primary = table.loc[table["role"] == "primary"]
    sensitivity = table.loc[table["role"] == "sensitivity"]
    summary: dict[str, Any] = {
        "version": "1.0",
        "readiness_gate": {
            "min_covered_species": int(min_covered_species),
            "min_covered_fraction": float(min_covered_fraction),
        },
        "primary_metric": "n_islands_analysis_ready under the direct-only tier",
        "primary_min_axis_ready": (
            int(primary["n_islands_analysis_ready"].min()) if len(primary) else 0
        ),
        "sensitivity_min_axis_ready": (
            int(sensitivity["n_islands_analysis_ready"].min()) if len(sensitivity) else 0
        ),
        "interpretation": (
            "Validated Low is reported only as a sensitivity tier. It amplifies a "
            "small number of direct facts into many cells and is withdrawn in bulk "
            "whenever a rebuild tightens, so it is never the primary numerator. "
            "Ready counts are not trait values and not biological absences."
        ),
    }
    return table, summary


@app.command("evaluate")
def evaluate(
    evidence_csv: Path = typer.Option(..., "--evidence-csv", exists=True),
    output_dir: Path = typer.Option(..., "--output-dir"),
    config_path: Path = typer.Option(
        Path("config/island_weighted_acquisition.yml"), "--config", exists=True
    ),
) -> None:
    """Score the real evidence ledger by analysis-ready islands, per axis and tier."""
    config = load_config(config_path)
    if "axes" not in config or "evaluation_tiers" not in config:
        raise typer.BadParameter("config must declare axes and evaluation_tiers")

    readiness_gate = config["readiness"]
    evidence = pd.read_csv(evidence_csv)
    island_species = pd.read_csv(Path(config["inputs"]["island_species"]))

    table, summary = evaluate_tiers(
        evidence,
        island_species,
        config["axes"],
        config["evaluation_tiers"],
        int(readiness_gate["min_covered_species"]),
        float(readiness_gate["min_covered_fraction"]),
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_dir / "island_readiness_by_tier_axis.csv", index=False)
    (output_dir / "island_readiness_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    typer.echo(table.to_string(index=False))
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


def source_yield(
    baseline: set[str],
    candidate: set[str],
    island_species: pd.DataFrame,
    species_value: pd.DataFrame,
    budget: int,
    min_covered_species: int,
    min_covered_fraction: float,
) -> dict[str, Any]:
    """Score one source by the analysis-ready islands it adds over a baseline.

    This is the measurement the strategy turns on. A source that adds thousands
    of species but no analysis-ready islands has not moved the analysis, and a
    source that adds few species inside the priority head can move it a lot.
    """
    target = set(
        species_value.loc[species_value["priority_rank"] <= int(budget), "accepted_species"]
    )

    before = island_readiness(
        island_species, baseline, min_covered_species, min_covered_fraction
    )
    after = island_readiness(
        island_species, baseline | candidate, min_covered_species, min_covered_fraction
    )

    ready_before = set(before.loc[before["analysis_ready"], "island_id"])
    ready_after = set(after.loc[after["analysis_ready"], "island_id"])

    new_species = candidate - baseline
    in_target = new_species & target
    mass = species_value.loc[
        species_value["accepted_species"].isin(new_species), "flora_mass_share"
    ].sum()

    return {
        "target_budget": int(budget),
        "n_candidate_species": len(candidate),
        "n_new_species": len(new_species),
        "n_new_species_in_priority_head": len(in_target),
        "priority_head_hit_rate": (len(in_target) / len(new_species)) if new_species else 0.0,
        "new_flora_mass_share": float(mass),
        "n_islands_analysis_ready_before": len(ready_before),
        "n_islands_analysis_ready_after": len(ready_after),
        "n_islands_analysis_ready_gained": len(ready_after - ready_before),
        "n_islands_analysis_ready_lost": len(ready_before - ready_after),
        "interpretation": (
            "Net analysis-ready island gain is the acceptance test for a source. "
            "Species counts and cell counts are inputs to it, not substitutes."
        ),
    }


@app.command("source-yield")
def source_yield_command(
    candidate_csv: Path = typer.Option(..., "--candidate-csv", exists=True),
    output_dir: Path = typer.Option(..., "--output-dir"),
    config_path: Path = typer.Option(
        Path("config/island_weighted_acquisition.yml"), "--config", exists=True
    ),
    baseline_csv: Path | None = typer.Option(None, "--baseline-csv"),
    budget: int = typer.Option(10000, "--budget"),
) -> None:
    """Report what one source adds, in analysis-ready islands, over a baseline."""
    config = load_config(config_path)
    readiness_gate = config["readiness"]

    master = pd.read_csv(Path(config["inputs"]["master_taxa"]))
    island_species = pd.read_csv(Path(config["inputs"]["island_species"]))
    species_value = build_species_value(master, config["scope_classes"])

    baseline = _read_covered(baseline_csv) if baseline_csv is not None else set()
    candidate = _read_covered(candidate_csv)

    summary = source_yield(
        baseline,
        candidate,
        island_species,
        species_value,
        budget,
        int(readiness_gate["min_covered_species"]),
        float(readiness_gate["min_covered_fraction"]),
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "source_yield_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


def _read_covered(path: Path) -> set[str]:
    frame = pd.read_csv(path)
    for column in ("accepted_species", "species", "accepted_name"):
        if column in frame.columns:
            return set(_text(frame[column]).loc[lambda s: s.ne("")])
    raise typer.BadParameter(
        f"{path} must carry an accepted_species, species or accepted_name column"
    )


@app.command("run")
def run(
    config_path: Path = typer.Option(
        Path("config/island_weighted_acquisition.yml"), "--config", exists=True
    ),
    output_dir: Path = typer.Option(..., "--output-dir"),
    covered_species_csv: Path | None = typer.Option(
        None,
        "--covered-species-csv",
        help=(
            "Optional list of species already carrying accepted evidence. When given, "
            "island readiness is also reported for the real current coverage."
        ),
    ),
) -> None:
    """Write the species value ledger, the planning frontier, and island readiness."""
    config = load_config(config_path)
    inputs = config["inputs"]
    readiness_gate = config["readiness"]
    min_covered_species = int(readiness_gate["min_covered_species"])
    min_covered_fraction = float(readiness_gate["min_covered_fraction"])

    master = pd.read_csv(Path(inputs["master_taxa"]))
    island_species = pd.read_csv(Path(inputs["island_species"]))

    species_value = build_species_value(master, config["scope_classes"])
    targets = [int(value) for value in config["frontier_targets"]]
    frontier, summary = build_frontier(
        species_value, island_species, targets, min_covered_species, min_covered_fraction
    )
    breakdown = pd.concat(
        [scope_breakdown(species_value, budget) for budget in targets], ignore_index=True
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    species_value.to_csv(output_dir / "species_acquisition_value.csv.gz", index=False)
    frontier.to_csv(output_dir / "acquisition_strategy_frontier.csv", index=False)
    breakdown.to_csv(output_dir / "target_scope_breakdown.csv", index=False)

    if covered_species_csv is not None:
        covered = _read_covered(covered_species_csv)
        current = island_readiness(
            island_species, covered, min_covered_species, min_covered_fraction
        )
        current.to_csv(output_dir / "current_island_readiness.csv.gz", index=False)
        summary["current_coverage"] = _readiness_summary(current, len(covered))

    (output_dir / "acquisition_value_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    app()
