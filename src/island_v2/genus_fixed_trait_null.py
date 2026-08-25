"""Genus-composition-preserving null model for Chapter 1 trait analyses.

Direct species-level trait states are permuted among species within the same genus.
Island occurrence and genus composition remain fixed. This is a lineage sensitivity,
not missing-data imputation and not evidence of within-lineage evolution.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.flora_status_support import STRATA, _text, attach_floristic_status, stratum_mask

app = typer.Typer(add_completion=False, no_args_is_help=True)


def classify_binary_value(value: object, positive: set[str], negative: set[str]) -> float:
    tokens = {token.strip() for token in str(value or "").split("|") if token.strip()}
    if not tokens:
        return np.nan
    if tokens <= positive:
        return 1.0
    if tokens <= negative:
        return 0.0
    return np.nan


def direct_binary_species(
    direct_evidence: pd.DataFrame,
    master_taxa: pd.DataFrame,
    spec: dict[str, Any],
) -> pd.DataFrame:
    required = {"accepted_species", "trait_name", "normalized_value"}
    missing = required - set(direct_evidence.columns)
    if missing:
        raise typer.BadParameter(f"direct evidence missing columns: {sorted(missing)}")
    if not {"accepted_species", "genus"}.issubset(master_taxa.columns):
        raise typer.BadParameter("master taxa must contain accepted_species and genus")

    trait_name = str(spec["trait_name"])
    positive = {str(v) for v in spec["positive_states"]}
    negative = {str(v) for v in spec["negative_states"]}
    frame = direct_evidence.loc[
        direct_evidence["trait_name"].astype(str).eq(trait_name),
        ["accepted_species", "normalized_value"],
    ].copy()
    frame["accepted_species"] = _text(frame["accepted_species"])
    frame["state"] = [classify_binary_value(v, positive, negative) for v in frame["normalized_value"]]
    frame = frame.dropna(subset=["state"])

    # Fail closed on species with contradictory direct records.
    collapsed = frame.groupby("accepted_species", as_index=False).agg(
        n_states=("state", "nunique"), state=("state", "first")
    )
    collapsed = collapsed.loc[collapsed["n_states"].eq(1), ["accepted_species", "state"]]

    taxa = master_taxa[["accepted_species", "genus"]].drop_duplicates("accepted_species").copy()
    taxa["accepted_species"] = _text(taxa["accepted_species"])
    taxa["genus"] = _text(taxa["genus"])
    return collapsed.merge(taxa, on="accepted_species", how="left", validate="one_to_one").loc[
        lambda x: x["genus"].ne("")
    ].reset_index(drop=True)


def _permuted_state_map(species_states: pd.DataFrame, rng: np.random.Generator) -> pd.Series:
    pieces: list[pd.Series] = []
    for _, group in species_states.groupby("genus", sort=True):
        species = group["accepted_species"].to_numpy(dtype=object)
        states = group["state"].to_numpy(dtype=float).copy()
        rng.shuffle(states)
        pieces.append(pd.Series(states, index=species))
    return pd.concat(pieces) if pieces else pd.Series(dtype=float)


def _aggregate(flora: pd.DataFrame, state_map: pd.Series, stratum: str) -> pd.DataFrame:
    subset = flora.loc[stratum_mask(flora, stratum), ["island_id", "accepted_species"]].copy()
    subset["state"] = subset["accepted_species"].map(state_map)
    subset = subset.dropna(subset=["state"])
    if subset.empty:
        return pd.DataFrame(columns=["island_id", "n_species", "successes", "share"])
    out = subset.groupby("island_id", as_index=False).agg(
        n_species=("state", "size"), successes=("state", "sum")
    )
    out["share"] = out["successes"] / out["n_species"]
    return out


def run_genus_fixed_null(
    island_species: pd.DataFrame,
    status_ledger: pd.DataFrame,
    direct_evidence: pd.DataFrame,
    master_taxa: pd.DataFrame,
    binary_outcomes: dict[str, Any],
    *,
    n_draws: int = 1000,
    seed: int = 20260825,
    min_species_per_island: int = 1,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if n_draws < 10:
        raise typer.BadParameter("n_draws must be at least 10")
    flora = attach_floristic_status(island_species, status_ledger)
    rng = np.random.default_rng(seed)
    result_rows: list[pd.DataFrame] = []
    audit_rows: list[dict[str, Any]] = []

    for outcome, spec in binary_outcomes.items():
        states = direct_binary_species(direct_evidence, master_taxa, spec)
        observed_map = states.set_index("accepted_species")["state"]
        genus_sizes = states.groupby("genus")["accepted_species"].nunique()
        informative_genera = set(genus_sizes.loc[genus_sizes >= 2].index)
        permutable = states.loc[states["genus"].isin(informative_genera)].copy()
        audit_rows.append({
            "outcome": outcome,
            "n_direct_binary_species": int(len(states)),
            "n_genera": int(states["genus"].nunique()),
            "n_permutable_species": int(len(permutable)),
            "n_permutable_genera": int(len(informative_genera)),
        })
        if permutable.empty:
            continue

        for stratum in STRATA:
            observed = _aggregate(flora, observed_map, stratum).rename(
                columns={"share": "observed_share", "n_species": "observed_n_species"}
            )
            observed = observed.loc[observed["observed_n_species"] >= min_species_per_island].copy()
            if observed.empty:
                continue

            draws: list[pd.DataFrame] = []
            for draw in range(1, n_draws + 1):
                permuted_map = observed_map.copy()
                permuted = _permuted_state_map(permutable, rng)
                permuted_map.loc[permuted.index] = permuted
                agg = _aggregate(flora, permuted_map, stratum)
                agg = agg.loc[agg["island_id"].isin(observed["island_id"]), ["island_id", "share"]]
                agg["draw"] = draw
                draws.append(agg)
            draw_table = pd.concat(draws, ignore_index=True)
            null = draw_table.groupby("island_id", as_index=False).agg(
                null_mean=("share", "mean"),
                null_sd=("share", "std"),
                null_q025=("share", lambda x: float(np.quantile(x, 0.025))),
                null_q975=("share", lambda x: float(np.quantile(x, 0.975))),
            )
            merged = observed.merge(null, on="island_id", how="inner", validate="one_to_one")
            merged["deviation_observed_minus_null"] = merged["observed_share"] - merged["null_mean"]
            merged["ses"] = merged["deviation_observed_minus_null"] / merged["null_sd"].replace(0, np.nan)
            merged["outcome"] = outcome
            merged["stratum"] = stratum
            result_rows.append(merged)

    results = pd.concat(result_rows, ignore_index=True) if result_rows else pd.DataFrame()
    return results, pd.DataFrame(audit_rows)


def _load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not payload.get("binary_outcomes"):
        raise typer.BadParameter("config must contain binary_outcomes")
    return payload


@app.command("run")
def run(
    island_species_csv: Path = typer.Option(..., exists=True),
    status_ledger_csv: Path = typer.Option(..., exists=True),
    direct_evidence_csv: Path = typer.Option(..., exists=True),
    master_taxa_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/flora_status_support.yml"), exists=True),
    n_draws: int = typer.Option(1000, min=10),
    seed: int = typer.Option(20260825),
    min_species_per_island: int = typer.Option(1, min=1),
) -> None:
    config = _load_config(config_path)
    results, audit = run_genus_fixed_null(
        pd.read_csv(island_species_csv), pd.read_csv(status_ledger_csv),
        pd.read_csv(direct_evidence_csv), pd.read_csv(master_taxa_csv),
        config["binary_outcomes"], n_draws=n_draws, seed=seed,
        min_species_per_island=min_species_per_island,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "genus_fixed_null_by_island.csv.gz", index=False)
    audit.to_csv(output_dir / "genus_fixed_null_audit.csv", index=False)
    manifest = {
        "contract": "genus_fixed_trait_null_v2",
        "n_draws": n_draws,
        "seed": seed,
        "min_species_per_island": min_species_per_island,
        "interpretation": "Lineage sensitivity only; no missing trait fill and no causal pollinator inference.",
    }
    (output_dir / "genus_fixed_null_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
