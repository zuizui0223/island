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

    # Species with contradictory direct states fail closed and are excluded.
    collapsed = frame.groupby("accepted_species", as_index=False).agg(
        n_states=("state", "nunique"), state=("state", "first")
    )
    collapsed = collapsed.loc[collapsed["n_states"].eq(1), ["accepted_species", "state"]]

    taxa = master_taxa[["accepted_species", "genus"]].drop_duplicates("accepted_species").copy()
    taxa["accepted_species"] = _text(taxa["accepted_species"])
    taxa["genus"] = _text(taxa["genus"])
    return (
        collapsed.merge(taxa, on="accepted_species", how="left", validate="one_to_one")
        .loc[lambda x: x["genus"].ne("")]
        .sort_values(["genus", "accepted_species"])
        .reset_index(drop=True)
    )


def _permuted_state_map(species_states: pd.DataFrame, rng: np.random.Generator) -> pd.Series:
    """One-draw helper retained for transparent tests and debugging."""
    pieces: list[pd.Series] = []
    for _, group in species_states.groupby("genus", sort=True):
        species = group["accepted_species"].to_numpy(dtype=object)
        states = group["state"].to_numpy(dtype=float).copy()
        rng.shuffle(states)
        pieces.append(pd.Series(states, index=species))
    return pd.concat(pieces) if pieces else pd.Series(dtype=float)


def _draw_state_matrix(
    states: pd.DataFrame,
    *,
    n_draws: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Return species x draw states, permuting only within genus.

    This replaces the historical draw-by-draw island groupby loop. Each genus gets
    an independent random rank permutation in every draw while single-species genera
    remain fixed. The statistical null is unchanged; only execution is vectorized.
    """
    observed = states["state"].to_numpy(dtype=np.float32)
    draws = np.repeat(observed[:, None], n_draws, axis=1)
    for _, group in states.groupby("genus", sort=True):
        idx = group.index.to_numpy(dtype=int)
        if len(idx) < 2:
            continue
        order = np.argsort(rng.random((len(idx), n_draws)), axis=0)
        draws[idx, :] = observed[idx][order]
    return draws


def _stratum_null_summary(
    flora: pd.DataFrame,
    states: pd.DataFrame,
    state_draws: np.ndarray,
    *,
    stratum: str,
    min_species_per_island: int,
) -> pd.DataFrame:
    species_to_index = pd.Series(
        np.arange(len(states), dtype=int), index=states["accepted_species"].to_numpy()
    )
    subset = flora.loc[
        stratum_mask(flora, stratum), ["island_id", "accepted_species"]
    ].drop_duplicates()
    subset["species_index"] = subset["accepted_species"].map(species_to_index)
    subset = subset.dropna(subset=["species_index"]).copy()
    if subset.empty:
        return pd.DataFrame()
    subset["species_index"] = subset["species_index"].astype(int)

    observed_states = states["state"].to_numpy(dtype=float)
    rows: list[dict[str, Any]] = []
    for island_id, group in subset.groupby("island_id", sort=True):
        indices = np.unique(group["species_index"].to_numpy(dtype=int))
        n_species = int(len(indices))
        if n_species < min_species_per_island:
            continue
        observed_share = float(observed_states[indices].mean())
        null_shares = state_draws[indices, :].mean(axis=0, dtype=np.float64)
        null_mean = float(np.mean(null_shares))
        null_sd = float(np.std(null_shares, ddof=1)) if len(null_shares) > 1 else np.nan
        q025, q975 = np.quantile(null_shares, [0.025, 0.975])
        deviation = observed_share - null_mean
        rows.append(
            {
                "island_id": island_id,
                "observed_n_species": n_species,
                "observed_share": observed_share,
                "null_mean": null_mean,
                "null_sd": null_sd,
                "null_q025": float(q025),
                "null_q975": float(q975),
                "deviation_observed_minus_null": deviation,
                "ses": deviation / null_sd if np.isfinite(null_sd) and null_sd > 0 else np.nan,
            }
        )
    return pd.DataFrame(rows)


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
        states = direct_binary_species(direct_evidence, master_taxa, spec).reset_index(drop=True)
        genus_sizes = states.groupby("genus")["accepted_species"].nunique()
        informative_genera = set(genus_sizes.loc[genus_sizes >= 2].index)
        n_permutable = int(states["genus"].isin(informative_genera).sum())
        audit_rows.append(
            {
                "outcome": outcome,
                "n_direct_binary_species": int(len(states)),
                "n_genera": int(states["genus"].nunique()),
                "n_permutable_species": n_permutable,
                "n_permutable_genera": int(len(informative_genera)),
            }
        )
        if states.empty or not informative_genera:
            continue

        state_draws = _draw_state_matrix(states, n_draws=n_draws, rng=rng)
        for stratum in STRATA:
            summary = _stratum_null_summary(
                flora,
                states,
                state_draws,
                stratum=stratum,
                min_species_per_island=min_species_per_island,
            )
            if summary.empty:
                continue
            summary["outcome"] = outcome
            summary["stratum"] = stratum
            result_rows.append(summary)

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
        pd.read_csv(island_species_csv),
        pd.read_csv(status_ledger_csv),
        pd.read_csv(direct_evidence_csv),
        pd.read_csv(master_taxa_csv),
        config["binary_outcomes"],
        n_draws=n_draws,
        seed=seed,
        min_species_per_island=min_species_per_island,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "genus_fixed_null_by_island.csv.gz", index=False)
    audit.to_csv(output_dir / "genus_fixed_null_audit.csv", index=False)
    manifest = {
        "contract": "genus_fixed_trait_null_v3_vectorized",
        "n_draws": n_draws,
        "seed": seed,
        "min_species_per_island": min_species_per_island,
        "execution": "vectorized species-by-draw genus permutations with island aggregation",
        "interpretation": "Lineage sensitivity only; no missing trait fill and no causal pollinator inference.",
    }
    (output_dir / "genus_fixed_null_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
