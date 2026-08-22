"""Direct-evidence trait analysis restricted to source-backed native flora.

This is the status-aware M0 checkpoint. It deliberately answers a narrow
question before adding lineage or Bombus terms: within each broad analysis
regime, is isolation associated with native floral/reproductive composition
after island area and climate are represented?

Only species-direct/synonym-direct resolved evidence is used. Multistate records
are mapped to a binary contrast only when every reported state falls on the same
side of the contrast; cross-side or otherwise ambiguous records are excluded.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2.m0_m4_full_analysis import fit_grouped_binomial_clustered

app = typer.Typer(add_completion=False, no_args_is_help=True)

DIRECT_SCOPES = {"species_direct", "synonym_direct"}

CONTRASTS = {
    "plain_colour": {
        "trait_name": "flower_primary_color",
        "positive": {"white", "green_brown_inconspicuous"},
        "negative": {"yellow_orange", "red_pink", "blue_purple"},
    },
    "generalized_form": {
        "trait_name": "floral_form",
        "positive": {"open_radial", "brush_puff", "composite_head"},
        "negative": {
            "tubular",
            "bell_campanulate",
            "funnel_trumpet",
            "urn_urceolate",
            "bilabiate",
            "salverform",
            "papilionaceous",
            "spurred",
        },
    },
    "self_compatibility": {
        "trait_name": "self_incompatibility",
        "positive": {"SC"},
        "negative": {"SI"},
    },
}

DEFAULT_BASELINE = [
    "log_distance_to_continent_km",
    "log_island_area_km2",
    "climate_pc1",
    "climate_pc2",
    "climate_pc3",
    "climate_pc4",
]


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def direct_binary_species(
    evidence: pd.DataFrame,
    *,
    trait_name: str,
    positive: set[str],
    negative: set[str],
) -> pd.DataFrame:
    required = {"accepted_species", "trait_name", "normalized_value"}
    missing = required - set(evidence.columns)
    if missing:
        raise typer.BadParameter(f"direct evidence missing columns: {sorted(missing)}")
    work = evidence.copy()
    work["accepted_species"] = _text(work["accepted_species"])
    work["trait_name"] = _text(work["trait_name"])
    if "resolution_status" in work.columns:
        work = work.loc[_text(work["resolution_status"]).str.lower().eq("resolved")].copy()
    if "evidence_scope" in work.columns:
        work = work.loc[
            _text(work["evidence_scope"]).str.lower().isin(DIRECT_SCOPES)
        ].copy()
    work = work.loc[work["trait_name"].eq(trait_name)].copy()

    rows: list[dict[str, Any]] = []
    for row in work.itertuples(index=False):
        states = {
            token.strip()
            for token in str(row.normalized_value).split("|")
            if token.strip() and token.strip().lower() != "nan"
        }
        if states and states.issubset(positive):
            state = 1
        elif states and states.issubset(negative):
            state = 0
        else:
            continue
        rows.append({"accepted_species": row.accepted_species, "state": state})

    result = pd.DataFrame(rows)
    if result.empty:
        return pd.DataFrame(columns=["accepted_species", "state"])
    counts = result.groupby("accepted_species")["state"].nunique()
    conflicts = set(counts.loc[counts.gt(1)].index)
    return (
        result.loc[~result["accepted_species"].isin(conflicts)]
        .drop_duplicates("accepted_species")
        .sort_values("accepted_species")
        .reset_index(drop=True)
    )


def confirmed_native_flora(status_ledger: pd.DataFrame) -> pd.DataFrame:
    required = {"island_id", "accepted_species", "origin_status"}
    missing = required - set(status_ledger.columns)
    if missing:
        raise typer.BadParameter(f"status ledger missing columns: {sorted(missing)}")
    frame = status_ledger.copy()
    frame["island_id"] = _text(frame["island_id"])
    frame["accepted_species"] = _text(frame["accepted_species"])
    frame["origin_status"] = _text(frame["origin_status"]).str.lower()
    frame = frame.loc[frame["origin_status"].eq("native")]
    return frame[["island_id", "accepted_species"]].drop_duplicates()


def build_native_counts(
    native_flora: pd.DataFrame,
    evidence: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for outcome, spec in CONTRASTS.items():
        species = direct_binary_species(evidence, **spec)
        joined = native_flora.merge(species, on="accepted_species", how="inner")
        if joined.empty:
            continue
        counts = (
            joined.groupby("island_id", as_index=False)
            .agg(successes=("state", "sum"), trials=("state", "size"))
        )
        counts["outcome"] = outcome
        rows.append(counts)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def fit_regime_models(
    counts: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    regimes: tuple[str, ...] = (
        "northern_midlatitude",
        "tropical",
        "southern_extratropical",
    ),
    baseline: list[str] | None = None,
    cluster_column: str = "spatial_block",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    predictors = list(baseline or DEFAULT_BASELINE)
    required_cov = {"island_id", "analysis_regime", cluster_column, *predictors}
    missing = required_cov - set(covariates.columns)
    if missing:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")
    data = counts.merge(
        covariates[list(required_cov)].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )

    coefficient_rows: list[pd.DataFrame] = []
    fit_rows: list[dict[str, Any]] = []
    for regime in regimes:
        for outcome in CONTRASTS:
            subset = data.loc[
                data["analysis_regime"].eq(regime) & data["outcome"].eq(outcome)
            ].copy()
            if subset.empty:
                continue
            coefficients, fit = fit_grouped_binomial_clustered(
                subset,
                successes_column="successes",
                trials_column="trials",
                predictors=predictors,
                cluster_column=cluster_column,
                z_threshold=1.96,
            )
            if not coefficients.empty:
                coefficients.insert(0, "outcome", outcome)
                coefficients.insert(0, "regime", regime)
                coefficient_rows.append(coefficients)
            fit_rows.append({"regime": regime, "outcome": outcome, **fit})
    coefficients = (
        pd.concat(coefficient_rows, ignore_index=True)
        if coefficient_rows
        else pd.DataFrame()
    )
    return coefficients, pd.DataFrame(fit_rows)


@app.command("run")
def run(
    status_ledger_csv: Path = typer.Option(..., exists=True),
    direct_evidence_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    status = pd.read_csv(status_ledger_csv)
    evidence = pd.read_csv(direct_evidence_csv)
    covariates = pd.read_csv(covariates_csv)
    native = confirmed_native_flora(status)
    counts = build_native_counts(native, evidence)
    coefficients, fits = fit_regime_models(counts, covariates)

    output_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(output_dir / "native_direct_binary_counts.csv", index=False)
    coefficients.to_csv(output_dir / "native_direct_m0_coefficients.csv", index=False)
    fits.to_csv(output_dir / "native_direct_m0_fit.csv", index=False)

    distance = coefficients.loc[
        coefficients["predictor"].eq("log_distance_to_continent_km")
    ].copy()
    distance.to_csv(output_dir / "native_direct_isolation_coefficients.csv", index=False)
    manifest = {
        "contract": "native_direct_regime_m0_v1",
        "evidence": "direct species evidence only",
        "flora": "source-backed confirmed native only",
        "outcomes": list(CONTRASTS),
        "baseline": DEFAULT_BASELINE,
        "n_native_island_species_rows": int(len(native)),
        "n_model_count_rows": int(len(counts)),
        "interpretation": (
            "This is the pre-lineage, pre-Bombus M0 checkpoint. A weak northern raw "
            "isolation slope does not reject the conditional Bombus hypothesis."
        ),
    }
    (output_dir / "native_direct_m0_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
