"""Multiplicity-controlled summary of Chapter 1 trait-centric probability curves."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _bh(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=values.index, dtype=float)
    ok = p.notna()
    if not ok.any():
        return out
    x = p.loc[ok].to_numpy(float)
    order = np.argsort(x)
    ranked = x[order]
    n = len(ranked)
    adjusted = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    restored = np.empty(n, dtype=float)
    restored[order] = np.clip(adjusted, 0, 1)
    out.loc[ok] = restored
    return out


def summarize_trait_probabilities(
    coefficients: pd.DataFrame,
    curves: pd.DataFrame,
    fits: pd.DataFrame,
    *,
    distance_predictor: str,
) -> pd.DataFrame:
    slope = coefficients.loc[coefficients["predictor"].astype(str).eq(distance_predictor)].copy()
    slope = slope.merge(
        fits[
            [
                "stratum",
                "trait_name",
                "trait_state",
                "context",
                "support_tier",
                "n_unique_islands",
                "n_clusters",
            ]
        ],
        on=["stratum", "trait_name", "trait_state", "context"],
        how="left",
        validate="one_to_one",
    )

    slope["q_within_stratum_confirmatory"] = np.nan
    confirmatory = slope["support_tier"].eq("confirmatory")
    for stratum, group in slope.loc[confirmatory].groupby("stratum"):
        slope.loc[group.index, "q_within_stratum_confirmatory"] = _bh(group["p_value"])
    slope["distance_fdr05"] = slope["q_within_stratum_confirmatory"].le(0.05)

    endpoints: list[dict[str, object]] = []
    for key, group in curves.groupby(
        ["stratum", "trait_name", "trait_state", "context"], sort=False
    ):
        group = group.sort_values("distance_quantile")
        if group.empty:
            continue
        low = group.iloc[0]
        high = group.iloc[-1]
        endpoints.append(
            {
                "stratum": key[0],
                "trait_name": key[1],
                "trait_state": key[2],
                "context": key[3],
                "distance_q_low": float(low["distance_quantile"]),
                "distance_low": float(low["distance_gradient_value"]),
                "probability_low": float(low["predicted_probability"]),
                "distance_q_high": float(high["distance_quantile"]),
                "distance_high": float(high["distance_gradient_value"]),
                "probability_high": float(high["predicted_probability"]),
                "probability_change_high_minus_low": float(
                    high["predicted_probability"] - low["predicted_probability"]
                ),
            }
        )
    endpoint_table = pd.DataFrame(endpoints)
    if not endpoint_table.empty:
        slope = slope.merge(
            endpoint_table,
            on=["stratum", "trait_name", "trait_state", "context"],
            how="left",
            validate="one_to_one",
        )
    return slope.sort_values(
        ["stratum", "support_tier", "q_within_stratum_confirmatory", "trait_name", "trait_state", "context"],
        na_position="last",
    ).reset_index(drop=True)


@app.command("run")
def run(
    coefficients_csv: Path = typer.Option(..., exists=True),
    curves_csv: Path = typer.Option(..., exists=True),
    fits_csv: Path = typer.Option(..., exists=True),
    output_csv: Path = typer.Option(...),
    distance_predictor: str = typer.Option("z_log_distance_to_continent_km"),
) -> None:
    summary = summarize_trait_probabilities(
        pd.read_csv(coefficients_csv),
        pd.read_csv(curves_csv),
        pd.read_csv(fits_csv),
        distance_predictor=distance_predictor,
    )
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_csv, index=False)
    typer.echo(
        f"trait probability slopes={len(summary)} confirmatory_fdr05="
        f"{int(summary.get('distance_fdr05', pd.Series(dtype=bool)).fillna(False).sum())}"
    )


if __name__ == "__main__":
    app()
