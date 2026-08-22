"""Secondary direct-evidence analysis of introduced/naturalized island flora.

Introduced taxa are not assumed to be a null control. Human transport, planting,
horticultural selection and propagule pressure can create their own geographic
structure. The stratum is therefore reported as a transport contrast: if its
isolation slopes resemble native slopes, native biological interpretations need
extra caution; if they differ, that difference is descriptive rather than a
causal control contrast.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

from island_v2.native_status_trait_analysis import build_native_counts, fit_regime_models

app = typer.Typer(add_completion=False, no_args_is_help=True)


def confirmed_introduced_flora(status_ledger: pd.DataFrame) -> pd.DataFrame:
    required = {"island_id", "accepted_species", "origin_status"}
    missing = required - set(status_ledger.columns)
    if missing:
        raise typer.BadParameter(f"status ledger missing columns: {sorted(missing)}")
    frame = status_ledger.copy()
    for column in ("island_id", "accepted_species", "origin_status"):
        frame[column] = frame[column].fillna("").astype(str).str.strip()
    frame["origin_status"] = frame["origin_status"].str.lower()
    frame = frame.loc[frame["origin_status"].eq("introduced")]
    return frame[["island_id", "accepted_species"]].drop_duplicates()


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
    introduced = confirmed_introduced_flora(status)
    counts = build_native_counts(introduced, evidence)
    coefficients, fits = fit_regime_models(counts, covariates)

    output_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(output_dir / "introduced_direct_binary_counts.csv", index=False)
    coefficients.to_csv(output_dir / "introduced_direct_m0_coefficients.csv", index=False)
    fits.to_csv(output_dir / "introduced_direct_m0_fit.csv", index=False)
    isolation = coefficients.loc[
        coefficients["predictor"].eq("log_distance_to_continent_km")
    ].copy()
    isolation.to_csv(output_dir / "introduced_direct_isolation_coefficients.csv", index=False)
    manifest = {
        "contract": "introduced_direct_regime_transport_contrast_v1",
        "evidence": "direct species evidence only",
        "flora": "source-backed introduced/naturalized only",
        "n_introduced_island_species_rows": int(len(introduced)),
        "interpretation": (
            "Introduced flora is a secondary human-transport contrast, not a clean "
            "negative control. Geographic slopes may reflect propagule pressure, "
            "horticultural selection, planting history, or environmental filtering."
        ),
    }
    (output_dir / "introduced_direct_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
