"""Category-preserving direct-trait decomposition after floristic-status and genus control.

This is the lightweight decision layer before any final INLA fit. It keeps the
predeclared four colour and four floral-form categories rather than collapsing
them to a single binary syndrome. For each category, the observed island share is
compared with the exact expectation implied by the direct trait distribution of
the same genera. Thus island species membership and genus composition stay fixed.

Only direct/synonym-direct, resolved species evidence is admitted. Ambiguous or
cross-category species are excluded rather than forced into a category.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer

from island_v2.status_stratified_lineage_analysis import (
    DEFAULT_BASELINE,
    fit_weighted_linear_clustered,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)
DIRECT_SCOPES = {"species_direct", "synonym_direct"}

CATEGORY_SPECS: dict[str, dict[str, Any]] = {
    "colour": {
        "trait_name": "flower_primary_color",
        "categories": {
            "plain": {"white", "green_brown_inconspicuous"},
            "yellow_orange": {"yellow_orange"},
            "red_pink": {"red_pink"},
            "blue_purple": {"blue_purple"},
        },
    },
    "form": {
        "trait_name": "floral_form",
        "categories": {
            "open_generalized": {"open_radial"},
            "tubular_trumpet": {
                "tubular",
                "bell_campanulate",
                "funnel_trumpet",
                "urn_urceolate",
                "salverform",
            },
            "zygomorphic_specialized": {"bilabiate", "papilionaceous", "spurred"},
            "composite_brush": {"composite_head", "brush_puff"},
        },
    },
}
DEFAULT_STRATA = ("all_native", "native_nonendemic", "endemic")
DEFAULT_REGIMES = (
    "northern_midlatitude",
    "tropical",
    "southern_extratropical",
)


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _json_safe_specs(specs: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Return a deterministic JSON-safe copy of the category contract."""
    safe: dict[str, dict[str, Any]] = {}
    for domain, spec in specs.items():
        safe[domain] = {
            "trait_name": str(spec["trait_name"]),
            "categories": {
                str(name): sorted(str(state) for state in states)
                for name, states in spec["categories"].items()
            },
        }
    return safe


def _classify_value(value: object, categories: dict[str, set[str]]) -> str | None:
    tokens = {token.strip() for token in str(value or "").split("|") if token.strip()}
    if not tokens:
        return None
    hits = [name for name, states in categories.items() if tokens <= states]
    return hits[0] if len(hits) == 1 else None


def direct_categorical_species(
    evidence: pd.DataFrame,
    master_taxa: pd.DataFrame,
    spec: dict[str, Any],
) -> pd.DataFrame:
    required = {"accepted_species", "trait_name", "normalized_value"}
    missing = required - set(evidence.columns)
    if missing:
        raise typer.BadParameter(f"direct evidence missing columns: {sorted(missing)}")
    if not {"accepted_species", "genus"}.issubset(master_taxa.columns):
        raise typer.BadParameter("master taxa must contain accepted_species and genus")

    work = evidence.copy()
    work["accepted_species"] = _text(work["accepted_species"])
    work["trait_name"] = _text(work["trait_name"])
    if "evidence_scope" in work.columns:
        work = work.loc[_text(work["evidence_scope"]).str.lower().isin(DIRECT_SCOPES)].copy()
    if "resolution_status" in work.columns:
        work = work.loc[_text(work["resolution_status"]).str.lower().eq("resolved")].copy()
    work = work.loc[work["trait_name"].eq(str(spec["trait_name"]))].copy()
    categories = {name: set(states) for name, states in spec["categories"].items()}
    work["category"] = [
        _classify_value(value, categories) for value in work["normalized_value"]
    ]
    work = work.dropna(subset=["category"])

    agreement = work.groupby("accepted_species")["category"].nunique()
    conflicts = set(agreement.loc[agreement.gt(1)].index)
    work = work.loc[~work["accepted_species"].isin(conflicts)]
    work = work[["accepted_species", "category"]].drop_duplicates("accepted_species")

    taxa = master_taxa[["accepted_species", "genus"]].drop_duplicates("accepted_species").copy()
    taxa["accepted_species"] = _text(taxa["accepted_species"])
    taxa["genus"] = _text(taxa["genus"])
    work = work.merge(taxa, on="accepted_species", how="left", validate="one_to_one")
    return work.loc[work["genus"].ne("")].reset_index(drop=True)


def _stratum_mask(status: pd.DataFrame, stratum: str) -> pd.Series:
    origin = _text(status["origin_status"]).str.lower()
    endemic = _text(status["endemic_status"]).str.lower()
    if stratum == "all_native":
        return origin.eq("native")
    if stratum == "native_nonendemic":
        return origin.eq("native") & endemic.eq("nonendemic")
    if stratum == "endemic":
        return origin.eq("native") & endemic.eq("endemic")
    raise typer.BadParameter(f"unsupported stratum: {stratum}")


def build_category_residuals(
    status_ledger: pd.DataFrame,
    evidence: pd.DataFrame,
    master_taxa: pd.DataFrame,
    *,
    specs: dict[str, dict[str, Any]] = CATEGORY_SPECS,
    strata: tuple[str, ...] = DEFAULT_STRATA,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required_status = {"island_id", "accepted_species", "origin_status", "endemic_status"}
    missing = required_status - set(status_ledger.columns)
    if missing:
        raise typer.BadParameter(f"status ledger missing columns: {sorted(missing)}")
    status = status_ledger.copy()
    status["island_id"] = _text(status["island_id"])
    status["accepted_species"] = _text(status["accepted_species"])
    status = status.drop_duplicates(["island_id", "accepted_species"])

    result_rows: list[pd.DataFrame] = []
    audit_rows: list[dict[str, Any]] = []
    for domain, spec in specs.items():
        species = direct_categorical_species(evidence, master_taxa, spec)
        categories = list(spec["categories"])
        genus_counts = (
            species.groupby(["genus", "category"], as_index=False)
            .agg(n_category=("accepted_species", "size"))
        )
        genus_total = species.groupby("genus")["accepted_species"].nunique().rename("n_genus")
        genus_probs = genus_counts.merge(genus_total, on="genus", how="left")
        genus_probs["expected_probability"] = genus_probs["n_category"] / genus_probs["n_genus"]
        probability_lookup = {
            category: genus_probs.loc[genus_probs["category"].eq(category)]
            .set_index("genus")["expected_probability"]
            for category in categories
        }
        audit_rows.append(
            {
                "domain": domain,
                "trait_name": spec["trait_name"],
                "n_direct_categorical_species": int(len(species)),
                "n_genera": int(species["genus"].nunique()),
                "n_categories": len(categories),
            }
        )

        for stratum in strata:
            flora = status.loc[
                _stratum_mask(status, stratum), ["island_id", "accepted_species"]
            ].copy()
            joined = flora.merge(species, on="accepted_species", how="inner")
            if joined.empty:
                continue
            for category in categories:
                work = joined.copy()
                work["observed"] = work["category"].eq(category).astype(float)
                work["expected"] = work["genus"].map(probability_lookup[category]).fillna(0.0)
                island = (
                    work.groupby("island_id", as_index=False)
                    .agg(
                        trials=("accepted_species", "size"),
                        successes=("observed", "sum"),
                        expected_successes=("expected", "sum"),
                    )
                )
                island["observed_share"] = island["successes"] / island["trials"]
                island["expected_genus_share"] = island["expected_successes"] / island["trials"]
                island["deviation_observed_minus_genus"] = (
                    island["observed_share"] - island["expected_genus_share"]
                )
                island["domain"] = domain
                island["category"] = category
                island["stratum"] = stratum
                result_rows.append(island)

    results = pd.concat(result_rows, ignore_index=True) if result_rows else pd.DataFrame()
    return results, pd.DataFrame(audit_rows)


def fit_category_isolation_models(
    residuals: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    regimes: tuple[str, ...] = DEFAULT_REGIMES,
    baseline: list[str] | None = None,
    cluster: str = "spatial_block",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    predictors = list(baseline or DEFAULT_BASELINE)
    required_cov = {"island_id", "analysis_regime", cluster, *predictors}
    missing = required_cov - set(covariates.columns)
    if missing:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")
    data = residuals.merge(
        covariates[list(required_cov)].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    coefficient_tables: list[pd.DataFrame] = []
    support_rows: list[dict[str, Any]] = []
    keys = data[["stratum", "domain", "category"]].drop_duplicates()
    for key in keys.itertuples(index=False):
        for regime in regimes:
            subset = data.loc[
                data["stratum"].eq(key.stratum)
                & data["domain"].eq(key.domain)
                & data["category"].eq(key.category)
                & data["analysis_regime"].eq(regime)
            ].copy()
            n = int(len(subset))
            support_rows.append(
                {
                    "stratum": key.stratum,
                    "domain": key.domain,
                    "category": key.category,
                    "regime": regime,
                    "n_islands": n,
                    "n_spatial_blocks": int(subset[cluster].dropna().nunique()),
                    "median_direct_trials": float(subset["trials"].median()) if n else np.nan,
                    "support_class": (
                        "confirmatory_count_met" if n >= 50
                        else "targeted_acquisition_zone" if n >= 30
                        else "below_pilot_support"
                    ),
                }
            )
            coefficients, _ = fit_weighted_linear_clustered(
                subset,
                response_column="deviation_observed_minus_genus",
                weight_column="trials",
                predictors=predictors,
                cluster_column=cluster,
            )
            if coefficients.empty:
                continue
            coefficients.insert(0, "category", key.category)
            coefficients.insert(0, "domain", key.domain)
            coefficients.insert(0, "regime", regime)
            coefficients.insert(0, "stratum", key.stratum)
            coefficient_tables.append(coefficients)
    return (
        pd.concat(coefficient_tables, ignore_index=True) if coefficient_tables else pd.DataFrame(),
        pd.DataFrame(support_rows),
    )


@app.command("run")
def run(
    status_ledger_csv: Path = typer.Option(..., exists=True),
    direct_evidence_csv: Path = typer.Option(..., exists=True),
    master_taxa_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    residuals, audit = build_category_residuals(
        pd.read_csv(status_ledger_csv),
        pd.read_csv(direct_evidence_csv),
        pd.read_csv(master_taxa_csv),
    )
    coefficients, support = fit_category_isolation_models(
        residuals, pd.read_csv(covariates_csv)
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    residuals.to_csv(output_dir / "status_category_residuals.csv.gz", index=False)
    audit.to_csv(output_dir / "status_category_audit.csv", index=False)
    coefficients.to_csv(output_dir / "status_category_coefficients.csv", index=False)
    support.to_csv(output_dir / "status_category_support.csv", index=False)
    isolation = (
        coefficients.loc[coefficients["predictor"].eq("log_distance_to_continent_km")].copy()
        if not coefficients.empty else coefficients.copy()
    )
    isolation.to_csv(output_dir / "status_category_isolation.csv", index=False)
    manifest = {
        "contract": "status_category_genus_expectation_v1",
        "evidence": "direct species evidence only",
        "status": "GIFT origin plus WCVP-corroborated regional endemism",
        "lineage_control": "exact category expectation from direct species in the same genus",
        "domains": _json_safe_specs(CATEGORY_SPECS),
        "guardrail": (
            "This is a category-preserving decision layer before any final INLA model. "
            "A category residual slope is not by itself a pollination mechanism."
        ),
    }
    (output_dir / "status_category_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
