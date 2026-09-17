"""Raw flower-colour audit for the Chapter 1 v13 island-syndrome synthesis.

The current v13 paper deliberately defines its primary floral response through
accessibility/generalization rather than colour.  This audit restores the final
flower-colour database as an explicit *secondary* response without collapsing it
into ``plain_colour`` or a pollination-syndrome score.

Five reported colour states are analysed separately:

- white
- red_pink
- yellow_orange
- blue_purple
- green_brown_inconspicuous

A multistate species contributes to every focal colour it is explicitly reported to
have.  Species with only non-focal states (for example ``other_described``) do not
enter the focal-colour denominator.  The resulting island responses therefore mean
"share of colour-resolved focal species reported with this colour", not mutually
exclusive multinomial probabilities.

For each colour the audit fits both:

1. colour ~ isolation + area + climate
2. colour ~ isolation + selfing_core + area + climate

The second model asks whether an isolation-associated colour shift is reducible to
the measured reproductive-assurance core.  It is a conditional decomposition, not
causal mediation and not evidence that a named pollinator selected the colour.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_context_analysis import (
    _chi_square_sf_integer_df,
    _fit_grouped_binomial_design,
)
from island_v2.flora_status_support import stratum_mask

app = typer.Typer(add_completion=False, no_args_is_help=True)

FOCAL_COLOURS = (
    "white",
    "red_pink",
    "yellow_orange",
    "blue_purple",
    "green_brown_inconspicuous",
)

_COLOUR_RE = re.compile(r"flower_primary_color\s*=\s*(\[[^\]]*\])")


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
    restored[order] = np.clip(adjusted, 0.0, 1.0)
    out.loc[ok] = restored
    return out


def parse_colour_states(value: object) -> set[str]:
    """Return explicit raw colour categories from one species-axis composition cell."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return set()
    text = str(value).strip()
    if not text:
        return set()
    match = _COLOUR_RE.search(text)
    if not match:
        return set()
    try:
        parsed = json.loads(match.group(1))
    except json.JSONDecodeError:
        return set()
    if not isinstance(parsed, list):
        return set()
    return {str(x).strip() for x in parsed if str(x).strip()}


def _eligible_colour_species(
    species_axis: pd.DataFrame,
    *,
    evidence_scope: str,
    colours: list[str],
) -> pd.DataFrame:
    required = {"accepted_species", "axis", "trait_composition", "quality"}
    missing = required - set(species_axis.columns)
    if missing:
        raise typer.BadParameter(f"species-axis table missing columns: {sorted(missing)}")
    work = species_axis.loc[species_axis["axis"].astype(str).eq("flower_colour")].copy()
    if evidence_scope == "direct":
        work = work.loc[
            work["quality"].fillna("").astype(str).str.lower().isin({"high", "medium"})
        ].copy()
    elif evidence_scope != "all":
        raise typer.BadParameter("evidence_scope must be 'all' or 'direct'")

    focal = set(colours)
    rows: list[dict[str, Any]] = []
    for row in work[["accepted_species", "trait_composition"]].itertuples(index=False):
        states = parse_colour_states(row.trait_composition)
        retained = states & focal
        if not retained:
            continue
        rows.append(
            {
                "accepted_species": str(row.accepted_species),
                "colour_states": retained,
            }
        )
    return pd.DataFrame(rows, columns=["accepted_species", "colour_states"])


def build_raw_colour_counts(
    species_axis: pd.DataFrame,
    status_flora: pd.DataFrame,
    *,
    evidence_scope: str,
    strata: list[str],
    colours: list[str],
) -> pd.DataFrame:
    """Build island x stratum x raw-colour presence counts.

    The denominator is the number of focal-colour-resolved species on the island.
    A multistate species can count as a success for more than one colour, by design.
    """
    required_flora = {
        "island_id",
        "accepted_species",
        "origin_status",
        "endemic_status",
        "floristic_status",
    }
    missing = required_flora - set(status_flora.columns)
    if missing:
        raise typer.BadParameter(f"status flora missing columns: {sorted(missing)}")

    states = _eligible_colour_species(
        species_axis,
        evidence_scope=evidence_scope,
        colours=colours,
    )
    if states.empty:
        return pd.DataFrame(
            columns=["island_id", "stratum", "colour", "successes", "trials", "share"]
        )
    flora = status_flora.copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    joined = flora.merge(states, on="accepted_species", how="inner", validate="many_to_one")

    rows: list[dict[str, Any]] = []
    for stratum in strata:
        subset = joined.loc[stratum_mask(joined, stratum)].copy()
        subset = subset.drop_duplicates(["island_id", "accepted_species"])
        for island_id, part in subset.groupby("island_id", sort=False):
            trials = int(part["accepted_species"].nunique())
            if trials <= 0:
                continue
            state_sets = part["colour_states"].tolist()
            for colour in colours:
                successes = int(sum(colour in state_set for state_set in state_sets))
                rows.append(
                    {
                        "island_id": str(island_id),
                        "stratum": str(stratum),
                        "colour": str(colour),
                        "successes": successes,
                        "trials": trials,
                        "share": successes / trials,
                    }
                )
    return pd.DataFrame(rows)


def _standardize(values: pd.Series) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def _fit_one_colour(
    work: pd.DataFrame,
    *,
    geography: str,
    baseline: list[str],
    cluster: str,
    conditional_selfing: bool,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    predictors = [geography]
    if conditional_selfing:
        predictors.append("selfing_core")
    predictors.extend(baseline)
    complete = work.dropna(subset=["successes", "trials", *predictors, cluster]).copy()
    if complete.empty:
        return pd.DataFrame(), {"status": "no_complete_rows"}
    names = ["intercept", *[f"z_{x}" for x in predictors]]
    columns = [np.ones(len(complete), dtype=float)]
    for predictor in predictors:
        columns.append(_standardize(complete[predictor]))
    coef, fit, covariance = _fit_grouped_binomial_design(
        pd.to_numeric(complete["successes"], errors="coerce").to_numpy(float),
        pd.to_numeric(complete["trials"], errors="coerce").to_numpy(float),
        np.column_stack(columns),
        names,
        complete[cluster].astype(str).to_numpy(),
    )
    return coef, {
        **fit,
        "covariance": covariance,
        "n_unique_islands": int(complete["island_id"].nunique()),
        "median_trials": float(pd.to_numeric(complete["trials"], errors="coerce").median()),
    }


def _selfing_table(island_scores: pd.DataFrame) -> pd.DataFrame:
    required = {"island_id", "stratum", "syndrome", "syndrome_score"}
    missing = required - set(island_scores.columns)
    if missing:
        raise typer.BadParameter(f"island score table missing columns: {sorted(missing)}")
    return (
        island_scores.loc[
            island_scores["syndrome"].astype(str).eq("selfing_core"),
            ["island_id", "stratum", "syndrome_score"],
        ]
        .rename(columns={"syndrome_score": "selfing_core"})
        .drop_duplicates(["island_id", "stratum"])
    )


def fit_raw_colour_models(
    colour_counts: pd.DataFrame,
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Fit per-colour isolation models with and without selfing-core adjustment."""
    geography = str(config["geography_column"])
    context_column = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    contexts = [str(x) for x in config["contexts"]]
    strata = [str(x) for x in config["strata"]]
    colours = [str(x) for x in config.get("colours", FOCAL_COLOURS)]
    tiers = {str(k): int(v) for k, v in config["support_tiers"].items()}

    required_counts = {"island_id", "stratum", "colour", "successes", "trials"}
    missing = required_counts - set(colour_counts.columns)
    if missing:
        raise typer.BadParameter(f"raw colour counts missing columns: {sorted(missing)}")
    needed_cov = ["island_id", geography, context_column, cluster, *baseline]
    missing = set(needed_cov) - set(covariates.columns)
    if missing:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")

    data = colour_counts.merge(
        covariates[needed_cov].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    data = data.merge(
        _selfing_table(island_scores),
        on=["island_id", "stratum"],
        how="left",
        validate="many_to_one",
    )
    for column in ["successes", "trials", "selfing_core", geography, *baseline]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data[context_column] = data[context_column].fillna("").astype(str)
    data[cluster] = data[cluster].fillna("").astype(str)

    rows: list[dict[str, Any]] = []
    for stratum in strata:
        for context in contexts:
            for colour in colours:
                base = data.loc[
                    data["stratum"].astype(str).eq(stratum)
                    & data[context_column].eq(context)
                    & data["colour"].astype(str).eq(colour)
                ].copy()
                for tier, threshold in tiers.items():
                    for model, conditional in (
                        ("unconditional", False),
                        ("conditional_selfing", True),
                    ):
                        required = ["successes", "trials", geography, *baseline, cluster]
                        if conditional:
                            required.append("selfing_core")
                        complete = base.dropna(subset=required).copy()
                        n_islands = int(complete["island_id"].nunique())
                        if n_islands < threshold:
                            rows.append(
                                {
                                    "stratum": stratum,
                                    "context": context,
                                    "colour": colour,
                                    "support_tier": tier,
                                    "threshold": threshold,
                                    "model": model,
                                    "status": "not_testable",
                                    "n_unique_islands": n_islands,
                                }
                            )
                            continue
                        coef, fit = _fit_one_colour(
                            complete,
                            geography=geography,
                            baseline=baseline,
                            cluster=cluster,
                            conditional_selfing=conditional,
                        )
                        if coef.empty:
                            rows.append(
                                {
                                    "stratum": stratum,
                                    "context": context,
                                    "colour": colour,
                                    "support_tier": tier,
                                    "threshold": threshold,
                                    "model": model,
                                    "status": str(fit.get("status", "fit_failed")),
                                    "n_unique_islands": n_islands,
                                }
                            )
                            continue
                        indexed = coef.set_index("predictor")
                        distance = indexed.loc[f"z_{geography}"]
                        row: dict[str, Any] = {
                            "stratum": stratum,
                            "context": context,
                            "colour": colour,
                            "support_tier": tier,
                            "threshold": threshold,
                            "model": model,
                            "status": "fit",
                            "n_unique_islands": int(fit["n_unique_islands"]),
                            "n_clusters": int(fit["n_clusters"]),
                            "median_trials": float(fit["median_trials"]),
                            "distance_estimate": float(distance["estimate_log_odds"]),
                            "distance_se": float(distance["cluster_robust_se"]),
                            "distance_p": float(distance["p_value"]),
                        }
                        if conditional:
                            selfing = indexed.loc["z_selfing_core"]
                            row.update(
                                {
                                    "selfing_core_estimate": float(selfing["estimate_log_odds"]),
                                    "selfing_core_se": float(selfing["cluster_robust_se"]),
                                    "selfing_core_p": float(selfing["p_value"]),
                                }
                            )
                        rows.append(row)

    result = pd.DataFrame(rows)
    if not result.empty and "distance_p" in result.columns:
        result["distance_q"] = np.nan
        fit_mask = result["status"].eq("fit")
        for _, index in result.loc[fit_mask].groupby(
            ["stratum", "context", "support_tier", "model"]
        ).groups.items():
            result.loc[index, "distance_q"] = _bh(result.loc[index, "distance_p"])
    return result


def _masked_z(work: pd.DataFrame, mask: np.ndarray, column: str) -> np.ndarray:
    out = np.zeros(len(work), dtype=float)
    out[mask] = _standardize(work.loc[mask, column])
    return out


def fit_within_context_joint_omnibus(
    colour_counts: pd.DataFrame,
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Joint Wald test that the raw five-colour isolation vector differs from zero."""
    geography = str(config["geography_column"])
    context_column = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    contexts = [str(x) for x in config["contexts"]]
    strata = [str(x) for x in config["strata"]]
    colours = [str(x) for x in config.get("colours", FOCAL_COLOURS)]
    tiers = {str(k): int(v) for k, v in config["support_tiers"].items()}

    data = colour_counts.merge(
        covariates[["island_id", geography, context_column, cluster, *baseline]].drop_duplicates(
            "island_id"
        ),
        on="island_id",
        how="left",
        validate="many_to_one",
    ).merge(
        _selfing_table(island_scores),
        on=["island_id", "stratum"],
        how="left",
        validate="many_to_one",
    )
    for column in ["successes", "trials", geography, "selfing_core", *baseline]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data[context_column] = data[context_column].fillna("").astype(str)
    data[cluster] = data[cluster].fillna("").astype(str)

    rows: list[dict[str, Any]] = []
    for stratum in strata:
        for context in contexts:
            base = data.loc[
                data["stratum"].astype(str).eq(stratum) & data[context_column].eq(context)
            ].copy()
            for tier, threshold in tiers.items():
                for model, conditional in (
                    ("unconditional", False),
                    ("conditional_selfing", True),
                ):
                    required = ["successes", "trials", geography, *baseline, cluster]
                    if conditional:
                        required.append("selfing_core")
                    work = base.dropna(subset=required).copy()
                    support = work.groupby("colour")["island_id"].nunique()
                    retained = [c for c in colours if int(support.get(c, 0)) >= threshold]
                    if not retained:
                        rows.append(
                            {
                                "stratum": stratum,
                                "context": context,
                                "support_tier": tier,
                                "threshold": threshold,
                                "model": model,
                                "status": "not_testable",
                                "n_retained_colours": 0,
                            }
                        )
                        continue
                    work = work.loc[work["colour"].isin(retained)].copy()
                    names: list[str] = []
                    columns: list[np.ndarray] = []
                    distance_names: list[str] = []
                    for colour in retained:
                        mask = work["colour"].astype(str).eq(colour).to_numpy()
                        indicator = mask.astype(float)
                        names.append(f"colour[{colour}]")
                        columns.append(indicator)
                        for predictor in baseline:
                            names.append(f"colour[{colour}]:z_{predictor}")
                            columns.append(_masked_z(work, mask, predictor))
                        if conditional:
                            names.append(f"colour[{colour}]:z_selfing_core")
                            columns.append(_masked_z(work, mask, "selfing_core"))
                        dname = f"colour[{colour}]:z_{geography}"
                        names.append(dname)
                        columns.append(_masked_z(work, mask, geography))
                        distance_names.append(dname)
                    coef, fit, covariance = _fit_grouped_binomial_design(
                        work["successes"].to_numpy(float),
                        work["trials"].to_numpy(float),
                        np.column_stack(columns),
                        names,
                        work[cluster].astype(str).to_numpy(),
                    )
                    if coef.empty:
                        rows.append(
                            {
                                "stratum": stratum,
                                "context": context,
                                "support_tier": tier,
                                "threshold": threshold,
                                "model": model,
                                "status": str(fit.get("status", "fit_failed")),
                                "n_retained_colours": len(retained),
                            }
                        )
                        continue
                    indexed = coef.set_index("predictor")
                    beta = np.array([float(indexed.loc[name, "estimate_log_odds"]) for name in distance_names])
                    indices = [names.index(name) for name in distance_names]
                    cov = covariance[np.ix_(indices, indices)]
                    df = int(np.linalg.matrix_rank(cov))
                    wald = float(beta @ np.linalg.pinv(cov) @ beta) if df > 0 else float("nan")
                    p = _chi_square_sf_integer_df(wald, df) if df > 0 else float("nan")
                    rows.append(
                        {
                            "stratum": stratum,
                            "context": context,
                            "support_tier": tier,
                            "threshold": threshold,
                            "model": model,
                            "status": "fit",
                            "n_retained_colours": len(retained),
                            "retained_colours": "|".join(retained),
                            "n_unique_islands": int(work["island_id"].nunique()),
                            "n_clusters": int(fit["n_clusters"]),
                            "wald_chisq": wald,
                            "df": df,
                            "p_value": p,
                            "distance_slopes": ";".join(
                                f"{colour}:{value:.8g}" for colour, value in zip(retained, beta, strict=True)
                            ),
                        }
                    )
    return pd.DataFrame(rows)


@app.command("run")
def run(
    species_axis_csv: Path = typer.Option(..., exists=True),
    status_flora_csv: Path = typer.Option(..., exists=True),
    island_scores_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    pattern_config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    evidence_scope: str = typer.Option("all"),
    stratum: list[str] = typer.Option(["all_observed"], "--stratum"),
) -> None:
    config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    config = dict(config)
    config["strata"] = [str(x) for x in stratum]
    config["colours"] = list(FOCAL_COLOURS)

    species_axis = pd.read_csv(species_axis_csv, dtype=str)
    status_flora = pd.read_csv(status_flora_csv)
    island_scores = pd.read_csv(island_scores_csv)
    covariates = pd.read_csv(covariates_csv)

    counts = build_raw_colour_counts(
        species_axis,
        status_flora,
        evidence_scope=evidence_scope,
        strata=config["strata"],
        colours=config["colours"],
    )
    models = fit_raw_colour_models(counts, island_scores, covariates, config)
    omnibus = fit_within_context_joint_omnibus(counts, island_scores, covariates, config)

    output_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(output_dir / "raw_colour_counts.csv.gz", index=False, compression="gzip")
    models.to_csv(output_dir / "raw_colour_model_results.csv", index=False)
    omnibus.to_csv(output_dir / "raw_colour_joint_omnibus.csv", index=False)
    manifest = {
        "contract": "chapter1_v13_raw_colour_audit_v1",
        "evidence_scope": evidence_scope,
        "strata": config["strata"],
        "colours": config["colours"],
        "multistate_species_policy": "presence_in_each_explicitly_reported_focal_colour",
        "nonfocal_only_species_in_denominator": False,
        "models": [
            "colour ~ isolation + area + climate",
            "colour ~ isolation + selfing_core + area + climate",
        ],
        "interpretation_ceiling": (
            "raw colour-composition response; not attraction intensity and not realized pollinator identity"
        ),
        "causal_mediation_claimed": False,
    }
    (output_dir / "raw_colour_audit_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
