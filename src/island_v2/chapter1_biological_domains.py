"""Post-hoc biological-domain sensitivity for Chapter 1.

This analysis never replaces the confirmatory final H1.  It re-expresses the frozen
component-trait evidence in three biologically interpretable domains:

A. reproductive assurance (frozen PR138 selfing_core);
B. accessibility/specialization (frozen PR138 generalized_accessible);
C. signal/display (colour, flower size, inflorescence display).

Domain C is not forced onto one attraction-intensity axis.  Its only directional
sensitivity is reduced flower size, reusing the pre-existing PR138 selfing-syndrome
flower-size recode. Colour hue and inflorescence type are analysed state by state.
"""
from __future__ import annotations

import json
import math
from functools import reduce
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml
from scipy.stats import t as student_t

from island_v2.chapter1_all_data_probability import (
    _assemble_cluster_covariance,
    _fit_single_beta_binomial,
    _standardize,
)
from island_v2.chapter1_final_inference_audit import (
    finite_cluster_one_sided_p,
    finite_cluster_two_sided_p,
)
from island_v2.chapter1_h1_final_directional import synthesize_regions
from island_v2.chapter1_h2_complete_selfing_sensitivity import (
    aggregate_island_score,
    build_species_score,
)
from island_v2.chapter1_h5_glopl_global_distance import (
    MEASUREMENT_COLUMNS,
    _site_key,
    assign_context,
    build_study_key,
)
from island_v2.chapter1_v14_h2_decomposition import _bh, _clustered_ols
from island_v2.chapter1_v14_h4_family_bridge import fit_global_family_score

app = typer.Typer(add_completion=False, no_args_is_help=True)
CONTRACT = "chapter1_biological_domains_v1"


def load_config(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict) or value.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected biological-domain contract")
    return value


def _truthy(series: pd.Series) -> pd.Series:
    return (
        series.fillna(False)
        .astype(str)
        .str.strip()
        .str.casefold()
        .isin({"true", "1", "yes", "y"})
    )


def _name_key(value: object) -> str:
    return " ".join(str(value).replace("_", " ").split()).casefold()


def domain_specs(
    syndrome_config: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    self_spec = syndrome_config["syndromes"]["selfing_core"]
    access_spec = syndrome_config["syndromes"]["generalized_accessible"]
    size_parent = syndrome_config["syndromes"]["selfing_syndrome"]["traits"][
        "flower_size_class"
    ]
    size_spec = {
        "role": "reduced_flower_size_directional_sensitivity",
        "traits": {
            "flower_size_class": {
                "weight": 1.0,
                "preferred": list(size_parent["preferred"]),
                "opposed": list(size_parent["opposed"]),
            }
        },
    }
    expected = config["domains"]["signal_display"]["directional_sensitivity"]
    if set(size_spec["traits"]["flower_size_class"]["preferred"]) != set(
        expected["preferred"]
    ):
        raise typer.BadParameter("flower-size preferred states drifted from contract")
    if set(size_spec["traits"]["flower_size_class"]["opposed"]) != set(
        expected["opposed"]
    ):
        raise typer.BadParameter("flower-size opposed states drifted from contract")
    return {
        "reproductive_assurance": self_spec,
        "accessibility_specialization": access_spec,
        "reduced_flower_size": size_spec,
    }


def build_domain_species_scores(
    state_audit: pd.DataFrame,
    syndrome_config: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, pd.DataFrame]:
    specs = domain_specs(syndrome_config, config)
    mins = config["species_score"]["minimum_informative_traits"]
    out: dict[str, pd.DataFrame] = {}
    for name, spec in specs.items():
        frame = build_species_score(
            state_audit,
            spec,
            minimum_informative_traits=int(mins[name]),
            require_all=False,
        ).rename(
            columns={
                "score": name,
                "n_informative_traits": f"{name}_n_informative_traits",
            }
        )
        out[name] = frame
    return out


def aggregate_domain_scores(
    flora: pd.DataFrame,
    species_scores: dict[str, pd.DataFrame],
    *,
    minimum_species: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    island_parts: list[pd.DataFrame] = []
    species_wide: pd.DataFrame | None = None
    for name, frame in species_scores.items():
        source = frame[["accepted_species", name]].rename(
            columns={name: "score"}
        )
        island = aggregate_island_score(
            flora,
            source,
            score_name=name,
        )
        island = island.loc[
            island[f"{name}_n_species"].ge(int(minimum_species))
        ].copy()
        island_parts.append(island)
        score_only = frame[["accepted_species", name]].dropna(subset=[name])
        species_wide = (
            score_only
            if species_wide is None
            else species_wide.merge(
                score_only,
                on="accepted_species",
                how="outer",
                validate="one_to_one",
            )
        )

    island_wide = reduce(
        lambda left, right: left.merge(right, on="island_id", how="outer"),
        island_parts,
    )
    assert species_wide is not None
    return island_wide, species_wide


def aggregate_common_species_scores(
    flora: pd.DataFrame,
    species_wide: pd.DataFrame,
    *,
    minimum_species: int,
) -> pd.DataFrame:
    domains = [
        "reproductive_assurance",
        "accessibility_specialization",
        "reduced_flower_size",
    ]
    common = species_wide.dropna(subset=domains).copy()
    base = flora[["island_id", "accepted_species"]].copy()
    base["island_id"] = base["island_id"].astype(str)
    base["accepted_species"] = base["accepted_species"].astype(str)
    joined = base.merge(
        common[["accepted_species", *domains]],
        on="accepted_species",
        how="inner",
        validate="many_to_many",
    ).drop_duplicates(["island_id", "accepted_species"])
    if joined.empty:
        return pd.DataFrame(columns=["island_id", *domains, "common_n_species"])
    result = (
        joined.groupby("island_id", as_index=False)
        .agg(
            reproductive_assurance=("reproductive_assurance", "mean"),
            accessibility_specialization=("accessibility_specialization", "mean"),
            reduced_flower_size=("reduced_flower_size", "mean"),
            common_n_species=("accepted_species", "nunique"),
        )
    )
    return result.loc[
        result["common_n_species"].ge(int(minimum_species))
    ].copy()


def _fit_domain_regional(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
    *,
    evidence_scope: str,
    support_scheme: str,
) -> pd.DataFrame:
    geography = str(config["geography_column"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    needed = ["island_id", geography, context, cluster, *baseline]
    work = island_scores.merge(
        covariates[needed].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="one_to_one",
    )
    rows: list[dict[str, Any]] = []
    domains = [
        "reproductive_assurance",
        "accessibility_specialization",
        "reduced_flower_size",
    ]
    for context_value in [str(x) for x in config["contexts"]]:
        part = work.loc[work[context].astype(str).eq(context_value)].copy()
        for domain in domains:
            result = _clustered_ols(
                part,
                response=domain,
                predictors=[geography, *baseline],
                cluster_column=cluster,
            )
            row: dict[str, Any] = {
                "evidence_scope": evidence_scope,
                "support_scheme": support_scheme,
                "context": context_value,
                "domain": domain,
                "status": result["status"],
                "n_islands": int(result.get("n_islands", 0)),
                "n_clusters": int(result.get("n_clusters", 0)),
            }
            if result["status"] == "fit":
                coef = result["coefficients"][f"z_{geography}"]
                estimate = float(coef["estimate"])
                se = float(coef["se"])
                g = int(result["n_clusters"])
                p_one = finite_cluster_one_sided_p(
                    estimate,
                    se,
                    g,
                    alternative="positive",
                )
                row.update(
                    estimate=estimate,
                    cluster_robust_se=se,
                    finite_cluster_df=g - 1,
                    p_one_sided_positive=p_one,
                    p_two_sided=finite_cluster_two_sided_p(
                        estimate,
                        se,
                        g,
                    ),
                    positive_direction=bool(estimate > 0),
                    supported_positive=bool(
                        estimate > 0 and p_one <= float(config["alpha"])
                    ),
                )
            rows.append(row)
    return pd.DataFrame(rows)


def domain_global_synthesis(
    regional: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    contexts = [str(x) for x in config["contexts"]]
    for (scope, scheme, domain), part in regional.groupby(
        ["evidence_scope", "support_scheme", "domain"],
        dropna=False,
    ):
        fitted = part.loc[part["status"].eq("fit")].copy()
        input_frame = pd.DataFrame(
            {
                "status": fitted["status"],
                "context": fitted["context"],
                "estimate": fitted["estimate"],
                "cluster_robust_se": fitted["cluster_robust_se"],
            }
        )
        summary = synthesize_regions(
            input_frame,
            contexts=contexts,
            alpha=float(config["alpha"]),
        )
        summary.update(
            evidence_scope=str(scope),
            support_scheme=str(scheme),
            domain=str(domain),
        )
        if summary.get("status") == "fit":
            strict = fitted.set_index("context").reindex(contexts)
            if strict["p_one_sided_positive"].notna().all():
                summary["strict_iut_p"] = float(
                    strict["p_one_sided_positive"].max()
                )
                summary["strict_four_region_supported"] = bool(
                    strict["positive_direction"].all()
                    and summary["strict_iut_p"] <= float(config["alpha"])
                )
        rows.append(summary)
    return pd.DataFrame(rows)


def _canonical_tokens(value: object) -> set[str]:
    return {
        x.strip()
        for x in str(value or "").split("|")
        if x.strip() and x.strip() != "unresolved"
    }


def signal_state_ledger(
    state_audit: pd.DataFrame,
) -> pd.DataFrame:
    traits = {
        "flower_primary_color",
        "flower_size_class",
        "inflorescence_display",
    }
    work = state_audit.loc[
        _truthy(state_audit["resolved_for_primary"])
        & state_audit["trait_name"].astype(str).isin(traits)
    ].copy()
    rows: list[dict[str, str]] = []
    for row in work[
        ["accepted_species", "trait_name", "canonical_signature"]
    ].drop_duplicates().itertuples(index=False):
        for state in sorted(_canonical_tokens(row.canonical_signature)):
            rows.append(
                {
                    "accepted_species": str(row.accepted_species),
                    "trait_name": str(row.trait_name),
                    "state": state,
                }
            )
    return pd.DataFrame(rows).drop_duplicates()


def signal_state_counts(
    flora: pd.DataFrame,
    ledger: pd.DataFrame,
) -> pd.DataFrame:
    base = flora[["island_id", "accepted_species"]].drop_duplicates().copy()
    base["island_id"] = base["island_id"].astype(str)
    base["accepted_species"] = base["accepted_species"].astype(str)

    trait_species = ledger[
        ["accepted_species", "trait_name"]
    ].drop_duplicates()
    state_species = ledger[
        ["accepted_species", "trait_name", "state"]
    ].drop_duplicates()
    trait_join = base.merge(
        trait_species,
        on="accepted_species",
        how="inner",
        validate="many_to_many",
    )
    denom = (
        trait_join.groupby(["island_id", "trait_name"], as_index=False)
        .agg(trials=("accepted_species", "nunique"))
    )
    state_join = base.merge(
        state_species,
        on="accepted_species",
        how="inner",
        validate="many_to_many",
    )
    succ = (
        state_join.groupby(
            ["island_id", "trait_name", "state"],
            as_index=False,
        )
        .agg(successes=("accepted_species", "nunique"))
    )
    state_list = ledger[
        ["trait_name", "state"]
    ].drop_duplicates()
    parts: list[pd.DataFrame] = []
    for trait, part in denom.groupby("trait_name", sort=False):
        states = state_list.loc[
            state_list["trait_name"].eq(trait),
            "state",
        ].tolist()
        expanded = part.assign(_k=1).merge(
            pd.DataFrame({"state": states, "_k": 1}),
            on="_k",
            how="inner",
        ).drop(columns="_k")
        expanded = expanded.merge(
            succ.loc[
                succ["trait_name"].eq(trait),
                ["island_id", "state", "successes"],
            ],
            on=["island_id", "state"],
            how="left",
            validate="one_to_one",
        )
        expanded["successes"] = expanded["successes"].fillna(0).astype(int)
        parts.append(expanded)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def fit_signal_states(
    counts: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
    *,
    evidence_scope: str,
) -> pd.DataFrame:
    geography = str(config["geography_column"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    needed = ["island_id", geography, context, cluster, *baseline]
    data = counts.merge(
        covariates[needed].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    for column in ["successes", "trials", geography, *baseline]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data = data.dropna(
        subset=["successes", "trials", geography, *baseline, cluster]
    )
    rows: list[dict[str, Any]] = []
    for context_value in [str(x) for x in config["contexts"]]:
        regional = data.loc[
            data[context].astype(str).eq(context_value)
        ].copy()
        for (trait, state), part in regional.groupby(
            ["trait_name", "state"],
            sort=False,
        ):
            n_islands = int(part["island_id"].nunique())
            positive = int(part.loc[part["successes"].gt(0), "island_id"].nunique())
            negative = int(
                part.loc[
                    part["successes"].lt(part["trials"]),
                    "island_id",
                ].nunique()
            )
            row: dict[str, Any] = {
                "evidence_scope": evidence_scope,
                "context": context_value,
                "trait_name": trait,
                "state": state,
                "n_islands": n_islands,
                "positive_islands": positive,
                "negative_islands": negative,
            }
            if (
                n_islands < int(config["minimum_islands"])
                or positive < 10
                or negative < 10
            ):
                row["status"] = "not_testable"
                rows.append(row)
                continue
            columns = [np.ones(len(part), dtype=float)]
            names = ["intercept"]
            for predictor in baseline:
                columns.append(_standardize(part[predictor]))
                names.append(f"z_{predictor}")
            columns.append(_standardize(part[geography]))
            names.append(f"z_{geography}")
            fit = _fit_single_beta_binomial(
                part["successes"].to_numpy(float),
                part["trials"].to_numpy(float),
                np.column_stack(columns),
                names,
                max_iter=1500,
            )
            cov, assembled, theta = _assemble_cluster_covariance(
                [fit],
                [part[cluster].astype(str).to_numpy()],
            )
            idx = assembled.index(f"z_{geography}")
            estimate = float(theta[idx])
            se = float(math.sqrt(max(float(cov[idx, idx]), 0.0)))
            g = int(part[cluster].nunique())
            row.update(
                status="fit",
                estimate=estimate,
                cluster_robust_se=se,
                n_clusters=g,
                finite_cluster_df=g - 1,
                p_two_sided=finite_cluster_two_sided_p(
                    estimate,
                    se,
                    g,
                ),
                optimizer_success=bool(fit["success"]),
            )
            rows.append(row)
    frame = pd.DataFrame(rows)
    frame["q_within_trait_region"] = np.nan
    fit_mask = frame["status"].eq("fit")
    for _, idx in frame.loc[fit_mask].groupby(
        ["evidence_scope", "context", "trait_name"]
    ).groups.items():
        frame.loc[idx, "q_within_trait_region"] = _bh(
            frame.loc[idx, "p_two_sided"]
        )
    return frame


def fit_h2_models(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
    *,
    evidence_scope: str,
    support_scheme: str,
) -> pd.DataFrame:
    geography = str(config["geography_column"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    needed = ["island_id", geography, context, cluster, *baseline]
    work = island_scores.merge(
        covariates[needed].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="one_to_one",
    )
    models = [
        (
            "accessibility_given_assurance",
            "accessibility_specialization",
            [geography, "reproductive_assurance", *baseline],
            True,
        ),
        (
            "reduced_size_given_assurance",
            "reduced_flower_size",
            [geography, "reproductive_assurance", *baseline],
            True,
        ),
        (
            "accessibility_given_assurance_and_size",
            "accessibility_specialization",
            [
                geography,
                "reproductive_assurance",
                "reduced_flower_size",
                *baseline,
            ],
            False,
        ),
        (
            "reduced_size_given_assurance_and_accessibility",
            "reduced_flower_size",
            [
                geography,
                "reproductive_assurance",
                "accessibility_specialization",
                *baseline,
            ],
            False,
        ),
    ]
    rows: list[dict[str, Any]] = []
    for context_value in [str(x) for x in config["contexts"]]:
        part = work.loc[work[context].astype(str).eq(context_value)].copy()
        for model_name, response, predictors, primary in models:
            result = _clustered_ols(
                part,
                response=response,
                predictors=predictors,
                cluster_column=cluster,
            )
            row: dict[str, Any] = {
                "evidence_scope": evidence_scope,
                "support_scheme": support_scheme,
                "context": context_value,
                "model": model_name,
                "response": response,
                "primary_posthoc": primary,
                "status": result["status"],
                "n_islands": int(result.get("n_islands", 0)),
                "n_clusters": int(result.get("n_clusters", 0)),
            }
            if result["status"] == "fit":
                coef = result["coefficients"][f"z_{geography}"]
                estimate = float(coef["estimate"])
                se = float(coef["se"])
                g = int(result["n_clusters"])
                row.update(
                    isolation_estimate=estimate,
                    isolation_se=se,
                    isolation_p_two_sided=finite_cluster_two_sided_p(
                        estimate,
                        se,
                        g,
                    ),
                    isolation_p_one_sided_positive=finite_cluster_one_sided_p(
                        estimate,
                        se,
                        g,
                        alternative="positive",
                    ),
                )
            rows.append(row)
    frame = pd.DataFrame(rows)
    frame["primary_posthoc_q"] = np.nan
    primary = frame["status"].eq("fit") & frame["primary_posthoc"].astype(bool)
    if bool(primary.any()):
        frame.loc[primary, "primary_posthoc_q"] = _bh(
            frame.loc[primary, "isolation_p_two_sided"]
        )
    return frame


def audit_h3(
    main_payload: dict[str, Any],
    offshore_payload: dict[str, Any],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    primary = main_payload["corrected"]["global_gradient"]
    for name, item in (
        ("corrected_global", primary),
        (
            "corrected_supplemental_only",
            main_payload["corrected"]["sensitivities"]["supplemental_only"],
        ),
        (
            "corrected_no_zero_constant",
            main_payload["corrected"]["sensitivities"]["no_zero_constant"],
        ),
        (
            "offshore_only",
            offshore_payload["offshore_continuous_gradient"],
        ),
    ):
        estimate = float(
            item.get("distance_slope", item.get("estimate"))
        )
        se = float(item.get("distance_slope_se", item.get("se")))
        g = int(item["n_publications"])
        rows.append(
            {
                "analysis": name,
                "estimate": estimate,
                "se": se,
                "n_publications": g,
                "finite_publication_p_two_sided": finite_cluster_two_sided_p(
                    estimate,
                    se,
                    g,
                ),
                "finite_publication_p_one_sided_positive":
                    finite_cluster_one_sided_p(
                        estimate,
                        se,
                        g,
                        alternative="positive",
                    ),
            }
        )
    return pd.DataFrame(rows)


def prepare_corrected_glopl(
    glopl: pd.DataFrame,
    corrected_sites: pd.DataFrame,
    parent_config: dict[str, Any],
    config: dict[str, Any],
) -> pd.DataFrame:
    needed = {
        "Latitude",
        "Longitude",
        "DOI",
        "Author",
        "Year",
        "Species_accepted_names",
        "PL_Effect_Size",
        *MEASUREMENT_COLUMNS,
    }
    if missing := needed - set(glopl.columns):
        raise typer.BadParameter(f"GloPL missing columns: {sorted(missing)}")
    work = glopl.copy()
    work["latitude"] = pd.to_numeric(work["Latitude"], errors="coerce")
    work["longitude"] = pd.to_numeric(work["Longitude"], errors="coerce")
    work["site_key"] = _site_key(work["latitude"], work["longitude"])
    work["study_key"] = build_study_key(work)
    work["analysis_regime"] = assign_context(
        work["latitude"],
        parent_config,
    )
    sites = corrected_sites[
        ["site_key", "spherical_distance_km"]
    ].drop_duplicates("site_key")
    work = work.merge(
        sites,
        on="site_key",
        how="left",
        validate="many_to_one",
    )
    work["PL_Effect_Size"] = pd.to_numeric(
        work["PL_Effect_Size"],
        errors="coerce",
    )
    work["log1p_distance"] = np.log1p(
        pd.to_numeric(
            work["spherical_distance_km"],
            errors="coerce",
        )
    )
    mean = float(
        config["H4"]["corrected_GloPL_distance"]["mean_log1p_distance"]
    )
    sd = float(
        config["H4"]["corrected_GloPL_distance"]["sd_log1p_distance"]
    )
    work["z_distance"] = (work["log1p_distance"] - mean) / sd
    work["species_key"] = work["Species_accepted_names"].map(_name_key)
    for column in MEASUREMENT_COLUMNS:
        work[column] = work[column].fillna("").astype(str).str.strip()
    valid = (
        np.isfinite(work["PL_Effect_Size"].to_numpy(float))
        & np.isfinite(work["z_distance"].to_numpy(float))
        & work["site_key"].astype(str).ne("")
        & work["study_key"].astype(str).ne("")
        & work["species_key"].astype(str).ne("")
        & work["analysis_regime"].isin(config["contexts"])
    )
    return work.loc[valid].copy()


def fit_h4_domains(
    glopl_rows: pd.DataFrame,
    species_scores: dict[str, pd.DataFrame],
    config: dict[str, Any],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for domain in [
        "reproductive_assurance",
        "accessibility_specialization",
        "reduced_flower_size",
    ]:
        scores = species_scores[domain][
            ["accepted_species", domain]
        ].dropna(subset=[domain]).copy()
        # Existing Chapter 1 H4 uses H2 soft membership on [0, 1], whereas
        # plant-side syndrome scores are stored as raw concordance on [-1, 1].
        # Reuse the frozen transform so reproductive-assurance and accessibility
        # must reproduce the existing exact-score bridge before interpreting
        # the newly added reduced-flower-size bridge.
        scores[domain] = (pd.to_numeric(scores[domain], errors="coerce") + 1.0) / 2.0
        scores["species_key"] = scores["accepted_species"].map(_name_key)
        scores = scores[
            ["species_key", domain]
        ].drop_duplicates("species_key")
        joined = glopl_rows.merge(
            scores,
            on="species_key",
            how="inner",
            validate="many_to_one",
        )
        group_cols = [
            "study_key",
            "site_key",
            "species_key",
            "analysis_regime",
            "z_distance",
            domain,
            *MEASUREMENT_COLUMNS,
        ]
        cells = (
            joined.groupby(
                group_cols,
                as_index=False,
                dropna=False,
            )
            .agg(
                PL_Effect_Size=("PL_Effect_Size", "mean"),
                n_effect_rows=("PL_Effect_Size", "size"),
            )
            .reset_index(drop=True)
        )
        counts = cells.groupby("study_key")["site_key"].transform("size")
        cells["analysis_weight"] = (
            float(config["H4"]["publication_total_weight"])
            / counts.astype(float)
        )
        for analysis, part in (
            ("primary", cells),
            (
                "supplemental_only",
                cells.loc[
                    cells["PL_Effect_Size_Type2"].astype(str).eq("Sup")
                ].copy(),
            ),
            (
                "no_zero_constant",
                cells.loc[
                    ~cells["Constant_added"]
                    .astype(str)
                    .str.casefold()
                    .isin({"true", "1", "yes", "y", "t"})
                ].copy(),
            ),
        ):
            result = fit_global_family_score(
                part,
                score_name=domain,
            )
            row: dict[str, Any] = {
                "domain": domain,
                "analysis": analysis,
                "evaluable": bool(result.get("evaluable", False)),
                "n_cells": int(result.get("n_cells", 0)),
                "n_publications": int(result.get("n_publications", 0)),
                "n_sites": int(result.get("n_sites", 0)),
                "n_species": int(result.get("n_species", 0)),
            }
            if result.get("evaluable"):
                estimate = float(result["estimate"])
                se = float(result["se"])
                g = int(result["n_publications"])
                row.update(
                    estimate=estimate,
                    se=se,
                    finite_publication_p_two_sided=
                        finite_cluster_two_sided_p(
                            estimate,
                            se,
                            g,
                        ),
                    finite_publication_p_one_sided_negative=
                        finite_cluster_one_sided_p(
                            estimate,
                            se,
                            g,
                            alternative="negative",
                        ),
                )
            rows.append(row)
    return pd.DataFrame(rows)


def run_analysis(
    flora: pd.DataFrame,
    state_audits: dict[str, pd.DataFrame],
    covariates: pd.DataFrame,
    syndrome_config: dict[str, Any],
    config: dict[str, Any],
    glopl: pd.DataFrame,
    corrected_glopl_sites: pd.DataFrame,
    parent_glopl_config: dict[str, Any],
    h3_main: dict[str, Any],
    h3_offshore: dict[str, Any],
) -> dict[str, Any]:
    minimum_species = int(
        config["species_score"]["minimum_scored_species_per_island"]
    )
    h1_parts: list[pd.DataFrame] = []
    h1_global_parts: list[pd.DataFrame] = []
    h2_parts: list[pd.DataFrame] = []
    state_parts: list[pd.DataFrame] = []
    score_support_rows: list[dict[str, Any]] = []
    direct_species_scores: dict[str, pd.DataFrame] | None = None

    for scope, audit in state_audits.items():
        species_scores = build_domain_species_scores(
            audit,
            syndrome_config,
            config,
        )
        if scope == "direct_only":
            direct_species_scores = species_scores
        island_scores, species_wide = aggregate_domain_scores(
            flora,
            species_scores,
            minimum_species=minimum_species,
        )
        common_scores = aggregate_common_species_scores(
            flora,
            species_wide,
            minimum_species=minimum_species,
        )

        for domain, frame in species_scores.items():
            score_support_rows.append(
                {
                    "evidence_scope": scope,
                    "domain": domain,
                    "n_species_scored": int(frame[domain].notna().sum()),
                }
            )
        score_support_rows.append(
            {
                "evidence_scope": scope,
                "domain": "common_three_domain_species",
                "n_species_scored": int(
                    species_wide.dropna(
                        subset=[
                            "reproductive_assurance",
                            "accessibility_specialization",
                            "reduced_flower_size",
                        ]
                    )["accepted_species"].nunique()
                ),
            }
        )

        for scheme, scores in (
            ("available_domain_support", island_scores),
            ("common_species_support", common_scores),
        ):
            h1 = _fit_domain_regional(
                scores,
                covariates,
                config,
                evidence_scope=scope,
                support_scheme=scheme,
            )
            h1_parts.append(h1)
            h1_global_parts.append(
                domain_global_synthesis(h1, config)
            )
            h2_parts.append(
                fit_h2_models(
                    scores,
                    covariates,
                    config,
                    evidence_scope=scope,
                    support_scheme=scheme,
                )
            )

        ledger = signal_state_ledger(audit)
        counts = signal_state_counts(flora, ledger)
        states = fit_signal_states(
            counts,
            covariates,
            config,
            evidence_scope=scope,
        )
        state_parts.append(states)

    if direct_species_scores is None:
        raise typer.BadParameter("direct_only evidence scope is required for H4")

    corrected_rows = prepare_corrected_glopl(
        glopl,
        corrected_glopl_sites,
        parent_glopl_config,
        config,
    )
    h4 = fit_h4_domains(
        corrected_rows,
        direct_species_scores,
        config,
    )
    h3 = audit_h3(h3_main, h3_offshore)

    return {
        "h1_regional": pd.concat(h1_parts, ignore_index=True),
        "h1_global": pd.concat(h1_global_parts, ignore_index=True),
        "signal_states": pd.concat(state_parts, ignore_index=True),
        "h2": pd.concat(h2_parts, ignore_index=True),
        "h3": h3,
        "h4": h4,
        "score_support": pd.DataFrame(score_support_rows),
    }


@app.command("run")
def run(
    status_flora_csv: Path = typer.Option(..., exists=True),
    state_audit_all_csv: Path = typer.Option(..., exists=True),
    state_audit_direct_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    syndrome_config_path: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    glopl_csv: Path = typer.Option(..., exists=True),
    corrected_glopl_sites_csv: Path = typer.Option(..., exists=True),
    parent_glopl_config_path: Path = typer.Option(..., exists=True),
    h3_json: Path = typer.Option(..., exists=True),
    h3_offshore_json: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    syndrome_config = yaml.safe_load(
        syndrome_config_path.read_text(encoding="utf-8")
    )
    parent_glopl_config = yaml.safe_load(
        parent_glopl_config_path.read_text(encoding="utf-8")
    )
    result = run_analysis(
        pd.read_csv(status_flora_csv),
        {
            "all_analysis_eligible": pd.read_csv(
                state_audit_all_csv
            ),
            "direct_only": pd.read_csv(
                state_audit_direct_csv
            ),
        },
        pd.read_csv(covariates_csv),
        syndrome_config,
        config,
        pd.read_csv(glopl_csv, encoding="latin-1"),
        pd.read_csv(corrected_glopl_sites_csv),
        parent_glopl_config,
        json.loads(h3_json.read_text(encoding="utf-8")),
        json.loads(h3_offshore_json.read_text(encoding="utf-8")),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, frame in result.items():
        assert isinstance(frame, pd.DataFrame)
        frame.to_csv(output_dir / f"{name}.csv", index=False)

    summary = {
        "contract": CONTRACT,
        "inferential_role": config["inferential_role"],
        "H3_changed_by_reclassification": False,
        "confirmatory_H1_replaced": False,
        "historical_selection_identified": False,
        "mediation_identified": False,
        "named_pollinator_identity_identified": False,
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary, indent=2))


if __name__ == "__main__":
    app()
