"""Global four-context GloBI source-breadth extension for Chapter 1 H5.

The predictor is the already frozen, outcome-blind GloBI genus breadth table. This
module expands only the island-side test: all observed floras are primary, source-backed
native floras are sensitivities, and all four predeclared analysis regimes are fitted.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml
from scipy import sparse

from island_v2.chapter1_all_data_probability import (
    _bh,
    _chi_square_sf_integer_df,
    _normal_two_sided_p,
)
from island_v2.chapter1_globi_source_breadth import sha256_file
from island_v2.chapter1_globi_source_breadth_v2 import (
    compute_effort_matched_enrichment,
    reference_effort_bins,
)
from island_v2.chapter1_pr138_lineage_representation_bridge import (
    _availability_matrices,
    _genus,
    _source_assignment_matrix,
    broad_source_availability,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _scope_mask(frame: pd.DataFrame, scope: str) -> pd.Series:
    if scope == "all_observed":
        return pd.Series(True, index=frame.index)
    if scope == "all_native":
        return frame["origin_status"].astype(str).eq("native")
    if scope == "native_nonendemic":
        return frame["floristic_status"].astype(str).eq("native_nonendemic")
    raise ValueError(f"unsupported flora scope: {scope}")


def _representation_count_matrix_scope(
    status_flora: pd.DataFrame,
    island_index: dict[str, int],
    genus_index: dict[str, int],
    *,
    scope: str,
) -> sparse.csr_matrix:
    required = {"island_id", "accepted_species", "origin_status", "floristic_status"}
    if missing := required - set(status_flora.columns):
        raise ValueError(f"status flora lacks columns: {sorted(missing)}")
    work = status_flora.loc[
        _scope_mask(status_flora, scope), ["island_id", "accepted_species"]
    ].drop_duplicates()
    work = work.copy()
    work["island_id"] = work["island_id"].astype(str)
    work["genus"] = work["accepted_species"].map(_genus)
    work = work.loc[
        work["island_id"].isin(island_index) & work["genus"].isin(genus_index)
    ].copy()
    grouped = work.groupby(["island_id", "genus"], as_index=False).agg(
        n_island_species=("accepted_species", "nunique")
    )
    rows = grouped["island_id"].map(island_index).to_numpy(int)
    cols = grouped["genus"].map(genus_index).to_numpy(int)
    return sparse.csr_matrix(
        (grouped["n_island_species"].to_numpy(float), (rows, cols)),
        shape=(len(island_index), len(genus_index)),
    )


def _richness_bins(values: np.ndarray) -> np.ndarray:
    bins = np.zeros_like(values, dtype=np.int8)
    bins[(values > 1) & (values <= 2)] = 1
    bins[(values > 2) & (values <= 4)] = 2
    bins[(values > 4) & (values <= 8)] = 3
    bins[values > 8] = 4
    return bins


def build_global_enrichment(
    genus_breadth: pd.DataFrame,
    gift_flora: pd.DataFrame,
    assignments: pd.DataFrame,
    status_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    predictor = config["predictor"]
    source = config["source_matching"]
    scopes = [str(config["flora_scopes"]["primary"]), *map(str, config["flora_scopes"]["sensitivities"])]
    thresholds = sorted(
        {
            int(predictor["primary_min_independent_references"]),
            *[int(x) for x in predictor["sensitivity_min_independent_references"]],
        }
    )
    matchings = [str(source["primary"]), str(source["sensitivity"])]
    source_modes = [str(x) for x in source["source_modes"]]
    minimum = int(source["minimum_represented_genera"])
    metric = str(predictor["metric"])

    required = {"genus", "n_independent_references", metric}
    if missing := required - set(genus_breadth.columns):
        raise ValueError(f"genus breadth lacks columns: {sorted(missing)}")
    islands = sorted(covariates["island_id"].astype(str).unique())
    island_index = {island: index for index, island in enumerate(islands)}
    parts: list[pd.DataFrame] = []

    for threshold in thresholds:
        eligible = genus_breadth.copy()
        eligible["n_independent_references"] = pd.to_numeric(
            eligible["n_independent_references"], errors="coerce"
        )
        eligible[metric] = pd.to_numeric(eligible[metric], errors="coerce")
        eligible = eligible.loc[
            eligible["n_independent_references"].ge(threshold) & eligible[metric].notna()
        ].drop_duplicates("genus")
        genera = sorted(set(eligible["genus"].astype(str)) - {""})
        if not genera:
            continue
        genus_index = {genus: index for index, genus in enumerate(genera)}
        ordered = eligible.set_index("genus").loc[genera]
        positions = ordered[metric].to_numpy(float)
        effort_bins = reference_effort_bins(
            ordered["n_independent_references"].to_numpy(float)
        )

        availability = broad_source_availability(gift_flora, set(genera))
        entities = sorted(
            set(pd.to_numeric(assignments["entity_ID"], errors="coerce").dropna().astype(int))
            | set(pd.to_numeric(availability["entity_ID"], errors="coerce").dropna().astype(int))
        )
        entity_index = {entity: index for index, entity in enumerate(entities)}
        presence, richness = _availability_matrices(availability, entity_index, genus_index)
        counts_by_scope = {
            scope: _representation_count_matrix_scope(
                status_flora, island_index, genus_index, scope=scope
            ).toarray()
            for scope in scopes
        }
        for source_mode in source_modes:
            assignment = _source_assignment_matrix(
                assignments, island_index, entity_index, source_mode=source_mode
            )
            prevalence = (assignment @ presence).toarray().astype(np.int16)
            source_richness = (assignment @ richness).toarray().astype(np.float32)
            richness_bin = _richness_bins(source_richness)
            for scope in scopes:
                counts = counts_by_scope[scope]
                for matching in matchings:
                    rows: list[dict[str, Any]] = []
                    for island_position, island_id in enumerate(islands):
                        result = compute_effort_matched_enrichment(
                            prevalence[island_position],
                            source_richness[island_position],
                            counts[island_position],
                            positions,
                            effort_bins,
                            matching=matching,
                            minimum_represented_genera=minimum,
                        )
                        if result is None:
                            continue
                        rows.append(
                            {
                                "island_id": island_id,
                                "flora_scope": scope,
                                "source_mode": source_mode,
                                "source_matching": matching,
                                "min_independent_references": threshold,
                                "metric": metric,
                                **result,
                            }
                        )
                    if rows:
                        parts.append(pd.DataFrame(rows))
            del prevalence, source_richness, richness_bin
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def _z(series: pd.Series) -> np.ndarray:
    x = pd.to_numeric(series, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def _fit_clustered_ols(
    y: np.ndarray, design: np.ndarray, clusters: np.ndarray
) -> tuple[np.ndarray, np.ndarray, int]:
    beta = np.linalg.pinv(design.T @ design) @ (design.T @ y)
    residual = y - design @ beta
    bread = np.linalg.pinv(design.T @ design)
    labels = np.asarray(clusters).astype(str)
    unique = np.unique(labels)
    meat = np.zeros((design.shape[1], design.shape[1]), dtype=float)
    for label in unique:
        mask = labels == label
        score = design[mask].T @ residual[mask]
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    n = len(y)
    k = design.shape[1]
    g = len(unique)
    if g > 1 and n > k:
        covariance *= (g / (g - 1.0)) * ((n - 1.0) / (n - k))
    return beta, covariance, g


def _prepare_model_data(
    enrichment: pd.DataFrame, covariates: pd.DataFrame, config: dict[str, Any]
) -> pd.DataFrame:
    model = config["model"]
    required = {
        "island_id",
        "analysis_regime",
        str(model["cluster"]),
        str(model["distance"]),
        str(model["area"]),
        *[str(x) for x in model["controls"]],
    }
    if missing := required - set(covariates.columns):
        raise ValueError(f"covariates lack columns: {sorted(missing)}")
    data = enrichment.merge(
        covariates[list(required)].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    numeric = [
        "entry_enrichment",
        str(model["distance"]),
        str(model["area"]),
        *[str(x) for x in model["controls"]],
    ]
    for column in numeric:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data = data.dropna(subset=numeric)
    data = data.loc[
        data["analysis_regime"].astype(str).isin([str(x) for x in config["contexts"]])
        & data[str(model["cluster"])].fillna("").astype(str).ne("")
    ].copy()
    return data


def fit_within_context(data: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    model = config["model"]
    minimum = int(config["support"]["confirmatory_min_islands_per_context"])
    distance = str(model["distance"])
    area = str(model["area"])
    controls = [str(x) for x in model["controls"]]
    cluster = str(model["cluster"])
    rows: list[dict[str, Any]] = []
    group_cols = [
        "flora_scope",
        "source_mode",
        "source_matching",
        "min_independent_references",
        "analysis_regime",
    ]
    for key, part in data.groupby(group_cols, sort=True):
        part = part.drop_duplicates("island_id").copy()
        meta = dict(zip(group_cols, key, strict=True))
        if len(part) < minimum:
            rows.append({**meta, "status": "not_testable", "n_islands": len(part)})
            continue
        zd = _z(part[distance])
        za = _z(part[area])
        columns = [np.ones(len(part), dtype=float), zd, za, zd * za]
        names = ["intercept", "distance", "area", "distance_by_area"]
        for control in controls:
            columns.append(_z(part[control]))
            names.append(control)
        beta, covariance, n_clusters = _fit_clustered_ols(
            part["entry_enrichment"].to_numpy(float),
            np.column_stack(columns),
            part[cluster].to_numpy(str),
        )
        record: dict[str, Any] = {
            **meta,
            "status": "fit",
            "n_islands": int(part["island_id"].nunique()),
            "n_clusters": int(n_clusters),
        }
        for term in ("distance", "distance_by_area"):
            index = names.index(term)
            estimate = float(beta[index])
            stderr = float(math.sqrt(max(float(covariance[index, index]), 0.0)))
            z_value = estimate / stderr if stderr > 0 else float("nan")
            record[f"{term}_estimate"] = estimate
            record[f"{term}_se"] = stderr
            record[f"{term}_p"] = _normal_two_sided_p(z_value)
        rows.append(record)
    result = pd.DataFrame(rows)
    fit = result["status"].eq("fit") if not result.empty else pd.Series(dtype=bool)
    for term in ("distance", "distance_by_area"):
        if not result.empty:
            result.loc[fit, f"{term}_q"] = result.loc[fit].groupby(
                [
                    "flora_scope",
                    "source_matching",
                    "min_independent_references",
                    "analysis_regime",
                ],
                group_keys=False,
            )[f"{term}_p"].transform(_bh)
    return result


def _joint(vector: np.ndarray, covariance: np.ndarray) -> tuple[float, int, float]:
    rank = int(np.linalg.matrix_rank(covariance))
    if rank <= 0:
        return float("nan"), 0, float("nan")
    statistic = float(vector @ np.linalg.pinv(covariance) @ vector)
    return statistic, rank, _chi_square_sf_integer_df(statistic, rank)


def _context_model(
    part: pd.DataFrame,
    contexts: list[str],
    config: dict[str, Any],
) -> dict[str, Any]:
    model = config["model"]
    reference = contexts[0]
    other = contexts[1:]
    cluster = str(model["cluster"])
    zd = _z(part[str(model["distance"])])
    za = _z(part[str(model["area"])])
    z_controls = {name: _z(part[name]) for name in map(str, model["controls"])}
    indicators = {
        context: part["analysis_regime"].astype(str).eq(context).to_numpy(float)
        for context in other
    }
    columns: list[np.ndarray] = [np.ones(len(part), dtype=float)]
    names = ["intercept"]
    for context in other:
        columns.append(indicators[context])
        names.append(f"context[{context}]")
    for name, values in z_controls.items():
        columns.append(values)
        names.append(name)
        for context in other:
            columns.append(values * indicators[context])
            names.append(f"{name}:context[{context}]")
    for name, values in (("area", za), ("distance", zd)):
        columns.append(values)
        names.append(name)
        for context in other:
            columns.append(values * indicators[context])
            names.append(f"{name}:context[{context}]")
    da = zd * za
    columns.append(da)
    names.append("distance_by_area")
    for context in other:
        columns.append(da * indicators[context])
        names.append(f"distance_by_area:context[{context}]")
    beta, covariance, n_clusters = _fit_clustered_ols(
        part["entry_enrichment"].to_numpy(float),
        np.column_stack(columns),
        part[cluster].to_numpy(str),
    )
    distance_indices = [names.index(f"distance:context[{context}]") for context in other]
    area_indices = [
        names.index(f"distance_by_area:context[{context}]") for context in other
    ]
    d_stat, d_df, d_p = _joint(
        beta[distance_indices], covariance[np.ix_(distance_indices, distance_indices)]
    )
    a_stat, a_df, a_p = _joint(
        beta[area_indices], covariance[np.ix_(area_indices, area_indices)]
    )
    return {
        "reference_context": reference,
        "n_islands": int(part["island_id"].nunique()),
        "n_clusters": int(n_clusters),
        "distance_heterogeneity_wald": d_stat,
        "distance_heterogeneity_df": d_df,
        "distance_heterogeneity_p": d_p,
        "area_heterogeneity_wald": a_stat,
        "area_heterogeneity_df": a_df,
        "area_heterogeneity_p": a_p,
    }


def fit_global_context_tests(data: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    contexts = [str(x) for x in config["contexts"]]
    minimum = int(config["support"]["confirmatory_min_islands_per_context"])
    rows: list[dict[str, Any]] = []
    group_cols = ["flora_scope", "source_mode", "source_matching", "min_independent_references"]
    for key, part in data.groupby(group_cols, sort=True):
        part = part.drop_duplicates(["island_id", "analysis_regime"]).copy()
        support = part.groupby("analysis_regime")["island_id"].nunique()
        meta = dict(zip(group_cols, key, strict=True))
        if any(int(support.get(context, 0)) < minimum for context in contexts):
            rows.append({**meta, "status": "not_testable"})
            continue
        rows.append({**meta, "status": "fit", **_context_model(part, contexts, config)})
    result = pd.DataFrame(rows)
    fit = result["status"].eq("fit") if not result.empty else pd.Series(dtype=bool)
    for family in ("distance_heterogeneity", "area_heterogeneity"):
        if not result.empty:
            result.loc[fit, f"{family}_q"] = result.loc[fit].groupby(
                ["flora_scope", "source_matching", "min_independent_references"],
                group_keys=False,
            )[f"{family}_p"].transform(_bh)
    return result


def fit_primary_pair(data: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    contexts = [str(x) for x in config["primary_pair"]]
    subset = data.loc[data["analysis_regime"].astype(str).isin(contexts)].copy()
    result = fit_global_context_tests(subset, {**config, "contexts": contexts})
    if result.empty:
        return result
    return result.rename(
        columns={
            "distance_heterogeneity_wald": "pair_distance_wald",
            "distance_heterogeneity_df": "pair_distance_df",
            "distance_heterogeneity_p": "pair_distance_p",
            "distance_heterogeneity_q": "pair_distance_q",
            "area_heterogeneity_wald": "pair_area_wald",
            "area_heterogeneity_df": "pair_area_df",
            "area_heterogeneity_p": "pair_area_p",
            "area_heterogeneity_q": "pair_area_q",
        }
    )


def classify_primary(global_tests: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    primary_scope = str(config["flora_scopes"]["primary"])
    primary_matching = str(config["source_matching"]["primary"])
    threshold = int(config["predictor"]["primary_min_independent_references"])
    source_modes = [str(x) for x in config["source_matching"]["source_modes"]]
    alpha = float(config["support"]["alpha"])
    target = global_tests.loc[
        global_tests["flora_scope"].astype(str).eq(primary_scope)
        & global_tests["source_matching"].astype(str).eq(primary_matching)
        & global_tests["min_independent_references"].eq(threshold)
        & global_tests["status"].eq("fit")
    ].copy()
    rows = []
    for family in ("distance_heterogeneity", "area_heterogeneity"):
        supported = target.loc[target[f"{family}_q"].le(alpha), "source_mode"].astype(str)
        n_supported = int(supported.nunique())
        complete = set(target["source_mode"].astype(str)) == set(source_modes)
        if complete and n_supported == len(source_modes):
            classification = "robust_global_context_heterogeneity"
        elif n_supported > 0:
            classification = "source_definition_sensitive_context_heterogeneity"
        else:
            classification = "global_context_heterogeneity_not_supported"
        rows.append(
            {
                "family": family,
                "flora_scope": primary_scope,
                "n_source_modes_expected": len(source_modes),
                "n_source_modes_supported": n_supported,
                "classification": classification,
            }
        )
    return pd.DataFrame(rows)


def run_global(
    genus_breadth: pd.DataFrame,
    gift_flora: pd.DataFrame,
    assignments: pd.DataFrame,
    status_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> dict[str, pd.DataFrame]:
    enrichment = build_global_enrichment(
        genus_breadth, gift_flora, assignments, status_flora, covariates, config
    )
    data = _prepare_model_data(enrichment, covariates, config)
    within = fit_within_context(data, config)
    global_tests = fit_global_context_tests(data, config)
    pair = fit_primary_pair(data, config)
    classification = classify_primary(global_tests, config)
    return {
        "enrichment": enrichment,
        "within": within,
        "global": global_tests,
        "pair": pair,
        "classification": classification,
    }


@app.command("run")
def run(
    genus_breadth_csv: Path = typer.Option(..., exists=True),
    gift_flora_csv: Path = typer.Option(..., exists=True),
    assignments_csv: Path = typer.Option(..., exists=True),
    status_flora_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    predictor_sha256: str = typer.Option(""),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if config.get("contract") != "chapter1_globi_global_v3":
        raise typer.BadParameter("unexpected global GloBI contract")
    observed_sha = sha256_file(genus_breadth_csv)
    if predictor_sha256 and observed_sha != predictor_sha256:
        raise typer.BadParameter("frozen GloBI predictor SHA mismatch")
    outputs = run_global(
        pd.read_csv(genus_breadth_csv),
        pd.read_csv(gift_flora_csv),
        pd.read_csv(assignments_csv),
        pd.read_csv(status_flora_csv),
        pd.read_csv(covariates_csv),
        config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs["enrichment"].to_csv(
        output_dir / "globi_global_island_enrichment.csv.gz",
        index=False,
        compression="gzip",
    )
    outputs["within"].to_csv(output_dir / "globi_global_within_context.csv", index=False)
    outputs["global"].to_csv(output_dir / "globi_global_context_tests.csv", index=False)
    outputs["pair"].to_csv(output_dir / "globi_global_north_tropical.csv", index=False)
    outputs["classification"].to_csv(
        output_dir / "globi_global_primary_classification.csv", index=False
    )
    manifest = {
        "contract": config["contract"],
        "predictor_sha256": observed_sha,
        "n_enrichment_rows": int(len(outputs["enrichment"])),
        "n_primary_all_observed_islands": int(
            outputs["enrichment"].loc[
                outputs["enrichment"]["flora_scope"].eq("all_observed"), "island_id"
            ].nunique()
        ),
        "contexts": [str(x) for x in config["contexts"]],
        "claim_boundary": config["claim_ceiling"],
    }
    (output_dir / "globi_global_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    lines = [
        "# Global GloBI H5 extension",
        "",
        f"- all-observed islands with GloBI enrichment: {manifest['n_primary_all_observed_islands']}",
        "- four contexts tested: " + ", ".join(manifest["contexts"]),
        "",
        "## Primary classifications",
        "",
        outputs["classification"].to_csv(index=False),
    ]
    (output_dir / "RESULT_SUMMARY.md").write_text("\n".join(lines), encoding="utf-8")
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
