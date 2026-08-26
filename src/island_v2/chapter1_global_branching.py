"""Global when/where test of plant-side branches predicted by the central hypothesis.

The central hypothesis is the Dispersal-filtered pollination-service branching
hypothesis. Chapter 1 does not observe pollination service, however, so this module
tests only its predeclared plant-side response geometry across broad analysis regimes
and formal biogeographic realms.

Branch scores remain continuous multivariate floral-trait concordance responses. They
must never be interpreted as observations of pollinator identity, loss, mobility, or
functional replacement. Those direct mechanism tests belong to ``izu-core``.
"""

from __future__ import annotations

import itertools
import json
from copy import deepcopy
from pathlib import Path
from typing import Annotated, Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_pr138_syndrome_analysis import (
    _between_contexts,
    _bh,
    _prepare,
    _within_context,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def build_branch_scores(
    island_scores: pd.DataFrame,
    branching_config: dict[str, Any],
) -> pd.DataFrame:
    """Build fail-closed island branch scores from predeclared syndrome components."""

    required = {
        "island_id",
        "stratum",
        "syndrome",
        "syndrome_score",
        "n_species",
    }
    missing = required - set(island_scores.columns)
    if missing:
        raise typer.BadParameter(f"island syndrome scores missing columns: {sorted(missing)}")

    work = island_scores[list(required)].copy()
    keys = ["island_id", "stratum", "syndrome"]
    duplicated = work.duplicated(keys, keep=False)
    if duplicated.any():
        examples = work.loc[duplicated, keys].head(5).to_dict("records")
        raise typer.BadParameter(f"duplicate island syndrome scores: {examples}")

    score_wide = work.pivot(
        index=["island_id", "stratum"],
        columns="syndrome",
        values="syndrome_score",
    )
    support_wide = work.pivot(
        index=["island_id", "stratum"],
        columns="syndrome",
        values="n_species",
    )

    rows: list[pd.DataFrame] = []
    for axis, spec in branching_config["branch_axes"].items():
        components = {str(k): float(v) for k, v in spec["components"].items()}
        component_names = list(components)
        if any(name not in score_wide.columns for name in component_names):
            continue
        scores = score_wide[component_names].apply(pd.to_numeric, errors="coerce")
        supports = support_wide[component_names].apply(pd.to_numeric, errors="coerce")
        complete = scores.notna().all(axis=1) & supports.notna().all(axis=1)
        if not complete.any():
            continue
        weights = pd.Series(components, dtype=float)
        axis_frame = scores.loc[complete].mul(weights, axis=1).sum(axis=1).rename(
            "syndrome_score"
        )
        out = axis_frame.reset_index()
        out["syndrome"] = str(axis)
        out["n_species"] = supports.loc[complete].min(axis=1).to_numpy(float)
        rows.append(out)

    if not rows:
        return pd.DataFrame(
            columns=["island_id", "stratum", "syndrome", "syndrome_score", "n_species"]
        )
    return pd.concat(rows, ignore_index=True)


def _run_context_layer(
    branch_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    pattern_config: dict[str, Any],
    *,
    layer_name: str,
    alpha: float,
    branching_config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    contexts = [str(x) for x in pattern_config["contexts"]]
    strata = [str(x) for x in pattern_config["strata"]]
    tiers = {str(k): int(v) for k, v in pattern_config["support_tiers"].items()}

    slope_parts: list[pd.DataFrame] = []
    within_rows: list[dict[str, Any]] = []
    between_rows: list[dict[str, Any]] = []
    for axis_set, axis_spec in branching_config["axis_sets"].items():
        axes = [str(x) for x in axis_spec["axes"]]
        selected = branch_scores.loc[branch_scores["syndrome"].isin(axes)].copy()
        data = _prepare(selected, covariates, pattern_config)
        for stratum in strata:
            for tier, threshold in tiers.items():
                for context in contexts:
                    slopes, result = _within_context(
                        data,
                        stratum=stratum,
                        context_value=context,
                        support_tier=tier,
                        threshold=threshold,
                        pattern_config=pattern_config,
                        syndrome_config={},
                    )
                    if not slopes.empty:
                        slopes.insert(0, "axis_set", str(axis_set))
                        slope_parts.append(slopes)
                    result["axis_set"] = str(axis_set)
                    result["axis_set_role"] = str(axis_spec["role"])
                    within_rows.append(result)
                for context_a, context_b in itertools.combinations(contexts, 2):
                    result = _between_contexts(
                        data,
                        stratum=stratum,
                        context_a=context_a,
                        context_b=context_b,
                        support_tier=tier,
                        threshold=threshold,
                        pattern_config=pattern_config,
                    )
                    result["axis_set"] = str(axis_set)
                    result["axis_set_role"] = str(axis_spec["role"])
                    between_rows.append(result)

    slopes = pd.concat(slope_parts, ignore_index=True) if slope_parts else pd.DataFrame()
    within = pd.DataFrame(within_rows)
    between = pd.DataFrame(between_rows)
    for frame in (slopes, within, between):
        frame.insert(0, "context_layer", layer_name)

    family = ["context_layer", "axis_set", "stratum", "support_tier"]
    if not slopes.empty and "p_value" in slopes.columns:
        slopes["q_axis_family"] = slopes.groupby(family, group_keys=False)["p_value"].transform(
            _bh
        )
        slopes["axis_supported"] = slopes["q_axis_family"].le(alpha).fillna(False)
    if "p_value" in within.columns:
        within["q_context_vector_family"] = within.groupby(
            family, group_keys=False
        )["p_value"].transform(_bh)
        within["context_vector_supported"] = within["q_context_vector_family"].le(
            alpha
        ).fillna(False)
    if "p_value" in between.columns:
        between["q_between_context_family"] = between.groupby(
            family, group_keys=False
        )["p_value"].transform(_bh)
        between["context_vector_difference_supported"] = between[
            "q_between_context_family"
        ].le(alpha).fillna(False)
    return slopes, within, between


def classify_branch_states(
    slopes: pd.DataFrame,
    within: pd.DataFrame,
    branching_config: dict[str, Any],
) -> pd.DataFrame:
    """Classify plant-side branches only after the multivariate context gate passes."""

    classified_sets = {
        str(name): spec
        for name, spec in branching_config["axis_sets"].items()
        if bool(spec.get("classify", False))
    }
    keys = ["context_layer", "axis_set", "stratum", "support_tier", "context"]
    rows: list[dict[str, Any]] = []
    for _, result in within.iterrows():
        axis_set = str(result.get("axis_set", ""))
        if axis_set not in classified_sets:
            continue
        axes = [str(x) for x in classified_sets[axis_set]["axes"]]
        base = {key: result.get(key) for key in keys}
        status = str(result.get("status", ""))
        vector_supported = bool(result.get("context_vector_supported", False))
        row: dict[str, Any] = {
            **base,
            "status": status,
            "threshold": result.get("threshold"),
            "n_unique_islands": result.get("n_unique_islands"),
            "n_clusters": result.get("n_clusters"),
            "context_vector_supported": vector_supported,
            "q_context_vector_family": result.get("q_context_vector_family"),
        }
        subset = slopes.copy()
        if not subset.empty:
            for key in keys:
                subset = subset.loc[
                    subset[key].astype(str).eq(str(result.get(key, "")))
                ]
        indexed = subset.set_index("syndrome") if not subset.empty else pd.DataFrame()

        positive: list[str] = []
        negative: list[str] = []
        for axis in axes:
            estimate = float("nan")
            q_value = float("nan")
            if not indexed.empty and axis in indexed.index:
                axis_row = indexed.loc[axis]
                if isinstance(axis_row, pd.DataFrame):
                    raise ValueError(f"duplicate fitted branch axis: {axis}")
                estimate = float(axis_row["distance_slope"])
                q_value = float(axis_row["q_axis_family"])
            supported = vector_supported and np.isfinite(q_value) and q_value <= float(
                branching_config["alpha"]
            )
            row[f"{axis}_slope"] = estimate
            row[f"{axis}_q"] = q_value
            row[f"{axis}_positive"] = bool(supported and estimate > 0)
            row[f"{axis}_negative"] = bool(supported and estimate < 0)
            if supported and estimate > 0:
                positive.append(axis)
            elif supported and estimate < 0:
                negative.append(axis)

        if status != "fit":
            branch_state = "not_testable"
        elif not vector_supported:
            branch_state = "no_supported_multivariate_branch"
        elif not positive and not negative:
            branch_state = "multivariate_branch_without_axiswise_resolution"
        else:
            parts: list[str] = []
            if "accessibility_generalization" in positive:
                parts.append("accessibility_generalization")
            if "accessibility_generalization" in negative:
                parts.append("specialized_access_maintenance_or_increase")
            if "reproductive_assurance" in positive:
                parts.append("reproductive_assurance")
            other = [
                axis
                for axis in [*positive, *negative]
                if axis not in {"accessibility_generalization", "reproductive_assurance"}
            ]
            parts.extend(other)
            branch_state = "|".join(parts) or "supported_axis_change"

        row["positive_axes"] = "|".join(positive)
        row["negative_axes"] = "|".join(negative)
        row["plant_side_branch_state"] = branch_state
        rows.append(row)
    return pd.DataFrame(rows)


def run_global_branching(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    branching_config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    branch_scores = build_branch_scores(island_scores, branching_config)
    all_slopes: list[pd.DataFrame] = []
    all_within: list[pd.DataFrame] = []
    all_between: list[pd.DataFrame] = []

    for layer_name, layer_spec in branching_config["context_layers"].items():
        column = str(layer_spec["column"])
        layer_covariates = covariates.copy()
        if column not in layer_covariates.columns:
            required = {"island_id", column}
            missing = required - set(realm_assignment.columns)
            if missing:
                raise typer.BadParameter(
                    f"context assignment missing columns: {sorted(missing)}"
                )
            assignment = realm_assignment[["island_id", column]].drop_duplicates(
                "island_id"
            )
            layer_covariates = layer_covariates.merge(
                assignment,
                on="island_id",
                how="left",
                validate="one_to_one",
            )
        layer_pattern = deepcopy(pattern_config)
        layer_pattern["context_column"] = column
        layer_pattern["contexts"] = [str(x) for x in layer_spec["contexts"]]
        slopes, within, between = _run_context_layer(
            branch_scores,
            layer_covariates,
            layer_pattern,
            layer_name=str(layer_name),
            alpha=float(branching_config["alpha"]),
            branching_config=branching_config,
        )
        all_slopes.append(slopes)
        all_within.append(within)
        all_between.append(between)

    slopes = pd.concat(all_slopes, ignore_index=True)
    within = pd.concat(all_within, ignore_index=True)
    between = pd.concat(all_between, ignore_index=True)
    states = classify_branch_states(slopes, within, branching_config)
    return branch_scores, slopes, within, between, states


@app.command("run")
def run(
    island_scores_csv: Annotated[Path, typer.Option(exists=True)],
    covariates_csv: Annotated[Path, typer.Option(exists=True)],
    realm_assignment_csv: Annotated[Path, typer.Option(exists=True)],
    pattern_config_path: Annotated[Path, typer.Option(exists=True)],
    branching_config_path: Annotated[Path, typer.Option(exists=True)],
    output_dir: Annotated[Path, typer.Option()],
) -> None:
    pattern_config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    branching_config = yaml.safe_load(branching_config_path.read_text(encoding="utf-8"))
    branch_scores, slopes, within, between, states = run_global_branching(
        pd.read_csv(island_scores_csv),
        pd.read_csv(covariates_csv),
        pd.read_csv(realm_assignment_csv),
        pattern_config,
        branching_config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    branch_scores.to_csv(output_dir / "island_plant_side_branch_scores.csv.gz", index=False)
    slopes.to_csv(output_dir / "global_branch_distance_slopes.csv", index=False)
    within.to_csv(output_dir / "global_branch_within_context_omnibus.csv", index=False)
    between.to_csv(output_dir / "global_branch_between_context_omnibus.csv", index=False)
    states.to_csv(output_dir / "global_plant_side_branch_states.csv", index=False)
    manifest = {
        "contract": str(branching_config["contract"]),
        "hypothesis_name": str(branching_config["hypothesis_name"]),
        "chapter_1_identifies": "global when/where plant-side branches",
        "chapter_2_repository": "https://github.com/zuizui0223/izu-core",
        "pollination_service_observed": False,
        "pollinator_identity_observed": False,
        "functional_replacement_identified": False,
        "primary_vector_contains_pollinator_names": False,
        "guild_concordance_is_non_exhaustive_secondary_evidence": True,
        "source_pool_standardized": False,
        "distance_interpretation": "composite isolation/connectivity/source-supply gradient",
        "n_branch_score_rows": len(branch_scores),
        "n_context_state_rows": len(states),
    }
    (output_dir / "global_branching_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
