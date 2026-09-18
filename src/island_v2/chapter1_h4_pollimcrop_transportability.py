"""Secondary PolLimCrop transportability test for Chapter 1 H4.

This module may read PolLimCrop pollen-limitation outcomes only after the frozen
outcome-blind preflight says at least one co-primary hypothesis is evaluable and
the response mapping has already been committed.

The crop-domain test is an independent transportability check. It cannot repair
the support-limited wild-plant temporal replication and cannot make the original
v13 H4 discovery confirmatory for wild floras.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_h4_prospective_temporal_replication import (
    prepare_frozen_trait_states,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT = "chapter1_h4_pollimcrop_response_mapping_v1"
HYPOTHESES = (
    "H4a_reproductive_assurance",
    "H4b_accessibility_generalization",
)
GRAIN = (
    "article_code",
    "species",
    "continent",
    "country",
    "locality",
    "experiment_year",
    "supplement_type",
    "scale",
    "crop_part",
)
ADJUSTMENT_CATEGORIES = ("continent", "supplement_type", "scale", "crop_part")
RAW_COLUMN_ALIASES = {
    "unicode": "record_unicode",
    "plant accession": "plant_accession",
    "crop part": "crop_part",
    "year_of_the_experiment": "experiment_year",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_mapping(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict) or value.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected PolLimCrop response-mapping contract")
    if value.get("status") != "frozen_before_row_level_PolLimCrop_outcome_read":
        raise typer.BadParameter("PolLimCrop response mapping is not frozen")
    return value


def load_preflight(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    expected = "chapter1_h4_pollimcrop_transportability_preflight_v1"
    lock_contract = "chapter1_h4_pollimcrop_transportability_preflight_result_lock_v1"
    if not isinstance(value, dict) or value.get("contract") != expected:
        raise typer.BadParameter("unexpected PolLimCrop preflight result")
    if value.get("lock_contract") != lock_contract:
        raise typer.BadParameter(
            "PolLimCrop outcome analysis requires the committed support-result lock"
        )
    if value.get("outcomes_read") is not False:
        raise typer.BadParameter("preflight result lock must remain outcome-blind")
    support = value.get("support", {})
    admitted = [
        name
        for name in HYPOTHESES
        if bool(support.get(name, {}).get("evaluable"))
    ]
    decision = value.get("decision", {})
    recorded = [str(x) for x in decision.get("admitted_hypotheses", [])]
    if recorded != admitted:
        raise typer.BadParameter(
            "PolLimCrop support-result lock admission list is inconsistent"
        )
    if bool(decision.get("outcome_extraction_authorized")) != bool(admitted):
        raise typer.BadParameter(
            "PolLimCrop support-result lock authorization is inconsistent"
        )
    return value


def evaluable_hypotheses(preflight: dict[str, Any]) -> list[str]:
    support = preflight.get("support", {})
    return [
        name
        for name in HYPOTHESES
        if bool(support.get(name, {}).get("evaluable"))
    ]


def require_support_before_outcome_read(preflight: dict[str, Any]) -> list[str]:
    admitted = evaluable_hypotheses(preflight)
    if not admitted:
        raise typer.BadParameter(
            "no PolLimCrop co-primary hypothesis passed the frozen support gate; "
            "outcome values must remain unread"
        )
    return admitted


def exact_binomial(value: object) -> str:
    tokens = str(value or "").strip().split()
    if len(tokens) != 2:
        return ""
    genus, epithet = tokens
    if not re.fullmatch(r"[A-Z][A-Za-z-]+", genus):
        return ""
    if not re.fullmatch(r"[a-z][A-Za-z-]+", epithet):
        return ""
    return f"{genus} {epithet}"


def admitted_species_for_outcome_read(
    trait_states: pd.DataFrame,
    parent_contract: dict[str, Any],
    preflight: dict[str, Any],
) -> set[str]:
    """Return only species belonging to hypotheses that passed frozen support gates."""

    admitted = require_support_before_outcome_read(preflight)
    frozen = prepare_frozen_trait_states(trait_states, parent_contract)
    mask = pd.Series(False, index=frozen.index)
    if "H4a_reproductive_assurance" in admitted:
        mask |= frozen["autonomous_selfing"].notna()
    if "H4b_accessibility_generalization" in admitted:
        mask |= frozen["accessibility_generalization_score"].notna()
    return set(frozen.loc[mask, "accepted_species"].astype(str))


def _unwrap_outer_record_text(
    dataset_csv: Path,
    *,
    encoding: str,
    wrapper: str,
) -> str:
    text = dataset_csv.read_bytes().decode(encoding)
    if wrapper == "none":
        return text
    if wrapper != "full_row_double_quoted_csv_field":
        raise typer.BadParameter(f"unexpected frozen PolLimCrop row wrapper: {wrapper}")
    reader = csv.reader(io.StringIO(text), delimiter=",", quotechar='"')
    records: list[str] = []
    for row in reader:
        if not row:
            continue
        if len(row) != 1:
            raise typer.BadParameter(
                "PolLimCrop outer wrapper did not decode to one field per record"
            )
        records.append(row[0])
    return "\n".join(records) + ("\n" if records else "")


def read_outcome_rows(
    dataset_csv: Path,
    preflight: dict[str, Any],
    *,
    allowed_species: set[str],
) -> pd.DataFrame:
    """Parse the frozen PL response only for species admitted by support gates."""

    require_support_before_outcome_read(preflight)
    if not allowed_species:
        raise typer.BadParameter("no admitted frozen species available for outcome read")
    expected_sha = str(preflight.get("figshare", {}).get("csv_sha256", ""))
    if not expected_sha:
        raise typer.BadParameter("preflight is missing the frozen PolLimCrop CSV SHA-256")
    observed_sha = _sha256(dataset_csv)
    if observed_sha != expected_sha:
        raise typer.BadParameter(
            f"PolLimCrop CSV SHA-256 mismatch: expected {expected_sha}, observed {observed_sha}"
        )

    source_format = preflight.get("source_format", {})
    delimiter = str(source_format.get("delimiter", ""))
    decimal_mark = str(source_format.get("decimal_mark", ""))
    encoding = str(source_format.get("encoding", ""))
    wrapper = str(source_format.get("outer_record_wrapping", ""))
    if (
        delimiter != ";"
        or decimal_mark != ","
        or encoding != "latin-1"
        or wrapper != "full_row_double_quoted_csv_field"
    ):
        raise typer.BadParameter("unexpected frozen PolLimCrop CSV format")

    response_column = str(
        preflight.get("response_mapping", {}).get("response_column", "")
    )
    if response_column != "PL_effect_size":
        raise typer.BadParameter("unexpected frozen PolLimCrop response column")

    normalized = _unwrap_outer_record_text(
        dataset_csv,
        encoding=encoding,
        wrapper=wrapper,
    )
    reader = csv.DictReader(io.StringIO(normalized), delimiter=delimiter)
    raw_fieldnames = [str(x) for x in (reader.fieldnames or [])]
    canonical_fieldnames = {
        RAW_COLUMN_ALIASES.get(column, column) for column in raw_fieldnames
    }
    required = set(GRAIN) | {response_column}
    if missing := required - canonical_fieldnames:
        raise typer.BadParameter(
            f"PolLimCrop outcome schema missing columns: {sorted(missing)}"
        )

    def parse_number(value: object) -> float:
        text = str(value or "").strip()
        if decimal_mark != ".":
            text = text.replace(decimal_mark, ".")
        return float(text)

    selected: list[dict[str, Any]] = []
    for raw in reader:
        canonical = {
            RAW_COLUMN_ALIASES.get(str(key), str(key)): value
            for key, value in raw.items()
            if key is not None
        }
        accepted_species = exact_binomial(canonical.get("species", ""))
        if accepted_species not in allowed_species:
            # Do not access the outcome field for support-failed species.
            continue
        try:
            effect = parse_number(canonical.get(response_column, ""))
            experiment_year = parse_number(canonical.get("experiment_year", ""))
        except ValueError:
            continue
        if not np.isfinite(effect) or not np.isfinite(experiment_year):
            continue
        article_code = str(canonical.get("article_code", "")).strip()
        if not article_code:
            continue
        row = {column: canonical.get(column, "") for column in GRAIN}
        row["PL_effectsize"] = effect
        row["experiment_year"] = experiment_year
        row["accepted_species"] = accepted_species
        selected.append(row)
    return pd.DataFrame(
        selected,
        columns=[*GRAIN, "PL_effectsize", "accepted_species"],
    )


def attach_frozen_traits(
    rows: pd.DataFrame,
    trait_states: pd.DataFrame,
    parent_contract: dict[str, Any],
) -> pd.DataFrame:
    frozen = prepare_frozen_trait_states(trait_states, parent_contract)
    return rows.merge(
        frozen,
        on="accepted_species",
        how="left",
        validate="many_to_one",
    )


def aggregate_analysis_cells(rows: pd.DataFrame) -> pd.DataFrame:
    """Aggregate exact frozen response-mapping cells and weight publications equally."""

    if rows.empty:
        return rows.copy()
    group_cols = [*GRAIN, "accepted_species"]
    out = (
        rows.groupby(group_cols, as_index=False, dropna=False)
        .agg(
            PL_effectsize=("PL_effectsize", "mean"),
            autonomous_selfing=("autonomous_selfing", "first"),
            accessibility_generalization_score=(
                "accessibility_generalization_score",
                "first",
            ),
            n_effect_rows=("PL_effectsize", "size"),
        )
        .reset_index(drop=True)
    )
    n_cells = out.groupby("article_code")["accepted_species"].transform("size").astype(float)
    out["analysis_weight"] = 1.0 / n_cells
    years = pd.to_numeric(out["experiment_year"], errors="coerce").to_numpy(float)
    sd = float(np.std(years, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        out["z_experiment_year"] = 0.0
    else:
        out["z_experiment_year"] = (years - float(np.mean(years))) / sd
    return out


def _categorical_columns(
    frame: pd.DataFrame,
    column: str,
) -> tuple[list[np.ndarray], list[str]]:
    values = frame[column].fillna("").astype(str)
    levels = sorted(value for value in values.unique() if value)
    if len(levels) <= 1:
        return [], []
    reference = levels[0]
    columns: list[np.ndarray] = []
    names: list[str] = []
    for level in levels[1:]:
        columns.append(values.eq(level).to_numpy(float))
        names.append(f"{column}={level}")
    if reference:
        return columns, names
    return columns, names


def _full_rank_design(
    frame: pd.DataFrame,
    predictor: str,
) -> tuple[np.ndarray, list[str]]:
    predictor_values = pd.to_numeric(frame[predictor], errors="coerce").to_numpy(float)
    year = pd.to_numeric(frame["z_experiment_year"], errors="coerce").to_numpy(float)
    columns: list[np.ndarray] = [
        np.ones(len(frame), dtype=float),
        predictor_values,
        year,
    ]
    names = ["intercept", predictor, "z_experiment_year"]
    for category in ADJUSTMENT_CATEGORIES:
        local_columns, local_names = _categorical_columns(frame, category)
        columns.extend(local_columns)
        names.extend(local_names)

    raw = np.column_stack(columns)
    if not np.isfinite(raw).all():
        raise typer.BadParameter("non-finite PolLimCrop design matrix")

    keep: list[int] = []
    current = np.empty((len(frame), 0), dtype=float)
    rank = 0
    for i in range(raw.shape[1]):
        candidate = np.column_stack([current, raw[:, i]])
        new_rank = int(np.linalg.matrix_rank(candidate))
        if new_rank > rank:
            keep.append(i)
            current = candidate
            rank = new_rank
    retained_names = [names[i] for i in keep]
    if predictor not in retained_names:
        raise typer.BadParameter("frozen trait predictor is not estimable in PolLimCrop design")
    return raw[:, keep], retained_names


def _cluster_meat(
    X: np.ndarray,
    weighted_residual: np.ndarray,
    labels: np.ndarray,
) -> np.ndarray:
    p = X.shape[1]
    meat = np.zeros((p, p), dtype=float)
    for label in pd.unique(labels):
        mask = labels == label
        score = X[mask].T @ weighted_residual[mask]
        meat += np.outer(score, score)
    return meat


def _cluster_correction(n: int, p: int, g: int) -> float:
    if g <= 1 or n <= p:
        return 1.0
    return (g / (g - 1.0)) * ((n - 1.0) / (n - p))


def _p2(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0)) if math.isfinite(z) else float("nan")


def _lower(z: float) -> float:
    return 0.5 * math.erfc(-z / math.sqrt(2.0)) if math.isfinite(z) else float("nan")


def fit_two_way_clustered_trait(
    cells: pd.DataFrame,
    predictor: str,
) -> dict[str, Any]:
    part = cells.loc[
        cells[predictor].notna()
        & np.isfinite(pd.to_numeric(cells[predictor], errors="coerce").to_numpy(float))
    ].copy()
    if part.empty or part[predictor].nunique() < 2:
        return {"evaluable": False, "reason": "predictor_support_not_estimable"}

    # The frozen mapping assigns each publication total weight 1.0 for each
    # co-primary trait analysis. Recompute after predictor-specific filtering.
    retained_per_publication = (
        part.groupby("article_code")["accepted_species"].transform("size").astype(float)
    )
    part["analysis_weight"] = 1.0 / retained_per_publication

    X, names = _full_rank_design(part, predictor)
    y = pd.to_numeric(part["PL_effectsize"], errors="coerce").to_numpy(float)
    w = pd.to_numeric(part["analysis_weight"], errors="coerce").to_numpy(float)
    n, p = X.shape
    if n <= p + 2 or np.linalg.matrix_rank(X) < p:
        return {"evaluable": False, "reason": "design_not_full_rank"}

    xtwx = X.T @ (w[:, None] * X)
    bread = np.linalg.inv(xtwx)
    beta = bread @ (X.T @ (w * y))
    residual = y - X @ beta
    weighted_residual = w * residual

    article = part["article_code"].astype(str).to_numpy()
    species = part["accepted_species"].astype(str).to_numpy()
    intersection = np.array(
        [f"{a}\x1f{s}" for a, s in zip(article, species, strict=True)],
        dtype=object,
    )

    g_article = int(pd.Series(article).nunique())
    g_species = int(pd.Series(species).nunique())
    g_intersection = int(pd.Series(intersection).nunique())
    if min(g_article, g_species) < 2:
        return {"evaluable": False, "reason": "insufficient_two_way_clusters"}

    meat_article = _cluster_meat(X, weighted_residual, article)
    meat_species = _cluster_meat(X, weighted_residual, species)
    meat_intersection = _cluster_meat(X, weighted_residual, intersection)
    meat = (
        _cluster_correction(n, p, g_article) * meat_article
        + _cluster_correction(n, p, g_species) * meat_species
        - _cluster_correction(n, p, g_intersection) * meat_intersection
    )
    cov = bread @ meat @ bread
    se = np.sqrt(np.clip(np.diag(cov), 0.0, None))

    index = {name: i for i, name in enumerate(names)}
    i = index[predictor]
    estimate = float(beta[i])
    predictor_se = float(se[i])
    z = estimate / predictor_se if predictor_se > 0 else float("nan")
    return {
        "evaluable": True,
        "predictor": predictor,
        "estimate": estimate,
        "se": predictor_se,
        "z": z,
        "two_sided_p": _p2(z),
        "one_sided_negative_p": _lower(z),
        "n_cells": int(len(part)),
        "n_publications": g_article,
        "n_species": g_species,
        "n_publication_species_clusters": g_intersection,
        "retained_design_columns": names,
    }


def run_transportability(
    cells: pd.DataFrame,
    preflight: dict[str, Any],
    mapping: dict[str, Any],
) -> dict[str, Any]:
    admitted = require_support_before_outcome_read(preflight)
    predictors = {
        "H4a_reproductive_assurance": "autonomous_selfing",
        "H4b_accessibility_generalization": "accessibility_generalization_score",
    }
    results: dict[str, Any] = {}
    for name in HYPOTHESES:
        if name not in admitted:
            results[name] = {
                "evaluable": False,
                "reason": "frozen_preflight_support_gate_failed",
                "outcomes_used": False,
            }
            continue
        predictor = predictors[name]
        fit = fit_two_way_clustered_trait(cells, predictor)
        alpha = float(mapping["primary_model"]["alpha_per_coprimary_hypothesis"])
        supported = bool(
            fit.get("evaluable")
            and float(fit.get("estimate", float("nan"))) < 0
            and float(fit.get("one_sided_negative_p", 1.0)) <= alpha
        )
        results[name] = {
            **fit,
            "alpha": alpha,
            "expected_direction": "negative",
            "supported": supported,
            "outcomes_used": bool(fit.get("evaluable")),
        }

    supported = [name for name, result in results.items() if result.get("supported")]
    return {
        "contract": CONTRACT,
        "inferential_role": "secondary_external_domain_transportability",
        "wild_temporal_replication_repaired": False,
        "primary_results": results,
        "supported_primary_hypotheses": supported,
        "n_supported_primary_hypotheses": int(len(supported)),
        "claim_ceiling": mapping["claim_ceiling"],
    }


@app.command("analyse")
def analyse(
    dataset_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    preflight_json: Path = typer.Option(..., exists=True, dir_okay=False),
    trait_states_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    parent_contract_path: Path = typer.Option(..., exists=True, dir_okay=False),
    response_mapping_path: Path = typer.Option(..., exists=True, dir_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    mapping = load_mapping(response_mapping_path)
    preflight = load_preflight(preflight_json)

    # Hard firewall: this check occurs before PL_effectsize is parsed.
    require_support_before_outcome_read(preflight)

    parent = yaml.safe_load(parent_contract_path.read_text(encoding="utf-8"))
    if not isinstance(parent, dict):
        raise typer.BadParameter("invalid parent H4 contract")
    trait_states = pd.read_csv(trait_states_csv)
    allowed_species = admitted_species_for_outcome_read(
        trait_states,
        parent,
        preflight,
    )

    rows = read_outcome_rows(
        dataset_csv,
        preflight,
        allowed_species=allowed_species,
    )
    rows = attach_frozen_traits(rows, trait_states, parent)
    cells = aggregate_analysis_cells(rows)
    result = run_transportability(cells, preflight, mapping)

    output_dir.mkdir(parents=True, exist_ok=True)
    cells[
        [
            "article_code",
            "accepted_species",
            "continent",
            "country",
            "locality",
            "experiment_year",
            "supplement_type",
            "scale",
            "crop_part",
            "autonomous_selfing",
            "accessibility_generalization_score",
            "analysis_weight",
        ]
    ].to_csv(output_dir / "ANALYSIS_CELLS_METADATA_AND_TRAITS.csv", index=False)
    (output_dir / "RESULT.json").write_text(
        json.dumps(result, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(result, indent=2))


if __name__ == "__main__":
    app()
