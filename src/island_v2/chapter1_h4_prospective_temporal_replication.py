"""Prospective temporal replication of the Chapter 1 H4 functional bridge.

The discovery result in v13 remains post-hoc. This module implements a new,
temporally non-overlapping validation cohort (publications from 2016 onward)
with a hard outcome-blind preflight gate.

Workflow:
1. preflight: metadata + frozen trait states only; outcome columns are rejected.
2. commit the generated support lock and replace the frozen_commit placeholder.
3. analyse: only then may PL effect columns be read.

Passing this replication can confirm a current trait--pollen-limitation
association. It cannot establish historical mediation or trait evolution.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_v13_functional_bridge import (
    aggregate_species_measurement_cells,
    fit_global_trait_level,
    fit_within_group_trait_level,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT = "chapter1_h4_prospective_temporal_replication_v1"
PLACEHOLDER_COMMIT = "REQUIRED_BEFORE_UNBLINDING"
ARCH_COMPONENTS = (
    "generalized_form",
    "actinomorphic_symmetry",
    "shallow_open_tube",
)
MEASUREMENT_COLUMNS = (
    "PL_Effect_Size_Type1",
    "PL_Effect_Size_Type2",
    "Constant_added",
    "Level_of_Supplementation",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected prospective H4 replication contract")
    if config.get("inferential_role") != "prospective_temporal_replication":
        raise typer.BadParameter("prospective H4 contract has wrong inferential role")
    return config


def _validate_binary(series: pd.Series, name: str) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    observed = set(values.dropna().unique().tolist())
    if not observed.issubset({0.0, 1.0}):
        raise typer.BadParameter(f"{name} must contain only 0/1/missing")
    return values


def validate_metadata_outcome_blind(
    metadata: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    spec = config["outcome_blind_preflight"]
    required = {str(x) for x in spec["required_metadata_columns"]}
    if missing := required - set(metadata.columns):
        raise typer.BadParameter(f"metadata missing required columns: {sorted(missing)}")

    forbidden = {str(x).casefold() for x in spec["forbidden_columns"]}
    observed_forbidden = [
        column for column in metadata.columns if str(column).casefold() in forbidden
    ]
    if observed_forbidden:
        raise typer.BadParameter(
            "outcome columns are forbidden before support lock: "
            f"{sorted(observed_forbidden)}"
        )

    work = metadata.copy()
    if work["experiment_key"].astype(str).duplicated().any():
        raise typer.BadParameter("experiment_key must be unique during preflight")

    dates = pd.to_datetime(work["publication_date"], errors="coerce", utc=True)
    if dates.isna().any():
        raise typer.BadParameter("publication_date contains missing or invalid dates")
    start = pd.Timestamp(config["temporal_holdout"]["validation_publication_window"]["start_date"], tz="UTC")
    end = pd.Timestamp(config["temporal_holdout"]["validation_publication_window"]["end_date"], tz="UTC")
    outside = (dates < start) | (dates > end)
    if bool(outside.any()):
        bad = work.loc[outside, ["experiment_key", "publication_date"]].head(5)
        raise typer.BadParameter(
            "metadata includes studies outside frozen temporal holdout: "
            + bad.to_dict("records").__repr__()
        )

    for column in ("experiment_key", "publication_id", "study_key", "site_key", "accepted_species"):
        work[column] = work[column].fillna("").astype(str).str.strip()
        if work[column].eq("").any():
            raise typer.BadParameter(f"{column} cannot be blank in admitted metadata")
    work["doi"] = work["doi"].fillna("").astype(str).str.strip().str.casefold()
    work["analysis_regime"] = work["analysis_regime"].fillna("").astype(str).str.strip()
    work["z_distance"] = pd.to_numeric(work["z_distance"], errors="coerce")
    if not np.isfinite(work["z_distance"].to_numpy(float)).all():
        raise typer.BadParameter("z_distance must be finite before support locking")
    return work.reset_index(drop=True)


def prepare_frozen_trait_states(traits: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    required = {str(x) for x in config["trait_snapshot"]["required_columns"]}
    if missing := required - set(traits.columns):
        raise typer.BadParameter(f"trait-state table missing columns: {sorted(missing)}")
    work = traits[list(required)].copy()
    work["accepted_species"] = work["accepted_species"].fillna("").astype(str).str.strip()
    if work["accepted_species"].eq("").any():
        raise typer.BadParameter("trait accepted_species cannot be blank")
    if work["accepted_species"].duplicated().any():
        raise typer.BadParameter("trait-state table must have one row per accepted species")
    for column in ("autonomous_selfing", *ARCH_COMPONENTS):
        work[column] = _validate_binary(work[column], column)

    component_count = work[list(ARCH_COMPONENTS)].notna().sum(axis=1)
    component_mean = work[list(ARCH_COMPONENTS)].mean(axis=1, skipna=True)
    minimum = int(
        config["primary_hypotheses"]["H4b_accessibility_generalization"][
            "minimum_nonmissing_components_per_species"
        ]
    )
    work["accessibility_component_count"] = component_count
    work["accessibility_generalization_score"] = component_mean.where(
        component_count.ge(minimum)
    )
    return work


def match_metadata_to_traits(
    metadata: pd.DataFrame,
    traits: pd.DataFrame,
) -> pd.DataFrame:
    return metadata.merge(
        traits,
        on="accepted_species",
        how="left",
        validate="many_to_one",
    )


def _h4a_support(matched: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    spec = config["primary_hypotheses"]["H4a_reproductive_assurance"]["support_gate"]
    rows = matched.loc[matched["autonomous_selfing"].notna()].copy()
    species = rows[["accepted_species", "autonomous_selfing"]].drop_duplicates()
    publications = rows[["study_key", "autonomous_selfing"]].drop_duplicates()

    species_state = species.groupby("autonomous_selfing")["accepted_species"].nunique()
    pub_state = publications.groupby("autonomous_selfing")["study_key"].nunique()
    counts = {
        "matched_species_total": int(species["accepted_species"].nunique()),
        "publications_total": int(rows["study_key"].nunique()),
        "species_state_0": int(species_state.get(0.0, 0)),
        "species_state_1": int(species_state.get(1.0, 0)),
        "publications_state_0": int(pub_state.get(0.0, 0)),
        "publications_state_1": int(pub_state.get(1.0, 0)),
    }
    evaluable = (
        counts["matched_species_total"] >= int(spec["minimum_matched_species_total"])
        and min(counts["species_state_0"], counts["species_state_1"])
        >= int(spec["minimum_species_per_state"])
        and counts["publications_total"] >= int(spec["minimum_publications_total"])
        and min(counts["publications_state_0"], counts["publications_state_1"])
        >= int(spec["minimum_publications_per_state"])
    )
    return {"evaluable": bool(evaluable), **counts}


def _h4b_support(matched: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    spec = config["primary_hypotheses"]["H4b_accessibility_generalization"]["support_gate"]
    rows = matched.loc[matched["accessibility_generalization_score"].notna()].copy()
    species = rows[
        ["accepted_species", "accessibility_generalization_score"]
    ].drop_duplicates()
    scores = pd.to_numeric(species["accessibility_generalization_score"], errors="coerce")
    low = scores.le(float(spec["low_score_max"]))
    high = scores.ge(float(spec["high_score_min"]))
    score_sd = float(scores.std(ddof=0)) if len(scores) else float("nan")
    counts = {
        "matched_species_total": int(species["accepted_species"].nunique()),
        "publications_total": int(rows["study_key"].nunique()),
        "score_sd": score_sd,
        "low_score_species": int(low.sum()),
        "high_score_species": int(high.sum()),
    }
    evaluable = (
        counts["matched_species_total"] >= int(spec["minimum_matched_species_total"])
        and counts["publications_total"] >= int(spec["minimum_publications_total"])
        and np.isfinite(score_sd)
        and score_sd >= float(spec["minimum_score_sd"])
        and counts["low_score_species"] >= int(spec["minimum_low_score_species"])
        and counts["high_score_species"] >= int(spec["minimum_high_score_species"])
    )
    return {"evaluable": bool(evaluable), **counts}


def build_support_lock(
    metadata_path: Path,
    trait_states_path: Path,
    config_path: Path,
    metadata: pd.DataFrame,
    traits: pd.DataFrame,
    matched: pd.DataFrame,
    config: dict[str, Any],
) -> dict[str, Any]:
    keys = sorted(metadata["experiment_key"].astype(str).tolist())
    support = {
        "H4a_reproductive_assurance": _h4a_support(matched, config),
        "H4b_accessibility_generalization": _h4b_support(matched, config),
    }
    return {
        "contract": CONTRACT,
        "inferential_role": "prospective_temporal_replication",
        "outcome_status": "unopened",
        "frozen_commit": PLACEHOLDER_COMMIT,
        "contract_sha256": _sha256(config_path),
        "metadata_sha256": _sha256(metadata_path),
        "trait_states_sha256": _sha256(trait_states_path),
        "admitted_experiment_keys_sha256": _sha256_text("\n".join(keys) + "\n"),
        "n_admitted_experiments": int(len(keys)),
        "n_unique_publications": int(metadata["study_key"].nunique()),
        "n_unique_species": int(metadata["accepted_species"].nunique()),
        "support_counts": support,
        "evaluable_primary_hypotheses": [
            name for name, result in support.items() if result["evaluable"]
        ],
        "rules": {
            "outcomes_read": False,
            "thresholds_frozen": True,
            "trait_components_frozen": True,
        },
    }


def _verify_unblinding_lock(
    lock: dict[str, Any],
    *,
    metadata_path: Path,
    trait_states_path: Path,
    config_path: Path,
) -> None:
    if lock.get("contract") != CONTRACT:
        raise typer.BadParameter("support lock contract mismatch")
    if lock.get("outcome_status") != "unopened":
        raise typer.BadParameter("support lock must record outcome_status=unopened")
    frozen_commit = str(lock.get("frozen_commit", ""))
    if not frozen_commit or frozen_commit == PLACEHOLDER_COMMIT:
        raise typer.BadParameter(
            "support lock must be committed and frozen_commit filled before unblinding"
        )
    expected = {
        "contract_sha256": _sha256(config_path),
        "metadata_sha256": _sha256(metadata_path),
        "trait_states_sha256": _sha256(trait_states_path),
    }
    for key, observed in expected.items():
        if str(lock.get(key, "")) != observed:
            raise typer.BadParameter(f"support lock digest mismatch: {key}")


def _prepare_outcome_rows(
    outcomes: pd.DataFrame,
    matched: pd.DataFrame,
    locked_keys: set[str],
) -> pd.DataFrame:
    required = {"experiment_key", "PL_Effect_Size", *MEASUREMENT_COLUMNS}
    if missing := required - set(outcomes.columns):
        raise typer.BadParameter(f"outcome table missing columns: {sorted(missing)}")
    work = outcomes.copy()
    work["experiment_key"] = work["experiment_key"].fillna("").astype(str).str.strip()
    new_keys = sorted(set(work["experiment_key"]) - locked_keys)
    if new_keys:
        raise typer.BadParameter(
            "outcome table contains experiments not frozen in support lock: "
            f"{new_keys[:5]}"
        )
    work["PL_Effect_Size"] = pd.to_numeric(work["PL_Effect_Size"], errors="coerce")
    work = work.loc[np.isfinite(work["PL_Effect_Size"].to_numpy(float))].copy()
    joined = work.merge(
        matched,
        on="experiment_key",
        how="left",
        validate="many_to_one",
        suffixes=("", "_locked"),
    )
    if joined["accepted_species"].isna().any():
        raise typer.BadParameter("outcome row failed locked metadata join")
    return joined


def _predictor_rows(
    joined: pd.DataFrame,
    predictor: str,
    trait_label: str,
) -> pd.DataFrame:
    rows = joined.loc[joined[predictor].notna()].copy()
    rows["species_key"] = rows["accepted_species"].astype(str)
    rows["trait"] = trait_label
    rows["trait_state"] = pd.to_numeric(rows[predictor], errors="coerce")
    return rows[
        [
            "study_key",
            "site_key",
            "species_key",
            "analysis_regime",
            "z_distance",
            "trait",
            "trait_state",
            "PL_Effect_Size",
            *MEASUREMENT_COLUMNS,
        ]
    ].copy()


def _fit_primary(
    rows: pd.DataFrame,
    *,
    publication_total_weight: float,
) -> dict[str, Any]:
    cells = aggregate_species_measurement_cells(
        rows,
        publication_total_weight=publication_total_weight,
    )
    result = fit_global_trait_level(cells)
    result["n_input_effect_rows"] = int(len(rows))
    return result


def _fit_secondary(
    rows: pd.DataFrame,
    *,
    publication_total_weight: float,
) -> dict[str, Any]:
    cells = aggregate_species_measurement_cells(
        rows,
        publication_total_weight=publication_total_weight,
    )
    within_publication = fit_within_group_trait_level(cells, ["study_key"])
    within_site = fit_within_group_trait_level(cells, ["study_key", "site_key"])
    return {
        "within_publication": within_publication,
        "within_publication_site": within_site,
    }


def run_replication(
    joined: pd.DataFrame,
    lock: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    publication_total_weight = float(config["analysis"]["publication_total_weight"])
    definitions = {
        "H4a_reproductive_assurance": (
            "autonomous_selfing",
            "autonomous_selfing",
        ),
        "H4b_accessibility_generalization": (
            "accessibility_generalization_score",
            "accessibility_generalization",
        ),
    }
    results: dict[str, Any] = {}
    for name, (predictor, label) in definitions.items():
        support = lock["support_counts"][name]
        alpha = float(config["primary_hypotheses"][name]["one_sided_alpha"])
        if not support["evaluable"]:
            results[name] = {
                "evaluable": False,
                "reason": "frozen_support_gate_failed",
                "support": support,
            }
            continue
        rows = _predictor_rows(joined, predictor, label)
        primary = _fit_primary(
            rows,
            publication_total_weight=publication_total_weight,
        )
        supported = bool(
            primary.get("evaluable")
            and float(primary.get("trait_state_estimate", float("nan"))) < 0
            and float(primary.get("trait_state_one_sided_negative_p", 1.0)) <= alpha
        )
        results[name] = {
            "evaluable": bool(primary.get("evaluable")),
            "alpha": alpha,
            "prediction": "negative_current_PL_association",
            "primary": primary,
            "supported": supported,
            "secondary_sensitivities": _fit_secondary(
                rows,
                publication_total_weight=publication_total_weight,
            ),
        }

    supported_names = [
        name for name, result in results.items() if bool(result.get("supported"))
    ]
    return {
        "contract": CONTRACT,
        "inferential_role": "prospective_temporal_replication",
        "primary_results": results,
        "supported_primary_hypotheses": supported_names,
        "n_supported_primary_hypotheses": int(len(supported_names)),
        "claim_ceiling": config["claim_ceiling"],
    }


@app.command("preflight")
def preflight(
    metadata_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    trait_states_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    metadata = validate_metadata_outcome_blind(pd.read_csv(metadata_csv), config)
    traits = prepare_frozen_trait_states(pd.read_csv(trait_states_csv), config)
    matched = match_metadata_to_traits(metadata, traits)
    lock = build_support_lock(
        metadata_csv,
        trait_states_csv,
        config_path,
        metadata,
        traits,
        matched,
        config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    matched.to_csv(output_dir / "OUTCOME_BLIND_MATCHED_PREFLIGHT.csv", index=False)
    (output_dir / "SUPPORT_LOCK_TEMPLATE.json").write_text(
        json.dumps(lock, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(lock, indent=2))


@app.command("analyse")
def analyse(
    metadata_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    trait_states_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    outcomes_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    support_lock_json: Path = typer.Option(..., exists=True, dir_okay=False),
    config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    lock = json.loads(support_lock_json.read_text(encoding="utf-8"))
    _verify_unblinding_lock(
        lock,
        metadata_path=metadata_csv,
        trait_states_path=trait_states_csv,
        config_path=config_path,
    )

    metadata = validate_metadata_outcome_blind(pd.read_csv(metadata_csv), config)
    keys = sorted(metadata["experiment_key"].astype(str).tolist())
    key_digest = _sha256_text("\n".join(keys) + "\n")
    if key_digest != str(lock.get("admitted_experiment_keys_sha256", "")):
        raise typer.BadParameter("locked experiment-key cohort changed")

    traits = prepare_frozen_trait_states(pd.read_csv(trait_states_csv), config)
    matched = match_metadata_to_traits(metadata, traits)
    outcomes = pd.read_csv(outcomes_csv)
    joined = _prepare_outcome_rows(outcomes, matched, set(keys))
    result = run_replication(joined, lock, config)

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "RESULT.json").write_text(
        json.dumps(result, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(result, indent=2))


if __name__ == "__main__":
    app()
