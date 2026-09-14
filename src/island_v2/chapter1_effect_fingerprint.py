"""Effect-size synthesis for the frozen Chapter 1 result.

This module does not refit biological models. It reads already frozen slope tables and
re-expresses them as (i) source-matched taxonomic attenuation profiles and (ii) atomic
trait-domain fingerprints. No new p-values or mechanism labels are created here.
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

app = typer.Typer(add_completion=False, no_args_is_help=True)

SCOPE_PATHS = {
    "all_analysis_eligible": "all",
    "direct_only": "direct",
}


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != "chapter1_effect_fingerprint_v1":
        raise typer.BadParameter("unexpected effect-fingerprint contract")
    return config


def _ci_columns(frame: pd.DataFrame, estimate: str, se: str) -> pd.DataFrame:
    out = frame.copy()
    out[estimate] = pd.to_numeric(out[estimate], errors="coerce")
    out[se] = pd.to_numeric(out[se], errors="coerce")
    out["ci_low"] = out[estimate] - 1.96 * out[se]
    out["ci_high"] = out[estimate] + 1.96 * out[se]
    out["direction"] = np.where(
        out["ci_low"].gt(0),
        "positive",
        np.where(out["ci_high"].lt(0), "negative", "uncertain"),
    )
    return out


def build_taxonomic_attenuation(root: Path, config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    axes = [str(x) for x in config["primary_taxonomic_axes"]]
    stages = [str(x) for x in config["stages"]]
    source_modes = [str(x) for x in config["source_modes"]]
    strata = [str(x) for x in config["strata"]]
    target = config["primary_taxonomic_context"]
    rows: list[pd.DataFrame] = []
    for evidence_scope in SCOPE_PATHS:
        path = root / "taxonomic-depth" / evidence_scope / "slopes.csv"
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_csv(path)
        required = {
            "source_mode", "context_layer", "context", "stratum", "support_tier",
            "axis_set", "syndrome", "distance_slope", "cluster_robust_se", "p_value",
            "n_islands",
        }
        if missing := required - set(frame.columns):
            raise ValueError(f"taxonomic slope table lacks columns: {sorted(missing)}")
        work = frame.loc[
            frame["source_mode"].astype(str).isin(source_modes)
            & frame["context_layer"].astype(str).eq(str(target["context_layer"]))
            & frame["context"].astype(str).eq(str(target["context"]))
            & frame["stratum"].astype(str).isin(strata)
            & frame["support_tier"].astype(str).eq("confirmatory")
        ].copy()
        work["stage"] = work["axis_set"].astype(str).str.replace(
            "taxonomic_stage__", "", regex=False
        )
        work["axis"] = work["syndrome"].astype(str)
        for prefix in ("observed_score__", "after_family_residual__", "after_genus_residual__"):
            work["axis"] = work["axis"].str.replace(prefix, "", regex=False)
        work = work.loc[work["stage"].isin(stages) & work["axis"].isin(axes)].copy()
        work.insert(0, "evidence_scope", evidence_scope)
        rows.append(_ci_columns(work, "distance_slope", "cluster_robust_se"))
    long = pd.concat(rows, ignore_index=True)

    expected_axes = set(axes)
    vector_rows: list[dict[str, Any]] = []
    group_cols = ["evidence_scope", "source_mode", "stratum"]
    for key, group in long.groupby(group_cols, sort=True):
        meta = dict(zip(group_cols, key, strict=True))
        rec: dict[str, Any] = {
            **meta,
            "context_layer": str(target["context_layer"]),
            "context": str(target["context"]),
        }
        complete = True
        for stage in stages:
            part = group.loc[group["stage"].eq(stage)].copy()
            found = set(part["axis"].astype(str))
            if found != expected_axes or part["distance_slope"].isna().any():
                complete = False
                rec[f"{stage}_vector_norm"] = np.nan
                continue
            rec[f"{stage}_vector_norm"] = float(
                np.linalg.norm(part["distance_slope"].to_numpy(float))
            )
        rec["complete_two_axis_vector"] = bool(complete)
        observed = float(rec.get("observed_score_vector_norm", np.nan))
        floor = float(
            config["taxonomic_attenuation"]["do_not_compute_ratio_when_observed_norm_below"]
        )
        if complete and math.isfinite(observed) and observed > floor:
            family = float(rec["after_family_residual_vector_norm"])
            genus = float(rec["after_genus_residual_vector_norm"])
            rec["family_retention_fraction"] = family / observed
            rec["genus_retention_fraction"] = genus / observed
            rec["family_attenuation_fraction"] = 1.0 - family / observed
            rec["genus_attenuation_fraction"] = 1.0 - genus / observed
            rec["conditional_genus_attenuation_fraction"] = (
                1.0 - genus / family if family > floor else np.nan
            )
        else:
            rec["family_retention_fraction"] = np.nan
            rec["genus_retention_fraction"] = np.nan
            rec["family_attenuation_fraction"] = np.nan
            rec["genus_attenuation_fraction"] = np.nan
            rec["conditional_genus_attenuation_fraction"] = np.nan
        vector_rows.append(rec)
    vectors = pd.DataFrame(vector_rows)
    return long, vectors


def _domain_map(config: dict[str, Any]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for domain, outcomes in config["atomic_fingerprint"]["outcome_domains"].items():
        for outcome in outcomes:
            outcome = str(outcome)
            if outcome in mapping:
                raise ValueError(f"duplicate atomic outcome across domains: {outcome}")
            mapping[outcome] = str(domain)
    return mapping


def build_atomic_fingerprint(root: Path, config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    mapping = _domain_map(config)
    contexts = {str(x) for x in config["atomic_fingerprint"]["contexts"]}
    strata = {str(x) for x in config["strata"]}
    tier = str(config["atomic_fingerprint"]["support_tier"])
    frames: list[pd.DataFrame] = []
    for evidence_scope, short_scope in SCOPE_PATHS.items():
        path = root / "atomic" / short_scope / "observed_within_outcome_slopes.csv"
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_csv(path)
        required = {
            "stratum", "support_tier", "context", "outcome",
            "geography_slope_log_odds_per_response_sd", "cluster_robust_se",
            "p_value", "n_islands",
        }
        if missing := required - set(frame.columns):
            raise ValueError(f"atomic slope table lacks columns: {sorted(missing)}")
        work = frame.loc[
            frame["support_tier"].astype(str).eq(tier)
            & frame["context"].astype(str).isin(contexts)
            & frame["stratum"].astype(str).isin(strata)
            & frame["outcome"].astype(str).isin(mapping)
        ].copy()
        work["domain"] = work["outcome"].astype(str).map(mapping)
        work.insert(0, "evidence_scope", evidence_scope)
        frames.append(
            _ci_columns(
                work,
                "geography_slope_log_odds_per_response_sd",
                "cluster_robust_se",
            )
        )
    long = pd.concat(frames, ignore_index=True)

    summary_rows: list[dict[str, Any]] = []
    group_cols = ["evidence_scope", "stratum", "context", "domain"]
    for key, group in long.groupby(group_cols, sort=True):
        slopes = group["geography_slope_log_odds_per_response_sd"].to_numpy(float)
        directions = group["direction"].astype(str)
        summary_rows.append(
            {
                **dict(zip(group_cols, key, strict=True)),
                "n_outcomes": int(len(group)),
                "mean_classic_signed_slope": float(np.mean(slopes)),
                "euclidean_slope_norm": float(np.linalg.norm(slopes)),
                "n_positive_ci": int(directions.eq("positive").sum()),
                "n_negative_ci": int(directions.eq("negative").sum()),
                "n_uncertain": int(directions.eq("uncertain").sum()),
                "median_n_islands": float(
                    pd.to_numeric(group["n_islands"], errors="coerce").median()
                ),
            }
        )
    summary = pd.DataFrame(summary_rows)
    return long, summary


def build_cross_context_geometry(long: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    order = [str(x) for x in config["atomic_fingerprint"]["outcome_order"]]
    contexts = [str(x) for x in config["atomic_fingerprint"]["contexts"]]
    if len(contexts) != 2:
        raise ValueError("cross-context angle requires exactly two frozen contexts")
    rows: list[dict[str, Any]] = []
    for (scope, stratum), group in long.groupby(["evidence_scope", "stratum"], sort=True):
        vectors: dict[str, np.ndarray] = {}
        complete = True
        for context in contexts:
            part = group.loc[group["context"].astype(str).eq(context)].set_index("outcome")
            if not set(order).issubset(part.index):
                complete = False
                break
            vectors[context] = np.asarray(
                [part.loc[outcome, "geography_slope_log_odds_per_response_sd"] for outcome in order],
                dtype=float,
            )
        if not complete or any(not np.isfinite(vector).all() for vector in vectors.values()):
            rows.append(
                {
                    "evidence_scope": scope,
                    "stratum": stratum,
                    "context_a": contexts[0],
                    "context_b": contexts[1],
                    "n_components": len(order),
                    "cosine_similarity": np.nan,
                    "vector_angle_degrees": np.nan,
                    "complete_vector": False,
                }
            )
            continue
        a = vectors[contexts[0]]
        b = vectors[contexts[1]]
        denom = float(np.linalg.norm(a) * np.linalg.norm(b))
        cosine = float(np.dot(a, b) / denom) if denom > 0 else np.nan
        angle = float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))) if np.isfinite(cosine) else np.nan
        rows.append(
            {
                "evidence_scope": scope,
                "stratum": stratum,
                "context_a": contexts[0],
                "context_b": contexts[1],
                "n_components": len(order),
                "cosine_similarity": cosine,
                "vector_angle_degrees": angle,
                "complete_vector": True,
            }
        )
    return pd.DataFrame(rows)


def build_manifest(
    config: dict[str, Any],
    taxonomic_vectors: pd.DataFrame,
    fingerprint_summary: pd.DataFrame,
    cross_context: pd.DataFrame,
) -> dict[str, Any]:
    attenuation = pd.to_numeric(
        taxonomic_vectors["genus_attenuation_fraction"], errors="coerce"
    ).dropna()
    conditional = pd.to_numeric(
        taxonomic_vectors["conditional_genus_attenuation_fraction"], errors="coerce"
    ).dropna()
    angles = pd.to_numeric(cross_context["vector_angle_degrees"], errors="coerce").dropna()
    return {
        "contract": config["contract"],
        "status": "descriptive_effect_fingerprint_complete",
        "pinned_workflow_run_id": int(config["pinned_input"]["workflow_run_id"]),
        "pinned_artifact_id": int(config["pinned_input"]["artifact_id"]),
        "pinned_artifact_digest": str(config["pinned_input"]["digest"]),
        "n_taxonomic_vector_profiles": int(len(taxonomic_vectors)),
        "genus_attenuation_fraction_range": (
            [float(attenuation.min()), float(attenuation.max())] if len(attenuation) else []
        ),
        "genus_attenuation_fraction_median": (
            float(attenuation.median()) if len(attenuation) else None
        ),
        "conditional_genus_attenuation_fraction_range": (
            [float(conditional.min()), float(conditional.max())] if len(conditional) else []
        ),
        "cross_context_angle_degree_range": (
            [float(angles.min()), float(angles.max())] if len(angles) else []
        ),
        "n_atomic_domain_summaries": int(len(fingerprint_summary)),
        "new_p_values_generated": False,
        "mechanism_promoted": False,
        "claim_boundary": str(config["claim_boundary"]),
    }


def run_synthesis(*, artifact_root: Path, config_path: Path, output_dir: Path) -> dict[str, Any]:
    config = load_config(config_path)
    tax_long, tax_vectors = build_taxonomic_attenuation(artifact_root, config)
    fp_long, fp_summary = build_atomic_fingerprint(artifact_root, config)
    cross_context = build_cross_context_geometry(fp_long, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    tax_long.to_csv(output_dir / "taxonomic_effects_long.csv", index=False)
    tax_vectors.to_csv(output_dir / "taxonomic_vector_attenuation.csv", index=False)
    fp_long.to_csv(output_dir / "atomic_response_fingerprint.csv", index=False)
    fp_summary.to_csv(output_dir / "atomic_domain_fingerprint_summary.csv", index=False)
    cross_context.to_csv(output_dir / "atomic_cross_context_vector_geometry.csv", index=False)
    manifest = build_manifest(config, tax_vectors, fp_summary, cross_context)
    (output_dir / "chapter1_effect_fingerprint_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("run")
def run_command(
    artifact_root: Path = typer.Option(..., exists=True, file_okay=False),
    config_path: Path = typer.Option(
        Path("config/chapter1_effect_fingerprint.yml"), exists=True, dir_okay=False
    ),
    output_dir: Path = typer.Option(...),
) -> None:
    typer.echo(
        json.dumps(
            run_synthesis(
                artifact_root=artifact_root,
                config_path=config_path,
                output_dir=output_dir,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    app()
