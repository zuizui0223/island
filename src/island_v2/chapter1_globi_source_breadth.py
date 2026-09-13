"""Source-side GloBI functional-channel breadth extension for Chapter 1 H3.

The workflow has two deliberately separated stages.

1. ``predictor`` reads only the pinned GloBI SUPPORTS export plus the frozen plant
   taxonomy and builds genus-level sampled functional-channel breadth. No island
   occurrence, distance, source assignment or trait outcome is available at this stage.
2. ``analyse`` takes that frozen genus predictor and asks whether represented island
   genera are enriched for broader sampled interaction channels relative to source-
   available genera matched on source prevalence and source species richness.

This is D3 association evidence. It cannot rescue N1, open N2, or identify pollinator
loss/effective service.
"""
from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_nee_channel_catalog import (
    _read_chunks,
    filter_interaction_chunk,
    load_config as load_channel_config,
)
from island_v2.chapter1_nee_channel_catalog_stream import verify_pinned_archives
from island_v2.chapter1_pr138_lineage_representation_bridge import (
    _availability_matrices,
    _bh,
    _binomial_key,
    _genus,
    _normalise_name,
    _representation_count_matrix,
    _source_assignment_matrix,
    broad_source_availability,
    compute_island_enrichment,
)
from island_v2.status_stratified_lineage_analysis import fit_weighted_linear_clustered

app = typer.Typer(add_completion=False, no_args_is_help=True)


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != "chapter1_globi_source_breadth_v1":
        raise typer.BadParameter("unexpected GloBI source-breadth contract")
    return config


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class PlantTaxonomyMatcher:
    def __init__(self, taxonomy: pd.DataFrame) -> None:
        required = {"accepted_species", "genus"}
        if missing := required - set(taxonomy.columns):
            raise ValueError(f"taxonomy lacks columns: {sorted(missing)}")
        frame = taxonomy[["accepted_species", "genus"]].copy().fillna("")
        frame["accepted_species"] = frame["accepted_species"].astype(str).str.strip()
        frame["genus"] = frame["genus"].astype(str).str.strip()
        frame = frame.loc[frame["accepted_species"].ne("")].drop_duplicates("accepted_species")
        frame["norm"] = frame["accepted_species"].map(_normalise_name)
        frame["binomial"] = frame["accepted_species"].map(_binomial_key)

        exact_groups = frame.groupby("norm")["accepted_species"].agg(lambda x: sorted(set(x)))
        self.exact = {key: values[0] for key, values in exact_groups.items() if key and len(values) == 1}
        binomial_groups = frame.groupby("binomial")["accepted_species"].agg(lambda x: sorted(set(x)))
        self.binomial = {
            key: values[0] for key, values in binomial_groups.items() if key and len(values) == 1
        }
        self.genus = frame.set_index("accepted_species")["genus"].to_dict()

    def match(self, value: object) -> tuple[str, str, str]:
        norm = _normalise_name(value)
        accepted = self.exact.get(norm, "")
        route = "exact" if accepted else ""
        if not accepted:
            key = _binomial_key(value)
            accepted = self.binomial.get(key, "") if key else ""
            route = "binomial" if accepted else ""
        if not accepted:
            return "", "", "unmatched"
        genus = str(self.genus.get(accepted, "")).strip() or _genus(accepted)
        if not genus:
            return "", "", "unmatched"
        return accepted, genus, route


class PlantBreadthAccumulator:
    def __init__(self, matcher: PlantTaxonomyMatcher, channels: list[str]) -> None:
        self.matcher = matcher
        self.channels = tuple(channels)
        self._stats: dict[str, dict[str, Any]] = {}
        self.n_input_evidence_rows = 0
        self.n_matched_evidence_rows = 0
        self.n_unmatched_evidence_rows = 0
        self.n_exact_rows = 0
        self.n_binomial_rows = 0

    def add(self, evidence: pd.DataFrame) -> None:
        if evidence.empty:
            return
        required = {"plant_taxon", "channel_id", "reference_key", "evidence_strength"}
        if missing := required - set(evidence.columns):
            raise ValueError(f"flower evidence lacks columns: {sorted(missing)}")
        for row in evidence.to_dict("records"):
            self.n_input_evidence_rows += 1
            channel = str(row["channel_id"])
            reference = str(row["reference_key"] or "").strip()
            if channel not in self.channels or not reference:
                continue
            species, genus, route = self.matcher.match(row["plant_taxon"])
            if not genus:
                self.n_unmatched_evidence_rows += 1
                continue
            self.n_matched_evidence_rows += 1
            if route == "exact":
                self.n_exact_rows += 1
            elif route == "binomial":
                self.n_binomial_rows += 1
            stats = self._stats.setdefault(
                genus,
                {
                    "species": set(),
                    "references": set(),
                    "channel_references": defaultdict(set),
                    "n_interaction_rows": 0,
                    "n_pollination_rows": 0,
                    "n_flower_visit_rows": 0,
                },
            )
            stats["species"].add(species)
            stats["references"].add(reference)
            stats["channel_references"][channel].add(reference)
            stats["n_interaction_rows"] += 1
            if str(row["evidence_strength"]) == "pollination":
                stats["n_pollination_rows"] += 1
            elif str(row["evidence_strength"]) == "flower_visit":
                stats["n_flower_visit_rows"] += 1

    def to_frame(self) -> pd.DataFrame:
        rows: list[dict[str, Any]] = []
        for genus, stats in sorted(self._stats.items()):
            channel_counts = {
                channel: len(stats["channel_references"].get(channel, set()))
                for channel in self.channels
            }
            total_channel_reference_units = int(sum(channel_counts.values()))
            positive = [value for value in channel_counts.values() if value > 0]
            if positive and total_channel_reference_units > 0:
                probabilities = np.asarray(positive, dtype=float) / total_channel_reference_units
                entropy = float(-np.sum(probabilities * np.log(probabilities)))
                effective = float(np.exp(entropy))
            else:
                effective = float("nan")
            record: dict[str, Any] = {
                "genus": genus,
                "n_matched_plant_species": len(stats["species"]),
                "n_independent_references": len(stats["references"]),
                "n_interaction_rows": int(stats["n_interaction_rows"]),
                "n_pollination_rows": int(stats["n_pollination_rows"]),
                "n_flower_visit_rows": int(stats["n_flower_visit_rows"]),
                "sampled_channel_count": int(sum(value > 0 for value in channel_counts.values())),
                "effective_channel_number": effective,
                "n_reference_channel_units": total_channel_reference_units,
            }
            for channel in self.channels:
                record[f"refs__{channel}"] = int(channel_counts[channel])
            rows.append(record)
        return pd.DataFrame(rows)


def build_genus_breadth_predictor(
    *,
    interactions_tsv_gz: Path,
    refuted_tsv_gz: Path,
    taxonomy_csv: Path,
    config_path: Path,
    channel_config_path: Path,
    output_dir: Path,
    chunksize: int = 200_000,
) -> dict[str, Any]:
    config = load_config(config_path)
    channel_config = load_channel_config(channel_config_path)
    policy = config["source_policy"]
    if str(channel_config["source_policy"]["resolved_version_doi"]) != str(policy["resolved_version_doi"]):
        raise ValueError("GloBI source version differs between frozen contracts")
    if str(channel_config["source_policy"].get("interactions_export_argument_type")) != "SUPPORTS":
        raise ValueError("source-breadth predictor requires the SUPPORTS export")
    if bool(channel_config["source_policy"].get("refuted_interactions_must_be_subtracted")):
        raise ValueError("source-breadth predictor prohibits SUPPORTS/REFUTES cross-subtraction")

    source_digest, refuted_digest = verify_pinned_archives(
        interactions_tsv_gz,
        refuted_tsv_gz,
        channel_config,
        str(policy["resolved_version_doi"]),
    )
    if source_digest != str(policy["interactions_tsv_sha256"]):
        raise ValueError("source digest differs from source-breadth freeze")
    if refuted_digest != str(policy["refuted_interactions_tsv_sha256"]):
        raise ValueError("refuted digest differs from source-breadth freeze")

    taxonomy = pd.read_csv(taxonomy_csv, dtype=str).fillna("")
    matcher = PlantTaxonomyMatcher(taxonomy)
    channels = [str(x) for x in config["functional_channels"]]
    accumulator = PlantBreadthAccumulator(matcher, channels)
    holdout_counts: defaultdict[str, int] = defaultdict(int)
    n_input_archive_rows = 0
    n_retained_flower_rows = 0
    required = list(channel_config["required_columns"])
    for chunk in _read_chunks(interactions_tsv_gz, required, chunksize):
        n_input_archive_rows += int(len(chunk))
        evidence, holdouts, _ = filter_interaction_chunk(chunk, channel_config, None)
        n_retained_flower_rows += int(len(evidence))
        accumulator.add(evidence)
        if not holdouts.empty:
            for reason, count in holdouts["reason"].value_counts().items():
                holdout_counts[str(reason)] += int(count)

    breadth = accumulator.to_frame()
    output_dir.mkdir(parents=True, exist_ok=True)
    breadth.to_csv(output_dir / "globi_source_genus_breadth.csv", index=False)
    receipt = {
        "contract": config["contract"],
        "status": "source_predictor_built_before_island_join",
        "source_version_doi": str(policy["resolved_version_doi"]),
        "source_sha256": source_digest,
        "refuted_sha256": refuted_digest,
        "support_refute_cross_subtraction_applied": False,
        "n_input_archive_rows": n_input_archive_rows,
        "n_retained_flower_interaction_rows": n_retained_flower_rows,
        "n_matched_evidence_rows": accumulator.n_matched_evidence_rows,
        "n_unmatched_evidence_rows": accumulator.n_unmatched_evidence_rows,
        "n_exact_match_rows": accumulator.n_exact_rows,
        "n_binomial_fallback_rows": accumulator.n_binomial_rows,
        "n_genus_breadth_rows": int(len(breadth)),
        "functional_channels": channels,
        "island_outcomes_or_distance_loaded": False,
        "no_record_interpretation": "missing_not_specialist",
        "breadth_file_sha256": sha256_file(output_dir / "globi_source_genus_breadth.csv"),
        "holdout_counts": dict(sorted(holdout_counts.items())),
    }
    (output_dir / "globi_source_genus_breadth_receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def _all_effort_thresholds(config: dict[str, Any]) -> list[int]:
    gate = config["sampling_effort_gate"]
    return sorted(
        {
            int(gate["primary_min_independent_references"]),
            *[int(x) for x in gate["sensitivity_min_independent_references"]],
        }
    )


def build_island_breadth_enrichment(
    *,
    genus_breadth: pd.DataFrame,
    gift_flora: pd.DataFrame,
    assignments: pd.DataFrame,
    status_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    required_breadth = {
        "genus", "n_independent_references", "sampled_channel_count", "effective_channel_number"
    }
    if missing := required_breadth - set(genus_breadth.columns):
        raise ValueError(f"genus breadth lacks columns: {sorted(missing)}")
    required_cov = {
        "island_id",
        "analysis_regime",
        *[str(x) for x in config["model"]["predictors"]],
        str(config["model"]["cluster_column"]),
    }
    if missing := required_cov - set(covariates.columns):
        raise ValueError(f"covariates lack source-breadth columns: {sorted(missing)}")

    islands = sorted(covariates["island_id"].astype(str).unique())
    island_index = {island: index for index, island in enumerate(islands)}
    source_modes = [str(x) for x in config["source_matching"]["source_modes"]]
    strata = [str(x) for x in config["source_matching"]["strata"]]
    matchings = [
        str(config["source_matching"]["primary_matching"]),
        str(config["source_matching"]["sensitivity_matching"]),
    ]
    min_repr = int(config["source_matching"]["minimum_represented_genera"])
    metrics = [str(x) for x in config["model"]["response_metrics"]]
    thresholds = _all_effort_thresholds(config)

    parts: list[pd.DataFrame] = []
    for threshold in thresholds:
        eligible = genus_breadth.copy()
        eligible["n_independent_references"] = pd.to_numeric(
            eligible["n_independent_references"], errors="coerce"
        )
        eligible = eligible.loc[eligible["n_independent_references"].ge(threshold)].copy()
        eligible = eligible.dropna(subset=metrics)
        genera = sorted(set(eligible["genus"].astype(str)) - {""})
        if not genera:
            continue
        genus_index = {genus: index for index, genus in enumerate(genera)}
        availability = broad_source_availability(gift_flora, set(genera))
        entities = sorted(
            set(pd.to_numeric(assignments["entity_ID"], errors="coerce").dropna().astype(int))
            | set(pd.to_numeric(availability["entity_ID"], errors="coerce").dropna().astype(int))
        )
        entity_index = {entity: index for index, entity in enumerate(entities)}
        presence, richness = _availability_matrices(availability, entity_index, genus_index)
        positions = {
            metric: eligible.drop_duplicates("genus").set_index("genus").loc[genera, metric].to_numpy(float)
            for metric in metrics
        }
        counts_by_stratum = {
            stratum: _representation_count_matrix(
                status_flora,
                island_index,
                genus_index,
                stratum=stratum,
                species_filter=None,
            ).toarray()
            for stratum in strata
        }
        for source_mode in source_modes:
            assignment_matrix = _source_assignment_matrix(
                assignments,
                island_index,
                entity_index,
                source_mode=source_mode,
            )
            prevalence = (assignment_matrix @ presence).toarray().astype(np.int16)
            source_richness = (assignment_matrix @ richness).toarray().astype(np.float32)
            for stratum in strata:
                island_counts = counts_by_stratum[stratum]
                for matching in matchings:
                    for metric in metrics:
                        rows: list[dict[str, Any]] = []
                        for island_position, island_id in enumerate(islands):
                            result = compute_island_enrichment(
                                prevalence[island_position],
                                source_richness[island_position],
                                island_counts[island_position],
                                positions[metric],
                                matching=matching,
                                minimum_represented_genera=min_repr,
                            )
                            if result is None:
                                continue
                            rows.append(
                                {
                                    "island_id": island_id,
                                    "metric": metric,
                                    "min_independent_references": threshold,
                                    "stratum": stratum,
                                    "source_mode": source_mode,
                                    "source_matching": matching,
                                    **result,
                                }
                            )
                        if rows:
                            parts.append(pd.DataFrame(rows))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def fit_breadth_models(
    enrichment: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    if enrichment.empty:
        return pd.DataFrame()
    predictors = [str(x) for x in config["model"]["predictors"]]
    cluster = str(config["model"]["cluster_column"])
    contexts = [str(x) for x in config["model"]["contexts"]]
    outcomes = [str(x) for x in config["island_enrichment"]["outcomes"]]
    weights = {str(k): str(v) for k, v in config["island_enrichment"]["weights"].items()}
    needed_cov = ["island_id", "analysis_regime", cluster, *predictors]
    data = enrichment.merge(
        covariates[needed_cov].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    group_cols = [
        "metric", "min_independent_references", "stratum", "source_mode", "source_matching"
    ]
    rows: list[dict[str, Any]] = []
    for key, group in data.groupby(group_cols, sort=True):
        meta = dict(zip(group_cols, key, strict=True))
        for context in contexts:
            context_data = group.loc[group["analysis_regime"].astype(str).eq(context)].copy()
            for outcome in outcomes:
                n_islands = int(context_data[["island_id", outcome]].dropna()["island_id"].nunique())
                if n_islands >= int(config["model"]["confirmatory_min_islands"]):
                    support = "confirmatory"
                elif n_islands >= int(config["model"]["pilot_min_islands"]):
                    support = "pilot"
                else:
                    support = "below_pilot"
                coefficients, fit = fit_weighted_linear_clustered(
                    context_data,
                    response_column=outcome,
                    weight_column=weights[outcome],
                    predictors=predictors,
                    cluster_column=cluster,
                )
                hit = coefficients.loc[
                    coefficients["predictor"].eq("log_distance_to_continent_km")
                ] if not coefficients.empty else pd.DataFrame()
                if hit.empty:
                    rows.append(
                        {
                            **meta,
                            "context": context,
                            "outcome": outcome,
                            "support_class": support,
                            "n_islands": n_islands,
                            "fit_status": str(fit.get("status", "unknown")),
                            "distance_slope": np.nan,
                            "cluster_robust_se": np.nan,
                            "p_value": np.nan,
                        }
                    )
                    continue
                row = hit.iloc[0]
                rows.append(
                    {
                        **meta,
                        "context": context,
                        "outcome": outcome,
                        "support_class": support,
                        "n_islands": n_islands,
                        "fit_status": str(fit.get("status", "fit")),
                        "distance_slope": float(row["estimate"]),
                        "cluster_robust_se": float(row["cluster_robust_se"]),
                        "p_value": float(row["p_value"]),
                    }
                )
    slopes = pd.DataFrame(rows)
    if slopes.empty:
        return slopes
    slopes["ci_low"] = slopes["distance_slope"] - 1.96 * slopes["cluster_robust_se"]
    slopes["ci_high"] = slopes["distance_slope"] + 1.96 * slopes["cluster_robust_se"]
    family = ["metric", "outcome", "context", "stratum", "source_matching"]
    slopes["q_value"] = slopes.groupby(family, group_keys=False)["p_value"].transform(_bh)
    return slopes


def classify_primary(slopes: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    primary = config["primary_prediction"]
    modes = [str(x) for x in config["source_matching"]["source_modes"]]
    primary_threshold = int(primary["effort_threshold"])
    sensitivities = [int(x) for x in config["sampling_effort_gate"]["sensitivity_min_independent_references"]]
    rows: list[dict[str, Any]] = []
    for context in [str(x) for x in config["model"]["contexts"]]:
        for stratum in [str(x) for x in config["source_matching"]["strata"]]:
            base = slopes.loc[
                slopes["metric"].astype(str).eq(str(primary["metric"]))
                & slopes["outcome"].astype(str).eq(str(primary["outcome"]))
                & slopes["source_matching"].astype(str).eq(str(primary["source_matching"]))
                & slopes["context"].astype(str).eq(context)
                & slopes["stratum"].astype(str).eq(stratum)
            ].copy()
            primary_rows = base.loc[base["min_independent_references"].eq(primary_threshold)]
            primary_modes = set(primary_rows["source_mode"].astype(str))
            primary_pass = (
                primary_modes == set(modes)
                and primary_rows["distance_slope"].gt(0).all()
                and primary_rows["q_value"].le(0.05).all()
                and primary_rows["support_class"].eq("confirmatory").all()
            )
            sensitivity_sign_ok = True
            sensitivity_cells = 0
            for threshold in sensitivities:
                part = base.loc[base["min_independent_references"].eq(threshold)]
                sensitivity_cells += len(part)
                if set(part["source_mode"].astype(str)) != set(modes) or not part["distance_slope"].gt(0).all():
                    sensitivity_sign_ok = False
            classification = (
                "robust_positive_source_breadth_enrichment"
                if primary_pass and sensitivity_sign_ok
                else "not_promoted"
            )
            rows.append(
                {
                    "context": context,
                    "stratum": stratum,
                    "primary_metric": str(primary["metric"]),
                    "primary_outcome": str(primary["outcome"]),
                    "primary_effort_threshold": primary_threshold,
                    "n_primary_source_modes": len(primary_modes),
                    "primary_all_positive": bool(len(primary_rows) and primary_rows["distance_slope"].gt(0).all()),
                    "primary_all_q_le_005": bool(len(primary_rows) and primary_rows["q_value"].le(0.05).all()),
                    "sensitivity_cells_checked": sensitivity_cells,
                    "sensitivity_no_sign_reversal": bool(sensitivity_sign_ok),
                    "classification": classification,
                }
            )
    return pd.DataFrame(rows)


def analyse_source_breadth(
    *,
    genus_breadth_csv: Path,
    gift_flora_csv: Path,
    assignments_csv: Path,
    status_flora_csv: Path,
    covariates_csv: Path,
    config_path: Path,
    output_dir: Path,
    predictor_sha256: str = "",
) -> dict[str, Any]:
    config = load_config(config_path)
    observed_sha = sha256_file(genus_breadth_csv)
    if predictor_sha256 and observed_sha != predictor_sha256:
        raise ValueError("genus-breadth predictor SHA changed after the outcome-blind build stage")
    breadth = pd.read_csv(genus_breadth_csv)
    gift = pd.read_csv(gift_flora_csv)
    assignments = pd.read_csv(assignments_csv)
    status = pd.read_csv(status_flora_csv)
    covariates = pd.read_csv(covariates_csv)
    enrichment = build_island_breadth_enrichment(
        genus_breadth=breadth,
        gift_flora=gift,
        assignments=assignments,
        status_flora=status,
        covariates=covariates,
        config=config,
    )
    slopes = fit_breadth_models(enrichment, covariates, config)
    classification = classify_primary(slopes, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    enrichment.to_csv(output_dir / "globi_source_breadth_island_enrichment.csv.gz", index=False, compression="gzip")
    slopes.to_csv(output_dir / "globi_source_breadth_context_slopes.csv", index=False)
    classification.to_csv(output_dir / "globi_source_breadth_primary_classification.csv", index=False)
    manifest = {
        "contract": config["contract"],
        "status": "secondary_H3_source_breadth_analysis_complete",
        "predictor_sha256": observed_sha,
        "n_genus_breadth_rows": int(len(breadth)),
        "n_island_enrichment_rows": int(len(enrichment)),
        "n_slope_rows": int(len(slopes)),
        "primary_promoted_cells": int(classification["classification"].eq("robust_positive_source_breadth_enrichment").sum()),
        "N1_rescued": False,
        "N2_opened": False,
        "claim_boundary": str(config["claim_ceiling"]),
    }
    (output_dir / "globi_source_breadth_analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("predictor")
def predictor_command(
    interactions_tsv_gz: Path = typer.Option(..., exists=True, dir_okay=False),
    refuted_tsv_gz: Path = typer.Option(..., exists=True, dir_okay=False),
    taxonomy_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    config_path: Path = typer.Option(Path("config/chapter1_globi_source_breadth.yml"), exists=True, dir_okay=False),
    channel_config_path: Path = typer.Option(Path("config/chapter1_nee_channel_taxon_catalog.yml"), exists=True, dir_okay=False),
    output_dir: Path = typer.Option(...),
    chunksize: int = typer.Option(200_000, min=1),
) -> None:
    typer.echo(
        json.dumps(
            build_genus_breadth_predictor(
                interactions_tsv_gz=interactions_tsv_gz,
                refuted_tsv_gz=refuted_tsv_gz,
                taxonomy_csv=taxonomy_csv,
                config_path=config_path,
                channel_config_path=channel_config_path,
                output_dir=output_dir,
                chunksize=chunksize,
            ),
            indent=2,
        )
    )


@app.command("analyse")
def analyse_command(
    genus_breadth_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    gift_flora_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    assignments_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    status_flora_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    covariates_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    config_path: Path = typer.Option(Path("config/chapter1_globi_source_breadth.yml"), exists=True, dir_okay=False),
    output_dir: Path = typer.Option(...),
    predictor_sha256: str = typer.Option(""),
) -> None:
    typer.echo(
        json.dumps(
            analyse_source_breadth(
                genus_breadth_csv=genus_breadth_csv,
                gift_flora_csv=gift_flora_csv,
                assignments_csv=assignments_csv,
                status_flora_csv=status_flora_csv,
                covariates_csv=covariates_csv,
                config_path=config_path,
                output_dir=output_dir,
                predictor_sha256=predictor_sha256,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    app()
