"""Partition the broad observed flora by source-backed and WCVP native compatibility.

This is a diagnostic for the all-data Chapter 1 route. It does not relabel WCVP-incompatible
or unclassifiable records as introduced. Instead it asks which evidence partition carries the
North--Tropical atomic-trait response difference.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.chapter1_all_data_probability import build_broad_counts, run_probability_analysis

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _tokens(value: object) -> set[str]:
    return {x.strip() for x in str(value or "").split("|") if x.strip()}


def classify_wcvp_partitions(
    status_flora: pd.DataFrame,
    wcvp_ranges: pd.DataFrame,
    island_tdwg: pd.DataFrame,
) -> pd.DataFrame:
    flora = status_flora.copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    flora["origin_status"] = flora["origin_status"].fillna("unresolved").astype(str)

    ranges = wcvp_ranges[["accepted_species", "native_l3_codes"]].drop_duplicates(
        "accepted_species"
    )
    ranges["accepted_species"] = ranges["accepted_species"].astype(str)
    mapping = island_tdwg[
        ["island_id", "tdwg_l3_code", "tdwg_match_status"]
    ].drop_duplicates("island_id")
    mapping["island_id"] = mapping["island_id"].astype(str)

    merged = flora.merge(ranges, on="accepted_species", how="left", validate="many_to_one")
    merged = merged.merge(mapping, on="island_id", how="left", validate="many_to_one")
    merged["native_l3_codes"] = merged["native_l3_codes"].fillna("").astype(str)
    merged["tdwg_l3_code"] = merged["tdwg_l3_code"].fillna("").astype(str)
    merged["tdwg_match_status"] = merged["tdwg_match_status"].fillna("").astype(str)

    partitions: list[str] = []
    for origin, code, match_status, native_codes in zip(
        merged["origin_status"],
        merged["tdwg_l3_code"],
        merged["tdwg_match_status"],
        merged["native_l3_codes"],
        strict=True,
    ):
        if origin == "native":
            partitions.append("source_native")
            continue
        if origin == "introduced":
            partitions.append("source_introduced")
            continue
        codes = _tokens(native_codes)
        if match_status == "accepted" and code and codes:
            if code in codes:
                partitions.append("wcvp_compatible_unresolved")
            else:
                partitions.append("wcvp_incompatible_unresolved")
        else:
            partitions.append("wcvp_unclassifiable_unresolved")
    merged["wcvp_partition"] = partitions
    return merged


def _run_one_partition(
    flora: pd.DataFrame,
    state_audit: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> dict[str, Any]:
    cfg = dict(config)
    cfg["strata"] = ["all_observed"]
    cfg["contexts"] = ["northern_midlatitude", "tropical"]
    cfg["between_contexts"] = [["northern_midlatitude", "tropical"]]
    counts = build_broad_counts(flora, state_audit, cfg)
    _, _, _, between = run_probability_analysis(counts, covariates, cfg)
    result: dict[str, Any] = {
        "n_flora_rows": int(len(flora)),
        "n_islands": int(flora["island_id"].nunique()),
    }
    if between.empty:
        result.update({"status": "not_testable", "p_value": float("nan")})
        return result
    row = between.iloc[0]
    result.update(
        {
            "status": str(row["status"]),
            "n_analysis_islands": row.get("n_unique_islands"),
            "n_clusters": row.get("n_clusters"),
            "n_retained_outcomes": row.get("n_retained_outcomes"),
            "p_value": row.get("p_value"),
        }
    )
    return result


def run_partition_diagnostic(
    status_flora: pd.DataFrame,
    state_audit: pd.DataFrame,
    covariates: pd.DataFrame,
    wcvp_ranges: pd.DataFrame,
    island_tdwg: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    classified = classify_wcvp_partitions(status_flora, wcvp_ranges, island_tdwg)
    masks = {
        "source_native": classified["wcvp_partition"].eq("source_native"),
        "source_introduced": classified["wcvp_partition"].eq("source_introduced"),
        "wcvp_compatible_unresolved": classified["wcvp_partition"].eq(
            "wcvp_compatible_unresolved"
        ),
        "wcvp_incompatible_unresolved": classified["wcvp_partition"].eq(
            "wcvp_incompatible_unresolved"
        ),
        "wcvp_unclassifiable_unresolved": classified["wcvp_partition"].eq(
            "wcvp_unclassifiable_unresolved"
        ),
        "regional_native_compatible": classified["wcvp_partition"].isin(
            {"source_native", "wcvp_compatible_unresolved"}
        ),
        "regionally_incompatible_or_introduced": classified["wcvp_partition"].isin(
            {"source_introduced", "wcvp_incompatible_unresolved"}
        ),
    }
    rows: list[dict[str, Any]] = []
    for label, mask in masks.items():
        part = classified.loc[mask].copy()
        result = _run_one_partition(part, state_audit, covariates, config)
        rows.append({"partition": label, **result})
    return pd.DataFrame(rows)


@app.command("run")
def run(
    status_flora_csv: Path = typer.Option(..., exists=True),
    state_audit_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    wcvp_ranges_csv: Path = typer.Option(..., exists=True),
    island_tdwg_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_csv: Path = typer.Option(...),
) -> None:
    result = run_partition_diagnostic(
        pd.read_csv(status_flora_csv),
        pd.read_csv(state_audit_csv),
        pd.read_csv(covariates_csv),
        pd.read_csv(wcvp_ranges_csv),
        pd.read_csv(island_tdwg_csv),
        yaml.safe_load(config_path.read_text(encoding="utf-8")),
    )
    result.insert(0, "evidence_scope", evidence_scope)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_csv, index=False)
    typer.echo(result.to_csv(index=False))


if __name__ == "__main__":
    app()
