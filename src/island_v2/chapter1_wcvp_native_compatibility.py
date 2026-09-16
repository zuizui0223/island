"""WCVP TDWG-L3 native-compatibility sensitivity for the all-data Chapter 1 route.

This does not claim exact island-native status. It upgrades unresolved island-species
records only when the focal island maps unambiguously to a TDWG level-3 unit that WCVP
lists in the accepted species' native range. Source-backed GIFT `introduced` rows remain
introduced and are never overwritten. Existing source-backed `native` rows remain native.
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


def build_regional_native_compatible_flora(
    status_flora: pd.DataFrame,
    wcvp_ranges: pd.DataFrame,
    island_tdwg: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    flora = status_flora.copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    flora["origin_status"] = flora["origin_status"].fillna("unresolved").astype(str)
    ranges = wcvp_ranges[["accepted_species", "native_l3_codes"]].drop_duplicates("accepted_species").copy()
    ranges["accepted_species"] = ranges["accepted_species"].astype(str)
    mapping = island_tdwg[["island_id", "tdwg_l3_code", "tdwg_match_status"]].drop_duplicates("island_id").copy()
    mapping["island_id"] = mapping["island_id"].astype(str)
    merged = flora.merge(ranges, on="accepted_species", how="left", validate="many_to_one").merge(
        mapping, on="island_id", how="left", validate="many_to_one"
    )
    merged["native_l3_codes"] = merged["native_l3_codes"].fillna("").astype(str)
    merged["tdwg_l3_code"] = merged["tdwg_l3_code"].fillna("").astype(str)
    merged["tdwg_match_status"] = merged["tdwg_match_status"].fillna("").astype(str)
    compatible = [
        bool(code) and status == "accepted" and code in _tokens(native_codes)
        for code, status, native_codes in zip(
            merged["tdwg_l3_code"], merged["tdwg_match_status"], merged["native_l3_codes"], strict=True
        )
    ]
    merged["wcvp_regional_native_compatible"] = compatible
    unresolved = merged["origin_status"].eq("unresolved")
    upgraded = unresolved & merged["wcvp_regional_native_compatible"]
    merged["origin_status_original"] = merged["origin_status"]
    merged.loc[upgraded, "origin_status"] = "native"
    # Preserve exact native/endemic labels where known. WCVP-compatible upgrades are only
    # native at regional resolution and therefore cannot be declared endemic/nonendemic.
    if "floristic_status" not in merged.columns:
        merged["floristic_status"] = "unresolved"
    merged.loc[upgraded, "floristic_status"] = "native_endemism_unresolved"
    audit = {
        "n_rows": int(len(merged)),
        "n_original_native": int(merged["origin_status_original"].eq("native").sum()),
        "n_original_introduced": int(merged["origin_status_original"].eq("introduced").sum()),
        "n_original_unresolved": int(merged["origin_status_original"].eq("unresolved").sum()),
        "n_wcvp_compatible_upgrades": int(upgraded.sum()),
        "n_islands_with_upgrade": int(merged.loc[upgraded, "island_id"].nunique()),
        "n_regional_native_rows": int(merged["origin_status"].eq("native").sum()),
        "n_regional_native_islands": int(merged.loc[merged["origin_status"].eq("native"), "island_id"].nunique()),
    }
    return merged, audit


def run_compatibility_analysis(
    status_flora: pd.DataFrame,
    state_audit: pd.DataFrame,
    covariates: pd.DataFrame,
    wcvp_ranges: pd.DataFrame,
    island_tdwg: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    upgraded, audit = build_regional_native_compatible_flora(status_flora, wcvp_ranges, island_tdwg)
    cfg = dict(config)
    cfg["strata"] = ["all_native"]
    counts = build_broad_counts(upgraded, state_audit, cfg)
    counts["stratum"] = "regional_native_compatible"
    cfg["strata"] = ["regional_native_compatible"]
    # run_probability_analysis only needs stratum labels in the count table/config.
    within_slopes, between_slopes, within, between = run_probability_analysis(counts, covariates, cfg)
    return within_slopes, between_slopes, within, between, audit


@app.command("run")
def run(
    status_flora_csv: Path = typer.Option(..., exists=True),
    state_audit_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    wcvp_ranges_csv: Path = typer.Option(..., exists=True),
    island_tdwg_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
) -> None:
    cfg = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    within_slopes, between_slopes, within, between, audit = run_compatibility_analysis(
        pd.read_csv(status_flora_csv), pd.read_csv(state_audit_csv), pd.read_csv(covariates_csv),
        pd.read_csv(wcvp_ranges_csv), pd.read_csv(island_tdwg_csv), cfg,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    for frame in (within_slopes, between_slopes, within, between):
        if not frame.empty:
            frame.insert(0, "evidence_scope", evidence_scope)
    within_slopes.to_csv(output_dir / "regional_native_within_slopes.csv", index=False)
    between_slopes.to_csv(output_dir / "regional_native_between_slopes.csv", index=False)
    within.to_csv(output_dir / "regional_native_within_omnibus.csv", index=False)
    between.to_csv(output_dir / "regional_native_between_omnibus.csv", index=False)
    pd.DataFrame([{"evidence_scope": evidence_scope, **audit}]).to_csv(output_dir / "regional_native_coverage.csv", index=False)
    typer.echo(pd.DataFrame([audit]).to_csv(index=False))
    typer.echo(between.to_csv(index=False))


if __name__ == "__main__":
    app()
