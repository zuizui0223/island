"""Fail-closed floristic-status support for Chapter 1 analyses.

Floristic status and trait evidence are separate evidence layers. Status is never
inferred from GBIF occurrence footprint. Missing or conflicting status remains
unresolved and is kept visible in downstream attrition reports.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

DIRECT_SCOPES = {"species_direct", "synonym_direct"}
ORIGIN_VALUES = {"native", "introduced", "unresolved"}
ENDEMISM_VALUES = {"endemic", "nonendemic", "unresolved"}
STRATA = ("all_native", "native_nonendemic", "endemic", "introduced", "unresolved")


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def normalize_island_species(island_species: pd.DataFrame) -> pd.DataFrame:
    frame = island_species.copy()
    if "accepted_species" not in frame.columns and "species" in frame.columns:
        frame = frame.rename(columns={"species": "accepted_species"})
    required = {"island_id", "accepted_species"}
    missing = required - set(frame.columns)
    if missing:
        raise typer.BadParameter(f"island species table missing columns: {sorted(missing)}")
    frame = frame[["island_id", "accepted_species"]].copy()
    frame["island_id"] = _text(frame["island_id"])
    frame["accepted_species"] = _text(frame["accepted_species"])
    frame = frame.loc[frame["island_id"].ne("") & frame["accepted_species"].ne("")]
    return frame.drop_duplicates().sort_values(["island_id", "accepted_species"]).reset_index(drop=True)


def collapse_status_ledger(status: pd.DataFrame) -> pd.DataFrame:
    required = {"island_id", "accepted_species", "origin_status", "endemic_status"}
    missing = required - set(status.columns)
    if missing:
        raise typer.BadParameter(f"status ledger missing columns: {sorted(missing)}")

    frame = status.copy()
    frame["island_id"] = _text(frame["island_id"])
    frame["accepted_species"] = _text(frame["accepted_species"])
    frame["origin_status"] = _text(frame["origin_status"]).str.lower().replace("", "unresolved")
    frame["endemic_status"] = _text(frame["endemic_status"]).str.lower().replace("", "unresolved")

    invalid_origin = sorted(set(frame["origin_status"]) - ORIGIN_VALUES)
    invalid_endemism = sorted(set(frame["endemic_status"]) - ENDEMISM_VALUES)
    if invalid_origin:
        raise typer.BadParameter(f"invalid origin_status values: {invalid_origin}")
    if invalid_endemism:
        raise typer.BadParameter(f"invalid endemic_status values: {invalid_endemism}")

    optional = [
        c for c in ("status_source", "status_reference", "status_evidence_id") if c in frame.columns
    ]
    rows: list[dict[str, Any]] = []
    for (island_id, species), group in frame.groupby(["island_id", "accepted_species"], sort=True):
        origins = {v for v in group["origin_status"] if v != "unresolved"}
        endemics = {v for v in group["endemic_status"] if v != "unresolved"}
        conflict = len(origins) > 1 or len(endemics) > 1
        origin = next(iter(origins)) if len(origins) == 1 and not conflict else "unresolved"
        endemic = next(iter(endemics)) if len(endemics) == 1 and not conflict else "unresolved"
        row: dict[str, Any] = {
            "island_id": island_id,
            "accepted_species": species,
            "origin_status": origin,
            "endemic_status": endemic,
            "status_conflict": conflict,
            "n_status_rows": int(len(group)),
        }
        for column in optional:
            row[column] = "|".join(sorted({v for v in _text(group[column]) if v}))
        rows.append(row)
    return pd.DataFrame(rows)


def attach_floristic_status(island_species: pd.DataFrame, status_ledger: pd.DataFrame) -> pd.DataFrame:
    flora = normalize_island_species(island_species)
    status = collapse_status_ledger(status_ledger)
    out = flora.merge(status, on=["island_id", "accepted_species"], how="left", validate="one_to_one")
    out["origin_status"] = _text(out["origin_status"]).replace("", "unresolved")
    out["endemic_status"] = _text(out["endemic_status"]).replace("", "unresolved")
    out["status_conflict"] = out["status_conflict"].eq(True)
    out["status_resolved"] = out["origin_status"].ne("unresolved")

    def combined(row: pd.Series) -> str:
        if row["origin_status"] == "introduced":
            return "introduced"
        if row["origin_status"] != "native":
            return "unresolved"
        if row["endemic_status"] == "endemic":
            return "endemic"
        if row["endemic_status"] == "nonendemic":
            return "native_nonendemic"
        return "native_endemism_unresolved"

    out["floristic_status"] = out.apply(combined, axis=1)
    return out.sort_values(["island_id", "accepted_species"]).reset_index(drop=True)


def stratum_mask(frame: pd.DataFrame, stratum: str) -> pd.Series:
    if stratum == "all_native":
        return frame["origin_status"].eq("native")
    if stratum == "native_nonendemic":
        return frame["floristic_status"].eq("native_nonendemic")
    if stratum == "endemic":
        return frame["floristic_status"].eq("endemic")
    if stratum == "introduced":
        return frame["floristic_status"].eq("introduced")
    if stratum == "unresolved":
        return frame["floristic_status"].isin({"unresolved", "native_endemism_unresolved"})
    raise typer.BadParameter(f"unknown stratum: {stratum}")


def direct_species_by_outcome(evidence: pd.DataFrame, outcomes: dict[str, Any]) -> dict[str, set[str]]:
    required = {"accepted_species", "trait_name"}
    missing = required - set(evidence.columns)
    if missing:
        raise typer.BadParameter(f"direct evidence table missing columns: {sorted(missing)}")
    frame = evidence.copy()
    frame["accepted_species"] = _text(frame["accepted_species"])
    frame["trait_name"] = _text(frame["trait_name"])
    if "evidence_scope" in frame.columns:
        frame = frame.loc[_text(frame["evidence_scope"]).str.lower().isin(DIRECT_SCOPES)].copy()
    if "resolution_status" in frame.columns:
        frame = frame.loc[_text(frame["resolution_status"]).str.lower().eq("resolved")].copy()

    result: dict[str, set[str]] = {}
    for outcome, spec in outcomes.items():
        traits = spec.get("trait_names") or [spec.get("trait_name")]
        names = {str(v) for v in traits if v}
        if not names:
            raise typer.BadParameter(f"outcome {outcome} declares no trait_name(s)")
        result[outcome] = set(frame.loc[frame["trait_name"].isin(names), "accepted_species"])
    return result


def build_direct_support(
    status_flora: pd.DataFrame,
    direct_by_outcome: dict[str, set[str]],
    covariates: pd.DataFrame,
    *,
    min_direct_species: int = 1,
) -> pd.DataFrame:
    if "island_id" not in covariates.columns:
        raise typer.BadParameter("covariates must contain island_id")
    cov = covariates.drop_duplicates("island_id").copy()
    cov["island_id"] = _text(cov["island_id"])
    rows: list[pd.DataFrame] = []
    for stratum in STRATA:
        subset = status_flora.loc[stratum_mask(status_flora, stratum)].copy()
        for outcome, covered in direct_by_outcome.items():
            work = subset[["island_id", "accepted_species"]].copy()
            work["has_direct"] = work["accepted_species"].isin(covered)
            counts = work.groupby("island_id", as_index=False).agg(
                n_stratum_species=("accepted_species", "size"),
                n_direct_species=("has_direct", "sum"),
            )
            counts["n_direct_species"] = counts["n_direct_species"].astype(int)
            counts["direct_fraction"] = counts["n_direct_species"] / counts["n_stratum_species"].replace(0, np.nan)
            counts["direct_eligible"] = counts["n_direct_species"].ge(min_direct_species)
            counts["stratum"] = stratum
            counts["outcome"] = outcome
            rows.append(counts)
    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if out.empty:
        return out
    return out.merge(cov, on="island_id", how="left", validate="many_to_one")


def _load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not payload.get("outcomes"):
        raise typer.BadParameter("config must contain outcomes")
    return payload


@app.command("run")
def run(
    island_species_csv: Path = typer.Option(..., exists=True),
    status_ledger_csv: Path = typer.Option(..., exists=True),
    direct_evidence_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/flora_status_support.yml"), exists=True),
) -> None:
    config = _load_config(config_path)
    status_flora = attach_floristic_status(pd.read_csv(island_species_csv), pd.read_csv(status_ledger_csv))
    direct = direct_species_by_outcome(pd.read_csv(direct_evidence_csv), config["outcomes"])
    support = build_direct_support(
        status_flora,
        direct,
        pd.read_csv(covariates_csv),
        min_direct_species=int(config.get("support", {}).get("min_direct_species_per_island", 1)),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    status_flora.to_csv(output_dir / "island_species_floristic_status.csv.gz", index=False)
    support.to_csv(output_dir / "direct_trait_support_by_island.csv.gz", index=False)
    manifest = {
        "contract": "flora_status_support_v2",
        "status_policy": "source-backed only; unresolved/conflicting status fails closed",
        "n_island_species_rows": int(len(status_flora)),
        "n_status_resolved_rows": int(status_flora["status_resolved"].sum()),
        "n_status_conflicts": int(status_flora["status_conflict"].sum()),
        "outcomes": sorted(direct),
    }
    (output_dir / "flora_status_support_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
