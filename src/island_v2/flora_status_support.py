"""Status-aware trait-support audit for island floral analyses.

This module deliberately separates two problems that were previously conflated:
(1) floristic status (native, endemic, introduced, unresolved) and
(2) whether a species has direct source-backed evidence for a focal trait.

It never infers floristic status from GBIF occurrence footprint. Status must be
provided by an explicit flora/checklist source. Missing or conflicting status is
retained as unresolved.
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


def _load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    required = {"outcomes", "support"}
    if not isinstance(payload, dict) or not required.issubset(payload):
        raise typer.BadParameter("config must contain outcomes and support")
    if not isinstance(payload["outcomes"], dict) or not payload["outcomes"]:
        raise typer.BadParameter("outcomes must be a non-empty mapping")
    return payload


def normalize_island_species(island_species: pd.DataFrame) -> pd.DataFrame:
    """Return one row per island x accepted species without inventing status."""
    frame = island_species.copy()
    if "accepted_species" not in frame.columns and "species" in frame.columns:
        frame = frame.rename(columns={"species": "accepted_species"})
    required = {"island_id", "accepted_species"}
    missing = required - set(frame.columns)
    if missing:
        raise typer.BadParameter(f"island species table missing columns: {sorted(missing)}")
    frame = frame.loc[:, ["island_id", "accepted_species"]].copy()
    frame["island_id"] = _text(frame["island_id"])
    frame["accepted_species"] = _text(frame["accepted_species"])
    frame = frame.loc[
        frame["island_id"].ne("") & frame["accepted_species"].ne("")
    ].drop_duplicates()
    return frame.sort_values(["island_id", "accepted_species"]).reset_index(drop=True)


def collapse_status_ledger(status: pd.DataFrame) -> pd.DataFrame:
    """Collapse source-backed status rows; any conflict fails closed to unresolved.

    Required semantic columns are origin_status and endemic_status. Sources may
    contribute multiple rows for the same island x species. Agreement is retained;
    disagreement is represented explicitly instead of resolved by hidden priority.
    """
    required = {"island_id", "accepted_species", "origin_status", "endemic_status"}
    missing = required - set(status.columns)
    if missing:
        raise typer.BadParameter(f"status ledger missing columns: {sorted(missing)}")

    frame = status.copy()
    for column in required:
        frame[column] = (
            _text(frame[column]).str.lower() if column.endswith("status") else _text(frame[column])
        )
    invalid_origin = sorted(set(frame["origin_status"]) - ORIGIN_VALUES)
    invalid_endemism = sorted(set(frame["endemic_status"]) - ENDEMISM_VALUES)
    if invalid_origin:
        raise typer.BadParameter(f"invalid origin_status values: {invalid_origin}")
    if invalid_endemism:
        raise typer.BadParameter(f"invalid endemic_status values: {invalid_endemism}")

    optional = [
        column
        for column in ("status_source", "status_reference", "status_evidence_id")
        if column in frame.columns
    ]
    rows: list[dict[str, Any]] = []
    for (island_id, species), group in frame.groupby(
        ["island_id", "accepted_species"], sort=True
    ):
        origins = {value for value in group["origin_status"] if value != "unresolved"}
        endemics = {value for value in group["endemic_status"] if value != "unresolved"}
        origin = next(iter(origins)) if len(origins) == 1 else "unresolved"
        endemic = next(iter(endemics)) if len(endemics) == 1 else "unresolved"
        conflict = len(origins) > 1 or len(endemics) > 1
        if conflict:
            origin = "unresolved"
            endemic = "unresolved"
        row: dict[str, Any] = {
            "island_id": island_id,
            "accepted_species": species,
            "origin_status": origin,
            "endemic_status": endemic,
            "status_conflict": conflict,
            "n_status_rows": int(len(group)),
        }
        for column in optional:
            values = sorted({value for value in _text(group[column]) if value})
            row[column] = "|".join(values)
        rows.append(row)
    return pd.DataFrame(rows)


def attach_floristic_status(
    island_species: pd.DataFrame, status_ledger: pd.DataFrame
) -> pd.DataFrame:
    """Attach status to the entire flora and retain unmatched rows as unresolved."""
    flora = normalize_island_species(island_species)
    status = collapse_status_ledger(status_ledger)
    out = flora.merge(
        status,
        on=["island_id", "accepted_species"],
        how="left",
        validate="one_to_one",
    )
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


def direct_species_by_outcome(
    evidence: pd.DataFrame, outcomes: dict[str, Any]
) -> dict[str, set[str]]:
    """Return species with direct evidence for each declared final outcome."""
    required = {"accepted_species", "trait_name"}
    missing = required - set(evidence.columns)
    if missing:
        raise typer.BadParameter(f"direct evidence table missing columns: {sorted(missing)}")
    frame = evidence.copy()
    frame["accepted_species"] = _text(frame["accepted_species"])
    frame["trait_name"] = _text(frame["trait_name"])
    if "evidence_scope" in frame.columns:
        scope = _text(frame["evidence_scope"]).str.lower()
        frame = frame.loc[scope.isin(DIRECT_SCOPES)].copy()
    if "resolution_status" in frame.columns:
        resolution = _text(frame["resolution_status"]).str.lower()
        frame = frame.loc[resolution.eq("resolved")].copy()

    result: dict[str, set[str]] = {}
    for outcome, spec in outcomes.items():
        traits = spec.get("trait_names") or [spec.get("trait_name")]
        traits = {str(value) for value in traits if value}
        if not traits:
            raise typer.BadParameter(f"outcome {outcome} declares no trait_name(s)")
        rows = frame.loc[frame["trait_name"].isin(traits)]
        result[outcome] = set(rows["accepted_species"])
    return result


def island_status_counts(status_flora: pd.DataFrame) -> pd.DataFrame:
    """Build island-level native/endemic/introduced counts for endemism models."""
    rows: list[dict[str, Any]] = []
    for island_id, group in status_flora.groupby("island_id", sort=True):
        native = group["origin_status"].eq("native")
        endemic = group["floristic_status"].eq("endemic")
        nonendemic = group["floristic_status"].eq("native_nonendemic")
        introduced = group["floristic_status"].eq("introduced")
        unresolved = ~group["floristic_status"].isin(
            {"endemic", "native_nonendemic", "introduced"}
        )
        n_native = int(native.sum())
        n_endemic = int(endemic.sum())
        rows.append(
            {
                "island_id": island_id,
                "n_flora_species": int(len(group)),
                "n_native_species": n_native,
                "n_native_nonendemic_species": int(nonendemic.sum()),
                "n_endemic_species": n_endemic,
                "n_introduced_species": int(introduced.sum()),
                "n_unresolved_status_species": int(unresolved.sum()),
                "native_status_fraction": float(native.mean()) if len(group) else np.nan,
                "endemic_share_of_native": (
                    float(n_endemic / n_native) if n_native > 0 else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def build_direct_support(
    status_flora: pd.DataFrame,
    direct_by_outcome: dict[str, set[str]],
    covariates: pd.DataFrame,
    *,
    min_direct_species: int,
    min_direct_fraction: float = 0.0,
) -> pd.DataFrame:
    """Count direct evidence per island x status stratum x final outcome."""
    if "island_id" not in covariates.columns:
        raise typer.BadParameter("covariates must contain island_id")
    cov = covariates.drop_duplicates("island_id").copy()
    cov["island_id"] = _text(cov["island_id"])
    rows: list[dict[str, Any]] = []
    for stratum in STRATA:
        subset = status_flora.loc[stratum_mask(status_flora, stratum)].copy()
        for outcome, covered in direct_by_outcome.items():
            work = subset.copy()
            work["has_direct"] = work["accepted_species"].isin(covered)
            counts = (
                work.groupby("island_id", as_index=False)
                .agg(
                    n_stratum_species=("accepted_species", "size"),
                    n_direct_species=("has_direct", "sum"),
                )
            )
            counts["n_direct_species"] = counts["n_direct_species"].astype(int)
            counts["direct_fraction"] = (
                counts["n_direct_species"] / counts["n_stratum_species"].replace(0, np.nan)
            )
            counts["direct_eligible"] = (
                counts["n_direct_species"].ge(min_direct_species)
                & counts["direct_fraction"].ge(min_direct_fraction)
            )
            counts["gap_to_min_direct_species"] = (
                min_direct_species - counts["n_direct_species"]
            ).clip(lower=0)
            counts["stratum"] = stratum
            counts["outcome"] = outcome
            rows.append(counts)
    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if out.empty:
        return out
    return out.merge(cov, on="island_id", how="left", validate="many_to_one")


def summarize_support(
    island_support: pd.DataFrame,
    *,
    regime_column: str,
    isolation_column: str,
    pilot_min_islands: int,
    confirmatory_min_islands: int,
) -> pd.DataFrame:
    """Summarize island counts and isolation support without inventing a gate."""
    required = {"stratum", "outcome", "direct_eligible", isolation_column}
    missing = required - set(island_support.columns)
    if missing:
        raise typer.BadParameter(f"island support missing columns: {sorted(missing)}")
    work = island_support.copy()
    if regime_column not in work.columns:
        work[regime_column] = "all"
    work[regime_column] = _text(work[regime_column]).replace("", "unresolved")
    work[isolation_column] = pd.to_numeric(work[isolation_column], errors="coerce")

    rows: list[dict[str, Any]] = []
    for (regime, stratum, outcome), group in work.groupby(
        [regime_column, "stratum", "outcome"], sort=True
    ):
        eligible = group.loc[group["direct_eligible"]].copy()
        iso = eligible[isolation_column].dropna().to_numpy(dtype=float)
        n_eligible = int(len(eligible))
        if n_eligible >= confirmatory_min_islands:
            class_ = "confirmatory_count_met"
        elif n_eligible >= pilot_min_islands:
            class_ = "targeted_acquisition_zone"
        else:
            class_ = "below_pilot_support"
        record: dict[str, Any] = {
            "regime": regime,
            "stratum": stratum,
            "outcome": outcome,
            "n_islands_with_stratum": int(len(group)),
            "n_direct_eligible_islands": n_eligible,
            "support_class": class_,
            "median_direct_fraction": float(group["direct_fraction"].median()),
            "median_direct_species": float(group["n_direct_species"].median()),
        }
        if len(iso):
            q05, q25, q50, q75, q95 = np.quantile(iso, [0.05, 0.25, 0.5, 0.75, 0.95])
            record.update(
                {
                    "isolation_min": float(np.min(iso)),
                    "isolation_q05": float(q05),
                    "isolation_q25": float(q25),
                    "isolation_median": float(q50),
                    "isolation_q75": float(q75),
                    "isolation_q95": float(q95),
                    "isolation_max": float(np.max(iso)),
                }
            )
        else:
            for name in (
                "isolation_min",
                "isolation_q05",
                "isolation_q25",
                "isolation_median",
                "isolation_q75",
                "isolation_q95",
                "isolation_max",
            ):
                record[name] = np.nan
        if "spatial_block" in eligible.columns:
            record["n_spatial_blocks"] = int(eligible["spatial_block"].dropna().nunique())
        rows.append(record)
    return pd.DataFrame(rows)


def acquisition_candidates(
    status_flora: pd.DataFrame,
    island_support: pd.DataFrame,
    direct_by_outcome: dict[str, set[str]],
    *,
    min_direct_species: int,
    confirmatory_min_islands: int,
    regime_column: str,
    isolation_column: str,
) -> pd.DataFrame:
    """Rank missing direct evidence only where model support is below target.

    A species receives credit for occurring on islands that are currently below
    the per-island direct-species minimum. `one_record_unlocks_islands` counts the
    highest-leverage case: islands sitting exactly one species below eligibility.
    """
    summaries = summarize_support(
        island_support,
        regime_column=regime_column,
        isolation_column=isolation_column,
        pilot_min_islands=0,
        confirmatory_min_islands=confirmatory_min_islands,
    )
    needs = summaries.loc[
        summaries["n_direct_eligible_islands"].lt(confirmatory_min_islands),
        ["regime", "stratum", "outcome"],
    ]
    if needs.empty:
        return pd.DataFrame()

    support_lookup = island_support.copy()
    support_lookup[regime_column] = _text(support_lookup[regime_column]).replace(
        "", "unresolved"
    )
    support_lookup[isolation_column] = pd.to_numeric(
        support_lookup[isolation_column], errors="coerce"
    )
    rows: list[dict[str, Any]] = []
    for need in needs.itertuples(index=False):
        stratum = str(need.stratum)
        outcome = str(need.outcome)
        regime = str(need.regime)
        covered = direct_by_outcome[outcome]
        target_support = support_lookup.loc[
            support_lookup["stratum"].eq(stratum)
            & support_lookup["outcome"].eq(outcome)
            & support_lookup[regime_column].eq(regime)
            & ~support_lookup["direct_eligible"]
            & support_lookup["gap_to_min_direct_species"].between(1, min_direct_species),
            ["island_id", "gap_to_min_direct_species", isolation_column],
        ].copy()
        if target_support.empty:
            continue
        candidates = status_flora.loc[
            stratum_mask(status_flora, stratum)
            & ~status_flora["accepted_species"].isin(covered)
            & status_flora["island_id"].isin(target_support["island_id"]),
            ["island_id", "accepted_species"],
        ].drop_duplicates()
        candidates = candidates.merge(target_support, on="island_id", how="inner")
        if candidates.empty:
            continue
        candidates["gap_weight"] = 1.0 / candidates["gap_to_min_direct_species"].astype(float)
        candidates["one_unlock"] = candidates["gap_to_min_direct_species"].eq(1)
        for species, group in candidates.groupby("accepted_species", sort=True):
            iso = group[isolation_column].dropna().to_numpy(dtype=float)
            rows.append(
                {
                    "regime": regime,
                    "stratum": stratum,
                    "outcome": outcome,
                    "accepted_species": species,
                    "n_target_islands": int(group["island_id"].nunique()),
                    "one_record_unlocks_islands": int(
                        group.loc[group["one_unlock"], "island_id"].nunique()
                    ),
                    "gap_weighted_support_gain": float(group["gap_weight"].sum()),
                    "target_isolation_min": float(np.min(iso)) if len(iso) else np.nan,
                    "target_isolation_max": float(np.max(iso)) if len(iso) else np.nan,
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(
        ["one_record_unlocks_islands", "gap_weighted_support_gain", "n_target_islands"],
        ascending=[False, False, False],
    ).reset_index(drop=True)


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
    island_species = pd.read_csv(island_species_csv)
    status_ledger = pd.read_csv(status_ledger_csv)
    direct_evidence = pd.read_csv(direct_evidence_csv)
    covariates = pd.read_csv(covariates_csv)

    status_flora = attach_floristic_status(island_species, status_ledger)
    direct_by_outcome = direct_species_by_outcome(direct_evidence, config["outcomes"])
    support_cfg = config["support"]
    island_support = build_direct_support(
        status_flora,
        direct_by_outcome,
        covariates,
        min_direct_species=int(support_cfg["min_direct_species_per_island"]),
        min_direct_fraction=float(support_cfg.get("min_direct_fraction", 0.0)),
    )
    summary = summarize_support(
        island_support,
        regime_column=str(support_cfg["regime_column"]),
        isolation_column=str(support_cfg["isolation_column"]),
        pilot_min_islands=int(support_cfg["pilot_min_islands"]),
        confirmatory_min_islands=int(support_cfg["confirmatory_min_islands"]),
    )
    candidates = acquisition_candidates(
        status_flora,
        island_support,
        direct_by_outcome,
        min_direct_species=int(support_cfg["min_direct_species_per_island"]),
        confirmatory_min_islands=int(support_cfg["confirmatory_min_islands"]),
        regime_column=str(support_cfg["regime_column"]),
        isolation_column=str(support_cfg["isolation_column"]),
    )
    endemism = island_status_counts(status_flora).merge(
        covariates.drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="one_to_one",
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    status_flora.to_csv(output_dir / "island_species_floristic_status.csv.gz", index=False)
    endemism.to_csv(output_dir / "endemism_response_input.csv", index=False)
    island_support.to_csv(output_dir / "direct_trait_support_by_island.csv.gz", index=False)
    summary.to_csv(output_dir / "direct_trait_support_summary.csv", index=False)
    candidates.to_csv(output_dir / "targeted_acquisition_candidates.csv", index=False)

    manifest = {
        "contract": "flora_status_support_v1",
        "status_policy": "source-backed only; no endemic/native inference from GBIF footprint",
        "n_island_species_rows": int(len(status_flora)),
        "n_status_resolved_rows": int(status_flora["status_resolved"].sum()),
        "n_status_conflicts": int(status_flora["status_conflict"].sum()),
        "outcomes": sorted(direct_by_outcome),
        "strata": list(STRATA),
        "support": support_cfg,
        "n_targeted_acquisition_rows": int(len(candidates)),
    }
    (output_dir / "flora_status_support_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
