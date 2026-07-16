"""Audit and connect the rectangular, evidence-grounded trait cascade."""

from __future__ import annotations

import glob
import json
from pathlib import Path

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_COLUMNS = {
    "accepted_species",
    "trait_name",
    "filled_value",
    "fill_tier",
    "evidence_scope",
    "confidence",
}
GROUNDED_ANALYSIS_TIERS = {
    "species_direct",
    "synonym_direct",
    "genus_inference",
    "family_inference",
}


def _analysis_eligibility_mask(frame: pd.DataFrame) -> pd.Series:
    """Return a fail-closed mask for both current and historical artifacts."""
    tier_supported = frame["fill_tier"].astype(str).isin(GROUNDED_ANALYSIS_TIERS)
    if "analysis_eligible" not in frame.columns:
        return tier_supported
    return (
        frame["analysis_eligible"].astype(str).str.lower().eq("true")
        & tier_supported
    )


def _load_shards(root: Path) -> pd.DataFrame:
    paths = sorted(Path(path) for path in glob.glob(str(root / "shards" / "shard_*" / "fills.csv.gz")))
    if len(paths) != 128:
        raise typer.BadParameter(f"expected 128 fill shards, found {len(paths)}")
    frames = [pd.read_csv(path, dtype=str).fillna("") for path in paths]
    frame = pd.concat(frames, ignore_index=True)
    missing = REQUIRED_COLUMNS.difference(frame.columns)
    if missing:
        raise typer.BadParameter(f"trait cascade is missing columns: {sorted(missing)}")
    return frame


def _load_tiers(path: Path) -> dict[str, object]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or "analysis_tiers" not in payload:
        raise typer.BadParameter("analysis-tier config must contain analysis_tiers")
    for tier_name, rule in dict(payload["analysis_tiers"]).items():
        allowed = {str(value) for value in rule.get("allowed_fill_tiers", [])}
        unsupported = allowed.difference(GROUNDED_ANALYSIS_TIERS)
        if unsupported:
            raise typer.BadParameter(
                f"analysis tier {tier_name!r} includes non-analysable fills: "
                f"{sorted(unsupported)}"
            )
    return payload


@app.command("audit")
def audit(
    cascade_root: Path = typer.Option(..., exists=True, file_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    """Verify shard completeness and emit trait/tier coverage diagnostics."""
    frame = _load_shards(cascade_root)
    output_dir.mkdir(parents=True, exist_ok=True)

    duplicates = int(frame.duplicated(["accepted_species", "trait_name"]).sum())
    n_species = int(frame["accepted_species"].nunique())
    n_traits = int(frame["trait_name"].nunique())
    expected_rows = n_species * n_traits
    if n_species != 106295:
        raise RuntimeError(f"expected 106295 species, found {n_species}")
    if n_traits < 1:
        raise RuntimeError("trait cascade contains no traits")
    if len(frame) != expected_rows or duplicates:
        raise RuntimeError(
            f"species-trait matrix is not rectangular: rows={len(frame)}, expected={expected_rows}, duplicates={duplicates}"
        )

    coverage = (
        frame.groupby(["trait_name", "fill_tier"], dropna=False)
        .size()
        .rename("n_rows")
        .reset_index()
    )
    totals = coverage.groupby("trait_name")["n_rows"].transform("sum")
    coverage["share_within_trait"] = coverage["n_rows"] / totals
    coverage.to_csv(output_dir / "trait_fill_tier_coverage.csv", index=False)

    summary = {
        "contract": "evidence_grounded_trait_cascade_audit_v2",
        "n_shards": 128,
        "n_rows": int(len(frame)),
        "n_species": n_species,
        "n_traits": n_traits,
        "duplicate_species_trait_rows": duplicates,
        "traits": sorted(frame["trait_name"].unique().tolist()),
        "fill_tier_counts": {
            str(key): int(value) for key, value in frame["fill_tier"].value_counts().items()
        },
        "n_analysis_eligible_rows": int(
            _analysis_eligibility_mask(frame).sum()
        ),
        "interpretation": (
            "Rectangular species-trait matrix. Genus/family rows are qualified inferred "
            "values; unresolved_no_evidence rows carry no biological value and are excluded."
        ),
    }
    (output_dir / "trait_cascade_audit.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(summary, indent=2))


@app.command("join-islands")
def join_islands(
    cascade_root: Path = typer.Option(..., exists=True, file_okay=False),
    island_species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    tiers_config: Path = typer.Option(
        Path("config/trait_cascade_analysis_tiers.yml"), exists=True
    ),
    species_column: str = typer.Option("accepted_species"),
) -> None:
    """Join the completed trait cascade to an island x species occurrence master."""
    frame = _load_shards(cascade_root)
    config = _load_tiers(tiers_config)
    occurrences = pd.read_csv(island_species_csv, dtype=str).fillna("")
    if species_column not in occurrences.columns:
        raise typer.BadParameter(f"island species table missing {species_column}")

    output_dir.mkdir(parents=True, exist_ok=True)
    all_species = set(frame["accepted_species"].astype(str))
    occurrence_species = set(occurrences[species_column].astype(str))
    missing_species = sorted(value for value in occurrence_species.difference(all_species) if value)
    pd.DataFrame({species_column: missing_species}).to_csv(
        output_dir / "island_species_missing_from_trait_cascade.csv", index=False
    )

    base = occurrences.rename(columns={species_column: "accepted_species"})
    # Compatibility with historical artifacts while retaining a hard safety
    # firewall: neither the removed global tier nor unresolved sentinels can
    # enter an analysis even if a config is accidentally broadened.
    analysable = frame.loc[_analysis_eligibility_mask(frame)].copy()
    for tier_name, rule in dict(config["analysis_tiers"]).items():
        allowed = {str(value) for value in rule.get("allowed_fill_tiers", [])}
        subset = analysable.loc[analysable["fill_tier"].isin(allowed)].copy()
        joined = base.merge(subset, on="accepted_species", how="left", validate="many_to_many")
        joined["analysis_tier"] = str(tier_name)
        joined.to_csv(output_dir / f"island_species_traits_{tier_name}.csv.gz", index=False, compression="gzip")

    manifest = {
        "contract": "island_trait_join_v1",
        "n_island_species_rows": int(len(base)),
        "n_unique_island_species": int(base["accepted_species"].nunique()),
        "n_species_missing_from_trait_cascade": int(len(missing_species)),
        "analysis_tiers": list(dict(config["analysis_tiers"]).keys()),
        "source_trait_species": int(frame["accepted_species"].nunique()),
        "source_trait_rows": int(len(frame)),
    }
    (output_dir / "island_trait_join_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
