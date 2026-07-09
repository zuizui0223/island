"""Generate unreviewed proxy candidates for gaps in reported trait evidence.

Reported evidence and proxies are deliberately separate. This module reads one or
more reported-candidate tables, keeps their explicit statements untouched, and
emits only new ``candidate_class=proxy`` rows for species lacking the relevant
reported statement. Proxies are phenotype/family-based candidates, never decided
values and never confirmed pollinator or breeding-system records.
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.callback()
def main() -> None:
    """Generate unreviewed reported-gap proxy candidates."""


PROXY_COLUMNS = [
    "accepted_species",
    "trait_name",
    "proxy_value",
    "candidate_class",
    "basis_traits",
    "basis_family",
    "source_type",
    "raw_description",
    "evidence_scope",
    "review_status",
]

POLLINATION_REPORTED_TRAITS = {"pollination_vector_reported"}
COMPAT_REPORTED_TRAITS = {"self_compatibility_reported", "self_incompatibility_reported"}
REPRO_REPORTED_TRAITS = {"autonomous_selfing_reported", "mating_system_reported"}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _normalise_reported(frame: pd.DataFrame) -> pd.DataFrame:
    table = frame.copy().fillna("")
    if "candidate_class" not in table.columns:
        table["candidate_class"] = "reported"
    if "raw_description" not in table.columns and "evidence_excerpt" in table.columns:
        table["raw_description"] = table["evidence_excerpt"]
    for column in ("accepted_species", "trait_name", "provisional_candidate_value"):
        if column not in table.columns:
            table[column] = ""
    return table


def load_reported(paths: Iterable[Path]) -> pd.DataFrame:
    frames = []
    for path in paths:
        if path.exists():
            frames.append(_normalise_reported(pd.read_csv(path, dtype=str).fillna("")))
    if not frames:
        return pd.DataFrame(columns=["accepted_species", "trait_name", "provisional_candidate_value"])
    return pd.concat(frames, ignore_index=True, sort=False).fillna("")


def _values_by_trait(reported: pd.DataFrame, species: str, trait_name: str) -> set[str]:
    if reported.empty:
        return set()
    subset = reported.loc[
        (reported["accepted_species"].astype(str) == species)
        & (reported["trait_name"].astype(str) == trait_name)
    ]
    return {_text(v).lower() for v in subset["provisional_candidate_value"] if _text(v)}


def _has_reported(reported: pd.DataFrame, species: str, traits: set[str]) -> bool:
    if reported.empty:
        return False
    subset = reported.loc[
        (reported["accepted_species"].astype(str) == species)
        & (reported["trait_name"].astype(str).isin(traits))
    ]
    return len(subset) > 0


def _row(
    species: str,
    trait_name: str,
    value: str,
    basis: list[str],
    family: str,
) -> dict[str, str]:
    return {
        "accepted_species": species,
        "trait_name": trait_name,
        "proxy_value": value,
        "candidate_class": "proxy",
        "basis_traits": "|".join(basis),
        "basis_family": family,
        "source_type": "rule_based_proxy",
        "raw_description": "Proxy generated from reported trait candidates and/or declared family priors; not an explicit source statement.",
        "evidence_scope": "proxy_not_reported",
        "review_status": "unreviewed_proxy",
    }


def floral_proxy_for(
    species: str,
    family: str,
    reported: pd.DataFrame,
    rules: dict[str, Any],
) -> dict[str, str] | None:
    if _has_reported(reported, species, POLLINATION_REPORTED_TRAITS):
        return None
    config = rules.get("floral_syndrome_proxy", {}) or {}
    if family in set(config.get("wind_like_families") or []):
        return _row(species, "floral_syndrome_proxy", "wind_like", ["family_prior:wind_like"], family)

    colours = _values_by_trait(reported, species, "flower_primary_color")
    symmetries = _values_by_trait(reported, species, "floral_symmetry")
    basis = [f"flower_primary_color={v}" for v in sorted(colours)] + [
        f"floral_symmetry={v}" for v in sorted(symmetries)
    ]
    if not basis:
        return None

    large_bee = config.get("large_bee_or_Bombus_like", {}) or {}
    if colours.intersection(set(large_bee.get("colours") or [])) and symmetries.intersection(
        set(large_bee.get("symmetries") or [])
    ):
        return _row(
            species,
            "floral_syndrome_proxy",
            "large_bee_or_Bombus_like_floral_phenotype_proxy",
            basis,
            family,
        )

    butterfly = config.get("butterfly_or_moth_like", {}) or {}
    if colours.intersection(set(butterfly.get("colours") or [])):
        return _row(species, "floral_syndrome_proxy", "butterfly_or_moth_like", basis, family)

    generalist = config.get("open_or_generalist_insect_like", {}) or {}
    if colours.intersection(set(generalist.get("colours") or [])) or symmetries.intersection(
        set(generalist.get("symmetries") or [])
    ):
        return _row(species, "floral_syndrome_proxy", "open_or_generalist_insect_like", basis, family)
    return None


def compatibility_proxy_for(
    species: str,
    family: str,
    reported: pd.DataFrame,
    rules: dict[str, Any],
) -> dict[str, str] | None:
    if _has_reported(reported, species, COMPAT_REPORTED_TRAITS):
        return None
    config = rules.get("compatibility_system_proxy", {}) or {}
    if family in set(config.get("likely_self_compatible_families") or []):
        return _row(species, "compatibility_system_proxy", "likely_self_compatible_proxy", ["family_prior:self_compatible"], family)
    if family in set(config.get("likely_self_incompatible_families") or []):
        return _row(species, "compatibility_system_proxy", "likely_self_incompatible_proxy", ["family_prior:self_incompatible"], family)
    return None


def reproductive_proxy_for(
    species: str,
    family: str,
    reported: pd.DataFrame,
    rules: dict[str, Any],
) -> dict[str, str] | None:
    if _has_reported(reported, species, REPRO_REPORTED_TRAITS):
        return None
    config = rules.get("reproductive_assurance_proxy", {}) or {}
    if family in set(config.get("reproductive_assurance_like_families") or []):
        return _row(
            species,
            "reproductive_assurance_proxy",
            "reproductive_assurance_like_proxy",
            ["family_prior:reproductive_assurance"],
            family,
        )
    return None


def generate_proxies(
    species_df: pd.DataFrame,
    reported: pd.DataFrame,
    rules: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if "accepted_species" not in species_df.columns:
        raise typer.BadParameter("species table must contain accepted_species")
    species_table = species_df.copy().fillna("")
    if "family" not in species_table.columns:
        species_table["family"] = ""

    rows: list[dict[str, str]] = []
    for item in species_table[["accepted_species", "family"]].drop_duplicates().itertuples(index=False):
        species = _text(item.accepted_species)
        family = _text(item.family)
        if not species:
            continue
        for maker in (floral_proxy_for, compatibility_proxy_for, reproductive_proxy_for):
            proxy = maker(species, family, reported, rules)
            if proxy is not None:
                rows.append(proxy)

    frame = pd.DataFrame(rows, columns=PROXY_COLUMNS).drop_duplicates()

    def _by(column: str) -> dict[str, int]:
        if frame.empty:
            return {}
        return {str(k): int(v) for k, v in frame[column].value_counts().sort_index().items()}

    queried = species_table["accepted_species"].astype(str).str.strip().ne("").sum()
    report = {
        "note": "Unreviewed *_proxy candidates only. Reported evidence is not modified, and proxies are not decided values.",
        "n_species_input": int(queried),
        "n_proxy_candidates": int(len(frame)),
        "n_species_with_proxy": int(frame["accepted_species"].nunique()) if len(frame) else 0,
        "n_by_trait_name": _by("trait_name"),
        "n_by_proxy_value": _by("proxy_value"),
    }
    return frame, report


def load_rules(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise typer.BadParameter("proxy rules config must be a mapping")
    return config


@app.command("generate")
def generate(
    species_csv: Path = typer.Option(..., exists=True, help="Species CSV with accepted_species and optional family."),
    reported_csv: list[Path] = typer.Option(..., "--reported-csv", exists=True, help="Reported candidate CSV; may be repeated."),
    output_dir: Path = typer.Option(...),
    rules_config: Path = typer.Option(Path("config/trait_proxy_rules.yml")),
) -> None:
    """Generate unreviewed proxy candidates for gaps in reported evidence."""
    species = pd.read_csv(species_csv, dtype=str).fillna("")
    reported = load_reported(reported_csv)
    frame, report = generate_proxies(species, reported, load_rules(rules_config))
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "trait_proxy_candidates.csv", index=False)
    (output_dir / "trait_proxy_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"{report['n_proxy_candidates']} unreviewed proxy candidate(s) for "
        f"{report['n_species_with_proxy']}/{report['n_species_input']} species. No value decided."
    )


if __name__ == "__main__":
    app()
