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
    "basis_source_urls",
    "basis_family",
    "source_type",
    "raw_description",
    "evidence_scope",
    "inference_rule",
    "inference_status",
    "confidence",
    "review_status",
]

POLLINATION_REPORTED_TRAITS = {
    "pollen_vector_mode",
    "pollination_functional_guild",
    "pollination_vector_reported",
}
COMPAT_REPORTED_TRAITS = {
    "self_incompatibility",
    "self_compatibility_reported",
    "self_incompatibility_reported",
}
REPRO_REPORTED_TRAITS = {
    "autonomous_selfing_capacity",
    "autonomous_selfing_reported",
    "cleistogamy",
    "mating_system",
    "mating_system_reported",
}


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
    for source_column in ("standardized_value", "candidate_value", "proxy_value"):
        if source_column in table.columns:
            missing = table["provisional_candidate_value"].astype(str).str.strip().eq("")
            table.loc[missing, "provisional_candidate_value"] = table.loc[missing, source_column]
    if "source_url" not in table.columns:
        table["source_url"] = ""
    return table


def load_reported(paths: Iterable[Path]) -> pd.DataFrame:
    frames = []
    for path in paths:
        if path.exists():
            frames.append(_normalise_reported(pd.read_csv(path, dtype=str).fillna("")))
    if not frames:
        return pd.DataFrame(
            columns=["accepted_species", "trait_name", "provisional_candidate_value"]
        )
    return pd.concat(frames, ignore_index=True, sort=False).fillna("")


def _reported_index(
    reported: pd.DataFrame,
) -> tuple[
    dict[tuple[str, str], set[str]],
    dict[str, set[str]],
    dict[tuple[str, str], set[str]],
]:
    """Index reported rows once so the proxy pass stays linear at global scale."""
    values: dict[tuple[str, str], set[str]] = {}
    traits: dict[str, set[str]] = {}
    urls: dict[tuple[str, str], set[str]] = {}
    for row in _normalise_reported(reported).itertuples(index=False):
        record = row._asdict()
        species = _text(record.get("accepted_species"))
        trait = _text(record.get("trait_name"))
        value = _text(record.get("provisional_candidate_value")).lower()
        if not species or not trait:
            continue
        key = (species, trait)
        traits.setdefault(species, set()).add(trait)
        if value:
            values.setdefault(key, set()).add(value)
        source_url = _text(record.get("source_url"))
        if source_url:
            urls.setdefault(key, set()).add(source_url)
    return values, traits, urls


def _has_reported(traits_by_species: dict[str, set[str]], species: str, traits: set[str]) -> bool:
    return bool(traits_by_species.get(species, set()).intersection(traits))


def _row(
    species: str,
    trait_name: str,
    value: str,
    basis: list[str],
    basis_urls: list[str],
    family: str,
    inference_rule: str,
) -> dict[str, str]:
    return {
        "accepted_species": species,
        "trait_name": trait_name,
        "proxy_value": value,
        "candidate_class": "proxy",
        "basis_traits": "|".join(basis),
        "basis_source_urls": "|".join(basis_urls),
        "basis_family": family,
        "source_type": "rule_based_proxy",
        "raw_description": "Proxy generated from reported trait candidates and/or declared family priors; not an explicit source statement.",
        "evidence_scope": "proxy_not_reported",
        "inference_rule": inference_rule,
        "inference_status": "likely",
        "confidence": "low",
        "review_status": "unreviewed_proxy",
    }


def floral_proxy_for(
    species: str,
    family: str,
    values: dict[tuple[str, str], set[str]],
    traits_by_species: dict[str, set[str]],
    urls: dict[tuple[str, str], set[str]],
    rules: dict[str, Any],
) -> dict[str, str] | None:
    if _has_reported(traits_by_species, species, POLLINATION_REPORTED_TRAITS):
        return None
    config = rules.get("floral_syndrome_proxy", {}) or {}
    if family in set(config.get("wind_like_families") or []):
        return _row(
            species,
            "floral_syndrome_proxy",
            "wind_like",
            ["family_prior:wind_like"],
            [],
            family,
            "likely_wind_like_from_declared_family_prior",
        )

    colours = values.get((species, "flower_primary_color"), set())
    symmetries = values.get((species, "floral_symmetry"), set())
    forms = values.get((species, "floral_form"), set())
    basis = (
        [f"flower_primary_color={v}" for v in sorted(colours)]
        + [f"floral_symmetry={v}" for v in sorted(symmetries)]
        + [f"floral_form={v}" for v in sorted(forms)]
    )
    basis_urls = sorted(
        urls.get((species, "flower_primary_color"), set())
        | urls.get((species, "floral_symmetry"), set())
        | urls.get((species, "floral_form"), set())
    )
    if not colours or not (symmetries or forms):
        return None

    warm_tubular = config.get("bird_or_butterfly_like", {}) or {}
    if colours.intersection(set(warm_tubular.get("colours") or [])) and forms.intersection(
        set(warm_tubular.get("forms") or [])
    ):
        return _row(
            species,
            "floral_syndrome_proxy",
            "bird_or_butterfly_like_floral_phenotype_proxy",
            basis,
            basis_urls,
            family,
            "likely_warm_colour_plus_restrictive_tube",
        )

    large_bee = config.get("large_bee_or_Bombus_like", {}) or {}
    bee_architecture = symmetries.intersection(
        set(large_bee.get("symmetries") or [])
    ) or forms.intersection(set(large_bee.get("forms") or []))
    if colours.intersection(set(large_bee.get("colours") or [])) and bee_architecture:
        return _row(
            species,
            "floral_syndrome_proxy",
            "large_bee_or_Bombus_like_floral_phenotype_proxy",
            basis,
            basis_urls,
            family,
            "likely_bee_like_colour_plus_restrictive_architecture",
        )

    generalist = config.get("open_or_generalist_insect_like", {}) or {}
    generalist_architecture = symmetries.intersection(
        set(generalist.get("symmetries") or [])
    ) or forms.intersection(set(generalist.get("forms") or []))
    if colours.intersection(set(generalist.get("colours") or [])) and generalist_architecture:
        return _row(
            species,
            "floral_syndrome_proxy",
            "open_or_generalist_insect_like",
            basis,
            basis_urls,
            family,
            "likely_open_generalist_from_colour_plus_architecture",
        )
    return None


def compatibility_proxy_for(
    species: str,
    family: str,
    values: dict[tuple[str, str], set[str]],
    traits_by_species: dict[str, set[str]],
    urls: dict[tuple[str, str], set[str]],
    rules: dict[str, Any],
) -> dict[str, str] | None:
    del values, urls
    if _has_reported(traits_by_species, species, COMPAT_REPORTED_TRAITS):
        return None
    config = rules.get("compatibility_system_proxy", {}) or {}
    if family in set(config.get("likely_self_compatible_families") or []):
        return _row(
            species,
            "compatibility_system_proxy",
            "likely_self_compatible_proxy",
            ["family_prior:self_compatible"],
            [],
            family,
            "likely_compatibility_from_declared_family_prior",
        )
    if family in set(config.get("likely_self_incompatible_families") or []):
        return _row(
            species,
            "compatibility_system_proxy",
            "likely_self_incompatible_proxy",
            ["family_prior:self_incompatible"],
            [],
            family,
            "likely_compatibility_from_declared_family_prior",
        )
    return None


def reproductive_proxy_for(
    species: str,
    family: str,
    values: dict[tuple[str, str], set[str]],
    traits_by_species: dict[str, set[str]],
    urls: dict[tuple[str, str], set[str]],
    rules: dict[str, Any],
) -> dict[str, str] | None:
    del values, urls
    if _has_reported(traits_by_species, species, REPRO_REPORTED_TRAITS):
        return None
    config = rules.get("reproductive_assurance_proxy", {}) or {}
    if family in set(config.get("reproductive_assurance_like_families") or []):
        return _row(
            species,
            "reproductive_assurance_proxy",
            "reproductive_assurance_like_proxy",
            ["family_prior:reproductive_assurance"],
            [],
            family,
            "likely_reproductive_assurance_from_declared_family_prior",
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

    values, traits_by_species, source_urls = _reported_index(reported)
    rows: list[dict[str, str]] = []
    for item in (
        species_table[["accepted_species", "family"]].drop_duplicates().itertuples(index=False)
    ):
        species = _text(item.accepted_species)
        family = _text(item.family)
        if not species:
            continue
        for maker in (floral_proxy_for, compatibility_proxy_for, reproductive_proxy_for):
            proxy = maker(
                species,
                family,
                values,
                traits_by_species,
                source_urls,
                rules,
            )
            if proxy is not None:
                rows.append(proxy)

    frame = pd.DataFrame(rows, columns=PROXY_COLUMNS).drop_duplicates()

    def _by(column: str) -> dict[str, int]:
        if frame.empty:
            return {}
        return {str(k): int(v) for k, v in frame[column].value_counts().sort_index().items()}

    queried = species_table["accepted_species"].astype(str).str.strip().ne("").sum()
    report = {
        "note": "Unreviewed likely *_proxy candidates only. Reported evidence is not modified, and proxies are not decided values.",
        "inference_status": "likely",
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
    species_csv: Path = typer.Option(
        ..., exists=True, help="Species CSV with accepted_species and optional family."
    ),
    reported_csv: list[Path] = typer.Option(
        ..., "--reported-csv", exists=True, help="Reported candidate CSV; may be repeated."
    ),
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
