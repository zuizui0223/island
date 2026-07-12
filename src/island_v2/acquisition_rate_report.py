"""Unified broad trait acquisition-rate report against the master species set.

The v2 design keeps many parallel acquisition channels (curated review, the
machine wave campaign, bulk source joins, and search-enabled LLM extraction).
Each channel writes its own staging output, so no single view answers the
operational question this report exists for: *of the whole master species
universe, what fraction now has a source-backed candidate for each trait, and
which channel supplied it.*

"Broad" here means source-backed coverage, not adjudicated truth. A species-trait
is broad-acquired when any channel whose config marks ``counts_as_broad: true``
supplied a non-empty value for it. Inference/proxy channels (for example the
search-enabled LLM track) are still counted and reported, but separately, and
never fold into the source-backed rate. Presence of a candidate is never a
biological absence for the missing remainder, and a broad-acquired value is
never an accepted trait value.
"""

from __future__ import annotations

import glob as globlib
import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.search_enabled_llm_campaign import (
    OUTPUT_COLUMNS as LLM_OUTPUT_COLUMNS,
    _parse_csv_row as parse_llm_result,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRIBUTION_COLUMNS = ["accepted_species", "trait_name", "source_key", "track", "counts_as_broad"]

# Empty / non-committal values never count as an acquired trait.
EMPTY_VALUES = {"", "unknown", "na", "n/a", "none", "unresolved"}

# search-enabled LLM CSV columns -> v2 canonical trait names.
LLM_TRAIT_MAP = {
    "flower_color": "flower_primary_color",
    "flower_shape": "floral_form",
    "pollination_guild": "pollination_functional_guild",
    "mating_system": "mating_system",
    "self_incompatibility": "self_incompatibility",
}


@app.callback()
def main() -> None:
    """Report the broad, all-channel trait acquisition rate over the master species."""


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _is_present(value: object) -> bool:
    return _text(value).lower() not in EMPTY_VALUES


def load_config(path: Path) -> dict[str, Any]:
    """Load and minimally validate the versioned acquisition-rate configuration."""
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    required = {"master_taxa_csv", "master_species_column", "core_traits", "reported_traits", "sources"}
    if not isinstance(config, dict) or not required.issubset(config):
        raise typer.BadParameter(
            "acquisition-rate config must contain master_taxa_csv, master_species_column, "
            "core_traits, reported_traits, and sources"
        )
    return config


def load_master_species(path: Path, column: str) -> list[str]:
    """Return the ordered, de-duplicated master species denominator."""
    frame = pd.read_csv(path, dtype=str).fillna("")
    if column not in frame.columns:
        raise typer.BadParameter(f"master taxa file missing column {column!r}")
    species = [_text(value) for value in frame[column]]
    seen: dict[str, None] = {}
    for name in species:
        if name and name not in seen:
            seen[name] = None
    return list(seen)


def _normalize_long(
    frame: pd.DataFrame, species_col: str, trait_col: str, value_col: str
) -> pd.DataFrame:
    """Reduce an arbitrary long trait table to present (species, trait) pairs."""
    if frame.empty:
        return pd.DataFrame(columns=["accepted_species", "trait_name"])
    out = pd.DataFrame(
        {
            "accepted_species": frame[species_col].map(_text),
            "trait_name": frame[trait_col].map(_text),
            "_present": frame[value_col].map(_is_present),
        }
    )
    out = out.loc[out["_present"] & out["accepted_species"].ne("") & out["trait_name"].ne("")]
    return out[["accepted_species", "trait_name"]].drop_duplicates()


def _first_column(frame: pd.DataFrame, candidates: list[str]) -> str:
    for name in candidates:
        if name in frame.columns:
            return name
    raise typer.BadParameter(f"table missing any of columns {candidates}")


def adapter_candidate_index(source: dict[str, Any]) -> pd.DataFrame:
    """Read one machine candidate index (deduplicated wave campaign output)."""
    path = Path(source["path"])
    if not path.exists():
        return pd.DataFrame(columns=["accepted_species", "trait_name"])
    frame = pd.read_csv(path, dtype=str).fillna("")
    return _normalize_long(frame, "accepted_species", "trait_name", "candidate_value")


def adapter_candidate_long(source: dict[str, Any]) -> pd.DataFrame:
    """Read generic long trait tables (bulk joins, other channels) by glob."""
    frames: list[pd.DataFrame] = []
    for path in sorted(globlib.glob(source["glob"], recursive=True)):
        frame = pd.read_csv(path, dtype=str).fillna("")
        if frame.empty:
            continue
        species_col = _first_column(frame, ["accepted_species", "scientific_name", "species"])
        value_col = _first_column(frame, ["candidate_value", "trait_value", "value"])
        frames.append(_normalize_long(frame, species_col, "trait_name", value_col))
    if not frames:
        return pd.DataFrame(columns=["accepted_species", "trait_name"])
    return pd.concat(frames, ignore_index=True).drop_duplicates()


def adapter_curated_evidence(source: dict[str, Any]) -> pd.DataFrame:
    """Read human-adjudicated accepted trait evidence by glob."""
    frames: list[pd.DataFrame] = []
    for path in sorted(globlib.glob(source["glob"], recursive=True)):
        frame = pd.read_csv(path, dtype=str).fillna("")
        if frame.empty:
            continue
        if "review_status" in frame.columns:
            frame = frame.loc[frame["review_status"].map(_text).str.lower().eq("accepted")]
        value_col = _first_column(frame, ["trait_value", "candidate_value", "value"])
        frames.append(_normalize_long(frame, "accepted_species", "trait_name", value_col))
    if not frames:
        return pd.DataFrame(columns=["accepted_species", "trait_name"])
    return pd.concat(frames, ignore_index=True).drop_duplicates()


def adapter_search_enabled_llm(source: dict[str, Any]) -> pd.DataFrame:
    """Parse search-enabled LLM result batches into present (species, trait) pairs."""
    rows: list[dict[str, str]] = []
    for path in sorted(globlib.glob(source["glob"], recursive=True)):
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                result = _text(json.loads(line).get("result"))
                parsed = parse_llm_result(result)
            except (ValueError, KeyError, json.JSONDecodeError):
                continue
            species = parsed.get("species", "")
            for llm_column, trait_name in LLM_TRAIT_MAP.items():
                if _is_present(parsed.get(llm_column, "")):
                    rows.append({"accepted_species": species, "trait_name": trait_name})
    if not rows:
        return pd.DataFrame(columns=["accepted_species", "trait_name"])
    frame = pd.DataFrame(rows)
    frame["accepted_species"] = frame["accepted_species"].map(_text)
    return frame.loc[frame["accepted_species"].ne("")].drop_duplicates()


ADAPTERS: dict[str, Callable[[dict[str, Any]], pd.DataFrame]] = {
    "candidate_index": adapter_candidate_index,
    "candidate_long": adapter_candidate_long,
    "curated_evidence": adapter_curated_evidence,
    "search_enabled_llm": adapter_search_enabled_llm,
}

# LLM_OUTPUT_COLUMNS is imported to fail loudly if the upstream contract changes.
assert set(LLM_TRAIT_MAP).issubset(LLM_OUTPUT_COLUMNS)


def collect_contributions(config: dict[str, Any], master: set[str]) -> tuple[pd.DataFrame, dict[str, int]]:
    """Run every configured source adapter and return master-scoped contributions."""
    frames: list[pd.DataFrame] = []
    off_master: dict[str, int] = {}
    for source_key, source in config["sources"].items():
        adapter_name = source.get("adapter")
        if adapter_name not in ADAPTERS:
            raise typer.BadParameter(f"source {source_key!r} has unknown adapter {adapter_name!r}")
        pairs = ADAPTERS[adapter_name](source)
        if pairs.empty:
            continue
        in_master = pairs["accepted_species"].isin(master)
        off_master[source_key] = int((~in_master).sum())
        pairs = pairs.loc[in_master].copy()
        pairs["source_key"] = source_key
        pairs["track"] = source.get("track", source_key)
        pairs["counts_as_broad"] = bool(source.get("counts_as_broad", True))
        frames.append(pairs[CONTRIBUTION_COLUMNS])
    if not frames:
        return pd.DataFrame(columns=CONTRIBUTION_COLUMNS), off_master
    return pd.concat(frames, ignore_index=True).drop_duplicates(), off_master


def _species_with_trait(contributions: pd.DataFrame, traits: list[str], *, broad_only: bool) -> set[str]:
    frame = contributions
    if broad_only:
        frame = frame.loc[frame["counts_as_broad"]]
    frame = frame.loc[frame["trait_name"].isin(traits)]
    return set(frame["accepted_species"])


def build_by_trait(contributions: pd.DataFrame, config: dict[str, Any], n_master: int) -> pd.DataFrame:
    """Per-trait broad vs. any-channel species coverage against the master."""
    rows: list[dict[str, Any]] = []
    for trait in config["reported_traits"]:
        broad = _species_with_trait(contributions, [trait], broad_only=True)
        any_channel = _species_with_trait(contributions, [trait], broad_only=False)
        rows.append(
            {
                "trait_name": trait,
                "n_species_broad": len(broad),
                "broad_rate": len(broad) / n_master if n_master else 0.0,
                "n_species_any_channel": len(any_channel),
                "any_channel_rate": len(any_channel) / n_master if n_master else 0.0,
            }
        )
    return pd.DataFrame(rows)


def build_by_track(contributions: pd.DataFrame, n_master: int) -> pd.DataFrame:
    """Per-channel species reach and species-trait pair counts."""
    if contributions.empty:
        return pd.DataFrame(
            columns=["track", "counts_as_broad", "n_species", "species_rate", "n_species_trait_pairs"]
        )
    rows: list[dict[str, Any]] = []
    for (track, counts_as_broad), group in contributions.groupby(["track", "counts_as_broad"], sort=True):
        species = set(group["accepted_species"])
        rows.append(
            {
                "track": track,
                "counts_as_broad": bool(counts_as_broad),
                "n_species": len(species),
                "species_rate": len(species) / n_master if n_master else 0.0,
                "n_species_trait_pairs": int(len(group.drop_duplicates(["accepted_species", "trait_name"]))),
            }
        )
    return pd.DataFrame(rows).sort_values(["counts_as_broad", "track"], ascending=[False, True]).reset_index(drop=True)


def build_summary(contributions: pd.DataFrame, config: dict[str, Any], n_master: int) -> dict[str, Any]:
    """Headline broad acquisition rates for the core minimum and any trait."""
    colour_traits = list(config["core_traits"].get("colour", []))
    form_traits = list(config["core_traits"].get("form", []))
    colour = _species_with_trait(contributions, colour_traits, broad_only=True)
    form = _species_with_trait(contributions, form_traits, broad_only=True)
    core_minimum = colour & form
    any_broad = set(contributions.loc[contributions["counts_as_broad"], "accepted_species"])
    any_channel = set(contributions["accepted_species"])

    def rate(count: int) -> float:
        return count / n_master if n_master else 0.0

    return {
        "version": "1.0",
        "n_master_species": n_master,
        "n_colour_broad": len(colour),
        "colour_broad_rate": rate(len(colour)),
        "n_form_broad": len(form),
        "form_broad_rate": rate(len(form)),
        "n_core_minimum_broad": len(core_minimum),
        "core_minimum_broad_rate": rate(len(core_minimum)),
        "n_any_trait_broad": len(any_broad),
        "any_trait_broad_rate": rate(len(any_broad)),
        "n_any_trait_any_channel": len(any_channel),
        "any_trait_any_channel_rate": rate(len(any_channel)),
        "interpretation": (
            "Broad, source-backed coverage over the master species universe. Rates are "
            "candidate presence, not accepted trait values or biological absences. "
            "Non-source-backed channels are reported by track but excluded from broad rates."
        ),
    }


@app.command("run")
def run(
    config_path: Path = typer.Option(Path("config/acquisition_rate.yml"), exists=True),
    output_dir: Path = typer.Option(..., help="Directory for the acquisition-rate report outputs."),
) -> None:
    """Write per-trait, per-track, and headline broad acquisition-rate outputs."""
    config = load_config(config_path)
    master_list = load_master_species(Path(config["master_taxa_csv"]), config["master_species_column"])
    master = set(master_list)
    n_master = len(master_list)

    contributions, off_master = collect_contributions(config, master)
    by_trait = build_by_trait(contributions, config, n_master)
    by_track = build_by_track(contributions, n_master)
    summary = build_summary(contributions, config, n_master)
    summary["off_master_rows_dropped_by_source"] = off_master

    output_dir.mkdir(parents=True, exist_ok=True)
    by_trait.to_csv(output_dir / "acquisition_rate_by_trait.csv", index=False)
    by_track.to_csv(output_dir / "acquisition_rate_by_track.csv", index=False)
    (output_dir / "acquisition_rate_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(
        f"Master {n_master} species | colour broad {summary['n_colour_broad']} "
        f"({summary['colour_broad_rate']:.2%}) | core-minimum {summary['n_core_minimum_broad']} "
        f"({summary['core_minimum_broad_rate']:.2%}) | any-trait broad {summary['n_any_trait_broad']} "
        f"({summary['any_trait_broad_rate']:.2%})."
    )


if __name__ == "__main__":
    app()
