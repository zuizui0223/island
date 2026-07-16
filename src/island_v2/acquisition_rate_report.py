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

from island_v2 import angiosperm_scope
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
    required = {
        "master_taxa_csv",
        "master_species_column",
        "core_traits",
        "reported_traits",
        "sources",
    }
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


def load_master_applicability(
    master_path: Path,
    species_column: str,
    config: dict[str, Any],
) -> pd.DataFrame | None:
    """Load the existing taxonomic scope gate without dropping master rows."""
    scope_path_value = _text(config.get("scope_config_path"))
    if not scope_path_value:
        return None
    scope_path = Path(scope_path_value)
    if not scope_path.exists():
        raise typer.BadParameter(f"scope config does not exist: {scope_path}")
    master = pd.read_csv(master_path, dtype=str).fillna("")
    if species_column != "accepted_species":
        master = master.rename(columns={species_column: "accepted_species"})
    scope_config = angiosperm_scope.load_config(scope_path)
    return angiosperm_scope.classify_scope(master, scope_config)


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
        value_col = _first_column(
            frame,
            ["candidate_value", "standardized_value", "trait_value", "value"],
        )
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
        value_col = _first_column(
            frame,
            ["trait_value", "standardized_value", "candidate_value", "value"],
        )
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


def collect_contributions(
    config: dict[str, Any], master: set[str]
) -> tuple[pd.DataFrame, dict[str, int]]:
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


def _species_with_trait(
    contributions: pd.DataFrame, traits: list[str], *, broad_only: bool
) -> set[str]:
    frame = contributions
    if broad_only:
        frame = frame.loc[frame["counts_as_broad"]]
    frame = frame.loc[frame["trait_name"].isin(traits)]
    return set(frame["accepted_species"])


def build_by_trait(
    contributions: pd.DataFrame, config: dict[str, Any], n_master: int
) -> pd.DataFrame:
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
            columns=[
                "track",
                "counts_as_broad",
                "n_species",
                "species_rate",
                "n_species_trait_pairs",
            ]
        )
    rows: list[dict[str, Any]] = []
    for (track, counts_as_broad), group in contributions.groupby(
        ["track", "counts_as_broad"], sort=True
    ):
        species = set(group["accepted_species"])
        rows.append(
            {
                "track": track,
                "counts_as_broad": bool(counts_as_broad),
                "n_species": len(species),
                "species_rate": len(species) / n_master if n_master else 0.0,
                "n_species_trait_pairs": int(
                    len(group.drop_duplicates(["accepted_species", "trait_name"]))
                ),
            }
        )
    return (
        pd.DataFrame(rows)
        .sort_values(["counts_as_broad", "track"], ascending=[False, True])
        .reset_index(drop=True)
    )


def _joined(values: pd.Series) -> str:
    return "|".join(sorted({_text(value) for value in values if _text(value)}))


def build_species_trait_ledger(
    master_species: list[str],
    contributions: pd.DataFrame,
    reported_traits: list[str],
    applicability: pd.DataFrame | None = None,
    angiosperm_only_traits: set[str] | None = None,
) -> pd.DataFrame:
    """Account for every master species x reported-trait cell explicitly."""
    base = pd.MultiIndex.from_product(
        [master_species, reported_traits],
        names=["accepted_species", "trait_name"],
    ).to_frame(index=False)
    if applicability is None:
        base["taxonomic_group"] = ""
        base["angiosperm_analysis_eligible"] = True
    else:
        base = base.merge(
            applicability[
                [
                    "accepted_species",
                    "taxonomic_group",
                    "angiosperm_analysis_eligible",
                ]
            ],
            on="accepted_species",
            how="left",
            validate="many_to_one",
        )
        base["taxonomic_group"] = base["taxonomic_group"].fillna("family_unresolved")
        base["angiosperm_analysis_eligible"] = (
            base["angiosperm_analysis_eligible"].fillna(False).astype(bool)
        )
    gated_traits = angiosperm_only_traits or set()
    base["trait_applicable"] = (
        ~base["trait_name"].isin(gated_traits) | base["angiosperm_analysis_eligible"]
    )
    selected = contributions.loc[contributions["trait_name"].isin(reported_traits)].copy()
    if selected.empty:
        grouped = pd.DataFrame(
            columns=[
                "accepted_species",
                "trait_name",
                "source_keys",
                "tracks",
                "source_backed",
            ]
        )
    else:
        grouped = (
            selected.groupby(["accepted_species", "trait_name"], sort=False)
            .agg(
                source_keys=("source_key", _joined),
                tracks=("track", _joined),
                source_backed=("counts_as_broad", "any"),
            )
            .reset_index()
        )
    ledger = base.merge(
        grouped,
        on=["accepted_species", "trait_name"],
        how="left",
        validate="one_to_one",
    )
    ledger["source_keys"] = ledger["source_keys"].fillna("")
    ledger["tracks"] = ledger["tracks"].fillna("")
    ledger["source_backed"] = ledger["source_backed"].fillna(False).astype(bool)
    ledger.loc[~ledger["trait_applicable"], "source_backed"] = False
    has_candidate = ledger["source_keys"].ne("")
    ledger["acquisition_status"] = "no_candidate_yet"
    ledger.loc[has_candidate & ~ledger["source_backed"], "acquisition_status"] = (
        "candidate_non_source_backed_only"
    )
    ledger.loc[ledger["source_backed"], "acquisition_status"] = "source_backed_candidate"
    ledger["unresolved_reason"] = "not_found_in_configured_sources_yet"
    ledger.loc[has_candidate & ~ledger["source_backed"], "unresolved_reason"] = (
        "only_non_source_backed_candidate"
    )
    ledger.loc[ledger["source_backed"], "unresolved_reason"] = ""
    ledger.loc[~ledger["trait_applicable"], "acquisition_status"] = "not_applicable_taxonomic_scope"
    ledger.loc[~ledger["trait_applicable"], "unresolved_reason"] = (
        "trait_not_applicable_to_non_angiosperm_scope"
    )
    return ledger


def build_species_ledger(
    master_species: list[str],
    species_trait_ledger: pd.DataFrame,
    reported_traits: list[str],
) -> pd.DataFrame:
    """Summarise all trait-cell states without dropping zero-hit species."""
    work = species_trait_ledger.copy()
    work["has_any_candidate"] = work["source_keys"].ne("")
    summary = (
        work.groupby("accepted_species", sort=False)
        .agg(
            taxonomic_group=("taxonomic_group", "first"),
            angiosperm_analysis_eligible=("angiosperm_analysis_eligible", "first"),
            n_reported_traits=("trait_name", "size"),
            n_applicable_traits=("trait_applicable", "sum"),
            n_traits_source_backed=("source_backed", "sum"),
            n_traits_any_channel=("has_any_candidate", "sum"),
        )
        .reindex(master_species)
        .reset_index()
    )
    broad = (
        work.loc[work["source_backed"]]
        .groupby("accepted_species", sort=False)["trait_name"]
        .agg("|".join)
    )
    missing = (
        work.loc[work["trait_applicable"] & ~work["source_backed"]]
        .groupby("accepted_species", sort=False)["trait_name"]
        .agg("|".join)
    )
    summary["source_backed_traits"] = summary["accepted_species"].map(broad).fillna("")
    summary["missing_source_backed_traits"] = summary["accepted_species"].map(missing).fillna("")
    summary["acquisition_status"] = "no_candidate_yet"
    summary["unresolved_reason"] = "not_found_in_configured_sources_yet"
    inference_only = summary["n_traits_any_channel"].gt(0) & summary["n_traits_source_backed"].eq(0)
    partial = summary["n_traits_source_backed"].gt(0) & summary["n_traits_source_backed"].lt(
        summary["n_applicable_traits"]
    )
    not_applicable = summary["n_applicable_traits"].eq(0)
    complete = summary["n_applicable_traits"].gt(0) & summary["n_traits_source_backed"].eq(
        summary["n_applicable_traits"]
    )
    summary.loc[inference_only, "acquisition_status"] = "candidate_non_source_backed_only"
    summary.loc[inference_only, "unresolved_reason"] = "only_non_source_backed_candidates"
    summary.loc[partial, "acquisition_status"] = "partial_source_backed"
    summary.loc[partial, "unresolved_reason"] = "source_backed_candidate_missing_for_some_traits"
    summary.loc[complete, "acquisition_status"] = "all_reported_traits_source_backed"
    summary.loc[complete, "unresolved_reason"] = ""
    summary.loc[not_applicable, "acquisition_status"] = "not_applicable_taxonomic_scope"
    summary.loc[not_applicable, "unresolved_reason"] = (
        "all_reported_traits_outside_angiosperm_scope"
    )
    return summary[
        [
            "accepted_species",
            "taxonomic_group",
            "angiosperm_analysis_eligible",
            "acquisition_status",
            "unresolved_reason",
            "n_reported_traits",
            "n_applicable_traits",
            "n_traits_source_backed",
            "n_traits_any_channel",
            "source_backed_traits",
            "missing_source_backed_traits",
        ]
    ]


def build_summary(
    contributions: pd.DataFrame, config: dict[str, Any], n_master: int
) -> dict[str, Any]:
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
    master_list = load_master_species(
        Path(config["master_taxa_csv"]), config["master_species_column"]
    )
    applicability = load_master_applicability(
        Path(config["master_taxa_csv"]),
        config["master_species_column"],
        config,
    )
    master = set(master_list)
    n_master = len(master_list)

    contributions, off_master = collect_contributions(config, master)
    by_trait = build_by_trait(contributions, config, n_master)
    by_track = build_by_track(contributions, n_master)
    species_trait_ledger = build_species_trait_ledger(
        master_list,
        contributions,
        list(config["reported_traits"]),
        applicability,
        set(config.get("angiosperm_only_traits") or []),
    )
    species_ledger = build_species_ledger(
        master_list,
        species_trait_ledger,
        list(config["reported_traits"]),
    )
    summary = build_summary(contributions, config, n_master)
    summary["off_master_rows_dropped_by_source"] = off_master
    summary["n_master_species_accounted"] = int(len(species_ledger))
    summary["n_species_trait_cells_accounted"] = int(len(species_trait_ledger))
    summary["species_acquisition_status_counts"] = {
        str(key): int(value)
        for key, value in species_ledger["acquisition_status"].value_counts().items()
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    by_trait.to_csv(output_dir / "acquisition_rate_by_trait.csv", index=False)
    by_track.to_csv(output_dir / "acquisition_rate_by_track.csv", index=False)
    species_ledger.to_csv(
        output_dir / "acquisition_status_by_species.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    species_trait_ledger.to_csv(
        output_dir / "acquisition_status_by_species_trait.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
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
