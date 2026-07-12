"""Extract the fixed 9-column trait table from frozen shard search packets.

This is deliberately conservative. Species-name-matched snippets are treated as direct
search evidence. Genus/family priors are only used for unresolved fields and are always
marked as low-confidence inference. The raw search packet remains the audit source.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

from island_v2.search_engine_trait_recovery import TRAITS, extract_traits

app = typer.Typer(add_completion=False, no_args_is_help=True)

TABLE_COLUMNS = [
    "species",
    "flower_color",
    "flower_shape",
    "pollination_guild",
    "pollination_notes",
    "mating_system",
    "self_incompatibility",
    "evidence_type",
    "confidence",
]

# Only broad, high-level priors are allowed here. They are sensitivity-layer values,
# never direct species reports.
FAMILY_PRIORS: dict[str, dict[str, str]] = {
    "Poaceae": {
        "flower_shape": "reduced_wind",
        "pollination_guild": "wind",
    },
    "Cyperaceae": {
        "flower_shape": "reduced_wind",
        "pollination_guild": "wind",
    },
    "Juncaceae": {
        "pollination_guild": "wind",
    },
    "Asteraceae": {
        "flower_shape": "composite_head",
        "pollination_guild": "mixed",
    },
    "Lamiaceae": {
        "flower_shape": "zygomorphic",
        "pollination_guild": "bees",
    },
    "Gesneriaceae": {
        "flower_shape": "tubular",
        "pollination_guild": "mixed",
    },
    "Orchidaceae": {
        "flower_shape": "zygomorphic",
        "pollination_guild": "mixed",
    },
}


def clean(value: object) -> str:
    return " ".join(str(value or "").split())


def domain_source_type(url: str, title: str) -> str:
    low = f"{url} {title}".casefold()
    if any(token in low for token in ("flora", "florabase", "efloras", "plants.jstor", "powo.science")):
        return "flora"
    if any(token in low for token in ("garden", "nursery", "hort", "rhs.org")):
        return "horticulture"
    return "review"


def direct_values(packet: dict[str, object]) -> tuple[dict[str, str], list[str], set[str]]:
    species = clean(packet.get("species"))
    species_low = species.casefold()
    values: dict[str, str] = {}
    notes: list[str] = []
    source_types: set[str] = set()

    for search in packet.get("searches", []):
        if not isinstance(search, dict):
            continue
        for result in search.get("results", []):
            if not isinstance(result, dict):
                continue
            title = clean(result.get("title"))
            snippet = clean(result.get("snippet"))
            combined = clean(f"{title}. {snippet}")
            if species_low not in combined.casefold():
                continue
            extracted = extract_traits(combined)
            for trait, value in extracted.items():
                values.setdefault(trait, value)
            if extracted:
                notes.append(combined[:240])
                source_types.add(domain_source_type(clean(result.get("url")), title))
    return values, notes, source_types


def build_row(packet: dict[str, object]) -> dict[str, str]:
    species = clean(packet.get("species"))
    family = clean(packet.get("family"))
    values, notes, source_types = direct_values(packet)
    inferred: list[str] = []

    for trait, value in FAMILY_PRIORS.get(family, {}).items():
        if trait not in values:
            values[trait] = value
            inferred.append(f"{trait} from {family} family-level prior")

    direct_count = sum(trait in values for trait in TRAITS) - len(inferred)
    if direct_count > 0:
        evidence_type = sorted(source_types)[0] if source_types else "review"
        confidence = "medium"
    elif inferred:
        evidence_type = "inference"
        confidence = "low"
    else:
        evidence_type = "inference"
        confidence = "low"

    note_parts: list[str] = []
    if notes:
        note_parts.append(notes[0])
    if inferred:
        note_parts.append("; ".join(inferred))
    if not note_parts:
        note_parts.append("No species-direct trait statement recovered from the two-query search packet.")

    return {
        "species": species,
        "flower_color": values.get("flower_color", "unknown"),
        "flower_shape": values.get("flower_shape", "unknown"),
        "pollination_guild": values.get("pollination_guild", "unknown"),
        "pollination_notes": " | ".join(note_parts),
        "mating_system": values.get("mating_system", "unknown"),
        "self_incompatibility": values.get("self_incompatibility", "unknown"),
        "evidence_type": evidence_type,
        "confidence": confidence,
    }


@app.command()
def extract(
    packets_jsonl: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    packets = [
        json.loads(line)
        for line in packets_jsonl.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    rows = [build_row(packet) for packet in packets]
    table = pd.DataFrame(rows, columns=TABLE_COLUMNS)
    output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_dir / "trait_table.csv", index=False)

    coverage = {
        trait: int(table[trait].ne("unknown").sum())
        for trait in TRAITS
    }
    report = {
        "n_species": len(table),
        "coverage": coverage,
        "n_species_with_any_non_unknown_trait": int(table[list(TRAITS)].ne("unknown").any(axis=1).sum()),
        "n_medium_confidence": int(table["confidence"].eq("medium").sum()),
        "n_low_confidence": int(table["confidence"].eq("low").sum()),
        "policy": "species-name-matched snippet extraction first; unresolved shape/guild may use explicit low-confidence family priors",
    }
    (output_dir / "trait_extraction_status.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
