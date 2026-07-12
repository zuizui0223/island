"""Build the compact v1-style reproductive/pollination trait table from frozen search evidence.

The output schema intentionally matches the manually prompted v1 workflow while preserving
an auditable distinction between explicit statements and inference. Search snippets are treated
as frozen evidence carriers; missing fields remain unknown.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

OUTPUT_COLUMNS = [
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

SI_TERMS = (
    r"\bself[- ]incompatible\b",
    r"\bself[- ]sterile\b",
    r"\bself incompatibility\b",
)
SC_TERMS = (
    r"\bself[- ]compatible\b",
    r"\bself[- ]fertile\b",
    r"\bself compatibility\b",
)
MATING_PATTERNS = {
    "obligate_outcrossing": (
        r"\bobligate(?:ly)? outcross(?:ing|er)?\b",
        r"\bstrict(?:ly)? outcross(?:ing|er)?\b",
    ),
    "mixed_mating": (
        r"\bmixed mating\b",
        r"\bmixed-mating\b",
        r"\bboth selfing and outcrossing\b",
    ),
    "mainly_selfing": (
        r"\bpredominantly self(?:ing|ed)\b",
        r"\bmainly self(?:ing|ed)\b",
        r"\bprimarily self(?:ing|ed)\b",
    ),
    "obligate_selfing": (
        r"\bobligate(?:ly)? self(?:ing|ed)\b",
        r"\bstrict(?:ly)? self(?:ing|ed)\b",
    ),
}

POLLINATION_ORDER = [
    "bumblebees",
    "bees",
    "flies",
    "butterflies",
    "moths",
    "birds",
    "wind",
    "self",
    "mixed",
]

GUILD_MAP = {
    "bumblebees": "bumblebees",
    "other_bees": "bees",
    "flies": "flies",
    "butterflies": "butterflies",
    "moths": "moths",
    "birds": "birds",
    "wind_water": "wind",
    "abiotic_wind": "wind",
    "self_only": "self",
    "mixed_or_generalist": "mixed",
    "mixed": "mixed",
    "biotic": "mixed",
}


def _text(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return " ".join(str(value).strip().split())


def _species_scoped(text: str, accepted: str, searched: str) -> bool:
    low = _text(text).casefold()
    names = []
    for name in (accepted, searched):
        parts = _text(name).casefold().split()
        if len(parts) >= 2:
            names.append(" ".join(parts[:2]))
    return any(name and name in low for name in names)


def _evidence_type(title: str, url: str, text: str) -> str:
    low = f"{title} {url} {text}".casefold()
    if any(term in low for term in ("pollination study", "field study", "experiment", "doi.org/", "journal")):
        return "field_study"
    if any(term in low for term in ("review", "meta-analysis", "meta analysis")):
        return "review"
    if any(term in low for term in ("flora", "florabase", "efloras", "plantnet", "herbarium")):
        return "flora"
    if any(term in low for term in ("nursery", "garden", "horticulture", "rhs.org", "plant care")):
        return "horticulture"
    return "inference"


def _best_value(rows: pd.DataFrame, trait_name: str) -> tuple[str, str, str]:
    subset = rows.loc[rows["trait_name"].eq(trait_name)].copy()
    if subset.empty:
        return "unknown", "", "low"
    rank_class = {"reported": 0, "likely": 1, "proxy": 1}
    rank_conf = {"high": 0, "medium": 1, "low": 2}
    subset["_rank_class"] = subset["candidate_class"].map(rank_class).fillna(2)
    subset["_rank_conf"] = subset["confidence"].map(rank_conf).fillna(3)
    subset = subset.sort_values(["_rank_class", "_rank_conf"], kind="stable")
    row = subset.iloc[0]
    return (
        _text(row.get("provisional_candidate_value")) or "unknown",
        _text(row.get("evidence_quote")),
        _text(row.get("confidence")) or "low",
    )


def _extract_si_mating(search_results: pd.DataFrame) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for record in search_results.fillna("").to_dict("records"):
        species = _text(record.get("accepted_species"))
        searched = _text(record.get("searched_name"))
        title = _text(record.get("result_title"))
        snippet = _text(record.get("result_snippet"))
        combined = f"{title}. {snippet}".strip()
        if not species or not _species_scoped(combined, species, searched):
            continue
        low = combined.casefold()
        state = result.setdefault(
            species,
            {
                "si": [],
                "mating": [],
                "notes": [],
                "evidence_types": [],
            },
        )
        for pattern in SI_TERMS:
            if re.search(pattern, low):
                state["si"].append(("SI", "high", combined))
                break
        for pattern in SC_TERMS:
            if re.search(pattern, low):
                state["si"].append(("SC", "high", combined))
                break
        for value, patterns in MATING_PATTERNS.items():
            if any(re.search(pattern, low) for pattern in patterns):
                state["mating"].append((value, "high", combined))
                break
        if "pollinat" in low or "self" in low or "outcross" in low:
            state["notes"].append(combined)
            state["evidence_types"].append(
                _evidence_type(title, _text(record.get("result_url")), combined)
            )
    return result


def _choose_direct(values: list[tuple[str, str, str]]) -> tuple[str, str, str]:
    if not values:
        return "unknown", "", "low"
    distinct = {value for value, _, _ in values}
    if len(distinct) > 1:
        notes = " | ".join(dict.fromkeys(note for _, _, note in values))
        return "unknown", notes, "low"
    return values[0]


def build_v1_style_table(
    species_table: pd.DataFrame,
    candidates: pd.DataFrame,
    search_results: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    direct = _extract_si_mating(search_results)
    candidate_groups = {
        species: frame.copy()
        for species, frame in candidates.fillna("").groupby("accepted_species", sort=False)
    }

    rows: list[dict[str, str]] = []
    for record in species_table.fillna("").to_dict("records"):
        species = _text(record.get("accepted_species"))
        if not species:
            continue
        cands = candidate_groups.get(species, pd.DataFrame(columns=candidates.columns))
        color, color_note, color_conf = _best_value(cands, "flower_primary_color")
        shape, shape_note, shape_conf = _best_value(cands, "floral_form")
        if shape == "unknown":
            shape, shape_note, shape_conf = _best_value(cands, "floral_symmetry")

        guild_candidates: list[tuple[str, str, str]] = []
        for trait in ("pollination_functional_guild", "pollen_vector_mode"):
            subset = cands.loc[cands.get("trait_name", pd.Series(dtype=str)).eq(trait)]
            for item in subset.to_dict("records"):
                mapped = GUILD_MAP.get(_text(item.get("provisional_candidate_value")), "")
                if mapped:
                    guild_candidates.append(
                        (
                            mapped,
                            _text(item.get("confidence")) or "low",
                            _text(item.get("evidence_quote")),
                        )
                    )
        guild = "unknown"
        guild_note = ""
        guild_conf = "low"
        if guild_candidates:
            by_order = {name: i for i, name in enumerate(POLLINATION_ORDER)}
            guild_candidates.sort(key=lambda x: by_order.get(x[0], 999))
            guild, guild_conf, guild_note = guild_candidates[0]

        state = direct.get(species, {"si": [], "mating": [], "notes": [], "evidence_types": []})
        si, si_note, si_conf = _choose_direct(state["si"])
        mating, mating_note, mating_conf = _choose_direct(state["mating"])

        notes = [note for note in (guild_note, si_note, mating_note, *state["notes"]) if note]
        notes = list(dict.fromkeys(notes))[:3]
        pollination_notes = " | ".join(notes) if notes else "unknown"

        evidence_types = [value for value in state["evidence_types"] if value]
        if not evidence_types and not cands.empty:
            evidence_types = ["inference"]
        evidence_type = evidence_types[0] if evidence_types else "inference"

        confidence_values = [color_conf, shape_conf, guild_conf, si_conf, mating_conf]
        non_low = [value for value in confidence_values if value in {"high", "medium"}]
        if "high" in non_low:
            confidence = "high"
        elif "medium" in non_low:
            confidence = "medium"
        else:
            confidence = "low"

        rows.append(
            {
                "species": species,
                "flower_color": color,
                "flower_shape": shape,
                "pollination_guild": guild,
                "pollination_notes": pollination_notes,
                "mating_system": mating,
                "self_incompatibility": si,
                "evidence_type": evidence_type,
                "confidence": confidence,
            }
        )

    frame = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    tracked = [
        "flower_color",
        "flower_shape",
        "pollination_guild",
        "mating_system",
        "self_incompatibility",
    ]
    coverage = {
        column: int(frame[column].ne("unknown").sum()) for column in tracked
    }
    report = {
        "n_species": int(len(frame)),
        "coverage": coverage,
        "n_species_with_any_non_unknown_trait": int(
            frame[tracked].ne("unknown").any(axis=1).sum()
        ),
        "policy": (
            "Compact v1-style output. Explicit search-result statements and likely/proxy "
            "candidates are retained upstream; unknown is used when no supported value exists."
        ),
    }
    return frame, report


@app.command()
def build(
    species_csv: Path = typer.Option(..., exists=True),
    candidates_csv: Path = typer.Option(..., exists=True),
    search_results_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    species = pd.read_csv(species_csv, dtype=str).fillna("")
    candidates = pd.read_csv(candidates_csv, dtype=str).fillna("")
    search = pd.read_csv(search_results_csv, dtype=str).fillna("")
    frame, report = build_v1_style_table(species, candidates, search)
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "v1_style_trait_table.csv", index=False)
    (output_dir / "v1_style_trait_table_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
