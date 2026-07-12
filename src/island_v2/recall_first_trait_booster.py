"""Fill nine-column traits with a simple accepted-species -> synonym -> genus -> family cascade."""

from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import typer

from island_v2.v1_category_traits import OUTPUT_COLUMNS, validate_result_table

app = typer.Typer(add_completion=False, no_args_is_help=True)

FIELDS = (
    "flower_color",
    "flower_shape",
    "pollination_guild",
    "mating_system",
    "self_incompatibility",
)

FAMILY_PRIORS: dict[str, dict[str, str]] = {
    "Poaceae": {"flower_shape": "reduced / grass-like", "pollination_guild": "wind", "self_incompatibility": "likely_SC"},
    "Cyperaceae": {"flower_shape": "reduced / sedge-like", "pollination_guild": "wind"},
    "Juncaceae": {"flower_shape": "reduced / rush-like", "pollination_guild": "wind"},
    "Asteraceae": {"flower_shape": "composite head / capitulum", "pollination_guild": "mixed"},
    "Lamiaceae": {"flower_shape": "zygomorphic / bilaterally symmetric", "pollination_guild": "bees"},
    "Fabaceae": {"flower_shape": "zygomorphic / bilaterally symmetric", "pollination_guild": "bees"},
    "Orchidaceae": {"flower_shape": "zygomorphic / bilaterally symmetric", "pollination_guild": "mixed"},
    "Gesneriaceae": {"flower_shape": "tubular", "pollination_guild": "mixed"},
    "Ericaceae": {"pollination_guild": "bees"},
    "Campanulaceae": {"pollination_guild": "bees"},
    "Boraginaceae": {"pollination_guild": "bees"},
    "Brassicaceae": {"pollination_guild": "bees", "self_incompatibility": "likely_SI"},
    "Solanaceae": {"pollination_guild": "bees", "self_incompatibility": "likely_SI"},
    "Rosaceae": {"pollination_guild": "bees", "self_incompatibility": "likely_SI"},
    "Papaveraceae": {"pollination_guild": "bees"},
    "Onagraceae": {"pollination_guild": "mixed"},
    "Rubiaceae": {"pollination_guild": "mixed"},
    "Apocynaceae": {"pollination_guild": "mixed"},
    "Malvaceae": {"pollination_guild": "bees"},
    "Myrtaceae": {"pollination_guild": "bees"},
}
SELFING_TO_MATING = {
    "autonomous_selfing": "mainly_selfing",
    "cleistogamy": "mainly_selfing",
    "likely_selfing": "mainly_selfing",
}
ALIAS_COLUMNS = ("synonym", "synonyms", "accepted_name_synonyms", "scientific_name", "original_name")
ALIAS_SPLIT = re.compile(r"\s*[|;]\s*")
PROVENANCE_COLUMNS = [
    "species", "field", "value", "mode", "source", "evidence_scope", "confidence", "basis"
]


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _known(value: object) -> bool:
    return bool(_text(value)) and _text(value).casefold() != "unknown"


def _load(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str).fillna("")


def _genus(species: str) -> str:
    return _text(species).split(" ", 1)[0] if _text(species) else ""


def _aliases(row: pd.Series) -> set[str]:
    values: set[str] = set()
    for column in ALIAS_COLUMNS:
        if column in row.index:
            values.update(value for value in ALIAS_SPLIT.split(_text(row[column])) if value)
    return values


def _mode(values: list[str]) -> str:
    known = [_text(value) for value in values if _known(value)]
    if not known:
        return ""
    counts = Counter(known)
    return sorted(counts, key=lambda value: (-counts[value], value))[0]


def _taxonomy(master: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, str]]:
    required = {"accepted_species", "family"}
    if missing := required.difference(master.columns):
        raise typer.BadParameter(f"master is missing columns: {sorted(missing)}")
    meta = master.copy().drop_duplicates("accepted_species")
    if "genus" not in meta.columns:
        meta["genus"] = meta["accepted_species"].map(_genus)
    else:
        meta["genus"] = meta.apply(
            lambda row: _text(row["genus"]) or _genus(row["accepted_species"]), axis=1
        )
    alias_to_accepted: dict[str, str] = {}
    for _, row in meta.iterrows():
        accepted = _text(row["accepted_species"])
        alias_to_accepted[accepted] = accepted
        for alias in _aliases(row):
            alias_to_accepted.setdefault(alias, accepted)
    return meta, alias_to_accepted


def _taxon_modes(
    meta: pd.DataFrame,
    alias_to_accepted: dict[str, str],
    tables: list[pd.DataFrame],
) -> tuple[dict[tuple[str, str, str], str], dict[tuple[str, str], str]]:
    tax = meta.set_index("accepted_species")[["genus", "family"]].to_dict("index")
    genus_values: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    family_values: dict[tuple[str, str], list[str]] = defaultdict(list)
    for table in tables:
        if "species" not in table.columns:
            continue
        for _, row in table.iterrows():
            reported = _text(row.get("species", ""))
            accepted = alias_to_accepted.get(reported, reported)
            info = tax.get(accepted)
            if not info:
                continue
            family, genus = _text(info["family"]), _text(info["genus"])
            for field in FIELDS:
                value = _text(row.get(field, ""))
                if _known(value):
                    genus_values[(family, genus, field)].append(value)
                    family_values[(family, field)].append(value)
    return (
        {key: _mode(values) for key, values in genus_values.items()},
        {key: _mode(values) for key, values in family_values.items()},
    )


def boost_table(
    master: pd.DataFrame,
    direct: pd.DataFrame,
    priority: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fill unresolved fields by accepted species, synonym, genus, then family."""
    if "species" not in direct.columns or "species" not in priority.columns:
        raise typer.BadParameter("direct and priority tables must contain species")
    meta, alias_to_accepted = _taxonomy(master)
    direct_i = direct.drop_duplicates("species").set_index("species")
    priority_i = priority.drop_duplicates("species").set_index("species")
    genus_values, family_values = _taxon_modes(meta, alias_to_accepted, [direct, priority])
    rows: list[dict[str, str]] = []
    provenance: list[dict[str, str]] = []

    for _, rec in meta.iterrows():
        species, genus, family = (_text(rec["accepted_species"]), _text(rec["genus"]), _text(rec["family"]))
        names = [species, *sorted(_aliases(rec))]
        direct_name = next((name for name in names if name in direct_i.index), "")
        priority_name = next((name for name in names if name in priority_i.index), "")
        d = direct_i.loc[direct_name] if direct_name else pd.Series(dtype=object)
        p = priority_i.loc[priority_name] if priority_name else pd.Series(dtype=object)
        out = {field: "unknown" for field in FIELDS}
        notes: list[str] = []
        used_direct = False
        used_inference = False

        for field in FIELDS:
            candidates = [
                (_text(d.get(field, "")), "species_direct" if direct_name == species else "synonym_direct", "medium", direct_name or species),
                (_text(p.get(field, "")), "species_likely" if priority_name == species else "synonym_likely", "low", priority_name or species),
            ]
            if field == "mating_system":
                selfing = _text(p.get("selfing_status", ""))
                candidates.append((SELFING_TO_MATING.get(selfing, ""), "species_likely", "low", f"selfing_status:{selfing}"))
            candidates.extend([
                (genus_values.get((family, genus, field), ""), "genus_fallback", "low", genus),
                (family_values.get((family, field), ""), "family_fallback", "low", family),
                (FAMILY_PRIORS.get(family, {}).get(field, ""), "family_prior", "low", family),
            ])
            for value, scope, confidence, source in candidates:
                if not _known(value):
                    continue
                out[field] = value
                is_direct = scope in {"species_direct", "synonym_direct"}
                used_direct = used_direct or is_direct
                used_inference = used_inference or not is_direct
                provenance.append({
                    "species": species,
                    "field": field,
                    "value": value,
                    "mode": "direct" if is_direct else "likely",
                    "source": source,
                    "evidence_scope": scope,
                    "confidence": confidence,
                    "basis": f"cascade:{scope}",
                })
                if scope in {"genus_fallback", "family_fallback", "family_prior"}:
                    notes.append(f"{field} filled by {scope}:{source}")
                break

        direct_notes = _text(d.get("pollination_notes", ""))
        if _known(direct_notes):
            notes.insert(0, direct_notes)
        evidence_type = _text(d.get("evidence_type", "")) if used_direct else "inference"
        if used_direct and used_inference and "inference" not in evidence_type.split("|"):
            evidence_type = f"{evidence_type}|inference" if evidence_type else "inference"
        if not evidence_type:
            evidence_type = "inference"
        rows.append({
            "species": species,
            "flower_color": out["flower_color"],
            "flower_shape": out["flower_shape"],
            "pollination_guild": out["pollination_guild"],
            "pollination_notes": " | ".join(notes[:6]) if notes else "unknown",
            "mating_system": out["mating_system"],
            "self_incompatibility": out["self_incompatibility"],
            "evidence_type": evidence_type,
            "confidence": "medium" if used_direct else "low",
        })

    table = validate_result_table(
        pd.DataFrame(rows, columns=OUTPUT_COLUMNS), expected_species=list(meta["accepted_species"])
    )
    return table, pd.DataFrame(provenance, columns=PROVENANCE_COLUMNS)


def coverage(table: pd.DataFrame) -> dict[str, int]:
    return {field: int(table[field].fillna("").ne("unknown").sum()) for field in FIELDS}


@app.command()
def run(
    master_csv: Path = typer.Option(..., exists=True),
    direct_csv: Path = typer.Option(..., exists=True),
    priority_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    direct = _load(direct_csv)
    table, provenance = boost_table(_load(master_csv), direct, _load(priority_csv))
    output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_dir / "recall_first_traits.csv", index=False)
    provenance.to_csv(output_dir / "recall_first_provenance.csv", index=False)
    report = {
        "n_species": int(len(table)),
        "direct_coverage": coverage(direct),
        "recall_first_coverage": coverage(table),
        "n_likely_provenance_rows": int(provenance["mode"].eq("likely").sum()),
        "scope_counts": provenance["evidence_scope"].value_counts().to_dict(),
        "policy": "accepted species -> synonym -> genus -> family; first known value wins; provenance keeps source/evidence_scope/confidence",
    }
    (output_dir / "recall_first_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    typer.echo(json.dumps(report))


if __name__ == "__main__":
    app()
