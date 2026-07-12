"""Boost nine-column trait recall without overwriting direct evidence.

The public-web campaign remains the source of direct species-level values. This module
adds an explicitly inferential fallback layer for unresolved fields, using the existing
priority/likely outputs first and a small set of auditable family priors last.

The resulting table is intended for the user's recall-first analysis layer. Every value
added here is marked through ``evidence_type=inference`` and low confidence unless the
row also contains direct evidence. Direct values always take precedence.
"""

from __future__ import annotations

import json
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

# These are deliberately coarse sensitivity-layer priors, never species reports.
# The goal is to recover likely ecological classes for large-scale analysis while
# keeping provenance explicit and confidence low.
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


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _known(value: object) -> bool:
    return bool(_text(value)) and _text(value) != "unknown"


def _load(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str).fillna("")


def boost_table(
    master: pd.DataFrame,
    direct: pd.DataFrame,
    priority: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return boosted nine-column rows and field-level provenance."""
    required_master = {"accepted_species", "family"}
    if missing := required_master.difference(master.columns):
        raise typer.BadParameter(f"master is missing columns: {sorted(missing)}")
    if "species" not in direct.columns or "species" not in priority.columns:
        raise typer.BadParameter("direct and priority tables must contain species")

    meta = master[["accepted_species", "family"]].drop_duplicates("accepted_species")
    direct_i = direct.drop_duplicates("species").set_index("species")
    priority_i = priority.drop_duplicates("species").set_index("species")

    rows: list[dict[str, str]] = []
    provenance: list[dict[str, str]] = []

    for rec in meta.itertuples(index=False):
        species = _text(rec.accepted_species)
        family = _text(rec.family)
        d = direct_i.loc[species] if species in direct_i.index else pd.Series(dtype=object)
        p = priority_i.loc[species] if species in priority_i.index else pd.Series(dtype=object)
        out = {field: "unknown" for field in FIELDS}
        notes: list[str] = []
        used_direct = False
        used_inference = False

        for field in FIELDS:
            direct_value = _text(d.get(field, ""))
            if _known(direct_value):
                out[field] = direct_value
                used_direct = True
                provenance.append({"species": species, "field": field, "value": direct_value, "mode": "direct", "basis": "species-level public-web evidence"})
                continue

            priority_value = _text(p.get(field, ""))
            if _known(priority_value):
                out[field] = priority_value
                used_inference = True
                provenance.append({"species": species, "field": field, "value": priority_value, "mode": "likely", "basis": "existing deterministic likely/priority layer"})
                continue

            if field == "mating_system":
                selfing = _text(p.get("selfing_status", ""))
                mapped = SELFING_TO_MATING.get(selfing)
                if mapped:
                    out[field] = mapped
                    used_inference = True
                    provenance.append({"species": species, "field": field, "value": mapped, "mode": "likely", "basis": f"selfing_status={selfing}"})
                    continue

            prior = FAMILY_PRIORS.get(family, {}).get(field)
            if prior:
                out[field] = prior
                used_inference = True
                notes.append(f"{field} inferred from {family} family prior")
                provenance.append({"species": species, "field": field, "value": prior, "mode": "likely", "basis": f"family_prior:{family}"})

        direct_notes = _text(d.get("pollination_notes", ""))
        if _known(direct_notes):
            notes.insert(0, direct_notes)

        evidence_type = _text(d.get("evidence_type", "")) if used_direct else "inference"
        if used_direct and used_inference and "inference" not in evidence_type.split("|"):
            evidence_type = f"{evidence_type}|inference" if evidence_type else "inference"
        if not evidence_type:
            evidence_type = "inference"

        confidence = "medium" if used_direct else "low"
        rows.append(
            {
                "species": species,
                "flower_color": out["flower_color"],
                "flower_shape": out["flower_shape"],
                "pollination_guild": out["pollination_guild"],
                "pollination_notes": " | ".join(notes[:4]) if notes else "unknown",
                "mating_system": out["mating_system"],
                "self_incompatibility": out["self_incompatibility"],
                "evidence_type": evidence_type,
                "confidence": confidence,
            }
        )

    table = validate_result_table(pd.DataFrame(rows, columns=OUTPUT_COLUMNS), expected_species=list(meta["accepted_species"]))
    prov = pd.DataFrame(provenance, columns=["species", "field", "value", "mode", "basis"])
    return table, prov


def coverage(table: pd.DataFrame) -> dict[str, int]:
    return {field: int(table[field].ne("unknown").sum()) for field in FIELDS}


@app.command()
def run(
    master_csv: Path = typer.Option(..., exists=True),
    direct_csv: Path = typer.Option(..., exists=True),
    priority_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    master = _load(master_csv)
    direct = _load(direct_csv)
    priority = _load(priority_csv)
    table, provenance = boost_table(master, direct, priority)
    output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_dir / "recall_first_traits.csv", index=False)
    provenance.to_csv(output_dir / "recall_first_provenance.csv", index=False)
    report = {
        "n_species": int(len(table)),
        "direct_coverage": coverage(direct),
        "recall_first_coverage": coverage(table),
        "n_likely_provenance_rows": int(provenance["mode"].eq("likely").sum()),
        "policy": "direct species evidence wins; existing likely layer second; auditable family priors last; inferred values are low-confidence",
    }
    (output_dir / "recall_first_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    typer.echo(json.dumps(report))


if __name__ == "__main__":
    app()
