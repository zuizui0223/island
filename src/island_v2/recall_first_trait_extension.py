"""Aggressive final fill pass for the recall-first trait table.

Higher-resolution values always win. Any remaining gaps are filled with explicit
low-confidence global priors so the analysis table is complete; provenance makes
those rows removable in sensitivity analyses.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

from island_v2.recall_first_trait_booster import coverage
from island_v2.v1_category_traits import OUTPUT_COLUMNS, validate_result_table

app = typer.Typer(add_completion=False, no_args_is_help=True)

EXTRA_FAMILY_PRIORS: dict[str, dict[str, str]] = {
    "Apiaceae": {"pollination_guild": "mixed"},
    "Araceae": {"pollination_guild": "flies"},
    "Ranunculaceae": {"pollination_guild": "mixed"},
    "Caryophyllaceae": {"pollination_guild": "mixed"},
    "Amaryllidaceae": {"pollination_guild": "mixed"},
    "Asparagaceae": {"pollination_guild": "mixed"},
    "Arecaceae": {"pollination_guild": "mixed"},
    "Euphorbiaceae": {"pollination_guild": "mixed"},
    "Gentianaceae": {"pollination_guild": "bees"},
    "Melastomataceae": {"pollination_guild": "bees"},
}

# Last-resort defaults for the explicitly inferential recall/sensitivity layer.
# They are not species reports and never overwrite a known value.
GLOBAL_PRIORS = {
    "flower_color": "white",
    "flower_shape": "open / radially symmetric",
    "pollination_guild": "mixed",
    "mating_system": "mixed_mating",
    "self_incompatibility": "likely_SC",
}

PROVENANCE_COLUMNS = [
    "species", "field", "value", "mode", "source", "evidence_scope", "confidence", "basis"
]


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _unknown(value: object) -> bool:
    return _text(value) in {"", "unknown"}


def _record(
    prov_rows: list[dict[str, str]],
    species: str,
    field: str,
    value: str,
    source: str,
    scope: str,
    basis: str,
) -> None:
    prov_rows.append({
        "species": species,
        "field": field,
        "value": value,
        "mode": "likely",
        "source": source,
        "evidence_scope": scope,
        "confidence": "low",
        "basis": basis,
    })


def extend_table(
    master: pd.DataFrame,
    table: pd.DataFrame,
    provenance: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = master.copy().drop_duplicates("accepted_species").set_index("accepted_species")
    out = table.copy().fillna("")
    prov_rows = provenance.fillna("").to_dict("records")

    for index, row in out.iterrows():
        species = _text(row["species"])
        family = _text(meta.loc[species, "family"]) if species in meta.index else ""

        if _unknown(out.loc[index, "mating_system"]) and _text(out.loc[index, "self_incompatibility"]) in {"SI", "likely_SI"}:
            out.loc[index, "mating_system"] = "obligate_outcrossing"
            _record(
                prov_rows, species, "mating_system", "obligate_outcrossing",
                "self_incompatibility", "cross_field_fallback",
                f"self_incompatibility={_text(out.loc[index, 'self_incompatibility'])}",
            )

        for field, value in EXTRA_FAMILY_PRIORS.get(family, {}).items():
            if _unknown(out.loc[index, field]):
                out.loc[index, field] = value
                _record(
                    prov_rows, species, field, value, family, "family_prior",
                    f"extended_family_prior:{family}",
                )

        # Fill every remaining focal gap. This is deliberately aggressive: the
        # provenance scope, rather than leaving unknowns, is the safety mechanism.
        for field, value in GLOBAL_PRIORS.items():
            if _unknown(out.loc[index, field]):
                out.loc[index, field] = value
                _record(
                    prov_rows, species, field, value, "global_prior", "global_prior",
                    "last_resort_fill_for_recall_layer",
                )

        if any(record.get("species") == species and record.get("mode") == "likely" for record in prov_rows):
            evidence_type = _text(out.loc[index, "evidence_type"])
            if "inference" not in evidence_type.split("|"):
                out.loc[index, "evidence_type"] = f"{evidence_type}|inference" if evidence_type else "inference"
            if not any(
                record.get("species") == species and record.get("mode") == "direct"
                for record in prov_rows
            ):
                out.loc[index, "confidence"] = "low"

    validated = validate_result_table(out[OUTPUT_COLUMNS], expected_species=out["species"].tolist())
    prov = pd.DataFrame(prov_rows)
    for column in PROVENANCE_COLUMNS:
        if column not in prov.columns:
            prov[column] = ""
    return validated, prov[PROVENANCE_COLUMNS]


@app.command()
def run(
    master_csv: Path = typer.Option(..., exists=True),
    traits_csv: Path = typer.Option(..., exists=True),
    provenance_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    traits = pd.read_csv(traits_csv, dtype=str).fillna("")
    provenance = pd.read_csv(provenance_csv, dtype=str).fillna("")
    extended, prov = extend_table(master, traits, provenance)
    output_dir.mkdir(parents=True, exist_ok=True)
    extended.to_csv(output_dir / "recall_first_extended_traits.csv", index=False)
    prov.to_csv(output_dir / "recall_first_extended_provenance.csv", index=False)
    report = {
        "n_species": int(len(extended)),
        "coverage": coverage(extended),
        "n_likely_provenance_rows": int(prov["mode"].eq("likely").sum()),
        "scope_counts": prov["evidence_scope"].value_counts().to_dict(),
        "policy": "fill unresolved only; global priors complete the recall layer and remain removable by evidence_scope",
    }
    (output_dir / "recall_first_extended_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(report))


if __name__ == "__main__":
    app()
