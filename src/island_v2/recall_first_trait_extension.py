"""Second-stage recall extension for unresolved nine-column traits.

This layer is intentionally inferential. It never overwrites a known value from the
species-level or first recall pass. Its purpose is to increase usable coverage in the
user's likely-enabled analysis layer while keeping every addition auditable and low
confidence.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

from island_v2.recall_first_trait_booster import FIELDS, coverage
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
    "Lamiaceae": {"self_incompatibility": "likely_SC"},
    "Fabaceae": {"self_incompatibility": "likely_SC"},
    "Gesneriaceae": {"self_incompatibility": "likely_SC"},
    "Orchidaceae": {"self_incompatibility": "likely_SC"},
}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _unknown(value: object) -> bool:
    return _text(value) in {"", "unknown"}


def extend_table(
    master: pd.DataFrame,
    table: pd.DataFrame,
    provenance: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = master[["accepted_species", "family"]].drop_duplicates("accepted_species").set_index(
        "accepted_species"
    )
    out = table.copy().fillna("")
    prov_rows = provenance.fillna("").to_dict("records")

    for index, row in out.iterrows():
        species = _text(row["species"])
        family = _text(meta.loc[species, "family"]) if species in meta.index else ""

        if _unknown(row["mating_system"]) and _text(row["self_incompatibility"]) in {
            "SI",
            "likely_SI",
        }:
            out.loc[index, "mating_system"] = "obligate_outcrossing"
            prov_rows.append(
                {
                    "species": species,
                    "field": "mating_system",
                    "value": "obligate_outcrossing",
                    "mode": "likely",
                    "basis": f"self_incompatibility={_text(row['self_incompatibility'])}",
                }
            )

        for field, value in EXTRA_FAMILY_PRIORS.get(family, {}).items():
            if _unknown(out.loc[index, field]):
                out.loc[index, field] = value
                prov_rows.append(
                    {
                        "species": species,
                        "field": field,
                        "value": value,
                        "mode": "likely",
                        "basis": f"extended_family_prior:{family}",
                    }
                )

        has_likely = any(
            record["species"] == species and record["mode"] == "likely" for record in prov_rows
        )
        if has_likely:
            evidence_type = _text(out.loc[index, "evidence_type"])
            if "inference" not in evidence_type.split("|"):
                out.loc[index, "evidence_type"] = (
                    f"{evidence_type}|inference" if evidence_type else "inference"
                )
            if _text(out.loc[index, "confidence"]) not in {"high", "medium"}:
                out.loc[index, "confidence"] = "low"

    validated = validate_result_table(out[OUTPUT_COLUMNS], expected_species=out["species"].tolist())
    prov = pd.DataFrame(prov_rows, columns=["species", "field", "value", "mode", "basis"])
    return validated, prov


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
        "policy": "fill unresolved only; SI/likely_SI may imply low-confidence obligate outcrossing; extended family priors remain inference",
    }
    (output_dir / "recall_first_extended_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(report))


if __name__ == "__main__":
    app()
