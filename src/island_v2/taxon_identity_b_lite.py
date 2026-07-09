"""B-lite taxon identity audit against a POWO/WCVP-style authority table.

This is deliberately an identity-stabilisation layer, not a trait source. It
matches pilot taxa against a downloaded/exported authority CSV and emits accepted
name, family, synonym/status, and review fields so downstream reported/proxy
candidates can stay attached to a stable taxon concept.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.callback()
def main() -> None:
    """Audit taxon identity against a downloaded POWO/WCVP-style authority table."""


OUTPUT_COLUMNS = [
    "submitted_name",
    "authority_matched_name",
    "authority_accepted_name",
    "authority_family",
    "authority_taxon_status",
    "authority_source",
    "authority_source_id",
    "authority_url",
    "match_status",
    "review_status",
]


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _first_existing(frame: pd.DataFrame, names: list[str]) -> str:
    for name in names:
        if name in frame.columns:
            return name
    return ""


def _normalise_authority(authority: pd.DataFrame, source_name: str) -> pd.DataFrame:
    table = authority.copy().fillna("")
    name_col = _first_existing(table, ["taxon_name", "scientificName", "scientific_name", "name"])
    accepted_col = _first_existing(
        table,
        ["accepted_name", "acceptedNameUsage", "accepted_scientific_name", "accepted_species"],
    )
    family_col = _first_existing(table, ["family", "Family"])
    status_col = _first_existing(table, ["taxonomic_status", "taxonomicStatus", "status"])
    id_col = _first_existing(table, ["taxon_id", "taxonID", "id", "kew_id", "powo_id"])
    url_col = _first_existing(table, ["source_url", "url", "references", "powo_url"])

    if not name_col:
        raise typer.BadParameter(
            "authority CSV must contain a taxon name column such as taxon_name, scientificName, or name"
        )
    if not accepted_col:
        accepted_col = name_col

    normalised = pd.DataFrame(
        {
            "authority_matched_name": table[name_col].map(_text),
            "authority_accepted_name": table[accepted_col].map(_text),
            "authority_family": table[family_col].map(_text) if family_col else "",
            "authority_taxon_status": table[status_col].map(_text) if status_col else "",
            "authority_source_id": table[id_col].map(_text) if id_col else "",
            "authority_url": table[url_col].map(_text) if url_col else "",
            "authority_source": source_name,
        }
    )
    normalised["match_key"] = normalised["authority_matched_name"].str.casefold()
    return normalised.loc[normalised["match_key"].ne("")].drop_duplicates("match_key")


def audit_identity(
    species: pd.DataFrame,
    authority: pd.DataFrame,
    source_name: str = "POWO/WCVP_authority_export",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if "accepted_species" not in species.columns:
        raise typer.BadParameter("species CSV must contain accepted_species")
    authority_norm = _normalise_authority(authority, source_name)
    authority_by_key = authority_norm.set_index("match_key", drop=False)

    rows: list[dict[str, str]] = []
    for name in species["accepted_species"].fillna("").astype(str).map(_text).drop_duplicates():
        if not name:
            continue
        key = name.casefold()
        if key in authority_by_key.index:
            hit = authority_by_key.loc[key]
            rows.append(
                {
                    "submitted_name": name,
                    "authority_matched_name": _text(hit["authority_matched_name"]),
                    "authority_accepted_name": _text(hit["authority_accepted_name"]),
                    "authority_family": _text(hit["authority_family"]),
                    "authority_taxon_status": _text(hit["authority_taxon_status"]),
                    "authority_source": _text(hit["authority_source"]),
                    "authority_source_id": _text(hit["authority_source_id"]),
                    "authority_url": _text(hit["authority_url"]),
                    "match_status": "matched_exact_name",
                    "review_status": "unreviewed_identity_match",
                }
            )
        else:
            rows.append(
                {
                    "submitted_name": name,
                    "authority_matched_name": "",
                    "authority_accepted_name": "",
                    "authority_family": "",
                    "authority_taxon_status": "",
                    "authority_source": source_name,
                    "authority_source_id": "",
                    "authority_url": "",
                    "match_status": "unmatched",
                    "review_status": "needs_identity_review",
                }
            )

    frame = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    report = {
        "note": "B-lite identity audit only. No trait value, establishment status, or analysis inclusion is decided.",
        "authority_source": source_name,
        "n_species_input": int(len(frame)),
        "n_matched_exact_name": int((frame["match_status"] == "matched_exact_name").sum()),
        "n_unmatched": int((frame["match_status"] == "unmatched").sum()),
        "n_name_changes_suggested": int(
            (
                frame["authority_accepted_name"].astype(str).str.strip().ne("")
                & (frame["authority_accepted_name"].astype(str) != frame["submitted_name"].astype(str))
            ).sum()
        ),
    }
    return frame, report


@app.command("audit")
def audit(
    species_csv: Path = typer.Option(..., exists=True, help="Species CSV with accepted_species."),
    authority_csv: Path = typer.Option(..., exists=True, help="Downloaded POWO/WCVP-style authority CSV."),
    output_dir: Path = typer.Option(...),
    source_name: str = typer.Option("POWO/WCVP_authority_export"),
) -> None:
    """Audit pilot species identity against a downloaded authority table."""
    species = pd.read_csv(species_csv, dtype=str).fillna("")
    authority = pd.read_csv(authority_csv, dtype=str).fillna("")
    frame, report = audit_identity(species, authority, source_name)
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "taxon_identity_b_lite_audit.csv", index=False)
    (output_dir / "taxon_identity_b_lite_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"Matched {report['n_matched_exact_name']}/{report['n_species_input']} taxon names "
        f"against {source_name}. No trait value decided."
    )


if __name__ == "__main__":
    app()
