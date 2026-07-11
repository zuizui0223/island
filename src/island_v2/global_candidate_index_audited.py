"""Build the global machine-candidate index with explicit exclusion audits.

Raw campaign waves are immutable acquisition records. When a later audit finds a
scope error, uncertainty statement, or parser artefact, the raw row is retained
in its wave and excluded only from the cumulative review index through a small,
version-controlled rule table. Exclusions never imply biological absence.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2.global_candidate_index import (
    _text,
    _write_gzip,
    build_candidate_index,
    build_counterfactual_screening,
    build_status,
    load_campaign_ledger,
    load_wave_candidates,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

EXCLUSION_KEY = [
    "accepted_species",
    "trait_name",
    "candidate_value",
    "source_url",
]
EXCLUSION_REQUIRED = {
    *EXCLUSION_KEY,
    "exclusion_reason",
    "exclusion_rule",
}


def load_candidate_exclusions(campaign_dir: Path) -> pd.DataFrame:
    """Load optional exact-key candidate exclusions from the campaign directory."""
    path = campaign_dir / "machine_candidate_exclusions.csv"
    columns = [*EXCLUSION_KEY, "exclusion_reason", "exclusion_rule", "added_at"]
    if not path.exists():
        return pd.DataFrame(columns=columns)
    table = pd.read_csv(path, dtype=str).fillna("")
    missing = EXCLUSION_REQUIRED.difference(table.columns)
    if missing:
        raise typer.BadParameter(
            f"candidate exclusion table {path} missing columns: {sorted(missing)}"
        )
    if "added_at" not in table.columns:
        table["added_at"] = ""
    for column in columns:
        table[column] = table[column].map(_text)
    duplicates = table.duplicated(EXCLUSION_KEY, keep=False)
    if duplicates.any():
        duplicated = table.loc[duplicates, EXCLUSION_KEY].drop_duplicates()
        raise typer.BadParameter(
            "candidate exclusion table contains duplicate exact keys: "
            f"{duplicated.to_dict('records')[:5]}"
        )
    return table[columns].copy()


def apply_candidate_exclusions(
    raw: pd.DataFrame,
    exclusions: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return indexable raw rows and the matched exclusion audit."""
    if raw.empty or exclusions.empty:
        empty_columns = [
            *raw.columns,
            "exclusion_reason",
            "exclusion_rule",
            "exclusion_added_at",
        ]
        return raw.copy(), pd.DataFrame(columns=empty_columns)

    frame = raw.copy()
    rules = exclusions.copy()
    for column in EXCLUSION_KEY:
        if column not in frame.columns:
            frame[column] = ""
        frame[column] = frame[column].map(_text)
        rules[column] = rules[column].map(_text)

    rules = rules.rename(columns={"added_at": "exclusion_added_at"})
    merged = frame.merge(
        rules,
        how="left",
        on=EXCLUSION_KEY,
        validate="many_to_one",
        indicator=True,
    )
    matched = merged.loc[merged["_merge"].eq("both")].drop(columns="_merge")
    retained = merged.loc[merged["_merge"].eq("left_only"), frame.columns]
    return retained.reset_index(drop=True), matched.reset_index(drop=True)


def audited_status(
    index: pd.DataFrame,
    screening: pd.DataFrame,
    wave_ids: list[str],
    exclusions: pd.DataFrame,
    matches: pd.DataFrame,
) -> dict[str, Any]:
    status = build_status(index, screening, wave_ids)
    status.update(
        {
            "n_candidate_exclusion_rules": int(len(exclusions)),
            "n_raw_candidate_rows_excluded": int(len(matches)),
            "n_species_with_excluded_candidates": int(
                matches["accepted_species"].nunique()
            )
            if not matches.empty
            else 0,
            "candidate_exclusion_policy": (
                "Raw wave evidence remains immutable. Exact-key audit rules remove "
                "known scope/uncertainty/parser artefacts only from cumulative machine "
                "indexes; exclusions are not biological absences."
            ),
        }
    )
    return status


@app.command("build")
def build(
    campaign_dir: Path = typer.Option(..., exists=True),
    output_dir: Path | None = typer.Option(None),
) -> None:
    """Write exclusion-aware candidate and counterfactual-screening indexes."""
    destination = output_dir or campaign_dir
    destination.mkdir(parents=True, exist_ok=True)

    raw, wave_ids = load_wave_candidates(campaign_dir)
    exclusions = load_candidate_exclusions(campaign_dir)
    retained, matches = apply_candidate_exclusions(raw, exclusions)
    index = build_candidate_index(retained)
    ledger = load_campaign_ledger(campaign_dir)
    screening = build_counterfactual_screening(index, ledger)
    status = audited_status(
        index,
        screening,
        wave_ids,
        exclusions,
        matches,
    )

    _write_gzip(index, destination / "machine_candidate_index.csv.gz")
    _write_gzip(screening, destination / "counterfactual_screening.csv.gz")
    _write_gzip(matches, destination / "machine_candidate_exclusion_matches.csv.gz")
    (destination / "candidate_index_status.json").write_text(
        json.dumps(status, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(
        f"Indexed {status['n_unique_candidate_rows']} machine candidate row(s) from "
        f"{status['n_waves_indexed']} wave(s); excluded "
        f"{status['n_raw_candidate_rows_excluded']} audited raw row(s); "
        f"{len(screening)} biotic screening species."
    )


if __name__ == "__main__":
    app()
