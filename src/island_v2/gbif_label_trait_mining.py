"""Mine floral colour from herbarium specimen label text in GBIF occurrences.

The unresolved tail of the master is island endemics with a median of four GBIF
records each. Almost none of them appear in a flora treatment or a pollination
study, so text mining over literature cannot reach them. Every one of those four
records, however, is a pressed specimen, and collectors routinely wrote the
flower colour on the label while the plant was still fresh -- precisely because
colour is the character a dried sheet loses.

That makes label text the one broad source that reaches the tail, and it is
species-direct by construction: the label belongs to a determined specimen.

The matching itself is not this module's. It lives in ``floral_text_matching``
and is shared with the protologue lane, because both lanes feed the same
candidate ledger and the same reviewer, and a rule applied in one and not the
other is not two policies but one policy applied inconsistently. Written on its
own this module had two measured defects the shared rules do not have: it read
"Flowers on reduced peduncles" as red, matching inside "reduced", and it read
"Leaves coriaceous, beneath brown. Flowers seen." as brown, a leaf clause
reaching across a full stop. Both produced a populated cell a reviewer would
read as an observation, which is worse than producing nothing.

What stays here is what is genuinely specific to a specimen label: which
occurrence fields carry label text, which basis-of-record counts as vouchered,
what a collector's shorthand voids, and the fact that a gathering split across
four herbaria is one observation and not four.

Output is a candidate ledger with the verbatim label quote retained. Nothing is
promoted to accepted evidence here.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2 import floral_text_matching as matching

# Outcome names shared with the protologue lane, imported rather than restated
# so a ledger holding rows from both can be read with one vocabulary.
from island_v2.floral_text_matching import (
    ACCEPTED,
    REJECT_COMPETING,
    REJECT_NEGATED,
    REJECT_NO_COLOUR,
    REJECT_NO_ORGAN,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

TRAIT_NAME = "flower_primary_color"
EVIDENCE_SCOPE = "species_direct"
SOURCE_TYPE = "herbarium_specimen_label"

# A record with no label text is a distinct failure from one whose label was
# read and simply says nothing about colour, so it keeps its own name.
REJECT_NO_LABEL = "no_label_text"
REJECT_BASIS = "basis_of_record_not_vouchered"
REJECT_RANK = "not_determined_to_species"
REJECT_OFF_MASTER = "species_not_in_island_master"

__all__ = [
    "ACCEPTED",
    "REJECT_BASIS",
    "REJECT_COMPETING",
    "REJECT_NEGATED",
    "REJECT_NO_COLOUR",
    "REJECT_NO_LABEL",
    "REJECT_NO_ORGAN",
    "REJECT_OFF_MASTER",
    "REJECT_RANK",
    "app",
]


@app.callback()
def main() -> None:
    """Extract species-direct floral colour candidates from specimen labels."""


def load_config(path: Path) -> dict[str, Any]:
    """Load the label-mining configuration and the vocabulary it inherits."""
    required = {
        "label_fields",
        "accepted_basis_of_record",
        "floral_organ_terms",
        "competing_organ_terms",
        "colour_terms",
        "negation_markers",
        "organ_proximity_chars",
        "plain_colour_values",
    }
    return matching.load_config(path, required)


def _text(value: object) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return " ".join(str(value).strip().split())


def join_label_text(row: pd.Series, fields: list[str]) -> str:
    """Concatenate the label-bearing fields present on this occurrence row."""
    parts = [_text(row.get(field)) for field in fields]
    return " ".join(part for part in parts if part)


def extract_colour(
    label_text: str, config: dict[str, Any]
) -> tuple[str, str, str, str]:
    """Return (outcome, normalized_value, matched_terms, verbatim_quote).

    The whole rule is the shared one, so a label and a protologue page reading
    the same sentence produce the same answer. Only the empty-text outcome is
    this lane's own.
    """
    return matching.extract_floral_colour(
        label_text, config, empty_text_outcome=REJECT_NO_LABEL
    )


def binary_plain_class(matched_terms: str, config: dict[str, Any]) -> str:
    """Return "plain", "nonplain" or "unresolved" for the terms that matched.

    This is the level the downstream model actually consumes, and it decides
    cases the five-value ontology cannot: "corolla pale reddish purple" is not
    uniquely red-pink or blue-purple, but it is certainly not plain.
    """
    return matching.binary_plain_class(matched_terms, config)


def mine_occurrences(
    occurrences: pd.DataFrame, master_species: set[str], config: dict[str, Any]
) -> pd.DataFrame:
    """Score every occurrence row and emit the full audit, accepted and not."""
    fields = [str(f) for f in config["label_fields"]]
    basis = {str(b).upper() for b in config["accepted_basis_of_record"]}
    required_rank = str(config.get("required_taxon_rank", "SPECIES")).upper()

    rows: list[dict[str, Any]] = []
    for _, row in occurrences.iterrows():
        species = _text(row.get("species") or row.get("acceptedScientificName"))
        label = join_label_text(row, fields)

        if _text(row.get("basisOfRecord")).upper() not in basis:
            outcome, value, matched, quote = REJECT_BASIS, "", "", ""
        elif _text(row.get("taxonRank")).upper() != required_rank:
            outcome, value, matched, quote = REJECT_RANK, "", "", ""
        elif species not in master_species:
            outcome, value, matched, quote = REJECT_OFF_MASTER, "", "", ""
        else:
            outcome, value, matched, quote = extract_colour(label, config)

        rows.append(
            {
                "accepted_species": species,
                "trait_name": TRAIT_NAME,
                "normalized_value": value,
                "binary_plain_class": binary_plain_class(matched, config),
                "outcome": outcome,
                "matched_colour_terms": matched,
                # The statement the colour was read from, and the whole label it
                # came from. A reviewer needs both: the sentence to judge the
                # claim, the label to judge the specimen.
                "exact_supporting_quote": (quote or label)[:500],
                "label_text": label[:500],
                "gbif_id": _text(row.get("gbifID")),
                "occurrence_id": _text(row.get("occurrenceID")),
                "dataset_key": _text(row.get("datasetKey")),
                "institution_code": _text(row.get("institutionCode")),
                "recorded_by": _text(row.get("recordedBy")),
                "event_date": _text(row.get("eventDate")),
                "locality": _text(row.get("locality")),
                "basis_of_record": _text(row.get("basisOfRecord")),
                "evidence_scope": EVIDENCE_SCOPE,
                "source_type": SOURCE_TYPE,
                "review_status": "unreviewed",
            }
        )
    return pd.DataFrame(rows)


def collapse_to_independent_events(
    audit: pd.DataFrame, config: dict[str, Any]
) -> pd.DataFrame:
    """One accepted candidate per collecting event, not per duplicate sheet.

    A single gathering is commonly split across many herbaria, so counting
    sheets would manufacture independent support that does not exist.
    """
    accepted = audit.loc[audit["outcome"].eq(ACCEPTED)].copy()
    if accepted.empty:
        return accepted.assign(n_duplicate_sheets=pd.Series(dtype=int))

    key_fields = [str(k) for k in (config.get("independence_key") or [])]
    alias = {"recordedBy": "recorded_by", "eventDate": "event_date", "locality": "locality"}
    columns = [
        "accepted_species",
        "trait_name",
        "normalized_value",
        "binary_plain_class",
    ] + [
        alias.get(k, k) for k in key_fields if alias.get(k, k) in accepted.columns
    ]

    accepted["n_duplicate_sheets"] = 1
    grouped = (
        accepted.groupby(columns, dropna=False, as_index=False)
        .agg(
            n_duplicate_sheets=("n_duplicate_sheets", "sum"),
            exact_supporting_quote=("exact_supporting_quote", "first"),
            label_text=("label_text", "first"),
            gbif_id=("gbif_id", "first"),
            dataset_key=("dataset_key", "first"),
            institution_code=("institution_code", "first"),
            matched_colour_terms=("matched_colour_terms", "first"),
        )
    )
    grouped["evidence_scope"] = EVIDENCE_SCOPE
    grouped["source_type"] = SOURCE_TYPE
    grouped["review_status"] = "unreviewed"
    return grouped


@app.command("run")
def run(
    occurrences_csv: Path = typer.Option(..., "--occurrences-csv", exists=True),
    master_taxa_csv: Path = typer.Option(
        Path("data/v2/staging/gbif/collected/island_taxa.csv"), "--master-taxa-csv", exists=True
    ),
    output_dir: Path = typer.Option(..., "--output-dir"),
    config_path: Path = typer.Option(
        Path("config/gbif_label_mining.yml"), "--config", exists=True
    ),
) -> None:
    """Mine label text into a reviewable floral-colour candidate ledger."""
    config = load_config(config_path)
    master = set(pd.read_csv(master_taxa_csv)["accepted_species"].map(_text))
    occurrences = pd.read_csv(occurrences_csv, dtype=str, keep_default_na=False)

    audit = mine_occurrences(occurrences, master, config)
    events = collapse_to_independent_events(audit, config)

    output_dir.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output_dir / "label_mining_audit.csv.gz", index=False)
    events.to_csv(output_dir / "label_colour_candidates.csv.gz", index=False)

    summary = {
        "version": "1.0",
        "n_occurrence_rows": int(len(audit)),
        "outcomes": {str(k): int(v) for k, v in audit["outcome"].value_counts().items()},
        "n_candidate_rows_before_dedup": int((audit["outcome"] == ACCEPTED).sum()),
        "n_independent_collecting_events": int(len(events)),
        "n_species_with_a_candidate": int(events["accepted_species"].nunique())
        if len(events)
        else 0,
        "binary_plain_class": {
            str(k): int(v)
            for k, v in (
                events["binary_plain_class"].value_counts().items() if len(events) else ()
            )
        },
        "interpretation": (
            "Candidates are unreviewed species-direct label statements with the "
            "verbatim quote retained. They are not accepted evidence, and a "
            "species without a candidate is unresolved, not colourless."
        ),
    }
    (output_dir / "label_mining_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    app()
