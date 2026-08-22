"""Mine floral colour from herbarium specimen label text in GBIF occurrences.

The unresolved tail of the master is island endemics with a median of four GBIF
records each. Almost none of them appear in a flora treatment or a pollination
study, so text mining over literature cannot reach them. Every one of those four
records, however, is a pressed specimen, and collectors routinely wrote the
flower colour on the label while the plant was still fresh -- precisely because
colour is the character a dried sheet loses.

That makes label text the one broad source that reaches the tail, and it is
species-direct by construction: the label belongs to a determined specimen.

The extraction is deliberately literal. A label is a handful of words written in
the field, not prose to be interpreted, so this module anchors every colour word
to a floral organ within a short window, rejects the row outright when the
nearest organ is a fruit or a leaf, and voids statements a collector negated.
That discipline is not theoretical: in the WFO rejection ledger 38.6% of
rejected candidates were a style, calyx or petal measurement being read as a
whole-flower value.

Output is a candidate ledger with the verbatim label quote retained. Nothing is
promoted to accepted evidence here.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

TRAIT_NAME = "flower_primary_color"
EVIDENCE_SCOPE = "species_direct"
SOURCE_TYPE = "herbarium_specimen_label"

REJECT_NO_LABEL = "no_label_text"
REJECT_BASIS = "basis_of_record_not_vouchered"
REJECT_RANK = "not_determined_to_species"
REJECT_OFF_MASTER = "species_not_in_island_master"
REJECT_NEGATED = "statement_negated_or_uncertain"
REJECT_NO_COLOUR = "no_colour_term"
REJECT_NO_ORGAN = "colour_not_anchored_to_floral_organ"
REJECT_COMPETING = "colour_belongs_to_non_floral_organ"
ACCEPTED = "candidate"


@app.callback()
def main() -> None:
    """Extract species-direct floral colour candidates from specimen labels."""


def load_config(path: Path) -> dict[str, Any]:
    """Load and minimally validate the versioned label-mining configuration."""
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    required = {
        "label_fields",
        "accepted_basis_of_record",
        "floral_organ_terms",
        "competing_organ_terms",
        "colour_terms",
        "negation_markers",
        "organ_proximity_chars",
    }
    if not isinstance(config, dict) or not required.issubset(config):
        raise typer.BadParameter(f"config must contain {sorted(required)}")
    return config


def _text(value: object) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return " ".join(str(value).strip().split())


def join_label_text(row: pd.Series, fields: list[str]) -> str:
    """Concatenate the label-bearing fields present on this occurrence row."""
    parts = [_text(row.get(field)) for field in fields]
    return " ".join(part for part in parts if part)


def _nearest_organ(
    text: str, position: int, floral: list[str], competing: list[str], window: int
) -> str:
    """Return "floral", "competing" or "" for the organ closest to a colour word.

    Distance is measured to the nearest occurrence of any organ term on either
    side. Ties go to the competing organ, so "fruits and flowers red" is
    rejected rather than accepted -- the conservative reading.
    """
    best_kind = ""
    best_distance = window + 1
    for kind, terms in (("competing", competing), ("floral", floral)):
        for term in terms:
            start = 0
            while True:
                index = text.find(term, start)
                if index == -1:
                    break
                distance = (
                    0
                    if index <= position <= index + len(term)
                    else min(abs(position - index), abs(position - (index + len(term))))
                )
                if distance <= window and (
                    distance < best_distance
                    or (distance == best_distance and kind == "competing")
                ):
                    best_distance = distance
                    best_kind = kind
                start = index + 1
    return best_kind


def extract_colour(label_text: str, config: dict[str, Any]) -> tuple[str, str, str]:
    """Return (outcome, normalized_value, matched_terms) for one label."""
    text = label_text.lower()
    if not text:
        return REJECT_NO_LABEL, "", ""

    for marker in config["negation_markers"]:
        if str(marker).lower() in text:
            return REJECT_NEGATED, "", ""

    colours: dict[str, str] = {str(k).lower(): str(v) for k, v in config["colour_terms"].items()}
    floral = [str(t).lower() for t in config["floral_organ_terms"]]
    competing = [str(t).lower() for t in config["competing_organ_terms"]]
    window = int(config["organ_proximity_chars"])

    values: list[str] = []
    matched: list[str] = []
    saw_colour = False
    saw_competing = False

    # Longer terms first so "yellowish" is not consumed by "yellow".
    for term in sorted(colours, key=len, reverse=True):
        start = 0
        while True:
            index = text.find(term, start)
            if index == -1:
                break
            start = index + len(term)
            # skip a hit that is merely the tail of a longer word already matched
            if any(term in seen for seen in matched):
                continue
            saw_colour = True
            kind = _nearest_organ(text, index, floral, competing, window)
            if kind == "competing":
                saw_competing = True
                continue
            if kind != "floral":
                continue
            value = colours[term]
            matched.append(term)
            if value not in values:
                values.append(value)

    if not saw_colour:
        return REJECT_NO_COLOUR, "", ""
    if not values:
        return (REJECT_COMPETING if saw_competing else REJECT_NO_ORGAN), "", ""
    if len(values) > 1:
        return ACCEPTED, str(config.get("multiple_value_result", "multicolored_variable")), "+".join(matched)
    return ACCEPTED, values[0], "+".join(matched)


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
            outcome, value, matched = REJECT_BASIS, "", ""
        elif _text(row.get("taxonRank")).upper() != required_rank:
            outcome, value, matched = REJECT_RANK, "", ""
        elif species not in master_species:
            outcome, value, matched = REJECT_OFF_MASTER, "", ""
        else:
            outcome, value, matched = extract_colour(label, config)

        rows.append(
            {
                "accepted_species": species,
                "trait_name": TRAIT_NAME,
                "normalized_value": value,
                "outcome": outcome,
                "matched_colour_terms": matched,
                "exact_supporting_quote": label[:500],
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
    columns = ["accepted_species", "trait_name", "normalized_value"] + [
        alias.get(k, k) for k in key_fields if alias.get(k, k) in accepted.columns
    ]

    accepted["n_duplicate_sheets"] = 1
    grouped = (
        accepted.groupby(columns, dropna=False, as_index=False)
        .agg(
            n_duplicate_sheets=("n_duplicate_sheets", "sum"),
            exact_supporting_quote=("exact_supporting_quote", "first"),
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
