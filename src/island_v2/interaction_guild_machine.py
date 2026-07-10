"""Build no-human-review pollinator-guild machine layers from interaction claims."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

GUILD_RULES: list[tuple[str, str, tuple[str, ...]]] = [
    ("bumblebees", "partner name contains Bombus", ("bombus",)),
    (
        "other_bees",
        "partner name is a non-Bombus bee taxon",
        (
            "apis",
            "lasioglossum",
            "halictus",
            "halictidae",
            "halictoides",
            "andrena",
            "xylocopa",
            "megachile",
            "colletes",
            "bee",
            "bees",
        ),
    ),
    (
        "flies",
        "partner name is a fly taxon",
        (
            "diptera",
            "calliphoridae",
            "platycheirus",
            "helina",
            "eristalis",
            "bibionidae",
            "musca",
            "eucalliphora",
            "syrphidae",
            "syrphid",
            "fly",
            "flies",
        ),
    ),
    ("butterflies", "partner name is a butterfly taxon", ("butterfly", "papilion", "nymphalid")),
    ("moths", "partner name is a moth taxon", ("moth", "sphingidae", "noctuidae")),
    ("birds", "partner name is a bird taxon", ("bird", "hummingbird", "aves", "nectariniidae")),
    ("bats", "partner name is a bat taxon", ("bat", "chiroptera", "pteropodidae")),
    (
        "beetles_wasps_ants",
        "partner name is a beetle, wasp, ant, or other small insect taxon",
        (
            "cantharidae",
            "coleoptera",
            "beetle",
            "vespidae",
            "wasp",
            "formicidae",
            "ant",
            "dermaptera",
            "earwig",
        ),
    ),
]

OUTPUT_COLUMNS = [
    "machine_guild_candidate_id",
    "accepted_species",
    "pollinator_guild",
    "evidence_type",
    "effectiveness_class",
    "interaction_type",
    "interaction_claim_class",
    "partner_taxon_name",
    "source_evidence_id",
    "source_citation",
    "source_url",
    "study_doi",
    "study_location",
    "event_date",
    "wild_or_cultivated",
    "machine_use_tier",
    "mapping_rule",
    "no_human_review_policy",
]

LONG_FORMAT_COLUMNS = [
    "accepted_species",
    "pollinator_guild",
    "evidence_type",
    "effectiveness_class",
    "source_citation",
    "source_url",
    "study_location",
    "season",
    "wild_or_cultivated",
    "review_status",
]

HOLDOUT_COLUMNS = [
    "source_evidence_id",
    "accepted_species",
    "partner_taxon_name",
    "interaction_type",
    "interaction_claim_class",
    "holdout_reason",
]


@app.callback()
def main() -> None:
    """Create no-human-review pollinator guild machine outputs."""


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _field(record: dict[str, Any], *names: str) -> str:
    for name in names:
        value = _text(record.get(name))
        if value:
            return value
    return ""


def _candidate_id(parts: list[str]) -> str:
    digest = hashlib.sha1("\x1f".join(parts).encode("utf-8")).hexdigest()[:16]
    return f"machine_guild:{digest}"


def classify_partner_guild(partner_taxon_name: str) -> tuple[str, str]:
    """Return (guild, mapping_rule) for a partner taxon name, or blanks."""
    name = _text(partner_taxon_name).casefold()
    if not name:
        return "", ""
    for guild, rule, terms in GUILD_RULES:
        if any(term in name for term in terms):
            return guild, rule
    return "", ""


def _season_from_date(event_date: str) -> str:
    text = _text(event_date)
    if len(text) < 7:
        return ""
    # Keep this coarse and hemisphere-neutral; exact phenology is not needed here.
    month = text[5:7]
    return month if month.isdigit() else ""


def build_machine_guild_candidates(globi: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Map interaction partner names to machine pollinator-guild candidates."""
    rows: list[dict[str, str]] = []
    holdouts: list[dict[str, str]] = []
    for record in globi.fillna("").to_dict("records"):
        species = _field(record, "query_taxon", "accepted_species")
        partner = _text(record.get("partner_taxon_name"))
        interaction_type = _text(record.get("interaction_type"))
        claim_class = _text(record.get("interaction_claim_class"))
        evidence_id = _field(record, "evidence_id", "source_evidence_id")
        if not species:
            continue
        guild, rule = classify_partner_guild(partner)
        if not guild:
            holdouts.append(
                {
                    "source_evidence_id": evidence_id,
                    "accepted_species": species,
                    "partner_taxon_name": partner,
                    "interaction_type": interaction_type,
                    "interaction_claim_class": claim_class,
                    "holdout_reason": "partner taxon name could not be mapped to a configured guild",
                }
            )
            continue
        study_url = _text(record.get("study_url"))
        study_doi = _text(record.get("study_doi"))
        event_date = _text(record.get("event_date"))
        rows.append(
            {
                "machine_guild_candidate_id": _candidate_id(
                    [species, guild, partner, evidence_id, interaction_type]
                ),
                "accepted_species": species,
                "pollinator_guild": guild,
                "evidence_type": "floral_visitor",
                "effectiveness_class": "unknown",
                "interaction_type": interaction_type,
                "interaction_claim_class": claim_class,
                "partner_taxon_name": partner,
                "source_evidence_id": evidence_id,
                "source_citation": _field(record, "study_citation", "study_source_citation"),
                "source_url": study_url or (f"https://doi.org/{study_doi}" if study_doi else ""),
                "study_doi": study_doi,
                "study_location": _text(record.get("locality")),
                "event_date": event_date,
                "wild_or_cultivated": "unknown",
                "machine_use_tier": "source_backed_pollinator_guild_candidate",
                "mapping_rule": rule,
                "no_human_review_policy": (
                    "Machine-only guild candidate from an explicit interaction claim. "
                    "Do not treat as confirmed effectiveness."
                ),
            }
        )
    candidates = pd.DataFrame(rows, columns=OUTPUT_COLUMNS).drop_duplicates(
        "machine_guild_candidate_id"
    )
    holdout = pd.DataFrame(holdouts, columns=HOLDOUT_COLUMNS).drop_duplicates()
    return candidates.reset_index(drop=True), holdout.reset_index(drop=True)


def to_long_pollinator_evidence(candidates: pd.DataFrame) -> pd.DataFrame:
    """Return schema-compatible long-format evidence with review_status=unreviewed."""
    if candidates.empty:
        return pd.DataFrame(columns=LONG_FORMAT_COLUMNS)
    table = pd.DataFrame(
        {
            "accepted_species": candidates["accepted_species"],
            "pollinator_guild": candidates["pollinator_guild"],
            "evidence_type": candidates["evidence_type"],
            "effectiveness_class": candidates["effectiveness_class"],
            "source_citation": candidates["source_citation"],
            "source_url": candidates["source_url"],
            "study_location": candidates["study_location"],
            "season": candidates["event_date"].map(_season_from_date),
            "wild_or_cultivated": candidates["wild_or_cultivated"],
            "review_status": "unreviewed",
        }
    )
    return table[LONG_FORMAT_COLUMNS].drop_duplicates().reset_index(drop=True)


def machine_guild_index(candidates: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "accepted_species",
        "n_machine_pollinator_guilds",
        "machine_pollinator_guilds",
        "n_interaction_rows",
        "has_bumblebee_claim",
        "has_non_bombus_bee_claim",
        "has_fly_claim",
        "machine_functional_replacement_signal",
    ]
    if candidates.empty:
        return pd.DataFrame(columns=columns)
    rows: list[dict[str, Any]] = []
    for species, group in candidates.groupby("accepted_species", sort=True):
        guilds = sorted(set(group["pollinator_guild"]))
        rows.append(
            {
                "accepted_species": species,
                "n_machine_pollinator_guilds": len(guilds),
                "machine_pollinator_guilds": "|".join(guilds),
                "n_interaction_rows": int(len(group)),
                "has_bumblebee_claim": "bumblebees" in guilds,
                "has_non_bombus_bee_claim": "other_bees" in guilds,
                "has_fly_claim": "flies" in guilds,
                "machine_functional_replacement_signal": len(guilds) >= 2,
            }
        )
    return pd.DataFrame(rows, columns=columns)


def summary_dict(
    globi: pd.DataFrame,
    candidates: pd.DataFrame,
    holdouts: pd.DataFrame,
    index: pd.DataFrame,
) -> dict[str, Any]:
    return {
        "n_input_interaction_claims": int(len(globi)),
        "n_machine_guild_candidates": int(len(candidates)),
        "n_holdouts": int(len(holdouts)),
        "n_species_with_machine_guilds": int(candidates["accepted_species"].nunique())
        if not candidates.empty
        else 0,
        "n_species_with_multi_guild_signal": int(
            index["machine_functional_replacement_signal"].sum()
        )
        if not index.empty
        else 0,
        "guild_counts": {
            str(k): int(v)
            for k, v in candidates["pollinator_guild"].value_counts().sort_index().items()
        }
        if not candidates.empty
        else {},
        "policy": (
            "No human review required. These are machine pollinator-guild candidates "
            "from source-backed interaction claims, not effectiveness decisions."
        ),
    }


def write_summary_markdown(
    path: Path,
    summary: dict[str, Any],
    index: pd.DataFrame,
) -> None:
    lines = [
        "# Machine Pollinator Guilds",
        "",
        "Human review was not required for this run. These are machine guild candidates, not confirmed effective pollinators.",
        "",
        "## Counts",
        "",
        f"- Input interaction claims: {summary['n_input_interaction_claims']}",
        f"- Machine guild candidates: {summary['n_machine_guild_candidates']}",
        f"- Holdouts: {summary['n_holdouts']}",
        f"- Species with machine guilds: {summary['n_species_with_machine_guilds']}",
        f"- Species with multi-guild signal: {summary['n_species_with_multi_guild_signal']}",
        "",
        "## Guild Counts",
        "",
    ]
    for guild, count in summary["guild_counts"].items():
        lines.append(f"- `{guild}`: {count}")
    lines.extend(
        [
            "",
            "## Species Index",
            "",
        ]
    )
    if index.empty:
        lines.append("(no mapped guilds)")
    else:
        columns = [
            "accepted_species",
            "n_machine_pollinator_guilds",
            "machine_pollinator_guilds",
            "n_interaction_rows",
            "machine_functional_replacement_signal",
        ]
        lines.append("| " + " | ".join(columns) + " |")
        lines.append("| " + " | ".join("---" for _ in columns) + " |")
        for record in index[columns].astype(str).to_dict("records"):
            values = [record[column].replace("|", "\\|") for column in columns]
            lines.append("| " + " | ".join(values) + " |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@app.command("build")
def build(
    globi_csv: Path = typer.Option(..., exists=True, help="GloBI interaction evidence CSV."),
    output_dir: Path = typer.Option(..., help="Directory for machine guild outputs."),
) -> None:
    """Build machine pollinator-guild candidates and index from GloBI claims."""
    globi = pd.read_csv(globi_csv, dtype=str).fillna("")
    candidates, holdouts = build_machine_guild_candidates(globi)
    long_format = to_long_pollinator_evidence(candidates)
    index = machine_guild_index(candidates)
    summary = summary_dict(globi, candidates, holdouts, index)

    output_dir.mkdir(parents=True, exist_ok=True)
    candidates.to_csv(output_dir / "machine_pollinator_guild_candidates.csv", index=False)
    holdouts.to_csv(output_dir / "machine_pollinator_guild_holdouts.csv", index=False)
    long_format.to_csv(output_dir / "machine_pollinator_evidence_long.csv", index=False)
    index.to_csv(output_dir / "machine_pollinator_guild_index.csv", index=False)
    (output_dir / "machine_pollinator_guild_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    write_summary_markdown(output_dir / "machine_pollinator_guild_summary.md", summary, index)
    typer.echo(
        f"Built {summary['n_machine_guild_candidates']} machine guild candidate(s) "
        f"for {summary['n_species_with_machine_guilds']} species; "
        f"{summary['n_holdouts']} holdout(s)."
    )


if __name__ == "__main__":
    app()
