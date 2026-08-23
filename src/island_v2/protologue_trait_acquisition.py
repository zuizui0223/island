"""Acquire floral descriptions from protologues, via IPNI and BHL.

The unresolved tail of the master is 74,846 species occurring on three islands
or fewer, with a median of four GBIF records each. Acquisition so far has been
source-driven toward large continental floras, which is why the tail was never
touched: no monograph, revision or pollination study covers most of these
species. One description does exist for every one of them, because the name was
validly published -- the protologue. For a species with four records it is
frequently the only floral description anywhere.

Two things separate this module from the other text-mining paths.

Multilingual by requirement, not preference. Protologues are Latin, French,
German or Spanish far more often than English, and this is exactly where an
English-only rule set has already been measured to fail: in the WFO rejection
ledger the one automatically recoverable slice was non-English statements the
rules never bound. Matching therefore runs over folded text -- lowercased,
diacritics stripped, ``ss`` for eszett -- against enumerated forms in all four
languages.

Word boundaries, not substrings. Label text is a handful of words; protologue
prose runs for paragraphs, and substring matching over it would read
``rotundifolia`` as German ``rot`` and ``albacea`` as Latin ``alba``. Those hits
would outnumber the real ones, so every term here is matched at a boundary.

Network access is confined to two thin clients, each taking an injected
``fetch`` so the extraction and selection logic is exercised without a network.
Items still in copyright are skipped, never scraped, and a record missing any
required lineage field is rejected rather than downgraded.

Output is a candidate ledger of unreviewed statements carrying the verbatim
quote. Nothing is promoted to accepted evidence here.
"""

from __future__ import annotations

import json
import re
import unicodedata
from datetime import date
from pathlib import Path
from typing import Any, Callable

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

TRAIT_NAME = "flower_primary_color"
EVIDENCE_SCOPE = "species_direct"
SOURCE_TYPE = "protologue"

REJECT_NO_CITATION = "no_protologue_citation"
REJECT_CITATION_INCOMPLETE = "citation_missing_required_fields"
REJECT_NO_SCAN = "no_public_domain_scan"
REJECT_IN_COPYRIGHT = "scan_in_copyright"
REJECT_NO_TEXT = "no_ocr_text"
REJECT_NEGATED = "statement_negated_or_uncertain"
REJECT_NO_COLOUR = "no_colour_term"
REJECT_NO_ORGAN = "colour_not_anchored_to_floral_organ"
REJECT_COMPETING = "colour_belongs_to_non_floral_organ"
REJECT_MISSING_LINEAGE = "incomplete_source_lineage"
ACCEPTED = "candidate"

Fetch = Callable[[str, dict[str, str]], Any]


@app.callback()
def main() -> None:
    """Acquire species-direct floral colour candidates from protologues."""


def load_config(path: Path) -> dict[str, Any]:
    """Load and minimally validate the versioned protologue configuration."""
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    required = {
        "target_selection",
        "citation_index",
        "scan_provider",
        "floral_organ_terms",
        "competing_organ_terms",
        "colour_terms",
        "latin_adjective_stems",
        "latin_adjective_endings",
        "negation_markers",
        "organ_proximity_chars",
        "required_lineage_fields",
    }
    if not isinstance(config, dict) or not required.issubset(config):
        raise typer.BadParameter(f"config must contain {sorted(required)}")
    return config


def expand_colour_terms(config: dict[str, Any]) -> dict[str, str]:
    """The full folded colour vocabulary, with Latin declensions applied.

    Latin adjectives agree with the organ noun, so which form appears is fixed
    by the sentence -- "corolla alba" but "flores albi" but "floribus albis" --
    and protologues use the plural at least as often as the dictionary form.
    Enumerating only the singular would drop most real statements, so stems are
    declined here instead. A generated non-word never matches, which makes
    over-generation free and under-generation the only real cost.
    """
    terms: dict[str, str] = {}
    for stems_key, endings_key in (
        ("latin_adjective_stems", "latin_adjective_endings"),
        ("latin_third_declension_stems", "latin_third_declension_endings"),
    ):
        endings = [str(e) for e in (config.get(endings_key) or [])]
        for stem, value in (config.get(stems_key) or {}).items():
            for ending in endings:
                terms[fold(f"{stem}{ending}")[0]] = str(value)
    # Literal entries win over generated ones: an irregular listed by hand is
    # there precisely because the stem rule cannot produce it correctly.
    for term, value in config["colour_terms"].items():
        terms[fold(str(term))[0]] = str(value)
    return terms


def _text(value: object) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return " ".join(str(value).strip().split())


def fold(text: str) -> tuple[str, list[int]]:
    """Fold text for matching and keep a map back to the original offsets.

    Folding lowercases, strips combining marks so ``blanchâtre`` matches
    ``blanchatre``, and expands eszett so ``weiß`` matches ``weiss``. Because
    eszett expands to two characters the offsets shift, so each folded character
    carries the index of the original character it came from -- that is what
    lets an accepted match be quoted verbatim rather than in folded form.
    """
    folded: list[str] = []
    origin: list[int] = []
    for index, char in enumerate(text):
        if char in ("ß", "ẞ"):  # eszett, capital eszett
            expanded = "ss"
        else:
            decomposed = unicodedata.normalize("NFKD", char)
            expanded = "".join(
                c for c in decomposed if not unicodedata.combining(c)
            ).lower()
        for produced in expanded:
            folded.append(produced)
            origin.append(index)
    return "".join(folded), origin


def _boundary_spans(text: str, term: str) -> list[tuple[int, int]]:
    """Every whole-word occurrence of ``term`` in already-folded ``text``."""
    if not term:
        return []
    pattern = re.compile(rf"(?<![a-z0-9]){re.escape(term)}(?![a-z0-9])")
    return [(m.start(), m.end()) for m in pattern.finditer(text)]


SENTENCE_DELIMITERS = ".;\n"


def sentence_spans(text: str) -> list[tuple[int, int]]:
    """Split text into the segments an organ description may not escape.

    Botanical descriptions are organised one organ per sentence or
    semicolon-delimited clause: "Folia coriacea, subtus fusca. Flores albi."
    Searching for the nearest organ across that boundary attaches the leaf's
    brown to the flower, which is the exact failure the WFO ledger measured at
    38.6% of its rejections. A period between them is a wall.

    Splitting on abbreviations ("2 m.") costs evidence rather than correctness,
    so the conservative boundary is the right one.
    """
    spans: list[tuple[int, int]] = []
    start = 0
    for index, char in enumerate(text):
        if char in SENTENCE_DELIMITERS:
            if index > start:
                spans.append((start, index))
            start = index + 1
    if start < len(text):
        spans.append((start, len(text)))
    return spans


def _containing_span(spans: list[tuple[int, int]], position: int) -> tuple[int, int]:
    for start, end in spans:
        if start <= position < end:
            return start, end
    return position, position


def _nearest_organ(
    span: tuple[int, int],
    organs: list[tuple[int, int, str]],
    bounds: tuple[int, int],
    window: int,
) -> str:
    """Return "floral", "competing" or "" for the organ closest to a colour word.

    Only organ terms inside ``bounds`` -- the colour word's own sentence -- are
    considered. Within it, distance is measured to the nearest whole-word organ
    on either side, and when nothing separates them ("flores rubri fructus", one
    character each way) the competing organ wins: the conservative reading is
    the one the WFO rejection ledger says matters.
    """
    start, end = span
    low, high = bounds
    best_kind = ""
    best_distance = window + 1
    for organ_start, organ_end, kind in organs:
        if organ_start < low or organ_end > high:
            continue
        if organ_start <= start and end <= organ_end:
            distance = 0
        else:
            distance = min(abs(start - organ_end), abs(organ_start - end))
        if distance <= window and (
            distance < best_distance
            or (distance == best_distance and kind == "competing")
        ):
            best_distance = distance
            best_kind = kind
    return best_kind


def _sentence_around(text: str, start: int, end: int, limit: int = 500) -> str:
    """The verbatim sentence containing an accepted match, capped at ``limit``."""
    left = max(
        (text.rfind(mark, 0, start) for mark in (".", ";", "\n")), default=-1
    )
    right_candidates = [
        position
        for position in (text.find(mark, end) for mark in (".", ";", "\n"))
        if position != -1
    ]
    right = min(right_candidates) + 1 if right_candidates else len(text)
    return " ".join(text[left + 1 : right].split())[:limit]


def extract_floral_colour(
    description: str, config: dict[str, Any]
) -> tuple[str, str, str, str]:
    """Return (outcome, normalized_value, matched_terms, verbatim_quote)."""
    if not description.strip():
        return REJECT_NO_TEXT, "", "", ""

    text, origin = fold(description)

    for marker in config["negation_markers"]:
        folded_marker, _ = fold(str(marker))
        if folded_marker and folded_marker in text:
            return REJECT_NEGATED, "", "", ""

    colours = expand_colour_terms(config)
    window = int(config["organ_proximity_chars"])
    spans = sentence_spans(text)

    # Organ positions are found once for the whole document rather than once per
    # colour word: an OCR'd page against two hundred organ terms and three
    # hundred colour forms is otherwise a five-figure number of regex scans.
    organs: list[tuple[int, int, str]] = []
    for kind, key in (
        ("competing", "competing_organ_terms"),
        ("floral", "floral_organ_terms"),
    ):
        for raw in config[key]:
            for organ_start, organ_end in _boundary_spans(text, fold(str(raw))[0]):
                organs.append((organ_start, organ_end, kind))

    values: list[str] = []
    matched: list[str] = []
    quote = ""
    saw_colour = False
    saw_competing = False

    for term in sorted(colours, key=len, reverse=True):
        for span in _boundary_spans(text, term):
            saw_colour = True
            kind = _nearest_organ(
                span, organs, _containing_span(spans, span[0]), window
            )
            if kind == "competing":
                saw_competing = True
                continue
            if kind != "floral":
                continue
            value = colours[term]
            if term not in matched:
                matched.append(term)
            if value not in values:
                values.append(value)
            if not quote:
                quote = _sentence_around(
                    description, origin[span[0]], origin[span[1] - 1] + 1
                )

    if not saw_colour:
        return REJECT_NO_COLOUR, "", "", ""
    if not values:
        return (REJECT_COMPETING if saw_competing else REJECT_NO_ORGAN), "", "", ""
    if len(values) > 1:
        return (
            ACCEPTED,
            str(config.get("multiple_value_result", "multicolored_variable")),
            "+".join(sorted(matched)),
            quote,
        )
    return ACCEPTED, values[0], "+".join(sorted(matched)), quote


def build_target_queue(
    master_taxa: pd.DataFrame,
    config: dict[str, Any],
    resolved_species: set[str] | None = None,
) -> pd.DataFrame:
    """Queue the tail species a protologue lookup should be attempted for.

    The tail is defined by island occupancy rather than record count: a species
    on many islands is not the acquisition problem even when poorly recorded,
    and it is the narrow endemics that carry the coverage-isolation confound.
    """
    selection = config["target_selection"]
    max_islands = int(selection["max_islands"])

    queue = master_taxa.loc[
        pd.to_numeric(master_taxa["n_islands"], errors="coerce").le(max_islands)
    ].copy()
    queue["accepted_species"] = queue["accepted_species"].map(_text)
    if resolved_species:
        queue = queue.loc[~queue["accepted_species"].isin(resolved_species)]

    queue = queue.loc[queue["accepted_species"].ne("")]
    queue["axis"] = str(selection.get("axis", TRAIT_NAME))
    columns = [c for c in ("accepted_species", "genus", "family", "n_islands", "n_records", "axis") if c in queue.columns]
    return (
        queue[columns]
        .sort_values(["n_islands", "n_records"], ascending=[True, True])
        .reset_index(drop=True)
    )


def parse_ipni_citation(payload: Any, species: str, config: dict[str, Any]) -> dict[str, str]:
    """Pull the protologue citation for ``species`` out of an IPNI response.

    IPNI returns matches for related names as well, so the record is only used
    when its name equals the queried species. A citation naming the work but not
    the page cannot locate a scan, so it is rejected rather than half-used.
    """
    results = (payload or {}).get("results") or []
    wanted = _text(species).lower()
    for record in results:
        if _text(record.get("name")).lower() != wanted:
            continue
        citation = {
            "publication": _text(record.get("publication")),
            "page": _text(record.get("referenceCollation") or record.get("page")),
            "year": _text(record.get("publicationYear")),
            "authors": _text(record.get("authors")),
            "ipni_id": _text(record.get("id") or record.get("fqId")),
            "source_citation": _text(record.get("reference")),
        }
        missing = [
            field
            for field in config["citation_index"]["fields_required"]
            if not citation.get(str(field))
        ]
        if missing:
            return {"outcome": REJECT_CITATION_INCOMPLETE, "missing": ",".join(missing)}
        if not citation["source_citation"]:
            citation["source_citation"] = (
                f"{citation['publication']} {citation['page']} ({citation['year']})"
            )
        citation["outcome"] = ""
        return citation
    return {"outcome": REJECT_NO_CITATION}


def _year(value: str) -> int | None:
    match = re.search(r"\b(1[5-9]\d{2}|20\d{2})\b", value or "")
    return int(match.group(1)) if match else None


def select_public_domain_page(
    payload: Any, citation: dict[str, str], config: dict[str, Any]
) -> dict[str, str]:
    """Choose the scanned page for a citation, or say why none may be used.

    Two gates, both required: BHL must declare the item public domain, and the
    year must be at or before the configured cutoff. An item in copyright is
    skipped, never scraped, whatever its scientific value.
    """
    provider = config["scan_provider"]
    permitted = {str(r).strip().lower() for r in provider["public_domain_rights"]}
    cutoff = int(provider["max_copyright_year"])

    items = (payload or {}).get("Result") or []
    if isinstance(items, dict):
        items = [items]

    saw_item = False
    saw_in_copyright = False
    for item in items:
        saw_item = True
        rights = _text(item.get("Rights") or item.get("CopyrightStatus"))
        year = _year(_text(item.get("PublicationDate"))) or _year(citation.get("year", ""))
        if rights.strip().lower() not in permitted or (year is not None and year > cutoff):
            saw_in_copyright = True
            continue

        wanted_page = _text(citation.get("page"))
        for page in item.get("Pages") or []:
            numbers = _text(page.get("PageNumbers") or page.get("PageNumber"))
            if wanted_page and wanted_page not in numbers:
                continue
            ocr = str(page.get("OcrText") or "")
            return {
                "outcome": "",
                "page_id": _text(page.get("PageID")),
                "item_id": _text(item.get("ItemID")),
                "source_url": _text(page.get("PageUrl") or page.get("ItemUrl")),
                "license": rights,
                "publication_year": str(year or ""),
                "ocr_text": ocr,
            }

    if saw_in_copyright and not saw_item:
        return {"outcome": REJECT_IN_COPYRIGHT}
    if saw_in_copyright:
        return {"outcome": REJECT_IN_COPYRIGHT}
    return {"outcome": REJECT_NO_SCAN}


def fetch_citation(species: str, fetch: Fetch, config: dict[str, Any]) -> dict[str, str]:
    """Ask IPNI for the protologue citation of one species."""
    payload = fetch(str(config["citation_index"]["base_url"]), {"q": species})
    return parse_ipni_citation(payload, species, config)


def fetch_page(
    citation: dict[str, str], fetch: Fetch, config: dict[str, Any], api_key: str
) -> dict[str, str]:
    """Ask BHL for a public-domain scan of the page the citation names."""
    payload = fetch(
        str(config["scan_provider"]["base_url"]),
        {
            "op": "PublicationSearch",
            "searchterm": citation.get("publication", ""),
            "searchtype": "F",
            "pages": "t",
            "ocr": "t",
            "apikey": api_key,
            "format": "json",
        },
    )
    return select_public_domain_page(payload, citation, config)


def require_lineage(record: dict[str, Any], config: dict[str, Any]) -> list[str]:
    """Lineage fields missing from a record, which reject it rather than warn."""
    return [
        str(field)
        for field in config["required_lineage_fields"]
        if not _text(record.get(str(field)))
    ]


def acquire(
    queue: pd.DataFrame,
    fetch: Fetch,
    config: dict[str, Any],
    api_key: str,
    retrieval_date: str | None = None,
) -> pd.DataFrame:
    """Walk the queue through IPNI and BHL and emit the full audit ledger."""
    retrieved = retrieval_date or date.today().isoformat()
    rows: list[dict[str, Any]] = []

    for species in queue["accepted_species"].map(_text):
        record: dict[str, Any] = {
            "accepted_species": species,
            "trait_name": TRAIT_NAME,
            "normalized_value": "",
            "matched_colour_terms": "",
            "exact_supporting_quote": "",
            "source_url": "",
            "source_citation": "",
            "license": "",
            "retrieval_date": retrieved,
            "ipni_id": "",
            "bhl_page_id": "",
            "evidence_scope": EVIDENCE_SCOPE,
            "source_type": SOURCE_TYPE,
            "review_status": "unreviewed",
        }

        citation = fetch_citation(species, fetch, config)
        if citation.get("outcome"):
            rows.append({**record, "outcome": citation["outcome"]})
            continue
        record["source_citation"] = citation["source_citation"]
        record["ipni_id"] = citation.get("ipni_id", "")

        page = fetch_page(citation, fetch, config, api_key)
        if page.get("outcome"):
            rows.append({**record, "outcome": page["outcome"]})
            continue
        record["source_url"] = page["source_url"]
        record["license"] = page["license"]
        record["bhl_page_id"] = page["page_id"]

        outcome, value, matched, quote = extract_floral_colour(page["ocr_text"], config)
        record["normalized_value"] = value
        record["matched_colour_terms"] = matched
        record["exact_supporting_quote"] = quote

        if outcome == ACCEPTED:
            missing = require_lineage(record, config)
            if missing:
                outcome = REJECT_MISSING_LINEAGE
                record["normalized_value"] = ""
        rows.append({**record, "outcome": outcome})

    return pd.DataFrame(rows)


@app.command("queue")
def queue_command(
    master_taxa_csv: Path = typer.Option(
        Path("data/v2/staging/gbif/collected/island_taxa.csv"),
        "--master-taxa-csv",
        exists=True,
    ),
    output_csv: Path = typer.Option(..., "--output-csv"),
    resolved_species_csv: Path | None = typer.Option(
        None, "--resolved-species-csv", exists=True
    ),
    config_path: Path = typer.Option(
        Path("config/protologue_acquisition.yml"), "--config", exists=True
    ),
) -> None:
    """Write the tail-species queue a protologue run should work through."""
    config = load_config(config_path)
    master = pd.read_csv(master_taxa_csv)
    resolved: set[str] | None = None
    if resolved_species_csv is not None:
        resolved = set(pd.read_csv(resolved_species_csv)["accepted_species"].map(_text))

    queue = build_target_queue(master, config, resolved)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    queue.to_csv(output_csv, index=False)
    typer.echo(
        json.dumps(
            {
                "n_master_species": int(len(master)),
                "n_queued": int(len(queue)),
                "max_islands": int(config["target_selection"]["max_islands"]),
                "n_already_resolved_excluded": int(len(resolved)) if resolved else 0,
            },
            indent=2,
            sort_keys=True,
        )
    )


@app.command("run")
def run(
    queue_csv: Path = typer.Option(..., "--queue-csv", exists=True),
    output_dir: Path = typer.Option(..., "--output-dir"),
    api_key_env: str = typer.Option("BHL_API_KEY", "--api-key-env"),
    limit: int = typer.Option(0, "--limit", help="0 walks the whole queue"),
    config_path: Path = typer.Option(
        Path("config/protologue_acquisition.yml"), "--config", exists=True
    ),
) -> None:
    """Resolve queued species to protologues and mine their floral statements."""
    import os

    import requests

    config = load_config(config_path)
    api_key = os.environ.get(api_key_env, "")
    if not api_key:
        raise typer.BadParameter(f"{api_key_env} is not set; BHL requires an API key")

    session = requests.Session()

    def fetch(url: str, params: dict[str, str]) -> Any:
        response = session.get(url, params=params, timeout=60)
        response.raise_for_status()
        return response.json()

    queue = pd.read_csv(queue_csv, dtype=str, keep_default_na=False)
    if limit > 0:
        queue = queue.head(limit)

    audit = acquire(queue, fetch, config, api_key)
    candidates = audit.loc[audit["outcome"].eq(ACCEPTED)]

    output_dir.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output_dir / "protologue_acquisition_audit.csv.gz", index=False)
    candidates.to_csv(output_dir / "protologue_colour_candidates.csv.gz", index=False)

    summary = {
        "version": "1.0",
        "n_species_attempted": int(len(audit)),
        "outcomes": {str(k): int(v) for k, v in audit["outcome"].value_counts().items()},
        "n_candidates": int(len(candidates)),
        "n_species_with_a_candidate": int(candidates["accepted_species"].nunique())
        if len(candidates)
        else 0,
        "interpretation": (
            "Candidates are unreviewed species-direct protologue statements "
            "carrying the verbatim quote and full source lineage. They are not "
            "accepted evidence, and a species without a candidate is unresolved "
            "rather than colourless -- most commonly because no public-domain "
            "scan of its protologue exists."
        ),
    }
    (output_dir / "protologue_acquisition_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    app()
