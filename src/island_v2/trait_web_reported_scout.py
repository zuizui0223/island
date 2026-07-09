"""Stage-1B (A): Wikipedia / Wikidata scout for `reported` trait candidates.

Primary literature is thin for obscure island endemics, but many have a Wikipedia
article or a Wikidata item. Per the evidence policy, an **explicit statement in any
source** — including web/DB — is kept as a ``reported`` candidate (not a proxy). This
scout resolves each species to its Wikidata item and English Wikipedia article, then
emits reported M0/M1 candidates wherever a controlled-vocabulary term appears verbatim,
retaining full provenance.

Fail-closed: every row is `candidate_class=reported`, unreviewed, with
`source_type`, `source_url`, `raw_description`, `evidence_scope`,
`wild_or_cultivated`, `review_status`. No value is decided; inference with no explicit
statement is out of scope here (that is a `*_proxy`, produced separately).
"""

from __future__ import annotations

import json
import re
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

WIKIDATA_API = "https://www.wikidata.org/w/api.php"
WIKIPEDIA_API = "https://en.wikipedia.org/w/api.php"
JsonGetter = Callable[[str, dict[str, Any]], dict[str, Any]]

REPORTED_COLUMNS = [
    "accepted_species",
    "trait_layer",
    "trait_name",
    "provisional_candidate_value",
    "matched_term",
    "candidate_class",
    "source",
    "source_type",
    "source_url",
    "raw_description",
    "evidence_scope",
    "wild_or_cultivated",
    "review_status",
]

PROHIBITED_OUTPUT_COLUMNS = {
    "flower_primary_color", "floral_symmetry", "floral_form", "trait", "trait_value",
    "trait_state", "native_status", "establishment_means", "is_native",
    "bombus_applicability", "applicability", "analysis_included", "review_status_accepted",
}


@app.callback()
def main() -> None:
    """Scout Wikipedia/Wikidata for unreviewed `reported` trait candidates."""


def load_keyword_layers(paths: dict[str, Path]) -> tuple[dict[str, Any], dict[str, str]]:
    """Merge per-layer keyword configs; return (traits map, trait_name -> layer)."""
    merged: dict[str, Any] = {}
    trait_layer: dict[str, str] = {}
    for layer, path in paths.items():
        config = yaml.safe_load(path.read_text(encoding="utf-8"))
        if not isinstance(config, dict) or "traits" not in config:
            raise typer.BadParameter(f"keyword config {path} must contain a traits mapping")
        for trait_name, canon_map in config["traits"].items():
            merged[trait_name] = canon_map
            trait_layer[trait_name] = layer
    return {"traits": merged}, trait_layer


def _excerpt(text: str, start: int, end: int, pad: int = 60) -> str:
    lo = max(0, start - pad)
    hi = min(len(text), end + pad)
    return re.sub(r"\s+", " ", text[lo:hi]).strip()


def extract_reported(
    species: str,
    text: str,
    source: str,
    source_type: str,
    source_url: str,
    evidence_scope: str,
    config: dict[str, Any],
    trait_layer: dict[str, str],
) -> list[dict[str, str]]:
    """Emit reported candidate rows for controlled terms found verbatim in text."""
    lowered = (text or "").lower()
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for trait_name, canon_map in config.get("traits", {}).items():
        for canonical, surface_terms in (canon_map or {}).items():
            for term in surface_terms or []:
                term_l = str(term).lower().strip()
                if not term_l:
                    continue
                m = re.search(rf"(?<!\w){re.escape(term_l)}(?!\w)", lowered)
                if not m:
                    continue
                key = (trait_name, canonical, term_l)
                if key in seen:
                    continue
                seen.add(key)
                rows.append(
                    {
                        "accepted_species": species,
                        "trait_layer": trait_layer.get(trait_name, ""),
                        "trait_name": str(trait_name),
                        "provisional_candidate_value": str(canonical),
                        "matched_term": term_l,
                        "candidate_class": "reported",
                        "source": source,
                        "source_type": source_type,
                        "source_url": source_url,
                        "raw_description": _excerpt(text, m.start(), m.end()),
                        "evidence_scope": evidence_scope,
                        "wild_or_cultivated": "unknown",
                        "review_status": "unreviewed",
                    }
                )
    return rows


def fetch_wikidata(getter: JsonGetter, name: str) -> dict[str, str]:
    """Resolve a species name to a Wikidata QID, English description, enwiki title."""
    search = getter(
        WIKIDATA_API,
        {"action": "wbsearchentities", "search": name, "language": "en",
         "uselang": "en", "type": "item", "format": "json", "limit": 1},
    )
    hits = search.get("search") or []
    if not hits:
        return {}
    qid = str(hits[0].get("id") or "")
    description = str(hits[0].get("description") or "")
    enwiki_title = ""
    if qid:
        ent = getter(
            WIKIDATA_API,
            {"action": "wbgetentities", "ids": qid, "props": "sitelinks",
             "sitefilter": "enwiki", "format": "json"},
        )
        sitelinks = (((ent.get("entities") or {}).get(qid) or {}).get("sitelinks") or {})
        enwiki_title = str((sitelinks.get("enwiki") or {}).get("title") or "")
    return {"qid": qid, "description": description, "enwiki_title": enwiki_title}


def fetch_wikipedia_extract(getter: JsonGetter, title: str) -> tuple[str, str, str]:
    """Return (plain-text extract, article URL, resolved title) for a Wikipedia page."""
    payload = getter(
        WIKIPEDIA_API,
        {"action": "query", "prop": "extracts", "explaintext": 1, "redirects": 1,
         "titles": title, "format": "json", "exlimit": 1},
    )
    pages = ((payload.get("query") or {}).get("pages") or {})
    for page in pages.values():
        if "missing" in page:
            continue
        resolved = str(page.get("title") or title)
        extract = str(page.get("extract") or "")
        url = "https://en.wikipedia.org/wiki/" + resolved.replace(" ", "_")
        return extract, url, resolved
    return "", "", ""


def _scope_for(species: str, resolved_title: str) -> str:
    """species_direct when the article is clearly the species page, else species_indirect."""
    tokens = [t for t in re.split(r"\s+", species.strip()) if t.isalpha()]
    epithet = tokens[1].lower() if len(tokens) >= 2 else ""
    title_l = resolved_title.lower()
    if species.lower() in title_l or (epithet and epithet in title_l):
        return "species_direct"
    return "species_indirect"


def scout_species(
    getter: JsonGetter, species: str, config: dict[str, Any], trait_layer: dict[str, str]
) -> tuple[list[dict[str, str]], list[str]]:
    rows: list[dict[str, str]] = []
    errors: list[str] = []
    enwiki_title = ""
    try:
        wd = fetch_wikidata(getter, species)
    except Exception as exc:  # noqa: BLE001 - one source failing must not abort the taxon
        wd = {}
        errors.append(f"wikidata:{species}:{exc}")
    if wd:
        enwiki_title = wd.get("enwiki_title", "")
        # The Wikidata item description itself can carry an explicit term.
        if wd.get("description"):
            rows.extend(
                extract_reported(
                    species, wd["description"], "wikidata", "wikidata",
                    f"https://www.wikidata.org/wiki/{wd.get('qid', '')}",
                    "species_direct", config, trait_layer,
                )
            )
    title = enwiki_title or species
    try:
        extract, url, resolved = fetch_wikipedia_extract(getter, title)
    except Exception as exc:  # noqa: BLE001
        return rows, errors + [f"wikipedia:{species}:{exc}"]
    if extract:
        rows.extend(
            extract_reported(
                species, extract, "wikipedia", "wikipedia", url,
                _scope_for(species, resolved), config, trait_layer,
            )
        )
    return rows, errors


def scout_web_reported(
    species_df: pd.DataFrame,
    getter: JsonGetter,
    config: dict[str, Any],
    trait_layer: dict[str, str],
    max_taxa: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if "accepted_species" not in species_df.columns:
        raise typer.BadParameter("species table must have an 'accepted_species' column")
    names = [n for n in dict.fromkeys(species_df["accepted_species"].fillna("").astype(str).str.strip()) if n][:max_taxa]

    rows: list[dict[str, str]] = []
    errors: list[str] = []
    species_with: set[str] = set()
    for name in names:
        species_rows, species_errors = scout_species(getter, name, config, trait_layer)
        rows.extend(species_rows)
        errors.extend(species_errors)
        if species_rows:
            species_with.add(name)

    frame = pd.DataFrame(rows, columns=REPORTED_COLUMNS)
    leaked = PROHIBITED_OUTPUT_COLUMNS.intersection(frame.columns)
    if leaked:
        raise typer.BadParameter(f"web reported scout leaked prohibited columns: {sorted(leaked)}")

    def _by(column: str) -> dict[str, int]:
        return {str(k): int(v) for k, v in frame[column].value_counts().sort_index().items()} if len(frame) else {}

    report = {
        "note": (
            "Unreviewed `reported` candidates from Wikipedia/Wikidata (explicit statements "
            "only). Web/DB-derived explicit statements are kept with full provenance; "
            "candidate_class=reported. Inference without an explicit statement is out of "
            "scope here (handled as *_proxy). No value is decided."
        ),
        "source_lane": "stage_1b_web_reported_scout",
        "sources": ["wikipedia", "wikidata"],
        "n_taxa_queried": len(names),
        "n_species_with_reported": len(species_with),
        "n_species_zero_reported": len(names) - len(species_with),
        "n_reported_candidates": int(len(frame)),
        "n_by_source": _by("source"),
        "n_by_trait_layer": _by("trait_layer"),
        "n_by_trait_name": _by("trait_name"),
        "n_lookup_errors": len(errors),
        "lookup_errors": errors[:20],
    }
    return frame, report


def _httpx_getter() -> JsonGetter:
    import httpx

    client = httpx.Client(
        timeout=45.0,
        follow_redirects=True,
        headers={"User-Agent": "island-floral-v2/0.1 (web-reported-scout; contact via repo)"},
    )

    def getter(url: str, params: dict[str, Any]) -> dict[str, Any]:
        response = client.get(url, params=params)
        response.raise_for_status()
        return response.json()

    return getter


@app.command("scout")
def scout(
    species_csv: Path = typer.Option(..., exists=True, help="Species CSV (needs accepted_species)."),
    output_dir: Path = typer.Option(...),
    m0_config: Path = typer.Option(Path("config/m0_descriptive_keywords.yml")),
    m1_config: Path = typer.Option(Path("config/m1_reported_keywords.yml")),
    max_taxa: int = typer.Option(50, min=1, max=500, help="Bounded species to query."),
) -> None:
    """Write unreviewed `reported` M0/M1 candidates from Wikipedia/Wikidata."""
    config, trait_layer = load_keyword_layers({"M0": m0_config, "M1": m1_config})
    species_df = pd.read_csv(species_csv, dtype=str).fillna("")
    frame, report = scout_web_reported(species_df, _httpx_getter(), config, trait_layer, max_taxa)
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "web_reported_candidates.csv", index=False)
    (output_dir / "web_reported_scout_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"{report['n_reported_candidates']} unreviewed reported candidate(s) for "
        f"{report['n_species_with_reported']}/{report['n_taxa_queried']} species "
        "(Wikipedia/Wikidata). No value decided."
    )


if __name__ == "__main__":
    app()
