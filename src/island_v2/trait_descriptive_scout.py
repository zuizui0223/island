"""Stage-1B: descriptive / flora scout for M0 floral-phenotype candidates.

Obscure island endemics often have little primary literature (Stage-1A returns
nothing on-species) but do have descriptive coverage. This scout retrieves free-text
species descriptions -- starting with the GBIF species descriptions API (free, no key)
-- and proposes **M0 candidate values only where a controlled-vocabulary term appears
verbatim**, each with a source excerpt.

Fail-closed: every row is an unreviewed candidate, not a decided value. Values live in
a long ``trait_name`` / ``provisional_candidate_value`` format (never a finalized wide
column), multiple colours yield multiple candidates (nothing is collapsed to a single
primary), and a blank is never biological absence. Additional sources (POWO/WCVP, flora
and horticulture databases, Wikidata/Wikipedia) can be added behind the same getter.
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

GBIF_API = "https://api.gbif.org/v1"
JsonGetter = Callable[[str, dict[str, Any]], dict[str, Any]]

CANDIDATE_COLUMNS = [
    "accepted_species",
    "trait_layer",
    "trait_name",
    "provisional_candidate_value",
    "matched_term",
    "candidate_kind",
    "evidence_excerpt",
    "source",
    "description_type",
    "source_url",
    "source_citation",
    "review_status",
]

# These finalized-value column names must never appear: values stay in the long
# trait_name / provisional_candidate_value form, keeping the output a candidate DB.
PROHIBITED_OUTPUT_COLUMNS = {
    "flower_primary_color", "floral_symmetry", "floral_form", "trait", "trait_value",
    "trait_state", "native_status", "establishment_means", "is_native",
    "bombus_applicability", "applicability", "analysis_included", "review_status_accepted",
}


@app.callback()
def main() -> None:
    """Scout descriptive sources for unreviewed M0 floral-phenotype candidates."""


def load_keywords(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or "traits" not in config:
        raise typer.BadParameter("M0 keyword config must contain a traits mapping")
    return config


def _excerpt(text: str, start: int, end: int, pad: int = 60) -> str:
    lo = max(0, start - pad)
    hi = min(len(text), end + pad)
    return re.sub(r"\s+", " ", text[lo:hi]).strip()


def extract_candidates(
    species: str,
    text: str,
    source_meta: dict[str, str],
    config: dict[str, Any],
) -> list[dict[str, str]]:
    """Return M0 candidate rows for every controlled term found verbatim in text."""
    lowered = (text or "").lower()
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for trait_name, canon_map in config.get("traits", {}).items():
        for canonical, surface_terms in (canon_map or {}).items():
            for term in surface_terms or []:
                term_l = str(term).lower().strip()
                if not term_l:
                    continue
                # whole-word (or phrase) match
                for m in re.finditer(rf"(?<!\w){re.escape(term_l)}(?!\w)", lowered):
                    key = (trait_name, canonical, term_l)
                    if key in seen:
                        break  # one candidate per (trait, value, term) per description
                    seen.add(key)
                    rows.append(
                        {
                            "accepted_species": species,
                            "trait_layer": "M0",
                            "trait_name": str(trait_name),
                            "provisional_candidate_value": str(canonical),
                            "matched_term": term_l,
                            "candidate_kind": "descriptive_keyword_match",
                            "evidence_excerpt": _excerpt(text, m.start(), m.end()),
                            "source": source_meta.get("source", "gbif_species_description"),
                            "description_type": source_meta.get("description_type", ""),
                            "source_url": source_meta.get("source_url", ""),
                            "source_citation": source_meta.get("source_citation", ""),
                            "review_status": "unreviewed",
                        }
                    )
                    break
    return rows


def fetch_species_key(getter: JsonGetter, name: str) -> str:
    payload = getter(f"{GBIF_API}/species/match", {"name": name})
    key = payload.get("usageKey") or payload.get("nubKey") or ""
    return str(key) if key else ""


def fetch_descriptions(getter: JsonGetter, key: str, limit: int) -> list[dict[str, Any]]:
    payload = getter(f"{GBIF_API}/species/{key}/descriptions", {"limit": limit})
    return list(payload.get("results") or [])[:limit]


def scout_species(
    getter: JsonGetter,
    species: str,
    config: dict[str, Any],
    max_descriptions: int,
) -> tuple[list[dict[str, str]], list[str]]:
    """Retrieve GBIF descriptions for one species and extract M0 candidates."""
    wanted_types = {str(t).lower() for t in config.get("description_types", [])}
    errors: list[str] = []
    try:
        key = fetch_species_key(getter, species)
    except Exception as exc:  # noqa: BLE001 - one source failing must not abort the run
        return [], [f"match:{species}:{exc}"]
    if not key:
        return [], []
    try:
        descriptions = fetch_descriptions(getter, key, max_descriptions)
    except Exception as exc:  # noqa: BLE001
        return [], [f"descriptions:{species}:{exc}"]

    rows: list[dict[str, str]] = []
    for doc in descriptions:
        dtype = str(doc.get("type") or "").lower()
        if wanted_types and dtype and dtype not in wanted_types:
            continue
        text = str(doc.get("description") or "")
        if not text.strip():
            continue
        meta = {
            "source": "gbif_species_description",
            "description_type": dtype,
            "source_url": f"https://www.gbif.org/species/{key}",
            "source_citation": str(doc.get("source") or ""),
        }
        rows.extend(extract_candidates(species, text, meta, config))
    return rows, errors


def scout_descriptive(
    species_df: pd.DataFrame,
    getter: JsonGetter,
    config: dict[str, Any],
    max_taxa: int,
    max_descriptions: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Build the M0 candidate table and a run report for the given species."""
    if "accepted_species" not in species_df.columns:
        raise typer.BadParameter("species table must have an 'accepted_species' column")
    names = (
        species_df["accepted_species"].fillna("").astype(str).str.strip()
    )
    names = [n for n in dict.fromkeys(names) if n][:max_taxa]

    rows: list[dict[str, str]] = []
    errors: list[str] = []
    species_with_candidates: set[str] = set()
    for name in names:
        species_rows, species_errors = scout_species(getter, name, config, max_descriptions)
        rows.extend(species_rows)
        errors.extend(species_errors)
        if species_rows:
            species_with_candidates.add(name)

    frame = pd.DataFrame(rows, columns=CANDIDATE_COLUMNS)
    leaked = PROHIBITED_OUTPUT_COLUMNS.intersection(frame.columns)
    if leaked:
        raise typer.BadParameter(f"descriptive scout leaked prohibited columns: {sorted(leaked)}")

    def _by(column: str) -> dict[str, int]:
        return {str(k): int(v) for k, v in frame[column].value_counts().sort_index().items()} if len(frame) else {}

    report = {
        "note": (
            "Unreviewed M0 descriptive candidates only. A candidate is a verbatim "
            "controlled-vocabulary match in a species description, with a source "
            "excerpt; it is not a decided trait value. Multiple colours yield multiple "
            "candidates (nothing collapsed); a blank is never biological absence."
        ),
        "source_lane": "stage_1b_descriptive_scout",
        "primary_source": "gbif_species_descriptions",
        "n_taxa_queried": len(names),
        "n_species_with_candidates": len(species_with_candidates),
        "n_species_zero_candidates": len(names) - len(species_with_candidates),
        "n_candidates": int(len(frame)),
        "n_by_trait_name": _by("trait_name"),
        "n_by_candidate_value": _by("provisional_candidate_value"),
        "n_lookup_errors": len(errors),
        "lookup_errors": errors[:20],
    }
    return frame, report


def _httpx_getter() -> JsonGetter:
    import httpx

    client = httpx.Client(
        timeout=45.0,
        follow_redirects=True,
        headers={"User-Agent": "island-floral-v2/0.1 (descriptive-scout)"},
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
    config_path: Path = typer.Option(Path("config/m0_descriptive_keywords.yml")),
    max_taxa: int = typer.Option(50, min=1, max=500, help="Bounded species to query."),
    max_descriptions: int = typer.Option(20, min=1, max=100, help="Bounded descriptions per species."),
) -> None:
    """Write unreviewed M0 descriptive candidates from GBIF species descriptions."""
    config = load_keywords(config_path)
    species_df = pd.read_csv(species_csv, dtype=str).fillna("")
    frame, report = scout_descriptive(species_df, _httpx_getter(), config, max_taxa, max_descriptions)
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "m0_descriptive_candidates.csv", index=False)
    (output_dir / "m0_descriptive_scout_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"{report['n_candidates']} unreviewed M0 candidate(s) for "
        f"{report['n_species_with_candidates']}/{report['n_taxa_queried']} species. "
        "No trait value was decided."
    )


if __name__ == "__main__":
    app()
