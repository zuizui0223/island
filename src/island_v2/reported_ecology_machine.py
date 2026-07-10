"""Extract reported reproductive and visitor ecology candidates from source text.

This module is the machine-only bridge for LLM or rule-based source-excerpt
extraction. It never decides curated truth. It normalizes source-backed
statements into the reported-ecology candidate CSV consumed by
``trait_multimodal_candidates``.
"""

from __future__ import annotations

import json
import re
import time
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.trait_source_discovery import openalex_abstract

app = typer.Typer(add_completion=False, no_args_is_help=True)

OPENALEX_WORKS_URL = "https://api.openalex.org/works"

OUTPUT_COLUMNS = [
    "accepted_species",
    "trait_name",
    "candidate_value",
    "source_url",
    "source_citation",
    "source_excerpt",
    "source_type",
    "extraction_model",
    "prompt_id",
    "evidence_scope",
    "matched_term",
    "confidence",
]

TEXT_SOURCE_COLUMNS = [
    "accepted_species",
    "source_text",
    "source_url",
    "source_citation",
    "source_type",
    "evidence_scope",
]

TRAIT_ALIASES = {
    "self_incompatibility": {
        "si": "SI",
        "self_incompatible": "SI",
        "self-incompatible": "SI",
        "self incompatible": "SI",
        "sc": "SC",
        "self_compatible": "SC",
        "self-compatible": "SC",
        "self compatible": "SC",
        "mixed": "mixed_or_variable",
        "mixed_or_variable": "mixed_or_variable",
        "unresolved": "unresolved",
    },
    "autonomous_selfing_capacity": {
        "absent": "absent",
        "no_autonomous_selfing": "absent",
        "delayed": "delayed",
        "delayed_selfing": "delayed",
        "facilitated": "facilitated",
        "autonomous": "autonomous",
        "autonomous_selfing": "autonomous",
        "spontaneous_selfing": "autonomous",
        "autogamy": "autonomous",
        "autogamous": "autonomous",
        "mixed": "mixed_or_variable",
        "mixed_or_variable": "mixed_or_variable",
        "unresolved": "unresolved",
    },
    "mating_system": {
        "outcrossing": "predominantly_outcrossing",
        "predominantly_outcrossing": "predominantly_outcrossing",
        "obligate_outcrossing": "predominantly_outcrossing",
        "mixed": "mixed_mating",
        "mixed_mating": "mixed_mating",
        "facultative_selfing": "mixed_mating",
        "selfing": "predominantly_selfing",
        "predominantly_selfing": "predominantly_selfing",
        "obligate_selfing": "obligate_selfing",
        "unresolved": "unresolved",
    },
    "pollen_vector_mode": {
        "biotic": "biotic",
        "animal": "biotic",
        "insect": "biotic",
        "entomophily": "biotic",
        "wind": "abiotic_wind",
        "wind_pollinated": "abiotic_wind",
        "anemophily": "abiotic_wind",
        "mixed": "mixed",
        "unresolved": "unresolved",
    },
    "pollination_functional_guild": {
        "bumblebee": "bumblebees",
        "bumblebees": "bumblebees",
        "bombus": "bumblebees",
        "bee": "other_bees",
        "bees": "other_bees",
        "fly": "flies",
        "flies": "flies",
        "butterfly": "butterflies",
        "butterflies": "butterflies",
        "moth": "moths",
        "moths": "moths",
        "bird": "birds",
        "birds": "birds",
        "bat": "bats",
        "bats": "bats",
        "beetle": "beetles_wasps_ants",
        "wasp": "beetles_wasps_ants",
        "ant": "beetles_wasps_ants",
        "mixed": "mixed_or_generalist",
        "generalist": "mixed_or_generalist",
        "unresolved": "unresolved",
    },
    "dichogamy": {
        "protandry": "protandry",
        "protandrous": "protandry",
        "protogyny": "protogyny",
        "protogynous": "protogyny",
        "absent": "absent",
        "variable": "variable",
        "unresolved": "unresolved",
    },
    "cleistogamy": {
        "absent": "absent",
        "facultative": "facultative",
        "obligate": "obligate",
        "unresolved": "unresolved",
    },
    "sex_system": {
        "hermaphroditic": "hermaphroditic",
        "bisexual": "hermaphroditic",
        "monoecious": "monoecious",
        "dioecious": "dioecious",
        "gynodioecious": "gynodioecious",
        "androdioecious": "androdioecious",
        "polygamous": "polygamous_or_other",
        "unresolved": "unresolved",
    },
}

PATTERNS: list[tuple[str, str, str]] = [
    ("self_incompatibility", "SI", r"\bself[- ]incompatib(?:le|ility)\b"),
    ("self_incompatibility", "SC", r"\bself[- ]compatib(?:le|ility)\b"),
    ("autonomous_selfing_capacity", "delayed", r"\bdelayed selfing\b"),
    ("autonomous_selfing_capacity", "autonomous", r"\b(?:autonomous|spontaneous) selfing\b"),
    ("autonomous_selfing_capacity", "autonomous", r"\bautogam(?:y|ous)\b"),
    ("autonomous_selfing_capacity", "absent", r"\bno autonomous selfing\b"),
    ("mating_system", "mixed_mating", r"\bmixed mating\b"),
    ("mating_system", "mixed_mating", r"\bfacultative selfing\b"),
    ("mating_system", "predominantly_outcrossing", r"\boutcrossing\b"),
    ("mating_system", "obligate_selfing", r"\bobligate selfing\b"),
    ("pollen_vector_mode", "abiotic_wind", r"\bwind[- ]pollinat(?:ed|ion)\b"),
    ("pollen_vector_mode", "biotic", r"\binsect[- ]pollinat(?:ed|ion)\b"),
    ("pollen_vector_mode", "biotic", r"\bentomophil(?:y|ous)\b"),
    ("pollination_functional_guild", "bumblebees", r"\b(?:Bombus|bumblebee|bumblebees)\b"),
    ("pollination_functional_guild", "other_bees", r"\bbee[- ]pollinated\b"),
    ("pollination_functional_guild", "other_bees", r"\bpollinated by bees\b"),
    ("pollination_functional_guild", "flies", r"\bfly[- ]pollinated\b"),
    ("pollination_functional_guild", "flies", r"\bmyophil(?:y|ous)\b"),
    ("pollination_functional_guild", "birds", r"\bbird[- ]pollinated\b"),
    ("pollination_functional_guild", "birds", r"\bornithophil(?:y|ous)\b"),
    ("pollination_functional_guild", "bats", r"\bbat[- ]pollinated\b"),
    ("pollination_functional_guild", "bats", r"\bchiropterophil(?:y|ous)\b"),
    ("dichogamy", "protandry", r"\bprotandr(?:y|ous)\b"),
    ("dichogamy", "protogyny", r"\bprotogyn(?:y|ous)\b"),
    ("cleistogamy", "facultative", r"\bfacultative cleistogam(?:y|ous)\b"),
    ("cleistogamy", "obligate", r"\bobligate cleistogam(?:y|ous)\b"),
    ("sex_system", "dioecious", r"\bdioecious\b"),
    ("sex_system", "monoecious", r"\bmonoecious\b"),
    ("sex_system", "gynodioecious", r"\bgynodioecious\b"),
]


@app.callback()
def main() -> None:
    """Extract source-backed ecology candidates for machine-only workflows."""


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _excerpt(text: str, start: int, end: int, pad: int = 90) -> str:
    lo = max(0, start - pad)
    hi = min(len(text), end + pad)
    return re.sub(r"\s+", " ", text[lo:hi]).strip()


def _strip_markup(text: object) -> str:
    return re.sub(r"<[^>]+>", " ", str(text or ""))


def _has_local_negation(text: str, start: int) -> bool:
    local = text[max(0, start - 35):start].casefold()
    return bool(re.search(r"\b(?:not|no|without|lack(?:s|ed|ing)?)\b", local))


def load_ontology(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    return data if isinstance(data, dict) else {}


def allowed_values(ontology: dict[str, Any], trait_name: str) -> set[str]:
    traits = ontology.get("traits") or {}
    trait = traits.get(trait_name) or {}
    return {str(value) for value in trait.get("allowed_values") or []}


def normalize_value(trait_name: str, raw_value: str, ontology: dict[str, Any]) -> str:
    value = _text(raw_value)
    allowed = allowed_values(ontology, trait_name)
    if value in allowed:
        return value
    alias = TRAIT_ALIASES.get(trait_name, {}).get(value.casefold().replace(" ", "_"))
    if alias in allowed:
        return str(alias)
    alias = TRAIT_ALIASES.get(trait_name, {}).get(value.casefold())
    if alias in allowed:
        return str(alias)
    return ""


def validate_candidate_row(row: dict[str, Any], ontology: dict[str, Any]) -> tuple[dict[str, str] | None, str]:
    trait_name = _text(row.get("trait_name"))
    candidate_value = normalize_value(trait_name, _text(row.get("candidate_value")), ontology)
    source_excerpt = _text(
        row.get("source_excerpt") or row.get("evidence_excerpt") or row.get("raw_description")
    )
    species = _text(row.get("accepted_species"))
    if not species:
        return None, "missing accepted_species"
    if not trait_name:
        return None, "missing trait_name"
    if not candidate_value:
        return None, f"invalid candidate_value for {trait_name}"
    if not source_excerpt:
        return None, "missing source_excerpt"
    output = {
        "accepted_species": species,
        "trait_name": trait_name,
        "candidate_value": candidate_value,
        "source_url": _text(row.get("source_url")),
        "source_citation": _text(row.get("source_citation")),
        "source_excerpt": source_excerpt,
        "source_type": _text(row.get("source_type")) or "source_text",
        "extraction_model": _text(row.get("extraction_model")) or "external_or_rule_extractor",
        "prompt_id": _text(row.get("prompt_id")) or "reported_ecology_v1",
        "evidence_scope": _text(row.get("evidence_scope")) or "species_direct",
        "matched_term": _text(row.get("matched_term")),
        "confidence": _text(row.get("confidence")) or "machine_extracted",
    }
    return output, ""


def extract_candidates_from_text_sources(
    sources: pd.DataFrame,
    ontology: dict[str, Any],
    *,
    extraction_model: str = "rule_based_reported_ecology_v1",
    prompt_id: str = "reported_ecology_rules_v1",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, str]] = []
    holdouts: list[dict[str, str]] = []
    for record in sources.fillna("").to_dict("records"):
        species = _text(record.get("accepted_species"))
        text = _text(
            record.get("source_text") or record.get("text") or record.get("abstract") or record.get("raw_description")
        )
        if not species or not text:
            continue
        seen: set[tuple[str, str, str]] = set()
        for trait_name, value, pattern in PATTERNS:
            for match in re.finditer(pattern, text, flags=re.IGNORECASE):
                if value not in {"absent"} and _has_local_negation(text, match.start()):
                    holdouts.append(
                        {
                            "accepted_species": species,
                            "trait_name": trait_name,
                            "matched_term": match.group(0),
                            "holdout_reason": "local negation before matched term",
                        }
                    )
                    continue
                final_value = normalize_value(trait_name, value, ontology)
                key = (species, trait_name, final_value)
                if not final_value or key in seen:
                    continue
                seen.add(key)
                rows.append(
                    {
                        "accepted_species": species,
                        "trait_name": trait_name,
                        "candidate_value": final_value,
                        "source_url": _text(record.get("source_url")),
                        "source_citation": _text(record.get("source_citation")),
                        "source_excerpt": _excerpt(text, match.start(), match.end()),
                        "source_type": _text(record.get("source_type")) or "source_text",
                        "extraction_model": extraction_model,
                        "prompt_id": prompt_id,
                        "evidence_scope": _text(record.get("evidence_scope")) or "species_direct",
                        "matched_term": match.group(0),
                        "confidence": "rule_match",
                    }
                )
    candidates = pd.DataFrame(rows, columns=OUTPUT_COLUMNS).drop_duplicates()
    holdout = pd.DataFrame(
        holdouts,
        columns=["accepted_species", "trait_name", "matched_term", "holdout_reason"],
    ).drop_duplicates()
    return candidates.reset_index(drop=True), holdout.reset_index(drop=True)


def read_jsonl_candidates(path: Path, ontology: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        try:
            payload = json.loads(line)
        except json.JSONDecodeError as exc:
            errors.append({"line": str(line_number), "error": f"invalid json: {exc}"})
            continue
        items = payload.get("candidates") if isinstance(payload, dict) else None
        if items is None:
            items = [payload]
        for item in items:
            if not isinstance(item, dict):
                errors.append({"line": str(line_number), "error": "candidate is not an object"})
                continue
            row, error = validate_candidate_row(item, ontology)
            if row is None:
                errors.append({"line": str(line_number), "error": error})
            else:
                rows.append(row)
    candidates = pd.DataFrame(rows, columns=OUTPUT_COLUMNS).drop_duplicates()
    error_table = pd.DataFrame(errors, columns=["line", "error"])
    return candidates.reset_index(drop=True), error_table


def openalex_text_sources(
    species_df: pd.DataFrame,
    max_taxa: int,
    per_page: int,
    pause_seconds: float,
) -> pd.DataFrame:
    import httpx

    terms = (
        "self-incompatibility self-compatible autonomous selfing mating system "
        "pollination pollinator floral visitor pollen vector"
    )
    selected = (
        species_df["accepted_species"].fillna("").astype(str).str.strip().drop_duplicates().head(max_taxa)
    )
    rows: list[dict[str, str]] = []
    with httpx.Client(
        timeout=45.0,
        follow_redirects=True,
        headers={"User-Agent": "island-floral-v2 reported-ecology-machine"},
    ) as client:
        for species in selected:
            if not species:
                continue
            response = client.get(
                OPENALEX_WORKS_URL,
                params={"search": f'"{species}" {terms}', "per_page": per_page},
            )
            response.raise_for_status()
            for work in response.json().get("results") or []:
                title = _text(_strip_markup(work.get("display_name") or work.get("title")))
                abstract = _text(openalex_abstract(work))
                source_text = _text(f"{title}. {abstract}")
                if not source_text:
                    continue
                lowered = source_text.casefold()
                scope = "species_direct" if species.casefold() in lowered else "species_indirect"
                doi = _text(work.get("doi"))
                rows.append(
                    {
                        "accepted_species": species,
                        "source_text": source_text,
                        "source_url": doi or _text(work.get("id")),
                        "source_citation": _text(f"{title} {work.get('publication_year') or ''}"),
                        "source_type": "openalex_title_abstract",
                        "evidence_scope": scope,
                    }
                )
            if pause_seconds > 0:
                time.sleep(pause_seconds)
    return pd.DataFrame(rows, columns=TEXT_SOURCE_COLUMNS)


def write_outputs(
    output_dir: Path,
    candidates: pd.DataFrame,
    holdouts: pd.DataFrame,
    summary_extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    candidates.to_csv(output_dir / "reported_ecology_candidates.csv", index=False)
    holdouts.to_csv(output_dir / "reported_ecology_holdouts.csv", index=False)
    summary = {
        "n_reported_ecology_candidates": int(len(candidates)),
        "n_species": int(candidates["accepted_species"].nunique()) if not candidates.empty else 0,
        "n_holdouts": int(len(holdouts)),
        "trait_counts": {
            str(k): int(v)
            for k, v in candidates["trait_name"].value_counts().sort_index().items()
        }
        if not candidates.empty
        else {},
        "human_review_required": False,
        "policy": "Machine-only source-excerpt candidates; preserve source excerpt and keep method tier downstream.",
    }
    if summary_extra:
        summary.update(summary_extra)
    (output_dir / "reported_ecology_extraction_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return summary


@app.command("extract-text")
def extract_text(
    source_csv: Path = typer.Option(..., exists=True, help="CSV with accepted_species and source_text."),
    output_dir: Path = typer.Option(...),
    ontology_path: Path = typer.Option(Path("config/trait_ontology.yml"), exists=True),
) -> None:
    """Extract rule-based reported ecology candidates from source-text CSV."""
    ontology = load_ontology(ontology_path)
    sources = pd.read_csv(source_csv, dtype=str).fillna("")
    candidates, holdouts = extract_candidates_from_text_sources(sources, ontology)
    summary = write_outputs(output_dir, candidates, holdouts, {"n_source_rows": int(len(sources))})
    typer.echo(
        f"Extracted {summary['n_reported_ecology_candidates']} reported ecology candidate(s) "
        f"from {summary['n_source_rows']} source row(s)."
    )


@app.command("import-jsonl")
def import_jsonl(
    jsonl_path: Path = typer.Option(..., exists=True, help="JSONL with LLM/external candidates."),
    output_dir: Path = typer.Option(...),
    ontology_path: Path = typer.Option(Path("config/trait_ontology.yml"), exists=True),
) -> None:
    """Import and validate LLM/external reported-ecology JSONL candidates."""
    ontology = load_ontology(ontology_path)
    candidates, errors = read_jsonl_candidates(jsonl_path, ontology)
    summary = write_outputs(output_dir, candidates, errors, {"n_import_errors": int(len(errors))})
    typer.echo(
        f"Imported {summary['n_reported_ecology_candidates']} reported ecology candidate(s); "
        f"{summary['n_import_errors']} import error(s)."
    )


@app.command("extract-openalex")
def extract_openalex(
    species_csv: Path = typer.Option(..., exists=True, help="Species CSV with accepted_species."),
    output_dir: Path = typer.Option(...),
    max_taxa: int = typer.Option(50, min=1, max=500),
    per_page: int = typer.Option(5, min=1, max=25),
    pause_seconds: float = typer.Option(0.05, min=0.0),
    ontology_path: Path = typer.Option(Path("config/trait_ontology.yml"), exists=True),
) -> None:
    """Fetch OpenAlex title/abstract text and run the rule-based ecology baseline."""
    species = pd.read_csv(species_csv, dtype=str).fillna("")
    if "accepted_species" not in species.columns:
        raise typer.BadParameter("species CSV must contain accepted_species")
    sources = openalex_text_sources(species, max_taxa, per_page, pause_seconds)
    ontology = load_ontology(ontology_path)
    candidates, holdouts = extract_candidates_from_text_sources(
        sources,
        ontology,
        extraction_model="rule_based_openalex_reported_ecology_v1",
        prompt_id="openalex_reported_ecology_rules_v1",
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    sources.to_csv(output_dir / "openalex_text_sources.csv", index=False)
    summary = write_outputs(
        output_dir,
        candidates,
        holdouts,
        {"n_source_rows": int(len(sources)), "source": "OpenAlex title/abstract"},
    )
    typer.echo(
        f"Extracted {summary['n_reported_ecology_candidates']} reported ecology candidate(s) "
        f"from {summary['n_source_rows']} OpenAlex source row(s)."
    )


if __name__ == "__main__":
    app()
