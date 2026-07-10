"""Build review queues for visual-trait and LLM-reported ecology candidates.

This module is intentionally lightweight. It does not download images, run an
image model, or call an LLM. It accepts the candidate outputs from those steps
as CSV files and converts them into the same fail-closed human review queue used
by the validation pilot.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2 import trait_candidate_review_queue as review_queue

app = typer.Typer(add_completion=False, no_args_is_help=True)

VISUAL_TRAITS = {
    "flower_primary_color",
    "floral_symmetry",
    "floral_form",
    "tube_depth_class",
    "flower_size_class",
    "inflorescence_display",
}

REPORTED_ECOLOGY_TRAITS = {
    "pollen_vector_mode",
    "pollination_functional_guild",
    "self_incompatibility",
    "autonomous_selfing_capacity",
    "mating_system",
    "herkogamy",
    "dichogamy",
    "cleistogamy",
    "sex_system",
}

VISUAL_TEMPLATE_COLUMNS = [
    "accepted_species",
    "trait_name",
    "candidate_value",
    "image_url",
    "image_page_url",
    "source_citation",
    "model_name",
    "model_version",
    "model_confidence",
    "prompt_id",
    "evidence_note",
]

REPORTED_ECOLOGY_TEMPLATE_COLUMNS = [
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
]


@app.callback()
def main() -> None:
    """Prepare multimodal candidate rows for human review."""


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _read_csv(path: Path | None) -> pd.DataFrame:
    if path is None or not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, dtype=str).fillna("")


def _pick(record: dict[str, Any], *names: str) -> str:
    for name in names:
        value = _text(record.get(name))
        if value:
            return value
    return ""


def _join_parts(*parts: str) -> str:
    return " | ".join(part for part in parts if part)


def _trait_layer(trait_name: str) -> str:
    if trait_name in VISUAL_TRAITS:
        return "M0"
    if trait_name == "pollination_functional_guild":
        return "M2"
    return "M1"


def _basis_metadata(record: dict[str, Any], *names: str) -> str:
    parts = []
    for name in names:
        value = _text(record.get(name))
        if value:
            parts.append(f"{name}={value}")
    return "|".join(parts)


def _sort_queue(queue: pd.DataFrame) -> pd.DataFrame:
    if queue.empty:
        return queue
    queue = queue.drop_duplicates("review_candidate_id", keep="first")
    queue["_priority_sort"] = queue["review_priority"].astype(int)
    queue = queue.sort_values(
        ["_priority_sort", "accepted_species", "trait_name", "candidate_value", "review_candidate_id"],
        kind="stable",
    )
    return queue.drop(columns=["_priority_sort"]).reset_index(drop=True)


def rows_from_visual_candidates(table: pd.DataFrame) -> list[dict[str, str]]:
    """Return queue rows for visible floral traits from an image-candidate table."""
    rows: list[dict[str, str]] = []
    for record in table.fillna("").to_dict("records"):
        species = _text(record.get("accepted_species"))
        trait_name = _text(record.get("trait_name"))
        value = _pick(record, "candidate_value", "visual_candidate_value", "provisional_candidate_value")
        image_url = _pick(record, "image_url", "identifier", "media_url", "source_url")
        if not species or not trait_name or not value or trait_name not in VISUAL_TRAITS:
            continue
        page_url = _pick(record, "image_page_url", "source_page_url", "source_url")
        note = _join_parts(
            "Image-derived visible-trait candidate; review morphology only.",
            _text(record.get("evidence_note")),
        )
        rows.append(
            review_queue._base_row(  # noqa: SLF001
                source_lane="visual_trait_candidate",
                accepted_species=species,
                trait_layer=_trait_layer(trait_name),
                trait_name=trait_name,
                candidate_value=value,
                candidate_class=_text(record.get("candidate_class")) or "visual_observed",
                candidate_kind=_text(record.get("candidate_kind")) or "image_model_candidate",
                evidence_scope=_text(record.get("evidence_scope")) or "image_visible_trait",
                source=_text(record.get("source")) or "image_manifest",
                source_type=_text(record.get("source_type")) or "image",
                source_url=image_url,
                source_citation=_join_parts(_text(record.get("source_citation")), page_url),
                matched_term=_text(record.get("matched_term")),
                basis_traits=_basis_metadata(
                    record,
                    "model_name",
                    "model_version",
                    "model_confidence",
                    "prompt_id",
                ),
                basis_family=_text(record.get("basis_family")),
                raw_description=note,
                source_review_status=_text(record.get("review_status")) or "unreviewed_visual",
            )
        )
    return rows


def rows_from_reported_ecology_candidates(table: pd.DataFrame) -> list[dict[str, str]]:
    """Return queue rows for source-excerpted reproductive and visitor evidence."""
    rows: list[dict[str, str]] = []
    for record in table.fillna("").to_dict("records"):
        species = _text(record.get("accepted_species"))
        trait_name = _text(record.get("trait_name"))
        value = _pick(record, "candidate_value", "reported_candidate_value", "provisional_candidate_value")
        excerpt = _pick(record, "source_excerpt", "raw_description", "evidence_excerpt")
        if not species or not trait_name or not value or trait_name not in REPORTED_ECOLOGY_TRAITS:
            continue
        rows.append(
            review_queue._base_row(  # noqa: SLF001
                source_lane="llm_reported_ecology",
                accepted_species=species,
                trait_layer=_trait_layer(trait_name),
                trait_name=trait_name,
                candidate_value=value,
                candidate_class=_text(record.get("candidate_class")) or "reported",
                candidate_kind=_text(record.get("candidate_kind")) or "llm_text_extraction",
                evidence_scope=_text(record.get("evidence_scope")) or "species_direct",
                source=_text(record.get("source")) or "llm_reported_ecology",
                source_type=_text(record.get("source_type")) or "literature_text",
                source_url=_text(record.get("source_url")),
                source_citation=_text(record.get("source_citation")),
                matched_term=_text(record.get("matched_term")),
                basis_traits=_basis_metadata(record, "extraction_model", "prompt_id"),
                basis_family=_text(record.get("basis_family")),
                raw_description=excerpt,
                source_review_status=_text(record.get("review_status")) or "unreviewed_llm",
            )
        )
    return rows


def build_multimodal_review_queue(
    visual_candidates: pd.DataFrame,
    reported_ecology_candidates: pd.DataFrame,
) -> pd.DataFrame:
    """Return a review queue without deciding any trait value."""
    rows = []
    if not visual_candidates.empty:
        rows.extend(rows_from_visual_candidates(visual_candidates))
    if not reported_ecology_candidates.empty:
        rows.extend(rows_from_reported_ecology_candidates(reported_ecology_candidates))
    queue = pd.DataFrame(rows, columns=review_queue.OUTPUT_COLUMNS)
    return _sort_queue(queue)


def template_tables() -> dict[str, pd.DataFrame]:
    """Return empty input templates with one illustrative row each."""
    visual = pd.DataFrame(
        [
            {
                "accepted_species": "Example species",
                "trait_name": "flower_primary_color",
                "candidate_value": "red_pink",
                "image_url": "https://example.test/image.jpg",
                "image_page_url": "https://example.test/specimen",
                "source_citation": "Image provider and license",
                "model_name": "manual_or_vision_model",
                "model_version": "",
                "model_confidence": "",
                "prompt_id": "",
                "evidence_note": "Visible petals appear red-pink; needs human review.",
            }
        ],
        columns=VISUAL_TEMPLATE_COLUMNS,
    )
    reported = pd.DataFrame(
        [
            {
                "accepted_species": "Example species",
                "trait_name": "self_incompatibility",
                "candidate_value": "SC",
                "source_url": "https://example.test/paper",
                "source_citation": "Author year",
                "source_excerpt": "The source text explicitly states self-compatibility.",
                "source_type": "paper_abstract_or_full_text",
                "extraction_model": "llm_name",
                "prompt_id": "reported_ecology_v1",
                "evidence_scope": "species_direct",
            }
        ],
        columns=REPORTED_ECOLOGY_TEMPLATE_COLUMNS,
    )
    return {
        "visual_trait_candidates_template.csv": visual,
        "reported_ecology_candidates_template.csv": reported,
    }


@app.command("templates")
def templates(output_dir: Path = typer.Option(..., help="Directory for input CSV templates.")) -> None:
    """Write CSV templates for visual and reported-ecology candidate inputs."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, table in template_tables().items():
        table.to_csv(output_dir / name, index=False)
    typer.echo(f"Wrote {len(template_tables())} template file(s) to {output_dir}.")


@app.command("build")
def build(
    output_dir: Path = typer.Option(..., help="Directory for multimodal review queue artifacts."),
    visual_candidates_csv: Path | None = typer.Option(
        None,
        exists=True,
        help="Optional image-derived visible-trait candidate CSV.",
    ),
    reported_ecology_candidates_csv: Path | None = typer.Option(
        None,
        exists=True,
        help="Optional LLM source-excerpted reproductive/visitor candidate CSV.",
    ),
    ontology_path: Path = typer.Option(
        Path("config/trait_ontology.yml"),
        exists=True,
        help="Trait ontology for validation.",
    ),
) -> None:
    """Write a review queue for visual-trait and reported-ecology candidates."""
    queue = build_multimodal_review_queue(
        _read_csv(visual_candidates_csv),
        _read_csv(reported_ecology_candidates_csv),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    queue_path = output_dir / "trait_candidate_review_queue.csv"
    queue.to_csv(queue_path, index=False)
    ontology = yaml.safe_load(ontology_path.read_text(encoding="utf-8"))
    summary = review_queue.write_review_outputs(queue, output_dir, ontology=ontology)
    summary["input_policy"] = {
        "visual": "Image rows are visible-trait candidates only; no reproductive or pollinator state is inferred from images.",
        "reported_ecology": "LLM rows must preserve source excerpts and remain review candidates until accepted by a reviewer.",
    }
    (output_dir / "trait_candidate_review_queue_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"Built {summary['n_review_candidates']} multimodal review candidate(s) "
        f"across {summary['n_species']} species."
    )


if __name__ == "__main__":
    app()
