"""Category-first trait extraction for the v1 island-floral analysis.

Collection keeps the richer nine-column LLM output requested by the user. The
analysis table is derived afterwards using the binary categories from v1:
plain/conspicuous colour, generalized/specialized floral form, and
self-compatible/self-incompatible reproductive state.

The module deliberately separates raw extraction from recoding. Low-confidence
or inference-only rows remain available for sensitivity analyses rather than
being silently promoted into the primary evidence layer.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

PROMPT_PLACEHOLDER = "{{SPECIES_LIST}}"
DEFAULT_PROMPT_PATH = Path("config/v1_category_trait_prompt.txt")

OUTPUT_COLUMNS = [
    "species",
    "flower_color",
    "flower_shape",
    "pollination_guild",
    "pollination_notes",
    "mating_system",
    "self_incompatibility",
    "evidence_type",
    "confidence",
]

POLLINATION_GUILDS = {
    "bees",
    "bumblebees",
    "flies",
    "butterflies",
    "moths",
    "birds",
    "wind",
    "self",
    "mixed",
    "unknown",
}
MATING_SYSTEMS = {
    "obligate_outcrossing",
    "mixed_mating",
    "mainly_selfing",
    "obligate_selfing",
    "unknown",
}
SELF_INCOMPATIBILITY = {"SI", "SC", "likely_SI", "likely_SC", "unknown"}
EVIDENCE_TYPES = {"field_study", "review", "flora", "horticulture", "inference"}
CONFIDENCE_VALUES = {"high", "medium", "low"}

DERIVED_COLUMNS = [
    "v1_color_binary",
    "v1_form_binary",
    "v1_self_compatibility_binary",
    "v1_animal_pollinated",
    "v1_mating_binary",
    "evidence_tier",
    "primary_analysis_eligible",
]

COLOUR_CONSPICUOUS = re.compile(
    r"\b(?:blue|purple|violet|lavender|lilac|red|pink|magenta|scarlet|crimson|"
    r"yellow|orange|gold|golden)\b",
    flags=re.IGNORECASE,
)
COLOUR_PLAIN = re.compile(
    r"\b(?:white|cream|creamy|ivory|green|greenish|brown|brownish|dull|"
    r"inconspicuous|colourless|colorless)\b",
    flags=re.IGNORECASE,
)

FORM_SPECIALIZED = re.compile(
    r"\b(?:zygomorphic|bilateral|bilabiate|two[- ]lipped|papilionaceous|"
    r"tubular|tube[- ]shaped|salverform|funnel|trumpet|urceolate|urn[- ]shaped|"
    r"spurred)\b",
    flags=re.IGNORECASE,
)
FORM_GENERALIZED = re.compile(
    r"\b(?:actinomorphic|radial|radially symmetric|open|bowl|dish|cup[- ]shaped|"
    r"rotate|brush|puff|composite head|capitulum|head[- ]like|reduced|wind[- ]type)\b",
    flags=re.IGNORECASE,
)


class CategoryTraitValidationError(ValueError):
    """Raised when an LLM result violates the category-first output contract."""


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _species_column(columns: list[str], requested: str | None = None) -> str:
    if requested:
        if requested not in columns:
            raise typer.BadParameter(
                f"species column {requested!r} not found; available={columns}"
            )
        return requested
    for candidate in ("accepted_species", "species", "scientificName", "canonicalName"):
        if candidate in columns:
            return candidate
    raise typer.BadParameter(
        "could not find a species column; expected one of accepted_species, species, "
        "scientificName, canonicalName"
    )


def load_prompt_template(path: Path = DEFAULT_PROMPT_PATH) -> str:
    if not path.exists():
        raise typer.BadParameter(f"prompt template does not exist: {path}")
    text = path.read_text(encoding="utf-8")
    if text.count(PROMPT_PLACEHOLDER) != 1:
        raise typer.BadParameter(
            f"prompt template must contain exactly one {PROMPT_PLACEHOLDER!r} placeholder"
        )
    return text


def render_prompt(species: list[str], template: str) -> str:
    cleaned = [_text(value) for value in species if _text(value)]
    if not cleaned:
        raise typer.BadParameter("cannot render an empty species batch")
    species_block = "\n".join(f"- {name}" for name in cleaned)
    return template.replace(PROMPT_PLACEHOLDER, species_block)


def prepare_prompt_batches(
    *,
    species_csv: Path,
    output_dir: Path,
    prompt_path: Path = DEFAULT_PROMPT_PATH,
    species_column: str | None = None,
    batch_size: int = 25,
    start_index: int = 0,
    max_species: int | None = None,
) -> dict[str, Any]:
    """Create deterministic prompt packets without calling an external LLM API."""
    if not species_csv.exists():
        raise typer.BadParameter(f"species CSV does not exist: {species_csv}")
    if batch_size < 1:
        raise typer.BadParameter("batch_size must be positive")
    if start_index < 0:
        raise typer.BadParameter("start_index cannot be negative")
    if max_species is not None and max_species < 1:
        raise typer.BadParameter("max_species must be positive when supplied")

    table = pd.read_csv(species_csv, dtype=str).fillna("")
    column = _species_column(list(table.columns), species_column)
    ordered: list[str] = []
    seen: set[str] = set()
    for value in table[column]:
        name = _text(value)
        if name and name not in seen:
            seen.add(name)
            ordered.append(name)
    selected = ordered[start_index:]
    if max_species is not None:
        selected = selected[:max_species]
    if not selected:
        raise typer.BadParameter("no species remain after applying selection limits")

    output_dir.mkdir(parents=True, exist_ok=True)
    template = load_prompt_template(prompt_path)
    packet_rows: list[dict[str, Any]] = []
    for batch_number, offset in enumerate(range(0, len(selected), batch_size), start=1):
        batch_species = selected[offset : offset + batch_size]
        batch_id = f"batch_{batch_number:05d}"
        batch_dir = output_dir / batch_id
        batch_dir.mkdir(parents=True, exist_ok=True)
        (batch_dir / "prompt.txt").write_text(
            render_prompt(batch_species, template),
            encoding="utf-8",
        )
        pd.DataFrame({"species": batch_species}).to_csv(
            batch_dir / "expected_species.csv",
            index=False,
        )
        pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(
            batch_dir / "result_template.csv",
            index=False,
        )
        packet_rows.append(
            {
                "batch_id": batch_id,
                "n_species": len(batch_species),
                "first_species": batch_species[0],
                "last_species": batch_species[-1],
                "prompt_path": str(batch_dir / "prompt.txt"),
                "expected_species_path": str(batch_dir / "expected_species.csv"),
                "result_template_path": str(batch_dir / "result_template.csv"),
            }
        )

    packet_index = pd.DataFrame(packet_rows)
    packet_index.to_csv(output_dir / "prompt_packet_index.csv", index=False)
    manifest = {
        "version": "1.0",
        "species_csv": str(species_csv),
        "species_column": column,
        "prompt_path": str(prompt_path),
        "batch_size": batch_size,
        "start_index": start_index,
        "n_species_available": len(ordered),
        "n_species_selected": len(selected),
        "n_batches": len(packet_rows),
        "output_columns": OUTPUT_COLUMNS,
        "policy": (
            "Prompt packets preserve category-rich collection output. Binary v1 "
            "analysis variables are derived only after result validation."
        ),
    }
    (output_dir / "prompt_packet_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest


def _normalise_evidence_types(value: object) -> str:
    raw = _text(value)
    if not raw:
        return ""
    parts = [
        _text(part)
        for part in re.split(r"[,;|/]", raw)
        if _text(part)
    ]
    unknown = sorted(set(parts).difference(EVIDENCE_TYPES))
    if unknown:
        raise CategoryTraitValidationError(
            f"invalid evidence_type value(s): {unknown}; allowed={sorted(EVIDENCE_TYPES)}"
        )
    ordered = [value for value in sorted(EVIDENCE_TYPES) if value in set(parts)]
    return "|".join(ordered)


def _validate_enum(
    table: pd.DataFrame,
    column: str,
    allowed: set[str],
) -> None:
    invalid = sorted(set(table[column]).difference(allowed))
    if invalid:
        raise CategoryTraitValidationError(
            f"invalid {column} value(s): {invalid}; allowed={sorted(allowed)}"
        )


def validate_result_table(
    table: pd.DataFrame,
    *,
    expected_species: list[str] | None = None,
) -> pd.DataFrame:
    """Validate and normalize one category-first LLM result table."""
    actual_columns = list(table.columns)
    if actual_columns != OUTPUT_COLUMNS:
        raise CategoryTraitValidationError(
            "result columns must match the prompt exactly and in order; "
            f"expected={OUTPUT_COLUMNS}, actual={actual_columns}"
        )
    result = table.copy().fillna("")
    for column in OUTPUT_COLUMNS:
        result[column] = result[column].map(_text)
    result["evidence_type"] = result["evidence_type"].map(_normalise_evidence_types)

    if result["species"].eq("").any():
        raise CategoryTraitValidationError("species cannot be blank")
    duplicates = result.loc[result["species"].duplicated(keep=False), "species"].unique()
    if len(duplicates):
        raise CategoryTraitValidationError(
            f"duplicate species rows: {sorted(duplicates.tolist())}"
        )

    _validate_enum(result, "pollination_guild", POLLINATION_GUILDS)
    _validate_enum(result, "mating_system", MATING_SYSTEMS)
    _validate_enum(result, "self_incompatibility", SELF_INCOMPATIBILITY)
    _validate_enum(result, "confidence", CONFIDENCE_VALUES)

    blank_evidence = result["evidence_type"].eq("")
    if blank_evidence.any():
        species = result.loc[blank_evidence, "species"].tolist()
        raise CategoryTraitValidationError(
            f"evidence_type cannot be blank: {species[:10]}"
        )

    insufficient = (
        result["flower_color"].str.casefold().eq("unknown")
        & result["flower_shape"].str.casefold().eq("unknown")
        & result["pollination_guild"].eq("unknown")
        & result["mating_system"].eq("unknown")
        & result["self_incompatibility"].eq("unknown")
    )
    invalid_confidence = insufficient & ~result["confidence"].eq("low")
    if invalid_confidence.any():
        species = result.loc[invalid_confidence, "species"].tolist()
        raise CategoryTraitValidationError(
            "rows with insufficient information must have low confidence: "
            f"{species[:10]}"
        )

    if expected_species is not None:
        expected = [_text(value) for value in expected_species if _text(value)]
        expected_set = set(expected)
        observed_set = set(result["species"])
        missing = sorted(expected_set.difference(observed_set))
        extra = sorted(observed_set.difference(expected_set))
        if missing or extra:
            raise CategoryTraitValidationError(
                f"species mismatch; missing={missing[:10]}, extra={extra[:10]}"
            )
        order = {species: index for index, species in enumerate(expected)}
        result = (
            result.assign(_order=result["species"].map(order))
            .sort_values("_order")
            .drop(columns="_order")
            .reset_index(drop=True)
        )
    return result


def _color_binary(value: object) -> str:
    text = _text(value)
    if not text or text.casefold() == "unknown":
        return "unknown"
    conspicuous = bool(COLOUR_CONSPICUOUS.search(text))
    plain = bool(COLOUR_PLAIN.search(text))
    if conspicuous:
        return "conspicuous"
    if plain:
        return "plain"
    return "unknown"


def _form_binary(value: object) -> str:
    text = _text(value)
    if not text or text.casefold() == "unknown":
        return "unknown"
    specialized = bool(FORM_SPECIALIZED.search(text))
    generalized = bool(FORM_GENERALIZED.search(text))
    if specialized and not generalized:
        return "specialized"
    if generalized and not specialized:
        return "generalized"
    if specialized and generalized:
        # A tubular or zygomorphic restriction takes precedence in the v1 split.
        return "specialized"
    return "unknown"


def _self_binary(value: object) -> str:
    mapping = {
        "SC": "self_compatible",
        "likely_SC": "self_compatible",
        "SI": "self_incompatible",
        "likely_SI": "self_incompatible",
        "unknown": "unknown",
    }
    return mapping.get(_text(value), "unknown")


def _animal_pollinated(value: object) -> str:
    guild = _text(value)
    if guild in {"bees", "bumblebees", "flies", "butterflies", "moths", "birds"}:
        return "yes"
    if guild in {"wind", "self"}:
        return "no"
    if guild == "mixed":
        return "mixed"
    return "unknown"


def _mating_binary(value: object) -> str:
    mapping = {
        "obligate_outcrossing": "outcrossing",
        "mixed_mating": "mixed",
        "mainly_selfing": "selfing",
        "obligate_selfing": "selfing",
        "unknown": "unknown",
    }
    return mapping.get(_text(value), "unknown")


def _evidence_tier(row: pd.Series) -> str:
    evidence = set(_text(row["evidence_type"]).split("|"))
    confidence = _text(row["confidence"])
    direct = evidence.intersection({"field_study", "review", "flora"})
    if direct and confidence in {"high", "medium"}:
        return "primary"
    if "horticulture" in evidence and confidence in {"high", "medium"}:
        return "secondary"
    return "inference_or_low"


def derive_v1_categories(validated: pd.DataFrame) -> pd.DataFrame:
    """Append v1 binary analysis variables without altering raw trait fields."""
    result = validated.copy()
    result["v1_color_binary"] = result["flower_color"].map(_color_binary)
    result["v1_form_binary"] = result["flower_shape"].map(_form_binary)
    result["v1_self_compatibility_binary"] = result["self_incompatibility"].map(
        _self_binary
    )
    result["v1_animal_pollinated"] = result["pollination_guild"].map(
        _animal_pollinated
    )
    result["v1_mating_binary"] = result["mating_system"].map(_mating_binary)
    result["evidence_tier"] = result.apply(_evidence_tier, axis=1)
    result["primary_analysis_eligible"] = (
        result["evidence_tier"].eq("primary")
        & ~result["v1_color_binary"].eq("unknown")
        & ~result["v1_form_binary"].eq("unknown")
        & ~result["v1_self_compatibility_binary"].eq("unknown")
    )
    return result[[*OUTPUT_COLUMNS, *DERIVED_COLUMNS]]


def ingest_result_csv(
    *,
    result_csv: Path,
    expected_species_csv: Path | None,
    output_csv: Path,
) -> dict[str, Any]:
    if not result_csv.exists():
        raise typer.BadParameter(f"result CSV does not exist: {result_csv}")
    expected_species: list[str] | None = None
    if expected_species_csv is not None:
        if not expected_species_csv.exists():
            raise typer.BadParameter(
                f"expected species CSV does not exist: {expected_species_csv}"
            )
        expected_table = pd.read_csv(expected_species_csv, dtype=str).fillna("")
        column = _species_column(list(expected_table.columns))
        expected_species = expected_table[column].map(_text).tolist()

    raw = pd.read_csv(result_csv, dtype=str, keep_default_na=False)
    validated = validate_result_table(raw, expected_species=expected_species)
    derived = derive_v1_categories(validated)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    derived.to_csv(output_csv, index=False)
    report = {
        "version": "1.0",
        "result_csv": str(result_csv),
        "expected_species_csv": (
            str(expected_species_csv) if expected_species_csv is not None else None
        ),
        "output_csv": str(output_csv),
        "n_species": len(derived),
        "coverage": {
            "flower_color_nonunknown": int(
                ~derived["v1_color_binary"].eq("unknown").sum()
            ),
            "flower_shape_nonunknown": int(
                ~derived["v1_form_binary"].eq("unknown").sum()
            ),
            "pollination_guild_nonunknown": int(
                ~derived["pollination_guild"].eq("unknown").sum()
            ),
            "mating_system_nonunknown": int(
                ~derived["mating_system"].eq("unknown").sum()
            ),
            "self_incompatibility_nonunknown": int(
                ~derived["self_incompatibility"].eq("unknown").sum()
            ),
            "primary_analysis_eligible": int(
                derived["primary_analysis_eligible"].astype(bool).sum()
            ),
        },
    }
    report_path = output_csv.with_suffix(".manifest.json")
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return report


@app.command("prepare")
def prepare_command(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    prompt_path: Path = typer.Option(DEFAULT_PROMPT_PATH, exists=True),
    species_column: str | None = typer.Option(None),
    batch_size: int = typer.Option(25, min=1, max=100),
    start_index: int = typer.Option(0, min=0),
    max_species: int | None = typer.Option(None, min=1),
) -> None:
    report = prepare_prompt_batches(
        species_csv=species_csv,
        output_dir=output_dir,
        prompt_path=prompt_path,
        species_column=species_column,
        batch_size=batch_size,
        start_index=start_index,
        max_species=max_species,
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


@app.command("ingest")
def ingest_command(
    result_csv: Path = typer.Option(..., exists=True),
    output_csv: Path = typer.Option(...),
    expected_species_csv: Path | None = typer.Option(None),
) -> None:
    try:
        report = ingest_result_csv(
            result_csv=result_csv,
            expected_species_csv=expected_species_csv,
            output_csv=output_csv,
        )
    except CategoryTraitValidationError as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
