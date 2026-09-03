"""Evidence-grounded taxonomic cascade for flower and reproductive traits.

The operational priority is trait yield, not measurement. Direct source-backed
evidence covers a small fraction of the master species, so this cascade fills
every remaining species-trait by descending a taxonomic evidence ladder:

    species_direct -> synonym_direct -> genus_inference -> family_inference
    -> unresolved_no_evidence

Every fill records ``fill_tier``, ``evidence_scope`` and ``confidence`` so a
low-resolution fill is always separable for sensitivity analysis and never
masquerades as accepted direct evidence. Fills are written to staging, never to
curated. A cascade fill is candidate coverage, never a biological absence, and
one trait is never derived from another (self-compatibility never fills
autonomous selfing). Inference tiers retain supporting taxa, support size,
winner share, and the full value distribution. A taxonomic inference is emitted
only when configured support and consensus thresholds are met; tied or otherwise
unqualified evidence remains explicitly unresolved. No dataset-wide modal value
is ever assigned to a species without a genus or family evidence link.
The mutually exclusive ``reported_value`` and ``inferred_value`` channels make
that separation explicit even for consumers that do not inspect ``fill_tier``.
"""

from __future__ import annotations

import glob as globlib
import hashlib
import json
import os
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.angiosperm_scope import classify_scope, load_config as load_scope_config
from island_v2.search_enabled_llm_campaign import _parse_csv_row as parse_llm_result

app = typer.Typer(add_completion=False, no_args_is_help=True)

EMPTY_VALUES = {"", "unknown", "na", "n/a", "none", "unresolved"}

LLM_TRAIT_MAP = {
    "flower_color": "flower_primary_color",
    "flower_shape": "floral_form",
    "pollination_guild": "pollination_functional_guild",
    "mating_system": "mating_system",
    "self_incompatibility": "self_incompatibility",
}

ALLMASTER_SHAPE_MAP = {
    "actinomorphic / radially symmetric": ("floral_symmetry", "actinomorphic"),
    "zygomorphic / bilaterally symmetric": ("floral_symmetry", "zygomorphic"),
    "asymmetric": ("floral_symmetry", "asymmetric"),
    "rotate / open": ("floral_form", "open_radial"),
    "open / bowl / cup-shaped": ("floral_form", "open_radial"),
    "bell-shaped": ("floral_form", "bell_campanulate"),
    "tubular": ("floral_form", "tubular"),
    "salverform": ("floral_form", "salverform"),
    "funnel-shaped": ("floral_form", "funnel_trumpet"),
    "funnel / trumpet-shaped": ("floral_form", "funnel_trumpet"),
    "urn-shaped": ("floral_form", "urn_urceolate"),
    "globose / spherical flower head": ("floral_form", "composite_head"),
    "composite head / capitulum": ("floral_form", "composite_head"),
    "papilionaceous": ("floral_form", "papilionaceous"),
    "spurred": ("floral_form", "spurred"),
    "spike / elongated inflorescence": (
        "inflorescence_display",
        "raceme_spike_panicle",
    ),
    # Bilabiate is reported morphology, but recoding it as zygomorphy would be
    # a derived inference. Retain only that a form was explicitly described.
    "bilabiate / two-lipped": ("floral_form", "other_described"),
}


def _map_allmaster_value(field: object, value: object) -> tuple[str, str]:
    """Map the legacy all-master ledger only where a v2 state is deterministic."""
    field_text = _text(field).casefold()
    value_text = _text(value).casefold()
    if field_text == "flower_color":
        collapsed = {
            "green/brown/inconspicuous": "green_brown_inconspicuous",
            "yellow/orange": "yellow_orange",
            "red/pink": "red_pink",
            "blue/purple": "blue_purple",
            "multicolored/variable": "multicolored_variable",
        }
        if value_text in collapsed:
            return "flower_primary_color", collapsed[value_text]
        palette = {
            "white": "white",
            "cream": "white",
            "green": "green_brown_inconspicuous",
            "brown": "green_brown_inconspicuous",
            "yellow": "yellow_orange",
            "orange": "yellow_orange",
            "red": "red_pink",
            "pink": "red_pink",
            "blue": "blue_purple",
            "purple": "blue_purple",
            "violet": "blue_purple",
            "black": "other_described",
        }
        categories = {palette[token.strip()] for token in value_text.split(",") if token.strip() in palette}
        if not categories:
            return "", ""
        if len(categories) == 1:
            return "flower_primary_color", next(iter(categories))
        return "flower_primary_color", "multicolored_variable"
    if field_text == "flower_shape":
        return ALLMASTER_SHAPE_MAP.get(value_text, ("", ""))
    if field_text == "mating_system":
        mapped = {
            "obligate_outcrossing": "predominantly_outcrossing",
            "mixed_mating": "mixed_mating",
            "mainly_selfing": "predominantly_selfing",
            "obligate_selfing": "obligate_selfing",
        }.get(value_text, "")
        return ("mating_system", mapped) if mapped else ("", "")
    if field_text == "self_incompatibility" and value_text.upper() in {"SI", "SC"}:
        return "self_incompatibility", value_text.upper()
    # Legacy pollinator guild is intentionally not imported into the exhaustive
    # flower/reproductive target matrix.
    return "", ""


def _normalize_direct_value(trait_name: object, value: object) -> str:
    """Normalize known source vocabularies without inventing a biological state."""
    trait = _text(trait_name)
    raw = _text(value)
    folded = raw.casefold()
    if trait in {"inflorescence_display"}:
        if raw.startswith("[") and raw.endswith("]"):
            try:
                values = json.loads(raw)
            except json.JSONDecodeError:
                return ""
            if isinstance(values, list):
                # Inflorescence mode may provide multi-value arrays.
                # Keep strict single-state ontology entries only.
                normalized = {_normalize_direct_value(trait, item) for item in values}
                normalized.discard("")
                if len(normalized) == 1:
                    return next(iter(normalized))
                return ""
            if isinstance(values, str):
                raw = values
                folded = raw.casefold()
        elif isinstance(value, str):
            raw = value.replace("|", ",")
            folded = raw.casefold()
    if trait == "flower_primary_color":
        if folded in {
            "white",
            "green_brown_inconspicuous",
            "yellow_orange",
            "red_pink",
            "blue_purple",
            "multicolored_variable",
            "other_described",
        }:
            return folded
        mapped_trait, mapped_value = _map_allmaster_value("flower_color", raw.replace("|", ","))
        return mapped_value if mapped_trait else ""
    if trait == "floral_form":
        canonical = {
            "open_radial",
            "bell_campanulate",
            "tubular",
            "salverform",
            "funnel_trumpet",
            "urn_urceolate",
            "brush_puff",
            "composite_head",
            "papilionaceous",
            "spurred",
            "reduced_wind",
            "other_described",
        }
        if folded in canonical:
            return folded
        aliases = {
            "rotate": "open_radial",
            "campanulate": "bell_campanulate",
            "funnelform": "funnel_trumpet",
            "urceolate": "urn_urceolate",
            "bilabiate": "other_described",
        }
        states = {aliases.get(token, token) for token in folded.split("|") if token}
        if states and states.issubset(canonical):
            return next(iter(states)) if len(states) == 1 else "other_described"
        return ""
    if trait == "floral_symmetry":
        aliases = {
            "actinomorphic": "actinomorphic",
            "radial": "actinomorphic",
            "zygomorphic": "zygomorphic",
            "bilateral": "zygomorphic",
            "asymmetric": "asymmetric",
        }
        states = {aliases[token] for token in folded.split("|") if token in aliases}
        return next(iter(states)) if len(states) == 1 else ""
    if trait == "self_incompatibility":
        return {
            "si": "SI",
            "self_incompatible": "SI",
            "sc": "SC",
            "self_compatible": "SC",
            "mixed_or_variable": "mixed_or_variable",
        }.get(folded, "")
    if trait == "mating_system":
        return {
            "predominantly_outcrossing": "predominantly_outcrossing",
            "obligate_outcrossing": "predominantly_outcrossing",
            "mixed_mating": "mixed_mating",
            "predominantly_selfing": "predominantly_selfing",
            "mainly_selfing": "predominantly_selfing",
            "obligate_selfing": "obligate_selfing",
        }.get(folded, "")
    if trait == "sex_system":
        canonical = {
            "hermaphroditic",
            "monoecious",
            "dioecious",
            "gynodioecious",
            "androdioecious",
            "polygamous_or_other",
        }
        aliases = {
            "bisexual": "hermaphroditic",
            "andromonoecious": "polygamous_or_other",
        }
        tokens = {token for token in folded.split("|") if token}
        states = {aliases.get(token, token) for token in tokens if token != "unisexual"}
        if not states or not states.issubset(canonical):
            return ""
        return next(iter(states)) if len(states) == 1 else "polygamous_or_other"
    return raw

OUTPUT_COLUMNS = [
    "accepted_species",
    "genus",
    "family",
    "trait_name",
    "reported_value",
    "inferred_value",
    "filled_value",
    "fill_tier",
    "evidence_scope",
    "confidence",
    "support_n",
    "total_support_n",
    "supporting_genera_n",
    "value_distribution",
    "supporting_taxa",
    "supporting_genera",
    "winner_share",
    "inference_basis",
    "unresolved_reason",
    "rejected_inference_diagnostics",
    "analysis_eligible",
    "direct_conflict",
    "direct_conflict_distribution",
]


@app.callback()
def main() -> None:
    """Fill every master species-trait by taxonomic cascade, tagged by resolution."""


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _is_present(value: object) -> bool:
    return _text(value).lower() not in EMPTY_VALUES


def _write_gzip(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})


def load_config(path: Path) -> dict[str, Any]:
    """Load and minimally validate the versioned cascade configuration."""
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    required = {"master_taxa_csv", "species_column", "genus_column", "family_column",
                "target_traits", "evidence_sources", "tier_labels", "angiosperm_scope_config",
                "trait_ontology_path"}
    if not isinstance(config, dict) or not required.issubset(config):
        raise typer.BadParameter(
            "cascade config must contain master_taxa_csv, species/genus/family columns, "
            "target_traits, evidence_sources, tier_labels, angiosperm_scope_config, "
            "and trait_ontology_path"
        )
    for key in (
        "min_genus_support",
        "min_family_support",
        "min_family_supporting_genera",
    ):
        if int(config.get(key, 1)) < 1:
            raise typer.BadParameter(f"{key} must be at least 1")
    for key in ("min_genus_consensus", "min_family_consensus"):
        value = float(config.get(key, 1.0))
        if not 0.0 < value <= 1.0:
            raise typer.BadParameter(f"{key} must be in (0, 1]")
    traits = [str(value) for value in config["target_traits"]]
    if len(traits) != len(set(traits)):
        raise typer.BadParameter("target_traits must be unique and ordered")
    policies = config.get("inference_policies") or {}
    if not isinstance(policies, dict):
        raise typer.BadParameter("inference_policies must be a trait-keyed mapping")
    unknown_policies = set(policies).difference(traits)
    if unknown_policies:
        raise typer.BadParameter(
            f"inference_policies contain traits outside target_traits: {sorted(unknown_policies)}"
        )
    for trait, policy in policies.items():
        if not isinstance(policy, dict):
            raise typer.BadParameter(f"inference policy for {trait} must be a mapping")
        for key in ("min_genus_support", "min_family_support", "min_family_supporting_genera"):
            if key in policy and int(policy[key]) < 1:
                raise typer.BadParameter(f"{trait}.{key} must be at least 1")
        for key in ("min_genus_consensus", "min_family_consensus"):
            if key in policy and not 0.0 < float(policy[key]) <= 1.0:
                raise typer.BadParameter(f"{trait}.{key} must be in (0, 1]")
        if "family_allowed" in policy and not isinstance(policy["family_allowed"], bool):
            raise typer.BadParameter(f"{trait}.family_allowed must be boolean")
    required_tiers = {
        "species_direct",
        "synonym_direct",
        "genus_inference",
        "family_inference",
        "unresolved_no_evidence",
    }
    missing_tiers = required_tiers.difference(config["tier_labels"])
    if missing_tiers:
        raise typer.BadParameter(
            f"cascade tier_labels missing required tiers: {sorted(missing_tiers)}"
        )
    return config


# --- evidence loading: return (accepted_species, trait_name, value, weight) ---

def _evidence_candidate_index(source: dict[str, Any]) -> pd.DataFrame:
    path = Path(source["path"])
    if not path.exists():
        return pd.DataFrame(columns=["accepted_species", "trait_name", "value", "weight"])
    frame = pd.read_csv(path, dtype=str).fillna("")
    required = {
        "accepted_species",
        "trait_name",
        "candidate_value",
        "evidence_scope",
        "source_url",
        "source_excerpt",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise typer.BadParameter(f"candidate index {path} lacks direct-evidence columns: {sorted(missing)}")
    frame = frame.loc[
        frame["evidence_scope"].map(_text).eq("species_direct")
        & frame["source_url"].map(_is_present)
        & frame["source_excerpt"].map(_is_present)
    ].copy()
    weight = pd.to_numeric(frame.get("n_wave_observations", 1), errors="coerce").fillna(1.0)
    return pd.DataFrame(
        {
            "accepted_species": frame["accepted_species"].map(_text),
            "trait_name": frame["trait_name"].map(_text),
            "value": frame["candidate_value"].map(_text),
            "weight": weight.to_numpy(),
        }
    )


def _evidence_candidate_long(source: dict[str, Any]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in sorted(globlib.glob(source["glob"], recursive=True)):
        frame = pd.read_csv(path, dtype=str).fillna("")
        if frame.empty:
            continue
        require_candidate_kind = bool(source.get("require_candidate_kind", True))
        required = {"trait_name", "evidence_scope", "source_url"}
        if require_candidate_kind:
            required.update({"candidate_kind", "evidence_excerpt"})
        missing = required.difference(frame.columns)
        if missing:
            raise typer.BadParameter(
                f"candidate table {path} lacks direct-evidence columns: {sorted(missing)}"
            )
        direct = frame["evidence_scope"].map(_text).eq("species_direct")
        if require_candidate_kind:
            direct &= frame["candidate_kind"].map(_text).eq("source_backed")
            direct &= frame["evidence_excerpt"].map(_is_present)
        direct &= frame["source_url"].map(_is_present)
        frame = frame.loc[direct].copy()
        if frame.empty:
            continue
        species_col = next((c for c in ("accepted_species", "scientific_name", "species") if c in frame), None)
        value_col = next(
            (
                c
                for c in (
                    "candidate_value",
                    "standardized_value",
                    "provisional_candidate_value",
                    "trait_value",
                    "value",
                )
                if c in frame
            ),
            None,
        )
        if species_col is None or value_col is None or "trait_name" not in frame:
            continue
        frames.append(
            pd.DataFrame(
                {
                    "accepted_species": frame[species_col].map(_text),
                    "trait_name": frame["trait_name"].map(_text),
                    "value": frame[value_col].map(_text),
                    "weight": 1.0,
                }
            )
        )
    if not frames:
        return pd.DataFrame(columns=["accepted_species", "trait_name", "value", "weight"])
    return pd.concat(frames, ignore_index=True)


def _evidence_validated_llm_bundle(source: dict[str, Any]) -> pd.DataFrame:
    """Deep-revalidate frozen packets, exact quotes, ontology, and all hashes."""
    from island_v2.all_species_trait_export import llm_extracted_evidence

    frames: list[pd.DataFrame] = []
    for candidate_name in sorted(globlib.glob(source["glob"], recursive=True)):
        candidate_path = Path(candidate_name)
        bundle = candidate_path.parent
        # This shared validator replays the full all-master acceptance contract:
        # packet/manifest hashes, run id, task/trait binding, allowed value,
        # species/source binding, source URL/text hash, and quote substring.
        llm_extracted_evidence([bundle])
        candidate_rows = len(pd.read_csv(candidate_path, dtype=str))
        frame = _evidence_candidate_long(
            {"glob": str(candidate_path), "require_candidate_kind": False}
        )
        if len(frame) != candidate_rows:
            raise typer.BadParameter(f"LLM evidence direct-scope row contract failed: {bundle}")
        frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=["accepted_species", "trait_name", "value", "weight"])
    return pd.concat(frames, ignore_index=True)


def _evidence_efloras_direct(source: dict[str, Any]) -> pd.DataFrame:
    """Read eFloras rows under its source-specific direct-evidence contract."""
    frames: list[pd.DataFrame] = []
    for path in sorted(globlib.glob(source["glob"], recursive=True)):
        frame = pd.read_csv(path, dtype=str).fillna("")
        if frame.empty:
            continue
        required = {
            "accepted_species",
            "trait_name",
            "candidate_value",
            "evidence_scope",
            "evidence_status",
            "source_url",
            "source_excerpt",
            "source_record_id",
            "extraction_method",
        }
        missing = required.difference(frame.columns)
        if missing:
            raise typer.BadParameter(f"eFloras candidate table {path} lacks: {sorted(missing)}")
        direct = (
            frame["evidence_scope"].map(_text).eq("species_direct")
            & frame["evidence_status"].map(_text).str.startswith("source_backed_")
            & frame["extraction_method"]
            .map(_text)
            .str.startswith("rule_based_floral_context_")
            & frame["source_url"].map(_is_present)
            & frame["source_excerpt"].map(_is_present)
            & frame["source_record_id"].map(_is_present)
        )
        frame = frame.loc[direct].copy()
        frames.append(
            pd.DataFrame(
                {
                    "accepted_species": frame["accepted_species"].map(_text),
                    "trait_name": frame["trait_name"].map(_text),
                    "value": frame["candidate_value"].map(_text),
                    "weight": 1.0,
                }
            )
        )
    if not frames:
        return pd.DataFrame(columns=["accepted_species", "trait_name", "value", "weight"])
    return pd.concat(frames, ignore_index=True)


def _evidence_allmaster_long(source: dict[str, Any]) -> pd.DataFrame:
    """Read the durable all-master reported-evidence ledger as direct evidence."""
    frames: list[pd.DataFrame] = []
    for path in sorted(globlib.glob(source["glob"], recursive=True)):
        frame = pd.read_csv(path, dtype=str).fillna("")
        required = {
            "species",
            "field",
            "value",
            "source_backed",
            "evidence_type",
            "source_url",
            "source_record_id",
            "source_excerpt",
        }
        if frame.empty:
            continue
        missing = required.difference(frame.columns)
        if missing:
            raise typer.BadParameter(f"all-master evidence {path} lacks columns: {sorted(missing)}")
        source_backed = frame["source_backed"].map(_text).str.casefold().isin(
            {"true", "1", "yes"}
        )
        directly_reported = frame["evidence_type"].map(_text).str.casefold().isin(
            {"field_study", "review", "flora", "horticulture"}
        )
        traceable = (
            frame["source_url"].map(_is_present)
            & frame["source_record_id"].map(_is_present)
            & frame["source_excerpt"].map(_is_present)
        )
        frame = frame.loc[source_backed & directly_reported & traceable].copy()
        source_kinds = {_text(value) for value in source.get("include_source_kinds", [])}
        if source_kinds:
            if "source_kind" not in frame.columns:
                raise typer.BadParameter(f"all-master evidence lacks source_kind: {path}")
            frame = frame.loc[frame["source_kind"].map(_text).isin(source_kinds)].copy()
        mapped = frame.apply(
            lambda row: _map_allmaster_value(row["field"], row["value"]), axis=1
        )
        frame["trait_name"] = mapped.map(lambda item: item[0])
        frame["mapped_value"] = mapped.map(lambda item: item[1])
        frame = frame.loc[
            frame["trait_name"].map(_is_present) & frame["mapped_value"].map(_is_present)
        ].copy()
        if frame.empty:
            continue
        frames.append(
            pd.DataFrame(
                {
                    "accepted_species": frame["species"].map(_text),
                    "trait_name": frame["trait_name"].map(_text),
                    "value": frame["mapped_value"].map(_text),
                    "weight": 1.0,
                }
            )
        )
    if not frames:
        return pd.DataFrame(columns=["accepted_species", "trait_name", "value", "weight"])
    return pd.concat(frames, ignore_index=True)


def _evidence_curated(source: dict[str, Any]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in sorted(globlib.glob(source["glob"], recursive=True)):
        frame = pd.read_csv(path, dtype=str).fillna("")
        if frame.empty:
            continue
        value_col = next((c for c in ("trait_value", "candidate_value", "value") if c in frame), None)
        required = {
            "accepted_species",
            "trait_name",
            "review_status",
            "evidence_scope",
            "source_url",
            "evidence_excerpt",
        }
        missing = required.difference(frame.columns)
        if value_col is None or missing:
            raise typer.BadParameter(
                f"curated evidence {path} lacks accepted-direct fields: {sorted(missing)}"
            )
        frame = frame.loc[
            frame["review_status"].map(_text).str.casefold().eq("accepted")
            & frame["evidence_scope"].map(_text).eq("species_direct")
            & frame["source_url"].map(_is_present)
            & frame["evidence_excerpt"].map(_is_present)
        ].copy()
        frames.append(
            pd.DataFrame(
                {
                    "accepted_species": frame["accepted_species"].map(_text),
                    "trait_name": frame["trait_name"].map(_text),
                    "value": frame[value_col].map(_text),
                    # curated evidence outweighs machine candidates in the modal vote.
                    "weight": 1000.0,
                }
            )
        )
    if not frames:
        return pd.DataFrame(columns=["accepted_species", "trait_name", "value", "weight"])
    return pd.concat(frames, ignore_index=True)


def _evidence_search_enabled_llm(source: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for path in sorted(globlib.glob(source["glob"], recursive=True)):
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                parsed = parse_llm_result(_text(json.loads(line).get("result")))
            except (ValueError, KeyError, json.JSONDecodeError):
                continue
            species = _text(parsed.get("species", ""))
            for llm_column, trait_name in LLM_TRAIT_MAP.items():
                value = _text(parsed.get(llm_column, ""))
                if species and _is_present(value):
                    rows.append(
                        {"accepted_species": species, "trait_name": trait_name, "value": value, "weight": 1.0}
                    )
    return pd.DataFrame(rows, columns=["accepted_species", "trait_name", "value", "weight"])


EVIDENCE_ADAPTERS = {
    "candidate_index": _evidence_candidate_index,
    "candidate_long": _evidence_candidate_long,
    "efloras_direct": _evidence_efloras_direct,
    "validated_llm_bundle": _evidence_validated_llm_bundle,
    "allmaster_long": _evidence_allmaster_long,
    "curated_evidence": _evidence_curated,
    "search_enabled_llm": _evidence_search_enabled_llm,
}


def load_direct_evidence(config: dict[str, Any], traits: set[str]) -> pd.DataFrame:
    """Union all evidence channels into (species, trait, value, weight) rows."""
    frames: list[pd.DataFrame] = []
    for key, source in config["evidence_sources"].items():
        adapter = source.get("adapter")
        if adapter not in EVIDENCE_ADAPTERS:
            raise typer.BadParameter(f"evidence source {key!r} has unknown adapter {adapter!r}")
        frame = EVIDENCE_ADAPTERS[adapter](source)
        if not frame.empty:
            frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=["accepted_species", "trait_name", "value", "weight"])
    evidence = pd.concat(frames, ignore_index=True)
    evidence = evidence.loc[
        evidence["trait_name"].isin(traits)
        & evidence["accepted_species"].map(_is_present)
        & evidence["value"].map(_is_present)
    ].reset_index(drop=True)
    normalized = [
        _normalize_direct_value(trait, value)
        for trait, value in zip(evidence["trait_name"], evidence["value"], strict=True)
    ]
    evidence["value"] = normalized
    unmapped = evidence["value"].map(_text).eq("")
    unmapped_by_trait = {
        str(key): int(value)
        for key, value in evidence.loc[unmapped, "trait_name"].value_counts().items()
    }

    ontology_path = Path(config["trait_ontology_path"])
    if not ontology_path.is_file():
        raise typer.BadParameter(f"trait ontology does not exist: {ontology_path}")
    ontology = yaml.safe_load(ontology_path.read_text(encoding="utf-8"))
    if not isinstance(ontology, dict) or not isinstance(ontology.get("traits"), dict):
        raise typer.BadParameter(f"invalid trait ontology: {ontology_path}")
    allowed_by_trait = {
        _text(trait): set(details.get("allowed_values", [])) - {"unresolved"}
        for trait, details in ontology["traits"].items()
    }
    missing_ontology = traits.difference(allowed_by_trait)
    if missing_ontology:
        raise typer.BadParameter(
            f"target traits missing from ontology: {sorted(missing_ontology)}"
        )
    ontology_valid = pd.Series(
        [
            not value or value in allowed_by_trait.get(trait, set())
            for trait, value in zip(evidence["trait_name"], evidence["value"], strict=True)
        ],
        index=evidence.index,
    )
    invalid = evidence.loc[~ontology_valid, ["accepted_species", "trait_name", "value"]]
    if not invalid.empty:
        sample = invalid.head(10).to_dict("records")
        raise typer.BadParameter(
            f"direct evidence contains {len(invalid)} out-of-ontology values; sample={sample}"
        )
    evidence = evidence.loc[
        evidence["value"].map(_is_present)
    ].reset_index(drop=True)
    evidence.attrs["direct_evidence_audit"] = {
        "n_normalization_unmapped": int(unmapped.sum()),
        "normalization_unmapped_by_trait": unmapped_by_trait,
        "n_ontology_invalid": 0,
        "ontology_path": str(ontology_path),
    }
    return evidence


def _modal(frame: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    """Return modal candidates with auditable support and consensus."""
    if frame.empty:
        return pd.DataFrame(
            columns=[
                *group_cols,
                "value",
                "weight",
                "value_distribution",
                "support_n",
                "total_support_n",
                "supporting_genera_n",
                "supporting_taxa",
                "supporting_genera",
                "total_weight",
                "winner_share",
                "top_tie_n",
            ]
        )
    totals = frame.groupby([*group_cols, "value"], as_index=False)["weight"].sum()
    totals = totals.sort_values(
        [*group_cols, "weight", "value"],
        ascending=[*([True] * len(group_cols)), False, True],
    )
    winners = totals.drop_duplicates(group_cols, keep="first")
    total_weight = totals.groupby(group_cols)["weight"].sum().rename("total_weight").reset_index()
    winners = winners.merge(total_weight, on=group_cols)
    winners["winner_share"] = winners["weight"] / winners["total_weight"]
    top_weight = totals.groupby(group_cols)["weight"].transform("max")
    ties = (
        totals.loc[totals["weight"].eq(top_weight)]
        .groupby(group_cols)
        .size()
        .rename("top_tie_n")
        .reset_index()
    )
    dist = (
        totals.groupby(group_cols)
        .apply(lambda g: json.dumps({str(v): round(float(w), 3) for v, w in zip(g["value"], g["weight"], strict=True)}), include_groups=False)
        .rename("value_distribution")
        .reset_index()
    )
    if "accepted_species" in frame.columns:
        support_frame = frame.assign(_support_taxon=frame["accepted_species"].map(_text))
        winner_keys = winners[[*group_cols, "value"]]
        winner_support_frame = support_frame.merge(
            winner_keys,
            on=[*group_cols, "value"],
            how="inner",
            validate="many_to_one",
        )
        support = (
            winner_support_frame.groupby(group_cols)["_support_taxon"]
            .nunique()
            .rename("support_n")
            .reset_index()
        )
        total_support = (
            support_frame.groupby(group_cols)["_support_taxon"]
            .nunique()
            .rename("total_support_n")
            .reset_index()
        )
        supporting_taxa = (
            winner_support_frame.groupby(group_cols)["_support_taxon"]
            .apply(
                lambda values: json.dumps(
                    sorted({value for value in values if value}), ensure_ascii=False
                )
            )
            .rename("supporting_taxa")
            .reset_index()
        )
        if "genus" in winner_support_frame.columns:
            winner_support_frame = winner_support_frame.assign(
                _support_genus=winner_support_frame["genus"].map(_text)
            )
            nonblank_genera = winner_support_frame.loc[
                winner_support_frame["_support_genus"].ne("")
            ]
            supporting_genera_n = (
                nonblank_genera.groupby(group_cols)["_support_genus"]
                .nunique()
                .rename("supporting_genera_n")
                .reindex(
                    pd.MultiIndex.from_frame(support[group_cols]),
                    fill_value=0,
                )
                .rename_axis(group_cols)
                .reset_index()
            )
            supporting_genera = (
                nonblank_genera.groupby(group_cols)["_support_genus"]
                .apply(
                    lambda values: json.dumps(
                        sorted(set(values)),
                        ensure_ascii=False,
                    )
                )
                .reindex(
                    pd.MultiIndex.from_frame(support[group_cols]),
                    fill_value="[]",
                )
                .rename("supporting_genera")
                .rename_axis(group_cols)
                .reset_index()
            )
        else:
            supporting_genera_n = support[group_cols].copy()
            supporting_genera_n["supporting_genera_n"] = 0
            supporting_genera = support[group_cols].copy()
            supporting_genera["supporting_genera"] = "[]"
    else:
        support = frame.groupby(group_cols)["value"].size().rename("support_n").reset_index()
        total_support = support.rename(columns={"support_n": "total_support_n"})
        supporting_taxa = support[group_cols].copy()
        supporting_taxa["supporting_taxa"] = "[]"
        supporting_genera_n = support[group_cols].copy()
        supporting_genera_n["supporting_genera_n"] = 0
        supporting_genera = support[group_cols].copy()
        supporting_genera["supporting_genera"] = "[]"
    return (
        winners.merge(dist, on=group_cols)
        .merge(support, on=group_cols)
        .merge(total_support, on=group_cols)
        .merge(supporting_genera_n, on=group_cols)
        .merge(supporting_taxa, on=group_cols)
        .merge(supporting_genera, on=group_cols)
        .merge(ties, on=group_cols)
    )


def _qualified_inference(
    candidate: Any,
    *,
    min_support: int,
    min_consensus: float,
    min_supporting_genera: int = 0,
) -> bool:
    """Require a unique winner plus explicit taxon support and consensus."""
    return bool(
        candidate is not None
        and int(candidate.support_n) >= min_support
        and int(candidate.top_tie_n) == 1
        and float(candidate.winner_share) >= min_consensus
        and int(candidate.supporting_genera_n) >= min_supporting_genera
    )


def _rejected_inference_diagnostic(
    candidate: Any,
    *,
    basis: str,
    min_support: int,
    min_consensus: float,
    min_supporting_genera: int = 0,
) -> dict[str, Any] | None:
    """Describe why an available taxonomic candidate did not qualify."""
    if candidate is None:
        return None
    failed_gates: list[str] = []
    if int(candidate.support_n) < min_support:
        failed_gates.append("winner_support_below_threshold")
    if int(candidate.top_tie_n) != 1:
        failed_gates.append("top_value_tied")
    if float(candidate.winner_share) < min_consensus:
        failed_gates.append("consensus_below_threshold")
    if int(candidate.supporting_genera_n) < min_supporting_genera:
        failed_gates.append("supporting_genera_below_threshold")
    return {
        "basis": basis,
        "support_n": int(candidate.support_n),
        "total_support_n": int(candidate.total_support_n),
        "supporting_genera_n": int(candidate.supporting_genera_n),
        "winner_share": round(float(candidate.winner_share), 6),
        "top_tie_n": int(candidate.top_tie_n),
        "value_distribution": json.loads(str(candidate.value_distribution)),
        "supporting_taxa": json.loads(str(candidate.supporting_taxa)),
        "supporting_genera": json.loads(str(candidate.supporting_genera)),
        "failed_gates": failed_gates,
    }


def build_model(master: pd.DataFrame, evidence: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    """Build the shared cascade model once from all evidence, for per-shard reuse.

    Genus/family distributions are computed over the full eligible master
    and evidence, so applying the model to any species subset yields identical
    fills to a single whole-universe pass.
    """
    # Multiple directly reported colours are a valid plural signal; retain them
    # as multicoloured/variable rather than choosing an arbitrary colour.
    species_evidence = evidence.copy()
    color = species_evidence["trait_name"].eq("flower_primary_color")
    multi_color = (
        species_evidence.loc[color]
        .groupby(["accepted_species", "trait_name"])["value"]
        .nunique()
        .loc[lambda values: values.gt(1)]
    )
    multi_color_keys = set(multi_color.index)
    if multi_color_keys:
        keys = list(zip(species_evidence["accepted_species"], species_evidence["trait_name"], strict=True))
        species_evidence.loc[
            [key in multi_color_keys for key in keys], "value"
        ] = "multicolored_variable"

    # Equal-weight exclusive states are unresolved direct conflicts. They stay
    # visible in the output but do not become reported values or train priors.
    species_candidates = _modal(species_evidence, ["accepted_species", "trait_name"])
    species_conflicts = species_candidates.loc[species_candidates["top_tie_n"].gt(1)].copy()
    species_direct = species_candidates.loc[species_candidates["top_tie_n"].eq(1)].copy()
    species_value = {
        (r.accepted_species, r.trait_name): r for r in species_direct.itertuples(index=False)
    }
    conflict_lookup = {
        (r.accepted_species, r.trait_name): r
        for r in species_conflicts.itertuples(index=False)
    }

    # Attach each species-direct value to its genus/family for inference tiers.
    taxo = master[["accepted_species", "genus", "family"]].drop_duplicates("accepted_species")
    direct_taxo = species_direct.merge(taxo, on="accepted_species", how="left")
    direct_taxo["weight"] = 1.0

    genus_modal = _modal(
        direct_taxo[["accepted_species", "genus", "trait_name", "value", "weight"]],
        ["genus", "trait_name"],
    )
    genus_lookup = {
        (r.genus, r.trait_name): r for r in genus_modal.itertuples(index=False)
    }
    family_modal = _modal(
        direct_taxo[
            ["accepted_species", "genus", "family", "trait_name", "value", "weight"]
        ],
        ["family", "trait_name"],
    )
    family_lookup = {
        (r.family, r.trait_name): r for r in family_modal.itertuples(index=False)
    }
    return {
        "species_value": species_value,
        "direct_conflicts": conflict_lookup,
        "genus_lookup": genus_lookup,
        "family_lookup": family_lookup,
    }


def fill_species_frame(
    species_frame: pd.DataFrame, model: dict[str, Any], config: dict[str, Any]
) -> pd.DataFrame:
    """Resolve one species subset without inventing a global biological value."""
    traits = list(config["target_traits"])
    labels = config["tier_labels"]
    inference_policies = config.get("inference_policies") or {}
    species_value = model["species_value"]
    direct_conflicts = model["direct_conflicts"]
    genus_lookup = model["genus_lookup"]
    family_lookup = model["family_lookup"]

    rows: list[dict[str, Any]] = []
    for sp in species_frame.itertuples(index=False):
        species = _text(sp.accepted_species)
        genus = _text(sp.genus)
        family = _text(sp.family)
        for trait in traits:
            policy = inference_policies.get(trait) or {}
            min_genus = int(policy.get("min_genus_support", config.get("min_genus_support", 1)))
            min_family = int(
                policy.get("min_family_support", config.get("min_family_support", 3))
            )
            min_family_genera = int(
                policy.get(
                    "min_family_supporting_genera",
                    config.get("min_family_supporting_genera", 2),
                )
            )
            min_genus_consensus = float(
                policy.get(
                    "min_genus_consensus",
                    config.get("min_genus_consensus", 1.0),
                )
            )
            min_family_consensus = float(
                policy.get(
                    "min_family_consensus",
                    config.get("min_family_consensus", 0.8),
                )
            )
            family_allowed = bool(policy.get("family_allowed", True))
            direct = species_value.get((species, trait))
            conflict = direct_conflicts.get((species, trait))
            if direct is not None:
                tier, value = "species_direct", direct.value
                candidate = direct
                inference_basis = "species_direct"
                unresolved_reason = ""
                rejected_diagnostics = ""
            else:
                gen = genus_lookup.get((genus, trait)) if genus else None
                fam = family_lookup.get((family, trait)) if family else None
                if _qualified_inference(
                    gen,
                    min_support=min_genus,
                    min_consensus=min_genus_consensus,
                ):
                    tier, value, candidate = "genus_inference", str(gen.value), gen
                    inference_basis = f"genus:{genus}"
                    unresolved_reason = ""
                    rejected_diagnostics = ""
                elif family_allowed and _qualified_inference(
                    fam,
                    min_support=min_family,
                    min_consensus=min_family_consensus,
                    min_supporting_genera=min_family_genera,
                ):
                    tier, value, candidate = "family_inference", str(fam.value), fam
                    inference_basis = f"family:{family}"
                    unresolved_reason = ""
                    rejected_diagnostics = ""
                else:
                    tier, value = "unresolved_no_evidence", "unresolved"
                    candidate = None
                    inference_basis = ""
                    diagnostics = [
                        diagnostic
                        for diagnostic in (
                            _rejected_inference_diagnostic(
                                gen,
                                basis=f"genus:{genus}",
                                min_support=min_genus,
                                min_consensus=min_genus_consensus,
                            ),
                            _rejected_inference_diagnostic(
                                fam,
                                basis=f"family:{family}",
                                min_support=min_family,
                                min_consensus=min_family_consensus,
                                min_supporting_genera=min_family_genera,
                            ),
                        )
                        if diagnostic is not None
                    ]
                    if fam is not None and not family_allowed:
                        family_diagnostic = _rejected_inference_diagnostic(
                            fam,
                            basis=f"family:{family}",
                            min_support=min_family,
                            min_consensus=min_family_consensus,
                            min_supporting_genera=min_family_genera,
                        )
                        if family_diagnostic is not None:
                            family_diagnostic["failed_gates"] = sorted(
                                set(family_diagnostic["failed_gates"])
                                | {"family_inference_disabled_for_trait"}
                            )
                            diagnostics[-1:] = [family_diagnostic]
                    rejected_diagnostics = json.dumps(
                        diagnostics, ensure_ascii=False, sort_keys=True
                    )
                    if conflict is not None:
                        unresolved_reason = (
                            "species_direct_conflict_no_qualified_taxonomic_inference"
                        )
                    elif not diagnostics:
                        unresolved_reason = "no_genus_or_family_direct_evidence"
                    elif not family_allowed:
                        unresolved_reason = "no_qualified_genus_inference_family_disabled"
                    else:
                        unresolved_reason = "no_qualified_genus_or_family_inference"
            label = labels[tier]
            is_reported = tier in {"species_direct", "synonym_direct"}
            is_inferred = tier in {"genus_inference", "family_inference"}
            support = int(candidate.support_n) if candidate is not None else 0
            total_support = int(candidate.total_support_n) if candidate is not None else 0
            supporting_genera_n = (
                int(candidate.supporting_genera_n) if candidate is not None else 0
            )
            dist = str(candidate.value_distribution) if candidate is not None else ""
            supporting_taxa = str(candidate.supporting_taxa) if candidate is not None else "[]"
            supporting_genera = (
                str(candidate.supporting_genera) if candidate is not None else "[]"
            )
            winner_share = (
                round(float(candidate.winner_share), 6) if candidate is not None else ""
            )
            rows.append(
                {
                    "accepted_species": species,
                    "genus": genus,
                    "family": family,
                    "trait_name": trait,
                    "reported_value": value if is_reported else "",
                    "inferred_value": value if is_inferred else "",
                    "filled_value": value,
                    "fill_tier": tier,
                    "evidence_scope": label["evidence_scope"],
                    "confidence": label["confidence"],
                    "support_n": support,
                    "total_support_n": total_support,
                    "supporting_genera_n": supporting_genera_n,
                    "value_distribution": dist if not is_reported else "",
                    "supporting_taxa": supporting_taxa,
                    "supporting_genera": supporting_genera,
                    "winner_share": winner_share,
                    "inference_basis": inference_basis,
                    "unresolved_reason": unresolved_reason,
                    "rejected_inference_diagnostics": rejected_diagnostics,
                    "analysis_eligible": "false" if tier == "unresolved_no_evidence" else "true",
                    "direct_conflict": "true" if conflict is not None else "false",
                    "direct_conflict_distribution": (
                        conflict.value_distribution if conflict is not None else ""
                    ),
                }
            )
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


def build_fills(master: pd.DataFrame, evidence: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Fill every master species-trait by descending the taxonomic cascade."""
    model = build_model(master, evidence, config)
    return fill_species_frame(master, model, config)


def build_coverage_summary(fills: pd.DataFrame, config: dict[str, Any], n_master: int) -> dict[str, Any]:
    """Per-trait grounded resolution and rectangular cell-presence coverage."""
    traits = list(config["target_traits"])
    by_trait: dict[str, Any] = {}
    for trait in traits:
        sub = fills.loc[fills["trait_name"].eq(trait)]
        n_cells_present = int(sub["accepted_species"].nunique())
        n_unresolved = int(
            sub.loc[sub["fill_tier"].eq("unresolved_no_evidence"), "accepted_species"].nunique()
        )
        n_filled = n_cells_present - n_unresolved
        tiers = {str(k): int(v) for k, v in sub["fill_tier"].value_counts().items()}
        direct = int(sub.loc[sub["fill_tier"].eq("species_direct"), "accepted_species"].nunique())
        conflicts = int(
            sub.loc[sub["direct_conflict"].eq("true"), "accepted_species"].nunique()
        )
        by_trait[trait] = {
            "n_filled": n_filled,
            "fill_rate": n_filled / n_master if n_master else 0.0,
            "n_cells_present": n_cells_present,
            "cell_presence_rate": n_cells_present / n_master if n_master else 0.0,
            "n_species_direct": direct,
            "n_species_direct_conflict": conflicts,
            "species_direct_rate": direct / n_master if n_master else 0.0,
            "n_unresolved_no_evidence": n_unresolved,
            "n_unknown_remaining": n_master - n_filled,
            "by_tier": tiers,
        }
    overall_tiers = {str(k): int(v) for k, v in fills["fill_tier"].value_counts().items()}
    return {
        "version": "1.0",
        "n_master_species": n_master,
        "n_target_traits": len(traits),
        "by_trait": by_trait,
        "fills_by_tier": overall_tiers,
        "interpretation": (
            "Evidence-grounded cascade coverage. Genus/family values require explicit "
            "support and consensus; cells without qualified taxonomic evidence remain "
            "unresolved_no_evidence and are not analysis-eligible."
        ),
    }


def build_benchmark(fills: pd.DataFrame, master: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Wide per-field fill for a deterministic species sample (the 100-sp benchmark)."""
    size = int(config.get("benchmark_sample_size", 100))
    species = sorted(master["accepted_species"].map(_text).loc[lambda s: s.ne("")].unique())[:size]
    sub = fills.loc[fills["accepted_species"].isin(species)].copy()
    sub["cell"] = sub["filled_value"] + " [" + sub["fill_tier"] + "]"
    wide = sub.pivot_table(
        index="accepted_species", columns="trait_name", values="cell", aggfunc="first"
    ).reindex(species)
    return wide.reset_index()


# --- sharding, checkpoint, resume ------------------------------------------

CONTRACT_VERSION = "trait_fill_cascade_shards_v4"


def _sha_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def stable_shard(species: str, shard_count: int) -> int:
    """Return a stable shard for a species, independent of input row order."""
    if shard_count < 1:
        raise ValueError("shard_count must be positive")
    return int(_sha_text(_text(species))[:16], 16) % shard_count


def _now() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _atomic_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def _atomic_write_gzip(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    _write_gzip(frame, tmp)
    os.replace(tmp, path)


def model_fingerprint(
    evidence: pd.DataFrame, config: dict[str, Any], master: pd.DataFrame | None = None
) -> str:
    """Fingerprint the shared model: any evidence or cascade-parameter change flips it."""
    ordered = (
        evidence[["accepted_species", "trait_name", "value", "weight"]]
        .astype(str)
        .sort_values(["accepted_species", "trait_name", "value", "weight"])
        if not evidence.empty
        else evidence
    )
    payload = "\n".join(
        "\t".join(row) for row in ordered.itertuples(index=False, name=None)
    ) if not evidence.empty else ""
    taxonomy_payload = ""
    if master is not None and not master.empty:
        taxonomy = (
            master[["accepted_species", "genus", "family"]]
            .astype(str)
            .sort_values(["accepted_species", "genus", "family"])
        )
        taxonomy_payload = "\n".join(
            "\t".join(row) for row in taxonomy.itertuples(index=False, name=None)
        )
    params = json.dumps(
        {
            "traits": list(config["target_traits"]),
            "min_genus_support": int(config.get("min_genus_support", 1)),
            "min_family_support": int(config.get("min_family_support", 3)),
            "min_family_supporting_genera": int(
                config.get("min_family_supporting_genera", 2)
            ),
            "min_genus_consensus": float(config.get("min_genus_consensus", 1.0)),
            "min_family_consensus": float(config.get("min_family_consensus", 0.8)),
            "inference_policies": dict(config.get("inference_policies") or {}),
            "tier_labels": config["tier_labels"],
            "contract": CONTRACT_VERSION,
        },
        sort_keys=True,
    )
    return _sha_text(params + "\n" + payload + "\n" + taxonomy_payload)


def shard_species_fingerprint(species_frame: pd.DataFrame) -> str:
    """Fingerprint exact shard species membership and its taxonomy."""
    taxonomy = (
        species_frame[["accepted_species", "genus", "family"]]
        .astype(str)
        .sort_values(["accepted_species", "genus", "family"])
    )
    return _sha_text(
        "\n".join("\t".join(row) for row in taxonomy.itertuples(index=False, name=None))
    )


def shard_summary_from_fills(shard_fills: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    """Compact per-shard trait counts, disjoint across shards so global sums are exact."""
    traits = list(config["target_traits"])
    by_trait: dict[str, Any] = {}
    for trait in traits:
        sub = shard_fills.loc[shard_fills["trait_name"].eq(trait)]
        n_cells_present = int(sub["accepted_species"].nunique())
        n_unresolved = int(
            sub.loc[sub["fill_tier"].eq("unresolved_no_evidence"), "accepted_species"].nunique()
        )
        by_trait[trait] = {
            "n_filled": n_cells_present - n_unresolved,
            "n_cells_present": n_cells_present,
            "n_unresolved_no_evidence": n_unresolved,
            "n_species_direct": int(
                sub.loc[sub["fill_tier"].eq("species_direct"), "accepted_species"].nunique()
            ),
            "n_species_direct_conflict": int(
                sub.loc[sub["direct_conflict"].eq("true"), "accepted_species"].nunique()
            ),
            "by_tier": {str(k): int(v) for k, v in sub["fill_tier"].value_counts().items()},
        }
    return {
        "n_species": int(shard_fills["accepted_species"].nunique()) if not shard_fills.empty else 0,
        "by_trait": by_trait,
        "fills_by_tier": {str(k): int(v) for k, v in shard_fills["fill_tier"].value_counts().items()},
    }


def aggregate_shard_summaries(
    summaries: list[dict[str, Any]], config: dict[str, Any], n_eligible: int
) -> dict[str, Any]:
    """Sum disjoint per-shard summaries into the global coverage summary."""
    traits = list(config["target_traits"])
    by_trait: dict[str, Any] = {}
    fills_by_tier: dict[str, int] = {}
    for trait in traits:
        n_filled = sum(s["by_trait"].get(trait, {}).get("n_filled", 0) for s in summaries)
        n_cells_present = sum(
            s["by_trait"].get(trait, {}).get("n_cells_present", 0) for s in summaries
        )
        n_unresolved = sum(
            s["by_trait"].get(trait, {}).get("n_unresolved_no_evidence", 0)
            for s in summaries
        )
        direct = sum(s["by_trait"].get(trait, {}).get("n_species_direct", 0) for s in summaries)
        conflicts = sum(
            s["by_trait"].get(trait, {}).get("n_species_direct_conflict", 0)
            for s in summaries
        )
        tiers: dict[str, int] = {}
        for s in summaries:
            for tier, count in s["by_trait"].get(trait, {}).get("by_tier", {}).items():
                tiers[tier] = tiers.get(tier, 0) + count
        by_trait[trait] = {
            "n_filled": n_filled,
            "fill_rate": n_filled / n_eligible if n_eligible else 0.0,
            "n_cells_present": n_cells_present,
            "cell_presence_rate": n_cells_present / n_eligible if n_eligible else 0.0,
            "n_species_direct": direct,
            "n_species_direct_conflict": conflicts,
            "species_direct_rate": direct / n_eligible if n_eligible else 0.0,
            "n_unresolved_no_evidence": n_unresolved,
            "n_unknown_remaining": n_eligible - n_filled,
            "by_tier": tiers,
        }
    for s in summaries:
        for tier, count in s.get("fills_by_tier", {}).items():
            fills_by_tier[tier] = fills_by_tier.get(tier, 0) + count
    return {
        "version": "1.0",
        "n_master_species": n_eligible,
        "n_target_traits": len(traits),
        "by_trait": by_trait,
        "fills_by_tier": fills_by_tier,
        "interpretation": (
            "Evidence-grounded cascade coverage. Genus/family values require explicit "
            "support and consensus; cells without qualified taxonomic evidence remain "
            "unresolved_no_evidence and are not analysis-eligible."
        ),
    }


def load_scoped_master_and_evidence(
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Load master, apply the angiosperm gate, and load eligible-only evidence."""
    master = pd.read_csv(config["master_taxa_csv"], dtype=str).fillna("").rename(
        columns={
            config["species_column"]: "accepted_species",
            config["genus_column"]: "genus",
            config["family_column"]: "family",
        }
    )
    master = master.loc[master["accepted_species"].map(_is_present)].drop_duplicates("accepted_species")

    scope_config = load_scope_config(Path(config["angiosperm_scope_config"]))
    scope = classify_scope(master, scope_config)
    eligible_species = set(scope.loc[scope["angiosperm_analysis_eligible"], "accepted_species"])
    eligible_master = (
        master.loc[master["accepted_species"].isin(eligible_species)]
        .sort_values("accepted_species")
        .reset_index(drop=True)
    )

    traits = set(config["target_traits"])
    evidence = load_direct_evidence(config, traits)
    direct_evidence_audit = evidence.attrs.get("direct_evidence_audit", {})
    evidence = evidence.loc[evidence["accepted_species"].isin(eligible_species)].reset_index(drop=True)
    evidence.attrs["direct_evidence_audit"] = direct_evidence_audit
    return master, eligible_master, evidence, {
        "scope": scope,
        "direct_evidence_audit": direct_evidence_audit,
    }


@app.command("run")
def run(
    config_path: Path = typer.Option(Path("config/trait_fill_cascade.yml"), exists=True),
    output_dir: Path = typer.Option(..., help="Directory for cascade fill outputs."),
    shard_count: int = typer.Option(128, min=1, help="Number of deterministic species shards."),
    max_shards: int = typer.Option(
        0, min=0, help="Max stale shards to process this run (0 = all remaining)."
    ),
    no_resume: bool = typer.Option(False, help="Recompute every shard, ignoring checkpoints."),
) -> None:
    """Fill eligible species over resumable shards, then finalize when all are current.

    Each shard is checkpointed with the model and species fingerprints. A shard is
    recomputed only when its inputs changed or it is missing, so repeat runs do not
    recompute the whole universe. --max-shards bounds work per invocation for
    scheduled lanes; rerun to resume and finalize.
    """
    config = load_config(config_path)
    master, eligible_master, evidence, extra = load_scoped_master_and_evidence(config)
    scope = extra["scope"]
    n_master = int(len(master))
    n_eligible = int(len(eligible_master))

    model = build_model(eligible_master, evidence, config)
    fingerprint = model_fingerprint(evidence, config, eligible_master)

    assignments = eligible_master["accepted_species"].map(lambda s: stable_shard(s, shard_count))
    shards_dir = output_dir / "shards"
    shards_dir.mkdir(parents=True, exist_ok=True)

    processed, skipped, stale = 0, 0, 0
    for shard_index in range(shard_count):
        shard_species = eligible_master.loc[assignments.eq(shard_index)].reset_index(drop=True)
        species_fp = shard_species_fingerprint(shard_species)
        shard_dir = shards_dir / f"shard_{shard_index:05d}"
        checkpoint_path = shard_dir / "shard_checkpoint.json"

        current = None
        if not no_resume and checkpoint_path.exists():
            try:
                current = json.loads(checkpoint_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                current = None
        up_to_date = (
            current is not None
            and current.get("model_fingerprint") == fingerprint
            and current.get("species_fingerprint") == species_fp
            and (shard_dir / "fills.csv.gz").exists()
        )
        if up_to_date:
            skipped += 1
            continue

        stale += 1
        if max_shards and processed >= max_shards:
            continue

        shard_fills = fill_species_frame(shard_species, model, config)
        summary = shard_summary_from_fills(shard_fills, config)
        _atomic_write_gzip(shard_fills, shard_dir / "fills.csv.gz")
        _atomic_json(
            {
                "contract_version": CONTRACT_VERSION,
                "shard_index": shard_index,
                "shard_count": shard_count,
                "model_fingerprint": fingerprint,
                "species_fingerprint": species_fp,
                "summary": summary,
                "updated_at": _now(),
            },
            checkpoint_path,
        )
        processed += 1

    remaining = stale - processed
    complete = remaining == 0
    _finalize(
        output_dir,
        shards_dir,
        shard_count,
        config,
        master,
        eligible_master,
        scope,
        fingerprint,
        n_master,
        n_eligible,
        complete,
        extra.get("direct_evidence_audit", {}),
    )

    typer.echo(
        f"Shards {shard_count}: processed {processed}, skipped up-to-date {skipped}, "
        f"remaining stale {remaining}. Eligible {n_eligible}/{n_master}. "
        f"{'Finalized (all shards current).' if complete else 'Partial; rerun to resume.'}"
    )


def _finalize(
    output_dir: Path,
    shards_dir: Path,
    shard_count: int,
    config: dict[str, Any],
    master: pd.DataFrame,
    eligible_master: pd.DataFrame,
    scope: pd.DataFrame,
    fingerprint: str,
    n_master: int,
    n_eligible: int,
    complete: bool,
    direct_evidence_audit: dict[str, Any],
) -> None:
    """Aggregate current shard summaries into the coverage summary and benchmark."""
    summaries: list[dict[str, Any]] = []
    shards_done = 0
    for shard_index in range(shard_count):
        checkpoint_path = shards_dir / f"shard_{shard_index:05d}" / "shard_checkpoint.json"
        if not checkpoint_path.exists():
            continue
        payload = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        if payload.get("model_fingerprint") != fingerprint:
            continue
        summaries.append(payload["summary"])
        shards_done += 1

    coverage = aggregate_shard_summaries(summaries, config, n_eligible)
    coverage["n_master_species_all"] = n_master
    coverage["n_angiosperm_eligible"] = n_eligible
    coverage["n_out_of_scope"] = n_master - n_eligible
    coverage["out_of_scope_by_group"] = {
        str(k): int(v)
        for k, v in scope.loc[~scope["angiosperm_analysis_eligible"], "taxonomic_group"].value_counts().items()
    }
    coverage["denominator"] = "angiosperm_eligible"
    coverage["shards_complete"] = complete
    coverage["n_shards_current"] = shards_done
    coverage["shard_count"] = shard_count
    coverage["direct_evidence_audit"] = direct_evidence_audit

    output_dir.mkdir(parents=True, exist_ok=True)
    scope.to_csv(output_dir / "angiosperm_scope_by_species.csv", index=False)
    (output_dir / "fill_coverage_summary.json").write_text(
        json.dumps(coverage, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    _atomic_json(
        {
            "contract_version": CONTRACT_VERSION,
            "shard_count": shard_count,
            "model_fingerprint": fingerprint,
            "n_master_species_all": n_master,
            "n_angiosperm_eligible": n_eligible,
            "n_shards_current": shards_done,
            "shards_complete": complete,
            "direct_evidence_audit": direct_evidence_audit,
            "upstream_allmaster_run_id": os.environ.get(
                "ISLAND_ALLMASTER_SOURCE_RUN_ID", ""
            ),
            "updated_at": _now(),
        },
        output_dir / "cascade_manifest.json",
    )

    # Benchmark reads only the shards holding the deterministic sample species.
    size = int(config.get("benchmark_sample_size", 100))
    sample = sorted(eligible_master["accepted_species"].map(_text).loc[lambda s: s.ne("")].unique())[:size]
    sample_shards = {stable_shard(s, shard_count) for s in sample}
    frames = []
    for shard_index in sorted(sample_shards):
        path = shards_dir / f"shard_{shard_index:05d}" / "fills.csv.gz"
        if path.exists():
            frames.append(pd.read_csv(path, dtype=str).fillna(""))
    if frames:
        benchmark = build_benchmark(pd.concat(frames, ignore_index=True), eligible_master, config)
    else:
        benchmark = pd.DataFrame(columns=["accepted_species"])
    benchmark.to_csv(output_dir / "benchmark_sample.csv", index=False)


@app.command("status")
def status(
    output_dir: Path = typer.Option(..., help="Cascade output directory to inspect."),
    config_path: Path = typer.Option(Path("config/trait_fill_cascade.yml"), exists=True),
    shard_count: int = typer.Option(128, min=1, help="Shard count used for the run."),
) -> None:
    """Report how many shards are current for the present inputs, without filling."""
    config = load_config(config_path)
    _, eligible_master, evidence, _ = load_scoped_master_and_evidence(config)
    fingerprint = model_fingerprint(evidence, config, eligible_master)
    assignments = eligible_master["accepted_species"].map(lambda s: stable_shard(s, shard_count))
    shards_dir = output_dir / "shards"

    current, stale, missing = 0, 0, 0
    for shard_index in range(shard_count):
        shard_species = eligible_master.loc[assignments.eq(shard_index)].reset_index(drop=True)
        species_fp = shard_species_fingerprint(shard_species)
        checkpoint_path = shards_dir / f"shard_{shard_index:05d}" / "shard_checkpoint.json"
        fills_path = shards_dir / f"shard_{shard_index:05d}" / "fills.csv.gz"
        if not checkpoint_path.exists() or not fills_path.exists():
            missing += 1
            continue
        payload = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        if payload.get("model_fingerprint") == fingerprint and payload.get("species_fingerprint") == species_fp:
            current += 1
        else:
            stale += 1
    typer.echo(
        f"Shards {shard_count}: current {current}, stale {stale}, missing {missing}. "
        f"{'All current.' if current == shard_count else 'Rerun run to resume.'}"
    )


if __name__ == "__main__":
    app()
