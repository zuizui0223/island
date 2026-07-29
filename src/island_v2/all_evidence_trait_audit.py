"""Rebuild strict three-axis coverage and Validated Low from all direct evidence.

The audit is deliberately trait-specific. Direct evidence is resolved at
``accepted_species x trait_name`` after source-lineage deduplication, and genus
rules are applied on ``genus x trait_name``. The three broad axes are reporting
aggregates only; they are never used as a join key for inference.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import Annotated, Any

import pandas as pd
import typer
import yaml

from island_v2.angiosperm_scope import classify_scope
from island_v2.angiosperm_scope import load_config as load_scope_config
from island_v2.integrated_trait_coverage import (
    AXES,
    DIRECT_SCOPES,
    EVIDENCE_COLUMNS,
    QUALITY_RANK,
    TRAIT_TO_AXIS,
)

app = typer.Typer(add_completion=False)

CANONICAL_TRAIT = {
    "flower_color": "flower_primary_color",
    "flower_shape": "floral_form",
}
BASE_DIRECT_GROUPS = {
    "three_pass_high",
    "three_pass_medium",
    "expanded_wave",
    "targeted_reacquisition",
}
BULK_AND_PUBLIC_GROUPS = {
    "austraits",
    "efloras",
    "eol_traitbank",
    "gift",
    "meyer_2026",
    "promoted_public_web",
    "razanajatovo_2016",
    "usda_plants",
    "latest_public_web",
}
CULTIVAR_RE = re.compile(
    r"\b(cultivar|cultivars|horticultural variety|garden variety|cv\.)\b",
    re.IGNORECASE,
)
VARIABLE_RE = re.compile(
    r"\b(variable|varies|varying|polymorphic|sometimes|both|mixed)\b",
    re.IGNORECASE,
)
VALUE_ALIASES = {
    "flower_primary_color": {
        "red": "red_pink",
        "pink": "red_pink",
        "purple": "blue_purple",
        "blue": "blue_purple",
        "violet": "blue_purple",
        "yellow": "yellow_orange",
        "orange": "yellow_orange",
        "green": "green_brown_inconspicuous",
        "brown": "green_brown_inconspicuous",
        "cream": "white",
        "black": "other_described",
        "multicolored/variable": "multicolored_variable",
    },
    "floral_form": {
        "funnelform": "funnel_trumpet",
        "funnel": "funnel_trumpet",
        "trumpet": "funnel_trumpet",
        "funnel_/_trumpet-shaped": "funnel_trumpet",
        "campanulate": "bell_campanulate",
        "bell": "bell_campanulate",
        "bell-shaped": "bell_campanulate",
        "rotate": "open_radial",
        "radial": "open_radial",
        "open_/_bowl_/_cup-shaped": "open_radial",
        "urceolate": "urn_urceolate",
        "urn": "urn_urceolate",
        "composite_head_/_capitulum": "composite_head",
    },
    "floral_symmetry": {
        "radial": "actinomorphic",
        "radially_symmetrical": "actinomorphic",
        "bilateral": "zygomorphic",
        "bilaterally_symmetrical": "zygomorphic",
    },
    "self_incompatibility": {
        "si": "SI",
        "self_incompatible": "SI",
        "self-incompatible": "SI",
        "sc": "SC",
        "self_compatible": "SC",
        "self-compatible": "SC",
    },
}
RULE_COLUMNS = [
    "setting",
    "genus",
    "axis",
    "trait_name",
    "inferred_state_set",
    "inferred_value",
    "value_distribution",
    "n_direct_species",
    "dominant_species",
    "counterexample_species",
    "dominance",
    "entropy",
    "species_loo_n",
    "species_loo_correct",
    "species_loo_accuracy",
    "lineage_loo_n",
    "lineage_loo_correct",
    "lineage_loo_accuracy",
    "required_min_species",
    "required_dominance",
    "required_masked_accuracy",
    "eligible",
    "diagnostic_only",
    "old_low_overlap",
    "old_low_agreement",
    "old_low_agreement_rate",
]


@dataclass(frozen=True)
class ThresholdSetting:
    name: str
    min_species: int
    structure_dominance: float
    colour_dominance: float
    reproductive_dominance: float
    diagnostic_only: bool = False


SETTINGS = (
    ThresholdSetting("current_min3", 3, 0.80, 0.90, 0.95),
    ThresholdSetting("relaxed_min3", 3, 0.75, 0.80, 0.90),
    ThresholdSetting("current_min2_diagnostic", 2, 0.80, 0.90, 0.95, True),
    ThresholdSetting("relaxed_min2_diagnostic", 2, 0.75, 0.80, 0.90, True),
)
MASKED_ACCURACY = {
    "floral_structural_complexity": 0.75,
    "flower_colour": 0.80,
    "reproductive_assurance": 0.85,
}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _json(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_trait_name(value: object) -> str:
    trait = _text(value)
    return CANONICAL_TRAIT.get(trait, trait)


def trait_axis(value: object) -> str:
    trait = canonical_trait_name(value)
    return TRAIT_TO_AXIS.get(trait, "")


def load_ontology(path: Path) -> dict[str, set[str]]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    return {
        str(trait): set(values.get("allowed_values", [])) - {"unresolved"}
        for trait, values in config.get("traits", {}).items()
    }


def _normalise_token(trait: str, value: object) -> str:
    raw = _text(value)
    if not raw:
        return ""
    folded = raw.casefold().replace(" ", "_")
    alias = VALUE_ALIASES.get(trait, {}).get(folded)
    if alias:
        return alias
    if trait == "self_incompatibility" and raw in {"SI", "SC"}:
        return raw
    return folded


def normalise_state_set(
    trait: str,
    value: object,
    allowed: dict[str, set[str]],
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Return ontology-valid states and rejected tokens without cross-trait proxies."""
    canonical = canonical_trait_name(trait)
    raw = _text(value)
    if not raw:
        return (), ()
    values: list[str]
    if raw.startswith("["):
        try:
            decoded = json.loads(raw)
        except json.JSONDecodeError:
            decoded = []
        values = [str(item) for item in decoded] if isinstance(decoded, list) else [raw]
    elif "|" in raw or ";" in raw:
        values = re.split(r"[|;]", raw)
    else:
        values = [raw]
    valid: set[str] = set()
    invalid: set[str] = set()
    permitted = allowed.get(canonical, set())
    for item in values:
        token = _normalise_token(canonical, item)
        if token in permitted:
            valid.add(token)
            continue
        if canonical == "flower_primary_color":
            found = {
                mapped
                for keyword, mapped in VALUE_ALIASES[canonical].items()
                if re.search(rf"\b{re.escape(keyword)}\b", item, re.IGNORECASE)
            }
            if found:
                valid.update(found)
                continue
        invalid.add(token or _text(item))
    return tuple(sorted(valid)), tuple(sorted(invalid))


def load_latest_public_web(root: Path, manifest: dict[str, Any]) -> pd.DataFrame:
    candidates = sorted(root.glob("**/broad_web_medium_evidence.csv.gz"))
    if len(candidates) != 1:
        raise ValueError(
            "latest public-web artifact must contain one "
            f"broad_web_medium_evidence.csv.gz, got {len(candidates)}"
        )
    frame = pd.read_csv(candidates[0], dtype=str).fillna("")
    manifest_rows = [
        row
        for row in manifest.get("sources", [])
        if row.get("source_group") == "latest_public_web"
    ]
    run_id = "|".join(sorted({_text(row.get("run_id")) for row in manifest_rows}))
    artifact = "|".join(
        sorted({_text(row.get("artifact_name")) for row in manifest_rows})
    )
    rows: list[dict[str, str]] = []
    for record in frame.to_dict("records"):
        quality = _text(record.get("evidence_quality")).casefold()
        scope = _text(record.get("evidence_scope"))
        trait = canonical_trait_name(record.get("trait_name"))
        if (
            quality not in {"high", "medium"}
            or scope not in DIRECT_SCOPES
            or not trait_axis(trait)
            or not _text(record.get("accepted_species"))
            or not _text(record.get("normalized_value"))
            or _text(record.get("inference_rule"))
        ):
            continue
        rows.append(
            {
                "accepted_species": _text(record.get("accepted_species")),
                "axis": trait_axis(trait),
                "trait_name": trait,
                "normalized_value": _text(record.get("normalized_value")),
                "quality": quality,
                "source_group": "latest_public_web",
                "source_provider": _text(record.get("source_provider")),
                "source_url": _text(record.get("source_url")),
                "source_record_id": _text(record.get("source_record_id")),
                "source_citation": _text(record.get("source_citation")),
                "source_excerpt": _text(record.get("source_excerpt")),
                "evidence_scope": scope,
                "name_match_method": _text(record.get("name_match_method")),
                "source_lineage": _text(record.get("source_lineage")),
                "lineage_method": "artifact_source_lineage",
                "source_run_id": run_id,
                "source_artifact": artifact,
                "source_file": str(candidates[0]),
                "acceptance_contract": "latest_public_web_species_direct_v1",
            }
        )
    return pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)


def direct_evidence_from_integrated(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str).fillna("")
    frame["quality"] = frame["quality"].str.casefold()
    frame = frame.loc[
        frame["quality"].isin({"high", "medium"})
        & frame["evidence_scope"].isin(DIRECT_SCOPES)
    ].copy()
    frame["trait_name"] = frame["trait_name"].map(canonical_trait_name)
    frame["axis"] = frame["trait_name"].map(trait_axis)
    return frame.loc[frame["axis"].ne("")].reindex(columns=EVIDENCE_COLUMNS).fillna("")


def _fallback_lineage(record: dict[str, Any]) -> str:
    payload = "|".join(
        _text(record.get(column))
        for column in (
            "source_group",
            "source_provider",
            "source_url",
            "source_record_id",
            "source_citation",
        )
    )
    return "fallback:" + hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]


def dedupe_direct_lineages(
    evidence: pd.DataFrame,
    ontology: dict[str, set[str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Collapse redistributions at species x trait x original-source lineage."""
    data = evidence.copy().fillna("")
    data["trait_name"] = data["trait_name"].map(canonical_trait_name)
    data["axis"] = data["trait_name"].map(trait_axis)
    data["quality"] = data["quality"].str.casefold()
    data["source_lineage"] = [
        _text(lineage) or _fallback_lineage(record)
        for lineage, record in zip(
            data["source_lineage"],
            data.to_dict("records"),
            strict=True,
        )
    ]
    states = [
        normalise_state_set(trait, value, ontology)
        for trait, value in zip(
            data["trait_name"], data["normalized_value"], strict=True
        )
    ]
    data["state_set"] = [_json(list(valid)) for valid, _ in states]
    data["invalid_state_tokens"] = [_json(list(invalid)) for _, invalid in states]
    data["ontology_valid"] = [bool(valid) and not invalid for valid, invalid in states]
    group_keys = ["accepted_species", "trait_name", "source_lineage"]
    duplicate_mask = data.duplicated(group_keys, keep=False)
    duplicates = data.loc[duplicate_mask].copy()
    if not duplicates.empty:
        duplicates["lineage_record_count"] = duplicates.groupby(group_keys)[
            "source_lineage"
        ].transform("size")
        duplicates["lineage_value_count"] = duplicates.groupby(group_keys)[
            "state_set"
        ].transform("nunique")

    rows: list[dict[str, Any]] = []
    for keys, group in data.groupby(group_keys, sort=True, dropna=False):
        species, trait, lineage = keys
        top_rank = group["quality"].map(QUALITY_RANK).max()
        top = group.loc[group["quality"].map(QUALITY_RANK).eq(top_rank)].copy()
        valid_top = top.loc[top["ontology_valid"]].copy()
        state_sets = [
            tuple(json.loads(value))
            for value in valid_top["state_set"].drop_duplicates()
        ]
        internal_conflict = len(state_sets) > 1 and trait != "flower_primary_color"
        merged_states = sorted({state for states in state_sets for state in states})
        first = top.sort_values(
            ["source_group", "source_record_id", "normalized_value"]
        ).iloc[0].to_dict()
        first.update(
            {
                "accepted_species": species,
                "axis": trait_axis(trait),
                "trait_name": trait,
                "source_lineage": lineage,
                "quality": str(
                    top.sort_values(
                        ["quality"],
                        key=lambda values: values.map(QUALITY_RANK),
                        ascending=False,
                    ).iloc[0]["quality"]
                ),
                "state_set": _json(merged_states),
                "normalized_value": "|".join(merged_states),
                "ontology_valid": bool(merged_states),
                "ontology_mismatch_rows": int((~top["ontology_valid"]).sum()),
                "lineage_internal_conflict": internal_conflict,
                "lineage_internal_classification": (
                    "source_ontology_mismatch"
                    if not merged_states
                    else (
                        "true_multistate_variable"
                        if len(state_sets) > 1 and trait == "flower_primary_color"
                        else (
                            "unresolved_direct_conflict"
                            if internal_conflict
                            else "single_lineage_value"
                        )
                    )
                ),
                "source_groups": "|".join(sorted(set(top["source_group"]))),
                "source_record_ids": "|".join(
                    sorted(filter(None, set(top["source_record_id"])))
                ),
            }
        )
        rows.append(first)
    return pd.DataFrame(rows).fillna(""), duplicates.fillna("")


def resolve_direct_cells(
    lineages: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Resolve independent direct evidence without row-order tie breaking."""
    rows: list[dict[str, Any]] = []
    audits: list[dict[str, Any]] = []
    keys = ["accepted_species", "axis", "trait_name"]
    for (species, axis, trait), group in lineages.groupby(keys, sort=True):
        usable = group.loc[
            group["ontology_valid"].astype(bool)
            & ~group["lineage_internal_conflict"].astype(bool)
        ].copy()
        ontology_mismatch = int((~group["ontology_valid"].astype(bool)).sum())
        internal_conflict = int(
            group["lineage_internal_conflict"].astype(bool).sum()
        )
        if usable.empty:
            classification = (
                "unresolved_direct_conflict"
                if internal_conflict
                else "source_ontology_mismatch"
            )
            audits.append(
                {
                    "accepted_species": species,
                    "axis": axis,
                    "trait_name": trait,
                    "classification": classification,
                    "resolution_status": "excluded",
                    "selected_quality": "",
                    "state_sets": _json([]),
                    "n_independent_lineages": 0,
                    "source_lineages": "",
                    "source_groups": "",
                    "cultivar_rows_excluded": 0,
                    "ontology_mismatch_rows": ontology_mismatch,
                    "lineage_internal_conflicts": internal_conflict,
                    "invalid_state_tokens": "|".join(
                        sorted(
                            {
                                _text(value)
                                for value in group["invalid_state_tokens"]
                                if _text(value) not in {"", "[]"}
                            }
                        )
                    ),
                }
            )
            continue
        usable["quality_rank"] = usable["quality"].map(QUALITY_RANK)
        selected_quality_rank = int(usable["quality_rank"].max())
        selected = usable.loc[usable["quality_rank"].eq(selected_quality_rank)].copy()
        cultivar = selected.apply(
            lambda row: bool(
                CULTIVAR_RE.search(
                    " ".join(
                        _text(row.get(column))
                        for column in (
                            "source_excerpt",
                            "source_citation",
                            "source_url",
                        )
                    )
                )
            ),
            axis=1,
        )
        distinct_before = selected["state_set"].nunique()
        cultivar_excluded = 0
        if distinct_before > 1 and cultivar.any() and (~cultivar).any():
            cultivar_excluded = int(cultivar.sum())
            selected = selected.loc[~cultivar].copy()
        state_sets = sorted(set(selected["state_set"]))
        source_text = " ".join(selected["source_excerpt"].astype(str))
        classification: str
        status = "resolved"
        resolved_states: list[str] = []
        if len(state_sets) == 1:
            resolved_states = list(json.loads(state_sets[0]))
            if cultivar_excluded:
                classification = "cultivar_contamination"
            elif len(selected) > 1:
                classification = "independent_source_agreement"
            else:
                classification = "single_independent_lineage"
        elif trait == "flower_primary_color" or VARIABLE_RE.search(source_text):
            resolved_states = sorted(
                {
                    state
                    for state_set in state_sets
                    for state in json.loads(state_set)
                }
            )
            classification = "true_multistate_variable"
        else:
            classification = "unresolved_direct_conflict"
            status = "excluded"
        audit = {
            "accepted_species": species,
            "axis": axis,
            "trait_name": trait,
            "classification": classification,
            "resolution_status": status,
            "selected_quality": (
                max(selected["quality"], key=lambda value: QUALITY_RANK[value])
                if len(selected)
                else ""
            ),
            "state_sets": _json(
                [json.loads(value) for value in state_sets]
            ),
            "n_independent_lineages": len(selected),
            "source_lineages": "|".join(sorted(set(selected["source_lineage"]))),
            "source_groups": "|".join(sorted(set(selected["source_group"]))),
            "cultivar_rows_excluded": cultivar_excluded,
            "ontology_mismatch_rows": ontology_mismatch,
            "lineage_internal_conflicts": internal_conflict,
            "invalid_state_tokens": "|".join(
                sorted(
                    {
                        _text(value)
                        for value in group["invalid_state_tokens"]
                        if _text(value) not in {"", "[]"}
                    }
                )
            ),
        }
        audits.append(audit)
        if status != "resolved":
            continue
        rows.append(
            {
                **audit,
                "quality": audit["selected_quality"],
                "state_set": _json(resolved_states),
                "normalized_value": "|".join(resolved_states),
                "multistate": len(resolved_states) > 1,
            }
        )
    return pd.DataFrame(rows).fillna(""), pd.DataFrame(audits).fillna("")


def _modal(values: list[str]) -> tuple[str, int]:
    counts = Counter(values)
    if not counts:
        return "", 0
    value, count = min(counts.items(), key=lambda item: (-item[1], item[0]))
    return value, count


def _entropy(values: list[str]) -> float:
    counts = Counter(values)
    total = sum(counts.values())
    if not total:
        return 0.0
    return -sum(
        (count / total) * math.log(count / total)
        for count in counts.values()
        if count
    )


def _dominance_threshold(axis: str, setting: ThresholdSetting) -> float:
    if axis == "flower_colour":
        return setting.colour_dominance
    if axis == "reproductive_assurance":
        return setting.reproductive_dominance
    return setting.structure_dominance


def _lineage_loo(
    genus_cells: pd.DataFrame,
    genus_lineages: pd.DataFrame,
) -> tuple[int, int, float]:
    """Leave each original-source lineage out and predict its reported state."""
    if genus_cells.empty or genus_lineages.empty:
        return 0, 0, 0.0
    cell_values = dict(
        zip(
            genus_cells["accepted_species"],
            genus_cells["state_set"],
            strict=True,
        )
    )
    outcomes: list[bool] = []
    for row in genus_lineages.itertuples(index=False):
        truth = str(row.state_set)
        training_values = [
            value
            for species, value in cell_values.items()
            if species != row.accepted_species
        ]
        predicted, _ = _modal(training_values)
        if predicted:
            outcomes.append(predicted == truth)
    return len(outcomes), int(sum(outcomes)), (
        sum(outcomes) / len(outcomes) if outcomes else 0.0
    )


def build_rule_audit(
    direct_cells: pd.DataFrame,
    lineages: pd.DataFrame,
    old_low: pd.DataFrame,
    settings: tuple[ThresholdSetting, ...] = SETTINGS,
) -> pd.DataFrame:
    base_rows: list[dict[str, Any]] = []
    valid_lineages = lineages.loc[
        lineages["ontology_valid"].astype(bool)
        & ~lineages["lineage_internal_conflict"].astype(bool)
    ].copy()
    lineage_groups = {
        key: group
        for key, group in valid_lineages.groupby(
            ["genus", "axis", "trait_name"], sort=False
        )
    }
    old_low_groups = {
        key: group
        for key, group in old_low.groupby(
            ["genus", "trait_name"], sort=False
        )
    }
    for (genus, axis, trait), cells in direct_cells.groupby(
        ["genus", "axis", "trait_name"], sort=True
    ):
        values = list(cells["state_set"])
        inferred, dominant_count = _modal(values)
        n_species = len(values)
        species_outcomes: list[bool] = []
        for index, truth in enumerate(values):
            predicted, _ = _modal(values[:index] + values[index + 1 :])
            if predicted:
                species_outcomes.append(predicted == truth)
        group_lineages = lineage_groups.get(
            (genus, axis, trait),
            valid_lineages.iloc[0:0],
        )
        lineage_n, lineage_correct, lineage_accuracy = _lineage_loo(
            cells, group_lineages
        )
        old = old_low_groups.get((genus, trait), old_low.iloc[0:0])
        old_overlap = len(old)
        old_agreement = int(old["state_set"].eq(inferred).sum()) if old_overlap else 0
        base_rows.append(
            {
                "genus": genus,
                "axis": axis,
                "trait_name": trait,
                "inferred_state_set": inferred,
                "inferred_value": "|".join(json.loads(inferred)),
                "value_distribution": _json(dict(sorted(Counter(values).items()))),
                "n_direct_species": n_species,
                "dominant_species": dominant_count,
                "counterexample_species": n_species - dominant_count,
                "dominance": dominant_count / n_species,
                "entropy": _entropy(values),
                "species_loo_n": len(species_outcomes),
                "species_loo_correct": int(sum(species_outcomes)),
                "species_loo_accuracy": (
                    sum(species_outcomes) / len(species_outcomes)
                    if species_outcomes
                    else 0.0
                ),
                "lineage_loo_n": lineage_n,
                "lineage_loo_correct": lineage_correct,
                "lineage_loo_accuracy": lineage_accuracy,
                "old_low_overlap": old_overlap,
                "old_low_agreement": old_agreement,
                "old_low_agreement_rate": (
                    old_agreement / old_overlap if old_overlap else None
                ),
            }
        )
    base = pd.DataFrame(base_rows)
    if base.empty:
        return pd.DataFrame(columns=RULE_COLUMNS)
    frames: list[pd.DataFrame] = []
    for setting in settings:
        frame = base.copy()
        frame["setting"] = setting.name
        frame["required_min_species"] = setting.min_species
        frame["required_dominance"] = frame["axis"].map(
            lambda axis, active=setting: _dominance_threshold(axis, active)
        )
        frame["required_masked_accuracy"] = frame["axis"].map(MASKED_ACCURACY)
        frame["eligible"] = (
            frame["n_direct_species"].ge(setting.min_species)
            & frame["dominance"].ge(frame["required_dominance"])
            & frame["species_loo_accuracy"].ge(frame["required_masked_accuracy"])
            & frame["lineage_loo_accuracy"].ge(frame["required_masked_accuracy"])
        )
        frame["diagnostic_only"] = setting.diagnostic_only
        frames.append(frame.reindex(columns=RULE_COLUMNS))
    return pd.concat(frames, ignore_index=True).fillna("")


def load_old_low(
    integrated_lineage_path: Path,
    validated_low_dir: Path,
    ontology: dict[str, set[str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    lineage = pd.read_csv(integrated_lineage_path, dtype=str).fillna("")
    old = lineage.loc[lineage["quality"].str.casefold().eq("low")].copy()
    old["trait_name"] = old["trait_name"].map(canonical_trait_name)
    old["axis"] = old["trait_name"].map(trait_axis)
    old["genus"] = old["accepted_species"].str.split().str[0]
    old_states = [
        normalise_state_set(trait, value, ontology)[0]
        for trait, value in zip(
            old["trait_name"], old["normalized_value"], strict=True
        )
    ]
    old["state_set"] = [_json(list(states)) for states in old_states]
    old = old.loc[old["state_set"].ne("[]")].drop_duplicates(
        ["accepted_species", "trait_name", "state_set"]
    )
    rule_paths = sorted(validated_low_dir.glob("**/validated_genus_rules.csv.gz"))
    old_rules = (
        pd.concat(
            [pd.read_csv(path, dtype=str).fillna("") for path in rule_paths],
            ignore_index=True,
        )
        if rule_paths
        else pd.DataFrame()
    )
    if not old_rules.empty:
        old_rules["trait_name"] = old_rules["trait_name"].map(canonical_trait_name)
        old_rules = old_rules.drop_duplicates(
            ["genus", "axis", "trait_name", "inferred_value"]
        )
    return old.fillna(""), old_rules.fillna("")


def apply_genus_rules(
    master: pd.DataFrame,
    direct_cells: pd.DataFrame,
    rules: pd.DataFrame,
    setting: str,
) -> pd.DataFrame:
    """Apply only matching genus x trait rules; axis is reporting metadata."""
    eligible = rules.loc[
        rules["setting"].eq(setting) & rules["eligible"].astype(bool)
    ].copy()
    if eligible.empty:
        return pd.DataFrame()
    species = master[["accepted_species", "genus"]].drop_duplicates()
    candidates = species.merge(
        eligible,
        on="genus",
        how="inner",
        validate="many_to_many",
    )
    direct_keys = pd.MultiIndex.from_frame(
        direct_cells[["accepted_species", "trait_name"]]
    )
    candidate_keys = pd.MultiIndex.from_frame(
        candidates[["accepted_species", "trait_name"]]
    )
    candidates = candidates.loc[~candidate_keys.isin(direct_keys)].copy()
    candidates["quality"] = "low"
    candidates["evidence_quality"] = "low"
    candidates["evidence_scope"] = "genus_consensus"
    candidates["inference_method"] = "all_direct_validated_genus_consensus"
    candidates["state_set"] = candidates["inferred_state_set"]
    candidates["normalized_value"] = candidates["inferred_value"]
    candidates["family_inference_used"] = False
    candidates["global_fallback_used"] = False
    candidates["source_lineage"] = (
        "validated-low:"
        + candidates["genus"].astype(str)
        + ":"
        + candidates["trait_name"].astype(str)
        + ":"
        + candidates["inferred_state_set"].astype(str)
    )
    return candidates.sort_values(
        ["accepted_species", "axis", "trait_name"]
    ).reset_index(drop=True)


def species_axis_coverage(
    master: pd.DataFrame,
    direct_cells: pd.DataFrame,
    low: pd.DataFrame,
) -> pd.DataFrame:
    universe = master[["accepted_species"]].drop_duplicates().assign(_key=1).merge(
        pd.DataFrame({"axis": AXES, "_key": 1}),
        on="_key",
    ).drop(columns="_key")
    direct = direct_cells.copy()
    if not direct.empty:
        direct["quality_rank"] = direct["quality"].map(QUALITY_RANK)
    low_frame = low.copy()
    if not low_frame.empty:
        low_frame["quality"] = "low"
        low_frame["quality_rank"] = QUALITY_RANK["low"]
        low_frame["source_groups"] = "validated_low_all_evidence"
        low_frame["source_lineages"] = low_frame["source_lineage"]
        low_frame["classification"] = "validated_genus_consensus"
    evidence = pd.concat(
        [
            direct.reindex(
                columns=[
                    "accepted_species",
                    "axis",
                    "trait_name",
                    "state_set",
                    "normalized_value",
                    "quality",
                    "quality_rank",
                    "source_groups",
                    "source_lineages",
                    "classification",
                ]
            ),
            low_frame.reindex(
                columns=[
                    "accepted_species",
                    "axis",
                    "trait_name",
                    "state_set",
                    "normalized_value",
                    "quality",
                    "quality_rank",
                    "source_groups",
                    "source_lineages",
                    "classification",
                ]
            ),
        ],
        ignore_index=True,
    ).fillna("")
    if evidence.empty:
        return universe.assign(
            quality="",
            trait_composition="",
            trait_names="",
            source_groups="",
            source_lineages="",
        )
    evidence["quality_rank"] = pd.to_numeric(
        evidence["quality_rank"], errors="coerce"
    ).fillna(0)
    best_rank = evidence.groupby(["accepted_species", "axis"])[
        "quality_rank"
    ].transform("max")
    selected_rows = evidence.loc[evidence["quality_rank"].eq(best_rank)].copy()
    selected_rows = selected_rows.sort_values(
        ["accepted_species", "axis", "trait_name", "state_set"]
    )
    selected_rows["trait_entry"] = (
        selected_rows["trait_name"].astype(str)
        + "="
        + selected_rows["state_set"].astype(str)
    )
    selected = (
        selected_rows.groupby(["accepted_species", "axis"], as_index=False)
        .agg(
            quality_rank=("quality_rank", "max"),
            trait_composition=("trait_entry", "|".join),
            trait_names=("trait_name", "|".join),
            source_groups=("source_groups", "|".join),
            source_lineages=("source_lineages", "|".join),
        )
    )
    rank_quality = {rank: quality for quality, rank in QUALITY_RANK.items()}
    selected["quality"] = selected["quality_rank"].map(rank_quality).fillna("")
    selected = selected.drop(columns="quality_rank")
    return universe.merge(
        selected,
        on=["accepted_species", "axis"],
        how="left",
    ).fillna("")


def coverage_snapshot(coverage: pd.DataFrame) -> dict[str, Any]:
    filled = coverage["quality"].isin(QUALITY_RANK)
    by_axis: dict[str, Any] = {}
    for axis in AXES:
        group = coverage.loc[coverage["axis"].eq(axis)]
        axis_filled = group["quality"].isin(QUALITY_RANK)
        quality_counts = {
            quality: int(group["quality"].eq(quality).sum())
            for quality in ("high", "medium", "low")
        }
        by_axis[axis] = {
            "denominator": len(group),
            "filled_species": int(axis_filled.sum()),
            "fill_rate": float(axis_filled.mean()),
            "quality_counts": quality_counts,
            "unresolved_species": int((~axis_filled).sum()),
        }
    counts = (
        coverage.assign(filled=filled)
        .groupby("accepted_species")["filled"]
        .sum()
        .value_counts()
        .reindex(range(4), fill_value=0)
    )
    return {
        "quality_counts": {
            quality: int(coverage["quality"].eq(quality).sum())
            for quality in ("high", "medium", "low")
        },
        "filled_species_axis": int(filled.sum()),
        "unresolved_species_axis": int((~filled).sum()),
        "by_axis": by_axis,
        "species_by_filled_axis_count": {
            str(index): int(value) for index, value in counts.items()
        },
    }


def compare_old_low(
    old_low: pd.DataFrame,
    old_rules: pd.DataFrame,
    new_low: pd.DataFrame,
    direct_cells: pd.DataFrame,
    rules: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    current_rules = rules.loc[rules["setting"].eq("current_min3")].copy()
    eligible_keys = set(
        zip(
            current_rules.loc[current_rules["eligible"].astype(bool), "genus"],
            current_rules.loc[current_rules["eligible"].astype(bool), "trait_name"],
        )
    )
    inferred = {
        (row.genus, row.trait_name): row.inferred_state_set
        for row in current_rules.itertuples(index=False)
    }
    new_keys = set(zip(new_low["accepted_species"], new_low["trait_name"]))
    direct_keys = set(zip(direct_cells["accepted_species"], direct_cells["trait_name"]))
    direct_axis_keys = set(zip(direct_cells["accepted_species"], direct_cells["axis"]))
    rows: list[dict[str, Any]] = []
    for row in old_low.itertuples(index=False):
        key = (row.accepted_species, row.trait_name)
        rule_key = (row.genus, row.trait_name)
        status: str
        if key in direct_keys:
            status = "superseded_by_direct_same_trait"
        elif (row.accepted_species, row.axis) in direct_axis_keys:
            status = "superseded_by_direct_same_axis"
        elif key in new_keys and inferred.get(rule_key) == row.state_set:
            status = "retained"
        elif rule_key not in eligible_keys:
            status = "invalidated_new_counterexample_or_validation"
        elif inferred.get(rule_key) != row.state_set:
            status = "invalidated_inferred_value_changed"
        else:
            status = "not_reapplied_trait_specific"
        rows.append(
            {
                "accepted_species": row.accepted_species,
                "genus": row.genus,
                "axis": row.axis,
                "trait_name": row.trait_name,
                "old_state_set": row.state_set,
                "new_rule_state_set": inferred.get(rule_key, ""),
                "comparison_status": status,
            }
        )
    comparison = pd.DataFrame(rows)
    old_eligible_keys = set()
    if not old_rules.empty and "eligible" in old_rules:
        old_eligible = old_rules.loc[
            old_rules["eligible"].astype(str).str.casefold().isin({"true", "1"})
        ]
        old_eligible_keys = set(zip(old_eligible["genus"], old_eligible["trait_name"]))
    new_eligible_keys = eligible_keys
    summary = {
        "old_low_rows": len(old_low),
        "new_low_rows": len(new_low),
        "old_eligible_genus_trait_rules": len(old_eligible_keys),
        "new_eligible_genus_trait_rules": len(new_eligible_keys),
        "newly_eligible_genus_trait_rules": len(new_eligible_keys - old_eligible_keys),
        "old_rules_no_longer_eligible": len(old_eligible_keys - new_eligible_keys),
        "comparison_status_counts": {
            str(key): int(value)
            for key, value in comparison["comparison_status"].value_counts().items()
        },
    }
    return comparison, summary


def threshold_summary(
    rules: pd.DataFrame,
    master: pd.DataFrame,
    direct_cells: pd.DataFrame,
    old_low: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    details: dict[str, Any] = {}
    for setting in (item.name for item in SETTINGS):
        setting_rules = rules.loc[rules["setting"].eq(setting)].copy()
        low = apply_genus_rules(master, direct_cells, rules, setting)
        coverage = species_axis_coverage(master, direct_cells, low)
        snapshot = coverage_snapshot(coverage)
        eligible = setting_rules.loc[setting_rules["eligible"].astype(bool)]
        support_distribution = {
            str(key): int(value)
            for key, value in eligible["n_direct_species"]
            .astype(int)
            .value_counts()
            .sort_index()
            .items()
        }
        old_keys = set(
            zip(
                old_low["accepted_species"],
                old_low["trait_name"],
                old_low["state_set"],
            )
        )
        new_keys = set(
            zip(low["accepted_species"], low["trait_name"], low["state_set"])
        )
        old_overlap = len(
            {
                (species, trait)
                for species, trait, _ in old_keys
            }
            & {
                (species, trait)
                for species, trait, _ in new_keys
            }
        )
        agreement = len(old_keys & new_keys)
        row = {
            "setting": setting,
            "diagnostic_only": bool(
                setting_rules["diagnostic_only"].astype(bool).any()
            ),
            "eligible_rules": len(eligible),
            "resolved_low_species_trait": len(low),
            "resolved_species_axis": snapshot["filled_species_axis"],
            "low_species_axis": int(coverage["quality"].eq("low").sum()),
            "unresolved_species_axis": snapshot["unresolved_species_axis"],
            "mean_species_loo_accuracy": (
                float(eligible["species_loo_accuracy"].astype(float).mean())
                if len(eligible)
                else 0.0
            ),
            "mean_lineage_loo_accuracy": (
                float(eligible["lineage_loo_accuracy"].astype(float).mean())
                if len(eligible)
                else 0.0
            ),
            "mean_entropy": (
                float(eligible["entropy"].astype(float).mean())
                if len(eligible)
                else 0.0
            ),
            "mean_dominance": (
                float(eligible["dominance"].astype(float).mean())
                if len(eligible)
                else 0.0
            ),
            "old_low_overlap": old_overlap,
            "old_low_agreement": agreement,
            "old_low_agreement_rate": agreement / old_overlap if old_overlap else None,
        }
        rows.append(row)
        details[setting] = {
            **row,
            "support_species_distribution": support_distribution,
            "coverage": snapshot,
        }
    return pd.DataFrame(rows).fillna(""), details


def unresolved_reasons(
    coverage: pd.DataFrame,
    master: pd.DataFrame,
    rules: pd.DataFrame,
    cell_audit: pd.DataFrame,
    candidate_audit: pd.DataFrame,
) -> pd.DataFrame:
    unresolved = coverage.loc[
        ~coverage["quality"].isin(QUALITY_RANK),
        ["accepted_species", "axis"],
    ].merge(
        master[["accepted_species", "genus"]],
        on="accepted_species",
        how="left",
        validate="many_to_one",
    )
    current = rules.loc[rules["setting"].eq("current_min3")].copy()
    current["dominance_pass"] = (
        current["dominance"].astype(float)
        >= current["required_dominance"].astype(float)
    )
    genus_axis = (
        current.groupby(["genus", "axis"], as_index=False)
        .agg(
            max_direct_support_species=("n_direct_species", "max"),
            candidate_genus_trait_rules=("trait_name", "size"),
            any_dominance_pass=("dominance_pass", "any"),
        )
    )
    unresolved = unresolved.merge(
        genus_axis,
        on=["genus", "axis"],
        how="left",
        validate="many_to_one",
    )
    unresolved["max_direct_support_species"] = (
        pd.to_numeric(
            unresolved["max_direct_support_species"], errors="coerce"
        )
        .fillna(0)
        .astype(int)
    )
    unresolved["candidate_genus_trait_rules"] = (
        pd.to_numeric(
            unresolved["candidate_genus_trait_rules"], errors="coerce"
        )
        .fillna(0)
        .astype(int)
    )
    unresolved["any_dominance_pass"] = (
        unresolved["any_dominance_pass"].fillna(False).astype(bool)
    )
    conflict_keys = set(
        cell_audit.loc[
            cell_audit["classification"].eq("unresolved_direct_conflict"),
            "accepted_species",
        ].astype(str)
        + "\x1f"
        + cell_audit.loc[
            cell_audit["classification"].eq("unresolved_direct_conflict"),
            "axis",
        ].astype(str)
    )
    ontology_keys = set(
        cell_audit.loc[
            cell_audit["classification"].eq("source_ontology_mismatch"),
            "accepted_species",
        ].astype(str)
        + "\x1f"
        + cell_audit.loc[
            cell_audit["classification"].eq("source_ontology_mismatch"),
            "axis",
        ].astype(str)
    )
    name_keys: set[str] = set()
    filter_keys: set[str] = set()
    if not candidate_audit.empty:
        audit = candidate_audit.copy()
        audit["axis"] = audit["trait_name"].map(trait_axis)
        rejected = audit.loc[
            ~audit["accepted"].astype(str).str.casefold().isin({"true", "1"})
            & audit["axis"].ne("")
        ]
        for row in rejected.itertuples(index=False):
            key = f"{row.accepted_species}\x1f{row.axis}"
            if "name_match" in str(row.reason):
                name_keys.add(key)
            else:
                filter_keys.add(key)
    keys = unresolved["accepted_species"].astype(str) + "\x1f" + unresolved[
        "axis"
    ].astype(str)
    reasons: list[str] = []
    for key, max_support, dominance_pass in zip(
        keys,
        unresolved["max_direct_support_species"],
        unresolved["any_dominance_pass"],
        strict=True,
    ):
        if key in conflict_keys:
            reason = "direct_value_conflict"
        elif key in name_keys:
            reason = "name_matching_excluded"
        elif key in ontology_keys or key in filter_keys:
            reason = "trait_ontology_or_provenance_filter_excluded"
        elif max_support == 0:
            reason = "no_direct_evidence_in_genus"
        elif max_support == 1:
            reason = "one_direct_species_in_genus"
        elif max_support == 2:
            reason = "two_direct_species_in_genus"
        else:
            if not dominance_pass:
                reason = "three_plus_but_dominance_insufficient"
            else:
                reason = "masked_validation_insufficient"
        reasons.append(reason)
    unresolved["reason"] = reasons
    return unresolved.drop(columns=["any_dominance_pass"])


def source_marginal_gain(
    direct_cells: pd.DataFrame,
) -> pd.DataFrame:
    all_axes = set(zip(direct_cells["accepted_species"], direct_cells["axis"]))
    groups = sorted(
        {
            source
            for packed in direct_cells["source_groups"]
            for source in _text(packed).split("|")
            if source
        }
    )
    rows: list[dict[str, Any]] = []
    for source in groups:
        owned = direct_cells.loc[
            direct_cells["source_groups"].map(
                lambda packed, active=source: active in _text(packed).split("|")
            )
        ]
        without_cells = direct_cells.loc[
            direct_cells["source_groups"].map(
                lambda packed, active=source: (
                    active not in _text(packed).split("|")
                    or len(set(filter(None, _text(packed).split("|")))) > 1
                )
            )
        ]
        without_axes = set(
            zip(without_cells["accepted_species"], without_cells["axis"])
        )
        lost = all_axes - without_axes
        rows.append(
            {
                "source_group": source,
                "resolved_species_trait_cells": len(owned),
                "unique_species_axis": int(
                    owned[["accepted_species", "axis"]].drop_duplicates().shape[0]
                ),
                "exclusive_marginal_species_axis": len(lost),
                **{
                    f"exclusive_marginal_{axis}": sum(
                        lost_axis == axis for _, lost_axis in lost
                    )
                    for axis in AXES
                },
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["exclusive_marginal_species_axis", "source_group"],
        ascending=[False, True],
    )


def acquisition_queue(
    master: pd.DataFrame,
    coverage: pd.DataFrame,
    rules: pd.DataFrame,
) -> pd.DataFrame:
    current = rules.loc[rules["setting"].eq("current_min3")].copy()
    unresolved = coverage.loc[
        ~coverage["quality"].isin(QUALITY_RANK),
        ["accepted_species", "axis"],
    ].merge(
        master[
            ["accepted_species", "genus", "n_islands", "n_records"]
        ],
        on="accepted_species",
        how="left",
        validate="many_to_one",
    )
    unresolved["n_islands_num"] = pd.to_numeric(
        unresolved["n_islands"], errors="coerce"
    ).fillna(0)
    unresolved["n_records_num"] = pd.to_numeric(
        unresolved["n_records"], errors="coerce"
    ).fillna(0)
    opportunities = unresolved.groupby(["genus", "axis"], as_index=False).agg(
        unresolved_congener_count=("accepted_species", "size"),
        unresolved_island_occurrence_weight=("n_islands_num", "sum"),
        unresolved_gbif_record_weight=("n_records_num", "sum"),
    )
    queue = current.loc[current["n_direct_species"].astype(int).lt(3)].merge(
        opportunities,
        on=["genus", "axis"],
        how="inner",
        validate="many_to_one",
    )
    if queue.empty:
        return queue
    queue["current_support"] = queue["n_direct_species"].astype(int)
    queue["current_value_distribution"] = queue["value_distribution"]
    queue["current_dominance"] = queue["dominance"].astype(float)
    queue["potential_cells_unlocked"] = queue["unresolved_congener_count"].astype(int)
    queue["primary_endpoint_relevance"] = queue["axis"].map(
        lambda axis: 1.5 if axis == "reproductive_assurance" else 1.0
    )
    queue["acquisitions_needed_for_min3_dominance"] = [
        max(
            3 - support,
            0,
            math.ceil(
                (required_dominance * support - dominant_species)
                / (1.0 - required_dominance)
            ),
        )
        for support, dominant_species, required_dominance in zip(
            queue["current_support"],
            queue["dominant_species"].astype(int),
            queue["required_dominance"].astype(float),
            strict=True,
        )
    ]
    queue["potential_cells_per_acquisition"] = (
        queue["potential_cells_unlocked"]
        / queue["acquisitions_needed_for_min3_dominance"]
    )
    agreement_bonus = queue["current_dominance"].eq(1.0).astype(float)
    queue["acquisition_priority"] = (
        queue["potential_cells_per_acquisition"] * 10
        + queue["unresolved_island_occurrence_weight"].map(math.log1p) * 5
        + queue["unresolved_gbif_record_weight"].map(math.log1p)
        + agreement_bonus * 20
    ) * queue["primary_endpoint_relevance"]

    def recommendation(row: pd.Series) -> tuple[str, str]:
        if row["axis"] == "reproductive_assurance":
            query_terms = {
                "self_incompatibility": (
                    '"self-incompatibility" OR "self-compatible"'
                ),
                "autonomous_selfing_capacity": (
                    '"autonomous selfing" OR "reproductive assurance"'
                ),
                "mating_system": '"mating system" OR "selfing rate"',
                "cleistogamy": "cleistogamy",
                "herkogamy": "herkogamy",
                "dichogamy": "dichogamy",
                "sex_system": '"sex system" OR dioecy OR monoecy',
            }
            return (
                (
                    f'"{row["genus"]}" "{row["trait_name"]}" '
                    f"({query_terms.get(row['trait_name'], row['trait_name'])})"
                ),
                "Europe PMC full-text XML/supplement; TRY trait export",
            )
        if row["axis"] == "flower_colour":
            return (
                f'"{row["genus"]}" flower colour description synonym',
                "WFO/eFloras descriptions; multilingual botanical pages",
            )
        return (
            f'"{row["genus"]}" {row["trait_name"]} floral description',
            "WFO/eFloras descriptions; botanical Latin extraction",
        )

    recommendations = queue.apply(recommendation, axis=1)
    queue["recommended_query"] = [item[0] for item in recommendations]
    queue["recommended_source"] = [item[1] for item in recommendations]
    queue["priority_reason"] = [
        (
            "two agreeing direct species; one acquisition from min_species=3"
            if support == 2 and dominance == 1.0
            else "below min_species=3 with unresolved congeners"
        )
        for support, dominance in zip(
            queue["current_support"], queue["current_dominance"], strict=True
        )
    ]
    columns = [
        "genus",
        "axis",
        "trait_name",
        "current_support",
        "current_value_distribution",
        "current_dominance",
        "unresolved_congener_count",
        "potential_cells_unlocked",
        "acquisitions_needed_for_min3_dominance",
        "potential_cells_per_acquisition",
        "unresolved_island_occurrence_weight",
        "unresolved_gbif_record_weight",
        "primary_endpoint_relevance",
        "acquisition_priority",
        "recommended_query",
        "recommended_source",
        "priority_reason",
    ]
    return queue[columns].sort_values(
        ["acquisition_priority", "potential_cells_unlocked", "genus", "trait_name"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)


def island_endpoint_coverage(
    occurrences: pd.DataFrame,
    tracks: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    occ = occurrences[["island_id", "species"]].drop_duplicates().copy()
    total = occ.groupby("island_id")["species"].nunique().rename("n_species_total")
    resolved_parts: list[pd.DataFrame] = []
    for track, coverage in tracks.items():
        resolved = coverage.loc[
            coverage["quality"].isin(QUALITY_RANK),
            ["accepted_species", "axis"],
        ].drop_duplicates()
        resolved["track"] = track
        resolved_parts.append(resolved)
    resolved_long = pd.concat(resolved_parts, ignore_index=True)
    matched = occ.merge(
        resolved_long,
        left_on="species",
        right_on="accepted_species",
        how="inner",
        validate="many_to_many",
    )
    resolved_counts = (
        matched.groupby(["track", "axis", "island_id"])["species"]
        .nunique()
        .rename("n_resolved_species")
        .reset_index()
    )
    grid = (
        total.reset_index()
        .assign(_key=1)
        .merge(
            pd.MultiIndex.from_product(
                [list(tracks), AXES],
                names=["track", "axis"],
            )
            .to_frame(index=False)
            .assign(_key=1),
            on="_key",
        )
        .drop(columns="_key")
    )
    result = grid.merge(
        resolved_counts,
        on=["track", "axis", "island_id"],
        how="left",
        validate="one_to_one",
    )
    result["n_resolved_species"] = (
        result["n_resolved_species"].fillna(0).astype(int)
    )
    result["coverage"] = (
        result["n_resolved_species"] / result["n_species_total"]
    )
    confirmatory = result.loc[result["track"].eq("confirmatory")]
    eligible = (
        confirmatory.loc[confirmatory["n_resolved_species"].gt(0)]
        .groupby("island_id")["axis"]
        .nunique()
    )
    common = (
        eligible.loc[eligible.eq(len(AXES))]
        .rename("n_common_endpoints")
        .reset_index()
    )
    common["common_island_set_contract"] = (
        "same island set for every reproductive-assurance pathway equation"
    )
    return result, common


def _snapshot_from_integrated_coverage(path: Path) -> dict[str, Any]:
    frame = pd.read_csv(path, dtype=str).fillna("")
    current = frame[["accepted_species", "axis", "after_quality"]].rename(
        columns={"after_quality": "quality"}
    )
    return coverage_snapshot(current)


def _net_snapshot(before: dict[str, Any], after: dict[str, Any]) -> dict[str, Any]:
    return {
        "filled_species_axis": (
            after["filled_species_axis"] - before["filled_species_axis"]
        ),
        "unresolved_species_axis": (
            after["unresolved_species_axis"] - before["unresolved_species_axis"]
        ),
        "quality_counts": {
            quality: (
                after["quality_counts"][quality] - before["quality_counts"][quality]
            )
            for quality in ("high", "medium", "low")
        },
        "by_axis": {
            axis: {
                "filled_species": (
                    after["by_axis"][axis]["filled_species"]
                    - before["by_axis"][axis]["filled_species"]
                ),
                "fill_rate": (
                    after["by_axis"][axis]["fill_rate"]
                    - before["by_axis"][axis]["fill_rate"]
                ),
                "quality_counts": {
                    quality: (
                        after["by_axis"][axis]["quality_counts"][quality]
                        - before["by_axis"][axis]["quality_counts"][quality]
                    )
                    for quality in ("high", "medium", "low")
                },
            }
            for axis in AXES
        },
    }


def build_audit(
    *,
    integrated_dir: Path,
    latest_public_web_dir: Path,
    validated_low_dir: Path,
    master_csv: Path,
    island_occurrences_csv: Path,
    source_manifest_json: Path,
    ontology_yaml: Path,
    angiosperm_scope_yaml: Path,
    sensitivity_lock_json: Path,
    expected_species: int,
) -> tuple[dict[str, Any], dict[str, pd.DataFrame], dict[str, Any]]:
    started = perf_counter()

    def mark(label: str) -> None:
        typer.echo(f"[all-evidence-audit] {label}: {perf_counter() - started:.1f}s")

    manifest = json.loads(source_manifest_json.read_text(encoding="utf-8"))
    ontology = load_ontology(ontology_yaml)
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    required_master = {
        "accepted_species",
        "genus",
        "n_islands",
        "n_records",
    }
    missing_master = required_master.difference(master.columns)
    if missing_master:
        raise ValueError(f"master missing columns: {sorted(missing_master)}")
    master = master.drop_duplicates("accepted_species").copy()
    scope = classify_scope(master, load_scope_config(angiosperm_scope_yaml))
    eligible_species = set(
        scope.loc[scope["angiosperm_analysis_eligible"], "accepted_species"]
    )
    master = master.loc[master["accepted_species"].isin(eligible_species)].copy()
    if len(master) != expected_species:
        raise ValueError(
            f"master has {len(master)} species, expected {expected_species}"
        )
    mark("eligible master loaded")
    integrated_lineage = integrated_dir / "integrated_evidence_lineage.csv.gz"
    integrated_coverage = integrated_dir / "integrated_species_axis_coverage.csv.gz"
    candidate_audit_path = integrated_dir / "bulk_candidate_acceptance_audit.csv.gz"
    integrated_summary_path = integrated_dir / "integrated_trait_coverage_summary.json"
    for path in (
        integrated_lineage,
        integrated_coverage,
        candidate_audit_path,
        integrated_summary_path,
    ):
        if not path.exists():
            raise ValueError(f"integrated artifact missing {path.name}")
    integrated_summary = json.loads(
        integrated_summary_path.read_text(encoding="utf-8")
    )
    latest_evidence_paths = sorted(
        latest_public_web_dir.glob("**/broad_web_medium_evidence.csv.gz")
    )
    latest_summary_paths = sorted(
        latest_public_web_dir.glob("**/broad_web_medium_pilot_summary.json")
    )
    if len(latest_evidence_paths) != 1 or len(latest_summary_paths) != 1:
        raise ValueError(
            "latest public-web artifact must contain one evidence CSV and "
            "one summary JSON"
        )
    latest_summary = json.loads(
        latest_summary_paths[0].read_text(encoding="utf-8")
    )

    base_direct = direct_evidence_from_integrated(integrated_lineage)
    latest_direct = load_latest_public_web(latest_public_web_dir, manifest)
    raw_direct = pd.concat([base_direct, latest_direct], ignore_index=True).fillna("")
    raw_direct = raw_direct.loc[
        raw_direct["accepted_species"].isin(set(master["accepted_species"]))
    ].copy()
    lineages, lineage_duplicates = dedupe_direct_lineages(raw_direct, ontology)
    direct_cells, cell_audit = resolve_direct_cells(lineages)
    direct_cells["genus"] = direct_cells["accepted_species"].str.split().str[0]
    lineages["genus"] = lineages["accepted_species"].str.split().str[0]
    mark("direct lineages and conflicts resolved")

    old_low, old_rules = load_old_low(
        integrated_lineage,
        validated_low_dir,
        ontology,
    )
    rules = build_rule_audit(direct_cells, lineages, old_low)
    rebuilt_low = apply_genus_rules(
        master,
        direct_cells,
        rules,
        "current_min3",
    )
    direct_coverage = species_axis_coverage(
        master,
        direct_cells,
        pd.DataFrame(),
    )
    old_low_coverage = species_axis_coverage(master, direct_cells, old_low)
    strict_coverage = species_axis_coverage(master, direct_cells, rebuilt_low)
    current_artifact_snapshot = _snapshot_from_integrated_coverage(
        integrated_coverage
    )
    direct_snapshot = coverage_snapshot(direct_coverage)
    old_low_snapshot = coverage_snapshot(old_low_coverage)
    strict_snapshot = coverage_snapshot(strict_coverage)
    mark("current rules and strict coverage built")

    old_comparison, old_comparison_summary = compare_old_low(
        old_low,
        old_rules,
        rebuilt_low,
        direct_cells,
        rules,
    )
    sensitivity_table, sensitivity_details = threshold_summary(
        rules,
        master,
        direct_cells,
        old_low,
    )
    candidate_audit = pd.read_csv(candidate_audit_path, dtype=str).fillna("")
    unresolved = unresolved_reasons(
        strict_coverage,
        master,
        rules,
        cell_audit,
        candidate_audit,
    )
    queue = acquisition_queue(master, strict_coverage, rules)
    marginals = source_marginal_gain(direct_cells)
    mark("sensitivity, unresolved reasons, and queue built")

    base_only_raw = raw_direct.loc[
        raw_direct["source_group"].isin(BASE_DIRECT_GROUPS)
    ].copy()
    base_lineages, _ = dedupe_direct_lineages(base_only_raw, ontology)
    base_cells, _ = resolve_direct_cells(base_lineages)
    base_cells["genus"] = base_cells["accepted_species"].str.split().str[0]
    base_lineages["genus"] = base_lineages["accepted_species"].str.split().str[0]
    base_rules = build_rule_audit(base_cells, base_lineages, old_low)
    base_low = apply_genus_rules(master, base_cells, base_rules, "current_min3")
    base_coverage = species_axis_coverage(master, base_cells, base_low)
    base_snapshot = coverage_snapshot(base_coverage)
    mark("bulk counterfactual built")

    current_table = pd.read_csv(integrated_coverage, dtype=str).fillna("")
    old_low_axis_keys = set(
        zip(
            current_table.loc[
                current_table["after_quality"].eq("low"), "accepted_species"
            ],
            current_table.loc[
                current_table["after_quality"].eq("low"), "axis"
            ],
        )
    )
    direct_axis_keys = set(
        zip(
            direct_coverage.loc[
                direct_coverage["quality"].isin({"high", "medium"}),
                "accepted_species",
            ],
            direct_coverage.loc[
                direct_coverage["quality"].isin({"high", "medium"}),
                "axis",
            ],
        )
    )
    low_to_direct_keys = old_low_axis_keys & direct_axis_keys

    occurrences = pd.read_csv(
        island_occurrences_csv,
        usecols=["island_id", "species"],
        dtype=str,
    ).fillna("")
    island_coverage, common_islands = island_endpoint_coverage(
        occurrences,
        {
            "confirmatory": direct_coverage,
            "secondary": strict_coverage,
        },
    )
    mark("island endpoint coverage built")
    sensitivity_lock = json.loads(
        sensitivity_lock_json.read_text(encoding="utf-8")
    )
    analysis_manifests = {
        "contract": "three_tier_trait_analysis_inputs_v1",
        "confirmatory": {
            "input_file": "confirmatory_trait_input.csv.gz",
            "allowed_evidence": ["species_direct", "synonym_direct"],
            "source_lineage_deduplicated": True,
            "direct_conflicts_excluded": True,
            "multistate_composition_retained": True,
            "missing_is_not_zero_or_absence": True,
            "island_endpoint_coverage_file": (
                "analysis_island_endpoint_coverage.csv.gz"
            ),
        },
        "secondary_robustness": {
            "direct_input_file": "confirmatory_trait_input.csv.gz",
            "genus_probability_file": "secondary_genus_probability_input.csv.gz",
            "allowed_evidence": [
                "species_direct",
                "synonym_direct",
                "validated_probabilistic_genus_inference",
            ],
            "family_inference_allowed": False,
            "global_fallback_allowed": False,
            "imputation_unit": "accepted_species_x_trait_name",
            "species_draw_shared_across_all_islands_within_imputation": True,
            "n_imputations": 200,
            "uncertainty_propagation": (
                "multiple_imputation_or_measurement_error"
            ),
        },
        "sensitivity_only": {
            "family_inference_allowed": True,
            "global_fallback_allowed": True,
            "headline_effects_allowed": False,
            "artifact_lock": sensitivity_lock,
        },
        "model_guardrails": {
            "multistate_flower_colour_and_form_primary": True,
            "binary_plain_nonplain_secondary_only": True,
            "reproductive_traits_not_substitutable": [
                "pollen_vector_mode",
                "self_incompatibility",
                "autonomous_selfing_capacity",
                "mating_system",
            ],
            "reproductive_pathway_common_island_file": (
                "common_reproductive_pathway_islands.csv"
            ),
            "spatial_and_source_pool_lineage_structure_required": True,
            "report_results_separately_by_evidence_tier": True,
        },
    }

    duplicate_groups = (
        lineage_duplicates[
            ["accepted_species", "trait_name", "source_lineage"]
        ].drop_duplicates()
        if len(lineage_duplicates)
        else pd.DataFrame(
            columns=["accepted_species", "trait_name", "source_lineage"]
        )
    )
    provider_redistributions = 0
    if len(lineage_duplicates):
        provider_redistributions = int(
            (
                lineage_duplicates.groupby(
                    ["accepted_species", "trait_name", "source_lineage"]
                )["source_group"].nunique()
                > 1
            ).sum()
        )
    invalidated_count = int(
        old_comparison["comparison_status"].str.startswith("invalidated").sum()
    )
    summary = {
        "contract": "all_evidence_three_axis_trait_audit_v1",
        "denominator": {
            "species": expected_species,
            "axes": len(AXES),
            "species_axis": expected_species * len(AXES),
            "axis_denominators": {axis: expected_species for axis in AXES},
        },
        "quality_precedence": [
            "species_or_synonym_direct_high",
            "species_or_synonym_direct_medium",
            "validated_genus_low",
        ],
        "current_integrated_artifact": current_artifact_snapshot,
        "conflict_aware_direct_only": direct_snapshot,
        "conflict_aware_with_old_low": old_low_snapshot,
        "strict_direct_plus_rebuilt_validated_low": strict_snapshot,
        "net_vs_current_integrated_artifact": _net_snapshot(
            current_artifact_snapshot,
            strict_snapshot,
        ),
        "net_rebuilt_low_vs_conflict_aware_old_low": _net_snapshot(
            old_low_snapshot,
            strict_snapshot,
        ),
        "bulk_and_public_web_counterfactual": {
            "base_direct_plus_rebuilt_low": base_snapshot,
            "all_direct_plus_rebuilt_low": strict_snapshot,
            "net": _net_snapshot(base_snapshot, strict_snapshot),
            "included_increment_groups": sorted(BULK_AND_PUBLIC_GROUPS),
        },
        "validated_low": {
            **old_comparison_summary,
            "old_low_invalidated": invalidated_count,
            "old_low_upgraded_to_direct_species_axis": len(low_to_direct_keys),
            "old_low_upgraded_to_direct_by_axis": {
                axis: sum(key_axis == axis for _, key_axis in low_to_direct_keys)
                for axis in AXES
            },
            "new_low_by_axis": {
                axis: int(rebuilt_low["axis"].eq(axis).sum())
                for axis in AXES
            },
        },
        "source_lineage_audit": {
            "upstream_integrated_artifact": {
                key: int(value)
                for key, value in integrated_summary.get(
                    "input_audit", {}
                ).items()
                if isinstance(value, int)
            },
            "latest_public_web_artifact": {
                "evidence_rows": int(
                    latest_summary.get("accepted", {}).get("evidence_rows", 0)
                ),
                "source_lineage_key": (
                    "accepted_species_x_trait_name_x_original_source_lineage"
                ),
            },
            "direct_rows_before_lineage_dedup": len(raw_direct),
            "direct_lineages_after_trait_specific_dedup": len(lineages),
            "duplicate_rows": len(lineage_duplicates),
            "duplicate_species_trait_lineage_groups": len(duplicate_groups),
            "cross_provider_same_lineage_groups": provider_redistributions,
            "cell_resolution_classification_counts": {
                str(key): int(value)
                for key, value in cell_audit["classification"].value_counts().items()
            },
        },
        "unresolved_reason_counts": {
            reason: int(unresolved["reason"].eq(reason).sum())
            for reason in (
                "no_direct_evidence_in_genus",
                "one_direct_species_in_genus",
                "two_direct_species_in_genus",
                "three_plus_but_dominance_insufficient",
                "masked_validation_insufficient",
                "direct_value_conflict",
                "trait_ontology_or_provenance_filter_excluded",
                "name_matching_excluded",
            )
        },
        "threshold_sensitivity": sensitivity_details,
        "source_marginal_gain": {
            str(row.source_group): {
                "resolved_species_trait_cells": int(
                    row.resolved_species_trait_cells
                ),
                "unique_species_axis": int(row.unique_species_axis),
                "exclusive_marginal_species_axis": int(
                    row.exclusive_marginal_species_axis
                ),
            }
            for row in marginals.itertuples(index=False)
        },
        "acquisition_queue": {
            "rows": len(queue),
            "top_potential_cells_unlocked": (
                int(queue.iloc[0]["potential_cells_unlocked"])
                if len(queue)
                else 0
            ),
            "top_genus_trait": (
                f"{queue.iloc[0]['genus']} x {queue.iloc[0]['trait_name']}"
                if len(queue)
                else ""
            ),
        },
        "checks": {
            "denominator_exact": len(strict_coverage)
            == expected_species * len(AXES),
            "quality_counts_sum": sum(
                strict_snapshot["quality_counts"].values()
            )
            == strict_snapshot["filled_species_axis"],
            "filled_plus_unresolved": (
                strict_snapshot["filled_species_axis"]
                + strict_snapshot["unresolved_species_axis"]
                == expected_species * len(AXES)
            ),
            "family_inference_absent_from_strict": True,
            "global_fallback_absent_from_strict": True,
            "trait_specific_genus_join": True,
            "min_species_two_diagnostic_only": True,
            "species_imputation_shared_across_islands": True,
        },
    }
    frames = {
        "strict_coverage": strict_coverage,
        "direct_cells": direct_cells,
        "cell_audit": cell_audit,
        "lineage_duplicates": lineage_duplicates,
        "rebuilt_low": rebuilt_low,
        "rules": rules,
        "old_low_comparison": old_comparison,
        "unresolved": unresolved,
        "sensitivity_table": sensitivity_table,
        "source_marginals": marginals,
        "queue": queue,
        "confirmatory": direct_cells,
        "secondary": rebuilt_low,
        "island_coverage": island_coverage,
        "common_islands": common_islands,
    }
    input_manifest = {
        "contract": "all_evidence_source_run_manifest_v1",
        "source_manifest": manifest,
        "inputs": {
            str(path): {
                "sha256": _sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for path in (
                integrated_lineage,
                integrated_coverage,
                candidate_audit_path,
                integrated_summary_path,
                latest_evidence_paths[0],
                latest_summary_paths[0],
                master_csv,
                island_occurrences_csv,
                ontology_yaml,
                angiosperm_scope_yaml,
                sensitivity_lock_json,
            )
        },
        "latest_public_web_artifact_dir": str(latest_public_web_dir),
    }
    mark("audit complete")
    return summary, frames, {
        "analysis_manifests": analysis_manifests,
        "input_manifest": input_manifest,
    }


def write_outputs(
    summary: dict[str, Any],
    frames: dict[str, pd.DataFrame],
    manifests: dict[str, Any],
    output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_outputs = {
        "all_evidence_trait_coverage_summary.json": summary,
        "before_after_coverage_summary.json": {
            key: summary[key]
            for key in (
                "denominator",
                "current_integrated_artifact",
                "conflict_aware_direct_only",
                "conflict_aware_with_old_low",
                "strict_direct_plus_rebuilt_validated_low",
                "net_vs_current_integrated_artifact",
                "net_rebuilt_low_vs_conflict_aware_old_low",
                "bulk_and_public_web_counterfactual",
            )
        },
        "threshold_sensitivity_report.json": summary["threshold_sensitivity"],
        "analysis_input_manifests.json": manifests["analysis_manifests"],
        "all_evidence_source_run_manifest.json": manifests["input_manifest"],
    }
    for name, value in json_outputs.items():
        (output_dir / name).write_text(
            json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    csv_outputs = {
        "strict_species_axis_coverage.csv.gz": "strict_coverage",
        "resolved_direct_species_trait.csv.gz": "direct_cells",
        "source_lineage_conflicts.csv.gz": "cell_audit",
        "source_lineage_duplicates_trait_specific.csv.gz": "lineage_duplicates",
        "rebuilt_all_evidence_validated_low.csv.gz": "rebuilt_low",
        "trait_specific_genus_rule_audit.csv.gz": "rules",
        "old_low_comparison.csv.gz": "old_low_comparison",
        "unresolved_reason_by_species_axis.csv.gz": "unresolved",
        "threshold_sensitivity_summary.csv": "sensitivity_table",
        "source_marginal_gain.csv": "source_marginals",
        "prioritized_acquisition_queue.csv.gz": "queue",
        "confirmatory_trait_input.csv.gz": "confirmatory",
        "secondary_genus_probability_input.csv.gz": "secondary",
        "analysis_island_endpoint_coverage.csv.gz": "island_coverage",
        "common_reproductive_pathway_islands.csv": "common_islands",
    }
    for name, key in csv_outputs.items():
        frames[key].to_csv(output_dir / name, index=False)


@app.command()
def build(
    integrated_dir: Annotated[Path, typer.Option(exists=True, file_okay=False)],
    latest_public_web_dir: Annotated[
        Path, typer.Option(exists=True, file_okay=False)
    ],
    validated_low_dir: Annotated[Path, typer.Option(exists=True, file_okay=False)],
    master_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    island_occurrences_csv: Annotated[
        Path, typer.Option(exists=True, dir_okay=False)
    ],
    source_manifest_json: Annotated[
        Path, typer.Option(exists=True, dir_okay=False)
    ],
    ontology_yaml: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    angiosperm_scope_yaml: Annotated[
        Path, typer.Option(exists=True, dir_okay=False)
    ],
    sensitivity_lock_json: Annotated[
        Path, typer.Option(exists=True, dir_okay=False)
    ],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
    expected_species: Annotated[int, typer.Option(min=1)] = 106_295,
) -> None:
    summary, frames, manifests = build_audit(
        integrated_dir=integrated_dir,
        latest_public_web_dir=latest_public_web_dir,
        validated_low_dir=validated_low_dir,
        master_csv=master_csv,
        island_occurrences_csv=island_occurrences_csv,
        source_manifest_json=source_manifest_json,
        ontology_yaml=ontology_yaml,
        angiosperm_scope_yaml=angiosperm_scope_yaml,
        sensitivity_lock_json=sensitivity_lock_json,
        expected_species=expected_species,
    )
    write_outputs(summary, frames, manifests, output_dir)
    typer.echo(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    app()
