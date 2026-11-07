"""Restore categorical floral profiles without changing the Chapter 1 snapshot.

The frozen species-by-axis snapshot suppresses lower-quality component traits when
another trait in the same axis has higher quality.  For this separate analysis we
retain its existing Direct states, recover missing component traits from the
reviewed species-by-trait ledger, and admit existing validated Low only when no
Direct key exists.  A multistate cell remains one set; sets are never multiplied
into unobserved colour-by-form combinations.
"""

from __future__ import annotations

from collections import Counter
import hashlib
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


TRAITS = (
    "flower_primary_color",
    "floral_form",
    "floral_symmetry",
    "tube_depth_class",
)
TRAIT_SETS = {
    "colour_form": TRAITS[:2],
    "colour_form_symmetry": TRAITS[:3],
    "colour_form_symmetry_depth": TRAITS,
}
KEY = ["accepted_species", "trait_name"]
QUALITY_RANK = {"low": 1, "medium": 2, "high": 3}
NONSPECIFIC = {"other_described", "multicolored_variable"}
METADATA = (
    "quality", "evidence_scope", "source_groups", "source_lineages", "record_source"
)


def _tokens(value: str) -> tuple[str, ...]:
    return tuple(sorted({part.strip() for part in value.split("|") if part.strip()}))


def _read_ledger(
    path: Path,
    *,
    label: str,
    allowed: dict[str, set[str]],
    raw: bool = False,
) -> tuple[pd.DataFrame, set[str], dict[str, Any]]:
    frame = pd.read_csv(path, dtype=str).fillna("")
    required = set(KEY + ["normalized_value", "quality", "source_groups", "source_lineages"])
    if raw:
        required.add("resolution_status")
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{label}: missing columns {sorted(missing)}")
    for column in required:
        frame[column] = frame[column].str.strip()
    if frame["accepted_species"].eq("").any():
        raise ValueError(f"{label}: blank accepted_species")
    if frame.duplicated(KEY).any():
        raise ValueError(f"{label}: duplicate accepted_species x trait_name key")

    species = set(frame["accepted_species"])
    frame = frame.loc[frame["trait_name"].isin(TRAITS)].copy()
    frame["quality"] = frame["quality"].str.lower()
    frame["_states"] = frame["normalized_value"].map(_tokens)
    frame["normalized_value"] = frame["_states"].map("|".join)
    qualities = {"high", "medium", "low"} if label == "snapshot_all" else {"high", "medium"}

    def exclusion(row: Any) -> str:
        if row.quality not in qualities:
            return "ineligible_quality"
        if raw and row.resolution_status.lower() != "resolved":
            return "nonresolved"
        states = set(_tokens(row.normalized_value))
        if not states or "unresolved" in states:
            return "missing_or_unresolved"
        if not states.issubset(allowed[row.trait_name]):
            return "invalid_ontology"
        if states & NONSPECIFIC:
            return "nonspecific_description"
        if not row.source_groups or not row.source_lineages:
            return "missing_provenance"
        return ""

    frame["_exclusion"] = [exclusion(row) for row in frame.itertuples(index=False)]
    frame["_valid"] = frame["_exclusion"].eq("")
    frame["evidence_scope"] = frame["quality"].map(
        lambda quality: "validated_low" if quality == "low" else "direct"
    )
    counts = Counter(frame.loc[~frame["_valid"], "_exclusion"])
    validation = {
        "focal_cells": int(len(frame)),
        "valid_focal_cells": int(frame["_valid"].sum()),
        "invalid_ontology_cells": int(counts.get("invalid_ontology", 0)),
        "nonresolved_cells": int(counts.get("nonresolved", 0)),
        "excluded_by_reason": dict(sorted(counts.items())),
        "by_trait": {
            trait: {
                "input_cells": int(frame["trait_name"].eq(trait).sum()),
                "valid_cells": int((frame["trait_name"].eq(trait) & frame["_valid"]).sum()),
            }
            for trait in TRAITS
        },
    }
    return frame, species, validation


def _key_index(frame: pd.DataFrame) -> pd.MultiIndex:
    return pd.MultiIndex.from_frame(frame[KEY])


def _wide(frame: pd.DataFrame, species: list[str]) -> pd.DataFrame:
    output = pd.DataFrame(index=pd.Index(species, name="accepted_species"))
    for trait in TRAITS:
        part = frame.loc[frame["trait_name"].eq(trait)].set_index("accepted_species")
        output[trait] = part["normalized_value"].reindex(output.index)
        for column in METADATA:
            values = part[column] if column in part else pd.Series(dtype=object)
            output[f"{trait}__{column}"] = values.reindex(output.index)
    return output.reset_index()


def _coverage(frame: pd.DataFrame) -> dict[str, Any]:
    return {
        "input_species_union": int(len(frame)),
        "by_trait": {trait: int(frame[trait].notna().sum()) for trait in TRAITS},
        **{
            label: int(frame[list(traits)].notna().all(axis=1).sum())
            for label, traits in TRAIT_SETS.items()
        },
    }


def load_profiles(
    direct_ledger_path: Path,
    snapshot_all_path: Path,
    snapshot_direct_path: Path,
    ontology_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Return (All, Direct, audit); profile columns are canonical pipe-delimited sets.

    Rows cover the union of species appearing in the three input ledgers, including
    rows with no eligible floral profile. Reindex to the separately frozen species
    universe before calculating population-wide coverage percentages.

    Unspecified descriptions (including an explicit colour plus unspecified
    variation) and ontology-invalid cells remain missing. Low never replaces a
    Direct key, including a Direct key excluded by these profile-specific rules.
    Existing snapshot Direct states can extend the raw states only when they also
    retain all raw source lineages at equal or higher quality; other disagreements
    stop the analysis for adjudication. No new state union is constructed here.
    """
    ontology = yaml.safe_load(Path(ontology_path).read_text(encoding="utf-8"))["traits"]
    allowed = {trait: set(ontology[trait]["allowed_values"]) for trait in TRAITS}
    raw, raw_species, raw_validation = _read_ledger(
        Path(direct_ledger_path), label="raw_direct", allowed=allowed, raw=True
    )
    snapshot_all, all_species, all_validation = _read_ledger(
        Path(snapshot_all_path), label="snapshot_all", allowed=allowed
    )
    snapshot_direct, direct_species, direct_validation = _read_ledger(
        Path(snapshot_direct_path), label="snapshot_direct", allowed=allowed
    )
    species = sorted(raw_species | all_species | direct_species)

    comparison = raw.merge(
        snapshot_direct, on=KEY, how="outer", suffixes=("_raw", "_snapshot"),
        validate="one_to_one", indicator=True,
    )
    shared = comparison.loc[comparison["_merge"].eq("both")].copy()
    differences = shared.loc[
        shared["normalized_value_raw"].ne(shared["normalized_value_snapshot"])
    ]
    supported_supersets = 0
    difference_records = []
    for row in differences.itertuples(index=False):
        states_extend = set(_tokens(row.normalized_value_raw)).issubset(
            _tokens(row.normalized_value_snapshot)
        )
        lineages_extend = set(_tokens(row.source_lineages_raw)).issubset(
            _tokens(row.source_lineages_snapshot)
        )
        quality_retained = QUALITY_RANK.get(row.quality_snapshot, 0) >= QUALITY_RANK.get(
            row.quality_raw, 0
        )
        if not (states_extend and lineages_extend and quality_retained):
            raise ValueError(
                "unexplained direct disagreement: "
                f"{row.accepted_species} / {row.trait_name}: "
                f"raw={row.normalized_value_raw!r}, snapshot={row.normalized_value_snapshot!r}"
            )
        supported_supersets += 1
        difference_records.append({
            "accepted_species": row.accepted_species,
            "trait_name": row.trait_name,
            "raw_states": row.normalized_value_raw,
            "snapshot_states": row.normalized_value_snapshot,
            "raw_source_lineages": row.source_lineages_raw,
            "snapshot_source_lineages": row.source_lineages_snapshot,
            "resolution": "retain_documented_snapshot_superset",
        })

    snapshot_keys = _key_index(snapshot_direct)
    recovered = raw.loc[raw["_valid"] & ~_key_index(raw).isin(snapshot_keys)].copy()
    recovered["record_source"] = "raw_direct_recovery"
    retained = snapshot_direct.loc[snapshot_direct["_valid"]].copy()
    retained["record_source"] = "snapshot_direct"
    selected_direct = pd.concat([retained, recovered], ignore_index=True)

    # A missing/ambiguous Direct value is not permission to substitute genus Low.
    direct_keys = _key_index(pd.concat([raw, snapshot_direct], ignore_index=True))
    low = snapshot_all.loc[snapshot_all["quality"].eq("low")].copy()
    low_blocked = _key_index(low).isin(direct_keys)
    selected_low = low.loc[low["_valid"] & ~low_blocked].copy()
    selected_low["record_source"] = "snapshot_validated_low"
    selected_all = pd.concat([selected_direct, selected_low], ignore_index=True)
    all_profiles = _wide(selected_all, species)
    direct_profiles = _wide(selected_direct, species)

    audit = {
        "contract": "floral_space_trait_profiles_v1",
        "traits": list(TRAITS),
        "input_validation": {
            "raw_direct": raw_validation,
            "snapshot_all": all_validation,
            "snapshot_direct": direct_validation,
        },
        "comparison": {
            "shared_focal_direct_keys": int(len(shared)),
            "raw_only_focal_keys": int(comparison["_merge"].eq("left_only").sum()),
            "snapshot_only_focal_keys": int(comparison["_merge"].eq("right_only").sum()),
            "state_disagreements": int(len(differences)),
            "snapshot_state_and_lineage_supersets": supported_supersets,
            "disagreement_records": difference_records,
        },
        "selection": {
            "snapshot_direct_cells_retained": int(len(retained)),
            "raw_direct_cells_restored": int(len(recovered)),
            "raw_direct_cells_restored_by_trait": {
                trait: int(recovered["trait_name"].eq(trait).sum()) for trait in TRAITS
            },
            "low_cells_blocked_by_direct_key": int(low_blocked.sum()),
            "low_cells_admitted": int(len(selected_low)),
        },
        "coverage": {
            "snapshot_all": _coverage(_wide(snapshot_all.loc[snapshot_all["_valid"]], species)),
            "snapshot_direct": _coverage(
                _wide(snapshot_direct.loc[snapshot_direct["_valid"]], species)
            ),
            "restored_all": _coverage(all_profiles),
            "restored_direct": _coverage(direct_profiles),
        },
        "policies": {
            "states": "sets retained; no Cartesian product or inferred combinations",
            "direct": "retain snapshot Direct; restore missing traits from reviewed raw Direct",
            "low": "existing snapshot validated Low only where no Direct key exists",
            "generic_and_invalid": "entire cell missing; no silent token deletion or relabelling",
            "denominator": "species union here; use frozen 106295-species universe for percentages",
        },
        "inputs": {
            label: {"path": str(path), "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest()}
            for label, path in (
                ("raw_direct", direct_ledger_path),
                ("snapshot_all", snapshot_all_path),
                ("snapshot_direct", snapshot_direct_path),
                ("ontology", ontology_path),
            )
        },
    }
    return all_profiles, direct_profiles, audit
