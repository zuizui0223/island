"""Materialize the Wave35 direct-evidence and trait-specific Low overlay.

Wave35 is deliberately additive to the frozen Wave34 secondary ledger.  New
species-direct High/Medium evidence may fill an unresolved cell or upgrade an
existing Low/Medium cell.  Validated Low may fill only a still-unresolved cell
and must retain the full ``genus x axis x trait_name`` key.  Family inference
and global fallback are not accepted.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

AXES = (
    "flower_colour",
    "floral_structural_complexity",
    "reproductive_assurance",
)
TRAIT_TO_AXIS = {
    "flower_primary_color": "flower_colour",
    "floral_form": "floral_structural_complexity",
    "floral_symmetry": "floral_structural_complexity",
    "tube_depth_class": "floral_structural_complexity",
    "flower_size_class": "floral_structural_complexity",
    "inflorescence_display": "floral_structural_complexity",
    "self_incompatibility": "reproductive_assurance",
    "autonomous_selfing_capacity": "reproductive_assurance",
    "mating_system": "reproductive_assurance",
    "cleistogamy": "reproductive_assurance",
}
EXPECTED_SPECIES = 106_295
EXPECTED_CELLS = EXPECTED_SPECIES * len(AXES)
QUALITY_RANK = {"": 0, "low": 1, "medium": 2, "high": 3}
GZIP = {"method": "gzip", "mtime": 0}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_states(value: object, *, label: str) -> str:
    try:
        states = json.loads(str(value))
    except json.JSONDecodeError as exc:
        raise ValueError(f"{label} is not valid JSON: {value!r}") from exc
    if not isinstance(states, list) or not states or not all(
        isinstance(state, str) and state for state in states
    ):
        raise ValueError(f"{label} must be a non-empty JSON string list")
    if len(states) != len(set(states)):
        raise ValueError(f"{label} contains duplicate states")
    return json.dumps(sorted(states), ensure_ascii=False, separators=(",", ":"))


def _split_pipe(value: object) -> set[str]:
    return {token for token in str(value).split("|") if token}


def _parse_composition(value: object, *, axis: str) -> dict[str, str]:
    result: dict[str, str] = {}
    if not str(value):
        return result
    for item in str(value).split("|"):
        trait, separator, states = item.partition("=")
        if not separator or TRAIT_TO_AXIS.get(trait) != axis:
            raise ValueError(f"invalid trait composition item for {axis}: {item!r}")
        result[trait] = _canonical_states(states, label=f"composition {trait}")
    return result


def _serialize_composition(composition: dict[str, str]) -> str:
    return "|".join(f"{trait}={composition[trait]}" for trait in sorted(composition))


def _coverage_summary(frame: pd.DataFrame) -> dict[str, Any]:
    resolved = frame["quality"].ne("")
    per_species = frame.assign(_resolved=resolved).groupby("accepted_species")["_resolved"].sum()
    by_axis: dict[str, Any] = {}
    for axis in AXES:
        part = frame.loc[frame["axis"].eq(axis)]
        filled = int(part["quality"].ne("").sum())
        by_axis[axis] = {
            "denominator": EXPECTED_SPECIES,
            "filled_species": filled,
            "fill_rate": filled / EXPECTED_SPECIES,
            "unresolved_species": EXPECTED_SPECIES - filled,
            "quality_counts": {
                quality: int(part["quality"].eq(quality).sum())
                for quality in ("high", "medium", "low")
            },
        }
    return {
        "filled_species_axis": int(resolved.sum()),
        "unresolved_species_axis": int((~resolved).sum()),
        "quality_counts": {
            quality: int(frame["quality"].eq(quality).sum())
            for quality in ("high", "medium", "low")
        },
        "by_axis": by_axis,
        "species_by_filled_axis_count": {
            str(count): int(per_species.eq(count).sum()) for count in range(4)
        },
    }


def _validate_universe(frame: pd.DataFrame, *, label: str) -> None:
    required = {
        "accepted_species",
        "axis",
        "trait_composition",
        "trait_names",
        "source_groups",
        "source_lineages",
        "quality",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{label} lacks columns: {sorted(missing)}")
    if len(frame) != EXPECTED_CELLS:
        raise ValueError(f"{label} has {len(frame)} rows, expected {EXPECTED_CELLS}")
    if frame[["accepted_species", "axis"]].duplicated().any():
        raise ValueError(f"{label} contains duplicate species x axis cells")
    if frame["accepted_species"].nunique() != EXPECTED_SPECIES:
        raise ValueError(f"{label} does not contain {EXPECTED_SPECIES} unique species")
    if set(frame["axis"]) != set(AXES):
        raise ValueError(f"{label} contains an unexpected axis")
    if set(frame["quality"]) - set(QUALITY_RANK):
        raise ValueError(f"{label} contains an unexpected quality")
    counts = frame.groupby("axis").size().to_dict()
    if any(counts.get(axis) != EXPECTED_SPECIES for axis in AXES):
        raise ValueError(f"{label} axis denominators are not exact: {counts}")


def _validate_direct(frame: pd.DataFrame) -> pd.DataFrame:
    required = {
        "accepted_species",
        "axis",
        "trait_name",
        "resolution_status",
        "quality",
        "state_set",
        "source_groups",
        "source_lineages",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"direct evidence lacks columns: {sorted(missing)}")
    if frame[["accepted_species", "trait_name"]].duplicated().any():
        raise ValueError("direct evidence contains duplicate species x trait rows")
    if not frame["resolution_status"].eq("resolved").all():
        raise ValueError("direct evidence contains an unresolved conflict")
    if set(frame["quality"]) - {"high", "medium"}:
        raise ValueError("direct evidence contains a non-direct quality")
    mapped = frame["trait_name"].map(TRAIT_TO_AXIS)
    if mapped.isna().any() or not mapped.eq(frame["axis"]).all():
        raise ValueError("direct evidence contains a cross-trait or cross-axis mapping")
    if frame["source_groups"].eq("").any() or frame["source_lineages"].eq("").any():
        raise ValueError("direct evidence lacks provenance")
    for row in frame.itertuples(index=False):
        _canonical_states(row.state_set, label=f"direct {row.accepted_species} {row.trait_name}")
    return frame


def _validate_low(low: pd.DataFrame, rules: pd.DataFrame) -> pd.DataFrame:
    low_required = {
        "accepted_species",
        "genus",
        "axis",
        "trait_name",
        "inferred_state_set",
        "quality",
        "family_inference_used",
        "global_fallback_used",
        "source_lineage",
    }
    rule_required = {
        "genus",
        "axis",
        "trait_name",
        "inferred_state_set",
        "eligible",
        "n_direct_species",
        "species_loo_accuracy",
        "lineage_loo_accuracy",
    }
    for label, frame, required in (("Low", low, low_required), ("rule", rules, rule_required)):
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"{label} input lacks columns: {sorted(missing)}")
    if low[["accepted_species", "trait_name"]].duplicated().any():
        raise ValueError("Low input contains duplicate species x trait rows")
    if not low["quality"].eq("low").all():
        raise ValueError("Low input contains a non-Low quality")
    if not low["family_inference_used"].astype(str).str.casefold().isin({"false", "0"}).all():
        raise ValueError("family inference is present")
    if not low["global_fallback_used"].astype(str).str.casefold().isin({"false", "0"}).all():
        raise ValueError("global fallback is present")
    if not low["accepted_species"].str.split().str[0].eq(low["genus"]).all():
        raise ValueError("Low input genus does not match accepted_species")
    mapped = low["trait_name"].map(TRAIT_TO_AXIS)
    if mapped.isna().any() or not mapped.eq(low["axis"]).all():
        raise ValueError("Low input contains a cross-trait or cross-axis mapping")
    eligible = rules.loc[rules["eligible"].astype(str).str.casefold().eq("true")].copy()
    if (pd.to_numeric(eligible["n_direct_species"]) < 3).any():
        raise ValueError("eligible rule has support from fewer than three species")
    for column in ("species_loo_accuracy", "lineage_loo_accuracy"):
        if (pd.to_numeric(eligible[column]) < 0.75).any():
            raise ValueError(f"eligible rule contains {column} below 0.75")
    rule_keys = set(
        eligible[["genus", "axis", "trait_name", "inferred_state_set"]]
        .astype(str)
        .itertuples(index=False, name=None)
    )
    low_keys = set(
        low[["genus", "axis", "trait_name", "inferred_state_set"]]
        .astype(str)
        .itertuples(index=False, name=None)
    )
    if not low_keys.issubset(rule_keys):
        raise ValueError("Low input refers to a rule that was not eligible")
    for row in low.itertuples(index=False):
        _canonical_states(
            row.inferred_state_set,
            label=f"Low {row.accepted_species} {row.trait_name}",
        )
    return eligible


def _aggregate_direct(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for (species, axis), group in frame.groupby(["accepted_species", "axis"], sort=True):
        composition = {
            row.trait_name: _canonical_states(row.state_set, label=row.trait_name)
            for row in group.sort_values("trait_name").itertuples(index=False)
        }
        quality = max(group["quality"], key=QUALITY_RANK.__getitem__)
        rows.append(
            {
                "accepted_species": species,
                "axis": axis,
                "trait_composition": _serialize_composition(composition),
                "trait_names": "|".join(sorted(composition)),
                "source_groups": "|".join(
                    sorted({token for value in group["source_groups"] for token in _split_pipe(value)})
                ),
                "source_lineages": "|".join(
                    sorted(
                        {token for value in group["source_lineages"] for token in _split_pipe(value)}
                    )
                ),
                "quality": quality,
            }
        )
    return pd.DataFrame(rows)


def _aggregate_low(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for (species, axis), group in frame.groupby(["accepted_species", "axis"], sort=True):
        composition = {
            row.trait_name: _canonical_states(row.inferred_state_set, label=row.trait_name)
            for row in group.sort_values("trait_name").itertuples(index=False)
        }
        rows.append(
            {
                "accepted_species": species,
                "axis": axis,
                "trait_composition": _serialize_composition(composition),
                "trait_names": "|".join(sorted(composition)),
                "source_groups": "wave35_trait_specific_validated_low",
                "source_lineages": "|".join(sorted(set(group["source_lineage"]))),
                "quality": "low",
            }
        )
    return pd.DataFrame(rows)


def build_wave35_overlay(
    *, baseline_csv: Path, input_dir: Path, output_dir: Path
) -> dict[str, Any]:
    """Build the lossless Wave35 overlay from audited, frozen inputs."""
    direct_path = input_dir / "wave35_resolved_direct_species_trait.csv.gz"
    low_path = input_dir / "wave35_candidate_validated_low_species_trait.csv.gz"
    rule_path = input_dir / "wave35_provider_touched_new_rule_audit.csv.gz"
    old_low_path = input_dir / "wave35_old_low_comparison.csv.gz"
    manifest_path = input_dir / "wave35_source_manifest.json"
    required_paths = (baseline_csv, direct_path, low_path, rule_path, old_low_path, manifest_path)
    missing = [str(path) for path in required_paths if not path.exists()]
    if missing:
        raise ValueError(f"missing inputs: {missing}")

    baseline = pd.read_csv(baseline_csv, dtype=str).fillna("")
    direct = _validate_direct(pd.read_csv(direct_path, dtype=str).fillna(""))
    low = pd.read_csv(low_path, dtype=str).fillna("")
    rules = pd.read_csv(rule_path, dtype=str).fillna("")
    eligible_rules = _validate_low(low, rules)
    old_low = pd.read_csv(old_low_path, dtype=str).fillna("")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    _validate_universe(baseline, label="Wave34 baseline")
    if manifest.get("contract") != "wave35_reproduction_morphology_source_manifest_v1":
        raise ValueError("unexpected source manifest contract")
    manifest_hashes = manifest.get("file_sha256")
    if not isinstance(manifest_hashes, dict) or not manifest_hashes:
        raise ValueError("source manifest lacks file hashes")
    for name, expected_hash in manifest_hashes.items():
        source_path = input_dir / str(name)
        if not source_path.is_file():
            raise ValueError(f"source manifest file is missing: {name}")
        observed_hash = _sha256(source_path)
        if observed_hash != str(expected_hash):
            raise ValueError(
                f"source manifest hash mismatch for {name}: {observed_hash} != {expected_hash}"
            )

    before = _coverage_summary(baseline)
    result = baseline.copy().set_index(["accepted_species", "axis"])
    changes: list[dict[str, str]] = []
    direct_axis = _aggregate_direct(direct)
    for row in direct_axis.itertuples(index=False):
        key = (row.accepted_species, row.axis)
        if key not in result.index:
            raise ValueError(f"direct evidence species x axis is outside the universe: {key}")
        before_quality = str(result.loc[key, "quality"])
        if before_quality in {"", "low"}:
            composition = _parse_composition(row.trait_composition, axis=row.axis)
            groups = _split_pipe(row.source_groups)
            lineages = _split_pipe(row.source_lineages)
        else:
            composition = _parse_composition(result.loc[key, "trait_composition"], axis=row.axis)
            composition.update(_parse_composition(row.trait_composition, axis=row.axis))
            groups = _split_pipe(result.loc[key, "source_groups"]) | _split_pipe(row.source_groups)
            lineages = _split_pipe(result.loc[key, "source_lineages"]) | _split_pipe(
                row.source_lineages
            )
        after_quality = max((before_quality, row.quality), key=QUALITY_RANK.__getitem__)
        result.loc[key, "trait_composition"] = _serialize_composition(composition)
        result.loc[key, "trait_names"] = "|".join(sorted(composition))
        result.loc[key, "source_groups"] = "|".join(sorted(groups))
        result.loc[key, "source_lineages"] = "|".join(sorted(lineages))
        result.loc[key, "quality"] = after_quality
        action = "direct_fill" if before_quality == "" else (
            "direct_upgrade" if QUALITY_RANK[after_quality] > QUALITY_RANK[before_quality] else "direct_enrichment"
        )
        changes.append(
            {
                "accepted_species": row.accepted_species,
                "axis": row.axis,
                "action": action,
                "quality_before": before_quality,
                "quality_after": after_quality,
                "trait_names": row.trait_names,
                "source_groups": row.source_groups,
            }
        )

    low_axis = _aggregate_low(low)
    for row in low_axis.itertuples(index=False):
        key = (row.accepted_species, row.axis)
        if key not in result.index:
            raise ValueError(f"Low evidence species x axis is outside the universe: {key}")
        if result.loc[key, "quality"] != "":
            raise ValueError(f"Low attempts to replace a resolved species x axis cell: {key}")
        for column in (
            "trait_composition",
            "trait_names",
            "source_groups",
            "source_lineages",
            "quality",
        ):
            result.loc[key, column] = getattr(row, column)
        changes.append(
            {
                "accepted_species": row.accepted_species,
                "axis": row.axis,
                "action": "validated_low_fill",
                "quality_before": "",
                "quality_after": "low",
                "trait_names": row.trait_names,
                "source_groups": row.source_groups,
            }
        )

    result = result.reset_index().sort_values(["accepted_species", "axis"]).reset_index(drop=True)
    _validate_universe(result, label="Wave35 coverage")
    after = _coverage_summary(result)
    changes_frame = pd.DataFrame(changes).sort_values(["accepted_species", "axis", "action"])
    comparison = baseline[["accepted_species", "axis", "quality"]].merge(
        result[["accepted_species", "axis", "quality"]],
        on=["accepted_species", "axis"],
        suffixes=("_before", "_after"),
        validate="one_to_one",
    )
    was_filled = comparison["quality_before"].ne("")
    now_filled = comparison["quality_after"].ne("")
    loss = int((was_filled & ~now_filled).sum())
    gain = int((~was_filled & now_filled).sum())
    if loss:
        raise ValueError(f"coverage loss must be zero, observed {loss}")
    if (comparison["quality_after"].map(QUALITY_RANK) < comparison["quality_before"].map(QUALITY_RANK)).any():
        raise ValueError("an existing quality was downgraded")

    action_counts = changes_frame["action"].value_counts().to_dict()
    by_axis_action = {
        axis: {
            action: int(
                ((changes_frame["axis"] == axis) & (changes_frame["action"] == action)).sum()
            )
            for action in ("direct_fill", "validated_low_fill", "direct_upgrade", "direct_enrichment")
        }
        for axis in AXES
    }
    invalidated = (
        int(old_low["comparison_status"].astype(str).str.startswith("invalidated_").sum())
        if "comparison_status" in old_low.columns
        else 0
    )
    summary: dict[str, Any] = {
        "contract": "wave35_lossless_direct_and_trait_specific_low_overlay_v1",
        "formal_wave33_run_id": 32932103226,
        "formal_wave34_run_id": 33065274749,
        "wave34_before": before,
        "wave35_after": after,
        "delta": {
            "gross_gain_species_axis": gain,
            "loss_species_axis": loss,
            "net_gain_species_axis": gain - loss,
            "by_axis_net_gain": {
                axis: after["by_axis"][axis]["filled_species"]
                - before["by_axis"][axis]["filled_species"]
                for axis in AXES
            },
            "action_counts": {key: int(value) for key, value in action_counts.items()},
            "by_axis_action_counts": by_axis_action,
        },
        "validated_low_audit": {
            "eligible_provider_touched_rules": len(eligible_rules),
            "applied_species_trait": len(low),
            "applied_species_axis": len(low_axis),
            "old_low_invalidated_in_strict_rebuild": invalidated,
            "strict_rebuild_is_diagnostic_not_subtractive": True,
        },
        "checks": {
            "fixed_denominator": True,
            "quality_precedence_high_medium_low": True,
            "direct_conflicts_excluded": True,
            "trait_specific_genus_join": True,
            "min_direct_species_three": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
            "baseline_loss_zero": loss == 0,
        },
        "input_hashes": {path.name: _sha256(path) for path in required_paths},
        "source_manifest": manifest,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    coverage_path = output_dir / "wave35_species_axis_coverage.csv.gz"
    direct_output = output_dir / "wave35_direct_species_trait.csv.gz"
    low_output = output_dir / "wave35_validated_low_species_trait.csv.gz"
    change_output = output_dir / "wave35_change_audit.csv.gz"
    result.to_csv(coverage_path, index=False, compression=GZIP)
    direct.to_csv(direct_output, index=False, compression=GZIP)
    low.to_csv(low_output, index=False, compression=GZIP)
    changes_frame.to_csv(change_output, index=False, compression=GZIP)
    summary["artifact_hashes"] = {
        path.name: _sha256(path)
        for path in (coverage_path, direct_output, low_output, change_output)
    }
    (output_dir / "wave35_coverage_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return summary


@app.command("build")
def build_command(
    baseline_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    input_dir: Annotated[Path, typer.Option(exists=True, file_okay=False)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
) -> None:
    """Validate frozen inputs and materialize the Wave35 overlay."""
    summary = build_wave35_overlay(
        baseline_csv=baseline_csv,
        input_dir=input_dir,
        output_dir=output_dir,
    )
    typer.echo(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    app()
