"""Validate and materialize trait-specific secondary-robustness coverage.

This module deliberately does not infer from family or from a global prior.  It
accepts only already-audited ``genus x trait_name`` rules and fills cells that
were unresolved in the frozen baseline.  Direct High/Medium cells are retained
byte-for-byte at the semantic column level and remain distinguishable from Low.
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
GZIP = {"method": "gzip", "mtime": 0}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_list(value: object, *, field: str) -> list[str]:
    try:
        parsed = json.loads(str(value))
    except json.JSONDecodeError as exc:
        raise ValueError(f"{field} is not valid JSON: {value!r}") from exc
    if not isinstance(parsed, list) or not parsed or not all(
        isinstance(item, str) and item for item in parsed
    ):
        raise ValueError(f"{field} must be a non-empty JSON string list")
    return parsed


def _bool_false(series: pd.Series) -> pd.Series:
    return series.astype(str).str.casefold().isin({"false", "0"})


def _coverage_summary(frame: pd.DataFrame) -> dict[str, Any]:
    resolved = frame["quality"].ne("")
    per_species = frame.assign(_resolved=resolved).groupby("accepted_species")["_resolved"].sum()
    quality_counts = {
        quality: int(frame["quality"].eq(quality).sum())
        for quality in ("high", "medium", "low")
    }
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
        "quality_counts": quality_counts,
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
    if set(frame["quality"]) - {"", "high", "medium", "low"}:
        raise ValueError(f"{label} contains an unexpected quality")
    axis_counts = frame.groupby("axis").size().to_dict()
    if any(axis_counts.get(axis) != EXPECTED_SPECIES for axis in AXES):
        raise ValueError(f"{label} axis denominators are not exact: {axis_counts}")


def _validate_frontier(
    baseline: pd.DataFrame,
    candidate_trait: pd.DataFrame,
    candidate_axis: pd.DataFrame,
    rules: pd.DataFrame,
) -> None:
    trait_required = {
        "accepted_species",
        "axis",
        "genus",
        "trait_name",
        "predicted_state_set",
        "n_direct_species",
        "n_source_lineages",
        "species_loo_accuracy",
        "lineage_loo_accuracy",
        "support_source_lineages",
    }
    axis_required = {
        "accepted_species",
        "axis",
        "genus",
        "trait_names",
        "predicted_state_sets",
        "quality",
        "analysis_tier",
        "family_inference",
        "global_fallback",
    }
    rule_required = {
        "genus",
        "axis",
        "trait_name",
        "eligible",
        "n_direct_species",
        "species_loo_accuracy",
        "lineage_loo_accuracy",
    }
    for label, frame, required in (
        ("candidate trait", candidate_trait, trait_required),
        ("candidate axis", candidate_axis, axis_required),
        ("rule frontier", rules, rule_required),
    ):
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"{label} lacks columns: {sorted(missing)}")

    if candidate_trait[["accepted_species", "trait_name"]].duplicated().any():
        raise ValueError("candidate trait contains duplicate species x trait rows")
    if candidate_axis[["accepted_species", "axis"]].duplicated().any():
        raise ValueError("candidate axis contains duplicate species x axis rows")

    mapped_axes = candidate_trait["trait_name"].map(TRAIT_TO_AXIS)
    if mapped_axes.isna().any() or not mapped_axes.eq(candidate_trait["axis"]).all():
        raise ValueError("candidate trait contains a cross-trait or cross-axis mapping")
    if (pd.to_numeric(candidate_trait["n_direct_species"]) < 3).any():
        raise ValueError("candidate trait contains support from fewer than three species")
    if (pd.to_numeric(candidate_trait["n_source_lineages"]) < 1).any():
        raise ValueError("candidate trait lacks an independent source lineage")
    for column in ("species_loo_accuracy", "lineage_loo_accuracy"):
        if (pd.to_numeric(candidate_trait[column]) < 0.8).any():
            raise ValueError(f"candidate trait contains {column} below 0.8")
    for value in candidate_trait["predicted_state_set"]:
        states = _json_list(value, field="predicted_state_set")
        if len(states) > 3 or len(set(states)) != len(states):
            raise ValueError("predicted_state_set must contain one to three unique states")
    for value in candidate_trait["support_source_lineages"]:
        _json_list(value, field="support_source_lineages")

    if not candidate_axis["quality"].eq("low").all():
        raise ValueError("candidate axis contains a non-Low fill")
    if not candidate_axis["analysis_tier"].eq("secondary_robustness").all():
        raise ValueError("candidate axis escaped the secondary robustness tier")
    if not _bool_false(candidate_axis["family_inference"]).all():
        raise ValueError("family inference is present")
    if not _bool_false(candidate_axis["global_fallback"]).all():
        raise ValueError("global fallback is present")

    trait_groups = (
        candidate_trait.groupby(["accepted_species", "axis"], sort=True)["trait_name"]
        .apply(lambda values: sorted(values))
        .to_dict()
    )
    axis_groups = {
        (row.accepted_species, row.axis): sorted(_json_list(row.trait_names, field="trait_names"))
        for row in candidate_axis.itertuples(index=False)
    }
    if trait_groups != axis_groups:
        raise ValueError("species x trait candidates do not reproduce species x axis candidates")

    eligible_rules = rules.loc[rules["eligible"].astype(str).str.casefold().eq("true")].copy()
    rule_keys = set(
        eligible_rules[["genus", "axis", "trait_name"]]
        .astype(str)
        .itertuples(index=False, name=None)
    )
    candidate_rule_keys = set(
        candidate_trait[["genus", "axis", "trait_name"]]
        .astype(str)
        .itertuples(index=False, name=None)
    )
    if not candidate_rule_keys.issubset(rule_keys):
        raise ValueError("candidate trait refers to a rule that was not eligible")

    baseline_quality = baseline.set_index(["accepted_species", "axis"])["quality"]
    candidate_keys = pd.MultiIndex.from_frame(candidate_axis[["accepted_species", "axis"]])
    if not baseline_quality.reindex(candidate_keys).fillna("__missing__").eq("").all():
        raise ValueError("frontier attempts to replace an already-resolved baseline cell")


def build_secondary_coverage(
    *, baseline_dir: Path, frontier_dir: Path, output_dir: Path
) -> dict[str, Any]:
    """Build an additive Low layer after validating every scientific gate."""
    baseline_path = baseline_dir / "secondary_species_axis_coverage.csv.gz"
    candidate_trait_path = frontier_dir / "eligible_candidate_species_trait.csv.gz"
    candidate_axis_path = frontier_dir / "eligible_candidate_species_axis.csv.gz"
    rules_path = frontier_dir / "trait_specific_rule_frontier.csv.gz"
    queue_path = frontier_dir / "prioritized_direct_acquisition_queue.csv.gz"
    required_paths = (
        baseline_path,
        candidate_trait_path,
        candidate_axis_path,
        rules_path,
        queue_path,
    )
    missing_paths = [str(path) for path in required_paths if not path.exists()]
    if missing_paths:
        raise ValueError(f"missing inputs: {missing_paths}")

    baseline = pd.read_csv(baseline_path, dtype=str).fillna("")
    candidate_trait = pd.read_csv(candidate_trait_path, dtype=str).fillna("")
    candidate_axis = pd.read_csv(candidate_axis_path, dtype=str).fillna("")
    rules = pd.read_csv(rules_path, dtype=str).fillna("")
    _validate_universe(baseline, label="Wave33 baseline")
    _validate_frontier(baseline, candidate_trait, candidate_axis, rules)

    before = _coverage_summary(baseline)
    result = baseline.copy().set_index(["accepted_species", "axis"])
    for row in candidate_axis.sort_values(["accepted_species", "axis"]).itertuples(index=False):
        key = (row.accepted_species, row.axis)
        trait_names = sorted(_json_list(row.trait_names, field="trait_names"))
        predicted_sets = [
            _json_list(serialized, field="predicted_state_sets item")
            for serialized in _json_list(row.predicted_state_sets, field="predicted_state_sets")
        ]
        if len(trait_names) != len(predicted_sets):
            raise ValueError(f"trait/state count mismatch for {key}")
        composition = "|".join(
            f"{trait}={json.dumps(states, ensure_ascii=False, separators=(',', ':'))}"
            for trait, states in zip(trait_names, predicted_sets, strict=True)
        )
        nested_lineages = [
            _json_list(serialized, field="support_source_lineages item")
            for serialized in _json_list(
                row.support_source_lineages, field="support_source_lineages"
            )
        ]
        lineages = sorted({lineage for group in nested_lineages for lineage in group})
        result.loc[key, "trait_composition"] = composition
        result.loc[key, "trait_names"] = "|".join(trait_names)
        result.loc[key, "source_groups"] = "wave34_trait_specific_validated_low"
        result.loc[key, "source_lineages"] = "|".join(lineages)
        result.loc[key, "quality"] = "low"

    result = result.reset_index().sort_values(["accepted_species", "axis"]).reset_index(drop=True)
    _validate_universe(result, label="Wave34 secondary coverage")
    after = _coverage_summary(result)
    keys = ["accepted_species", "axis"]
    comparison = baseline[keys + ["quality"]].merge(
        result[keys + ["quality"]], on=keys, suffixes=("_before", "_after"), validate="one_to_one"
    )
    old_resolved = comparison["quality_before"].ne("")
    new_resolved = comparison["quality_after"].ne("")
    loss = int((old_resolved & ~new_resolved).sum())
    gross_gain = int((~old_resolved & new_resolved).sum())
    if loss:
        raise ValueError(f"coverage loss must be zero, observed {loss}")
    retained = comparison.loc[old_resolved]
    if not retained["quality_before"].eq(retained["quality_after"]).all():
        raise ValueError("an existing High/Medium/Low quality was changed")

    by_axis_gain = {
        axis: after["by_axis"][axis]["filled_species"]
        - before["by_axis"][axis]["filled_species"]
        for axis in AXES
    }
    summary = {
        "contract": "wave34_additive_trait_specific_secondary_robustness_v1",
        "formal_baseline_run_id": 32932103226,
        "wave33_baseline": before,
        "wave34_secondary": after,
        "delta": {
            "gross_gain_species_axis": gross_gain,
            "loss_species_axis": loss,
            "net_gain_species_axis": gross_gain - loss,
            "by_axis_net_gain": by_axis_gain,
        },
        "frontier": {
            "candidate_species_trait": len(candidate_trait),
            "candidate_species_axis": len(candidate_axis),
            "eligible_genus_trait_rules": int(
                rules["eligible"].astype(str).str.casefold().eq("true").sum()
            ),
            "queue_rows": len(pd.read_csv(queue_path, dtype=str)),
        },
        "checks": {
            "fixed_denominator": True,
            "candidate_cells_were_unresolved": True,
            "trait_specific_genus_join": True,
            "min_direct_species_three": True,
            "species_loo_at_least_0_8": True,
            "source_lineage_loo_at_least_0_8": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
            "direct_high_medium_unchanged": True,
            "baseline_loss_zero": loss == 0,
        },
        "input_hashes": {path.name: _sha256(path) for path in required_paths},
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    coverage_path = output_dir / "wave34_secondary_species_axis_coverage.csv.gz"
    result.to_csv(coverage_path, index=False, compression=GZIP)
    candidate_trait.to_csv(
        output_dir / "wave34_validated_low_species_trait.csv.gz",
        index=False,
        compression=GZIP,
    )
    summary["artifact_hashes"] = {
        coverage_path.name: _sha256(coverage_path),
        "wave34_validated_low_species_trait.csv.gz": _sha256(
            output_dir / "wave34_validated_low_species_trait.csv.gz"
        ),
    }
    (output_dir / "wave34_secondary_coverage_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return summary


@app.command("build")
def build_command(
    baseline_dir: Annotated[Path, typer.Option(exists=True, file_okay=False)],
    frontier_dir: Annotated[Path, typer.Option(exists=True, file_okay=False)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
) -> None:
    """Validate the frozen inputs and materialize the additive secondary ledger."""
    summary = build_secondary_coverage(
        baseline_dir=baseline_dir,
        frontier_dir=frontier_dir,
        output_dir=output_dir,
    )
    typer.echo(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    app()
