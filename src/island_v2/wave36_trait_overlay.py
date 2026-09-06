"""Build the lossless Wave36 NHM and reproductive-remine coverage overlay."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

from island_v2.wave35_trait_overlay import (
    AXES,
    QUALITY_RANK,
    _aggregate_direct,
    _aggregate_low,
    _parse_composition,
    _serialize_composition,
    _split_pipe,
    _validate_direct,
    _validate_low,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

EXPECTED_SPECIES = 106_295
EXPECTED_CELLS = EXPECTED_SPECIES * len(AXES)
GZIP = {"method": "gzip", "mtime": 0}
TEXT_SUFFIXES = {".csv", ".json", ".jsonl", ".md", ".toml", ".yml", ".yaml"}


def _sha256(path: Path) -> str:
    payload = path.read_bytes()
    if path.suffix.casefold() in TEXT_SUFFIXES:
        payload = payload.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def _local_validate_universe(frame: pd.DataFrame, *, label: str) -> None:
    # The shared validator uses the production constants.  Keep a local gate so
    # the test fixture can exercise the overlay on a tiny exact universe.
    required = {
        "accepted_species",
        "axis",
        "trait_composition",
        "trait_names",
        "source_groups",
        "source_lineages",
        "quality",
    }
    if missing := required.difference(frame.columns):
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


def _local_coverage_summary(frame: pd.DataFrame) -> dict[str, Any]:
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


def build_wave36_overlay(
    *, baseline_csv: Path, input_dir: Path, output_dir: Path
) -> dict[str, Any]:
    direct_path = input_dir / "wave36_resolved_direct_species_trait.csv.gz"
    low_path = input_dir / "wave36_candidate_validated_low_species_trait.csv.gz"
    rule_path = input_dir / "wave36_provider_touched_new_rule_audit.csv.gz"
    manifest_path = input_dir / "wave36_source_manifest.json"
    required = (baseline_csv, direct_path, low_path, rule_path, manifest_path)
    if missing := [str(path) for path in required if not path.is_file()]:
        raise ValueError(f"missing inputs: {missing}")

    baseline = pd.read_csv(baseline_csv, dtype=str).fillna("")
    direct = _validate_direct(pd.read_csv(direct_path, dtype=str).fillna(""))
    low = pd.read_csv(low_path, dtype=str).fillna("")
    rules = pd.read_csv(rule_path, dtype=str).fillna("")
    eligible_rules = _validate_low(low, rules)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    _local_validate_universe(baseline, label="Wave35 baseline")
    if manifest.get("contract") != "wave36_nhm_reproductive_remine_source_manifest_v1":
        raise ValueError("unexpected source manifest contract")
    hashes = manifest.get("file_sha256")
    if not isinstance(hashes, dict) or not hashes:
        raise ValueError("source manifest lacks file hashes")
    for name, expected in hashes.items():
        path = input_dir / str(name)
        if not path.is_file():
            raise ValueError(f"source manifest file is missing: {name}")
        if _sha256(path) != str(expected):
            raise ValueError(f"source manifest hash mismatch for {name}")

    before = _local_coverage_summary(baseline)
    result = baseline.copy().set_index(["accepted_species", "axis"])
    changes: list[dict[str, str]] = []
    for row in _aggregate_direct(direct).itertuples(index=False):
        key = (row.accepted_species, row.axis)
        if key not in result.index:
            raise ValueError(f"direct evidence outside denominator: {key}")
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
        changes.append(
            {
                "accepted_species": row.accepted_species,
                "axis": row.axis,
                "action": (
                    "direct_fill"
                    if before_quality == ""
                    else (
                        "direct_upgrade"
                        if QUALITY_RANK[after_quality] > QUALITY_RANK[before_quality]
                        else "direct_enrichment"
                    )
                ),
                "quality_before": before_quality,
                "quality_after": after_quality,
                "trait_names": row.trait_names,
                "source_groups": row.source_groups,
            }
        )

    for row in _aggregate_low(low).itertuples(index=False):
        key = (row.accepted_species, row.axis)
        if key not in result.index:
            raise ValueError(f"Low evidence outside denominator: {key}")
        if str(result.loc[key, "quality"]):
            raise ValueError(f"Low attempts to replace a resolved cell: {key}")
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

    result = result.reset_index().sort_values(["accepted_species", "axis"])
    _local_validate_universe(result, label="Wave36 coverage")
    after = _local_coverage_summary(result)
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
    if (
        comparison["quality_after"].map(QUALITY_RANK)
        < comparison["quality_before"].map(QUALITY_RANK)
    ).any():
        raise ValueError("an existing quality was downgraded")

    summary: dict[str, Any] = {
        "contract": "wave36_lossless_nhm_reproductive_remine_overlay_v1",
        "formal_wave33_run_id": 32932103226,
        "formal_wave35_run_id": 33072865120,
        "wave35_before": before,
        "wave36_after": after,
        "delta": {
            "gross_gain_species_axis": gain,
            "loss_species_axis": loss,
            "net_gain_species_axis": gain - loss,
            "by_axis_net_gain": {
                axis: after["by_axis"][axis]["filled_species"]
                - before["by_axis"][axis]["filled_species"]
                for axis in AXES
            },
            "action_counts": {
                str(key): int(value)
                for key, value in changes_frame["action"].value_counts().items()
            },
        },
        "validated_low_audit": {
            "eligible_provider_touched_rules": len(eligible_rules),
            "applied_species_trait": len(low),
            "strict_rebuild_is_diagnostic_not_subtractive": True,
        },
        "checks": {
            "fixed_denominator": True,
            "quality_precedence_high_medium_low": True,
            "direct_conflicts_excluded": True,
            "trait_specific_genus_join": True,
            "lineage_leave_one_out_required": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
            "baseline_loss_zero": loss == 0,
        },
        "input_hashes": {path.name: _sha256(path) for path in required},
        "source_manifest": manifest,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "wave36_species_axis_coverage.csv.gz": result,
        "wave36_direct_species_trait.csv.gz": direct,
        "wave36_validated_low_species_trait.csv.gz": low,
        "wave36_change_audit.csv.gz": changes_frame,
    }
    for name, frame in outputs.items():
        frame.to_csv(output_dir / name, index=False, compression=GZIP)
    summary["artifact_hashes"] = {name: _sha256(output_dir / name) for name in outputs}
    (output_dir / "wave36_coverage_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


@app.command("build")
def build(
    baseline_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    input_dir: Annotated[Path, typer.Option(exists=True, file_okay=False)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
) -> None:
    typer.echo(
        json.dumps(
            build_wave36_overlay(
                baseline_csv=baseline_csv,
                input_dir=input_dir,
                output_dir=output_dir,
            ),
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    app()
