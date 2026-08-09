"""Build uncertainty-aware genus predictions without changing strict coverage.

The strict ledger continues to use the shared current-min3 Validated Low rules.
This module materializes two separate analysis inputs from the same trait-specific
rule audit:

* relaxed-min3 predictions for secondary multiple-imputation analyses; and
* min-species-2 predictions for sensitivity diagnostics only.

No family rule, global fallback, or genus-by-axis join is performed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.all_evidence_trait_audit import apply_genus_rules


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _pairs(frame: pd.DataFrame) -> set[tuple[str, str]]:
    return set(zip(frame["accepted_species"], frame["trait_name"], strict=True))


def _eligible_rule_metrics(rules: pd.DataFrame, setting: str) -> dict[str, Any]:
    selected = rules.loc[
        rules["setting"].eq(setting) & rules["eligible"].astype(bool)
    ].copy()

    def summarize(group: pd.DataFrame) -> dict[str, Any]:
        species_n = group["species_loo_n"].sum()
        lineage_n = group["lineage_loo_n"].sum()
        return {
            "eligible_rules": len(group),
            "median_support_species": float(group["n_direct_species"].median()),
            "mean_dominance": float(group["dominance"].mean()),
            "mean_entropy": float(group["entropy"].mean()),
            "species_loo_accuracy": (
                float(group["species_loo_correct"].sum() / species_n)
                if species_n
                else None
            ),
            "lineage_loo_accuracy": (
                float(group["lineage_loo_correct"].sum() / lineage_n)
                if lineage_n
                else None
            ),
        }

    return {
        **summarize(selected),
        "by_trait": {
            str(trait): summarize(group)
            for trait, group in selected.groupby("trait_name", sort=True)
        },
    }


def _incremental_predictions(
    predictions: pd.DataFrame,
    baseline_pairs: set[tuple[str, str]],
    strict_filled_axes: set[tuple[str, str]],
    master: pd.DataFrame,
    *,
    role: str,
) -> pd.DataFrame:
    incremental = predictions.loc[
        [pair not in baseline_pairs for pair in _row_pairs(predictions)]
    ].copy()
    incremental["empirical_modal_probability"] = pd.to_numeric(
        incremental["dominance"], errors="coerce"
    )
    incremental["potential_new_species_axis"] = [
        (species, axis) not in strict_filled_axes
        for species, axis in zip(
            incremental["accepted_species"], incremental["axis"], strict=True
        )
    ]
    islands = master[["accepted_species", "n_islands"]].drop_duplicates(
        "accepted_species"
    )
    incremental = incremental.merge(
        islands, on="accepted_species", how="left", validate="many_to_one"
    )
    incremental["analysis_role"] = role
    incremental["quality"] = "probabilistic_low"
    incremental["strict_coverage_eligible"] = False
    incremental["family_inference_used"] = False
    incremental["global_fallback_used"] = False
    incremental["join_key"] = "genus_x_trait_name"
    return incremental.sort_values(
        [
            "potential_new_species_axis",
            "n_islands",
            "empirical_modal_probability",
            "accepted_species",
            "trait_name",
        ],
        ascending=[False, False, False, True, True],
    ).reset_index(drop=True)


def _row_pairs(frame: pd.DataFrame) -> list[tuple[str, str]]:
    return list(zip(frame["accepted_species"], frame["trait_name"], strict=True))


def _prediction_metrics(frame: pd.DataFrame) -> dict[str, Any]:
    potential = frame.loc[frame["potential_new_species_axis"].astype(bool)]
    return {
        "additional_species_trait": len(
            frame.drop_duplicates(["accepted_species", "trait_name"])
        ),
        "additional_species": frame["accepted_species"].nunique(),
        "potential_new_species_axis": len(
            potential.drop_duplicates(["accepted_species", "axis"])
        ),
        "by_axis": {
            str(axis): len(group.drop_duplicates(["accepted_species", "axis"]))
            for axis, group in potential.groupby("axis", sort=True)
        },
        "by_trait": {
            str(trait): len(group.drop_duplicates(["accepted_species", "trait_name"]))
            for trait, group in frame.groupby("trait_name", sort=True)
        },
    }


def build(
    *,
    master_csv: Path,
    direct_csv: Path,
    rules_csv: Path,
    strict_coverage_csv: Path,
    strict_low_csv: Path,
    output_dir: Path,
) -> dict[str, Any]:
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    direct = pd.read_csv(direct_csv, dtype=str).fillna("")
    rules = pd.read_csv(rules_csv, low_memory=False)
    strict_coverage = pd.read_csv(strict_coverage_csv, dtype=str).fillna("")
    strict_low = pd.read_csv(strict_low_csv, dtype=str).fillna("")

    if len(strict_coverage) != 318_885:
        raise ValueError("strict coverage must contain exactly 318885 species-axis rows")
    if strict_coverage.duplicated(["accepted_species", "axis"]).any():
        raise ValueError("strict coverage has duplicate species-axis rows")
    strict_species = set(strict_coverage["accepted_species"])
    if len(strict_species) != 106_295:
        raise ValueError("probabilistic Low requires the fixed 106295-species universe")
    master = master.loc[master["accepted_species"].isin(strict_species)].copy()
    if master["accepted_species"].nunique() != 106_295:
        raise ValueError("master does not cover every species in strict coverage")

    strict_pairs = _pairs(strict_low)
    strict_filled_axes = set(
        zip(
            strict_coverage.loc[
                strict_coverage["quality"].ne(""), "accepted_species"
            ],
            strict_coverage.loc[strict_coverage["quality"].ne(""), "axis"],
            strict=True,
        )
    )
    relaxed = apply_genus_rules(master, direct, rules, "relaxed_min3")
    min2 = apply_genus_rules(master, direct, rules, "relaxed_min2_diagnostic")
    relaxed_incremental = _incremental_predictions(
        relaxed,
        strict_pairs,
        strict_filled_axes,
        master,
        role="secondary_multiple_imputation_relaxed_min3",
    )
    min2_incremental = _incremental_predictions(
        min2,
        strict_pairs,
        strict_filled_axes,
        master,
        role="sensitivity_only_min_species_2",
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    relaxed_path = output_dir / "secondary_probabilistic_genus_low.csv.gz"
    min2_path = output_dir / "sensitivity_min2_genus_low.csv.gz"
    relaxed_incremental.to_csv(relaxed_path, index=False)
    min2_incremental.to_csv(min2_path, index=False)

    strict_filled = int(strict_coverage["quality"].ne("").sum())
    summary: dict[str, Any] = {
        "contract": "secondary_probabilistic_genus_low_v1",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "denominator": {"species": 106_295, "species_axis": 318_885},
        "strict_coverage": {
            "filled_species_axis": strict_filled,
            "unresolved_species_axis": 318_885 - strict_filled,
            "changed_by_this_artifact": False,
            "strict_validated_low_species_trait": len(strict_pairs),
        },
        "relaxed_min3_secondary": {
            **_prediction_metrics(relaxed_incremental),
            "rule_validation": _eligible_rule_metrics(rules, "relaxed_min3"),
            "analysis_role": "secondary robustness with species-level multiple imputation",
        },
        "min2_diagnostic": {
            **_prediction_metrics(min2_incremental),
            "rule_validation": _eligible_rule_metrics(
                rules, "relaxed_min2_diagnostic"
            ),
            "analysis_role": "sensitivity only; not recommended for main analysis",
        },
        "guardrails": {
            "trait_specific_join": "genus_x_trait_name",
            "axis_only_join": False,
            "family_inference": False,
            "global_fallback": False,
            "direct_cells_overridden": False,
            "strict_coverage_modified": False,
            "same_species_draw_shared_across_islands": True,
        },
    }
    summary_path = output_dir / "secondary_probabilistic_genus_low_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    inputs = [master_csv, direct_csv, rules_csv, strict_coverage_csv, strict_low_csv]
    outputs = [relaxed_path, min2_path, summary_path]
    manifest = {
        "contract": "secondary_probabilistic_genus_low_manifest_v1",
        "inputs": {
            str(path): {"sha256": _sha256(path), "bytes": path.stat().st_size}
            for path in inputs
        },
        "outputs": {
            path.name: {"sha256": _sha256(path), "bytes": path.stat().st_size}
            for path in outputs
        },
    }
    (output_dir / "secondary_probabilistic_genus_low_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--direct-csv", type=Path, required=True)
    parser.add_argument("--rules-csv", type=Path, required=True)
    parser.add_argument("--strict-coverage-csv", type=Path, required=True)
    parser.add_argument("--strict-low-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(**vars(args)), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
