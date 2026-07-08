"""Build a stratified gold-standard validation pilot and evaluate its metrics.

Two commands:

- ``build`` selects a small pilot (target_min..target_max) from a candidate
  species table by spreading selections **evenly across biological strata**
  (native/introduced x growth form x family x flower conspicuousness) — never by
  GBIF record count, which would bias toward human-associated / cultivated /
  easily observed species and hide the failure modes we want to measure.
- ``evaluate`` computes the pilot's review metrics, including
  ``irrelevant_literature_rate`` as a **primary** metric.

Fail-closed: the builder only nominates species for human review and decides no
trait value; the evaluator counts only explicit judgements — a blank is never a
judgement, and a rate with a zero denominator is reported as ``null`` (undefined),
never fabricated.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.callback()
def main() -> None:
    """Build a stratified validation pilot and evaluate its metrics."""


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or "stratification_dimensions" not in config:
        raise typer.BadParameter("pilot config must contain stratification_dimensions")
    return config


def _stratum_key(row: pd.Series, columns: list[str]) -> tuple[str, ...]:
    return tuple(str(row.get(c, "")).strip().lower() or "unknown" for c in columns)


def build_stratified_pilot(
    candidates: pd.DataFrame, config: dict[str, Any]
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Select up to target_max species spread evenly across occupied strata.

    Deterministic: within each stratum species are taken in sorted order, and
    strata are visited round-robin (largest strata first, ties broken by key) so
    the sample stays reproducible and balanced rather than record-count driven.
    """
    species_col = next(
        (c for c in ("accepted_species", "species", "scientific_name") if c in candidates.columns),
        None,
    )
    if species_col is None:
        raise typer.BadParameter("candidate table must have an 'accepted_species' column")

    columns = [
        dim["column"]
        for dim in config["stratification_dimensions"].values()
        if dim.get("column") in candidates.columns
    ]
    if not columns:
        raise typer.BadParameter("none of the configured stratification columns are present")

    target_max = int(config.get("target_max", 50))
    frame = candidates.copy()
    frame["_stratum"] = [tuple(_stratum_key(row, columns)) for _, row in frame.iterrows()]

    # Bucket sorted species indices per stratum.
    buckets: dict[tuple[str, ...], list[int]] = {}
    for stratum, group in frame.groupby("_stratum", sort=True):
        ordered = group.sort_values(species_col).index.tolist()
        buckets[stratum] = ordered

    # Round-robin across strata, largest first, so no single stratum dominates.
    order = sorted(buckets, key=lambda k: (-len(buckets[k]), k))
    selected: list[int] = []
    exhausted = False
    while len(selected) < target_max and not exhausted:
        exhausted = True
        for stratum in order:
            if buckets[stratum]:
                selected.append(buckets[stratum].pop(0))
                exhausted = False
                if len(selected) >= target_max:
                    break

    pilot = frame.loc[selected].drop(columns=["_stratum"]).reset_index(drop=True)

    target_min = int(config.get("target_min", 30))
    strata_counts = {
        "|".join(k): int(sum(1 for i in selected if frame.loc[i, "_stratum"] == k))
        for k in order
    }
    summary = {
        "note": (
            "Species nominated for the validation pilot only; no trait value is "
            "decided. Selection is stratified across biological axes, never by GBIF "
            "record count."
        ),
        "n_candidates": int(len(candidates)),
        "n_strata_occupied": int(len(buckets)),
        "n_selected": int(len(pilot)),
        "target_min": target_min,
        "target_max": target_max,
        "meets_target_min": bool(len(pilot) >= min(target_min, len(candidates))),
        "stratification_columns": columns,
        "selected_per_stratum": {k: v for k, v in strata_counts.items() if v},
    }
    return pilot, summary


def _rate(numerator: int, denominator: int) -> float | None:
    """Return numerator/denominator, or None when the denominator is zero."""
    return round(numerator / denominator, 4) if denominator else None


def evaluate_pilot(review: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    """Compute pilot review metrics from explicit judgement columns only.

    Every rate counts only non-blank judgements in the relevant column; a rate
    with no judgements is ``None`` (undefined). ``irrelevant_literature_rate`` is a
    primary metric, not a secondary agreement statistic.
    """

    def _col(name: str) -> pd.Series:
        if name in review.columns:
            return review[name].astype(str).str.strip().str.lower()
        return pd.Series([""] * len(review), dtype=str)

    literature = _col("literature_relevant")
    n_lit = int(literature.isin({"relevant", "irrelevant"}).sum())
    n_irrelevant = int((literature == "irrelevant").sum())

    scope = _col("evidence_scope")
    n_scope = int((scope != "").sum())
    n_species_direct = int((scope == "species_direct").sum())

    src_to_trait = _col("source_to_trait_relevant")
    n_src = int(src_to_trait.isin({"yes", "no"}).sum())
    n_src_relevant = int((src_to_trait == "yes").sum())

    resolved = _col("trait_resolved")
    n_resolved_judged = int(resolved.isin({"resolved", "unresolved"}).sum())
    n_unresolved = int((resolved == "unresolved").sum())

    transfer = _col("hierarchical_transfer_correct")
    n_transfer = int(transfer.isin({"correct", "wrong"}).sum())
    n_wrong_transfer = int((transfer == "wrong").sum())

    metrics = {
        # PRIMARY: how much off-target literature the retrieval pulled in.
        "irrelevant_literature_rate": {
            "rate": _rate(n_irrelevant, n_lit),
            "numerator": n_irrelevant,
            "denominator": n_lit,
        },
        "species_direct_evidence_proportion": {
            "rate": _rate(n_species_direct, n_scope),
            "numerator": n_species_direct,
            "denominator": n_scope,
        },
        "source_to_trait_relevance_rate": {
            "rate": _rate(n_src_relevant, n_src),
            "numerator": n_src_relevant,
            "denominator": n_src,
        },
        "trait_specific_unresolved_rate": {
            "rate": _rate(n_unresolved, n_resolved_judged),
            "numerator": n_unresolved,
            "denominator": n_resolved_judged,
        },
        "wrong_hierarchical_transfer_rate": {
            "rate": _rate(n_wrong_transfer, n_transfer),
            "numerator": n_wrong_transfer,
            "denominator": n_transfer,
        },
    }
    return {
        "note": (
            "Rates count only explicit judgements; a blank is never a judgement and "
            "a zero-denominator rate is null (undefined), never fabricated."
        ),
        "n_review_rows": int(len(review)),
        "primary_metrics": list(config.get("primary_metrics", [])),
        "metrics": metrics,
    }


@app.command("build")
def build(
    candidates_csv: Path = typer.Option(..., exists=True, help="Candidate species with stratification columns."),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/pilot_stratification.yml")),
) -> None:
    """Select a stratified validation pilot (for human review) and write a summary."""
    config = load_config(config_path)
    candidates = pd.read_csv(candidates_csv, dtype=str).fillna("")
    pilot, summary = build_stratified_pilot(candidates, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    pilot.to_csv(output_dir / "validation_pilot_species.csv", index=False)
    (output_dir / "validation_pilot_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    typer.echo(
        f"Selected {summary['n_selected']} pilot species across "
        f"{summary['n_strata_occupied']} strata (target {summary['target_min']}-{summary['target_max']})."
    )


@app.command("evaluate")
def evaluate(
    review_csv: Path = typer.Option(..., exists=True, help="Human-review record table for the pilot."),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/pilot_stratification.yml")),
) -> None:
    """Compute the pilot review metrics (irrelevant_literature_rate is primary)."""
    config = load_config(config_path)
    review = pd.read_csv(review_csv, dtype=str).fillna("")
    report = evaluate_pilot(review, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "validation_pilot_metrics.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    irr = report["metrics"]["irrelevant_literature_rate"]["rate"]
    typer.echo(
        f"Evaluated {report['n_review_rows']} review row(s); "
        f"irrelevant_literature_rate={irr}."
    )


if __name__ == "__main__":
    app()
