"""Build a conservative multi-source genus completion layer for plant SC/SI.

Every configured species-direct source is combined before a genus can qualify.
A single mixed/variable species, a within-species source conflict, or both SI and
SC within the genus disqualifies transfer.  Qualifying targets are emitted as
low-confidence hierarchical candidates and never as species-direct evidence.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
import typer
import yaml

from island_v2.angiosperm_scope import classify_scope, load_config as load_scope_config
from island_v2.meyer_2026_reproductive_traits import (
    CANDIDATE_COLUMNS,
    _text,
    _write_csv_gzip,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.callback()
def main() -> None:
    """Build labelled likely_SC/likely_SI genus candidates from all direct sources."""


def load_config(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise typer.BadParameter(f"Expected a mapping in {path}")
    return value


def load_direct_candidates(paths: Iterable[Path], trait_name: str) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    required = {
        "accepted_species",
        "genus",
        "family",
        "trait_name",
        "standardized_value",
        "candidate_kind",
        "evidence_scope",
        "source_name",
    }
    for path in paths:
        frame = pd.read_csv(path, dtype=str).fillna("")
        if missing := required.difference(frame.columns):
            raise typer.BadParameter(f"{path} lacks candidate columns: {sorted(missing)}")
        frame = frame.loc[
            frame["trait_name"].eq(trait_name)
            & frame["candidate_kind"].eq("source_backed")
            & frame["evidence_scope"].eq("species_direct")
        ].copy()
        frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=sorted(required))
    return pd.concat(frames, ignore_index=True).drop_duplicates()


def _species_summary(direct: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    explicit = set(config["explicit_values"])
    disqualifying = set(config.get("disqualifying_values") or [])
    rows: list[dict[str, Any]] = []
    for (species, genus, family), group in direct.groupby(
        ["accepted_species", "genus", "family"], sort=True
    ):
        values = sorted({_text(value) for value in group["standardized_value"] if _text(value)})
        explicit_values = sorted(set(values).intersection(explicit))
        ambiguous = bool(set(values).intersection(disqualifying)) or len(explicit_values) != 1
        rows.append(
            {
                "accepted_species": species,
                "genus": genus,
                "family": family,
                "state": "" if ambiguous else explicit_values[0],
                "ambiguous": ambiguous,
                "reported_values": json.dumps(values, ensure_ascii=False),
                "source_names": sorted(set(group["source_name"].map(_text))),
            }
        )
    return pd.DataFrame(rows)


def _qualify_genus(group: pd.DataFrame, config: dict[str, Any]) -> tuple[str, int] | None:
    minimum = int(config["minimum_direct_species"])
    if len(group) < minimum or group["ambiguous"].astype(bool).any():
        return None
    states = group["state"].value_counts()
    if len(states) != 1:
        return None
    purity = float(states.iloc[0] / states.sum())
    if purity < float(config["minimum_modal_purity"]):
        return None
    return _text(states.index[0]), int(states.sum())


def build_combined_genus_inferences(
    master: pd.DataFrame,
    direct: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    species_summary = _species_summary(direct, config)
    qualifying_rows: list[dict[str, Any]] = []
    for genus, group in species_summary.groupby("genus", sort=True):
        qualified = _qualify_genus(group, config)
        if qualified is None:
            continue
        value, support = qualified
        qualifying_rows.append(
            {
                "genus": genus,
                "standardized_value": value,
                "support_n": support,
                "supporting_taxa": sorted(group["accepted_species"].tolist()),
                "source_names": sorted(
                    {source for sources in group["source_names"] for source in sources}
                ),
            }
        )
    qualifying = pd.DataFrame(qualifying_rows)
    lookup = {
        row["genus"]: row for row in qualifying.to_dict("records")
    } if not qualifying.empty else {}

    loo_total = 0
    loo_correct = 0
    for _genus, group in species_summary.groupby("genus", sort=False):
        for held in group.itertuples(index=False):
            if not _text(held.state):
                continue
            remaining = group.loc[group["accepted_species"].ne(held.accepted_species)]
            qualified = _qualify_genus(remaining, config)
            if qualified is None:
                continue
            loo_total += 1
            loo_correct += qualified[0] == _text(held.state)

    directly_reported = set(species_summary["accepted_species"])
    rows: list[dict[str, Any]] = []
    for target in master.itertuples(index=False):
        species = _text(target.accepted_species)
        genus = _text(target.genus)
        if species in directly_reported or genus not in lookup:
            continue
        support = lookup[genus]
        value = _text(support["standardized_value"])
        support_n = int(support["support_n"])
        sources = support["source_names"]
        rows.append(
            {
                "evidence_id": f"combined_si_genus_v1:{genus}:{species}",
                "source_name": "combined_species_direct_si_sc_v1",
                "source_record_id": f"genus:{genus}:{value}:{species}",
                "accepted_species": species,
                "genus": genus,
                "family": _text(target.family),
                "name_match_method": "accepted_name_target_for_genus_inference",
                "trait_name": config["trait_name"],
                "raw_value": value,
                "standardized_value": value,
                "candidate_kind": "hierarchical_inference",
                "evidence_scope": "genus_inference",
                "source_type": "review",
                "source_reliability": "B_curated_database_or_institution",
                "source_citation": "Combined species-direct sources: " + ", ".join(sources),
                "source_url": "",
                "evidence_excerpt": (
                    f"All {support_n} directly reported {genus} species across "
                    f"{', '.join(sources)} have state {value}; the target itself is unreported."
                ),
                "supporting_taxa": json.dumps(support["supporting_taxa"], ensure_ascii=False),
                "inference_rule": (
                    f"Multi-source unanimous genus transfer; support_n={support_n}, "
                    f"minimum={config['minimum_direct_species']}, "
                    f"purity={config['minimum_modal_purity']}."
                ),
                "confidence": config["confidence"],
                "inference_note": (
                    "Completion/sensitivity only. Any direct mixed state or source conflict "
                    "disqualifies its genus; any target species record supersedes this row."
                ),
                "needs_human_review": True,
                "review_status": "pending",
                "analysis_role": config["analysis_role"],
            }
        )
    inference = pd.DataFrame(rows, columns=CANDIDATE_COLUMNS)
    report = {
        "n_direct_candidate_rows": int(len(direct)),
        "n_direct_species_any_state": int(len(species_summary)),
        "n_ambiguous_or_mixed_direct_species": int(species_summary["ambiguous"].sum()),
        "n_qualifying_genera": len(qualifying),
        "n_inferred_species": len(inference),
        "leave_one_out_n": loo_total,
        "leave_one_out_correct": loo_correct,
        "leave_one_out_accuracy": loo_correct / loo_total if loo_total else None,
    }
    return inference, qualifying, report


def run_inference(*, config_path: Path, output_dir: Path) -> dict[str, Any]:
    config = load_config(config_path)
    direct_paths = [Path(path) for path in config["direct_candidate_csvs"]]
    direct = load_direct_candidates(direct_paths, config["trait_name"])
    master = pd.read_csv(config["master_csv"], dtype=str).fillna("")
    scope = classify_scope(master, load_scope_config(Path(config["scope_config"])))
    eligible = set(
        scope.loc[scope["angiosperm_analysis_eligible"], "accepted_species"].map(_text)
    )
    master = master.loc[master["accepted_species"].isin(eligible)].drop_duplicates(
        "accepted_species"
    )
    inference, qualifying, report = build_combined_genus_inferences(master, direct, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    candidate_path = output_dir / "genus_inference_candidates.csv.gz"
    genus_path = output_dir / "qualifying_genera.csv"
    _write_csv_gzip(inference, candidate_path)
    if qualifying.empty:
        qualifying = pd.DataFrame(
            columns=[
                "genus",
                "standardized_value",
                "support_n",
                "supporting_taxa",
                "source_names",
            ]
        )
    qualifying.assign(
        supporting_taxa=qualifying["supporting_taxa"].map(
            lambda value: json.dumps(value, ensure_ascii=False)
        ),
        source_names=qualifying["source_names"].map(
            lambda value: json.dumps(value, ensure_ascii=False)
        ),
    ).to_csv(genus_path, index=False)
    manifest = {
        "version": config["version"],
        "config_path": str(config_path),
        "direct_candidate_csvs": [str(path) for path in direct_paths],
        **report,
        "outputs": {
            "genus_inference_candidates": str(candidate_path),
            "qualifying_genera": str(genus_path),
        },
        "policy": (
            "Low-confidence multi-source genus completion only. Mixed/variable species, "
            "within-species conflicts and mixed genera disqualify transfer."
        ),
    }
    manifest_path = output_dir / "inference_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    manifest["manifest"] = str(manifest_path)
    return manifest


@app.command("build")
def build(
    config_path: Path = typer.Option(Path("config/si_genus_inference.yml"), exists=True),
    output_dir: Path = typer.Option(
        Path("data/v2/staging/traits/inference/si_genus_v1")
    ),
) -> None:
    report = run_inference(config_path=config_path, output_dir=output_dir)
    typer.echo(
        f"Built {report['n_inferred_species']:,} low-confidence genus SI/SC candidates "
        f"from {report['n_qualifying_genera']:,} unanimous genera; "
        f"LOO accuracy={report['leave_one_out_accuracy']:.3f}."
    )


if __name__ == "__main__":
    app()
