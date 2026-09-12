"""Memory-bounded runner for the frozen N1 GloBI channel catalog.

The scientific rules live in ``chapter1_nee_channel_catalog``. This module only
adds archive-pin verification and chunk-wise aggregation so the multi-GB GloBI
0.9 archive can be processed on a standard GitHub Actions runner without keeping
all retained interaction rows in memory.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2.chapter1_nee_channel_catalog import (
    CATALOG_COLUMNS,
    _read_chunks,
    _refuted_fingerprints,
    filter_interaction_chunk,
    load_config,
    sha256_file,
)

app = typer.Typer(add_completion=False)


class CatalogAccumulator:
    """Aggregate species-level catalog statistics without retaining raw rows."""

    def __init__(self) -> None:
        self._stats: dict[tuple[str, str], dict[str, Any]] = {}

    def add(self, evidence: pd.DataFrame) -> None:
        if evidence.empty:
            return
        for row in evidence.to_dict("records"):
            key = (str(row["channel_id"]), str(row["pollinator_species"]))
            if key not in self._stats:
                self._stats[key] = {
                    "pollinator_genus": str(row["pollinator_genus"]),
                    "pollinator_family": str(row["pollinator_family"]),
                    "pollinator_order": str(row["pollinator_order"]),
                    "pollinator_class": str(row["pollinator_class"]),
                    "n_interaction_rows": 0,
                    "n_pollination_rows": 0,
                    "n_flower_visit_rows": 0,
                    "references": set(),
                    "plants": set(),
                }
            stats = self._stats[key]
            stats["n_interaction_rows"] += 1
            if str(row["evidence_strength"]) == "pollination":
                stats["n_pollination_rows"] += 1
            elif str(row["evidence_strength"]) == "flower_visit":
                stats["n_flower_visit_rows"] += 1
            reference = str(row["reference_key"])
            if reference:
                stats["references"].add(reference)
            plant = str(row["plant_taxon"])
            if plant:
                stats["plants"].add(plant)

    def to_frame(
        self,
        config: dict[str, Any],
        *,
        source_version_doi: str,
        source_sha256: str,
    ) -> pd.DataFrame:
        threshold = config["qualification"]["confirmatory_catalog_taxon"]
        rows: list[dict[str, Any]] = []
        for (channel, species), stats in sorted(self._stats.items()):
            references = sorted(stats["references"])
            plants = sorted(stats["plants"])
            n_pollination = int(stats["n_pollination_rows"])
            confirmatory = n_pollination >= 1 or (
                len(references)
                >= int(threshold["or_min_independent_flower_visit_references"])
                and len(plants) >= int(threshold["min_distinct_plant_taxa_for_visit_route"])
            )
            rows.append(
                {
                    "channel_id": channel,
                    "pollinator_species": species,
                    "pollinator_genus": stats["pollinator_genus"],
                    "pollinator_family": stats["pollinator_family"],
                    "pollinator_order": stats["pollinator_order"],
                    "pollinator_class": stats["pollinator_class"],
                    "n_interaction_rows": int(stats["n_interaction_rows"]),
                    "n_pollination_rows": n_pollination,
                    "n_flower_visit_rows": int(stats["n_flower_visit_rows"]),
                    "n_independent_references": len(references),
                    "n_distinct_plant_taxa": len(plants),
                    "strongest_evidence": "pollination" if n_pollination else "flower_visit",
                    "catalog_tier": "confirmatory" if confirmatory else "sensitivity",
                    "reference_examples": " || ".join(references[:5]),
                    "source_version_doi": source_version_doi,
                    "source_sha256": source_sha256,
                }
            )
        return pd.DataFrame(rows, columns=CATALOG_COLUMNS)


def verify_pinned_archives(
    interactions_tsv_gz: Path,
    refuted_tsv_gz: Path,
    config: dict[str, Any],
    source_version_doi: str,
) -> tuple[str, str]:
    """Hard-stop if the supplied archives differ from the pre-outcome GloBI pin."""
    policy = config["source_policy"]
    expected_doi = str(policy["resolved_version_doi"])
    if source_version_doi.strip().lower() != expected_doi.lower():
        raise ValueError(
            f"source version DOI mismatch: expected {expected_doi}, got {source_version_doi}"
        )
    source_digest = sha256_file(interactions_tsv_gz)
    refuted_digest = sha256_file(refuted_tsv_gz)
    if source_digest != str(policy["interactions_tsv_sha256"]):
        raise ValueError("interactions.tsv.gz SHA-256 does not match frozen GloBI 0.9 pin")
    if refuted_digest != str(policy["refuted_interactions_tsv_sha256"]):
        raise ValueError("refuted-interactions.tsv.gz SHA-256 does not match frozen GloBI 0.9 pin")
    return source_digest, refuted_digest


def run_streaming_catalog(
    *,
    interactions_tsv_gz: Path,
    refuted_tsv_gz: Path,
    source_version_doi: str,
    output_dir: Path,
    chunksize: int,
    config_path: Path,
    skip_digest_check: bool = False,
) -> dict[str, Any]:
    config = load_config(config_path)
    if skip_digest_check:
        source_digest = sha256_file(interactions_tsv_gz)
        refuted_digest = sha256_file(refuted_tsv_gz)
    else:
        source_digest, refuted_digest = verify_pinned_archives(
            interactions_tsv_gz,
            refuted_tsv_gz,
            config,
            source_version_doi,
        )

    refuted = _refuted_fingerprints(refuted_tsv_gz, config, chunksize)
    accumulator = CatalogAccumulator()
    holdout_counts: defaultdict[str, int] = defaultdict(int)
    n_input = 0
    n_retained = 0
    n_refuted = 0
    required = list(config["required_columns"])
    for chunk in _read_chunks(interactions_tsv_gz, required, chunksize):
        n_input += int(len(chunk))
        evidence, holdouts, removed = filter_interaction_chunk(chunk, config, refuted)
        n_refuted += int(removed)
        n_retained += int(len(evidence))
        accumulator.add(evidence)
        if not holdouts.empty:
            for reason, count in holdouts["reason"].value_counts().items():
                holdout_counts[str(reason)] += int(count)

    catalog = accumulator.to_frame(
        config,
        source_version_doi=source_version_doi,
        source_sha256=source_digest,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    catalog.to_csv(output_dir / "channel_taxon_catalog.csv", index=False)
    confirmatory = catalog.loc[catalog["catalog_tier"].eq("confirmatory")].copy()
    confirmatory.to_csv(output_dir / "channel_taxon_catalog_confirmatory.csv", index=False)

    per_channel = []
    for channel in config["channel_rules"]:
        all_n = int(catalog["channel_id"].eq(channel).sum()) if len(catalog) else 0
        conf_n = int(confirmatory["channel_id"].eq(channel).sum()) if len(confirmatory) else 0
        per_channel.append(
            {
                "channel_id": channel,
                "n_catalog_taxa": all_n,
                "n_confirmatory_taxa": conf_n,
            }
        )

    receipt = {
        "contract": config["contract"],
        "status": "catalog_built_pre_N1",
        "concept_doi": config["source_policy"]["concept_doi"],
        "source_version_doi": source_version_doi,
        "source_version": config["source_policy"]["version"],
        "source_published_date": config["source_policy"]["published_date"],
        "source_sha256": source_digest,
        "refuted_sha256": refuted_digest,
        "digest_check_skipped": bool(skip_digest_check),
        "n_input_interaction_rows": n_input,
        "n_refuted_flower_interactions_removed": n_refuted,
        "n_retained_flower_interaction_rows": n_retained,
        "n_catalog_taxa": int(len(catalog)),
        "n_confirmatory_catalog_taxa": int(len(confirmatory)),
        "channels": per_channel,
        "holdout_counts": dict(sorted(holdout_counts.items())),
        "primary_catalog_rule": (
            "species-resolved explicit flower interaction; >=1 pollination claim OR "
            ">=2 independent flower-visit references across >=2 plant taxa"
        ),
        "uses_focal_plant_traits": False,
        "N1_fitted": False,
    }
    (output_dir / "channel_taxon_catalog_receipt.json").write_text(
        json.dumps(receipt, indent=2), encoding="utf-8"
    )
    return receipt


@app.command("build")
def build_command(
    interactions_tsv_gz: Path = typer.Option(..., exists=True, dir_okay=False),
    refuted_tsv_gz: Path = typer.Option(..., exists=True, dir_okay=False),
    output_dir: Path = typer.Option(...),
    source_version_doi: str = typer.Option("10.5281/zenodo.20546682"),
    chunksize: int = typer.Option(200_000, min=1),
    config_path: Path = typer.Option(Path("config/chapter1_nee_channel_taxon_catalog.yml")),
) -> None:
    receipt = run_streaming_catalog(
        interactions_tsv_gz=interactions_tsv_gz,
        refuted_tsv_gz=refuted_tsv_gz,
        source_version_doi=source_version_doi,
        output_dir=output_dir,
        chunksize=chunksize,
        config_path=config_path,
    )
    typer.echo(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    app()
