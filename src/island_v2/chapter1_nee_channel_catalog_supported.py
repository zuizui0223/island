"""Canonical streaming GloBI catalog runner using provider argument semantics.

GloBI's standard ``interactions.tsv.gz`` export is generated from SUPPORTS argument
edges, while ``refuted-interactions.tsv.gz`` is generated from REFUTES argument
edges. The refuted product is therefore a contradiction-audit surface, not a bag of
rows to subtract from the supported export. This runner preserves that boundary.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import typer

from island_v2.chapter1_nee_channel_catalog import (
    _read_chunks,
    filter_interaction_chunk,
    load_config,
    sha256_file,
)
from island_v2.chapter1_nee_channel_catalog_stream import (
    CatalogAccumulator,
    verify_pinned_archives,
)

app = typer.Typer(add_completion=False)


def run_supported_catalog(
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
    policy = config["source_policy"]
    if str(policy.get("interactions_export_argument_type")) != "SUPPORTS":
        raise ValueError("canonical catalog requires GloBI interactions export = SUPPORTS")
    if str(policy.get("refuted_export_argument_type")) != "REFUTES":
        raise ValueError("canonical catalog requires GloBI refuted export = REFUTES")
    if bool(policy.get("refuted_interactions_must_be_subtracted")):
        raise ValueError("canonical provider semantics prohibit support/refute cross-subtraction")

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

    accumulator = CatalogAccumulator()
    holdout_counts: defaultdict[str, int] = defaultdict(int)
    n_input = 0
    n_retained = 0
    required = list(config["required_columns"])
    for chunk in _read_chunks(interactions_tsv_gz, required, chunksize):
        n_input += int(len(chunk))
        # GloBI interactions.tsv is already the SUPPORTS export. Refuted rows live
        # in a separate REFUTES export and are audited separately; no subtraction
        # is performed here, even if a taxon pair has contradictory arguments.
        evidence, holdouts, _ = filter_interaction_chunk(chunk, config, None)
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
        per_channel.append(
            {
                "channel_id": channel,
                "n_catalog_taxa": int(catalog["channel_id"].eq(channel).sum()) if len(catalog) else 0,
                "n_confirmatory_taxa": int(confirmatory["channel_id"].eq(channel).sum()) if len(confirmatory) else 0,
            }
        )

    receipt = {
        "contract": config["contract"],
        "status": "catalog_built_pre_N1",
        "source_semantics_revision": policy["source_semantics_revision"],
        "concept_doi": policy["concept_doi"],
        "source_version_doi": source_version_doi,
        "source_version": policy["version"],
        "source_published_date": policy["published_date"],
        "source_sha256": source_digest,
        "refuted_sha256": refuted_digest,
        "digest_check_skipped": bool(skip_digest_check),
        "interactions_export_argument_type": "SUPPORTS",
        "refuted_export_argument_type": "REFUTES",
        "refuted_archive_role": "contradiction_audit_only",
        "support_refute_cross_subtraction_applied": False,
        "n_input_interaction_rows": n_input,
        "n_refuted_flower_interactions_removed": 0,
        "n_retained_flower_interaction_rows": n_retained,
        "n_catalog_taxa": int(len(catalog)),
        "n_confirmatory_catalog_taxa": int(len(confirmatory)),
        "channels": per_channel,
        "holdout_counts": dict(sorted(holdout_counts.items())),
        "primary_catalog_rule": (
            "species-resolved explicit SUPPORTS flower interaction; >=1 pollination claim OR "
            ">=2 independent flower-visit references across >=2 plant taxa"
        ),
        "uses_focal_plant_traits": False,
        "N1_fitted": False,
    }
    (output_dir / "channel_taxon_catalog_receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
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
    receipt = run_supported_catalog(
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
