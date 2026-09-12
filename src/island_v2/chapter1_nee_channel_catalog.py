"""Build the outcome-independent pollinator taxon catalog used by N1.

Primary use is a pinned, versioned GloBI integrated interaction archive. The live
API is intentionally not used here: the catalog must be reproducible and frozen
before island-retention outcomes are inspected.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

EVIDENCE_COLUMNS = [
    "channel_id",
    "pollinator_species",
    "pollinator_genus",
    "pollinator_family",
    "pollinator_order",
    "pollinator_class",
    "plant_taxon",
    "interaction_type_id",
    "interaction_type_name",
    "evidence_strength",
    "reference_key",
    "reference_citation",
    "reference_doi",
    "reference_url",
    "source_namespace",
    "source_archive_uri",
    "source_doi",
]
CATALOG_COLUMNS = [
    "channel_id",
    "pollinator_species",
    "pollinator_genus",
    "pollinator_family",
    "pollinator_order",
    "pollinator_class",
    "n_interaction_rows",
    "n_pollination_rows",
    "n_flower_visit_rows",
    "n_independent_references",
    "n_distinct_plant_taxa",
    "strongest_evidence",
    "catalog_tier",
    "reference_examples",
    "source_version_doi",
    "source_sha256",
]
HOLDOUT_COLUMNS = ["reason", "interaction_type_id", "pollinator_name", "plant_name"]


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise typer.BadParameter("channel catalog config must be a mapping")
    if config.get("contract") != "chapter1_nee_channel_taxon_catalog_v1":
        raise typer.BadParameter("unexpected channel catalog contract")
    return config


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _row_value(row: dict[str, object], prefix: str, suffix: str) -> str:
    return _text(row.get(f"{prefix}Taxon{suffix}"))


def _reference_key(row: dict[str, object]) -> str:
    doi = _text(row.get("referenceDoi"))
    if doi:
        return "doi:" + doi.lower()
    citation = _text(row.get("referenceCitation"))
    if citation:
        return "citation:" + citation
    url = _text(row.get("referenceUrl"))
    if url:
        return "url:" + url
    source_doi = _text(row.get("sourceDOI"))
    if source_doi:
        return "source_doi:" + source_doi.lower()
    archive = _text(row.get("sourceArchiveURI"))
    if archive:
        return "archive:" + archive
    namespace = _text(row.get("sourceNamespace"))
    source_citation = _text(row.get("sourceCitation"))
    if namespace or source_citation:
        return "source:" + namespace + "|" + source_citation
    return ""


def interaction_fingerprint(row: dict[str, object]) -> str:
    """Stable key used only to subtract the separately published refuted archive."""
    fields = [
        _text(row.get("sourceTaxonName")),
        _text(row.get("interactionTypeId")),
        _text(row.get("targetTaxonName")),
        _reference_key(row),
    ]
    return hashlib.sha256("\x1f".join(fields).encode("utf-8")).hexdigest()


def _channel_for_pollinator(
    *, genus: str, family: str, order: str, taxon_class: str, config: dict[str, Any]
) -> str:
    rules = config["channel_rules"]
    if genus == str(rules["bombus"]["genus_exact"]):
        return "bombus"
    bee = rules["non_bombus_bees"]
    if family in set(bee["families"]) and genus not in set(bee["exclude_genera"]):
        return "non_bombus_bees"
    if order == str(rules["lepidoptera"]["order_exact"]):
        return "lepidoptera"
    if taxon_class == str(rules["flower_visiting_birds"]["class_exact"]):
        return "flower_visiting_birds"
    if order == str(rules["diptera"]["order_exact"]):
        return "diptera"
    return ""


def _extract_record(row: dict[str, object], config: dict[str, Any]) -> tuple[dict[str, str] | None, dict[str, str] | None]:
    interaction_id = _text(row.get("interactionTypeId"))
    interaction = config["accepted_interactions"].get(interaction_id)
    if not interaction:
        return None, None

    pollinator_prefix = "source" if interaction["pollinator_side"] == "source" else "target"
    plant_prefix = "source" if interaction["plant_side"] == "source" else "target"
    plant_kingdom = _row_value(row, plant_prefix, "KingdomName")
    plant_name = _row_value(row, plant_prefix, "SpeciesName") or _row_value(row, plant_prefix, "Name")
    pollinator_name = _row_value(row, pollinator_prefix, "SpeciesName") or _row_value(
        row, pollinator_prefix, "Name"
    )

    if plant_kingdom != str(config["plant_side_requirement"]["kingdom"]):
        return None, {
            "reason": "plant_side_not_resolved_as_Plantae",
            "interaction_type_id": interaction_id,
            "pollinator_name": pollinator_name,
            "plant_name": plant_name,
        }

    species = _row_value(row, pollinator_prefix, "SpeciesName")
    genus = _row_value(row, pollinator_prefix, "GenusName")
    family = _row_value(row, pollinator_prefix, "FamilyName")
    order = _row_value(row, pollinator_prefix, "OrderName")
    taxon_class = _row_value(row, pollinator_prefix, "ClassName")
    channel = _channel_for_pollinator(
        genus=genus, family=family, order=order, taxon_class=taxon_class, config=config
    )
    if not channel:
        return None, None
    if not species:
        return None, {
            "reason": "pollinator_species_unresolved",
            "interaction_type_id": interaction_id,
            "pollinator_name": pollinator_name,
            "plant_name": plant_name,
        }
    reference_key = _reference_key(row)
    if not reference_key:
        return None, {
            "reason": "interaction_reference_unresolved",
            "interaction_type_id": interaction_id,
            "pollinator_name": species,
            "plant_name": plant_name,
        }
    return {
        "channel_id": channel,
        "pollinator_species": species,
        "pollinator_genus": genus,
        "pollinator_family": family,
        "pollinator_order": order,
        "pollinator_class": taxon_class,
        "plant_taxon": plant_name,
        "interaction_type_id": interaction_id,
        "interaction_type_name": _text(row.get("interactionTypeName")),
        "evidence_strength": str(interaction["strength"]),
        "reference_key": reference_key,
        "reference_citation": _text(row.get("referenceCitation")),
        "reference_doi": _text(row.get("referenceDoi")),
        "reference_url": _text(row.get("referenceUrl")),
        "source_namespace": _text(row.get("sourceNamespace")),
        "source_archive_uri": _text(row.get("sourceArchiveURI")),
        "source_doi": _text(row.get("sourceDOI")),
    }, None


def filter_interaction_chunk(
    chunk: pd.DataFrame,
    config: dict[str, Any],
    refuted_fingerprints: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    """Keep explicit flower interactions for N0 channels and return holdouts."""
    missing = set(config["required_columns"]).difference(chunk.columns)
    if missing:
        raise ValueError(f"GloBI interaction table missing required columns: {sorted(missing)}")
    refuted = refuted_fingerprints or set()
    evidence_rows: list[dict[str, str]] = []
    holdout_rows: list[dict[str, str]] = []
    n_refuted = 0
    accepted_ids = set(config["accepted_interactions"])
    subset = chunk.loc[chunk["interactionTypeId"].fillna("").astype(str).isin(accepted_ids)]
    for row in subset.fillna("").to_dict("records"):
        if interaction_fingerprint(row) in refuted:
            n_refuted += 1
            continue
        evidence, holdout = _extract_record(row, config)
        if evidence is not None:
            evidence_rows.append(evidence)
        if holdout is not None:
            holdout_rows.append(holdout)
    return (
        pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS),
        pd.DataFrame(holdout_rows, columns=HOLDOUT_COLUMNS),
        n_refuted,
    )


def build_catalog(
    evidence: pd.DataFrame,
    config: dict[str, Any],
    *,
    source_version_doi: str,
    source_sha256: str,
) -> pd.DataFrame:
    """Aggregate explicit interactions to a species-level frozen channel catalog."""
    missing = set(EVIDENCE_COLUMNS).difference(evidence.columns)
    if missing:
        raise ValueError(f"channel interaction evidence missing columns: {sorted(missing)}")
    thresholds = config["qualification"]
    rows: list[dict[str, Any]] = []
    for (channel, species), group in evidence.groupby(
        ["channel_id", "pollinator_species"], sort=True
    ):
        strengths = group["evidence_strength"].astype(str)
        references = sorted({value for value in group["reference_key"].astype(str) if value})
        plants = sorted({value for value in group["plant_taxon"].astype(str) if value})
        n_pollination = int(strengths.eq("pollination").sum())
        n_visits = int(strengths.eq("flower_visit").sum())
        confirmatory = n_pollination >= 1 or (
            len(references)
            >= int(thresholds["confirmatory_catalog_taxon"]["or_min_independent_flower_visit_references"])
            and len(plants)
            >= int(thresholds["confirmatory_catalog_taxon"]["min_distinct_plant_taxa_for_visit_route"])
        )
        tier = "confirmatory" if confirmatory else "sensitivity"
        first = group.iloc[0]
        rows.append(
            {
                "channel_id": channel,
                "pollinator_species": species,
                "pollinator_genus": first["pollinator_genus"],
                "pollinator_family": first["pollinator_family"],
                "pollinator_order": first["pollinator_order"],
                "pollinator_class": first["pollinator_class"],
                "n_interaction_rows": int(len(group)),
                "n_pollination_rows": n_pollination,
                "n_flower_visit_rows": n_visits,
                "n_independent_references": len(references),
                "n_distinct_plant_taxa": len(plants),
                "strongest_evidence": "pollination" if n_pollination else "flower_visit",
                "catalog_tier": tier,
                "reference_examples": " || ".join(references[:5]),
                "source_version_doi": source_version_doi,
                "source_sha256": source_sha256,
            }
        )
    return pd.DataFrame(rows, columns=CATALOG_COLUMNS)


def _read_chunks(path: Path, columns: list[str], chunksize: int) -> Any:
    return pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        dtype=str,
        usecols=columns,
        chunksize=chunksize,
        low_memory=False,
    )


def _refuted_fingerprints(path: Path, config: dict[str, Any], chunksize: int) -> set[str]:
    fingerprints: set[str] = set()
    columns = list(config["required_columns"])
    for chunk in _read_chunks(path, columns, chunksize):
        accepted_ids = set(config["accepted_interactions"])
        subset = chunk.loc[chunk["interactionTypeId"].fillna("").astype(str).isin(accepted_ids)]
        for row in subset.fillna("").to_dict("records"):
            fingerprints.add(interaction_fingerprint(row))
    return fingerprints


def _append_gzip_csv(path: Path, table: pd.DataFrame, *, header: bool) -> None:
    if table.empty:
        return
    mode = "wt" if header else "at"
    with gzip.open(path, mode, encoding="utf-8", newline="") as handle:
        table.to_csv(handle, index=False, header=header)


@app.command("build")
def build_command(
    interactions_tsv_gz: Path = typer.Option(..., exists=True, dir_okay=False),
    refuted_tsv_gz: Path = typer.Option(..., exists=True, dir_okay=False),
    source_version_doi: str = typer.Option(..., help="Resolved version DOI, not only the concept DOI."),
    output_dir: Path = typer.Option(...),
    chunksize: int = typer.Option(200_000, min=1),
    config_path: Path = typer.Option(Path("config/chapter1_nee_channel_taxon_catalog.yml")),
) -> None:
    """Build a reproducible species-level channel catalog from a pinned GloBI archive."""
    config = load_config(config_path)
    concept_doi = str(config["source_policy"]["concept_doi"])
    if source_version_doi.strip().lower() == concept_doi.lower():
        raise typer.BadParameter("source_version_doi must be a resolved version DOI, not the concept DOI")
    source_digest = sha256_file(interactions_tsv_gz)
    refuted_digest = sha256_file(refuted_tsv_gz)
    refuted = _refuted_fingerprints(refuted_tsv_gz, config, chunksize)

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "channel_flower_interaction_evidence.csv.gz"
    if evidence_path.exists():
        evidence_path.unlink()
    evidence_frames: list[pd.DataFrame] = []
    holdout_counts: defaultdict[str, int] = defaultdict(int)
    n_input = 0
    n_refuted = 0
    wrote_header = False
    columns = list(config["required_columns"])
    for chunk in _read_chunks(interactions_tsv_gz, columns, chunksize):
        n_input += int(len(chunk))
        evidence, holdouts, removed = filter_interaction_chunk(chunk, config, refuted)
        n_refuted += removed
        if not evidence.empty:
            _append_gzip_csv(evidence_path, evidence, header=not wrote_header)
            wrote_header = True
            evidence_frames.append(evidence)
        if not holdouts.empty:
            for reason, count in holdouts["reason"].value_counts().items():
                holdout_counts[str(reason)] += int(count)

    all_evidence = (
        pd.concat(evidence_frames, ignore_index=True)
        if evidence_frames
        else pd.DataFrame(columns=EVIDENCE_COLUMNS)
    )
    catalog = build_catalog(
        all_evidence,
        config,
        source_version_doi=source_version_doi,
        source_sha256=source_digest,
    )
    catalog.to_csv(output_dir / "channel_taxon_catalog.csv", index=False)
    catalog.loc[catalog["catalog_tier"].eq("confirmatory")].to_csv(
        output_dir / "channel_taxon_catalog_confirmatory.csv", index=False
    )
    summary = {
        "contract": config["contract"],
        "concept_doi": concept_doi,
        "source_version_doi": source_version_doi,
        "source_sha256": source_digest,
        "refuted_sha256": refuted_digest,
        "n_input_interaction_rows": n_input,
        "n_refuted_flower_interactions_removed": n_refuted,
        "n_retained_flower_interaction_rows": int(len(all_evidence)),
        "n_catalog_taxa": int(len(catalog)),
        "n_confirmatory_catalog_taxa": int(catalog["catalog_tier"].eq("confirmatory").sum())
        if len(catalog)
        else 0,
        "channel_counts_confirmatory": {
            str(key): int(value)
            for key, value in catalog.loc[
                catalog["catalog_tier"].eq("confirmatory"), "channel_id"
            ].value_counts().sort_index().items()
        }
        if len(catalog)
        else {},
        "holdout_counts": dict(sorted(holdout_counts.items())),
        "primary_rule": (
            "species-resolved explicit GloBI flower interactions only; pollination claim OR "
            ">=2 independent flower-visit references across >=2 plant taxa"
        ),
    }
    (output_dir / "channel_taxon_catalog_receipt.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    typer.echo(json.dumps(summary, indent=2))


if __name__ == "__main__":
    app()
