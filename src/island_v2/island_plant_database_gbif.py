from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, help="Build the GBIF candidate-flora layer for Island Plant Database 2.0.")

PAIR_COLUMNS = {
    "island_id",
    "species",
    "n_records",
    "n_unique_gbif_ids",
    "basis_of_record_set",
    "review_status",
}
TAXA_COLUMNS = {"accepted_species", "genus", "family", "n_islands", "n_records"}
BLOCK_MEMBER_COLUMNS = {"block_id", "island_id"}


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{label} missing columns: {missing}")


def _normalise_name(value: object) -> str:
    return " ".join(str(value or "").strip().split())


def taxon_id_for_name(name: str) -> str:
    normalised = _normalise_name(name)
    if not normalised:
        raise ValueError("cannot mint a taxon_id from a blank species name")
    digest = hashlib.sha256(normalised.encode("utf-8")).hexdigest()[:20]
    return f"gbif_name_{digest}"


def evidence_id_for_pair(campaign_id: str, island_id: str, taxon_id: str) -> str:
    payload = f"{campaign_id}|{island_id}|{taxon_id}".encode("utf-8")
    return f"gbif_pair_{hashlib.sha256(payload).hexdigest()[:24]}"


def _write_csv_gzip_deterministic(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            frame.to_csv(zipped, index=False)


def _campaign_downloads(campaign: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for entry in campaign.get("ledger", []):
        status = str(entry.get("request_status", ""))
        download_key = str(entry.get("download_key", "")).strip()
        doi = str(entry.get("doi", "")).strip()
        if status != "succeeded":
            continue
        if not download_key or not doi:
            raise ValueError(f"succeeded GBIF block lacks download_key/doi: {entry}")
        current = str(entry.get("block_id", "")).strip()
        migrated = str(entry.get("migrated_from_block_id", "")).strip()
        for block_id in {current, migrated} - {""}:
            rows.append(
                {
                    "block_id": block_id,
                    "download_key": download_key,
                    "doi": doi,
                }
            )
    table = pd.DataFrame(rows).drop_duplicates()
    if table.empty:
        raise ValueError("campaign contains no succeeded GBIF downloads")
    conflicting = table.groupby("block_id").size()
    if (conflicting > 1).any():
        examples = conflicting[conflicting > 1].index.tolist()[:10]
        raise ValueError(f"GBIF block IDs map to multiple downloads: {examples}")
    return table


def build_island_download_map(
    block_members: pd.DataFrame,
    campaign: dict[str, Any],
) -> pd.DataFrame:
    _require_columns(block_members, BLOCK_MEMBER_COLUMNS, "GBIF block-members table")
    if block_members["island_id"].duplicated().any():
        examples = block_members.loc[block_members["island_id"].duplicated(), "island_id"].head(10).tolist()
        raise ValueError(f"GBIF block partition contains duplicate island IDs: {examples}")
    downloads = _campaign_downloads(campaign)
    mapped = block_members[["block_id", "island_id"]].merge(
        downloads,
        on="block_id",
        how="left",
        validate="many_to_one",
    )
    if mapped[["download_key", "doi"]].isna().any(axis=None):
        missing = mapped.loc[mapped["download_key"].isna(), "block_id"].drop_duplicates().head(10).tolist()
        raise ValueError(f"block-members rows lack a succeeded campaign download mapping: {missing}")
    mapped = mapped.fillna("")
    return mapped.sort_values("island_id", kind="stable").reset_index(drop=True)


def _build_taxa(taxa_summary: pd.DataFrame, campaign_id: str) -> pd.DataFrame:
    _require_columns(taxa_summary, TAXA_COLUMNS, "GBIF taxa summary")
    if taxa_summary["accepted_species"].duplicated().any():
        raise ValueError("GBIF taxa summary contains duplicate accepted_species rows")
    table = taxa_summary.copy().fillna("")
    table["accepted_species"] = table["accepted_species"].map(_normalise_name)
    if table["accepted_species"].eq("").any():
        raise ValueError("GBIF taxa summary contains blank accepted_species")
    taxon_ids = table["accepted_species"].map(taxon_id_for_name)
    if taxon_ids.duplicated().any():
        raise ValueError("name-derived provisional taxon IDs are not unique")
    return pd.DataFrame(
        {
            "taxon_id": taxon_ids,
            "accepted_name": table["accepted_species"],
            "authorship": "",
            "taxonomic_rank": "species",
            "family": table["family"].astype(str),
            "genus": table["genus"].astype(str),
            "backbone_name": "GBIF interpreted occurrence species field",
            "backbone_key": "",
            "backbone_version": campaign_id,
            "taxonomic_status": "unresolved_taxonomy",
            "release_status": "review_required",
            "source_license": "",
        }
    ).sort_values("accepted_name", kind="stable").reset_index(drop=True)


def build_gbif_candidate_core(
    pair_path: Path,
    taxa_summary_path: Path,
    block_members_path: Path,
    campaign_path: Path,
    islands_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    pairs = pd.read_csv(pair_path, dtype=str).fillna("")
    taxa_summary = pd.read_csv(taxa_summary_path, dtype=str).fillna("")
    block_members = pd.read_csv(block_members_path, dtype=str).fillna("")
    islands = pd.read_csv(islands_path, dtype=str).fillna("")
    campaign = json.loads(campaign_path.read_text(encoding="utf-8"))
    campaign_id = str(campaign.get("campaign_id", "")).strip()
    if not campaign_id:
        raise ValueError("GBIF campaign has no campaign_id")

    _require_columns(pairs, PAIR_COLUMNS, "GBIF island-species snapshot")
    _require_columns(islands, {"island_id"}, "Database 2.0 islands")
    if pairs[["island_id", "species"]].duplicated().any():
        raise ValueError("GBIF island-species snapshot contains duplicate island × species rows")

    island_ids = set(islands["island_id"])
    pair_islands = set(pairs["island_id"])
    unknown_islands = sorted(pair_islands - island_ids)
    if unknown_islands:
        raise ValueError(f"GBIF pairs reference islands outside Database 2.0: {unknown_islands[:10]}")

    island_downloads = build_island_download_map(block_members, campaign)
    mapped_islands = set(island_downloads["island_id"])
    if mapped_islands != island_ids:
        missing = sorted(island_ids - mapped_islands)
        extra = sorted(mapped_islands - island_ids)
        raise ValueError(
            "GBIF frozen block universe differs from Database 2.0 island universe: "
            f"missing={missing[:5]} extra={extra[:5]}"
        )

    taxa = _build_taxa(taxa_summary, campaign_id)
    name_to_taxon = dict(zip(taxa["accepted_name"], taxa["taxon_id"], strict=True))
    pairs["species"] = pairs["species"].map(_normalise_name)
    missing_taxa = sorted(set(pairs["species"]) - set(name_to_taxon))
    if missing_taxa:
        raise ValueError(f"GBIF pair snapshot contains species missing from taxa summary: {missing_taxa[:10]}")
    pairs["taxon_id"] = pairs["species"].map(name_to_taxon)
    pairs = pairs.merge(
        island_downloads[["island_id", "block_id", "download_key", "doi"]],
        on="island_id",
        how="left",
        validate="many_to_one",
    )

    evidence_ids = [
        evidence_id_for_pair(campaign_id, island_id, taxon_id)
        for island_id, taxon_id in zip(pairs["island_id"], pairs["taxon_id"], strict=True)
    ]
    if len(set(evidence_ids)) != len(evidence_ids):
        raise ValueError("GBIF aggregate evidence IDs are not unique")

    n_records = pd.to_numeric(pairs["n_records"], errors="raise").astype(int)
    if (n_records < 1).any():
        raise ValueError("GBIF island-species pairs must contain at least one occurrence record")

    island_taxa = pd.DataFrame(
        {
            "island_id": pairs["island_id"],
            "taxon_id": pairs["taxon_id"],
            "membership_status": "candidate",
            "establishment_status": "unknown",
            "endemic_status": "unknown",
            "first_record_year": "",
            "last_record_year": "",
            "occurrence_record_count": n_records,
            "specimen_record_count": "",
            "evidence_count": 1,
            "review_status": "candidate_occurrence_only_unresolved_taxonomy",
            "release_status": "review_required",
            "provenance_id": evidence_ids,
        }
    ).sort_values(["island_id", "taxon_id"], kind="stable").reset_index(drop=True)

    evidence = pd.DataFrame(
        {
            "evidence_id": evidence_ids,
            "island_id": pairs["island_id"],
            "taxon_id": pairs["taxon_id"],
            "source_type": "gbif_download_aggregate",
            "source_record_id": pairs["download_key"],
            "source_url": "https://www.gbif.org/occurrence/download/" + pairs["download_key"].astype(str),
            "dataset_key_or_doi": pairs["doi"],
            "basis_of_record": pairs["basis_of_record_set"],
            "event_date": "",
            "coordinate_uncertainty_m": "",
            "evidence_role": "supports_presence",
            "source_license": "",
            "rights_status": "review_required",
            "review_status": "aggregate_pair_requires_dataset_license_audit",
        }
    ).sort_values("evidence_id", kind="stable").reset_index(drop=True)

    download_rights = island_downloads[["download_key", "doi"]].drop_duplicates().copy()
    download_rights["source_url"] = (
        "https://www.gbif.org/occurrence/download/" + download_rights["download_key"].astype(str)
    )
    rights = pd.DataFrame(
        {
            "object_type": "gbif_download",
            "object_id": download_rights["download_key"],
            "source_id": download_rights["doi"],
            "source_url": download_rights["source_url"],
            "source_license": "",
            "rights_status": "review_required",
            "rights_evidence": "https://www.gbif.org/terms",
            "redistributable_value": "false",
            "redistributable_provenance": "true",
            "review_status": "pending_constituent_dataset_license_audit",
        }
    ).sort_values("object_id", kind="stable").reset_index(drop=True)

    manifest = {
        "schema_version": 1,
        "database_id": "global_island_plant_database",
        "version": "2.0.0-alpha1",
        "build_stage": "gbif_candidate_flora_core",
        "campaign_id": campaign_id,
        "scientific_interpretation": {
            "gbif_pair_status": "candidate_not_verified_flora",
            "native_status_inferred": False,
            "endemic_status_inferred": False,
            "absence_inferred_from_missing_records": False,
        },
        "counts": {
            "islands": int(len(islands)),
            "taxa": int(len(taxa)),
            "island_taxa": int(len(island_taxa)),
            "evidence": int(len(evidence)),
            "gbif_downloads": int(len(rights)),
        },
        "rights": {
            "gbif_candidate_rows_release_status": "review_required",
            "reason": "pair-level aggregate currently lacks constituent dataset license lineage",
        },
        "inputs": {
            "islands_sha256": _sha256_file(islands_path),
            "gbif_pairs_sha256": _sha256_file(pair_path),
            "gbif_taxa_summary_sha256": _sha256_file(taxa_summary_path),
            "gbif_block_members_sha256": _sha256_file(block_members_path),
            "gbif_campaign_sha256": _sha256_file(campaign_path),
        },
    }
    return taxa, island_taxa, evidence, rights, manifest


@app.command("build")
def build_command(
    pair_path: Path = typer.Option(..., "--pairs", exists=True, dir_okay=False),
    taxa_summary_path: Path = typer.Option(..., "--taxa-summary", exists=True, dir_okay=False),
    block_members_path: Path = typer.Option(..., "--block-members", exists=True, dir_okay=False),
    campaign_path: Path = typer.Option(..., "--campaign", exists=True, dir_okay=False),
    islands_path: Path = typer.Option(..., "--islands", exists=True, dir_okay=False),
    output_dir: Path = typer.Option(..., "--output-dir"),
) -> None:
    taxa, island_taxa, evidence, rights, manifest = build_gbif_candidate_core(
        pair_path,
        taxa_summary_path,
        block_members_path,
        campaign_path,
        islands_path,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv_gzip_deterministic(taxa, output_dir / "taxa.csv.gz")
    _write_csv_gzip_deterministic(island_taxa, output_dir / "island_taxa.csv.gz")
    _write_csv_gzip_deterministic(evidence, output_dir / "evidence.csv.gz")
    rights.to_csv(output_dir / "RIGHTS_LEDGER.csv", index=False)
    (output_dir / "DATABASE_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    typer.echo(
        f"Built Database 2.0 alpha1 candidate flora: {len(taxa):,} taxa, "
        f"{len(island_taxa):,} island-taxon pairs"
    )


if __name__ == "__main__":
    app()
