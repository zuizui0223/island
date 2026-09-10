from __future__ import annotations

import hashlib
import json
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd
import typer
import yaml

from island_v2.chapter1_database_manifest import load_manifest, validate_species_axis

app = typer.Typer(help="Build a fail-closed Zenodo candidate for a Chapter 1 database snapshot.")

PUBLIC_FILES = (
    "species_axis_coverage.csv.gz",
    "direct_species_trait_ledger.csv.gz",
    "integration_summary.json",
    "baseline_recovery_manifest.json",
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def tokenize_lineages(series: pd.Series) -> Counter[str]:
    counter: Counter[str] = Counter()
    for value in series.fillna("").astype(str):
        for token in value.split("|"):
            token = token.strip()
            if token:
                counter[token] += 1
    return counter


def source_family(lineage: str) -> str:
    """Collapse row-level provenance IDs to a licensable provider/source family.

    The goal is not to erase provenance: the raw lineage is retained separately. The
    family is only the unit on which redistribution policy is evaluated.
    """
    token = lineage.strip()
    lower = token.lower()

    if lower.startswith("dataset:dryad."):
        return "dataset:dryad"
    if lower.startswith("doi:10.5061/dryad."):
        return "dataset:dryad"
    if lower.startswith("origin:austraits:"):
        return "dataset:austraits"
    if lower.startswith("pladias:"):
        return "database:pladias"
    if lower.startswith("origin:usda_plants"):
        return "database:usda_plants"
    if lower.startswith("bhl:"):
        return "provider:bhl"
    if lower.startswith("baseflor:"):
        return "database:baseflor"
    if lower.startswith("ecoflora:"):
        return "database:ecoflora"
    if lower.startswith("efloras:") or lower.startswith("efloras-treatment:"):
        return "provider:efloras"
    if lower.startswith("floraweb:") or lower.startswith("biolflor:floraweb"):
        return "database:floraweb_biolflor"
    if lower.startswith("citation:"):
        return "unresolved:citation_hash"
    if lower.startswith("validated-low:"):
        return "derived:validated_low_without_direct_source_lineage"
    if lower.startswith("url:"):
        host = urlparse(token[4:]).netloc.lower().removeprefix("www.")
        return f"domain:{host or 'unknown'}"
    if lower.startswith("provider_treatment:"):
        parts = token.split(":")
        return ":".join(parts[:2]).lower()
    if lower.startswith("provider_compilation:"):
        parts = token.split(":")
        return ":".join(parts[:2]).lower()
    if lower.startswith("dataset:"):
        parts = token.split(":")
        return ":".join(parts[:2]).lower()
    if lower.startswith("origin:"):
        parts = token.split(":")
        return ":".join(parts[:2]).lower()
    if lower.startswith("doi:"):
        return "unresolved:other_doi"

    parts = token.split(":")
    if len(parts) >= 2:
        return ":".join(parts[:2]).lower()
    return f"unresolved:{lower or 'empty'}"


def source_status(source: str, policy: dict[str, object]) -> tuple[str, str]:
    default = str(policy.get("default_status", "review_required"))
    default_note = str(policy.get("default_note", "No redistribution decision recorded."))
    for rule in policy.get("rules", []):
        pattern = str(rule["pattern"])
        if re.search(pattern, source, flags=re.IGNORECASE):
            return str(rule["status"]), str(rule.get("note", ""))
    return default, default_note


def build_inventory(
    frame: pd.DataFrame,
    policy: dict[str, object],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return family-level policy inventory plus raw-lineage provenance details."""
    lineage_counts = tokenize_lineages(frame["source_lineages"])
    family_lineages: dict[str, set[str]] = defaultdict(set)
    family_lineage_mentions: Counter[str] = Counter()
    for lineage, count in lineage_counts.items():
        family = source_family(lineage)
        family_lineages[family].add(lineage)
        family_lineage_mentions[family] += count

    family_cell_mentions: Counter[str] = Counter()
    for value in frame["source_lineages"].fillna("").astype(str):
        families = {
            source_family(token.strip())
            for token in value.split("|")
            if token.strip()
        }
        family_cell_mentions.update(families)

    family_rows: list[dict[str, object]] = []
    for family, cell_count in family_cell_mentions.most_common():
        status, note = source_status(family, policy)
        examples = sorted(family_lineages[family])
        family_rows.append(
            {
                "source_family": family,
                "resolved_cell_mentions": cell_count,
                "lineage_mentions": family_lineage_mentions[family],
                "distinct_source_lineages": len(examples),
                "example_lineage": examples[0] if examples else "",
                "redistribution_status": status,
                "policy_note": note,
            }
        )
    family_inventory = pd.DataFrame(family_rows)

    lineage_rows: list[dict[str, object]] = []
    for lineage, count in lineage_counts.most_common():
        family = source_family(lineage)
        status, note = source_status(family, policy)
        lineage_rows.append(
            {
                "source_lineage": lineage,
                "source_family": family,
                "lineage_mentions": count,
                "redistribution_status": status,
                "policy_note": note,
            }
        )
    lineage_inventory = pd.DataFrame(lineage_rows)
    return family_inventory, lineage_inventory


@app.command("build")
def build(
    source_dir: Path = typer.Option(..., "--source-dir", exists=True, file_okay=False),
    manifest_path: Path = typer.Option(..., "--manifest", exists=True, dir_okay=False),
    source_policy_path: Path = typer.Option(..., "--source-policy", exists=True, dir_okay=False),
    output_dir: Path = typer.Option(..., "--output-dir"),
    public: bool = typer.Option(
        False,
        "--public",
        help="Require every source family to be explicitly redistributable before copying data files.",
    ),
) -> None:
    manifest = load_manifest(manifest_path)
    policy = yaml.safe_load(source_policy_path.read_text(encoding="utf-8")) or {}
    species_axis = source_dir / "species_axis_coverage.csv.gz"
    validation = validate_species_axis(species_axis, manifest)

    frame = pd.read_csv(species_axis, dtype=str).fillna("")
    resolved = frame[frame["quality"].str.lower().isin({"high", "medium", "low"})].copy()
    inventory, lineage_details = build_inventory(resolved, policy)
    blockers = inventory[inventory["redistribution_status"] != "redistributable"].copy()

    output_dir.mkdir(parents=True, exist_ok=True)
    inventory.to_csv(output_dir / "SOURCE_LICENSE_INVENTORY.csv", index=False)
    lineage_details.to_csv(output_dir / "SOURCE_LINEAGE_DETAILS.csv", index=False)
    blockers.to_csv(output_dir / "RELEASE_BLOCKERS.csv", index=False)

    release_ready = blockers.empty
    if public and not release_ready:
        raise typer.BadParameter(
            f"public release blocked: {len(blockers)} source families require redistribution review"
        )

    copied: list[dict[str, object]] = []
    if not public or release_ready:
        for name in PUBLIC_FILES:
            src = source_dir / name
            if not src.exists():
                raise FileNotFoundError(f"required release file missing: {src}")
            dst = output_dir / name
            shutil.copy2(src, dst)
            copied.append({"file": name, "sha256": sha256_file(dst), "bytes": dst.stat().st_size})

    release_manifest = {
        "schema_version": 1,
        "database_id": manifest.database.database_id,
        "database_version": manifest.database.version,
        "analysis_contract": manifest.database.analysis_contract,
        "source_workflow_run_id": manifest.database.source.run_id,
        "source_artifact_name": manifest.database.source.artifact_name,
        "source_species_axis_sha256": manifest.database.source.sha256,
        "validation": validation,
        "release_ready_for_public_zenodo": release_ready,
        "public_mode_requested": public,
        "files": copied,
        "license_blocker_count": int(len(blockers)),
        "distinct_source_families": int(len(inventory)),
        "distinct_source_lineages": int(len(lineage_details)),
    }
    (output_dir / "RELEASE_MANIFEST.json").write_text(
        json.dumps(release_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    data_dictionary = "# Chapter 1 database data dictionary\n\n"
    data_dictionary += (
        "`species_axis_coverage.csv.gz` contains one row per accepted species × analysis axis.\n\n"
    )
    data_dictionary += "| column | meaning |\n| --- | --- |\n"
    data_dictionary += "| accepted_species | accepted analysis taxon name |\n"
    data_dictionary += (
        "| axis | one of flower_colour, floral_structural_complexity, reproductive_assurance |\n"
    )
    data_dictionary += (
        "| trait_composition | normalized composition/value payload used to build analysis traits |\n"
    )
    data_dictionary += "| trait_names | contributing normalized trait names |\n"
    data_dictionary += "| source_groups | acquisition/integration route labels; not a licensing unit |\n"
    data_dictionary += "| source_lineages | row-level provenance retained through integration |\n"
    data_dictionary += "| quality | high, medium, low, or unresolved/blank |\n"
    data_dictionary += (
        "\n`SOURCE_LICENSE_INVENTORY.csv` evaluates normalized source families; "
        "`SOURCE_LINEAGE_DETAILS.csv` retains the exact lineage tokens behind those families.\n"
    )
    (output_dir / "DATA_DICTIONARY.md").write_text(data_dictionary, encoding="utf-8")

    status = "READY" if release_ready else "BLOCKED_PENDING_SOURCE_LICENSE_REVIEW"
    typer.echo(f"{status}: {output_dir}")


if __name__ == "__main__":
    app()
