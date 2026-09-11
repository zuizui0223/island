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
    """Collapse row-level provenance IDs to a source family used for rights review.

    Exact lineage strings are retained in a separate audit table. The family is only
    the unit on which a common redistribution decision can be applied safely.
    """
    token = lineage.strip()
    lower = token.lower()

    if lower.startswith("dataset:dryad.") or lower.startswith("doi:10.5061/dryad."):
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


def load_lineage_family_map(path: Path | None) -> dict[str, str]:
    """Load an audit-only exact-lineage -> source-family override map.

    This map changes only the rights-review grouping. It never mutates the scientific
    Database 1.0 ledger or its stored source_lineages.
    """
    if path is None:
        return {}
    frame = pd.read_csv(path, dtype=str).fillna("")
    required = {"source_lineage", "source_family"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"lineage-family map missing columns: {sorted(missing)}")
    if frame["source_lineage"].eq("").any() or frame["source_family"].eq("").any():
        raise ValueError("lineage-family map contains blank source_lineage/source_family")
    conflicts = frame.groupby("source_lineage")["source_family"].nunique()
    conflicts = conflicts[conflicts > 1]
    if not conflicts.empty:
        raise ValueError(f"lineage-family map has {len(conflicts)} conflicting lineages")
    deduped = frame.drop_duplicates(["source_lineage", "source_family"])
    return dict(zip(deduped["source_lineage"], deduped["source_family"], strict=True))


def family_for_lineage(lineage: str, overrides: dict[str, str] | None = None) -> str:
    if overrides and lineage in overrides:
        return overrides[lineage]
    return source_family(lineage)


def audit_bucket(source: str) -> str:
    """Return a non-binding triage bucket; never use this field to grant rights."""
    prefix = source.split(":", 1)[0].lower() if source else "unresolved"
    known = {
        "database",
        "dataset",
        "derived",
        "domain",
        "origin",
        "provider",
        "provider_compilation",
        "provider_treatment",
        "unresolved",
    }
    if prefix in known:
        return prefix
    return f"prefix:{prefix or 'empty'}"


def source_decision(source: str, policy: dict[str, object]) -> tuple[str, str, str]:
    default_status = str(policy.get("default_status", "review_required"))
    default_license = policy.get("default_license")
    default_note = str(policy.get("default_note", "No redistribution decision recorded."))
    for rule in policy.get("rules", []):
        pattern = str(rule["pattern"])
        if re.search(pattern, source, flags=re.IGNORECASE):
            license_id = rule.get("license")
            return (
                str(rule["status"]),
                "" if license_id is None else str(license_id),
                str(rule.get("note", "")),
            )
    return (
        default_status,
        "" if default_license is None else str(default_license),
        default_note,
    )


def validate_public_release_license(policy: dict[str, object]) -> str:
    """Require a release licence compatible with any admitted ShareAlike source.

    The source policy remains authoritative per source family. The public release
    licence is an umbrella for this compilation and never replaces upstream terms.
    """
    release_license = str(policy.get("public_release_license", "")).strip()
    if not release_license:
        raise ValueError("source policy must declare public_release_license")
    wiki_status, wiki_license, _ = source_decision("domain:en.wikipedia.org", policy)
    if (
        wiki_status == "redistributable"
        and wiki_license == "CC-BY-SA-4.0"
        and release_license != "CC-BY-SA-4.0"
    ):
        raise ValueError(
            "Wikipedia redistribution requires public_release_license=CC-BY-SA-4.0"
        )
    return release_license


def build_inventory(
    frame: pd.DataFrame,
    policy: dict[str, object],
    lineage_family_overrides: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return source-family rights inventory and exact-lineage provenance details."""
    overrides = lineage_family_overrides or {}
    lineage_counts = tokenize_lineages(frame["source_lineages"])
    family_lineages: dict[str, set[str]] = defaultdict(set)
    family_lineage_mentions: Counter[str] = Counter()
    for lineage, count in lineage_counts.items():
        family = family_for_lineage(lineage, overrides)
        family_lineages[family].add(lineage)
        family_lineage_mentions[family] += count

    family_cell_mentions: Counter[str] = Counter()
    for value in frame["source_lineages"].fillna("").astype(str):
        families = {
            family_for_lineage(token.strip(), overrides)
            for token in value.split("|")
            if token.strip()
        }
        family_cell_mentions.update(families)

    family_rows: list[dict[str, object]] = []
    for family, cell_count in family_cell_mentions.most_common():
        status, license_id, note = source_decision(family, policy)
        examples = sorted(family_lineages[family])
        family_rows.append(
            {
                "source_family": family,
                "audit_bucket": audit_bucket(family),
                "resolved_cell_mentions": cell_count,
                "lineage_mentions": family_lineage_mentions[family],
                "distinct_source_lineages": len(examples),
                "example_lineage": examples[0] if examples else "",
                "redistribution_status": status,
                "source_license": license_id,
                "policy_note": note,
            }
        )
    family_inventory = pd.DataFrame(family_rows)

    lineage_rows: list[dict[str, object]] = []
    for lineage, count in lineage_counts.most_common():
        family = family_for_lineage(lineage, overrides)
        status, license_id, note = source_decision(family, policy)
        lineage_rows.append(
            {
                "source_lineage": lineage,
                "source_family": family,
                "audit_bucket": audit_bucket(family),
                "lineage_mentions": count,
                "family_override_applied": lineage in overrides,
                "redistribution_status": status,
                "source_license": license_id,
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
    lineage_family_map: Path | None = typer.Option(
        None,
        "--lineage-family-map",
        exists=True,
        dir_okay=False,
        help="Audit-only exact source_lineage -> normalized source_family override CSV.",
    ),
    public: bool = typer.Option(
        False,
        "--public",
        help="Require every source family to have explicit redistribution permission and license.",
    ),
) -> None:
    manifest = load_manifest(manifest_path)
    policy = yaml.safe_load(source_policy_path.read_text(encoding="utf-8")) or {}
    release_license = validate_public_release_license(policy)
    overrides = load_lineage_family_map(lineage_family_map)
    species_axis = source_dir / "species_axis_coverage.csv.gz"
    validation = validate_species_axis(species_axis, manifest)

    frame = pd.read_csv(species_axis, dtype=str).fillna("")
    resolved = frame[frame["quality"].str.lower().isin({"high", "medium", "low"})].copy()
    inventory, lineage_details = build_inventory(resolved, policy, overrides)
    blockers = inventory[
        (inventory["redistribution_status"] != "redistributable")
        | inventory["source_license"].eq("")
    ].copy()

    output_dir.mkdir(parents=True, exist_ok=True)
    inventory.to_csv(output_dir / "SOURCE_LICENSE_INVENTORY.csv", index=False)
    lineage_details.to_csv(output_dir / "SOURCE_LINEAGE_DETAILS.csv", index=False)
    blockers.to_csv(output_dir / "RELEASE_BLOCKERS.csv", index=False)

    bucket_summary = (
        blockers.groupby("audit_bucket", dropna=False)
        .agg(
            blocker_families=("source_family", "nunique"),
            resolved_cell_mentions=("resolved_cell_mentions", "sum"),
            exact_lineages=("distinct_source_lineages", "sum"),
        )
        .reset_index()
        .sort_values("resolved_cell_mentions", ascending=False)
    )
    bucket_summary.to_csv(output_dir / "RIGHTS_TRIAGE_BUCKETS.csv", index=False)

    release_ready = blockers.empty
    if public and not release_ready:
        raise typer.BadParameter(
            f"public release blocked: {len(blockers)} source families require rights review"
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

    licenses_observed = sorted(
        set(inventory.loc[inventory["source_license"].ne(""), "source_license"].astype(str))
    )
    applied_overrides = int(lineage_details["family_override_applied"].sum())
    release_manifest = {
        "schema_version": 3,
        "database_id": manifest.database.database_id,
        "database_version": manifest.database.version,
        "analysis_contract": manifest.database.analysis_contract,
        "source_workflow_run_id": manifest.database.source.run_id,
        "source_artifact_name": manifest.database.source.artifact_name,
        "source_species_axis_sha256": manifest.database.source.sha256,
        "validation": validation,
        "release_ready_for_public_zenodo": release_ready,
        "public_mode_requested": public,
        "public_release_license": release_license,
        "files": copied,
        "license_blocker_count": int(len(blockers)),
        "distinct_source_families": int(len(inventory)),
        "distinct_source_lineages": int(len(lineage_details)),
        "lineage_family_override_rows": int(len(overrides)),
        "lineage_family_overrides_applied": applied_overrides,
        "source_licenses_with_explicit_permission": licenses_observed,
    }
    (output_dir / "RELEASE_MANIFEST.json").write_text(
        json.dumps(release_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    license_notice = (
        "# License and attribution notice\n\n"
        f"The public Chapter 1 Database 1.0 compilation is distributed under **{release_license}** "
        "to the extent copyright or database rights apply to the compilation.\n\n"
        "This umbrella licence does **not** replace, broaden, or relicense upstream source terms. "
        "Source-specific licences, rights decisions, and provenance are retained in "
        "`SOURCE_LICENSE_INVENTORY.csv` and `SOURCE_LINEAGE_DETAILS.csv`; downstream users must "
        "preserve the applicable source attribution and licence obligations.\n\n"
        "Wikipedia-derived rows in Database 1.0 contain normalized trait facts and exact article-URL "
        "provenance, not copied article prose. Their source family is recorded as CC BY-SA 4.0 and "
        "the compilation-level CC BY-SA 4.0 licence preserves the ShareAlike boundary.\n\n"
        "Rows whose source family remains `review_required` are not authorized for public release "
        "until that source-specific rights gate is closed.\n"
    )
    (output_dir / "LICENSE_NOTICE.md").write_text(license_notice, encoding="utf-8")

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
        "`SOURCE_LINEAGE_DETAILS.csv` retains the exact lineage tokens behind those families. "
        "When supplied, a lineage-family map changes only audit grouping and is recorded by "
        "`family_override_applied`; it never changes Database 1.0 source lineage bytes. "
        "`RIGHTS_TRIAGE_BUCKETS.csv` is only a workload summary and never grants redistribution rights. "
        "`LICENSE_NOTICE.md` records the compilation-level public licence while preserving all "
        "source-specific licence obligations.\n"
    )
    (output_dir / "DATA_DICTIONARY.md").write_text(data_dictionary, encoding="utf-8")

    status = "READY" if release_ready else "BLOCKED_PENDING_SOURCE_LICENSE_REVIEW"
    typer.echo(f"{status}: {output_dir}")


if __name__ == "__main__":
    app()
