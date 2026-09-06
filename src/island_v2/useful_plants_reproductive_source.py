"""Promote audited Ken Fern/PFAF self-fertility pages as Medium evidence.

The three websites redistribute closely related Useful Plants content.  Every
row for the same species therefore shares one source lineage.  The structured
``Self-fertile: Yes`` field can establish self compatibility for a wild,
sexually reproducing species, but it never establishes autonomous selfing.
Apomictic/asexual, cultivar-only, identity-mismatched, and family-mismatched
pages are rejected upstream and checked again here.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS

app = typer.Typer(add_completion=False, no_args_is_help=True)

ACCEPTED_KEN_STATUS = "accepted_wild_self_fertile"
ACCEPTED_PFAF_STATUS = "accepted_self_fertile_statement"
PROHIBITED_RE = r"(?i)\b(?:apomictic|agamosperm\w*|asexual|cultivar[- ]only)\b"

AUDIT_COLUMNS = [
    "accepted_species",
    "provider",
    "source_url",
    "supporting_excerpt",
    "page_sha256",
    "content_fingerprint",
    "decision",
    "reason",
    "source_lineage",
]


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _id(*parts: object) -> str:
    payload = "\n".join(_text(part) for part in parts)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]


def _bool(value: object) -> bool:
    return _text(value).casefold() in {"true", "1", "yes"}


def _lineage(species: str) -> str:
    return f"ken-fern-useful-plants:species:{species.casefold()}"


def _candidate_rows(frame: pd.DataFrame, provider: str) -> list[dict[str, Any]]:
    required = {
        "accepted_species",
        "source_url",
        "supporting_excerpt",
        "normalized_value",
        "trait_name",
        "status",
        "page_sha256",
        "content_fingerprint",
    }
    if missing := required.difference(frame.columns):
        raise ValueError(f"{provider} candidate file missing columns: {sorted(missing)}")
    rows: list[dict[str, Any]] = []
    for record in frame.fillna("").to_dict("records"):
        species = _text(record["accepted_species"])
        status = _text(record["status"])
        expected = ACCEPTED_KEN_STATUS if provider == "ken_fern" else ACCEPTED_PFAF_STATUS
        reason = "selected"
        if status != expected:
            reason = "not_upstream_accepted"
        elif _text(record["trait_name"]) != "self_incompatibility":
            reason = "trait_category_mismatch"
        elif _text(record["normalized_value"]) != "SC":
            reason = "value_category_mismatch"
        elif not _text(record["supporting_excerpt"]):
            reason = "missing_exact_supporting_excerpt"
        elif pd.Series([_text(record["supporting_excerpt"])]).str.contains(
            PROHIBITED_RE, regex=True
        ).iloc[0]:
            reason = "nonsexual_or_cultivar_category"
        elif provider == "ken_fern" and _bool(record.get("nonsexual_qualified")):
            reason = "nonsexual_reproduction_category"
        elif provider == "pfaf" and any(
            _bool(record.get(column))
            for column in (
                "cultivar_qualified",
                "nonsexual_qualified",
                "dioecious_conflict",
            )
        ):
            reason = "upstream_category_exclusion"
        rows.append(
            {
                **record,
                "provider": provider,
                "source_lineage": _lineage(species),
                "decision": "accept" if reason == "selected" else "exclude",
                "reason": reason,
            }
        )
    return rows


def build_package(
    ken_candidates: pd.DataFrame,
    pfaf_candidates: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return Medium direct evidence and the provider/lineage audit."""

    rows = _candidate_rows(ken_candidates, "ken_fern") + _candidate_rows(
        pfaf_candidates, "pfaf"
    )
    candidates = pd.DataFrame(rows).fillna("")
    audit = candidates.reindex(columns=AUDIT_COLUMNS).copy()
    accepted = candidates.loc[candidates["decision"].eq("accept")].copy()
    conflicts = (
        accepted.groupby("accepted_species")["normalized_value"].nunique().gt(1)
    )
    conflict_species = set(conflicts.loc[conflicts].index)
    if conflict_species:
        accepted = accepted.loc[
            ~accepted["accepted_species"].isin(conflict_species)
        ].copy()
        mask = audit["accepted_species"].isin(conflict_species)
        audit.loc[mask, "decision"] = "exclude"
        audit.loc[mask, "reason"] = "within_lineage_direct_conflict"

    evidence_rows: list[dict[str, str]] = []
    for record in accepted.to_dict("records"):
        species = _text(record["accepted_species"])
        provider = _text(record["provider"])
        provider_name = (
            "Useful Tropical/Temperate Plants"
            if provider == "ken_fern"
            else "Plants For A Future"
        )
        evidence_rows.append(
            {
                "accepted_species": species,
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "normalized_value": "SC",
                "quality": "medium",
                "source_group": "useful_plants_public_web",
                "source_provider": provider_name,
                "source_url": _text(record["source_url"]),
                "source_record_id": f"{provider}:" + _id(
                    species, record["source_url"], record["page_sha256"]
                ),
                "source_citation": (
                    f"{provider_name} species page, accessed "
                    f"{_text(record.get('retrieved_at_utc'))[:10]}."
                ),
                "source_excerpt": _text(record["supporting_excerpt"]),
                "evidence_scope": "species_direct",
                "name_match_method": "accepted_name_exact",
                "source_lineage": _lineage(species),
                "lineage_method": "ken_fern_provider_redistribution_collapse",
                "source_run_id": "local-wave35-useful-plants-20260827",
                "source_artifact": (
                    "ken_fern_self_fertile_candidates.csv.gz"
                    if provider == "ken_fern"
                    else "pfaf_self_fertile_candidates.csv.gz"
                ),
                "source_file": _text(record.get("cache_file")),
                "acceptance_contract": (
                    "exact_wild_sexual_self_fertile_maps_to_sc_medium_v1"
                ),
            }
        )
    evidence = pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS)
    return evidence, audit


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(
        path,
        index=False,
        compression={"method": "gzip", "mtime": 0}
        if path.suffix == ".gz"
        else None,
    )


@app.command("build")
def build(
    ken_candidates_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    pfaf_candidates_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
) -> None:
    ken = pd.read_csv(ken_candidates_csv, dtype=str).fillna("")
    pfaf = pd.read_csv(pfaf_candidates_csv, dtype=str).fillna("")
    evidence, audit = build_package(ken, pfaf)
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "useful_plants_reviewed_direct_evidence.csv.gz"
    audit_path = output_dir / "useful_plants_provider_lineage_audit.csv.gz"
    _write_csv(evidence, evidence_path)
    _write_csv(audit, audit_path)
    provider_counts = (
        audit.loc[audit["decision"].eq("accept"), "provider"]
        .value_counts()
        .sort_index()
        .to_dict()
    )
    summary = {
        "contract": "useful_plants_reproductive_source_package_v1",
        "strict_direct": {
            "rows_before_lineage_dedup": len(evidence),
            "species_trait_after_lineage_dedup": int(
                evidence[["accepted_species", "trait_name", "source_lineage"]]
                .drop_duplicates()
                .shape[0]
            ),
            "species": int(evidence["accepted_species"].nunique()),
            "quality": "medium",
            "trait_name": "self_incompatibility",
            "normalized_value": "SC",
        },
        "provider_rows": provider_counts,
        "excluded": audit.loc[audit["decision"].eq("exclude"), "reason"]
        .value_counts()
        .sort_index()
        .to_dict(),
        "checks": {
            "autonomous_selfing_not_inferred": True,
            "apomixis_excluded": True,
            "cultivar_only_excluded": True,
            "provider_redistributions_share_lineage": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
        },
        "inputs": {
            str(ken_candidates_csv): _sha256(ken_candidates_csv),
            str(pfaf_candidates_csv): _sha256(pfaf_candidates_csv),
        },
        "artifacts": {
            evidence_path.name: _sha256(evidence_path),
            audit_path.name: _sha256(audit_path),
        },
    }
    summary_path = output_dir / "useful_plants_source_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    app()
