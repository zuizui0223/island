"""Recover normalized source families behind Database 1.0 citation lineages.

Rights/provenance audit only. This script never changes scientific trait values and
never grants redistribution rights. It resolves opaque ``citation:*`` provenance
against one or more immutable reviewed evidence tables, then maps only recovered
lineages to stable source families suitable for the existing rights policy.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd


MANUAL_RECEIPTS: dict[str, dict[str, str]] = {
    "citation:Chung-Chung-Oh-Epperson-2000-Heredity-85-490-497": {
        "source_family": "publication:chung_etal_2000_heredity",
        "source_provider": "reviewed_manual_literature",
        "source_url": "https://doi.org/10.1046/j.1365-2540.2000.00781.x",
        "source_citation": (
            "Chung, Chung, Oh & Epperson (2000). Spatial genetic structure in a "
            "Neolitsea sericea population (Lauraceae). Heredity 85:490-497."
        ),
    },
    "citation:Fischer-Rahelivololona-2007-Adansonia-29-269-315": {
        "source_family": "publication:fischer_rahelivololona_2007_adansonia",
        "source_provider": "MNHN Science Press",
        "source_url": "https://sciencepress.mnhn.fr/sites/default/files/articles/pdf/a2007n2a8.pdf",
        "source_citation": "Fischer & Rahelivololona (2007), Adansonia 29(2):269-315",
    },
    "citation:keil1975:pectis-cylindrica-autogamy": {
        "source_family": "publication:keil_1975_pectis_cylindrica",
        "source_provider": "reviewed_manual_literature",
        "source_url": "",
        "source_citation": (
            "Keil, D.J. (1975). Pectis cylindrica (Compositae) established as a member "
            "of the Texas flora and confirmed as a distinct species. Southwestern "
            "Naturalist 20:286-287."
        ),
    },
}


def _tokens(value: object) -> list[str]:
    return [token.strip() for token in str(value or "").split("|") if token.strip()]


def _slug(value: object) -> str:
    text = re.sub(r"[^a-z0-9]+", "_", str(value or "").casefold()).strip("_")
    return text or "unknown"


def _citation_shape(lineage: str) -> str:
    suffix = lineage.removeprefix("citation:")
    if re.fullmatch(r"[0-9a-f]{24}", suffix):
        return "hex24"
    if re.fullmatch(r"[0-9a-f]{20}", suffix):
        return "hex20"
    return "named"


def _family(provider: object, url: object, lineage: str) -> str:
    """Normalize recovered provenance without inventing redistribution rights.

    For GIFT, ``citation:*`` identifies the underlying cited flora/publication rather
    than GIFT itself, so keep it as a publication-citation family. The rights policy
    has no automatic grant for that family and therefore remains fail-closed.
    """
    p = str(provider or "").strip()
    folded = p.casefold()
    u = str(url or "").strip()
    host = urlparse(u).netloc.casefold().removeprefix("www.") if u else ""
    if "dryad" in folded or "dryad" in host or "doi.org/10.5061/dryad" in u.casefold():
        return "dataset:dryad"
    if folded == "flora_of_australia_official":
        return "provider_treatment:flora_of_australia"
    if folded.startswith("sanbi_eflora_south_africa"):
        return "provider_treatment:sanbi_eflora_south_africa"
    if folded == "proteus":
        return "database:proteus"
    if folded == "gift_v3_2_direct":
        return f"publication_citation:{lineage.removeprefix('citation:').casefold()}"
    if host:
        return f"domain:{host}"
    return f"provider:{_slug(p)}"


def _load_reviewed(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str).fillna("")
    required = {"source_lineage", "source_provider", "source_url", "source_citation"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"reviewed provenance table {path} missing columns: {sorted(missing)}")
    frame = frame.loc[frame["source_lineage"].str.startswith("citation:")].copy()
    frame["source_family"] = [
        _family(provider, url, lineage)
        for provider, url, lineage in zip(
            frame["source_provider"],
            frame["source_url"],
            frame["source_lineage"],
            strict=True,
        )
    ]
    frame["provenance_source"] = path.name
    return frame


def recover(
    species_axis_path: Path,
    reviewed_provenance_paths: tuple[Path, ...],
    output_dir: Path,
) -> dict[str, object]:
    if not reviewed_provenance_paths:
        raise ValueError("at least one reviewed provenance table is required")

    coverage = pd.read_csv(species_axis_path, dtype=str).fillna("")
    if "source_lineages" not in coverage.columns:
        raise ValueError("species-axis ledger missing source_lineages")
    target = sorted(
        {
            token
            for value in coverage["source_lineages"]
            for token in _tokens(value)
            if token.startswith("citation:")
        }
    )
    if not target:
        raise ValueError("no citation lineages in Database 1.0")

    candidates = pd.concat(
        [_load_reviewed(path) for path in reviewed_provenance_paths],
        ignore_index=True,
        sort=False,
    ).fillna("")
    candidates = candidates.loc[candidates["source_lineage"].isin(target)].copy()

    rows: list[dict[str, object]] = []
    for lineage in target:
        found = candidates.loc[candidates["source_lineage"].eq(lineage)].copy()
        families = sorted(set(found["source_family"]))
        providers = sorted(set(found["source_provider"]))
        urls = sorted({x for x in found["source_url"] if x})
        citations = sorted({x for x in found["source_citation"] if x})
        sources = sorted(set(found["provenance_source"]))
        provenance_method = "reviewed_provenance_table" if len(found) else ""

        receipt = MANUAL_RECEIPTS.get(lineage)
        if receipt is not None:
            manual_family = receipt["source_family"]
            if families and families != [manual_family]:
                raise ValueError(
                    f"manual receipt conflicts with reviewed provenance: {lineage} -> "
                    f"{families} vs {manual_family}"
                )
            families = [manual_family]
            providers = sorted(set(providers) | {receipt["source_provider"]})
            if receipt["source_url"]:
                urls = sorted(set(urls) | {receipt["source_url"]})
            citations = sorted(set(citations) | {receipt["source_citation"]})
            sources = sorted(set(sources) | {"committed_or_verified_manual_receipt"})
            provenance_method = "manual_receipt"

        if len(families) > 1:
            raise ValueError(
                f"citation lineage maps to conflicting source families: {lineage} -> {families}"
            )
        recovered = len(families) == 1
        rows.append(
            {
                "source_lineage": lineage,
                "citation_shape": _citation_shape(lineage),
                "source_family": families[0] if recovered else "",
                "source_provenance_recovered": recovered,
                "source_providers": "|".join(providers),
                "source_urls": "|".join(urls),
                "source_citations": " || ".join(citations),
                "provenance_sources": "|".join(sources),
                "provenance_method": provenance_method,
                "unresolved_reason": (
                    "" if recovered else "no_reviewed_or_manual_provenance_match"
                ),
            }
        )

    audit = pd.DataFrame(rows)
    recovered_rows = audit.loc[audit["source_provenance_recovered"]].copy()
    unresolved = audit.loc[~audit["source_provenance_recovered"]].copy()
    output_dir.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output_dir / "CITATION_PROVENANCE_AUDIT.csv", index=False)
    recovered_rows[["source_lineage", "source_family"]].drop_duplicates().to_csv(
        output_dir / "CITATION_LINEAGE_FAMILY_MAP.csv", index=False
    )
    unresolved.to_csv(output_dir / "CITATION_UNRESOLVED.csv", index=False)

    provider_summary = (
        recovered_rows.groupby("source_family", dropna=False)
        .agg(
            citation_lineages=("source_lineage", "nunique"),
            example_lineage=("source_lineage", "first"),
            example_provider=("source_providers", "first"),
        )
        .reset_index()
        .sort_values(["citation_lineages", "source_family"], ascending=[False, True])
    )
    provider_summary.to_csv(output_dir / "CITATION_SOURCE_FAMILY_SUMMARY.csv", index=False)

    unresolved_shapes = unresolved["citation_shape"].value_counts().to_dict()
    summary = {
        "contract": "chapter1_database_v1_citation_lineage_provenance_audit_v2",
        "distinct_citation_lineages": int(len(audit)),
        "source_provenance_recovered": int(len(recovered_rows)),
        "source_provenance_recovery_rate": len(recovered_rows) / len(audit),
        "unresolved_citation_lineages": int(len(unresolved)),
        "unresolved_hex24_lineages": int(unresolved_shapes.get("hex24", 0)),
        "unresolved_hex20_lineages": int(unresolved_shapes.get("hex20", 0)),
        "unresolved_named_lineages": int(unresolved_shapes.get("named", 0)),
        "reviewed_provenance_tables": [str(path) for path in reviewed_provenance_paths],
        "distinct_recovered_source_families": int(
            provider_summary["source_family"].nunique()
        ),
        "scientific_database_modified": False,
        "rights_granted": False,
    }
    (output_dir / "CITATION_PROVENANCE_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-axis", type=Path, required=True)
    parser.add_argument(
        "--reviewed-provenance",
        type=Path,
        action="append",
        required=True,
        help="Repeat for each immutable reviewed provenance CSV/CSV.GZ.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    recover(args.species_axis, tuple(args.reviewed_provenance), args.output_dir)


if __name__ == "__main__":
    main()
